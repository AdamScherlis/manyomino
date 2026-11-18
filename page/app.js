// Page wiring.  DATA is injected by analysis/build_page.py.
(function () {
  const $ = id => document.getElementById(id);
  const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();

  // ---------------------------------------------------------------- hero viewer
  // the sample ships as base64 text (one nibble per cell underneath)
  fetch(DATA.heroUrl).then(r => r.text()).then(txt => {
    const bin = atob(txt.trim()), u8 = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
    const poly = decodePoly(u8.buffer);
    $('hero-loading').hidden = true;
    PolyViewer($('hero'), poly);
  }).catch(() => { $('hero-loading').textContent = 'could not load the sample'; });

  // ---------------------------------------------------------------- move diagrams
  function diagram(el, W, H, layers, extra) {
    const u = 22, pad = 6;
    let s = `<svg viewBox="0 0 ${W * u + 2 * pad} ${H * u + 2 * pad}" role="img">`;
    const rect = (c, cls) => `<rect class="cell ${cls}" x="${pad + c[0] * u}" y="${pad + (H - 1 - c[1]) * u}" width="${u}" height="${u}" rx="2"/>`;
    layers.forEach(([cells, cls]) => cells.forEach(c => { s += rect(c, cls); }));
    const mid = c => [pad + (c[0] + .5) * u, pad + (H - .5 - c[1]) * u];
    s += (extra ? extra(mid, u) : '') + '</svg>';
    el.innerHTML = s;
  }
  const arrowDef = '<defs><marker id="ah" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="var(--mut)"/></marker></defs>';
  (function () {
    const E = DATA.example, B = new Set(E.B);
    // relocate: one tip cell moves
    const from = E.reloc.from, to = E.reloc.to, cut = E.reloc.cut;
    const rest = E.before.filter(c => !(c[0] === from[0] && c[1] === from[1]) && !(c[0] === cut[0] && c[1] === cut[1]));
    diagram($('dg-relocate'), E.w, E.h, [[rest, ''], [[cut], 'c'], [[from], 'gone'], [[to], 'b']], (mid) => {
      const a = mid(from), b = mid(to);
      return arrowDef + `<path class="arrow" marker-end="url(#ah)" d="M${a[0] - 8} ${a[1] - 8} Q ${(a[0] + b[0]) / 2 - 50} ${(a[1] + b[1]) / 2 + 10} ${b[0] + 6} ${b[1] + 13}"/>`;
    });
    // transplant: branch B moves from c to c'
    const A = E.before.filter((c, i) => !B.has(i) && i !== E.c && i !== E.a);
    const oldB = E.before.filter((c, i) => B.has(i)), newB = E.after.filter((c, i) => B.has(i));
    diagram($('dg-transplant'), E.w, E.h, [[A, ''], [[E.before[E.c], E.before[E.a]], 'c'], [oldB, 'gone'], [newB, 'b']], (mid, u) => {
      const c = mid(E.before[E.c]), a = mid(E.before[E.a]);
      const p = mid(oldB[0]), q = mid(newB[0]);
      return arrowDef + `<text class="dgl" x="${c[0]}" y="${c[1] + 4.5}" text-anchor="middle">c</text><text class="dgl" x="${a[0]}" y="${a[1] + 4.5}" text-anchor="middle">c′</text>` +
        `<path class="arrow" marker-end="url(#ah)" d="M${c[0] + u * .2} ${c[1] - u * 1.8} Q ${(c[0] + a[0]) / 2} ${c[1] - u * 5.2} ${a[0] + u * .3} ${a[1] - u * 1.6}"/>`;
    });
  })();

  // ---------------------------------------------------------------- endless zoom
  // Levels are G x G pictures, each twice the resolution of the last and centred on a cell of
  // it (page/zoom.js).  Every point on screen is drawn from the finest level whose window
  // contains it, so the picture is a nest of squares, finest in the middle, and no level is ever
  // exchanged for another at a size where the eye could tell: a new one comes in at the centre
  // with cells a tenth of a pixel wide, and at its edge a level is identical to the next coarser
  // one.  Finer levels have less ink per unit area (the shape's dimension is below 2), so each is
  // brightened by the matching amount, which keeps the brightness even across the nest and
  // steady in time.  Drawing is one WebGL2 shader, with a plain canvas fallback; levels are made
  // in a worker where the browser allows one.
  (function () {
    const Z0 = DATA.zoom, cv = $('zoom-canvas'); if (!Z0 || !cv) return;
    const small = cv.clientWidth <= 420;
    let gl = null; try { gl = cv.getContext('webgl2', { alpha: false, antialias: false }); } catch (e) {}
    if (gl && (gl.getParameter(gl.MAX_TEXTURE_SIZE) < 2048 || gl.getParameter(gl.MAX_ARRAY_TEXTURE_LAYERS) < 8)) gl = null;
    const ctx = gl ? null : cv.getContext('2d');
    const G = gl && !small ? 2048 : 1024, C = G / 2, LIMIT = zoomLimit(G), LAYERS = 8, MAXDRAW = 6, AHEAD = 4;
    const P0 = 1.15 * (Math.max(240, cv.clientWidth) / 2) / C;          // px per cell of level k when the zoom s equals k: its window then just covers the view
    const PMIN = gl ? 0.1 : 0.45, DIM = 2 - 1 / 0.6412;               // a level comes in when its cells reach PMIN px; ink per area goes as (cell size)^DIM
    const bin = atob(Z0.bits), base = new Uint8Array(G * G);
    { const ox = (G - Z0.w) >> 1, oy = (G - Z0.h) >> 1;
      for (let y = 0; y < Z0.h; y++) for (let x = 0; x < Z0.w; x++) { const i = y * Z0.w + x; if (bin.charCodeAt(i >> 3) >> (i & 7) & 1) base[(x + ox) + (y + oy) * G] = 255; } }
    let seed = (Date.now() ^ (Math.random() * 1e9)) >>> 0, L = [], s = 0, s0 = 0, kref = 0, cam = [C, C], running = false, userPaused = false, raf = 0, last = 0;
    let R = null, building = -1;                                  // main-thread refinement, used only without a worker
    let worker = null, useWorker = true, busy = false, job = 0, sent = 0;
    const opts = { sweeps: Z0.sweeps };
    const workerSrc = 'const ZOOM=' + JSON.stringify(ZOOM) + ';' + zoomRate.toString() + ';' + zoomLimit.toString() + ';' + Refiner.toString() + ';' +
      'let R=null,key=-1;onmessage=e=>{const m=e.data;if(!R||key!==m.seed){R=Refiner(m.G,m.seed,m.opts);key=m.seed;}' +
      'R.start(m.pic,m.fx,m.fy);while(!R.step(1e9));const pic=R.result();postMessage({id:m.id,pic},[pic.buffer]);};';
    function dropWorker(forGood) { if (worker) worker.terminate(); worker = null; busy = false; job++; if (forGood) useWorker = false; }
    function startWorker() {
      try {
        worker = new Worker(URL.createObjectURL(new Blob([workerSrc], { type: 'text/javascript' })));
        worker.onmessage = e => { if (e.data.id !== job) return; busy = false; add(e.data.pic); };
        worker.onerror = () => dropWorker(true);
      } catch (e) { worker = null; useWorker = false; }
    }
    const hex = (v, d) => { const m = /#?([0-9a-f]{6})/i.exec(css(v)) || [0, d], n = parseInt(m[1], 16); return [n >> 16, n >> 8 & 255, n & 255]; };

    // ---- drawing with WebGL2: one pass over the screen, walking the nest from the finest level outwards
    let prog = null, U = null, lutTex = null;
    if (gl) {
      const sh = (type, src) => { const o = gl.createShader(type); gl.shaderSource(o, src); gl.compileShader(o); return o; };
      prog = gl.createProgram();
      gl.attachShader(prog, sh(gl.VERTEX_SHADER, '#version 300 es\nvoid main(){vec2 p=vec2(float((gl_VertexID<<1)&2),float(gl_VertexID&2));gl_Position=vec4(p*2.0-1.0,0.0,1.0);}'));
      gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, `#version 300 es
precision highp float; precision highp sampler2DArray;
uniform sampler2DArray u_tex; uniform sampler2D u_lut;
uniform vec4 u_lv[${MAXDRAW}];        // camera in the level's cells (x, y), cells per device pixel, opacity
uniform vec4 u_meta[${MAXDRAW}];      // texture layer, brightness, whether the level thins out towards its edge
uniform int u_n; uniform vec2 u_res; uniform float u_G; uniform vec3 u_bg, u_ink;
out vec4 o;
void main() {
  vec2 d = gl_FragCoord.xy - 0.5 * u_res; d.y = -d.y;
  float lum = 0.0, rem = 1.0;
  for (int i = 0; i < ${MAXDRAW}; i++) {
    if (i >= u_n) break;
    float z = u_lv[i].z; vec2 uv = (u_lv[i].xy + d * z) / u_G;
    float r = 2.0 * max(abs(uv.x - 0.5), abs(uv.y - 0.5));
    float m = (1.0 - smoothstep(0.985, 1.0, r)) * u_lv[i].w;
    float cov = textureGrad(u_tex, vec3(uv, u_meta[i].x), vec2(z / u_G, 0.0), vec2(0.0, z / u_G)).r;
    float gf = mix(1.0, textureLod(u_lut, vec2(r, 0.5), 0.0).r, u_meta[i].z);
    lum += rem * m * u_meta[i].y * gf * cov; rem *= 1.0 - m;
  }
  o = vec4(pow(mix(u_bg, u_ink, clamp(lum, 0.0, 1.0)), vec3(1.0 / 2.2)), 1.0);
}`));
      gl.linkProgram(prog);
      if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) { console.error('zoom shader: ' + gl.getProgramInfoLog(prog)); }
      gl.useProgram(prog); gl.bindVertexArray(gl.createVertexArray());
      U = {}; ['u_tex', 'u_lut', 'u_lv', 'u_meta', 'u_n', 'u_res', 'u_G', 'u_bg', 'u_ink'].forEach(nm => { U[nm] = gl.getUniformLocation(prog, nm); });
      gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
      gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D_ARRAY, gl.createTexture());
      gl.texStorage3D(gl.TEXTURE_2D_ARRAY, 7, gl.R8, G, G, LAYERS);
      gl.texParameteri(gl.TEXTURE_2D_ARRAY, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR); gl.texParameteri(gl.TEXTURE_2D_ARRAY, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D_ARRAY, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE); gl.texParameteri(gl.TEXTURE_2D_ARRAY, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      // how much of the doubled shape is left at each distance from the middle of a level, as a brightness correction:
      // where a level has had few sweeps it still has more cells than it will end with
      const mt = Z0.mt, left = t => { let k = 0; while (k < mt.length - 2 && mt[k + 1][0] < t) k++; const a = mt[k], b = mt[k + 1]; return t >= b[0] ? b[1] : a[1] + (b[1] - a[1]) * (t - a[0]) / (b[0] - a[0]); };
      const lut = new Uint8Array(256); for (let i = 0; i < 256; i++) lut[i] = Math.round(255 * left(Z0.sweeps) / left(Z0.sweeps * zoomRate((i + 0.5) / 256)));
      gl.activeTexture(gl.TEXTURE1); lutTex = gl.createTexture(); gl.bindTexture(gl.TEXTURE_2D, lutTex);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, 256, 1, 0, gl.RED, gl.UNSIGNED_BYTE, lut);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      gl.uniform1i(U.u_tex, 0); gl.uniform1i(U.u_lut, 1); gl.uniform1f(U.u_G, G);
    }
    function show(lv, j) {                                        // hand a level's picture to whatever draws it
      if (gl) { gl.activeTexture(gl.TEXTURE0); gl.texSubImage3D(gl.TEXTURE_2D_ARRAY, 0, 0, 0, j % LAYERS, G, G, 1, gl.RED, gl.UNSIGNED_BYTE, lv.pic); gl.generateMipmap(gl.TEXTURE_2D_ARRAY); return; }
      const c = document.createElement('canvas'); c.width = c.height = G;
      const x = c.getContext('2d'), im = x.createImageData(G, G), d = new Uint32Array(im.data.buffer), k = hex('--plate-ink', 'E9EBF2'), ink = (255 << 24 | k[2] << 16 | k[1] << 8 | k[0]) >>> 0;
      for (let i = 0; i < G * G; i++) if (lv.pic[i]) d[i] = ink;
      x.putImageData(im, 0, 0); lv.cv = c;
    }
    function add(pic) {
      const lv = { pic, focus: nearestCell(pic, G, C, C, 300) || [C, C], born: running ? performance.now() : -1e9, taper: 1 };
      show(lv, L.length); L.push(lv); if (!running) draw();
    }
    function forget(from) { for (let i = from; i < L.length; i++) if (L[i] && L[i].cv) L[i].cv.width = 0; L.length = from; building = -1; if (busy) dropWorker(false); }
    function reset() {
      seed = (seed + 0x9E3779B9) >>> 0; forget(0); R = null;
      // dive towards a randomly chosen cell of the first shape
      let f = null;
      for (let t = 0; t < 400 && !f; t++) { const x = C + ((Math.random() * 2 - 1) * LIMIT | 0), y = C + ((Math.random() * 2 - 1) * LIMIT | 0); if (base[x + y * G]) f = [x, y]; }
      const lv = { pic: base, focus: f || nearestCell(base, G, C, C, LIMIT), born: -1e9, taper: 0 };
      show(lv, 0); L = [lv]; kref = 0; cam = [C, C];
      s = s0 = Math.log2(0.92 * Math.max(240, cv.clientWidth) / (Math.max(Z0.w, Z0.h) * P0)); draw();
    }
    // make levels up to AHEAD beyond the zoom
    function build(ms) {
      if (L.length - 1 >= Math.floor(s) + AHEAD) return;
      const prev = L[L.length - 1];
      if (useWorker) {
        if (busy) { if (performance.now() - sent > 30000) dropWorker(true); return; }
        if (!worker) startWorker();
        if (worker) { busy = true; sent = performance.now(); worker.postMessage({ id: ++job, G, seed, opts, pic: prev.pic, fx: prev.focus[0], fy: prev.focus[1] }); return; }
      }
      if (building < 0) { if (!R) R = Refiner(G, seed, opts); R.start(prev.pic, prev.focus[0], prev.focus[1]); building = L.length; }
      if (R.step(ms)) { building = -1; add(R.result()); }
    }
    // a point given in the cells of level a, in the cells of level b
    function conv(p, a, b) {
      let x = p[0], y = p[1];
      for (; a < b; a++) { x = 2 * (x - L[a].focus[0]) + C; y = 2 * (y - L[a].focus[1]) + C; }
      for (; a > b; a--) { x = L[a - 1].focus[0] + (x - C) / 2; y = L[a - 1].focus[1] + (y - C) / 2; }
      return [x, y];
    }
    // keep the camera in the cells of the level the zoom has reached (or the finest there is)
    function rebase() {
      const k = Math.max(L.findIndex(v => v), Math.min(Math.floor(Math.max(0, s)), L.length - 1));
      if (k !== kref) { cam = conv(cam, kref, k); kref = k; }
    }
    // where the zoom is heading: the middle of the finest level's chosen cell
    function target() { const D = L.length - 1, f = L[D].focus; return conv([f[0] + 0.5, f[1] + 0.5], D, kref); }
    const smooth = (x, a, b) => { const t = Math.min(1, Math.max(0, (x - a) / (b - a))); return t * t * (3 - 2 * t); };
    function size() {
      const dpr = Math.min(2, window.devicePixelRatio || 1), w = Math.round(cv.clientWidth * dpr);
      if (cv.width !== w) { cv.width = cv.height = w; }
      return dpr;
    }
    const px = j => P0 * Math.pow(2, s - j);                                  // CSS px per cell of level j
    function nest() {                                                          // the levels to draw, finest first
      rebase(); const out = [], now = performance.now();
      for (let j = L.length - 1; j >= Math.max(0, kref - 1) && out.length < MAXDRAW; j--) {
        if (!L[j]) break;
        const p = px(j), a = (j === 0 ? 1 : smooth(p, 0.8 * PMIN, PMIN)) * Math.min(1, (now - L[j].born) / 500);
        if (a > 0) out.push({ j, p, a, c: conv(cam, kref, j) });
      }
      return out;
    }
    function draw() {
      const dpr = size(), V = cv.width, N = nest();
      if (gl) {
        const lv = new Float32Array(4 * MAXDRAW), meta = new Float32Array(4 * MAXDRAW), bg = hex('--plate', '0D0F16'), ink = hex('--plate-ink', 'E9EBF2'), lin = c => c.map(v => Math.pow(v / 255, 2.2));
        N.forEach((e, i) => { lv.set([e.c[0], e.c[1], 1 / (e.p * dpr), e.a], 4 * i); meta.set([e.j % LAYERS, Math.min(3, Z0.gain * Math.pow(e.p, -DIM)), L[e.j].taper, 0], 4 * i); });
        gl.viewport(0, 0, V, V); gl.uniform4fv(U.u_lv, lv); gl.uniform4fv(U.u_meta, meta); gl.uniform1i(U.u_n, N.length); gl.uniform2f(U.u_res, V, V);
        gl.uniform3fv(U.u_bg, lin(bg)); gl.uniform3fv(U.u_ink, lin(ink)); gl.drawArrays(gl.TRIANGLES, 0, 3);
      } else {
        ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.clearRect(0, 0, V, V); ctx.imageSmoothingEnabled = true; ctx.imageSmoothingQuality = 'high';
        for (let i = N.length - 1; i >= 0; i--) {                                // coarsest first; a finer level replaces what is under its window
          const e = N[i], q = e.p * dpr, x = V / 2 - e.c[0] * q, y = V / 2 - e.c[1] * q;
          if (i < N.length - 1) ctx.clearRect(x + q, y + q, (G - 2) * q, (G - 2) * q);
          ctx.globalAlpha = e.a; ctx.drawImage(L[e.j].cv, x, y, G * q, G * q);
        }
        ctx.globalAlpha = 1;
      }
      const mag = Math.pow(2, Math.max(0, s - s0)), cells = Z0.n * Math.pow(Z0.growth, Math.max(0, s - s0));
      const big = v => v < 1e6 ? fmtNum(+v.toPrecision(2)) : (v / Math.pow(10, Math.floor(Math.log10(v)))).toFixed(1) + '×10' + sup(Math.floor(Math.log10(v)));
      $('zoom-readout').innerHTML = s <= s0 ? 'the whole shape: <b>' + fmtNum(Z0.n) + '</b> cells' :
        'magnified <b>' + big(mag) + '×</b> · at this resolution the whole shape would have about <b>' + big(cells) + '</b> cells';
    }
    const sup = n => String(n).split('').map(ch => '⁰¹²³⁴⁵⁶⁷⁸⁹'[+ch]).join('');
    // head for the cell nearest a point on screen: the finest level that can be re-centred there
    // gets a new focus, and everything finer is thrown away and made again
    function steer(x, y) {
      rebase(); const V = cv.clientWidth;
      for (let j = L.length - 1; j >= kref; j--) {
        const p = px(j), c = conv(cam, kref, j); let u = c[0] + (x - V / 2) / p, v = c[1] + (y - V / 2) / p;
        if (j > kref && Math.max(Math.abs(u - C), Math.abs(v - C)) > LIMIT) continue;
        u = Math.min(C + LIMIT, Math.max(C - LIMIT, u)); v = Math.min(C + LIMIT, Math.max(C - LIMIT, v));
        const hit = nearestCell(L[j].pic, G, Math.round(u), Math.round(v), Math.min(400, Math.ceil(40 / p)));
        if (!hit || Math.max(Math.abs(hit[0] - C), Math.abs(hit[1] - C)) > LIMIT) continue;
        L[j].focus = hit; forget(j + 1); if (!running) draw();
        return;
      }
    }
    function frame(t) {
      const dt = Math.min(0.1, (t - last) / 1000); last = t;
      // full speed while levels are well ahead of the zoom, slowing to a stop if they fall behind
      const lead = L.length - 1 - s;
      s += Math.pow(10, +$('zoom-speed').value) * dt * Math.min(1, Math.max(0, (lead - 0.4) / 1.2));
      rebase(); const T = target(), w = 1 - Math.exp(-dt / 0.9); cam = [cam[0] + (T[0] - cam[0]) * w, cam[1] + (T[1] - cam[1]) * w];
      for (let k = 0; k < kref - 1; k++) if (L[k]) { if (L[k].cv) L[k].cv.width = 0; L[k] = null; }          // forget what is behind
      build(7); draw();
      if (running) raf = requestAnimationFrame(frame);
    }
    function setRunning(r) { running = r; $('zoom-run').textContent = r ? 'Pause' : 'Zoom'; cancelAnimationFrame(raf); if (r) { last = performance.now(); raf = requestAnimationFrame(frame); } }
    $('zoom-run').addEventListener('click', () => { userPaused = running; setRunning(!running); });
    $('zoom-new').addEventListener('click', () => { reset(); if (!running && !userPaused) setRunning(true); });
    cv.addEventListener('click', ev => { const r = cv.getBoundingClientRect(); steer(ev.clientX - r.left, ev.clientY - r.top); });
    new ResizeObserver(() => { if (!running) draw(); }).observe(cv);
    reset();
    const still = window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches;
    new IntersectionObserver(es => {
      const visible = es[0].isIntersecting;
      if (!visible && running) setRunning(false); else if (visible && !running && !userPaused && !still) setRunning(true);
    }, { threshold: 0.3 }).observe($('zoomer'));
    window.__zoom = { get s() { return s; }, set s(v) { s = v; }, get levels() { return L.length; }, get worker() { return useWorker && !!worker; }, G, gl: !!gl, draw, build, setRunning };
  })();

  // ---------------------------------------------------------------- live race
  (function () {
    const N = DATA.demo.n, lanes = [];
    const mk = (canvasId, statId, pcp, color) => ({ canvas: $(canvasId), stat: $(statId), pcp, color, S: null, moves: 0, trace: [], next: 1000, view: null });
    lanes.push(mk('race-a', 'race-a-stat', 0, '--s2'), mk('race-b', 'race-b-stat', 0.4, '--s1'));
    let running = false, raf = 0, seed = 20261114;
    function reset() {
      const cells = []; for (let i = 0; i < N; i++) cells.push([i % 40, (i / 40) | 0]);
      lanes.forEach((L, k) => { L.S = Sampler(cells, 9, seed + 17 * k); L.moves = 0; L.trace = [[1, Math.sqrt(L.S.rg2())]]; L.next = 1000; L.view = null; });
      seed += 101; drawAll(); drawTrace();
    }
    function drawLane(L) {
      const c = L.canvas, dpr = Math.min(2, devicePixelRatio || 1), r = c.getBoundingClientRect();
      const w = Math.round(r.width * dpr), h = Math.round(r.height * dpr);
      if (c.width !== w || c.height !== h) { c.width = w; c.height = h; }
      const ctx = c.getContext('2d'), S = L.S, n = S.n, now = performance.now();
      let x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9;
      for (let i = 0; i < n; i++) { const x = S.X[i], y = S.Y[i]; if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y; }
      const tgt = { cx: (x0 + x1 + 1) / 2, cy: (y0 + y1 + 1) / 2, span: Math.max(x1 - x0, y1 - y0) + 5 };
      if (!L.view) L.view = tgt; else for (const k in tgt) L.view[k] += (tgt[k] - L.view[k]) * 0.12;
      const s = Math.min(w, h) / L.view.span, D = S.depths();
      let md = 1; for (let i = 0; i < n; i++) if (D[i] > md) md = D[i];
      ctx.fillStyle = 'rgb(' + PLATE.join(',') + ')'; ctx.fillRect(0, 0, w, h);
      const sz = Math.max(1, s - (s > 5 ? 1 : 0));
      for (let i = 0; i < n; i++) {
        const col = rampColor(D[i] / md), f = S.flash[i] ? Math.max(0, 1 - (now - S.flash[i]) / 450) : 0;
        ctx.fillStyle = f > 0 ? `rgb(${Math.round(col[0] + (255 - col[0]) * f)},${Math.round(col[1] + (255 - col[1]) * f)},${Math.round(col[2] + (255 - col[2]) * f)})` : `rgb(${col[0]},${col[1]},${col[2]})`;
        ctx.fillRect(w / 2 + (S.X[i] - L.view.cx) * s, h / 2 - (S.Y[i] + 1 - L.view.cy) * s, sz, sz);
      }
      const rg = Math.sqrt(S.rg2());
      L.stat.innerHTML = `<b>${L.moves.toLocaleString('en-US')}</b> moves &middot; radius <b>${rg.toFixed(1)}</b>`;
    }
    function drawAll() { lanes.forEach(drawLane); }
    function drawTrace() {
      const el = $('race-trace'), W = Math.max(280, el.clientWidth), H = 190, m = { l: 40, r: 14, t: 8, b: 36 };
      const xmax = Math.max(5, Math.ceil(Math.log10(Math.max(1e5, lanes[0].moves)) + 0.02)), xmin = 3;
      const eq = DATA.demo, ymax = Math.max(46, eq.rg + 2.2 * eq.rg_sd), ymin = 10;
      const X = v => m.l + (Math.log10(Math.max(v, 1e3)) - xmin) / (xmax - xmin) * (W - m.l - m.r), Y = v => H - m.b - (v - ymin) / (ymax - ymin) * (H - m.t - m.b);
      let s = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" role="img" aria-label="Radius of gyration against moves for the two live walks">`;
      s += `<rect x="${m.l}" width="${W - m.l - m.r}" y="${Y(eq.rg + eq.rg_sd)}" height="${Y(eq.rg - eq.rg_sd) - Y(eq.rg + eq.rg_sd)}" fill="var(--s3)" opacity=".16"/>`;
      for (let y = 10; y <= ymax; y += 10) s += `<line class="grid" x1="${m.l}" x2="${W - m.r}" y1="${Y(y)}" y2="${Y(y)}"/><text class="tick" x="${m.l - 7}" y="${Y(y) + 4}" text-anchor="end">${y}</text>`;
      for (let e = xmin; e <= xmax; e++) s += `<line class="grid" x1="${X(10 ** e)}" x2="${X(10 ** e)}" y1="${m.t}" y2="${H - m.b}"/><text class="tick" x="${X(10 ** e)}" y="${H - m.b + 16}" text-anchor="middle">${e >= 6 ? `10<tspan dy="-5" font-size="8.5">${e}</tspan>` : fmtNum(10 ** e)}</text>`;
      s += `<text class="alabel" x="${(m.l + W - m.r) / 2}" y="${H - 4}" text-anchor="middle">moves made</text>`;
      s += `<text class="dlabel" x="${m.l + 8}" y="${Y(eq.rg + eq.rg_sd) - 5}">uniform samples</text>`;
      lanes.forEach(L => { s += `<polyline fill="none" stroke="var(${L.color})" stroke-width="2" stroke-linejoin="round" points="${L.trace.map(p => X(p[0]).toFixed(1) + ',' + Y(Math.min(p[1], ymax)).toFixed(1)).join(' ')}"/>`; });
      el.innerHTML = s + '</svg>';
    }
    let last = 0, lastTrace = 0;
    function frame(t) {
      if (!running) return;
      const dt = Math.min(0.1, (t - last) / 1000 || 0.016); last = t;
      const rate = Math.pow(10, +$('race-speed').value), k = Math.max(1, Math.round(rate * dt));
      lanes.forEach(L => {
        const S = L.S, t0 = performance.now();
        S.setClock(t0);
        let i = 0;
        for (; i < k; i++) {
          if (L.pcp && S.rnd() < L.pcp) S.transplant(); else S.relocate();
          if ((i & 1023) === 1023 && performance.now() - t0 > 7) { i++; break; }
        }
        L.moves += i;
        if (L.moves >= L.next) { L.trace.push([L.moves, Math.sqrt(S.rg2())]); L.next = Math.max(L.moves + 200, L.moves * 1.02); if (L.trace.length > 1500) L.trace = L.trace.filter((_, j) => j % 2 === 0); }
      });
      drawAll();
      if (t - lastTrace > 250) { drawTrace(); lastTrace = t; }
      raf = requestAnimationFrame(frame);
    }
    function setRunning(r) { running = r; $('race-run').textContent = r ? 'Pause' : 'Run'; if (r) { last = performance.now(); raf = requestAnimationFrame(frame); } else cancelAnimationFrame(raf); }
    $('race-run').addEventListener('click', () => { userPaused = running; setRunning(!running); });
    $('race-reset').addEventListener('click', reset);
    const so = () => { $('race-speed-out').textContent = fmtNum(Math.round(Math.pow(10, +$('race-speed').value) / 1000) * 1000) + ' moves/s'; };
    $('race-speed').addEventListener('input', so); so();
    let userPaused = false;
    reset();
    new ResizeObserver(() => { drawAll(); drawTrace(); }).observe($('race-trace'));
    // run only while on screen; never auto-start for reduced motion
    const reduce = matchMedia('(prefers-reduced-motion: reduce)').matches;
    new IntersectionObserver(es => {
      const vis = es[0].isIntersecting;
      if (!vis && running) setRunning(false);
      else if (vis && !running && !userPaused && !reduce) setRunning(true);
    }, { threshold: 0.25 }).observe($('race'));
  })();

  // ---------------------------------------------------------------- charts
  const S1 = 'var(--s1)', S2 = 'var(--s2)', S3 = 'var(--s3)';
  const pow = (a, z) => n => a * Math.pow(n, z);
  (function () {
    const fo = DATA.fit_old, fn = DATA.fit_new, fO = pow(fo.a, fo.z), fN = pow(fn.a, fn.z);
    const nO = DATA.reloc.map(r => r.n), nN = DATA.sizes.map(r => r.n);
    chart($('ch-tau'), {
      xlog: true, ylog: true, xlabel: 'cells in the polyomino', ylabel: 'autocorrelation time (moves)', xname: 'cells', yname: 'τ (moves)',
      xdomain: [40, 2.2e6], ydomain: [1e3, 1e16], rightPad: 18, ysig: 3,
      aria: 'Autocorrelation time against size for the two walks, log-log',
      series: [
        { color: S2, dash: true, line: [[Math.max(...nO), fO(Math.max(...nO))], [2e6, fO(2e6)]] },
        { color: S1, dash: true, line: [[100, fN(100)], [Math.min(...nN), fN(Math.min(...nN))]] },
        { name: 'one cell at a time', color: S2, pts: DATA.reloc.map(r => [r.n, r.tau_rg2, r.tau_rg2_err]), line: [[Math.min(...nO), fO(Math.min(...nO))], [Math.max(...nO), fO(Math.max(...nO))]], label: 'n' + sup(fo.z.toFixed(2)), labelAt: [3e4, fO(3e4)], labelDx: -10, labelDy: -6, labelAnchor: 'end' },
        { name: 'with transplants', color: S1, pts: DATA.sizes.map(r => [r.n, r.tau_rg2, r.tau_rg2_err]), line: [[Math.min(...nN), fN(Math.min(...nN))], [Math.max(...nN), fN(Math.max(...nN))]], label: 'n' + sup(fn.z.toFixed(2)), labelAt: [3e4, fN(3e4)], labelDx: 6, labelDy: 18 },
      ],
    });
  })();
  function sup(s) { return `<tspan dy="-5" font-size="9">${s}</tspan>`; }
  (function () {
    const eq = DATA.square_eq;
    chart($('ch-square'), {
      xlog: true, ylog: false, xlabel: 'moves made', ylabel: 'radius of gyration', xname: 'moves', yname: 'radius',
      xdomain: [1e5, Math.max(...DATA.square_old.map(p => p[0])) * 1.15], ydomain: [0, Math.max(420, eq.rg + 2 * eq.rg_sd)], noTable: true, rightPad: 16,
      aria: 'Radius of gyration against moves, starting from a square of 30,000 cells',
      bands: [{ y0: eq.rg - eq.rg_sd, y1: eq.rg + eq.rg_sd, color: S3, label: 'uniform samples', labelAt: eq.rg + eq.rg_sd, left: true }],
      series: [
        { name: 'one cell at a time', color: S2, line: DATA.square_old.filter(p => p[0] >= 1e5).map(p => [p[0], Math.sqrt(p[1])]) },
        { name: 'with transplants', color: S1, line: DATA.square_new.filter(p => p[0] >= 1e5).map(p => [p[0], Math.sqrt(p[1])]) },
      ],
    });
  })();
  (function () {
    const f = DATA.fit_nu, nu0 = 0.6412;
    const mine = DATA.sizes.map(r => [r.n, r.rg, r.rg_err]);
    chart($('ch-rg'), {
      xlog: true, ylog: true, xlabel: 'cells in the polyomino', ylabel: 'radius of gyration', xname: 'cells', yname: 'mean radius', ysig: 5,
      aria: 'Mean radius of gyration against size, log-log, with power-law fit',
      series: [
        { color: S1, line: [[mine[0][0], f.a * Math.pow(mine[0][0], f.nu)], [mine[mine.length - 1][0], f.a * Math.pow(mine[mine.length - 1][0], f.nu)]] },
        { color: S1, pts: mine, label: 'slope ' + f.nu.toFixed(4), labelAt: [3000, f.a * Math.pow(3000, f.nu)], labelDx: 12, labelDy: 14 },
      ],
    });
    chart($('ch-resid'), {
      xlog: true, xlabel: 'cells in the polyomino', ylabel: 'radius ÷ n^0.6412', xname: 'cells', yname: 'radius ÷ n^0.6412', ysig: 4,
      aria: 'Mean radius divided by n to the 0.6412 against size, this work and the earlier gallery',
      series: [
        { name: 'this sampler', color: S1, pts: DATA.sizes.filter(r => r.n >= 250).map(r => [r.n, r.rg / Math.pow(r.n, nu0), r.rg_err / Math.pow(r.n, nu0)]) },
        { name: 'earlier gallery (one-cell walk)', color: S2, pts: DATA.gallery.map(g => [g[0], g[1] / Math.pow(g[0], nu0)]) },
      ],
    });
  })();
  (function () {
    const pts = DATA.sizes.filter(r => r.n >= 250).map(r => [r.n, r.lam, r.lam_err]);
    chart($('ch-lam'), {
      xlog: true, xlabel: 'cells in the polyomino', ylabel: 'estimate of λ', xname: 'cells', yname: 'λ estimate', ysig: 6, height: 300,
      yfmt: v => v.toFixed(3), aria: 'Estimate of Klarner’s constant from each size, with the accepted value',
      bands: [{ y0: DATA.lambda_ref - 1e-7, y1: DATA.lambda_ref + 1e-7, mid: DATA.lambda_ref, color: S3, label: 'accepted value 4.0625696', labelAt: DATA.lambda_ref, left: true }],
      series: [{ color: S1, pts }],
    });
  })();
  // ---------------------------------------------------------------- convergence checks
  if (DATA.conv && $('ch-z')) (function () {
    const C = DATA.conv, rows = C.tests, lim = 6;
    // one dot per observable, one row per test; a little vertical jitter so dots do not hide each other
    const jit = i => ((i * 0.6180339887) % 1 - 0.5) * 0.56;
    const mk = kind => rows.flatMap((t, r) => (t.kind || 'chain') === kind ? t.z.map((z, i) => [Math.max(-lim, Math.min(lim, z[1])), rows.length - r + jit(i), 0, z[0], z[1], t.label.replace('|', ', ')]) : []);
    chart($('ch-z'), {
      xlabel: 'difference, in standard errors', ylabel: '', xdomain: [-lim - 0.3, lim + 0.3], ydomain: [0.4, rows.length + 1.05],
      height: 64 + 40 * rows.length, leftPad: 124, rightPad: 14, noTable: true,
      yticks: rows.map((t, r) => [rows.length - r, t.label]),
      vbands: [{ x0: -2, x1: 2, color: S3, label: 'within 2' }],
      aria: 'For each test, the difference between the chain and the reference for every measurement, in standard errors',
      tip: p => `<b>${p[5]}</b><br>${p[3]}: ${p[4] >= 0 ? '+' : ''}${p[4].toFixed(2)}`,
      series: [{ name: 'the walk against a trusted reference', color: S1, r: 3, pts: mk('chain') },
               { name: 'a walk made wrong on purpose', color: S2, r: 3, pts: mk('control') },
               { name: 'reference that fails its own check', color: 'var(--mut)', r: 3, pts: mk('strained') }].filter(se => se.pts.length),
    });
    const sl = C.slow;
    chart($('ch-slow'), {
      xlog: true, xlabel: 'cells in the polyomino', ylabel: 'moves per cell to forget', xname: 'cells', yname: 'autocorrelation time ÷ n', ysig: 3,
      ydomain: [0, Math.max(...sl.map(r => r.mode + r.mode_err), ...sl.map(r => r.single)) * 1.15], height: 320,
      aria: 'Autocorrelation time in moves per cell against size: the radius, the slowest single observable, and the slowest combination found',
      series: [
        { name: 'squared radius', color: S1, pts: sl.map(r => [r.n, r.rg2, r.rg2_err]), line: sl.map(r => [r.n, r.rg2]) },
        { name: 'slowest single measurement', color: S2, pts: sl.map(r => [r.n, r.single, r.single_err]), line: sl.map(r => [r.n, r.single]) },
        { name: 'slowest combination found', color: S3, pts: sl.map(r => [r.n, r.mode, r.mode_err]), line: sl.map(r => [r.n, r.mode]) },
      ],
    });
    if ($('fig-ring') && DATA.clump) (function () {
      const R = DATA.clump.ring, c = 10, d = (cells, cls) => `<path class="${cls}" d="${cells.map(p => `M${p[0] * c + 1} ${(R.h - 1 - p[1]) * c + 1}h${c - 2}v${c - 2}h${2 - c}z`).join('')}"/>`;
      $('fig-ring').innerHTML = `<svg viewBox="0 0 ${R.w * c} ${R.h * c}" role="img" aria-label="A clump of ${R.cells} cells: a thin closed ring about ${R.w - 12} cells across, drawn in colour, with nearby cells of the shape in grey">${d(R.near, 'near')}${d(R.ring, 'ring')}</svg>`;
    })();
    if ($('ch-clump-tau') && DATA.clump) (function () {
      const K = DATA.clump, ramp = ['color-mix(in srgb, var(--s1) 50%, var(--panel))', 'color-mix(in srgb, var(--s1) 75%, var(--panel))', S1];
      const nm = r => r.n.toLocaleString('en-US') + '-cell shapes';
      chart($('ch-clump-dens'), {
        xlog: true, ylog: true, xlabel: 'cells in the clump, at least', ylabel: 'clumps per cell', xname: 'clump size at least', yname: 'per cell', ysig: 3, height: 320,
        aria: 'Number of clumps of at least a given size per cell, against clump size, at three shape sizes; the three curves coincide',
        series: K.base.map((r, i) => ({ name: nm(r), color: ramp[i], r: 3, pts: r.dens, line: r.dens.map(p => [p[0], p[1]]) })),
      });
      chart($('ch-clump-tau'), {
        xlog: true, ylog: true, xlabel: 'cells in the clump, at least', ylabel: 'moves per cell to forget', xname: 'clump size at least', yname: 'moves per cell', ysig: 3, height: 320,
        aria: 'Autocorrelation time of the number of clumps of at least a given size, against clump size, at three shape sizes; it rises with clump size',
        series: K.base.map((r, i) => ({ name: nm(r), color: ramp[i], r: 3, pts: r.tau })),
      });
      if ($('ch-clump-fix') && K.fix) chart($('ch-clump-fix'), {
        xlog: true, ylog: true, xlabel: 'cells in the clump, at least', ylabel: 'moves per cell to forget', xname: 'clump size at least', yname: 'moves per cell', ysig: 3, height: 320,
        aria: 'Autocorrelation time of the number of clumps of at least a given size, for the walk as used and with larger boxes drawn more often',
        series: [{ name: 'the walk as used', color: S1, r: 3, pts: K.pooled, line: K.pooled.map(p => [p[0], p[1]]) },
                 { name: 'large boxes drawn more often', color: S2, r: 3, pts: K.fix, line: K.fix.map(p => [p[0], p[1]]) }],
      });
    })();
    if ($('ch-zoom-mass') && DATA.zoom && DATA.zoom.mass) (function () {
      const Z = DATA.zoom, ramp = ['color-mix(in srgb, var(--s1) 60%, var(--panel))', S1];
      chart($('ch-zoom-mass'), {
        xlog: true, xlabel: 'rounds of the walk', ylabel: 'cells ÷ original', xname: 'rounds', yname: 'cells ÷ original', ysig: 5, height: 320, ydomain: [2.9, 3.3],
        aria: 'Number of cells in a doubled shape relative to the original, against rounds of the local walk; it falls from 4 and settles near 2.95',
        bands: [{ y0: Z.pred[0] - Z.pred[1], y1: Z.pred[0] + Z.pred[1], color: S3, label: '2 to the power 1/ν', labelAt: Z.pred[0] + Z.pred[1] }],
        series: Z.mass.map((m, i) => ({ name: m.name, color: ramp[i % 2], r: 3, pts: m.pts, line: m.pts.map(p => [p[0], p[1]]) })),
      });
      const clip = z => Math.max(-6, Math.min(6, z));
      if ($('ch-zoom-z') && Z.test) chart($('ch-zoom-z'), {
        xlog: true, xlabel: 'rounds per level', ylabel: 'standard errors from sampled', xname: 'rounds per level', yname: 'standard errors', ysig: 3, height: 320, noTable: true,
        xdomain: [Math.min(...Z.test.map(t => t.T)) / 1.6, Math.max(...Z.test.map(t => t.T)) * 1.6], ydomain: [-6, 6],
        aria: 'Sixty measurements of shapes refined four times, as differences from sampled shapes in standard errors, for several numbers of rounds per level',
        bands: [{ y0: -2, y1: 2, color: S3, label: 'within two', labelAt: 2 }],
        series: [{ name: 'the shape at large (48 measurements)', color: S1, r: 2.6, pts: [].concat(...Z.test.map(t => t.rest.map((z, i) => [t.T * Math.pow(1.22, (i % 9 - 4) / 4), clip(z)]))) },
                 { name: 'counts of small features per cell (12)', color: S2, r: 2.6, pts: [].concat(...Z.test.map(t => t.local.map((z, i) => [t.T * Math.pow(1.22, (i % 5 - 2) / 2), clip(z)]))) }],
      });
    })();
    if ($('ch-fossil') && DATA.fossil) (function () {
      const F = DATA.fossil, ramp = ['color-mix(in srgb, var(--s1) 55%, var(--panel))', 'color-mix(in srgb, var(--s1) 78%, var(--panel))', S1];
      const cs = F.curves.slice(-3);
      chart($('ch-fossil'), {
        ylog: true, xlabel: 'moves made per cell', ylabel: 'cells in largest intact piece', xname: 'moves per cell', yname: 'cells', ysig: 4,
        xdomain: [0, 615], ydomain: [1.5, Math.max(...cs.map(c => c.pts[0][1])) * 1.6], height: 340, noTable: true,
        aria: 'Size of the largest piece of the starting shape that is still intact, against moves per cell, for several sizes',
        bands: [{ y0: F.floor[0] - 0.4, y1: F.floor[1] + 0.4, color: S3, label: 'coincidence', labelAt: F.floor[1] + 0.4, left: true }],
        tip: (p, se) => `<b>${se.name}</b><br>after ${p[0]} moves per cell<br>largest intact piece: ${p[1].toLocaleString('en-US')} cells`,
        series: cs.map((c, i) => ({ name: c.n.toLocaleString('en-US') + ' cells', color: ramp[i + 3 - cs.length], line: c.pts }))
          .concat(F.hero.length ? [{ name: F.hero_n.toLocaleString('en-US') + ' cells', color: S2, pts: F.hero.map(p => [p[0], p[1]]) }] : []),
      });
    })();
  })();
})();
