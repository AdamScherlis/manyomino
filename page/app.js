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
  // it (page/zoom.js).  The camera sits at zoom s: level floor(s) is on screen at 2 to 4 px per
  // cell, and the next one fades in over it at half that.
  (function () {
    const Z0 = DATA.zoom, cv = $('zoom-canvas'); if (!Z0 || !cv) return;
    const G = 768, P0 = 2, AHEAD = 3, ctx = cv.getContext('2d');
    const bin = atob(Z0.bits), base = new Uint8Array(G * G);
    { const ox = (G - Z0.w) >> 1, oy = (G - Z0.h) >> 1;
      for (let y = 0; y < Z0.h; y++) for (let x = 0; x < Z0.w; x++) { const i = y * Z0.w + x; if (bin.charCodeAt(i >> 3) >> (i & 7) & 1) base[(x + ox) + (y + oy) * G] = 1; } }
    let seed = (Date.now() ^ (Math.random() * 1e9)) >>> 0, R = null, L = [], building = -1, s = -0.75, running = false, userPaused = false, visible = false, raf = 0, last = 0, ink = 0;
    function inkColour() {
      const m = /#?([0-9a-f]{6})/i.exec(css('--plate-ink')) || [0, 'E9EBF2'], v = parseInt(m[1], 16);
      return (255 << 24 | (v & 255) << 16 | (v >> 8 & 255) << 8 | v >> 16) >>> 0;              // ABGR
    }
    function picture(pic) {
      const c = document.createElement('canvas'); c.width = c.height = G;
      const x = c.getContext('2d'), im = x.createImageData(G, G), d = new Uint32Array(im.data.buffer);
      for (let i = 0; i < G * G; i++) if (pic[i]) d[i] = ink;
      x.putImageData(im, 0, 0); return c;
    }
    function level(pic, focus) { return { pic, focus, cv: picture(pic), from: [G / 2 + 1, G / 2 + 1], f0: 0 }; }   // the camera runs from `from` (at f0) to the focus (at 1)
    function reset() {
      ink = inkColour(); R = Refiner(G, seed, { sweeps: Z0.sweeps }); seed = (seed + 0x9E3779B9) >>> 0;
      // dive towards a randomly chosen cell of the first shape, not too near its edge
      const lim = R.limit - 8, c = G >> 1; let f = null;
      for (let t = 0; t < 400 && !f; t++) { const x = c + ((Math.random() * 2 - 1) * lim | 0), y = c + ((Math.random() * 2 - 1) * lim | 0); if (base[x + y * G]) f = [x, y]; }
      L = [level(base, f || nearestCell(base, G, c, c, lim))]; building = -1; s = -0.75; draw();
    }
    // make levels up to AHEAD beyond the one on screen, a few milliseconds at a time
    function build(ms) {
      const want = Math.max(0, Math.floor(s)) + AHEAD;
      if (building < 0) {
        if (L.length > want) return;
        const prev = L[L.length - 1]; R.start(prev.pic, prev.focus[0], prev.focus[1]); building = L.length;
      }
      if (R.step(ms)) {
        const pic = R.result(); L.push(level(pic, nearestCell(pic, G, G >> 1, G >> 1, 200) || [G >> 1, G >> 1])); building = -1;
        for (let k = 0; k < L.length - AHEAD - 3; k++) L[k] = null;                              // forget what is far behind
      }
    }
    // where the camera is, in the coordinates of level A, at fraction f of the way to the next level:
    // it closes on the focus so that the focus drifts to the middle of the screen at a steady pace
    function camera(A, f) {
      const w = f <= A.f0 ? 1 : ((1 - f) * Math.pow(2, -f)) / ((1 - A.f0) * Math.pow(2, -A.f0));
      return [A.focus[0] + 0.5 + (A.from[0] - A.focus[0] - 0.5) * w, A.focus[1] + 0.5 + (A.from[1] - A.focus[1] - 0.5) * w];
    }
    // head for the cell nearest a point on screen: the level on screen (or the one fading in) gets
    // a new focus and everything finer is thrown away
    function steer(px, py) {
      const dpr = size(), V = cv.width, k = Math.max(0, Math.floor(s)), f = Math.max(0, s - k), A = L[k]; if (!A) return;
      const [camx, camy] = camera(A, f), p = P0 * dpr * Math.pow(2, f);
      let x = camx + (px * dpr - V / 2) / p, y = camy + (py * dpr - V / 2) / p, j = k, lv = A;
      if (f > 0.5 && L[k + 1]) { x = (x - A.focus[0]) * 2 + G / 2; y = (y - A.focus[1]) * 2 + G / 2; j = k + 1; lv = L[j]; }
      const lim = R.limit - 4, c = G >> 1, hit = nearestCell(lv.pic, G, Math.round(x), Math.round(y), 60);
      if (!hit || Math.abs(hit[0] - c) > lim || Math.abs(hit[1] - c) > lim) return;
      if (j === k) { lv.from = [camx, camy]; lv.f0 = Math.min(f, 0.95); }
      lv.focus = hit; L.length = j + 1; building = -1;
      if (!running) draw();
    }
    const smooth = (x, a, b) => { const t = Math.min(1, Math.max(0, (x - a) / (b - a))); return t * t * (3 - 2 * t); };
    function size() {
      const dpr = Math.min(2, window.devicePixelRatio || 1), w = Math.round(cv.clientWidth * dpr);
      if (cv.width !== w) { cv.width = cv.height = w; }
      return dpr;
    }
    function draw() {
      const dpr = size(), V = cv.width, k = Math.max(0, Math.floor(s)), f = s - k, A = L[k];
      ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.clearRect(0, 0, V, V); ctx.imageSmoothingEnabled = false;
      if (!A) return;
      const c = G / 2, [camx, camy] = camera(A, f), p = P0 * dpr * Math.pow(2, f);
      const put = (lv, x, y, q, alpha) => { if (alpha <= 0) return; ctx.globalAlpha = alpha; ctx.drawImage(lv.cv, V / 2 - x * q, V / 2 - y * q, G * q, G * q); };
      const B = L[k + 1], fade = B ? smooth(f, 0.3, 0.65) : 0;
      put(A, camx, camy, p, 1 - (B ? smooth(f, 0.5, 0.9) : 0));
      if (B) put(B, (camx - A.focus[0]) * 2 + c, (camy - A.focus[1]) * 2 + c, p / 2, fade);
      ctx.globalAlpha = 1;
      const mag = Math.pow(2, Math.max(0, s)), cells = Z0.n * Math.pow(Z0.growth, Math.max(0, s));
      const big = v => v < 1e6 ? fmtNum(+v.toPrecision(2)) : (v / Math.pow(10, Math.floor(Math.log10(v)))).toFixed(1) + '\u00d710' + sup(Math.floor(Math.log10(v)));
      $('zoom-readout').innerHTML = s <= 0 ? 'the whole shape: <b>' + fmtNum(Z0.n) + '</b> cells' :
        'magnified <b>' + big(mag) + '\u00d7</b> \u00b7 at this resolution the whole shape would have about <b>' + big(cells) + '</b> cells';
    }
    const sup = n => String(n).split('').map(ch => '⁰¹²³⁴⁵⁶⁷⁸⁹'[+ch]).join('');
    function frame(t) {
      const dt = Math.min(0.1, (t - last) / 1000); last = t;
      s = Math.min(s + Math.pow(10, +$('zoom-speed').value) * dt, L.length - 1 + 0.3);    // wait at the fade-in point for a level not made yet
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
      visible = es[0].isIntersecting;
      if (!visible && running) setRunning(false); else if (visible && !running && !userPaused && !still) setRunning(true);
    }, { threshold: 0.3 }).observe($('zoomer'));
    window.__zoom = { get s() { return s; }, set s(v) { s = v; }, get levels() { return L.length; }, draw, build, setRunning };
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
