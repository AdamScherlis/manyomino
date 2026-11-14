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
})();
