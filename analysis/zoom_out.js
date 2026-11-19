// The outward generator (page/zoom.js, Outward) against the truth, and run a long way out.
//
// 1. Truth.  In real samples of 1,024,000 cells (prod/snap_*_n1024000_*.bin), take the window
//    of G sites around a cell, hand it to the generator, and compare what it puts in the windows
//    twice and four times as wide with what the sample really has there.  A generated cell of
//    level j stands for GROW^j real cells where the level is fully relaxed, and for fewer where
//    it has had fewer sweeps (the same correction the page uses for brightness).
// 2. Stationarity.  From one window, many batches outward: the number of cells in the fully
//    relaxed ring of each new level, and its per-cell statistics, should not drift.
//
// The truth test uses a small window (256 sites by default), so that the windows two and four
// times as wide are still well inside a sample some 7,000 cells across: around cells near the
// middle of a sample the count grows faster than around a typical cell of a much larger shape,
// and near its edge more slowly.  Roots are any cells of the samples.
//
//   node analysis/zoom_out.js [ROOTS_PER_SAMPLE=25] [BATCHES=30] [G=256] [G_FAR=1024] [OUT=zoom_out.json]
const fs = require('fs'), path = require('path');
const { Outward, zoomRate } = require(path.join(__dirname, '..', 'page', 'zoom.js'));
const NR = +process.argv[2] || 25, NB = +process.argv[3] || 30, GFAR = +process.argv[5] || 1024, T = 16, GROW = 2.9605;
let G = +process.argv[4] || 256, C = G >> 1;
const MT = [[0, 4.0], [0.25, 3.6757], [0.5, 3.4798], [0.75, 3.3576], [1, 3.2733], [1.5, 3.1704], [2, 3.1114], [3, 3.0533], [4, 3.023], [6, 2.9953], [8, 2.9796], [12, 2.9708], [16, 2.964]];
const left = t => { let k = 0; while (k < MT.length - 2 && MT[k + 1][0] < t) k++; const a = MT[k], b = MT[k + 1]; return t >= b[0] ? b[1] : a[1] + (b[1] - a[1]) * (t - a[0]) / (b[0] - a[0]); };
let wring = null; const setG = g => { G = g; C = G >> 1; wring = new Float64Array(C + 1); for (let r = 0; r <= C; r++) wring[r] = left(T) / left(T * zoomRate((r + 0.5) / C)); }; setG(G);
let seed = 12345; const rnd = () => { seed = (seed + 0x6D2B79F5) | 0; let t = Math.imul(seed ^ (seed >>> 15), 1 | seed); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
const ring = (x, y) => Math.max(x < C ? C - 1 - x : x - C, y < C ? C - 1 - y : y - C);
// weighted cells of a level outside its held middle half, in units of that level's fully relaxed cells
function outer(pic) { let m = 0; for (let y = 0; y < G; y++) for (let x = 0; x < G; x++) if (pic[x + y * G]) { const r = ring(x, y); if (r >= C >> 1) m += wring[r]; } return m; }
function ringStats(pic) {      // the fully relaxed ring, 0.5 < r < 0.72
  let n = 0, per = 0, sq = 0, d1 = 0; const o = (x, y) => pic[x + y * G] ? 1 : 0, lo = (C >> 1) + 4, hi = Math.floor(0.72 * C) - 6;
  for (let y = 1; y < G - 1; y++) for (let x = 1; x < G - 1; x++) {
    const r = ring(x, y); if (r < lo || r > hi) continue;
    if (!o(x, y)) { if (o(x + 1, y) || o(x - 1, y) || o(x, y + 1) || o(x, y - 1)) per++; continue; }
    n++; const d = o(x + 1, y) + o(x - 1, y) + o(x, y + 1) + o(x, y - 1); if (d === 1) d1++; if (o(x + 1, y) && o(x, y + 1) && o(x + 1, y + 1)) sq++;
  }
  return { n, per: per / n, sq: sq / n, d1: d1 / n };
}
const files = fs.readdirSync(path.join(__dirname, '..', 'prod')).filter(f => /^snap_._n1024000_.*\.bin$/.test(f)).sort();
const rows = [], fits = []; let last = null; const t0 = Date.now();
// rings of the first new level, as shares of the way from its middle to its edge: generated against real, cell for cell
const EDGES = [0.5, 0.55, 0.6, 0.72, 0.85, 1.0], ringGen = new Float64Array(EDGES.length - 1), ringReal = new Float64Array(EDGES.length - 1);
const which = r => { for (let q = 0; q < EDGES.length - 1; q++) if (r >= EDGES[q] * C && r < EDGES[q + 1] * C) return q; return -1; };
for (const f of files.filter((_, i) => i % 2 === 0)) {
  setG(+process.argv[4] || 256);
  const buf = fs.readFileSync(path.join(__dirname, '..', 'prod', f)); const a = new Int32Array(buf.buffer, buf.byteOffset, buf.length / 4), n = a[0];
  let mx = 0, my = 0; for (let i = 0; i < n; i++) { mx += a[1 + 2 * i]; my += a[2 + 2 * i]; } mx /= n; my /= n;
  for (let k = 0; k < NR; k++) {
    const i = (rnd() * n) | 0;
    const rx = a[1 + 2 * i], ry = a[2 + 2 * i], F = new Uint8Array(G * G), real = [0, 0, 0];
    for (let j = 0; j < n; j++) {
      const dx = a[1 + 2 * j] - rx, dy = a[2 + 2 * j] - ry, d = Math.max(dx < 0 ? -dx - 1 : dx, dy < 0 ? -dy - 1 : dy);
      if (d < 2 * G) real[2]++; if (d < G) real[1]++; if (d < C) { real[0]++; F[(dx + C) + (dy + C) * G] = 255; } else if (d < G) { const q = which(d / 2); if (q >= 0) ringReal[q] += 1 / GROW; }
    }
    const O = Outward(G, (rnd() * 4294967296) >>> 0, { sweeps: T }), lv = []; let top = F, fitSum = 0;
    while (lv.length < 2) { O.begin(top); let r; while (!(r = O.step(1e9))); lv.push(...r); top = r[r.length - 1]; if (lv.length === r.length) fits.push(O.fit[0]); }
    const g1 = real[0] + GROW * outer(lv[0]), g2 = g1 + GROW * GROW * outer(lv[1]);
    for (let y = 0; y < G; y++) for (let x = 0; x < G; x++) if (lv[0][x + y * G]) { const r = ring(x, y), q = which(r); if (q >= 0) ringGen[q] += wring[r]; }
    rows.push({ file: f, real, gen: [real[0], Math.round(g1), Math.round(g2)] }); last = { rx, ry, a, n };
    if (k === 0) console.log(f, 'real', real.join(' '), 'generated', Math.round(g1), Math.round(g2), ((Date.now() - t0) / 1000).toFixed(0) + ' s');
  }
}
const lg = (v, i) => Math.log(v[i] / v[i - 1]), ms = v => { const m = v.reduce((p, q) => p + q, 0) / v.length; return [m, Math.sqrt(v.reduce((p, q) => p + (q - m) ** 2, 0) / (v.length - 1))]; };
const corr = (x, y) => { const [mx, sx] = ms(x), [my, sy] = ms(y); return x.reduce((p, q, i) => p + (q - mx) * (y[i] - my), 0) / ((x.length - 1) * sx * sy); };
const sum = {}, se = (v) => ms(v)[1] / Math.sqrt(v.length);
const rings = EDGES.slice(0, -1).map((e, q) => [e, EDGES[q + 1], ringGen[q] / ringReal[q]]);
console.log('generated over real by ring of the first new level:', rings.map(r => `${r[0]}-${r[1]}: ${r[2].toFixed(2)}`).join('  '));
for (const i of [1, 2]) {
  const r = rows.map(q => lg(q.real, i)), g = rows.map(q => lg(q.gen, i)), dens = rows.map(q => Math.log(q.real[0]));
  sum['step' + i] = { real: ms(r).concat(se(r)).map(v => v / Math.LN2), gen: ms(g).concat(se(g)).map(v => v / Math.LN2), corr_with_inside_real: corr(dens, r), corr_with_inside_gen: corr(dens, g), corr_real_gen: corr(r, g) };
  console.log(`cells in the window ${2 ** i}x as wide over the one inside it, log2: real ${(ms(r)[0] / Math.LN2).toFixed(3)} ± ${(se(r) / Math.LN2).toFixed(3)} sd ${(ms(r)[1] / Math.LN2).toFixed(3)}; generated ${(ms(g)[0] / Math.LN2).toFixed(3)} ± ${(se(g) / Math.LN2).toFixed(3)} sd ${(ms(g)[1] / Math.LN2).toFixed(3)}; ` +
              `correlation with the density inside: real ${corr(dens, r).toFixed(2)}, generated ${corr(dens, g).toFixed(2)}; real with generated ${corr(r, g).toFixed(2)}`);
}
// stationarity: a long way out from a window of the page's size around the last root
setG(GFAR); let F = new Uint8Array(G * G);
for (let j = 0; j < last.n; j++) { const x = last.a[1 + 2 * j] - last.rx + C, y = last.a[2 + 2 * j] - last.ry + C; if (x >= 0 && y >= 0 && x < G && y < G) F[x + y * G] = 255; }
const O = Outward(G, 99, { sweeps: T }), far = [], tf = Date.now();
let links = 0;
for (let b = 0; b < NB; b++) { O.begin(F); let lv; while (!(lv = O.step(1e9))); for (const p of lv) far.push(ringStats(p)); F = lv[lv.length - 1]; }
console.log('links', O.joined.join(' cells '), 'in', far.length, 'levels; blocks of 64 differing at the best match, mean over roots:', (fits.reduce((p, q) => p + q, 0) / fits.length).toFixed(1));
const perLevel = (Date.now() - tf) / 1000 / far.length;
const half = v => ['n', 'per', 'sq', 'd1'].map(k => v.reduce((p, q) => p + q[k], 0) / v.length);
console.log(`${far.length} levels outward: cells in the relaxed ring, perimeter, 2x2 blocks and tips per cell: first half ${half(far.slice(0, far.length >> 1)).map(v => v.toFixed(4)).join(' ')}; second half ${half(far.slice(far.length >> 1)).map(v => v.toFixed(4)).join(' ')}`);
const OUT = process.argv[6] || 'zoom_out.json';
fs.writeFileSync(path.join(__dirname, OUT), JSON.stringify({ fit: fits.reduce((p, q) => p + q, 0) / fits.length, rings, links: O.joined, per_batch: O.J, G: +process.argv[4] || 256, G_far: GFAR, seconds_per_level: perLevel, sweeps: T, roots: rows.length, rows, summary: sum, levels: far.length, first_half: half(far.slice(0, far.length >> 1)), second_half: half(far.slice(far.length >> 1)),
  ring: far.map(q => [q.n, +q.per.toFixed(4), +q.sq.toFixed(4), +q.d1.toFixed(4)]) }));
