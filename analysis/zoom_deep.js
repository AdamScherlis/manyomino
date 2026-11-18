// The browser's zoom engine (page/zoom.js) run many levels deep under node: does anything drift?
// At each level: per-cell statistics of the picture, and the growth of the cell count in the
// region that is kept.  Writes analysis/zoom_deep.json.
//   node analysis/zoom_deep.js [LEVELS=400] [SWEEPS=32] [SEED=7] [G=1536]
const fs = require('fs'), path = require('path'); const { Refiner, nearestCell } = require(path.join(__dirname, '..', 'page', 'zoom.js'));
const G = +process.argv[5] || 1536, LEVELS = +process.argv[2] || 400, T = +process.argv[3] || 32, seed = +process.argv[4] || 7;
const buf = fs.readFileSync(path.join(__dirname, '..', 'prod', 'snap_A_n16000_p0_0001.bin')); const a = new Int32Array(buf.buffer, buf.byteOffset, buf.length / 4); const n0 = a[0];
let xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9;
for (let i = 0; i < n0; i++) { xmin = Math.min(xmin, a[1 + 2 * i]); xmax = Math.max(xmax, a[1 + 2 * i]); ymin = Math.min(ymin, a[2 + 2 * i]); ymax = Math.max(ymax, a[2 + 2 * i]); }
let pic = new Uint8Array(G * G); const ox = (G - (xmax - xmin + 1)) >> 1, oy = (G - (ymax - ymin + 1)) >> 1;
for (let i = 0; i < n0; i++) pic[(a[1 + 2 * i] - xmin + ox) + (a[2 + 2 * i] - ymin + oy) * G] = 1;
function stats(pic, lo, hi) {           // cells with lo <= x,y < hi
  let n = 0, per = 0, sq = 0, d1 = 0, d3 = 0, st = 0; const o = (x, y) => pic[x + y * G];
  for (let y = lo; y < hi; y++) for (let x = lo; x < hi; x++) {
    if (!o(x, y)) { if (o(x + 1, y) || o(x - 1, y) || o(x, y + 1) || o(x, y - 1)) per++; continue; }
    n++; const E = o(x + 1, y), Wn = o(x - 1, y), N = o(x, y + 1), Sn = o(x, y - 1), d = E + Wn + N + Sn;
    if (d === 1) d1++; if (d === 3) d3++; if (d === 2 && ((E && Wn) || (N && Sn))) st++; if (E && N && o(x + 1, y + 1)) sq++;
  }
  return { n, per: per / n, sq: sq / n, d1: d1 / n, d3: d3 / n, st: st / n };
}
const R = Refiner(G, seed, { sweeps: T });
let f = nearestCell(pic, G, G >> 1, G >> 1, 300); let tot = { per: 0, sq: 0, d1: 0, d3: 0, st: 0, w: 0 }, ratios = [];
const t0 = Date.now(); let moves = 0; const rows = [];
for (let k = 1; k <= LEVELS; k++) {
  // cells of the coarse picture in the square of side G/4 about the focus, which becomes the middle half of the next
  let nc = 0; for (let y = f[1] - G / 8; y < f[1] + G / 8; y++) for (let x = f[0] - G / 8; x < f[0] + G / 8; x++) nc += pic[x + y * G];
  R.start(pic, f[0], f[1]); const nstart = R.n; while (!R.step(1000)); moves += R.moves;
  const fine = R.result(); const st = stats(fine, 64, G - 64);
  let nf = 0; for (let y = G / 4; y < 3 * G / 4; y++) for (let x = G / 4; x < 3 * G / 4; x++) nf += fine[x + y * G];
  ratios.push([nc, nf]); rows.push([k, st.n, st.per, st.sq, st.d1, st.d3, st.st, nf / nc]);
  if (k <= 3 || k % 50 === 0) console.log(`level ${k}: ${st.n} cells in view (movable ${nstart} -> ${R.n}); perim ${st.per.toFixed(4)} squares ${st.sq.toFixed(4)} deg1 ${st.d1.toFixed(4)} deg3 ${st.d3.toFixed(4)} straight ${st.st.toFixed(4)}; growth ${(nf / nc).toFixed(3)}`);
  if (k > 1) { for (const key of ['per', 'sq', 'd1', 'd3', 'st']) tot[key] += st[key] * st.n; tot.w += st.n; }
  pic = fine; f = nearestCell(pic, G, G >> 1, G >> 1, 300);
  if (!f) { console.log('lost the shape at level', k); break; }
}
const el = (Date.now() - t0) / 1000;
console.log(`T=${T}: ${LEVELS} levels in ${el.toFixed(1)} s (${(moves / el / 1e6).toFixed(2)} M moves/s, ${(el / LEVELS).toFixed(2)} s per level)`);
console.log(`weighted over levels 2+: perim ${(tot.per / tot.w).toFixed(4)} (real 1.1950)  squares ${(tot.sq / tot.w).toFixed(4)} (0.0706)  deg1 ${(tot.d1 / tot.w).toFixed(4)} (0.1692)  deg3 ${(tot.d3 / tot.w).toFixed(4)} (0.2585)  straight ${(tot.st / tot.w).toFixed(4)} (0.2957)`);
const sc = ratios.reduce((s, r) => s + r[0], 0), sf = ratios.reduce((s, r) => s + r[1], 0);
console.log(`growth of the cell count in the same region, summed over levels: ${(sf / sc).toFixed(4)}; mean of ratios ${(ratios.reduce((s, r) => s + r[1] / r[0], 0) / ratios.length).toFixed(4)}`);
const half = a => { const m = a.reduce((p, r) => p + r[1], 0); return ['per', 'sq', 'd1', 'd3', 'st'].map((_, j) => a.reduce((p, r) => p + r[2 + j] * r[1], 0) / m); };
const body = rows.slice(1), mid = body.length >> 1;
fs.writeFileSync(path.join(__dirname, 'zoom_deep.json'), JSON.stringify({ levels: LEVELS, sweeps: T, seed, G, columns: ['level', 'cells', 'perim', 'squares', 'deg1', 'deg3', 'straight', 'growth'],
  mean: half(body), first_half: half(body.slice(0, mid)), second_half: half(body.slice(mid)), growth: sf / sc, seconds_per_level: el / LEVELS, rows: rows.map(r => r.map((v, i) => i < 2 ? v : +v.toFixed(4))) }));
