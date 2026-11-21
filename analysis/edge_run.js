// For analysis/edge_test.py: halve the world in BALL.bin in several ways, with the page's own code, and
// describe what lands in the box.
//   node analysis/edge_run.js BALL.bin SEED [VARIANTS]  ->  one JSON line: { variant: [numbers...] }
const fs = require('fs'), path = require('path');
const { outRad, Outgoer } = require(path.join(__dirname, '..', 'page', 'outward.js'));
const E = new WebAssembly.Instance(new WebAssembly.Module(fs.readFileSync(path.join(__dirname, '..', 'page', 'outward.wasm'))), {}).exports;
E.init(12);
const raw = new Int32Array(fs.readFileSync(process.argv[2]).buffer.slice(0)), n = raw[0], x = new Int32Array(n), y = new Int32Array(n), B = 8;
for (let i = 0; i < n; i++) { x[i] = raw[1 + 2 * i]; y[i] = raw[2 + 2 * i]; }
// name: the engine's settings for the halving
const VAR = process.argv[4] ? JSON.parse(process.argv[4]) : { 'as on the page': {}, 'free 256': { RULE: 0, COVER: 0, W: 1, TC: 256, PARTS: 0 } };
const key = (a, b) => (b + 4096) * 8192 + a + 4096, out = {};
for (const v in VAR) {
  const O = Outgoer(E, Object.assign({ J: 1, seed: +process.argv[3] }, VAR[v])); O.begin({ x, y });
  const or = new Set(); for (let i = 0; i < n; i++) or.add(key(x[i] >> 1, y[i] >> 1));
  O.halve({});
  const m = E.count(), X = new Int32Array(E.memory.buffer, E.xs(), m), Y = new Int32Array(E.memory.buffer, E.ys(), m), Lb = new Int32Array(E.memory.buffer, E.partlabels(), m), S = new Map();
  let inside = 0, edge = 0, out1 = 0, bare = 0, far = 0;
  for (let i = 0; i < m; i++) { const r = outRad(X[i], Y[i]); S.set(key(X[i], Y[i]), i); if (r < B) inside++; if (r === B) out1++; }
  // the edge ring in order, for runs; edge cells where the finer shape had none (bare), or none within one site (far)
  const ring = []; for (let t = -B; t < B - 1; t++) ring.push([t, -B]); for (let t = -B; t < B - 1; t++) ring.push([B - 1, t]); for (let t = B - 1; t > -B; t--) ring.push([t, B - 1]); for (let t = B - 1; t > -B; t--) ring.push([-B, t]);
  const o = ring.map(([a, b]) => S.has(key(a, b)) ? 1 : 0); let runs = 0;
  ring.forEach(([a, b], i) => { if (!o[i]) return; edge++; if (!o[(i + ring.length - 1) % ring.length]) runs++;
    if (!or.has(key(a, b))) { bare++; let near = false; for (let da = -1; da <= 1; da++) for (let db = -1; db <= 1; db++) if (or.has(key(a + da, b + db))) near = true; if (!near) far++; } });
  if (edge === ring.length) runs = 1;
  let lost = 0; ring.forEach(([a, b]) => { if (!or.has(key(a, b))) return; for (let da = -1; da <= 1; da++) for (let db = -1; db <= 1; db++) if (S.has(key(a + da, b + db))) return; lost++; });
  // pieces of what is in the box: touching cells are joined unless both carry labels and the labels differ
  const seen = new Set(); let parts = 0, touching = 0;
  for (let i = 0; i < m; i++) { if (outRad(X[i], Y[i]) >= B || seen.has(i)) continue; parts++; const st = [i]; seen.add(i);
    while (st.length) { const a = st.pop(); for (const [da, db] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) { const j = S.get(key(X[a] + da, Y[a] + db)); if (j === undefined || outRad(X[j], Y[j]) >= B) continue; if (VAR[v].PARTS !== 0 && Lb[a] !== Lb[j]) { touching++; continue; } if (!seen.has(j)) { seen.add(j); st.push(j); } } } }
  out[v] = [n / m, inside, edge, runs, parts, out1, bare, far, lost, touching / 2];
}
console.log(JSON.stringify(out));
