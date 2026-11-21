// For analysis/part_test.py: make the level-1 world from BALL.bin with the page's code and count how many of
// the labelled parts of its box are still separate within each radius (in world cells).
//   node analysis/part_run.js BALL.bin SEED [key=value ...]  ->  one JSON line
const fs = require('fs'), path = require('path');
const { Outgoer, outRad, OUT } = require(path.join(__dirname, '..', 'page', 'outward.js'));
const E = new WebAssembly.Instance(new WebAssembly.Module(fs.readFileSync(path.join(__dirname, '..', 'page', 'outward.wasm'))), {}).exports; E.init(12);
const raw = new Int32Array(fs.readFileSync(process.argv[2]).buffer.slice(0)), n = raw[0], x = new Int32Array(n), y = new Int32Array(n);
for (let i = 0; i < n; i++) { x[i] = raw[1 + 2 * i]; y[i] = raw[2 + 2 * i]; }
const opt = { J: 0, seed: +process.argv[3] }; process.argv.slice(4).forEach(a => { const [k, v] = a.split('='); opt[k] = +v; });
const O = Outgoer(E, opt); O.begin({ x, y }); O.next(); const w = O.world, B = OUT.B, key = (a, b) => (b + 4096) * 8192 + a + 4096, at = new Map();
for (let i = 0; i < w.x.length; i++) at.set(key(w.x[i], w.y[i]), i);
const out = [];
for (const R of [B, B + 1, B + 2, 3 * B / 2, 2 * B, 3 * B, 4 * B]) {
  const up = new Map(), find = a => { while (up.get(a) !== a) { up.set(a, up.get(up.get(a))); a = up.get(a); } return a; };
  for (let i = 0; i < w.x.length; i++) if (outRad(w.x[i], w.y[i]) < R) up.set(i, i);
  for (const i of up.keys()) for (const [dx, dy] of [[1, 0], [0, 1]]) { const j = at.get(key(w.x[i] + dx, w.y[i] + dy)); if (j === undefined || !up.has(j)) continue; if (w.lab[i] && w.lab[j] && w.lab[i] !== w.lab[j]) continue; up.set(find(i), find(j)); }
  const cls = new Set(); for (let i = 0; i < w.x.length; i++) if (w.lab[i]) cls.add(find(i)); out.push(cls.size);
}
console.log(JSON.stringify(out));
