// Worlds only, a long way out, several times over: does the amount of shape near the root drift from level to level?
//   node analysis/out_drift.js BALL.bin LEVELS SEED [key=value ...]   ->  one JSON line per level: cells of the world in its box, in its ring, by layer, parts kept apart
const fs = require('fs'), path = require('path');
const { Outgoer, OUT, outRad } = require(path.join(__dirname, '..', 'page', 'outward.js'));
const E = new WebAssembly.Instance(new WebAssembly.Module(fs.readFileSync(path.join(__dirname, '..', 'page', 'outward.wasm'))), {}).exports; E.init(12);
const raw = new Int32Array(fs.readFileSync(process.argv[2]).buffer.slice(0)), n = raw[0], x = new Int32Array(n), y = new Int32Array(n);
for (let i = 0; i < n; i++) { x[i] = raw[1 + 2 * i]; y[i] = raw[2 + 2 * i]; }
const opt = { J: 0, seed: +process.argv[4] }; process.argv.slice(5).forEach(a => { const [k, v] = a.split('='); opt[k] = +v; });
const O = Outgoer(E, opt), B = OUT.B; O.begin({ x, y });
for (let L = 1; L <= +process.argv[3]; L++) {
  const p = O.next(), w = O.world, lay = new Array(4 * B).fill(0); let box = 0, ring = 0, far = 0;
  for (let i = 0; i < w.x.length; i++) { const r = outRad(w.x[i], w.y[i]); if (r < 4 * B) lay[r]++; if (r < B) box++; else if (r < 2 * B) ring++; else if (r < 4 * B) far++; }
  console.log(JSON.stringify({ level: L, box, ring, far, parts: p.info.parts, lay }));
}
