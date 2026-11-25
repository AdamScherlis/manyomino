// For analysis/mass_model.py: the real cells of one shape and the worlds 1, 2 and 3 levels down (halved and thinned
// only, no sampler), written as int32: for each of the four, a count, then the x's, then the y's.
//   node analysis/mass_dump.js BALL.bin SEED OUT.bin
const fs = require('fs'), path = require('path');
const { Outgoer } = require(path.join(__dirname, '..', 'page', 'outward.js'));
const E = new WebAssembly.Instance(new WebAssembly.Module(fs.readFileSync(path.join(__dirname, '..', 'page', 'outward.wasm'))), {}).exports; E.init(12);
const raw = new Int32Array(fs.readFileSync(process.argv[2]).buffer.slice(0)), n = raw[0], x = new Int32Array(n), y = new Int32Array(n);
for (let i = 0; i < n; i++) { x[i] = raw[1 + 2 * i]; y[i] = raw[2 + 2 * i]; }
const O = Outgoer(E, { J: 1, seed: +process.argv[3], N: 0, S: 0 }); O.begin({ x, y }); const ws = [{ x, y }];
for (let L = 1; L <= 3; L++) { O.next(); ws.push({ x: O.world.x.slice(), y: O.world.y.slice() }); }
const out = new Int32Array(ws.reduce((a, w) => a + 1 + 2 * w.x.length, 0)); let p = 0;
for (const w of ws) { out[p++] = w.x.length; out.set(w.x, p); p += w.x.length; out.set(w.y, p); p += w.y.length; }
fs.writeFileSync(process.argv[4], Buffer.from(out.buffer));
