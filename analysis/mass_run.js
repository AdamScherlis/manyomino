// For analysis/mass_test.py: a real shape and the worlds made from it by halving and thinning alone (no sampler),
// counted in pixels 8 real cells wide: the real cells, and the cells of the world 1, 2 and 3 levels down (which are
// 2, 4 and 8 real cells wide, so 2, 1 and 0 stages finer than the pixel).
//   node analysis/mass_run.js BALL.bin SEED
const fs = require('fs'), path = require('path');
const { Outgoer } = require(path.join(__dirname, '..', 'page', 'outward.js'));
const E = new WebAssembly.Instance(new WebAssembly.Module(fs.readFileSync(path.join(__dirname, '..', 'page', 'outward.wasm'))), {}).exports; E.init(12);
const raw = new Int32Array(fs.readFileSync(process.argv[2]).buffer.slice(0)), n = raw[0], x = new Int32Array(n), y = new Int32Array(n);
for (let i = 0; i < n; i++) { x[i] = raw[1 + 2 * i]; y[i] = raw[2 + 2 * i]; }
const O = Outgoer(E, { J: 1, seed: +process.argv[3], N: 0, S: 0 }); O.begin({ x, y });
const worlds = [{ x, y }]; for (let L = 1; L <= 3; L++) { O.next(); worlds.push({ x: O.world.x.slice(), y: O.world.y.slice() }); }
const P = 16, W = 2 * P;                                 // pixels out to 16 from the root: 128 real cells
const field = k => { const g = new Array(W * W).fill(0), w = worlds[k], sh = 3 - k; for (let i = 0; i < w.x.length; i++) { const u = (w.x[i] >> sh) + P, v = (w.y[i] >> sh) + P; if (u >= 0 && v >= 0 && u < W && v < W) g[v * W + u]++; } return g; };
console.log(JSON.stringify([3, 2, 1, 0].map(field)));    // [s]: the field counted s stages below the pixel (s = 3: real cells)
