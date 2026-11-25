// Runs the page's zoom out (page/outward.js + page/outward.wasm) outside the browser, for tests.
//   node analysis/out_run.js BALL.bin OUT_PREFIX LEVELS SEED [J=5] [key=value ...]
// (OUTWARD_WASM names another build of the engine: outward_big.wasm takes worlds of up to 262,144 cells, which the
// runs that only halve and refine a large real shape need.)
// BALL.bin: the world at level 0 (int32 count, then x y pairs, root at the origin; src/ball.cpp makes one).
// Writes OUT_PREFIX_l<L>.bin for each level: the generated picture at the resolution of real cells,
// r < 2 B 2^L (the same format), and prints one line of facts per level.
const fs = require('fs'), path = require('path');
const { Outgoer, OUT } = require(path.join(__dirname, '..', 'page', 'outward.js'));
const [ball, prefix, levels, seed] = process.argv.slice(2, 6), J = +(process.argv[6] || 5), opt = { J, seed: +seed };
process.argv.slice(7).forEach(a => { const [k, v] = a.split('='); opt[k] = +v; });
const E = new WebAssembly.Instance(new WebAssembly.Module(fs.readFileSync(process.env.OUTWARD_WASM || path.join(__dirname, '..', 'page', 'outward.wasm'))), {}).exports;
E.init(12);
const raw = new Int32Array(fs.readFileSync(ball).buffer.slice(0)), n = raw[0], x = new Int32Array(n), y = new Int32Array(n);
for (let i = 0; i < n; i++) { x[i] = raw[1 + 2 * i]; y[i] = raw[2 + 2 * i]; }
const O = Outgoer(E, opt); O.begin({ x, y });
for (let L = 1; L <= +levels; L++) {
  const pic = O.next(), c = O.stages[Math.min(L, J)], out = new Int32Array(1 + 2 * c.x.length);
  out[0] = c.x.length; for (let i = 0; i < c.x.length; i++) { out[1 + 2 * i] = c.x[i]; out[2 + 2 * i] = c.y[i]; }
  if (prefix !== '-') fs.writeFileSync(prefix + '_l' + L + '.bin', Buffer.from(out.buffer));
  console.log(JSON.stringify(Object.assign({ cells: c.x.length }, pic.info)));
}
