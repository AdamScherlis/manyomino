// For analysis/world_test.py: the level-1 world made from BALL.bin by the page's code, beside the halved real shape
// it starts from (same box, same held layers), by layer of the world.
//   node analysis/world_run.js BALL.bin SEED [key=value ...]  ->  one JSON line
const fs = require('fs'), path = require('path');
const { Outgoer, outRad, OUT } = require(path.join(__dirname, '..', 'page', 'outward.js'));
const E = new WebAssembly.Instance(new WebAssembly.Module(fs.readFileSync(path.join(__dirname, '..', 'page', 'outward.wasm'))), {}).exports; E.init(12);
const raw = new Int32Array(fs.readFileSync(process.argv[2]).buffer.slice(0)), n = raw[0], x = new Int32Array(n), y = new Int32Array(n), B = OUT.B, real = new Array(2 * B).fill(0); let fine = 0;
for (let i = 0; i < n; i++) { x[i] = raw[1 + 2 * i]; y[i] = raw[2 + 2 * i]; const r = outRad(x[i], y[i]); if (r < 2 * B) real[r]++; if (r < 4 * B) fine++; }
const key = (a, b) => (b + 4096) * 8192 + a + 4096;
function facts(c) { const lay = new Array(2 * B).fill(0), s = new Set(); let joins = 0;
  for (let i = 0; i < c.x.length; i++) { const r = outRad(c.x[i], c.y[i]); if (r < 2 * B) lay[r]++; s.add(key(c.x[i], c.y[i])); }
  for (let i = 0; i < c.x.length; i++) if (outRad(c.x[i], c.y[i]) === B) for (const [dx, dy] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) if (outRad(c.x[i] + dx, c.y[i] + dy) === B - 1 && s.has(key(c.x[i] + dx, c.y[i] + dy))) joins++;
  return { lay, joins }; }
let halved = null; const opt = { J: 0, seed: +process.argv[3], onHalved: c => { halved = facts(c); } }; process.argv.slice(4).forEach(a => { const [k, v] = a.split('='); opt[k] = +v; });
const O = Outgoer(E, opt); O.begin({ x, y }); O.next();
console.log(JSON.stringify({ real, fine, halved, world: facts(O.world), crossings: O.marks.length / 2 }));   // fine: real cells within twice the picture's reach, which the halved shape's 2B layers stand for
