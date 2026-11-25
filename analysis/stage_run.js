// For analysis/stage_test.py: the zoom out on one real shape, only halving and refining (no sampler), with the cells
// counted by layer at the end of every stage of every level's picture, and the true cells at that resolution
// beside them (the world as many levels down, which here is the real shape halved and thinned that many times).
//   node analysis/stage_run.js BALL.bin SEED LEVELS
const fs = require('fs'), path = require('path');
const { Outgoer, OUT, outRad } = require(path.join(__dirname, '..', 'page', 'outward.js'));
const [ball, seed, LV] = [process.argv[2], +process.argv[3], +process.argv[4]], B = OUT.B;
const wasm = new WebAssembly.Module(fs.readFileSync(process.env.OUTWARD_WASM || path.join(__dirname, '..', 'page', 'outward.wasm')));
const raw = new Int32Array(fs.readFileSync(ball).buffer.slice(0)), n = raw[0], x = new Int32Array(n), y = new Int32Array(n);
for (let i = 0; i < n; i++) { x[i] = raw[1 + 2 * i]; y[i] = raw[2 + 2 * i]; }
const layers = (c, m) => { const r = new Array(m).fill(0); for (let i = 0; i < c.x.length; i++) { const q = outRad(c.x[i], c.y[i]); if (q < m) r[q]++; } return r; };
const engine = () => { const E = new WebAssembly.Instance(wasm, {}).exports; E.init(12); return E; };
// the worlds first (a run that makes no picture worth keeping), then the pictures
const worlds = [{ x, y }]; { const O = Outgoer(engine(), { J: 1, seed, N: 0, S: 0 }); O.begin({ x, y }); for (let L = 1; L <= LV; L++) { O.next(); worlds.push({ x: O.world.x.slice(), y: O.world.y.slice() }); } }
const out = [];
const O = Outgoer(engine(), { J: LV, seed, N: 0, S: 0, onStage: (L, j, what, c) => { if (what !== 'end') return; const m = (2 * B + OUT.MARGIN) << j; out.push({ level: L, stage: j, made: layers(c, m), real: layers(worlds[L - j], m) }); } });
O.begin({ x, y }); for (let L = 1; L <= LV; L++) O.next();
console.log(JSON.stringify(out));
