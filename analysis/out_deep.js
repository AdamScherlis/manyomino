// The zoom out run a long way, to see that nothing drifts: the page's own code (page/outward.js and
// page/outward.wasm) at the resolution the page uses, from the real cells around one cell of a sample.
//   node analysis/out_deep.js BALL.bin LEVELS SEED [J=6]      # writes analysis/out_deep.json
// For each level: seconds to make it; cells of the world in its ring; parts kept apart in its box; cells added as
// joins; and, in the finished picture of the ring, cells, empty sites next to the shape per cell and cells with
// one neighbour per cell (two numbers that describe the fine texture).
const fs = require('fs'), path = require('path');
const { Outgoer, OUT, outRad } = require(path.join(__dirname, '..', 'page', 'outward.js'));
const [ball, levels, seed] = process.argv.slice(2, 5), J = +(process.argv[5] || 6);
const E = new WebAssembly.Instance(new WebAssembly.Module(fs.readFileSync(path.join(__dirname, '..', 'page', 'outward.wasm'))), {}).exports; E.init(12);
const raw = new Int32Array(fs.readFileSync(ball).buffer.slice(0)), n = raw[0], x = new Int32Array(n), y = new Int32Array(n);
for (let i = 0; i < n; i++) { x[i] = raw[1 + 2 * i]; y[i] = raw[2 + 2 * i]; }
const O = Outgoer(E, { J, seed: +seed }), G = O.G, B = OUT.B; O.begin({ x, y }); const rows = [];
for (let L = 1; L <= +levels; L++) {
  const pic = O.next(), w = O.world, h = G >> 1; let ringWorld = 0, cells = 0, empty = 0, tips = 0;
  for (let i = 0; i < w.x.length; i++) { const r = outRad(w.x[i], w.y[i]); if (r >= B && r < 2 * B) ringWorld++; }
  for (let v = 1; v < G - 1; v++) for (let u = 1; u < G - 1; u++) {
    const r = outRad(u - h, v - h); if (r < h / 2 + 8 || r >= h - 8) continue;
    const k = v * G + u, nb = (pic[k - 1] ? 1 : 0) + (pic[k + 1] ? 1 : 0) + (pic[k - G] ? 1 : 0) + (pic[k + G] ? 1 : 0);
    if (pic[k]) { cells++; if (nb === 1) tips++; } else if (nb) empty++;
  }
  rows.push({ level: L, s: pic.info.ms / 1000, sWorld: pic.info.msWorld / 1000, ringWorld, parts: pic.info.parts, links: pic.info.links, loose: pic.info.loose, cells, empty: empty / cells, tips: tips / cells });
  console.log(JSON.stringify(rows[rows.length - 1]));
}
fs.writeFileSync(path.join(__dirname, 'out_deep.json'), JSON.stringify({ J, G, levels: +levels, settings: OUT, rows }, null, 1));
