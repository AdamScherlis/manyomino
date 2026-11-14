// Dry run of the video's chain: prints what each bar does, writes events.json for the synth.
const { Sampler } = require('../page/sampler.js'), SIM = require('./sim.js'), fs = require('fs');
const sim = SIM.makeSim(Sampler); const t0 = Date.now();
let acc = 0, maxm = 0, rel = 0, top = [];
const EV = [];
for (let f = 0; f < SIM.NF; f++) {
  const r = sim.step(f); rel += r.relAcc;
  EV.push({ f, name: r.name, R: r.R, T: r.T, rel: r.relAcc, nacc: r.burst ? r.burst.length : 0, rg2: +sim.s.rg2().toFixed(1), top: r.burst ? r.burst.slice(0, 4).map(e => ({ m: e.m, pan: +e.pan.toFixed(3) })) : [] });
  if (r.burst) { acc += r.burst.length; if (r.burst.length) { maxm = Math.max(maxm, r.burst[0].m); top.push(r.burst[0].m); } }
  if (f % SIM.BAR === SIM.BAR - 1) {
    console.log(`bar ${String(f / SIM.BAR | 0).padStart(2)} ${r.name.padEnd(6)} moves ${r.moves} rg2 ${sim.s.rg2().toFixed(0).padStart(5)} relAcc ${String(rel).padStart(6)} cpAcc ${String(acc).padStart(5)} max ${String(maxm).padStart(4)}  tops ${top.join(' ')}`);
    acc = 0; maxm = 0; rel = 0; top = [];
  }
}
console.log((Date.now() - t0) / 1000, 's');
if (process.argv[2]) fs.writeFileSync(process.argv[2], JSON.stringify(EV));
