#!/bin/sh
# Eight long chains of worlds (no pictures), each from the real cells around a cell of a different sample, to see
# whether the amount of shape near the root drifts from level to level.  Writes analysis/out_drift.json.
#   scripts/out_drift.sh [LEVELS=240] [key=value ...]     (needs ./ball and the samples in prod/; four chains at a time;
#   chains already finished in zoomtest/drift are kept, so a run that was cut short can be started again)
cd "$(dirname "$0")/.."
L=${1:-240}; [ $# -gt 0 ] && shift; mkdir -p zoomtest/drift
python3 - <<'PY'
import numpy as np, os, subprocess
snaps = sorted(f for f in os.listdir("prod") if f.startswith("snap_") and "n1024000" in f and f.endswith(".bin"))[:8]
for k, s in enumerate(snaps):
    if os.path.exists(f"zoomtest/drift/ball_{k}.bin"): continue          # (a run that was cut short is taken up again)
    a = np.fromfile("prod/" + s, dtype="<i4"); T = a[1:1 + 2 * a[0]].reshape(-1, 2); rng = np.random.default_rng(900 + k)
    while True:
        c = T[rng.integers(len(T))]
        if subprocess.run(["./ball", "prod/" + s, f"zoomtest/drift/ball_{k}.bin", str(c[0]), str(c[1]), "32", "16000"], stderr=subprocess.DEVNULL).returncode == 0: break
PY
for half in "0 1 2 3" "4 5 6 7"; do
  for k in $half; do
    [ "$(grep -c . zoomtest/drift/chain_$k.txt 2>/dev/null)" = "$L" ] || node analysis/out_drift.js zoomtest/drift/ball_$k.bin "$L" $((500 + k)) "$@" > zoomtest/drift/chain_$k.txt &
  done; wait
done
python3 - <<'PY'
import json, glob, numpy as np
ch = [[json.loads(l) for l in open(f) if l.startswith("{")] for f in sorted(glob.glob("zoomtest/drift/chain_*.txt"))]
n = min(len(c) for c in ch); x = np.arange(n); out = {"chains": len(ch), "levels": n, "what": {}}
for k in ("box", "ring", "far", "parts"):
    a = np.array([[r[k] for r in c[:n]] for c in ch], float); sl = np.array([np.polyfit(x, r, 1)[0] * 100 for r in a]); q = n // 4
    out["what"][k] = {"mean": float(a.mean()), "se": float(a.mean(1).std(ddof=1) / np.sqrt(len(ch))), "slope per 100 levels": float(sl.mean()), "slope se": float(sl.std(ddof=1) / np.sqrt(len(ch))),
                      "quarters": [float(a[:, i * q:(i + 1) * q].mean()) for i in range(4)]}
lay = np.array([[r["lay"] for r in c[:n]] for c in ch], float); out["layers"] = [float(v) for v in lay.mean((0, 1))]
json.dump(out, open("analysis/out_drift.json", "w"), indent=1); print(json.dumps(out["what"], indent=1))
PY
