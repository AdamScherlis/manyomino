#!/bin/sh
# Four long chains of worlds (no pictures) from the real cells around one cell, to see whether the amount of shape near
# the root drifts from level to level.  Writes analysis/out_drift.json.
#   scripts/out_drift.sh BALL.bin [LEVELS=120]        (BALL.bin from ./ball, as for analysis/out_run.js)
cd "$(dirname "$0")/.."
B=$1; L=${2:-120}; mkdir -p zoomtest/e2e
for k in 1 2 3 4; do node analysis/out_drift.js "$B" "$L" $((100 + k)) > zoomtest/e2e/drift_$k.txt & done; wait
python3 - <<'PY'
import json, glob, numpy as np
ch = [[json.loads(l) for l in open(f) if l.startswith("{")] for f in sorted(glob.glob("zoomtest/e2e/drift_*.txt"))]
n = min(len(c) for c in ch); q = n // 4; out = {"chains": len(ch), "levels": n, "quarters": {}}
for k in ("box", "ring", "far", "parts"):
    a = np.array([[r[k] for r in c[:n]] for c in ch], float)
    out["quarters"][k] = [[float(a[:, i * q:(i + 1) * q].mean()), float(a[:, i * q:(i + 1) * q].mean(1).std(ddof=1) / np.sqrt(len(ch)))] for i in range(4)]
json.dump(out, open("analysis/out_drift.json", "w"), indent=1); print(json.dumps(out["quarters"]))
PY
