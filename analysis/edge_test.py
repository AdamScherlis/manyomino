"""What should the halving of a world be?

Halving (every 2x2 block holding a cell becomes a cell) gives a shape that is too thick; a local
walk then thins it.  Run freely for long, the walk gives a typical polyomino (analysis/coarsen_test.py)
but also moves the tips of branches about, so the result no longer lines up with the finer shape.
Kept to the sites where the finer shape has cells ("on"), or to those and the sites round them
("near"), it stays lined up.  This compares what each leaves in the box round the root with the
long free walk, on the same real cells:

  cells per cell       finer cells for each cell of the halved world
  inside, edge, runs   cells in the box; cells on its outermost layer; runs of those
  parts                pieces of what is in the box
  bare, far            edge cells where the finer shape has no cell; none even next door
  lost                 edge cells of the plain halving with no cell left within one site
  touching             pairs of cells in the box that touch without being joined (different parts of the finer shape)

  python3 analysis/edge_test.py ROOTS [SEED=1] [JOBS=4] [VARIANTS]
      # VARIANTS: JSON, name -> settings for page/outward.js (RULE 0 free, 1 on, 3 near; TC sweeps; COVER 0/1; W weight;
      # PARTS 0/1); must include "free 256": {"RULE":0,"COVER":0,"W":1,"TC":256,"PARTS":0}
"""
import sys, os, json, subprocess
import numpy as np
from multiprocessing import Pool
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
NAMES = ["cells per cell", "inside", "edge", "runs", "parts", "next layer out", "bare", "far", "lost", "touching"]


def one(job):
    k, fn, c, seed, var = job; tmp = os.path.join(ROOT, "zoomtest", "e2e", f"edge_{os.getpid()}.bin")
    if subprocess.run([os.path.join(ROOT, "ball"), fn, tmp, str(c[0]), str(c[1]), "32", "16000"], stderr=subprocess.DEVNULL).returncode: return None
    o = subprocess.run(["node", os.path.join(ROOT, "analysis", "edge_run.js"), tmp, str(seed + k)] + ([var] if var else []), capture_output=True, text=True)
    return json.loads(o.stdout) if o.returncode == 0 else None


def main():
    R = int(sys.argv[1]); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1; J = int(sys.argv[3]) if len(sys.argv) > 3 else 4; var = sys.argv[4] if len(sys.argv) > 4 else ""
    rng = np.random.default_rng(seed); os.makedirs(os.path.join(ROOT, "zoomtest", "e2e"), exist_ok=True)
    snaps = sorted(f for f in os.listdir(os.path.join(ROOT, "prod")) if f.startswith("snap_") and "n1024000" in f and f.endswith(".bin")); jobs = []
    for s in snaps:
        fn = os.path.join(ROOT, "prod", s); a = np.fromfile(fn, dtype="<i4"); T = a[1:1 + 2 * a[0]].reshape(-1, 2)
        for _ in range(-(-R // len(snaps))): c = T[rng.integers(len(T))]; jobs.append((len(jobs), fn, (int(c[0]), int(c[1])), seed * 100000, var))
    with Pool(J) as p: res = [r for r in p.map(one, jobs, chunksize=1) if r]
    ref = np.array([r["free 256"] for r in res]); print(f"{len(res)} roots; mean (difference from the long free walk in standard errors)")
    print(f"  {'':12s}" + "".join(f"{nm:>22s}" for nm in NAMES))
    for v in res[0]:
        a = np.array([r[v] for r in res]); d = a - ref; se = d.std(0, ddof=1) / np.sqrt(len(res)) + 1e-12
        print(f"  {v:12s}" + "".join(f"{m:14.3f} ({z:+5.1f})" for m, z in zip(a.mean(0), d.mean(0) / se)))
    json.dump(res, open(os.path.join(ROOT, "zoomtest", "e2e", "edge.json"), "w"))


if __name__ == "__main__":
    main()
