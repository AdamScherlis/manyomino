"""Far out, does a picture made of coarse cells look like the real shape seen from that far?

Past the level at which a texel of the zoom out is one real cell, a texel stands for a block of them, and the
screen averages a few texels into each pixel.  What a real shape would show in a pixel is the number of its
cells there.  For real shapes this compares the two at one pixel size (8 real cells): the count of real cells
in each pixel, and the count of cells of the shape halved and thinned (as the worlds of the zoom out are) until
its cells are 1, 2 and 4 to a pixel side, which is 0, 1 and 2 stages finer than the pixel.  The page draws its
far levels between 1 and 2 stages finer than a pixel of the screen.

  python3 analysis/mass_test.py ROOTS [SEED=1] [JOBS=4]       # writes analysis/mass_test.json
"""
import sys, os, json, subprocess
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e2e_test import read, ROOT


def one(job):
    k, fn, c, seed = job; tmp = os.path.join(ROOT, "zoomtest", "e2e", f"mass_{os.getpid()}.bin")
    if subprocess.run([os.path.join(ROOT, "ball"), fn, tmp, str(c[0]), str(c[1]), "256", "60000"], stderr=subprocess.DEVNULL).returncode: return None
    o = subprocess.run(["node", os.path.join(ROOT, "analysis", "mass_run.js"), tmp, str(seed + k)], capture_output=True, text=True)
    return json.loads(o.stdout) if o.returncode == 0 else None


def main():
    R = int(sys.argv[1]); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1; J = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    rng = np.random.default_rng(seed); os.makedirs(os.path.join(ROOT, "zoomtest", "e2e"), exist_ok=True)
    snaps = sorted(f for f in os.listdir(os.path.join(ROOT, "prod")) if f.startswith("snap_") and "n1024000" in f and f.endswith(".bin")); jobs = []
    for s in snaps:
        fn = os.path.join(ROOT, "prod", s); T = read(fn)
        for _ in range(-(-R // len(snaps))): c = T[rng.integers(len(T))]; jobs.append((len(jobs), fn, (int(c[0]), int(c[1])), seed * 100000))
    with Pool(J) as p: res = [r for r in p.map(one, jobs, chunksize=1) if r]
    F = np.array(res, float); n = len(F); out = {"roots": n, "samples": len(snaps), "pixel": 8, "stages": []}      # F[root, s, pixel]
    real = F[:, 3]; rb = real / real[real > 0].mean()
    for s in range(4):
        f = F[:, s]; lit = (f > 0).sum(1); b = f / f[f > 0].mean(); both = (f > 0) & (real > 0)
        row = {"stages finer than the pixel": s, "real cells": s == 3, "lit pixels": [float(lit.mean()), float(lit.std(ddof=1) / np.sqrt(n))], "lit, of real": float(lit.mean() / (real > 0).sum(1).mean()),
               "spread of brightness": float(b[f > 0].std() / b[f > 0].mean()), "lit here only": float(((f > 0) & (real == 0)).sum(1).mean()), "lit in the real only": float(((f == 0) & (real > 0)).sum(1).mean()),
               "correlation with real": 1.0 if s == 3 else float(np.corrcoef(b[both], rb[both])[0, 1]) if b[both].std() > 0 else None}
        out["stages"].append(row)
    json.dump(out, open(os.path.join(ROOT, "analysis", "mass_test.json"), "w"), indent=1)
    print(f"{n} roots; pixels of 8 real cells out to 128 cells from the root")
    for r in out["stages"]: print(f"  {r['stages finer than the pixel']} stages finer{' (real cells)' if r['real cells'] else ''}: lit pixels {r['lit pixels'][0]:.1f} ± {r['lit pixels'][1]:.1f} ({100 * r['lit, of real']:.0f}% of real), "
                                  f"spread of brightness among lit pixels {r['spread of brightness']:.3f}, brightness against real over pixels lit in both r = {r['correlation with real'] if r['correlation with real'] is not None else float('nan'):.3f}, "
                                  f"lit here only {r['lit here only']:.1f}, in the real only {r['lit in the real only']:.1f}")


if __name__ == "__main__": main()
