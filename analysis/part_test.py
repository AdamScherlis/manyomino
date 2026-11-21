"""Do the worlds of the zoom out rejoin their strands where real shapes do?

In a real sample, the cells within 16 of a given cell fall into a few separate pieces: strands that
are joined to each other only further out.  Widen the window and they join up.  The zoom out keeps
such pieces apart inside the box of each world (see page/outward.js) and leaves the sampler to
join them outside it.  This counts, for real samples and for the first world made from them, how
many of the pieces are still separate within each wider window.

  python3 analysis/part_test.py ROOTS [SEED=1] [JOBS=4]       # writes analysis/part_test.json
"""
import sys, os, json, subprocess
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e2e_test import read, raster, pieces, rad, ROOT, B
RS = [2 * B, 2 * B + 2, 2 * B + 4, 3 * B, 4 * B, 6 * B, 8 * B]          # in real cells; half of each in world cells


def one(job):
    k, fn, c, seed = job; T = read(fn); R = RS[-1]; g = raster(T - np.array(c), R); r = rad(R)
    lab0, s0 = pieces(g * (r < RS[0]), True); real = []
    for Rk in RS:
        labk, _ = pieces(g * (r < Rk), True); real.append(len({int(labk[lab0 == i + 1][0]) for i in range(len(s0))}))
    tmp = os.path.join(ROOT, "zoomtest", "e2e", f"part_{os.getpid()}.bin")
    if subprocess.run([os.path.join(ROOT, "ball"), fn, tmp, str(c[0]), str(c[1]), "32", "16000"], stderr=subprocess.DEVNULL).returncode: return None
    o = subprocess.run(["node", os.path.join(ROOT, "analysis", "part_run.js"), tmp, str(seed + k)], capture_output=True, text=True)
    return (real, json.loads(o.stdout)) if o.returncode == 0 else None


def main():
    R = int(sys.argv[1]); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1; J = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    rng = np.random.default_rng(seed); os.makedirs(os.path.join(ROOT, "zoomtest", "e2e"), exist_ok=True)
    snaps = sorted(f for f in os.listdir(os.path.join(ROOT, "prod")) if f.startswith("snap_") and "n1024000" in f and f.endswith(".bin")); jobs = []
    for s in snaps:
        fn = os.path.join(ROOT, "prod", s); T = read(fn)
        for _ in range(-(-R // len(snaps))): c = T[rng.integers(len(T))]; jobs.append((len(jobs), fn, (int(c[0]), int(c[1])), seed * 100000))
    with Pool(J) as p: res = [r for r in p.map(one, jobs, chunksize=1) if r]
    a = np.array([r[0] for r in res], float); b = np.array([r[1] for r in res], float); n = len(res)
    out = {"roots": n, "radii": RS, "real": [[float(m), float(s)] for m, s in zip(a.mean(0), a.std(0, ddof=1) / np.sqrt(n))], "made": [[float(m), float(s)] for m, s in zip(b.mean(0), b.std(0, ddof=1) / np.sqrt(n))]}
    json.dump(out, open(os.path.join(ROOT, "analysis", "part_test.json"), "w"), indent=1)
    print(f"{n} roots; pieces of what lies within {RS[0]} cells that are still separate within R cells")
    print("  R        " + "".join(f"{r:>14d}" for r in RS)); print("  real     " + "".join(f"{m:8.2f}±{s:.2f} " for m, s in out["real"])); print("  made     " + "".join(f"{m:8.2f}±{s:.2f} " for m, s in out["made"]))


if __name__ == "__main__":
    main()
