"""Does a stage of refining leave the cells it should, row by row?

The zoom out (page/outward.js) refines each ring in stages, and a fault at one stage is carried into all the
later ones.  When the sampler is left out, every world is a real shape halved and thinned, so the picture of
level L at stage j has a true counterpart: the world L - j levels down.  This counts both by layer, in zones
named by layers of the level's own world (the ring is layers 8 to 15, the strip refined with it 16 and 17):
the two layers just outside the seam, the rest of the ring, its last two layers, and the strip.

  python3 analysis/stage_test.py ROOTS [SEED=1] [JOBS=4]       # writes analysis/stage_test.json
"""
import sys, os, json, subprocess
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e2e_test import read, ROOT, B
LEVELS = 3                                                                           # a real shape out to 256 cells gives three levels


def one(job):
    k, fn, c, seed = job; tmp = os.path.join(ROOT, "zoomtest", "e2e", f"stage_{os.getpid()}.bin")
    if subprocess.run([os.path.join(ROOT, "ball"), fn, tmp, str(c[0]), str(c[1]), "256", "60000"], stderr=subprocess.DEVNULL).returncode: return None
    o = subprocess.run(["node", os.path.join(ROOT, "analysis", "stage_run.js"), tmp, str(seed + k), str(LEVELS)], capture_output=True, text=True)
    return json.loads(o.stdout) if o.returncode == 0 else None


def main():
    R = int(sys.argv[1]); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1; J = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    rng = np.random.default_rng(seed); os.makedirs(os.path.join(ROOT, "zoomtest", "e2e"), exist_ok=True)
    snaps = sorted(f for f in os.listdir(os.path.join(ROOT, "prod")) if f.startswith("snap_") and "n1024000" in f and f.endswith(".bin")); jobs = []
    for s in snaps:
        fn = os.path.join(ROOT, "prod", s); T = read(fn)
        for _ in range(-(-R // len(snaps))): c = T[rng.integers(len(T))]; jobs.append((len(jobs), fn, (int(c[0]), int(c[1])), seed * 100000))
    with Pool(J) as p: res = [r for r in p.map(one, jobs, chunksize=1) if r]
    n = len(res); zones = {"first two layers of the ring": (8, 10), "middle of the ring": (10, 14), "last two layers of the ring": (14, 16), "strip": (16, 18)}; rows = []
    for i in range(len(res[0])):
        L, j = res[0][i]["level"], res[0][i]["stage"]; made = np.array([r[i]["made"] for r in res], float); real = np.array([r[i]["real"] for r in res], float); u = 1 << j
        row = {"level": L, "stage": j, "last": j == L}
        for name, (a, b) in zones.items():
            m, t = made[:, a * u:b * u].sum(1), real[:, a * u:b * u].sum(1); row[name] = [float(m.mean() / t.mean()), float((m - t).std(ddof=1) / np.sqrt(n) / t.mean())]
        rows.append(row)
    json.dump({"roots": n, "samples": len(snaps), "levels": LEVELS, "zones": list(zones), "rows": rows}, open(os.path.join(ROOT, "analysis", "stage_test.json"), "w"), indent=1)
    print(f"{n} roots; cells made / true cells at that resolution, by zone")
    for r in rows: print(f"  level {r['level']} stage {r['stage']}{' (last)' if r['last'] else '       '}: " + "; ".join(f"{k} {r[k][0]:.3f} ± {r[k][1]:.3f}" for k in zones))


if __name__ == "__main__": main()
