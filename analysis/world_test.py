"""Is a world of the zoom out the world it should be?

Two things are asked of the step from one world to the next (page/outward.js).  Halving must leave a shape
with as many cells, layer by layer round the root, as a real polyomino has round one of its cells: that is
what "the same thing at twice the scale" means, and it fixes the weight of the thinning walk.  And the
sampler, which draws afresh everything beyond the box and the layers held with it, must give what the
halved real shape has there.  For cells of real samples this counts, by layer: the real shape at its own
scale, the real shape halved, and the first world made from it.

  python3 analysis/world_test.py ROOTS [SEED=1] [JOBS=4]       # writes analysis/world_test.json
"""
import sys, os, json, subprocess
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e2e_test import read, ROOT, B


def one(job):
    k, fn, c, seed = job; tmp = os.path.join(ROOT, "zoomtest", "e2e", f"world_{os.getpid()}.bin")
    if subprocess.run([os.path.join(ROOT, "ball"), fn, tmp, str(c[0]), str(c[1]), "32", "16000"], stderr=subprocess.DEVNULL).returncode: return None
    o = subprocess.run(["node", os.path.join(ROOT, "analysis", "world_run.js"), tmp, str(seed + k)], capture_output=True, text=True)
    return json.loads(o.stdout) if o.returncode == 0 else None


def main():
    R = int(sys.argv[1]); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1; J = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    rng = np.random.default_rng(seed); os.makedirs(os.path.join(ROOT, "zoomtest", "e2e"), exist_ok=True)
    snaps = sorted(f for f in os.listdir(os.path.join(ROOT, "prod")) if f.startswith("snap_") and "n1024000" in f and f.endswith(".bin")); jobs = []
    for s in snaps:
        fn = os.path.join(ROOT, "prod", s); T = read(fn)
        for _ in range(-(-R // len(snaps))): c = T[rng.integers(len(T))]; jobs.append((len(jobs), fn, (int(c[0]), int(c[1])), seed * 100000))
    with Pool(J) as p: res = [r for r in p.map(one, jobs, chunksize=1) if r]
    n = len(res); real = np.array([r["real"] for r in res], float); T = np.array([r["halved"]["lay"] for r in res], float); W = np.array([r["world"]["lay"] for r in res], float)
    se = lambda a: float(a.std(ddof=1) / np.sqrt(n)); ms = lambda a: [float(a.mean()), se(a)]
    bands = {"box": (0, B), "ring": (B, 2 * B), "first two layers drawn afresh": None}
    hold = int(np.argmax(np.abs(W - T).sum(0) > 0)) if np.abs(W - T).sum() > 0 else 2 * B     # the first layer the sampler changed
    bands["first two layers drawn afresh"] = (hold, hold + 2)
    out = {"roots": n, "samples": len(snaps), "held to": hold, "layers": {"real": [float(v) for v in real.mean(0)], "halved": [float(v) for v in T.mean(0)], "world": [float(v) for v in W.mean(0)]},
           "bands": {k: {"real": ms(real[:, a:b].sum(1)), "halved": ms(T[:, a:b].sum(1)), "world": ms(W[:, a:b].sum(1)), "world - halved": ms((W - T)[:, a:b].sum(1))} for k, (a, b) in bands.items()},
           "joins across the edge of the box": {"halved": ms(np.array([r["halved"]["joins"] for r in res], float)), "world": ms(np.array([r["world"]["joins"] for r in res], float))},
           "crossings": ms(np.array([r["crossings"] for r in res], float)),
           "real cells per halved cell": float(sum(r["fine"] for r in res) / T.sum())}
    json.dump(out, open(os.path.join(ROOT, "analysis", "world_test.json"), "w"), indent=1)
    print(f"{n} roots; the sampler holds layers below {hold}; {out['real cells per halved cell']:.3f} real cells for each cell of the halved shape")
    print("layer    " + " ".join(f"{k:6d}" for k in range(B - 3, 2 * B))); 
    for nm, a in (("real", real), ("halved", T), ("world", W)): print(f"{nm:8s} " + " ".join(f"{v:6.2f}" for v in a.mean(0)[B - 3:]))
    for k, v in out["bands"].items(): print(f"  {k}: real {v['real'][0]:.2f} ± {v['real'][1]:.2f}, halved {v['halved'][0]:.2f} ± {v['halved'][1]:.2f}, world {v['world'][0]:.2f} ± {v['world'][1]:.2f}, world - halved {v['world - halved'][0]:+.2f} ± {v['world - halved'][1]:.2f}")


if __name__ == "__main__":
    main()
