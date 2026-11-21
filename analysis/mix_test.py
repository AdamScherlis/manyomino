"""How long must the sampler run, with the box held, before the padding is forgotten?

A piece of M cells is cut from a real coarsened sample around a root (src/ball.cpp), padded to
N cells in one of two ways (cells at random sites next to the shape; or turned copies of the
piece hung on its tips), and the sampler is run with the box of half-width B held.  The cells in
rings round the root are counted as it goes and compared with the real shape's.

  python3 analysis/mix_test.py N M ROOTS [B=8] [LEVEL=2] [SEED=1] [JOBS=4]   # writes zoomtest/mix/mix_N_M_B_LEVEL_SEED.json
  S_AT="25 100 400" MODES="1" python3 analysis/mix_test.py ...               # other run lengths; one kind of padding only
"""
import sys, os, json, subprocess
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ext_test import RINGS, ROOT

S_AT = [int(v) for v in os.environ.get("S_AT", "25 50 100 200 400 800 1600").split()]
MODES = [int(v) for v in os.environ.get("MODES", "0 1").split()]


def read(fn): a = np.fromfile(fn, dtype="<i4"); return a[1:1 + 2 * a[0]].reshape(-1, 2).astype(np.int64)


def rings(xy, c, B):
    d = xy - c; r = np.maximum(np.where(d[:, 0] < 0, -d[:, 0] - 1, d[:, 0]), np.where(d[:, 1] < 0, -d[:, 1] - 1, d[:, 1]))
    return [int((r < B).sum())] + [int(((r >= a * B) & (r < b * B)).sum()) for a, b in zip(RINGS[:-1], RINGS[1:])]


def one(job):
    k, fn, c, N, M, B, seed = job; tmp = os.path.join(ROOT, "zoomtest", "mix", f"tmp_{os.getpid()}.bin")
    if subprocess.run([os.path.join(ROOT, "ball"), fn, tmp, str(c[0]), str(c[1]), str(B), str(M)], stderr=subprocess.DEVNULL).returncode: return None
    out = {}
    for mode in MODES:
        o = subprocess.run([os.path.join(ROOT, "outward_native"), "mix", tmp, str(N), str(B), str(seed + k), str(mode)] + [str(s) for s in S_AT], capture_output=True, text=True, check=True).stdout
        out[mode] = [[float(v) for v in l.split()] for l in o.strip().split("\n")]
    return out


def main():
    N, M, R = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); B = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    L = int(sys.argv[5]) if len(sys.argv) > 5 else 2; seed = int(sys.argv[6]) if len(sys.argv) > 6 else 1; J = int(sys.argv[7]) if len(sys.argv) > 7 else 4
    rng = np.random.default_rng(seed); os.makedirs(os.path.join(ROOT, "zoomtest", "mix"), exist_ok=True)
    tags = sorted(f[:-7] for f in os.listdir(os.path.join(ROOT, "zoomtest", "truth")) if f.endswith(f"_l{L}.bin"))
    jobs = []; truth = []
    per = -(-R // len(tags))
    for t in tags:
        fn = os.path.join(ROOT, "zoomtest", "truth", f"{t}_l{L}.bin"); T = read(fn)
        for _ in range(per):
            c = T[rng.integers(len(T))]; jobs.append((len(jobs), fn, (int(c[0]), int(c[1])), N, M, B, seed * 100000)); truth.append(rings(T, c, B))
    with Pool(J) as p: res = p.map(one, jobs, chunksize=1)
    keep = [i for i, r in enumerate(res) if r is not None]; T = np.array([truth[i] for i in keep], float)
    json.dump({"N": N, "M": M, "B": B, "level": L, "S": [0] + S_AT, "rings": RINGS, "truth": T.tolist(), "runs": [res[i] for i in keep]}, open(os.path.join(ROOT, "zoomtest", "mix", f"mix_{N}_{M}_{B}_{L}_{seed}.json"), "w"))
    print(f"N={N} M={M} B={B} level {L}: {len(keep)} roots; generated over real, by ring (inside the box first), and Rg^2")
    for mode, name in ((0, "random sites"), (1, "copies on the tips")):
        if mode not in MODES: continue
        print(" padding:", name)
        for si, s in enumerate([0] + S_AT):
            G = np.array([res[i][mode][si][2:10] for i in keep]); r2 = np.mean([res[i][mode][si][10] for i in keep])
            boot = np.array([G[j].sum(0) / T[j].sum(0) for j in (rng.integers(len(keep), size=len(keep)) for _ in range(300))])
            print(f"  S={s:5d}: " + " ".join(f"{g / t:.2f}±{e:.2f}" for g, t, e in zip(G.sum(0), T.sum(0), boot.std(0))) + f"   {r2:.0f}")


if __name__ == "__main__":
    main()
