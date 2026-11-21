"""How well does a finite polyomino stand in for the surroundings of a window?

The exact way to extend a known window outwards is to sample the rest of the polyomino given
it: the real sampler run with the window's cells held (manyomino -freeze).  That is exact for a
polyomino of N cells; the question is how large N must be, and how long the chain must run,
before the surroundings of the window are indistinguishable from those in a very large one.

Here the window is the box of half-width B around a cell of a real sample of 1,024,000 cells.
The starting state is a connected piece of N cells of the real sample around the box
(src/ball.cpp); the sampler then runs for S moves per cell with the box held; and
the cells it leaves in rings around the box are counted against what the sample really has.

  python3 analysis/ext_test.py N S ROOTS [B=32] [SEED=1]      # appends to zoomtest/ext/ext_N_S.json
"""
import sys, os, json, subprocess, collections
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
RINGS = [1, 1.5, 2, 2.5, 3, 4, 6, 8]           # ring edges, in units of the box half-width


def start_state(sample, root, B, N, tmp):
    """a connected piece of N cells of the sample around the box (src/ball.cpp); None if there is none"""
    r = subprocess.run([os.path.join(ROOT, "ball"), sample, tmp + "_b.bin", str(int(root[0])), str(int(root[1])), str(B), str(N)], stderr=subprocess.DEVNULL)
    if r.returncode: return None
    return np.fromfile(tmp + "_b.bin", dtype="<i4")[1:].reshape(-1, 2).astype(np.int64)


def rings(d, B):
    r = np.maximum(np.where(d[:, 0] < 0, -d[:, 0] - 1, d[:, 0]), np.where(d[:, 1] < 0, -d[:, 1] - 1, d[:, 1]))
    return [int(((r >= a * B) & (r < b * B)).sum()) for a, b in zip(RINGS[:-1], RINGS[1:])], int((r < B).sum())


def main():
    N, S, R = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); B = int(sys.argv[4]) if len(sys.argv) > 4 else 32; seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    rng = np.random.default_rng(seed); files = sorted(f for f in os.listdir(os.path.join(ROOT, "prod")) if f.startswith("snap_") and "_n1024000_" in f)
    out = os.path.join(ROOT, "zoomtest", "ext", f"ext_{N}_{S}_{B}_{seed}.json"); rows = []; tmp = os.path.join(ROOT, "zoomtest", "ext", f"tmp_{N}_{S}_{B}_{seed}")
    k = 0; per = max(1, R // len(files) + 1)
    while len(rows) < R:
        f = files[(k * 7 + seed) % len(files)]; k += 1; fn = os.path.join(ROOT, "prod", f)
        a = np.fromfile(fn, dtype="<i4"); xy = a[1:].reshape(-1, 2).astype(np.int64)
        for _ in range(per):
            root = xy[rng.integers(len(xy))]; st = start_state(fn, root, B, N, tmp)
            if st is None: continue
            real, inside = rings(xy - root, B); init, _ = rings(st, B)
            if S > 0:
                subprocess.run([os.path.join(ROOT, "manyomino"), "-n", str(N), "-load", tmp + "_b.bin", "-freeze", str(-B), str(-B), str(B), str(B), "-pcp", "0.4", "-seed", str(int(rng.integers(1 << 30))),
                                "-moves", str(S * N), "-every", str(S * N), "-out", tmp + ".txt", "-snap", tmp + "_o.bin"], check=True, stderr=subprocess.DEVNULL)
                o = np.fromfile(tmp + "_o.bin", dtype="<i4")[1:].reshape(-1, 2).astype(np.int64); gen, ins2 = rings(o, B); assert ins2 == inside
            else: gen = init
            rows.append({"file": f, "inside": inside, "real": real, "start": init, "gen": gen})
            if len(rows) % 25 == 0: print(len(rows), f, inside, real, gen, flush=True)
            if len(rows) >= R: break
    json.dump({"N": N, "S": S, "B": B, "rings": RINGS, "rows": rows}, open(out, "w"))
    R_ = np.array([r["real"] for r in rows], float); G_ = np.array([r["gen"] for r in rows], float)
    boot = np.array([(G_[i].sum(0) / R_[i].sum(0)) for i in (rng.integers(len(rows), size=len(rows)) for _ in range(400))])
    print(f"N={N} S={S} B={B} roots={len(rows)}: generated over real by ring " + "  ".join(f"{a}-{b}: {g / r:.3f}±{e:.3f}" for a, b, g, r, e in zip(RINGS[:-1], RINGS[1:], G_.sum(0), R_.sum(0), boot.std(0))))


if __name__ == "__main__":
    main()
