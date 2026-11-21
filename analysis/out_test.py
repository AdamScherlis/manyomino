"""The outward step, end to end, against the truth.

One step outward takes a polyomino of N cells that stands for the shape around a root at one
scale and returns one that stands for it at twice the scale:

  1. coarsen it (src/zoom.cpp -coarsen: every 2x2 block holding a cell becomes a cell, and the
     local walk thins the result for TC sweeps);
  2. hold the box of half-width B around the root, add cells outside it until there are N
     again, and run the real sampler with the box held for S moves per cell (manyomino -freeze):
     that samples the rest of a polyomino of N cells given what is in the box.

The truth is at hand for a few steps.  A real sample of 1,024,000 cells, coarsened whole, is the
real shape at each coarser scale (zoomtest/truth/TAG_lK.bin, K = 1..3, from
scripts/coarsen_chain.sh with KEEP=1).  The test starts from a polyomino of N cells cut from the
level-1 truth around a root, with its surroundings resampled, steps it outward twice, and
compares the cells in rings around the root with the truth at levels 2 and 3.

  python3 analysis/out_test.py N S ROOTS [B=32] [TC=256] [SEED=1]      # writes zoomtest/out/out_N_S_B_TC_SEED.json
"""
import sys, os, json, subprocess
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ext_test import start_state, RINGS, ROOT

BIN = os.path.join(ROOT, "manyomino"); ZOOM = os.path.join(ROOT, "zoom")


def read(fn): a = np.fromfile(fn, dtype="<i4"); return a[1:1 + 2 * a[0]].reshape(-1, 2).astype(np.int64)
def write(fn, xy): np.concatenate([[len(xy)], xy.astype("<i4").ravel()]).astype("<i4").tofile(fn)


def ring_counts(xy, c, B):
    d = xy - c; r = np.maximum(np.where(d[:, 0] < 0, -d[:, 0] - 1, d[:, 0]), np.where(d[:, 1] < 0, -d[:, 1] - 1, d[:, 1]))
    return [int(((r >= a * B) & (r < b * B)).sum()) for a, b in zip(RINGS[:-1], RINGS[1:])], int((r < B).sum())


def pad(xy, c, B, N, rng):
    """add cells at random empty sites next to the shape, outside the box, until there are N"""
    occ = {(int(a), int(b)) for a, b in xy.tolist()}; cells = list(occ); out = lambda x, y: not (c[0] - B <= x < c[0] + B and c[1] - B <= y < c[1] + B)
    while len(cells) < N:
        x, y = cells[rng.integers(len(cells))]; dx, dy = ((1, 0), (-1, 0), (0, 1), (0, -1))[rng.integers(4)]; s = (x + dx, y + dy)
        if s not in occ and out(*s): occ.add(s); cells.append(s)
    return np.array(cells, dtype=np.int64)


def resample(xy, c, B, S, tmp, seed):
    write(tmp + "_i.bin", xy); n = len(xy)
    subprocess.run([BIN, "-n", str(n), "-load", tmp + "_i.bin", "-freeze", str(c[0] - B), str(c[1] - B), str(c[0] + B), str(c[1] + B), "-pcp", "0.4", "-seed", str(seed),
                    "-moves", str(S * n), "-every", str(S * n), "-out", tmp + ".txt", "-snap", tmp + "_o.bin"], check=True, stderr=subprocess.DEVNULL)
    return read(tmp + "_o.bin")


def coarsen(xy, TC, tmp, seed):
    write(tmp + "_c.bin", xy)
    subprocess.run([ZOOM, "-load", tmp + "_c.bin", "-coarsen", "-sweeps", str(TC), "-seed", str(seed), "-out", tmp + "_d.bin"], check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    return read(tmp + "_d.bin")


def main():
    N, S, R = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); B = int(sys.argv[4]) if len(sys.argv) > 4 else 32
    TC = int(sys.argv[5]) if len(sys.argv) > 5 else 256; seed = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    rng = np.random.default_rng(seed); os.makedirs(os.path.join(ROOT, "zoomtest", "out"), exist_ok=True)
    tags = sorted(f[:-7] for f in os.listdir(os.path.join(ROOT, "zoomtest", "truth")) if f.endswith("_l3.bin"))
    tmp = os.path.join(ROOT, "zoomtest", "out", f"tmp_{N}_{S}_{B}_{TC}_{seed}"); rows = []; k = 0
    while len(rows) < R:
        tag = tags[(k * 5 + seed) % len(tags)]; k += 1
        T = [None] + [read(os.path.join(ROOT, "zoomtest", "truth", f"{tag}_l{j}.bin")) for j in (1, 2, 3)]
        for _ in range(2):
            c1 = T[1][rng.integers(len(T[1]))]; c = [None, c1, c1 >> 1, c1 >> 2]
            st = start_state(os.path.join(ROOT, 'zoomtest', 'truth', f'{tag}_l1.bin'), c1, B, N, tmp)
            if st is None: continue
            Q = resample(st + c1, c1, B, S, tmp, int(rng.integers(1 << 30))); row = {"tag": tag, "truth": [], "gen": [], "inside": []}
            for j in (1, 2, 3):
                if j > 1:
                    C = coarsen(Q, TC, tmp, int(rng.integers(1 << 30))); row.setdefault("kept", []).append(len(C))
                    if len(C) > N: break
                    Q = resample(pad(C, c[j], B, N, rng), c[j], B, S, tmp, int(rng.integers(1 << 30)))
                t, ti = ring_counts(T[j], c[j], B); g, gi = ring_counts(Q, c[j], B); row["truth"].append(t); row["gen"].append(g); row["inside"].append([ti, gi])
            if len(row["gen"]) < 3: continue
            rows.append(row); print(len(rows), tag, row["inside"], row["truth"][1], row["gen"][1], row["truth"][2], row["gen"][2], flush=True)
            if len(rows) >= R: break
    json.dump({"N": N, "S": S, "B": B, "TC": TC, "rings": RINGS, "rows": rows}, open(os.path.join(ROOT, "zoomtest", "out", f"out_{N}_{S}_{B}_{TC}_{seed}.json"), "w"))
    for j, name in ((0, "level 1 (surroundings resampled)"), (1, "level 2 (one step out)"), (2, "level 3 (two steps out)")):
        Tm = np.array([r["truth"][j] for r in rows], float); Gm = np.array([r["gen"][j] for r in rows], float)
        boot = np.array([Gm[i].sum(0) / Tm[i].sum(0) for i in (rng.integers(len(rows), size=len(rows)) for _ in range(400))])
        print(f"N={N} S={S} B={B} TC={TC} {name}: generated over truth by ring " + "  ".join(f"{a}-{b}: {g / t:.3f}±{e:.3f}" for a, b, g, t, e in zip(RINGS[:-1], RINGS[1:], Gm.sum(0), Tm.sum(0), boot.std(0))))


if __name__ == "__main__":
    main()
