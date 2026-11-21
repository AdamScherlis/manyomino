"""The zoom out, end to end, against real samples.

For a root cell of a real sample of 1,024,000 cells, the page's own code (page/outward.js and
page/outward.wasm, run by analysis/out_run.js) is given the real cells near the root and makes
LEVELS levels outward.  At the resolution of real cells the result is a picture of everything
within 8 * 2^(LEVELS+1) cells of the root, of which only the middle (within 16 cells) is real.
It is compared with the real cells in the same square around the same root.

  python3 analysis/e2e_test.py ROOTS [LEVELS=4] [SEED=1] [JOBS=4] [TAG=base] [key=value ...]
      # key=value pairs go to the engine (N, S, TC, T, ...); writes zoomtest/e2e/TAG.json and prints a table.
      # With TAG "final" it also writes analysis/e2e_test.json and analysis/e2e_pics.npz, which the page is built from.
"""
import sys, os, json, subprocess
import numpy as np
from multiprocessing import Pool
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
B = 8


def read(fn): a = np.fromfile(fn, dtype="<i4"); return a[1:1 + 2 * a[0]].reshape(-1, 2).astype(np.int64)


def raster(xy, R):
    """the square r < R about the root as a 2R x 2R array of 0 and 1"""
    k = (xy[:, 0] >= -R) & (xy[:, 0] < R) & (xy[:, 1] >= -R) & (xy[:, 1] < R); g = np.zeros((2 * R, 2 * R), np.uint8); g[xy[k, 1] + R, xy[k, 0] + R] = 1; return g


def rad(R): a = np.arange(-R, R); a = np.where(a < 0, -a - 1, a); return np.maximum(a[:, None], a[None, :])


def pieces(g, labels=False):
    """sizes of the connected pieces of a 0/1 array (largest first), or with labels=True the label array and sizes by label"""
    H, Wd = g.shape; lab = np.zeros(g.shape, np.int32); sizes = []; ys, xs = np.nonzero(g)
    for y0, x0 in zip(ys.tolist(), xs.tolist()):
        if lab[y0, x0]: continue
        st = [(y0, x0)]; lab[y0, x0] = len(sizes) + 1; m = 0
        while st:
            y, x = st.pop(); m += 1
            for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                a, b = y + dy, x + dx
                if 0 <= a < H and 0 <= b < Wd and g[a, b] and not lab[a, b]: lab[a, b] = len(sizes) + 1; st.append((a, b))
        sizes.append(m)
    return (lab, sizes) if labels else np.array(sorted(sizes, reverse=True))


def features(g, levels):
    """numbers describing a picture: cells by distance, local patterns in each level's ring, pieces"""
    R = g.shape[0] // 2; r = rad(R); f = {}; gi = g.astype(np.int32)
    p = np.pad(gi, 1); nb = p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2]
    full = np.zeros_like(gi); full[:-1, :-1] = gi[:-1, :-1] * gi[1:, :-1] * gi[:-1, 1:] * gi[1:, 1:]
    edges = [B * 2 ** (L + k / 4) for L in range(1, levels + 1) for k in range(4)] + [B * 2 ** (levels + 1)]
    for a, b in zip(edges[:-1], edges[1:]): f[f"cells {a:.0f}-{b:.0f}"] = float(gi[(r >= a) & (r < b)].sum())
    for L in range(1, levels + 1):
        a = B * 2 ** L; m = (r >= a) & (r < 2 * a) & (r < R - 1); c = gi * m; n = max(1, c.sum())
        f[f"L{L} cells"] = float(c.sum())
        for k in (1, 2, 3, 4): f[f"L{L} with {k} neighbours"] = float(((nb == k) * c).sum())
        f[f"L{L} empty sites next to the shape"] = float(((nb > 0) * (1 - gi) * m).sum())
        f[f"L{L} full 2x2"] = float((full * m).sum())
        for s in (4, 16):
            if a // s >= 2: q = (gi * m).reshape(2 * R // s, s, 2 * R // s, s).max((1, 3)); f[f"L{L} boxes of {s} holding cells"] = float(q.sum())
        w = 2 if L > 1 else 1                                                      # a thin shell each side of the seam
        f[f"L{L} seam: cells within {w} inside"] = float(gi[(r >= a - w) & (r < a)].sum()); f[f"L{L} seam: cells within {w} outside"] = float(gi[(r >= a) & (r < a + w)].sum())
        f[f"L{L} seam: joins across"] = float(((r == a - 1) * gi * np.maximum.reduce([np.roll(gi * (r == a), s_, ax) for s_, ax in ((1, 0), (-1, 0), (1, 1), (-1, 1))])).sum())
    # pieces are counted in the inner three quarters of the picture, so that its edge cuts through the middle of the last ring
    Rf = R; R = 3 * R // 4; g = g[Rf - R:Rf + R, Rf - R:Rf + R]; r = rad(R)
    lab, s = pieces(g, True); rim = np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]])); rim = set(rim[rim > 0].tolist())
    fl = np.array([m for k, m in enumerate(s) if k + 1 not in rim and lab[R, R] != k + 1] or [0])       # pieces that neither hold the root nor reach the edge of the picture
    f["floating pieces"] = float((fl > 0).sum()); f["floating pieces of 10 or more"] = float((fl >= 10).sum()); f["cells in floating pieces"] = float(fl.sum())
    for L in range(1, levels + 1):
        a = B * 2 ** L; ids = np.unique(lab[(r >= a) & (r < 2 * a)]); inn = set(np.unique(lab[r < a]).tolist())
        if a >= R: break
        f[f"L{L} floating pieces first met in this ring"] = float(sum(1 for k in ids.tolist() if k > 0 and k not in rim and k != lab[R, R] and k not in inn and s[k - 1] >= 3))
    s = np.array(sorted(s, reverse=True)); f["pieces"] = float(len(s)); f["pieces of 2 to 9 cells"] = float(((s >= 2) & (s < 10)).sum()); f["pieces of 10 to 99"] = float(((s >= 10) & (s < 100)).sum())
    f["pieces of 100 or more"] = float((s >= 100).sum()); f["single cells"] = float((s == 1).sum()); f["share in the largest piece"] = float(s[0] / s.sum()); f["cells"] = float(gi.sum())
    return f


def one(job):
    k, fn, c, levels, seed, extra = job; tmp = os.path.join(ROOT, "zoomtest", "e2e", f"tmp_{os.getpid()}")
    if subprocess.run([os.path.join(ROOT, "ball"), fn, tmp + "_b.bin", str(c[0]), str(c[1]), "32", "16000"], stderr=subprocess.DEVNULL).returncode: return None
    o = subprocess.run(["node", os.path.join(ROOT, "analysis", "out_run.js"), tmp + "_b.bin", tmp, str(levels), str(seed + k), str(levels)] + extra, capture_output=True, text=True)
    if o.returncode: return None
    R = B * 2 ** (levels + 1); gen = raster(read(f"{tmp}_l{levels}.bin"), R)
    return features(gen, levels), [json.loads(l) for l in o.stdout.strip().split("\n")], gen if k < 12 else None


def main():
    R = int(sys.argv[1]); levels = int(sys.argv[2]) if len(sys.argv) > 2 else 4; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    J = int(sys.argv[4]) if len(sys.argv) > 4 else 4; tag = sys.argv[5] if len(sys.argv) > 5 else "base"; extra = sys.argv[6:]
    rng = np.random.default_rng(seed); os.makedirs(os.path.join(ROOT, "zoomtest", "e2e"), exist_ok=True)
    snaps = sorted(f for f in os.listdir(os.path.join(ROOT, "prod")) if f.startswith("snap_") and "n1024000" in f and f.endswith(".bin"))
    jobs = []; real = []; per = -(-R // len(snaps)); Rr = B * 2 ** (levels + 1); pics = []
    for s in snaps:
        fn = os.path.join(ROOT, "prod", s); T = read(fn)
        for _ in range(per):
            c = T[rng.integers(len(T))]; g = raster(T - c, Rr); real.append(features(g, levels)); jobs.append((len(jobs), fn, (int(c[0]), int(c[1])), levels, seed * 100000, extra))
            if len(jobs) <= 12: pics.append(g)
    with Pool(J) as p: res = p.map(one, jobs, chunksize=1)
    keep = [i for i, r in enumerate(res) if r is not None]; names = list(real[0].keys())
    Rm = np.array([[real[i][k] for k in names] for i in keep]); Gm = np.array([[res[i][0][k] for k in names] for i in keep]); n = len(keep)
    info = [res[i][1] for i in keep]
    np.savez_compressed(os.path.join(ROOT, "zoomtest", "e2e", tag + "_pics.npz"), real=np.array(pics), gen=np.array([res[i][2] for i in range(min(12, len(res))) if res[i] is not None and res[i][2] is not None]))
    rows = []
    for j, k in enumerate(names):
        per_cell = k.startswith("L") and not k.endswith(" cells") and "seam" not in k and "boxes" not in k
        if per_cell:                                                                 # a ratio to the cells in that ring
            L = k.split()[0]; d = names.index(L + " cells")
            boot = np.array([Gm[ix, j].sum() / Gm[ix, d].sum() - Rm[ix, j].sum() / Rm[ix, d].sum() for ix in (rng.integers(n, size=n) for _ in range(300))])
            g, r_ = Gm[:, j].sum() / Gm[:, d].sum(), Rm[:, j].sum() / Rm[:, d].sum(); rows.append((k + " (per cell)", r_, g, g / r_ - 1, (g - r_) / boot.std()))
        else:
            d = Gm[:, j] - Rm[:, j]; se = d.std(ddof=1) / np.sqrt(n); r_, g = Rm[:, j].mean(), Gm[:, j].mean(); rows.append((k, r_, g, g / r_ - 1 if r_ else 0.0, (g - r_) / se if se > 0 else 0.0))
    json.dump({"roots": n, "levels": levels, "extra": extra, "rows": rows, "info": info}, open(os.path.join(ROOT, "zoomtest", "e2e", tag + ".json"), "w"))
    if tag == "final":                                                               # the run the page and the README quote
        mean = lambda k: [float(np.mean([i[L][k] for i in info])) for L in range(levels)]
        json.dump({"roots": n, "samples": len(snaps), "levels": levels, "seed": seed, "rows": rows, "seconds": [m / 1000 for m in mean("ms")], "links": mean("links"), "unjoined": mean("unjoined"),
                   "fineParts": mean("fineParts"), "parts": mean("parts")}, open(os.path.join(ROOT, "analysis", "e2e_test.json"), "w"), indent=1)
        np.savez_compressed(os.path.join(ROOT, "analysis", "e2e_pics.npz"), real=np.packbits(np.array(pics, bool), axis=-1), made=np.packbits(np.array([res[i][2] for i in range(min(12, len(res))) if res[i] is not None and res[i][2] is not None], bool), axis=-1))
    print(f"{tag}: {n} roots, {levels} levels, {' '.join(extra)}")
    print(f"  {'':46s} {'real':>10s} {'made':>10s} {'diff':>8s} {'s.e.':>6s}")
    for k, r_, g, d, z in rows: print(f"  {k:46s} {r_:10.4g} {g:10.4g} {100 * d:+7.1f}% {z:+6.1f}")
    z = np.array([r[4] for r in rows]); print(f"  beyond 2 s.e.: {int((np.abs(z) > 2).sum())} of {len(z)}; largest {np.abs(z).max():.1f}")
    for L in range(levels): print(f"  level {L + 1}: {np.mean([i[L]['ms'] for i in info]):.0f} ms, links {np.mean([i[L]['links'] for i in info]):.1f} cells, unjoined {np.mean([i[L]['unjoined'] for i in info]):.2f}, parts in the finer box {np.mean([i[L]['fineParts'] for i in info]):.2f}, after halving {np.mean([i[L]['partsHalved'] for i in info]):.2f}, after thinning {np.mean([i[L]['partsThinned'] for i in info]):.2f}, in the new world {np.mean([i[L]['parts'] for i in info]):.2f}")


if __name__ == "__main__":
    main()
