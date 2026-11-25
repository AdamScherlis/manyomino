"""Could each dot of a far level be given a number of real cells drawn from a table?

analysis/mass_test.py shows that far out the worlds' cells ("dots") light too few pixels and light them evenly,
where the real cells a dot stands for would light more and unevenly.  The cheapest repair would be to draw,
for every dot and every empty site beside one, a number of real cells from a table of what real shapes have
there, given the dots around it (the eight neighbours exactly, and how many of the sixteen sites beyond).  This
learns such tables from half the roots, for dots 2, 4 and 8 real cells wide, draws from them for the other
half, and counts lit pixels and the spread of brightness in pixels 1, 2 and 4 dots wide: for the real cells,
for the draws, and for the dots as the page shows them now.

  python3 analysis/mass_model.py ROOTS [SEED=1] [JOBS=4]       # writes analysis/mass_model.json
"""
import sys, os, json, subprocess
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e2e_test import read, ROOT


def one(job):
    k, fn, c, seed = job; d = os.path.join(ROOT, "zoomtest", "mass"); tmp = os.path.join(d, f"ball_{os.getpid()}.bin"); out = os.path.join(d, f"w_{k}.bin")
    if subprocess.run([os.path.join(ROOT, "ball"), fn, tmp, str(c[0]), str(c[1]), "256", "60000"], stderr=subprocess.DEVNULL).returncode: return None
    return out if subprocess.run(["node", os.path.join(ROOT, "analysis", "mass_dump.js"), tmp, str(seed + k), out], capture_output=True).returncode == 0 else None


def load(fn):
    a = np.fromfile(fn, dtype=np.int32); ws = []; p = 0
    for _ in range(4): n = a[p]; ws.append((a[p + 1:p + 1 + n], a[p + 1 + n:p + 1 + 2 * n])); p += 1 + 2 * n
    return ws


def grids(ws, d):
    """Dots, the key of each site (what is around it), and real cells per site, out to 192 real cells from the root."""
    R = 192 >> d; G = 2 * R + 4; occ = np.zeros((G, G), np.int32); m = np.zeros((G, G), np.int32)
    x, y = ws[d]; u, v = x + R + 2, y + R + 2; ok = (u >= 0) & (v >= 0) & (u < G) & (v < G); occ[v[ok], u[ok]] = 1
    x, y = ws[0]; u, v = (x >> d) + R + 2, (y >> d) + R + 2; ok = (u >= 0) & (v >= 0) & (u < G) & (v < G); np.add.at(m, (v[ok], u[ok]), 1)
    nb = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if a or b]; k = sum(np.roll(np.roll(occ, a, 0), b, 1) << i for i, (a, b) in enumerate(nb))
    k2 = sum(np.roll(np.roll(occ, a, 0), b, 1) for a in range(-2, 3) for b in range(-2, 3) if max(abs(a), abs(b)) == 2)
    key = occ * 1024 + k * 4 + np.minimum(k2, 6) // 2
    return occ[2:-2, 2:-2], key[2:-2, 2:-2], m[2:-2, 2:-2]


def pool(f, p): G = f.shape[0] // p * p; return f[:G, :G].reshape(G // p, p, G // p, p).sum((1, 3))


def main():
    R = int(sys.argv[1]); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1; J = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    rng = np.random.default_rng(seed); os.makedirs(os.path.join(ROOT, "zoomtest", "mass"), exist_ok=True)
    snaps = sorted(f for f in os.listdir(os.path.join(ROOT, "prod")) if f.startswith("snap_") and "n1024000" in f and f.endswith(".bin")); jobs = []
    for s in snaps:
        fn = os.path.join(ROOT, "prod", s); T = read(fn)
        for _ in range(-(-R // len(snaps))): c = T[rng.integers(len(T))]; jobs.append((len(jobs), fn, (int(c[0]), int(c[1])), seed * 100000))
    with Pool(J) as p: files = [r for r in p.map(one, jobs, chunksize=1) if r]
    learn, test = files[0::2], files[1::2]; out = {"roots learnt from": len(learn), "roots drawn for": len(test), "samples": len(snaps), "rows": []}
    for d in (1, 2, 3):
        tab = {}
        for fn in learn:
            occ, key, m = grids(load(fn), d); o = np.argsort(key.ravel(), kind="stable"); ks = key.ravel()[o]; ms = m.ravel()[o]; cut = np.flatnonzero(np.diff(ks)) + 1
            for kk, mm in zip(ks[np.r_[0, cut]], np.split(ms, cut)): tab.setdefault(int(kk), []).append(mm)
        tab = {k: np.concatenate(v) for k, v in tab.items()}; res = {}
        for fn in test:
            occ, key, m = grids(load(fn), d); mod = np.zeros(m.size, np.int32); kf = key.ravel()
            for kv in np.unique(kf):
                t = tab.get(int(kv))
                if t is not None: sel = kf == kv; mod[sel] = rng.choice(t, int(sel.sum()))
            for name, f in (("real", m), ("drawn", mod.reshape(m.shape)), ("dots", occ)):
                for p in (1, 2, 4):
                    g = pool(f.astype(float), p); lit = g > 0; b = g[lit] / g[lit].mean(); res.setdefault((name, p), []).append((lit.sum(), b.std()))
        for p in (1, 2, 4):
            row = {"dot, in real cells": 2 ** d, "pixel, in dots": p}
            for name in ("real", "drawn", "dots"): a = np.array(res[(name, p)], float); row[name] = {"lit pixels": float(a[:, 0].mean()), "spread of brightness": float(a[:, 1].mean())}
            for name in ("drawn", "dots"): row[name]["lit, against real"] = row[name]["lit pixels"] / row["real"]["lit pixels"] - 1
            out["rows"].append(row)
    json.dump(out, open(os.path.join(ROOT, "analysis", "mass_model.json"), "w"), indent=1)
    print(f"tables from {len(learn)} roots, drawn for {len(test)} others")
    for r in out["rows"]: print(f"  dot {r['dot, in real cells']} cells, pixel {r['pixel, in dots']} dots: lit pixels real {r['real']['lit pixels']:.0f}, drawn {100 * r['drawn']['lit, against real']:+.1f}%, dots {100 * r['dots']['lit, against real']:+.1f}%; "
                               f"spread of brightness real {r['real']['spread of brightness']:.2f}, drawn {r['drawn']['spread of brightness']:.2f}, dots {r['dots']['spread of brightness']:.2f}")


if __name__ == "__main__": main()
