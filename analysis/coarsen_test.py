"""Can a coarsened sample be told from a sampled one?

scripts/coarsen_chain.sh coarsens real samples of 1,024,000 cells four times over (each 2x2
block holding any cell becomes one cell, and the local walk of src/zoom.cpp thins the result),
which leaves about 13,500 cells.  Here the 60 observables of those shapes are set against the
same observables in equilibrium at 13,650 cells (zoomtest/ref, from scripts/zoom_test.sh), as
analysis/zoom_test.py does in the other direction.

  python3 analysis/coarsen_test.py [T ...]      # reads zoomtest/coarse/c{T}_*.feat, writes analysis/coarsen_test.json
"""
import glob, json, math, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zoom_test as zt

ROOT = zt.ROOT; N0 = 13650


def reference():
    B = []; A = []
    for fn in sorted(glob.glob(os.path.join(ROOT, "zoomtest", "ref", "r?.feat"))):
        a = np.fromfile(fn); a = a[: len(a) // 61 * 61].reshape(-1, 61)[:, 1:]; L = len(a) // 8; A.append(a)
        B += [a[i * L:(i + 1) * L].mean(0) for i in range(8)]
    B = np.array(B); A = np.vstack(A)
    return B.mean(0), B.std(0, ddof=1) / math.sqrt(len(B)), A.std(0, ddof=1), len(A)


def main():
    Ts = [int(a) for a in sys.argv[1:]] or sorted({int(re.search(r"c(\d+)_", os.path.basename(f)).group(1)) for f in glob.glob(os.path.join(ROOT, "zoomtest", "coarse", "c*_*.feat"))})
    names, mu1, err1, sd1, nb, nr = zt.reference(); expo = zt.exponents(names, mu1, err1)
    mu, mu_err, sd, readings = reference(); loc = np.array([nm in zt.LOCAL for nm in names])
    out = {"names": names, "local": list(zt.LOCAL), "readings": readings, "T": {}}
    for T in Ts:
        rows = []; ns = []; growth = {}
        for fn in sorted(glob.glob(os.path.join(ROOT, "zoomtest", "coarse", f"c{T}_*.feat"))):
            a = np.fromfile(fn)
            if a.size != 61: continue
            n = int(np.fromfile(fn[:-5] + ".bin", dtype="<i4", count=1)[0]); ns.append(n); rows.append(a[1:] * (N0 / n) ** expo)
            k = 0
            for line in open(fn[:-5] + ".log"):
                g = re.match(r"zoom: (\d+) -> (\d+) cells", line)
                if g: k += 1; growth.setdefault(k, []).append(int(g.group(2)))
        if len(rows) < 4: continue
        Z = np.array(rows); ns = np.array(ns); zm, zs = Z.mean(0), Z.std(0, ddof=1); se = np.sqrt(sd ** 2 / len(Z) + mu_err ** 2)
        with np.errstate(divide="ignore", invalid="ignore"): z = np.where(se > 0, (zm - mu) / se, 0.0)
        # cells kept per level: level k's count over level k-1's (1,024,000 before the first)
        prev = np.full(len(growth[1]), 1024000.0); keep = []
        for k in sorted(growth): cur = np.array(growth[k], float); r = prev / cur; keep.append([float(r.mean()), float(r.std(ddof=1) / math.sqrt(len(r)))]); prev = cur
        print(f"\nT={T}: {len(Z)} coarsened shapes of {ns.mean():,.0f} ± {ns.std(ddof=1):,.0f} cells; cells before over cells after, per level: " + ", ".join(f"{a:.4f}±{b:.4f}" for a, b in keep))
        print(f"  beyond 2 standard errors: {int((np.abs(z) > 2).sum())} of {len(names)}; largest {np.abs(z).max():.1f}; per-cell counts ({int(loc.sum())}): rms {math.sqrt((z[loc] ** 2).mean()):.2f}, largest {np.abs(z[loc]).max():.1f}; "
              f"the rest: rms {math.sqrt((z[~loc] ** 2).mean()):.2f}, largest {np.abs(z[~loc]).max():.1f}, beyond 2: {int((np.abs(z[~loc]) > 2).sum())}")
        for j in np.argsort(-np.abs(z))[:8]:
            print(f"    {names[j]:>10s}  z {z[j]:+5.1f}   coarsened {zm[j]:.6g} ± {zs[j] / math.sqrt(len(Z)):.3g}   equilibrium {mu[j]:.6g} ± {mu_err[j]:.3g}   ({100 * (zm[j] / mu[j] - 1) if abs(mu[j]) > 8 * mu_err[j] else float('nan'):+.2f}%)")
        out["T"][T] = {"shapes": len(Z), "cells": [float(ns.mean()), float(ns.std(ddof=1))], "z": [float(v) for v in z], "keep": keep, "mean": [float(v) for v in zm], "eq": [float(v) for v in mu], "eq_err": [float(v) for v in mu_err]}
    json.dump(out, open(os.path.join(ROOT, "analysis", "coarsen_test.json"), "w"))


if __name__ == "__main__":
    main()
