"""Can a refined shape be told from a sampled one?

scripts/zoom_chain.sh refines samples of 13,650 cells four times over, to about 1,030,000
cells (the count is left wherever the refinement puts it).  Here the 60 observables of those
shapes are set against the same observables along the equilibrium runs at 1,024,000 cells
(prod/ts_?_n1024000.part*.feat).  An observable that grows with the size of the shape is first
brought to 1,024,000 cells by its own power law, measured between the equilibrium runs at 64,000
and at 1,024,000 cells; the refined shapes are within a few per cent of that size, so the
correction is small beside their spread.

For each observable: the mean over the refined shapes minus the equilibrium mean, in units of
the standard error that difference would have if the refined shapes were equilibrium samples
(the equilibrium spread over the square root of their number, combined with the error of the
equilibrium mean from the scatter of blocks of the runs); and the spread of the refined shapes
relative to the equilibrium spread.

  python3 analysis/zoom_test.py [T ...]      # reads zoomtest/z{T}_*.feat, writes analysis/zoom_test.json
"""
import glob, json, math, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hero_slow import chains, N

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
BLOCK = 150          # moves per cell per block of the equilibrium runs, several relaxation times
LOCAL = ("perim", "holes", "cycles", "leaves", "deg2", "deg3", "deg4", "straight", "blocks22", "noncut", "blob_cells", "blob_count")


def reference():
    names, ch = chains(); blocks = []; allrows = []
    for t, c in ch.items():
        X = c["X"]; L = max(2, int(round(BLOCK / c["step"]))); allrows.append(X)
        blocks += [X[i:i + L].mean(0) for i in range(0, len(X) - L + 1, L)]
    B = np.array(blocks); A = np.vstack(allrows)
    return names, B.mean(0), B.std(0, ddof=1) / math.sqrt(len(B)), A.std(0, ddof=1), len(B), len(A)


def exponents(names, mu, mu_err):
    """power of the cell count that each observable's mean follows, from 64,000 to 1,024,000 cells
    (0 where the mean is not clearly away from zero at both sizes)"""
    rows = []
    for fn in sorted(glob.glob(os.path.join(ROOT, "conv", "u64000_?.feat"))):
        nm = open(fn + ".names").read().split(); a = np.memmap(fn, dtype=np.float64, mode="r"); a = a[: len(a) // len(nm) * len(nm)].reshape(-1, len(nm))
        assert nm[1:] == names; rows.append(np.array(a[len(a) // 10::7, 1:]))
    m64 = np.vstack(rows).mean(0); out = np.zeros(len(names))
    for j in range(len(names)):
        if abs(mu[j]) > 8 * mu_err[j] and m64[j] * mu[j] > 0: out[j] = math.log(mu[j] / m64[j]) / math.log(16)
    return out


def refined(T, expo):
    rows = []; ns = []
    for fn in sorted(glob.glob(os.path.join(ROOT, "zoomtest", f"z{T}_*.feat"))):
        a = np.fromfile(fn)
        if a.size != 61: continue
        n = int(np.fromfile(fn[:-5] + ".bin", dtype="<i4", count=1)[0]); ns.append(n)
        rows.append(a[1:] * (N / n) ** expo)
    return np.array(rows), np.array(ns)


def levels(T):
    """growth of the cell count at each level, from the logs"""
    m = {}
    for fn in glob.glob(os.path.join(ROOT, "zoomtest", f"z{T}_*.log")):
        k = 0
        for line in open(fn):
            g = re.match(r"zoom: (\d+) -> (\d+) cells \(x([\d.]+) of the original (\d+)\)", line)
            if g: k += 1; m.setdefault(k, []).append(float(g.group(3)))
    return m


def main():
    Ts = [int(a) for a in sys.argv[1:]] or sorted({int(re.search(r"z(\d+)_", f).group(1)) for f in glob.glob(os.path.join(ROOT, "zoomtest", "z*_*.feat"))})
    names, mu, mu_err, sd, nb, nr = reference(); expo = exponents(names, mu, mu_err)
    out = {"names": names, "local": list(LOCAL), "blocks": nb, "readings": nr, "exponent": [float(v) for v in expo], "T": {}}
    print(f"equilibrium: {nr} readings in {nb} blocks of {BLOCK} moves per cell")
    print("grow with size: " + " ".join(f"{nm}^{e:.2f}" for nm, e in zip(names, expo) if abs(e) > 0.02))
    for T in Ts:
        Z, ns = refined(T, expo)
        if len(Z) < 4: continue
        zm, zs = Z.mean(0), Z.std(0, ddof=1); se = np.sqrt(sd ** 2 / len(Z) + mu_err ** 2)
        with np.errstate(divide="ignore", invalid="ignore"):
            z = np.where(se > 0, (zm - mu) / se, 0.0); spread = np.where(sd > 0, zs / sd, 1.0); shift = np.where(sd > 0, (zm - mu) / sd, 0.0)
        lv = levels(T); order = np.argsort(-np.abs(z))
        print(f"\nT={T}: {len(Z)} refined shapes; growth per level " + ", ".join(f"{np.mean(lv[k]):.4f}±{np.std(lv[k], ddof=1) / math.sqrt(len(lv[k])):.4f}" for k in (1, 2, 3, 4) if k in lv)
              + f"; cells {ns.mean():,.0f} ± {ns.std(ddof=1):,.0f}")
        print(f"  beyond 2 standard errors: {int((np.abs(z) > 2).sum())} of {len(names)}; beyond 3: {int((np.abs(z) > 3).sum())}; largest {np.abs(z).max():.1f}; rms {math.sqrt((z ** 2).mean()):.2f}")
        loc = np.array([nm in LOCAL for nm in names])
        print(f"  counts of small features per cell ({int(loc.sum())}): rms {math.sqrt((z[loc] ** 2).mean()):.2f}, largest {np.abs(z[loc]).max():.1f}; the rest ({int((~loc).sum())}): rms {math.sqrt((z[~loc] ** 2).mean()):.2f}, largest {np.abs(z[~loc]).max():.1f}")
        for j in order[:12]: print(f"    {names[j]:>10s}  z {z[j]:+5.1f}   shift {shift[j]:+.2f} sd   refined {zm[j]:.6g} ± {zs[j] / math.sqrt(len(Z)):.3g}   equilibrium {mu[j]:.6g} ± {mu_err[j]:.3g}   ({100 * (zm[j] / mu[j] - 1) if abs(mu[j]) > 8 * mu_err[j] else float('nan'):+.2f}%)   spread ratio {spread[j]:.2f}")
        out["T"][T] = {"shapes": len(Z), "z": [float(v) for v in z], "shift": [float(v) for v in shift], "spread": [float(v) for v in spread],
                       "growth": {str(k): [float(np.mean(lv[k])), float(np.std(lv[k], ddof=1) / math.sqrt(len(lv[k])))] for k in (1, 2, 3, 4) if k in lv},
                       "cells": [float(ns.mean()), float(ns.std(ddof=1))], "mean": [float(v) for v in zm], "sem": [float(v) for v in zs / math.sqrt(len(Z))],
                       "eq": [float(v) for v in mu], "eq_err": [float(v) for v in mu_err]}
    json.dump(out, open(os.path.join(ROOT, "analysis", "zoom_test.json"), "w"))


if __name__ == "__main__":
    main()
