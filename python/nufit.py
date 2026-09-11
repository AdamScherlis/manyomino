"""Fit the radius-of-gyration exponent nu: <Rg> ~ n^nu.

Takes (n, csv[, csv...]) groups; pools chains at the same n (they must be
independent runs), then does a weighted least-squares fit of log<Rg> vs
log n.  Expect nu ~ 0.64.

Usage: python3 nufit.py n1:file1[,file1b] n2:file2 ...
"""

import math
import sys

from analyze import analyze


def pooled_rg(paths, burn_frac=0.5):
    means, ses = [], []
    for p in paths:
        r = analyze(p, burn_frac)
        means.append(r["rg_mean"])
        ses.append(r["rg_stderr"])
    wts = [1.0 / s**2 for s in ses]
    wsum = sum(wts)
    mean = sum(m * w for m, w in zip(means, wts)) / wsum
    se = math.sqrt(1.0 / wsum)
    return mean, se


def wls(xs, ys, ss):
    """Weighted least squares y = a + b x; returns (a, b, sb)."""
    wts = [1.0 / s**2 for s in ss]
    sw = sum(wts)
    sx = sum(w * x for w, x in zip(wts, xs))
    sy = sum(w * y for w, y in zip(wts, ys))
    sxx = sum(w * x * x for w, x in zip(wts, xs))
    sxy = sum(w * x * y for w, x, y in zip(wts, xs, ys))
    d = sw * sxx - sx * sx
    b = (sw * sxy - sx * sy) / d
    a = (sxx * sy - sx * sxy) / d
    sb = math.sqrt(sw / d)
    return a, b, sb


def run(groups, burn_frac=0.5):
    xs, ys, ss = [], [], []
    for n, paths in groups:
        rg, se = pooled_rg(paths, burn_frac)
        print(f"n={n:6d}  <Rg> = {rg:9.3f} +- {se:.3f}   ({len(paths)} chain(s))")
        xs.append(math.log(n))
        ys.append(math.log(rg))
        ss.append(se / rg)
    a, b, sb = wls(xs, ys, ss)
    print(f"\nnu = {b:.4f} +- {sb:.4f}   (amplitude {math.exp(a):.4f})")
    return b, sb


if __name__ == "__main__":
    groups = []
    for arg in sys.argv[1:]:
        n, files = arg.split(":")
        groups.append((int(n), files.split(",")))
    run(groups)
