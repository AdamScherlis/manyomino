"""Validate the tilted sampler against the EXACT reweighted distribution:
pi(shape) ~ exp(-beta * Rg^2) over the enumerated shape list.

Usage: python3 validate_beta.py counts.txt n beta"""
import math
import sys

from enumerate_fixed import enumerate_fixed
from validate_fast import pack
from stats import chi2_sf


def rg2(s):
    n = len(s)
    sx = sum(x for x, y in s); sy = sum(y for x, y in s)
    sx2 = sum(x * x for x, y in s); sy2 = sum(y * y for x, y in s)
    return (sx2 + sy2) / n - (sx / n) ** 2 - (sy / n) ** 2


def run(path, n, beta):
    shapes = sorted(enumerate_fixed(n)[n])
    logw = [-beta * rg2(s) for s in shapes]
    m = max(logw)
    w = [math.exp(l - m) for l in logw]
    z = sum(w)
    probs = [x / z for x in w]
    counts = {}
    total = 0
    for line in open(path):
        k, v = line.split(); counts[int(k)] = int(v); total += int(v)
    chi2 = 0.0
    for s, p in zip(shapes, probs):
        e = p * total
        o = counts.get(pack(s), 0)
        if e > 0:
            chi2 += (o - e) ** 2 / e
    dof = len(shapes) - 1
    pv = chi2_sf(chi2, dof)
    print(f"n={n} beta={beta}: chi2={chi2:.1f} dof={dof} p={pv:.4f} "
          f"(min/max prob ratio {min(probs)/max(probs):.2e})")
    return pv


if __name__ == "__main__":
    p = run(sys.argv[1], int(sys.argv[2]), float(sys.argv[3]))
    print("PASS" if p > 0.001 else "FAIL")
    sys.exit(0 if p > 0.001 else 1)
