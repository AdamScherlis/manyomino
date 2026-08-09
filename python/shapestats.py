"""Universal shape statistics from converged snapshots: gyration-tensor
asphericity, cycle density, perimeter density.

Asphericity A = <(l1 - l2)^2 / (l1 + l2)^2> over samples, where l1 >= l2
are the gyration tensor eigenvalues -- a universal ratio for the branched
polymer class.  Usage: python3 shapestats.py"""
import glob
import math
import re
from collections import defaultdict


def load(p):
    return [tuple(map(int, l.split())) for l in open(p)]


def stats(cells):
    n = len(cells)
    cs = set(cells)
    mx = sum(x for x, y in cells) / n
    my = sum(y for x, y in cells) / n
    sxx = sum((x - mx) ** 2 for x, y in cells) / n
    syy = sum((y - my) ** 2 for x, y in cells) / n
    sxy = sum((x - mx) * (y - my) for x, y in cells) / n
    tr, det = sxx + syy, sxx * syy - sxy * sxy
    disc = math.sqrt(max(tr * tr / 4 - det, 0.0))
    l1, l2 = tr / 2 + disc, tr / 2 - disc
    asph = (l1 - l2) ** 2 / (l1 + l2) ** 2
    e = sum(1 for (x, y) in cells for nb in ((x + 1, y), (x, y + 1)) if nb in cs)
    cyc = e - n + 1
    per = len({nb for (x, y) in cells
               for nb in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1))
               if nb not in cs})
    return asph, cyc / n, per / n, l1 / l2


by_n = defaultdict(list)
for f in glob.glob("../gallery/raw/n*_*.txt"):
    m = re.search(r"n(\d+)_(bar|rect|bar2)_(\d+)\.txt$", f)
    if not m:
        continue
    n = int(m.group(1))
    by_n[n].append(f)

print(f"{'n':>7} {'#':>3} {'asphericity':>12} {'l1/l2':>7} {'cycles/n':>9} {'perim/n':>8}")
for n in sorted(by_n):
    if n < 1000:
        continue
    vals = [stats(load(f)) for f in by_n[n]]
    k = len(vals)
    def m_se(i):
        xs = [v[i] for v in vals]
        m = sum(xs) / k
        se = (sum((x - m) ** 2 for x in xs) / (k * max(k - 1, 1))) ** 0.5
        return m, se
    a, ase = m_se(0)
    c, cse = m_se(1)
    p, pse = m_se(2)
    r, _ = m_se(3)
    print(f"{n:>7} {k:>3} {a:8.4f}±{ase:.4f} {r:7.2f} {c:7.4f}±{cse:.4f} {p:6.4f}±{pse:.4f}")
