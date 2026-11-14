"""Exact averages over ALL fixed polyominoes of size n (brute force), for testing the sampler:
site perimeter t, removable (non-cut) cells k, holes, squared radius of gyration."""
import sys
from fractions import Fraction

def canon(cells):
    mx = min(x for x, y in cells); my = min(y for x, y in cells)
    return frozenset((x - mx, y - my) for x, y in cells)

def grow(shapes):
    out = set()
    for s in shapes:
        for (x, y) in s:
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                c = (x + dx, y + dy)
                if c not in s:
                    out.add(canon(s | {c}))
    return out

def connected(s):
    s = set(s); start = next(iter(s)); seen = {start}; st = [start]
    while st:
        x, y = st.pop()
        for c in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
            if c in s and c not in seen:
                seen.add(c); st.append(c)
    return len(seen) == len(s)

def stats(s):
    per = set()
    for (x, y) in s:
        for c in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
            if c not in s: per.add(c)
    k = sum(1 for c in s if len(s) == 1 or connected(s - {c}))
    n = len(s)
    sx = sum(x for x, y in s); sy = sum(y for x, y in s); sq = sum(x * x + y * y for x, y in s)
    rg2 = Fraction(sq, n) - Fraction(sx * sx + sy * sy, n * n)
    return len(per), k, rg2

if __name__ == "__main__":
    N = int(sys.argv[1])
    shapes = {frozenset({(0, 0)})}
    for n in range(2, N + 1):
        shapes = grow(shapes)
        T = K = 0; R = Fraction(0)
        for s in shapes:
            t, k, r = stats(s); T += t; K += k; R += r
        A = len(shapes)
        print(f"n={n} A={A} <t>={T / A:.8f} <k>={K / A:.8f} <rg2>={float(R / A):.8f}", flush=True)
