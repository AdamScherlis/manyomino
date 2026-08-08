"""Exact census of fixed n-polyominoes as packed ints (memory-light), then
compare cycle-count and diagonal-near-miss distributions against sampler
output.  Usage: python3 cyclestats12.py counts.txt n"""
import sys
from collections import Counter

def decode(key, n):
    cells = []
    for _ in range(n):
        b = key & 0xFF
        key >>= 8
        cells.append((b >> 4, b & 0xF))
    return cells

def pack(cells):
    bs = sorted((x << 4) | y for x, y in cells)
    k = 0
    for b in bs:
        k = (k << 8) | b
    return k

def canon(cells):
    mx = min(x for x, y in cells); my = min(y for x, y in cells)
    return pack([(x - mx, y - my) for x, y in cells])

def grow(nmax):
    shapes = {pack([(0, 0)])}
    for n in range(2, nmax + 1):
        nxt = set()
        for key in shapes:
            cells = decode(key, n - 1)
            cs = set(cells)
            per = set()
            for (x, y) in cells:
                for nb in ((x+1,y),(x-1,y),(x,y+1),(x,y-1)):
                    if nb not in cs:
                        per.add(nb)
            for p in per:
                nxt.add(canon(cells + [p]))
        shapes = nxt
        print(f"n={n}: {len(shapes)}", file=sys.stderr, flush=True)
    return shapes

def cycles(cells):
    cs = set(cells)
    e = sum(1 for (x,y) in cells for nb in ((x+1,y),(x,y+1)) if nb in cs)
    return e - len(cells) + 1

def near_misses(cells):
    cs = set(cells)
    k = 0
    for (x, y) in cells:
        for dx, dy in ((1,1),(1,-1)):
            if (x+dx,y+dy) in cs and (x+dx,y) not in cs and (x,y+dy) not in cs:
                k += 1
    return k

if __name__ == "__main__":
    path, n = sys.argv[1], int(sys.argv[2])
    shapes = grow(n)
    exp_cyc, exp_nm = Counter(), Counter()
    for key in shapes:
        cells = decode(key, n)
        exp_cyc[cycles(cells)] += 1
        exp_nm[near_misses(cells)] += 1
    obs_cyc, obs_nm = Counter(), Counter()
    total = 0
    unknown = 0
    with open(path) as f:
        for line in f:
            k, v = line.split(); k, v = int(k), int(v)
            if k not in shapes:
                unknown += v
                continue
            cells = decode(k, n)
            obs_cyc[cycles(cells)] += v
            obs_nm[near_misses(cells)] += v
            total += v
    print(f"shapes={len(shapes)} observations={total} unknown_shapes={unknown}")
    print("\ncycle count: exact freq vs sampled")
    for c in sorted(exp_cyc):
        e = exp_cyc[c] / len(shapes); o = obs_cyc[c] / total
        se = (e * (1 - e) / total) ** 0.5
        print(f"  c={c}: exact {e:.6f}  sampled {o:.6f}  z={(o-e)/se:+.2f}  ({exp_cyc[c]} shapes)")
    print("\ndiagonal near-miss count: exact freq vs sampled")
    for m in sorted(exp_nm):
        e = exp_nm[m] / len(shapes); o = obs_nm[m] / total
        se = (e * (1 - e) / total) ** 0.5
        print(f"  m={m}: exact {e:.6f}  sampled {o:.6f}  z={(o-e)/se:+.2f}  ({exp_nm[m]} shapes)")
