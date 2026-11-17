"""Find the loop-connected clumps of a saved shape (a .bin snapshot: int32 n, then n pairs of
int32 x, y) and return their cells.  Used to draw what a large clump looks like."""
import sys
import numpy as np


def read(fn):
    a = np.fromfile(fn, dtype="<i4"); n = a[0]; return a[1:1 + 2 * n].reshape(n, 2)


def clumps(xy):
    """list of arrays of cell indices, one per clump of at least 4 cells, largest first"""
    n = len(xy); idx = {(int(x), int(y)): i for i, (x, y) in enumerate(xy)}
    nbr = [[idx[(x + dx, y + dy)] for dx, dy in ((1, 0), (0, 1), (-1, 0), (0, -1)) if (x + dx, y + dy) in idx] for x, y in map(tuple, xy.tolist())]
    disc = [0] * n; low = [0] * n; par = [-1] * n; pit = [0] * n; stk = [0]; cstk = [0]; t = 1; disc[0] = low[0] = 1; out = []
    while stk:
        v = stk[-1]
        if pit[v] < len(nbr[v]):
            w = nbr[v][pit[v]]; pit[v] += 1
            if not disc[w]:
                t += 1; par[w] = v; disc[w] = low[w] = t; stk.append(w); cstk.append(w)
            elif w != par[v]:
                low[v] = min(low[v], disc[w])
        else:
            stk.pop(); p = par[v]
            if p >= 0: low[p] = min(low[p], low[v])
            if p < 0 or low[v] > disc[p]:
                c = []
                while True:
                    u = cstk.pop(); c.append(u)
                    if u == v: break
                if len(c) >= 4: out.append(np.array(c))
    return sorted(out, key=len, reverse=True)


if __name__ == "__main__":
    xy = read(sys.argv[1]); cl = clumps(xy)
    print(len(xy), "cells;", len(cl), "clumps of 4 or more; largest:", [len(c) for c in cl[:10]])
