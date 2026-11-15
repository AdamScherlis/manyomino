"""Machine check, on small sizes, of every step of the proof in proof/README.md.

For all ordered pairs (A, B) of fixed polyominoes with n cells it builds the canonical path,
and checks that
  1. every step is a legal move of the relocation chain and every intermediate shape is connected;
  2. the cells already removed from A (and those of B not yet placed) always form at most
     2*log2(n) + 4 connected pieces, each touching the rest;
  3. the pair (A, B) can be recovered from the intermediate shape, the step number, the dividing
     column and those pieces (the map used to count paths through a transition is injective);
  4. the number of paths through any transition is below the bound of Lemma 3;
and it compares the resulting bound on the relaxation time with the exact one, obtained from
the eigenvalues of the chain's transition matrix.

  python3 proof/check_paths.py 7        # sizes 2..7; size 8 takes about an hour
"""
import math, sys, itertools
from collections import defaultdict, deque
import numpy as np

DIRS = ((1, 0), (0, 1), (-1, 0), (0, -1))


def canon(cells):
    mx = min(x for x, y in cells); my = min(y for x, y in cells)
    return tuple(sorted((x - mx, y - my) for x, y in cells))


def enumerate_polyominoes(n):
    cur = {((0, 0),)}
    for _ in range(n - 1):
        nxt = set()
        for p in cur:
            s = set(p)
            for (x, y) in p:
                for dx, dy in DIRS:
                    q = (x + dx, y + dy)
                    if q not in s: nxt.add(canon(p + (q,)))
        cur = nxt
    return sorted(cur)


def connected(cells):
    cells = set(cells)
    if not cells: return True
    start = next(iter(cells)); seen = {start}; dq = deque([start])
    while dq:
        x, y = dq.popleft()
        for dx, dy in DIRS:
            q = (x + dx, y + dy)
            if q in cells and q not in seen: seen.add(q); dq.append(q)
    return len(seen) == len(cells)


def order(cells, root):
    """The canonical order of Definition 2: a breadth-first spanning tree from the root, then a
    depth-first preorder visiting the children of each cell in increasing order of subtree size
    (heaviest last).  Returns (order, parent, subtree size)."""
    s = set(cells); parent = {root: None}; bfs = [root]; dq = deque([root])
    while dq:
        v = dq.popleft()
        for dx, dy in DIRS:
            q = (v[0] + dx, v[1] + dy)
            if q in s and q not in parent: parent[q] = v; bfs.append(q); dq.append(q)
    size = {v: 1 for v in bfs}; kids = defaultdict(list)
    for v in reversed(bfs):
        if parent[v] is not None: size[parent[v]] += size[v]; kids[parent[v]].append(v)
    out = []; stack = [root]
    while stack:
        v = stack.pop(); out.append(v)
        ch = sorted(kids[v], key=lambda c: (size[c], DIRS.index((c[0] - v[0], c[1] - v[1]))))
        stack.extend(reversed(ch))                       # lightest child is visited first
    return out, parent, size


def pieces(cells):
    """connected components of a set of cells"""
    cells = set(cells); out = []
    while cells:
        start = cells.pop(); comp = {start}; dq = deque([start])
        while dq:
            x, y = dq.popleft()
            for dx, dy in DIRS:
                q = (x + dx, y + dy)
                if q in cells: cells.remove(q); comp.add(q); dq.append(q)
        out.append(comp)
    return out


def adjacent(a, b):
    b = set(b)
    return any((x + dx, y + dy) in b for (x, y) in a for dx, dy in DIRS)


def main(nmax):
    for n in range(2, nmax + 1):
        P = enumerate_polyominoes(n); N = len(P); idx = {p: i for i, p in enumerate(P)}
        L = 2 * math.log2(n) + 4
        # canonical orders: A is taken apart from its rightmost (then topmost) cell, B is built from its leftmost (then lowest) cell
        ordA = [order(p, max(p)) [0] for p in P]
        ordB = [order(p, min(p)) [0] for p in P]
        # 2. piece counts of every suffix
        worst = 0
        for p, oa, ob in zip(P, ordA, ordB):
            for o in (oa, ob):
                for m in range(1, n):
                    pc = pieces(o[m:]); worst = max(worst, len(pc))
                    assert connected(o[:m]) and all(adjacent(c, o[:m]) for c in pc)
        assert worst <= L, (worst, L)
        # the chain: remove a cell that is not a cut cell, add one at an empty site next to the rest; each (cell, site) has probability 1/(4 n^2)
        T = np.zeros((N, N))
        for i, p in enumerate(P):
            s = set(p)
            for c in p:
                rest = s - {c}
                if not connected(rest): continue
                sites = {(x + dx, y + dy) for (x, y) in rest for dx, dy in DIRS} - rest
                for q in sites:
                    T[i, idx[canon(rest | {q})]] += 1.0 / (4 * n * n)
        off = T.sum(1) - np.diag(T); assert np.allclose(T, T.T)
        M = T - np.diag(np.diag(T)) + np.diag(1 - off)
        ev = np.sort(np.linalg.eigvalsh(M)); gap = 1 - ev[-2] if N > 1 else 1.0
        # 1, 3, 4. all canonical paths
        load = defaultdict(int); keys = set(); steps = 0; maxlen = 0
        for ia, A in enumerate(P):
            oa = ordA[ia]; ra = oa[0]
            for ib, B in enumerate(P):
                if ia == ib: continue
                ob0 = ordB[ib]; sh = (ra[0] + 1 - ob0[0][0], ra[1] - ob0[0][1])       # put B's root just right of A's root
                ob = [(x + sh[0], y + sh[1]) for x, y in ob0]
                Z = set(oa); assert not (Z & set(ob))
                zi = ia
                for k in range(1, n + 1):
                    c = oa[n - k]; q = ob[k - 1]
                    rest = Z - {c}
                    assert c in Z and q not in Z and connected(rest) and adjacent([q], rest)      # a legal relocation
                    Z2 = rest | {q}; z2 = idx[canon(Z2)]
                    if z2 != zi:
                        load[(zi, z2)] += 1; steps += 1
                        # the pair is determined by: the transition, the step, the dividing column, the removed part of A and the missing part of B,
                        # all written relative to the shape Z (so that they do not depend on where A was placed)
                        mx = min(x for x, y in Z); my = min(y for x, y in Z)
                        rel = lambda S: tuple(sorted((x - mx, y - my) for x, y in S))
                        key = (zi, z2, k, ra[0] - mx, rel(oa[n - k + 1:]), rel(ob[k - 1:]))
                        assert key not in keys; keys.add(key)
                    Z = Z2; zi = z2
                assert zi == ib
        Nmax = max(load.values())
        rho = max(cnt * n / (N * T[e[0], e[1]]) for e, cnt in load.items())       # Sinclair's congestion, with every path counted at full length n
        lemma3 = n * n * (L + 1) ** 2 * (2 * n + 2) ** (6 * L) * N
        thm = 4 * n ** 5 * (L + 1) ** 2 * (2 * n + 2) ** (6 * L)
        assert Nmax <= lemma3 and 1 / gap <= rho * (1 + 1e-9) <= thm
        print(f"n={n}: {N} shapes, {N * (N - 1)} paths, all steps legal; pieces <= {worst} (bound {L:.1f}); "
              f"most paths through one transition {Nmax} = {Nmax / N:.2f} x number of shapes; "
              f"exact relaxation time {1 / gap:.1f} moves <= path bound {rho:.0f} <= theorem {thm:.3g}", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
