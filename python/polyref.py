"""Slow, obviously-correct reference implementation of the fixed-n single-cell chain.

State: a set A of n cells in Z^2, 4-connected. One move:
  1. Pick a cell c uniformly from A.
  2. If A - {c} is disconnected: stay put (still counts as a step).
  3. Else pick s uniformly from the perimeter sites of A - {c} (empty sites
     4-adjacent to A - {c}; c itself is always one of them, so s = c is legal).
  4. New state: (A - {c}) | {s}.

Everything here is recomputed from scratch each move so there is nothing
incremental to get wrong.  This is the ground truth the fast implementation
is differential-tested against.
"""

import random


def neighbors(c):
    x, y = c
    return ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1))


def is_connected(cells):
    if not cells:
        return True
    it = iter(cells)
    start = next(it)
    seen = {start}
    stack = [start]
    while stack:
        c = stack.pop()
        for nb in neighbors(c):
            if nb in cells and nb not in seen:
                seen.add(nb)
                stack.append(nb)
    return len(seen) == len(cells)


def perimeter(cells):
    """Empty sites 4-adjacent to at least one cell (a set of *sites*, not pairs)."""
    per = set()
    for c in cells:
        for nb in neighbors(c):
            if nb not in cells:
                per.add(nb)
    return per


def canonical(cells):
    """Translation-canonical form: shift min corner to origin, sorted tuple."""
    mx = min(x for x, y in cells)
    my = min(y for x, y in cells)
    return tuple(sorted((x - mx, y - my) for x, y in cells))


def bar(n):
    return {(i, 0) for i in range(n)}


def rect(n):
    """Roughly sqrt(n) x sqrt(n) rectangle: full rows plus a partial last row."""
    w = max(1, round(n ** 0.5))
    return {(i % w, i // w) for i in range(n)}


def rg2(cells):
    n = len(cells)
    sx = sum(x for x, y in cells)
    sy = sum(y for x, y in cells)
    sx2 = sum(x * x for x, y in cells)
    sy2 = sum(y * y for x, y in cells)
    return (sx2 + sy2) / n - (sx / n) ** 2 - (sy / n) ** 2


class ReferenceChain:
    def __init__(self, n, seed=0, init="bar"):
        assert n >= 2
        self.n = n
        self.rng = random.Random(seed)
        self.cells = bar(n) if init == "bar" else rect(n)
        assert len(self.cells) == n and is_connected(self.cells)

    def step(self):
        """One move of the chain. Returns True if the state changed."""
        c = self.rng.choice(sorted(self.cells))
        rest = self.cells - {c}
        if not is_connected(rest):
            return False
        s = self.rng.choice(sorted(perimeter(rest)))
        self.cells = rest | {s}
        return s != c

    def check_invariants(self):
        assert len(self.cells) == self.n
        assert is_connected(self.cells)
