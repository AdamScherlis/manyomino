"""Verify detailed balance for a recorded single-cell transition.

Input: the post-move dump and the move (c removed, s added), as written by
`manyomino cyclewatch`.  Reconstructs the predecessor A = A' - {s} + {c} and
checks, entirely from scratch (no sampler code):
  - both states are valid n-cell polyominoes (connected, right size);
  - A - c = A' - s (the shared intermediate);
  - c is removable in A, s is removable in A' (both directions legal);
  - s is in the perimeter of A - c, c is in the perimeter of A' - s;
  - T(A -> A') = 1/(n |P(A - c)|) equals T(A' -> A) = 1/(n |P(A' - s)|),
    with the two perimeters computed independently from each state.

Usage: python3 verify_db.py after.txt meta.txt
"""
import sys
from collections import deque


def load(path):
    return set(tuple(map(int, l.split())) for l in open(path))


def meta(path):
    m = {}
    for line in open(path):
        parts = line.split()
        if parts[0] in ("c", "s"):
            m[parts[0]] = (int(parts[1]), int(parts[2]))
        else:
            m[parts[0]] = " ".join(parts[1:])
    return m


def connected(cs):
    it = iter(cs)
    start = next(it)
    seen = {start}
    dq = deque([start])
    while dq:
        x, y = dq.popleft()
        for nb in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
            if nb in cs and nb not in seen:
                seen.add(nb)
                dq.append(nb)
    return len(seen) == len(cs)


def perimeter(cs):
    per = set()
    for (x, y) in cs:
        for nb in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
            if nb not in cs:
                per.add(nb)
    return per


def run(after_path, meta_path):
    A2 = load(after_path)
    m = meta(meta_path)
    c, s = m["c"], m["s"]
    n = len(A2)
    assert s in A2 and c not in A2, "move endpoints inconsistent with dump"
    A1 = (A2 - {s}) | {c}
    assert len(A1) == len(A2) == n
    assert connected(A1), "predecessor disconnected"
    assert connected(A2), "successor disconnected"
    mid_f = A1 - {c}
    mid_r = A2 - {s}
    assert mid_f == mid_r, "A - c != A' - s"
    assert connected(mid_f), "intermediate disconnected (c not removable)"
    Pf = perimeter(mid_f)
    Pr = perimeter(mid_r)
    assert s in Pf, "s not selectable forward"
    assert c in Pr, "c not selectable backward"
    assert Pf == Pr and len(Pf) == len(Pr)
    pf = 1.0 / (n * len(Pf))
    pr = 1.0 / (n * len(Pr))
    print(f"{after_path}:")
    print(f"  n = {n}; move c={c} -> s={s}   ({m.get('cause', 'birth')})")
    print(f"  |P(A-c)| = {len(Pf)}  (computed from predecessor)")
    print(f"  |P(A'-s)| = {len(Pr)}  (computed from successor)")
    print(f"  T(A->A') = 1/(n*|P|) = {pf:.6e}")
    print(f"  T(A'->A) = 1/(n*|P|) = {pr:.6e}")
    assert pf == pr
    print("  DETAILED BALANCE HOLDS (exact equality)")
    return A1, A2, c, s


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2])
