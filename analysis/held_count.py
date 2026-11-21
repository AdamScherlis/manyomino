"""Count the polyominoes of n cells that contain (0,0) and (1,1) and neither (1,0) nor (0,1):
the shapes `outward_native held n ...` can visit, each of which it should visit equally often.

    python3 analysis/held_count.py 9 10 12 13
    python3 analysis/held_count.py cut 6 8 10      # for `outward_native cut`: shapes that contain (0,0) and (1,0)
                                                   # and are in one piece even if those two do not count as joined

Every fixed polyomino is listed once (Redelmeier's method) and counted once for each of its
cells c with c + (1,1) in it and c + (1,0), c + (0,1) not.

Below 12 cells the sampler reaches exactly half of them: the two held cells are joined by a path
round one side of the box or the other, and changing sides needs a ring of 12 cells round both.
"""
import sys
sys.setrecursionlimit(10000)


def count(n, cutmode=False):
    total = 0
    cells = set()
    seen = set()

    def ok(c):
        return c[1] > 0 or (c[1] == 0 and c[0] >= 0)

    def joined_without(a, b):
        """is the shape still in one piece if cells a and b, which touch, do not count as joined?"""
        seen_ = {a}; st = [a]
        while st:
            c = st.pop()
            for d in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                e = (c[0] + d[0], c[1] + d[1])
                if e in cells and e not in seen_ and {c, e} != {a, b}: seen_.add(e); st.append(e)
        return len(seen_) == len(cells)

    def tally():
        t = 0
        if cutmode: return sum(1 for (x, y) in cells if (x + 1, y) in cells and joined_without((x, y), (x + 1, y)))
        for (x, y) in cells:
            if (x + 1, y + 1) in cells and (x + 1, y) not in cells and (x, y + 1) not in cells:
                t += 1
        return t

    def grow(untried):
        nonlocal total
        untried = list(untried)
        while untried:
            c = untried.pop()
            cells.add(c)
            if len(cells) == n:
                total += tally()
            else:
                new = []
                for d in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    e = (c[0] + d[0], c[1] + d[1])
                    if ok(e) and e not in seen:
                        seen.add(e)
                        new.append(e)
                grow(untried + new)
                for e in new:
                    seen.discard(e)
            cells.discard(c)

    seen.add((0, 0))
    grow([(0, 0)])
    return total


if __name__ == "__main__":
    cutmode = len(sys.argv) > 1 and sys.argv[1] == "cut"
    for a in sys.argv[2 if cutmode else 1:]:
        print(a, count(int(a), cutmode), flush=True)
