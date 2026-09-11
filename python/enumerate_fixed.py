"""Brute-force enumerator of fixed polyominoes (translation classes), OEIS A001168.

Grows size-n shapes from size-(n-1) shapes by adding perimeter sites and
deduplicating translation-canonical forms.  Independent of the sampler code
paths (only shares the trivial neighbor/perimeter/canonical helpers).
"""

import sys

from polyref import canonical, perimeter

A001168 = [None, 1, 2, 6, 19, 63, 216, 760, 2725, 9910, 36446]


def enumerate_fixed(nmax):
    """Return {n: set of canonical shapes} for n = 1..nmax."""
    shapes = {1: {((0, 0),)}}
    for n in range(2, nmax + 1):
        grown = set()
        for sh in shapes[n - 1]:
            cs = set(sh)
            for p in perimeter(cs):
                grown.add(canonical(cs | {p}))
        shapes[n] = grown
    return shapes


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    shapes = enumerate_fixed(nmax)
    ok = True
    for n in range(1, nmax + 1):
        got = len(shapes[n])
        want = A001168[n] if n < len(A001168) else "?"
        status = "OK" if got == want else "MISMATCH"
        ok = ok and got == want
        print(f"n={n:2d}  count={got:6d}  A001168={want}  {status}")
    sys.exit(0 if ok else 1)
