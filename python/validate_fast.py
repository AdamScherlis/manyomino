"""Chi-square the Rust sampler's shape counts against uniform over the exact
A001168 shape list, and cross-check that every enumerated shape was seen.

Reads the "canonical_u64 count" lines produced by `manyomino chi`.

Usage: python3 validate_fast.py counts.txt n
"""

import sys

from enumerate_fixed import enumerate_fixed
from stats import chisquare_uniform


def pack(shape):
    """Match State::canonical_u64: byte (x<<4)|y per cell, sorted, folded."""
    bs = sorted((x << 4) | y for x, y in shape)
    k = 0
    for b in bs:
        k = (k << 8) | b
    return k


def run(path, n):
    shapes = enumerate_fixed(n)[n]
    keys = {pack(s) for s in shapes}
    counts = {}
    with open(path) as f:
        for line in f:
            k, v = line.split()
            counts[int(k)] = int(v)
    unknown = set(counts) - keys
    assert not unknown, f"sampler produced {len(unknown)} shapes not in the exact list"
    missing = keys - set(counts)
    assert not missing, f"{len(missing)} shapes never observed"
    obs = [counts[k] for k in sorted(keys)]
    total = sum(obs)
    chi2, dof, p = chisquare_uniform(obs, total)
    print(f"{path}: n={n} shapes={len(keys)} obs={total}")
    print(f"chi2={chi2:.1f} dof={dof} p={p:.4f}")
    return p


if __name__ == "__main__":
    p = run(sys.argv[1], int(sys.argv[2]))
    ok = p > 0.001
    print("PASS" if ok else "FAIL (p < 0.001)")
    sys.exit(0 if ok else 1)
