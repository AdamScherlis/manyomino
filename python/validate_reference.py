"""Validation of the reference chain: chi-square uniformity at small n.

Runs the reference chain, records the translation-canonical shape every
`thin` moves, and chi-squares the counts against uniform over the exact
A001168 shape list from the enumerator.  Catches proposal-distribution bugs
(perimeter pair-vs-site weighting, missing s=c no-op, wrong perimeter set).

Usage: python3 validate_reference.py [n] [observations] [thin] [seed]
"""

import sys
import time

from enumerate_fixed import enumerate_fixed
from polyref import ReferenceChain, canonical
from stats import chisquare_uniform


def run(n=5, obs=200_000, thin=10, seed=1, burn=10_000):
    shapes = sorted(enumerate_fixed(n)[n])
    index = {s: i for i, s in enumerate(shapes)}
    counts = [0] * len(shapes)

    chain = ReferenceChain(n, seed=seed, init="bar")
    for _ in range(burn):
        chain.step()
    t0 = time.time()
    for i in range(obs):
        for _ in range(thin):
            chain.step()
        counts[index[canonical(chain.cells)]] += 1
        if i % 500 == 0:
            chain.check_invariants()
    dt = time.time() - t0

    assert min(counts) > 0, "some shape never observed"
    chi2, dof, p = chisquare_uniform(counts, obs)
    print(
        f"n={n} shapes={len(shapes)} obs={obs} thin={thin} seed={seed} "
        f"({obs * thin / dt:.0f} moves/s)"
    )
    print(f"chi2={chi2:.1f} dof={dof} p={p:.4f}")
    return p


if __name__ == "__main__":
    args = sys.argv[1:]
    n = int(args[0]) if len(args) > 0 else 5
    obs = int(args[1]) if len(args) > 1 else 200_000
    thin = int(args[2]) if len(args) > 2 else 10
    seed = int(args[3]) if len(args) > 3 else 1
    p = run(n, obs, thin, seed)
    ok = p > 0.001
    print("PASS" if ok else "FAIL (p < 0.001)")
    sys.exit(0 if ok else 1)
