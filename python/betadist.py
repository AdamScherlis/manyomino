"""beta=0 distributions of the tilted-study observables vs n.

For each results/beta/n<N>_b0.csv: mean, stddev, rel. sigma, skewness of
Rg2, perimeter/cell, cycles/cell, asphericity, movable fraction; stderr
of the mean from 8-batch means on the post-burn series.
Usage: python3 betadist.py [results/beta]
"""
import glob
import math
import os
import re
import sys


def moments(xs):
    n = len(xs)
    m = sum(xs) / n
    v = sum((x - m) ** 2 for x in xs) / (n - 1)
    sd = math.sqrt(v)
    sk = (sum((x - m) ** 3 for x in xs) / n) / sd ** 3 if sd > 0 else 0.0
    nb = 8
    b = n // nb
    ms = [sum(xs[i * b:(i + 1) * b]) / b for i in range(nb)]
    mm = sum(ms) / nb
    se = math.sqrt(sum((x - mm) ** 2 for x in ms) / (nb - 1) / nb)
    return m, se, sd, sk


def main(d):
    groups = {}
    for f in glob.glob(os.path.join(d, "n*_b0*.csv")):
        m = re.match(r"n(\d+)_b0(_g\d+)?\.csv", os.path.basename(f))
        if m:
            groups.setdefault(int(m.group(1)), []).append(f)
    files = sorted(groups.items())
    print(f"{'n':>6} {'obs':>10} {'mean':>12} {'stderr':>10} {'stddev':>10} "
          f"{'rel sd':>7} {'skew':>6}")
    for n, fs in files:
        rows = set()
        for f in fs:
            with open(f) as fh:
                fh.readline()
                for line in fh:
                    p = line.rstrip().split(",")
                    if len(p) < 6:
                        continue
                    try:
                        rows.add(tuple(float(x) for x in p))
                    except ValueError:
                        continue
        cols = [[] for _ in range(6)]
        for r in sorted(rows):
            for i, v in enumerate(r):
                cols[i].append(v)
        burn = int(len(cols[0]) * 0.2)
        for ci, name, scale in ((1, "Rg2", 1.0), (2, "perim/n", 1.0 / n),
                                (3, "cyc/n", 1.0 / n), (4, "asph", 1.0),
                                (5, "movfrac", 1.0)):
            xs = [v * scale for v in cols[ci][burn:]]
            if len(xs) < 64:
                continue
            m, se, sd, sk = moments(xs)
            print(f"{n:>6} {name:>10} {m:>12.5g} {se:>10.2g} {sd:>10.3g} "
                  f"{sd/m:>7.3f} {sk:>6.2f}")
        print()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "../results/beta")
