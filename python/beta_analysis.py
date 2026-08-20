"""Analyze the beta sweep: response curves vs b = beta*<Rg2>_0, and entropy
via thermodynamic integration  d logZ/d beta = -<Rg^2>_beta.

S(beta) = logZ(beta) + beta*<Rg2>_beta;  reported as [S(beta)-S(0)]/n.
Usage: python3 beta_analysis.py results/beta n R2REF"""
import glob
import math
import os
import sys


def series(path, col, burn=0.4):
    rows = []
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.rstrip().split(",")
            try:
                rows.append(float(p[col]))
            except (ValueError, IndexError):
                pass
    xs = rows[int(len(rows) * burn):]
    if len(xs) < 32:
        return None
    n = len(xs)
    m = sum(xs) / n
    nb = 8
    b = n // nb
    ms = [sum(xs[i*b:(i+1)*b]) / b for i in range(nb)]
    mm = sum(ms) / nb
    se = math.sqrt(sum((x - mm) ** 2 for x in ms) / (nb - 1) / nb)
    sd = math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1))
    return m, se, sd


def run(d, n, r2ref):
    rows = []
    for f in sorted(glob.glob(os.path.join(d, f"n{n}_b*.csv"))):
        tag = os.path.basename(f)[len(f"n{n}_b"):-4]
        b = float(tag.replace("p", ".").replace("m", "-"))
        cols = {}
        for ci, name in ((1, "rg2"), (2, "perim"), (3, "cyc"), (4, "asph"), (5, "mov")):
            r = series(f, ci)
            if r:
                cols[name] = r
        if "rg2" in cols:
            rows.append((b, cols))
    rows.sort()
    print(f"{'b':>6} {'beta':>10} {'<Rg2>/R2ref':>11} {'perim/n':>8} {'cyc/n':>7} "
          f"{'movable':>8} {'asph':>6}")
    for b, c in rows:
        beta = b / r2ref
        print(f"{b:>6.2f} {beta:>10.2e} {c['rg2'][0]/r2ref:>11.4f} "
              f"{c['perim'][0]/n:>8.4f} {c['cyc'][0]/n:>7.4f} "
              f"{c.get('mov',(float('nan'),))[0]:>8.4f} {c.get('asph',(float('nan'),))[0]:>6.3f}")
    # thermodynamic integration over beta (sorted)
    bs = [b / r2ref for b, _ in rows]
    r2s = [c["rg2"][0] for _, c in rows]
    # logZ relative to beta=0: integrate -<R2> d beta from 0
    i0 = min(range(len(bs)), key=lambda i: abs(bs[i]))
    logz = {i0: 0.0}
    for i in range(i0 + 1, len(bs)):
        logz[i] = logz[i - 1] - 0.5 * (r2s[i] + r2s[i - 1]) * (bs[i] - bs[i - 1])
    for i in range(i0 - 1, -1, -1):
        logz[i] = logz[i + 1] + 0.5 * (r2s[i] + r2s[i + 1]) * (bs[i + 1] - bs[i])
    print(f"\n{'b':>6} {'dlogZ/n':>10} {'dS/n':>10}")
    for i, (b, c) in enumerate(rows):
        beta = b / r2ref
        ds = (logz[i] + beta * r2s[i] - (logz[i0] + bs[i0] * r2s[i0])) / n
        print(f"{b:>6.2f} {logz[i]/n:>10.5f} {ds:>10.5f}")


if __name__ == "__main__":
    run(sys.argv[1], int(sys.argv[2]), float(sys.argv[3]))
