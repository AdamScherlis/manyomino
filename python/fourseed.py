"""Four-chain two-seed test: pool chains by init lineage (stringy vs
compact), compare pooled means, and also report the chain-to-chain scatter.

Usage: python3 fourseed.py <stringy1.csv,stringy2.csv> <compact1.csv,compact2.csv> [burn]"""
import math
import sys

from analyze import analyze


def pool(paths, burn):
    rs = []
    for p in paths:
        try:
            rs.append(analyze(p, burn))
        except Exception:
            pass
    if not rs:
        return None
    w = [1 / r["rg2_stderr"] ** 2 for r in rs]
    m = sum(r["rg2_mean"] * wi for r, wi in zip(rs, w)) / sum(w)
    return m, math.sqrt(1 / sum(w)), rs


def run(stringy, compact, burn=0.5):
    a = pool(stringy, burn)
    b = pool(compact, burn)
    if a is None or b is None:
        print("insufficient data")
        return False
    for name, (m, se, rs) in (("stringy", a), ("compact", b)):
        per = "  ".join(f"{r['rg2_mean']:.4g}±{r['rg2_stderr']:.2g}(ESS {r['ess']:.0f})" for r in rs)
        print(f"{name:>8}: pooled {m:.6g} ± {se:.3g}   [{per}]")
    z = abs(a[0] - b[0]) / math.sqrt(a[1] ** 2 + b[1] ** 2)
    ess_ok = all(r["ess"] >= 8 for r in a[2] + b[2])
    print(f"|z| = {z:.2f}   {'PASS' if z < 3 and ess_ok else 'FAIL'}")
    return z < 3 and ess_ok


if __name__ == "__main__":
    stringy = sys.argv[1].split(",")
    compact = sys.argv[2].split(",")
    burn = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
    sys.exit(0 if run(stringy, compact, burn) else 1)
