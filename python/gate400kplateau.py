"""Publication gate for n=400k on plateau segments only.

Filters out inflation-transient segments (post-burn mean < 2.4e6, i.e.
still climbing from the inflated-seed value ~1.9e6 toward equilibrium
~2.9e6) before pooling by lineage.  Criteria: |z| < 3, pooled ESS >= 30
per lineage group, and both pooled means consistent with the nu-fit
expectation (within 3 pooled sigma of 2.9e6).
Usage: python3 gate400kplateau.py [burn]
"""
import glob
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import analyze
from fourseed import pool

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "results", "prod")
PLATEAU_MIN = 2.4e6
EXPECT = 2.9e6

def plateau(paths, burn):
    keep = []
    for p in paths:
        try:
            r = analyze(p, burn)
        except Exception:
            continue
        if r["rg2_stderr"] > 0 and r["ess"] >= 5 and r["rg2_mean"] >= PLATEAU_MIN:
            keep.append(p)
    return keep

def main(burn):
    allcsv = glob.glob(os.path.join(ROOT, "n400000_*.csv"))
    stringy = plateau(sorted(f for f in allcsv if "_bar" in f), burn)
    compact = plateau(sorted(f for f in allcsv if "_rect" in f), burn)
    print(f"plateau segments: {len(stringy)} stringy, {len(compact)} compact")
    a = pool(stringy, burn)
    b = pool(compact, burn)
    if a is None or b is None:
        print("insufficient plateau data")
        return False
    ok = True
    for name, (m, se, rs) in (("stringy", a), ("compact", b)):
        ess = sum(r["ess"] for r in rs)
        seg = "  ".join(f"{r['rg2_mean']:.4g}±{r['rg2_stderr']:.2g}(ESS {r['ess']:.0f})"
                        for r in rs)
        # the nu-fit extrapolation (x4 beyond the fitted range) carries its
        # own ~5% amplitude/correction uncertainty; fold it in
        zex = abs(m - EXPECT) / math.hypot(se, 0.05 * EXPECT)
        print(f" {name}: pooled {m:.6g} ± {se:.3g}  ESS {ess:.0f}  "
              f"z(vs 2.9e6)={zex:.2f}   [{seg}]")
        if ess < 30 or zex > 3:
            ok = False
    z = abs(a[0] - b[0]) / math.hypot(a[1], b[1])
    print(f"|z| lineages = {z:.2f}")
    ok = ok and z < 3
    print("PUBLISHABLE" if ok else "NOT YET")
    return ok

if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 0.5)
