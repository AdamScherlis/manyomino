"""Two-seed convergence test: compare equilibrium R_g^2 statistics of a
bar-seeded and a rect-seeded chain at the same n.

PASS iff the autocorrelation-corrected means agree within 3 combined
standard errors (and each chain has a sane ESS).

Usage: python3 twoseed.py bar.csv rect.csv [burn_frac]
"""

import math
import sys

from analyze import analyze


def run(bar_path, rect_path, burn_frac=0.5):
    a = analyze(bar_path, burn_frac)
    b = analyze(rect_path, burn_frac)
    dz = abs(a["rg2_mean"] - b["rg2_mean"]) / math.sqrt(
        a["rg2_stderr"] ** 2 + b["rg2_stderr"] ** 2
    )
    for r, name in ((a, "bar"), (b, "rect")):
        print(
            f"{name:>4}: <Rg^2> = {r['rg2_mean']:9.1f} +- {r['rg2_stderr']:7.1f}  "
            f"tau = {r['tau_steps']:.2e} steps  ESS = {r['ess']:.0f}"
        )
    print(f"|z| = {dz:.2f}")
    ok = dz < 3.0 and a["ess"] >= 10 and b["ess"] >= 10
    print("PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    burn = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
    sys.exit(0 if run(sys.argv[1], sys.argv[2], burn) else 1)
