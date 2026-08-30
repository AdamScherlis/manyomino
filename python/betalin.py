"""Linear-in-b response figure for the tilted ensemble, wide range with
the deep-negative branch, plus ABSOLUTE entropy per cell.

Four panels vs b = beta*<Rg2>_0 on a LINEAR axis:
 1. <Rg2>/<Rg2>_0            (log y: spans ~50x)
 2. S/n absolute (nats/cell, linear y).  Anchor: S(0)/n from the animal
    growth constant, ln(lambda) - ln(n)/n  (lambda = 4.0626, theta = 1);
    the anchor carries ~0.005 nats/cell systematic.  Tilted values by
    thermodynamic integration d logZ/d beta = -<Rg2>_beta over the sweep
    grid.
 3. local observables (perim/cell, cycles/cell, movable) linear y
 4. asphericity, linear y
Usage: python3 betalin.py [out.png] [BMIN BMAX]
"""
import math
import sys

from render import write_png
from plotpng import Canvas
from betafig import Panel, INK, GRID, COLS, sweep, entropy
import os

LAM = 4.0626


def s0_per_cell(n):
    return math.log(LAM) - math.log(n) / n


class LinPanel(Panel):
    def frame_lin(self, title, xticks):
        cv = self.cv
        for r in xticks:
            x = self.X(r)
            cv.line(x, self.py, x, self.py + self.ph, GRID)
            lab = str(int(r))
            cv.text(x - 6 * len(lab), self.py + self.ph + 10, lab, INK)
        cv.line(self.px, self.py + self.ph, self.px + self.pw,
                self.py + self.ph, INK)
        cv.line(self.px, self.py, self.px, self.py + self.ph, INK)
        cv.text(self.px + 8, self.py - 22, title, INK, 2)


def main(out, bmin=-32.0, bmax=32.0):
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                     "results", "beta")
    data = {n: sweep(d, n) for n in (1000, 3000)}
    r2ref = {1000: 1379.0, 3000: 5385.0}
    W, H = 1250, 980
    cv = Canvas(W, H)
    PW, PH = 480, 330
    X0, Y0, GX, GY = 95, 70, 130, 130
    xticks = [x for x in (-30, -20, -10, 0, 10, 20, 30)
              if bmin <= x <= bmax]

    def inrange(rows, key):
        return [(b, c[key]) for b, c in rows if bmin <= b <= bmax]

    # panel 1: Rg2 ratio, log y
    p1 = LinPanel(cv, X0, Y0, PW, PH, (bmin, bmax), (0.08, 80), False, True)
    p1.frame_lin("<Rg2>/<Rg2>(0)  (log y)", xticks)
    for v, lab in ((0.1, "0.1"), (1.0, "1"), (10.0, "10")):
        p1.ytick(v, lab)
    for n, rows in data.items():
        r20 = [c["rg2"] for b, c in rows if b == 0][0]
        p1.curve([(b, v / r20) for b, v in inrange(rows, "rg2")], COLS[n])

    # panel 2: absolute entropy per cell, linear
    p2 = LinPanel(cv, X0 + PW + GX, Y0, PW, PH, (bmin, bmax), (0.0, 1.5),
                  False, False)
    p2.frame_lin("S/n  (nats/cell, absolute)", xticks)
    for v in (0.0, 0.5, 1.0, 1.4):
        p2.ytick(v, f"{v:.1f}")
    for n, rows in data.items():
        ent = entropy(rows, r2ref[n], n)  # [(b, -dS/n)]
        s0 = s0_per_cell(n)
        p2.curve([(b, s0 - v) for b, v in ent if bmin <= b <= bmax], COLS[n])
    cv.text(p2.px + 12, p2.py + p2.ph - 40,
            "anchor S(0)/n = ln(4.0626) - ln(n)/n", INK)

    # panel 3: locals, linear
    p3 = LinPanel(cv, X0, Y0 + PH + GY, PW, PH, (bmin, bmax), (0.0, 2.1),
                  False, False)
    p3.frame_lin("locals (per cell)", xticks)
    for v in (0.0, 0.5, 1.0, 1.5, 2.0):
        p3.ytick(v, f"{v:.1f}")
    styles = {"perim": (59, 91, 240), "cyc": (210, 80, 15),
              "mov": (90, 160, 90)}
    for n, rows in data.items():
        r = 4 if n == 3000 else 2
        for k, col in styles.items():
            pts = inrange(rows, k)
            if k in ("perim", "cyc"):
                pts = [(b, v / n) for b, v in pts]
            if k == "cyc":
                pts = [(b, v * 10) for b, v in pts]  # x10 for visibility
            p3.curve(pts, col, r=r)
    cv.text(p3.px + 14, p3.py + 12, "perim/cell", styles["perim"])
    cv.text(p3.px + 14, p3.py + 32, "cycles/cell x10", styles["cyc"])
    cv.text(p3.px + 14, p3.py + 52, "movable frac", styles["mov"])

    # panel 4: asphericity, linear
    p4 = LinPanel(cv, X0 + PW + GX, Y0 + PH + GY, PW, PH, (bmin, bmax),
                  (0.0, 1.05), False, False)
    p4.frame_lin("asphericity", xticks)
    for v in (0.0, 0.25, 0.5, 0.75, 1.0):
        p4.ytick(v, f"{v:.2f}")
    for n, rows in data.items():
        p4.curve(inrange(rows, "asph"), COLS[n])

    for i, (n, col) in enumerate(sorted(COLS.items())):
        cv.text(X0 + PW - 150, Y0 + 14 + 20 * i, f"n={n}", col)
    cv.text(X0 + 150, H - 34,
            "b = beta <Rg2>(0)   (negative b rewards spread)", INK, 2)
    write_png(out, W, H, cv.rows)
    print(out)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "../results/betalin.png"
    if len(sys.argv) > 3:
        main(out, float(sys.argv[2]), float(sys.argv[3]))
    else:
        main(out)
