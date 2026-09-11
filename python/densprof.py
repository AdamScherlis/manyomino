"""Radial density profile: occupancy fraction vs distance from centroid,
averaged over equilibrated snapshots.  Compares the b=30000 tilted
n=100k ensemble against beta=0 at the same size.
Usage: python3 densprof.py [out.png]
"""
import glob
import math
import sys

from render import write_png
from plotpng import Canvas
import betafig  # extends the bitmap FONT with full lowercase
from betafig import Panel, INK, GRID


def profile(paths, dr=3.0):
    """Mean occupancy fraction per annulus, averaged over snapshots."""
    prof_sum = {}
    prof_cnt = {}
    for p in paths:
        cells = []
        with open(p) as f:
            for line in f:
                x, y = line.split()
                cells.append((int(x), int(y)))
        n = len(cells)
        cx = sum(c[0] for c in cells) / n
        cy = sum(c[1] for c in cells) / n
        counts = {}
        for x, y in cells:
            k = int(math.hypot(x - cx, y - cy) / dr)
            counts[k] = counts.get(k, 0) + 1
        kmax = max(counts)
        for k in range(kmax + 1):
            area = math.pi * ((k + 1) ** 2 - k ** 2) * dr * dr
            rho = counts.get(k, 0) / area
            prof_sum[k] = prof_sum.get(k, 0.0) + rho
            prof_cnt[k] = prof_cnt.get(k, 0) + 1
    # average only over snapshots that reach the bin
    return {k: prof_sum[k] / prof_cnt[k] for k in prof_sum}, dr


def main(out):
    tilt = sorted(
        glob.glob("../gallery/raw/tilt_b30000_*.txt"),
        key=lambda p: int(p.split("_")[-1][:-4]),
    )
    tilt = [p for p in tilt if int(p.split("_")[-1][:-4]) >= 60000000]
    ref = sorted(glob.glob("../gallery/raw/n100000_rect_*.txt"))[-3:]
    pt, drt = profile(tilt)
    pr, drr = profile(ref, dr=8.0)
    print(f"{len(tilt)} tilted snapshots, {len(ref)} beta=0 snapshots")

    W, H = 1100, 700
    cv = Canvas(W, H)
    p = Panel(cv, 100, 60, 900, 540, (0, 900), (0, 1.05), False, False)
    # frame + ticks (linear x, so draw manually)
    for r in range(0, 901, 150):
        x = p.X(r)
        cv.line(x, p.py, x, p.py + p.ph, GRID)
        cv.text(x - 12, p.py + p.ph + 10, str(r), INK)
    for v in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        p.ytick(v, f"{v:.1f}")
    cv.line(p.px, p.py + p.ph, p.px + p.pw, p.py + p.ph, INK)
    cv.line(p.px, p.py, p.px, p.py + p.ph, INK)

    COL_T = (210, 80, 15)
    COL_R = (59, 91, 240)
    p.curve([((k + 0.5) * drt, v) for k, v in sorted(pt.items())
             if (k + 0.5) * drt <= 900], COL_T, r=0)
    p.curve([((k + 0.5) * drr, v) for k, v in sorted(pr.items())
             if (k + 0.5) * drr <= 900], COL_R, r=0)
    cv.text(p.px + 220, p.py + 30, "b=30000 (collapsed)", COL_T)
    cv.text(p.px + 440, p.py + 300, "b=0 (branched polymer)", COL_R)
    cv.text(p.px + 8, p.py - 30,
            "occupied fraction vs distance from centroid, n=100000", INK, 2)
    cv.text(p.px + 330, H - 26, "r (lattice units)", INK, 2)
    write_png(out, W, H, cv.rows)
    print(out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "../results/densprof_b30000.png")
