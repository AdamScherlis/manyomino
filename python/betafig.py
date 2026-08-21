"""Response-curve figure for the tilted ensemble pi ~ exp(-beta*Rg^2).

Four panels vs the dimensionless coupling b = beta*<Rg^2>_0 (positive
branch, log axis): Rg^2 squeeze, entropy cost -dS/n, local observables
relative to b=0, asphericity.  Data: results/beta sweep CSVs at n=1000
and n=3000.  Stdlib-only PNG via the plotpng canvas.
Usage: python3 betafig.py [out.png]
"""
import math
import sys

from render import write_png
from plotpng import Canvas, FONT
from beta_analysis import series
import glob
import os

FONT.update({
    "a": [0, 0, 0b01110, 0b00001, 0b01111, 0b10001, 0b01111],
    "b": [0b10000, 0b10000, 0b10110, 0b11001, 0b10001, 0b11001, 0b10110],
    "d": [0b00001, 0b00001, 0b01101, 0b10011, 0b10001, 0b10011, 0b01101],
    "h": [0b10000, 0b10000, 0b10110, 0b11001, 0b10001, 0b10001, 0b10001],
    "i": [0b00100, 0, 0b01100, 0b00100, 0b00100, 0b00100, 0b01110],
    "m": [0, 0, 0b11010, 0b10101, 0b10101, 0b10101, 0b10101],
    "o": [0, 0, 0b01110, 0b10001, 0b10001, 0b10001, 0b01110],
    "p": [0, 0, 0b10110, 0b11001, 0b11001, 0b10110, 0b10000],
    "r": [0, 0, 0b10110, 0b11001, 0b10000, 0b10000, 0b10000],
    "t": [0b01000, 0b01000, 0b11110, 0b01000, 0b01000, 0b01001, 0b00110],
    "v": [0, 0, 0b10001, 0b10001, 0b10001, 0b01010, 0b00100],
    "y": [0, 0, 0b10001, 0b10001, 0b01111, 0b00001, 0b01110],
    "f": [0b00110, 0b01000, 0b11110, 0b01000, 0b01000, 0b01000, 0b01000],
    "k": [0b10000, 0b10000, 0b10010, 0b10100, 0b11000, 0b10100, 0b10010],
    "x": [0, 0, 0b10001, 0b01010, 0b00100, 0b01010, 0b10001],
    "w": [0, 0, 0b10001, 0b10001, 0b10101, 0b10101, 0b01010],
    "z": [0, 0, 0b11111, 0b00010, 0b00100, 0b01000, 0b11111],
    "S": [0b01111, 0b10000, 0b10000, 0b01110, 0b00001, 0b00001, 0b11110],
    "/": [0b00001, 0b00010, 0b00010, 0b00100, 0b01000, 0b01000, 0b10000],
    "<": [0b00010, 0b00100, 0b01000, 0b10000, 0b01000, 0b00100, 0b00010],
    ">": [0b01000, 0b00100, 0b00010, 0b00001, 0b00010, 0b00100, 0b01000],
})

INK = (24, 26, 32)
GRID = (215, 218, 226)
COLS = {1000: (59, 91, 240), 3000: (210, 80, 15)}
GUIDE = (120, 190, 120)


def sweep(d, n):
    """[(b, {rg2,perim,cyc,asph,mov} -> mean)] for positive b, sorted."""
    rows = []
    for f in sorted(glob.glob(os.path.join(d, f"n{n}_b*.csv"))):
        tag = os.path.basename(f)[len(f"n{n}_b"):-4]
        try:
            b = float(tag.replace("p", ".").replace("m", "-"))
        except ValueError:
            continue
        cols = {}
        for ci, name in ((1, "rg2"), (2, "perim"), (3, "cyc"), (4, "asph"),
                         (5, "mov")):
            r = series(f, ci)
            if r:
                cols[name] = r[0]
        if "rg2" in cols:
            rows.append((b, cols))
    rows.sort()
    return rows


def entropy(rows, r2ref, n):
    """[(b, -dS/n)] by thermodynamic integration over the full sweep."""
    bs = [b / r2ref for b, _ in rows]
    r2 = [c["rg2"] for _, c in rows]
    i0 = min(range(len(bs)), key=lambda i: abs(bs[i]))
    lz = {i0: 0.0}
    for i in range(i0 + 1, len(bs)):
        lz[i] = lz[i - 1] - 0.5 * (r2[i] + r2[i - 1]) * (bs[i] - bs[i - 1])
    for i in range(i0 - 1, -1, -1):
        lz[i] = lz[i + 1] + 0.5 * (r2[i] + r2[i + 1]) * (bs[i + 1] - bs[i])
    out = []
    for i, (b, c) in enumerate(rows):
        ds = (lz[i] + bs[i] * r2[i] - (lz[i0] + bs[i0] * r2[i0])) / n
        out.append((b, -ds))
    return out


class Panel:
    def __init__(self, cv, px, py, pw, ph, xlim, ylim, logx, logy):
        self.cv, self.px, self.py, self.pw, self.ph = cv, px, py, pw, ph
        self.xlim, self.ylim, self.logx, self.logy = xlim, ylim, logx, logy

    def X(self, x):
        a, b = self.xlim
        t = ((math.log10(x) - math.log10(a)) / (math.log10(b) - math.log10(a))
             if self.logx else (x - a) / (b - a))
        return self.px + t * self.pw

    def Y(self, y):
        a, b = self.ylim
        t = ((math.log10(y) - math.log10(a)) / (math.log10(b) - math.log10(a))
             if self.logy else (y - a) / (b - a))
        return self.py + self.ph - t * self.ph

    def frame(self, title):
        cv = self.cv
        # x grid at decades of b
        e0 = 0 if not self.logx else math.ceil(math.log10(self.xlim[0]))
        if self.logx:
            e = e0
            while 10 ** e <= self.xlim[1]:
                x = self.X(10 ** e)
                cv.line(x, self.py, x, self.py + self.ph, GRID)
                lab = ("1" if e == 0 else "10" if e == 1 else
                       "0.1" if e == -1 else f"10^{e}")
                cv.text(x - 6 * len(lab), self.py + self.ph + 8, lab, INK)
                e += 1
        cv.line(self.px, self.py + self.ph, self.px + self.pw,
                self.py + self.ph, INK)
        cv.line(self.px, self.py, self.px, self.py + self.ph, INK)
        cv.text(self.px + 8, self.py - 22, title, INK, 2)

    def ytick(self, v, lab):
        y = self.Y(v)
        self.cv.line(self.px, y, self.px + self.pw, y, GRID)
        self.cv.text(self.px - 12 * len(lab) - 6, y - 7, lab, INK)

    def curve(self, pts, c, r=4):
        pp = None
        for x, y in pts:
            P = (self.X(x), self.Y(y))
            if pp:
                self.cv.line(pp[0], pp[1], P[0], P[1], c, 2)
            pp = P
        for x, y in pts:
            self.cv.disc(self.X(x), self.Y(y), r, c)


def main(out):
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                     "results", "beta")
    data = {n: sweep(d, n) for n in (1000, 3000)}
    r2ref = {1000: 1379.0, 3000: 5385.0}
    W, H = 1250, 980
    cv = Canvas(W, H)
    PW, PH = 480, 330
    X0, Y0, GX, GY = 95, 70, 130, 130

    def positive(rows, key):
        return [(b, c[key]) for b, c in rows if b > 0]

    # --- panel 1: Rg^2 squeeze (log-log) ---
    p1 = Panel(cv, X0, Y0, PW, PH, (0.1, 30000), (0.15, 1.3), True, True)
    p1.frame("<Rg2> / <Rg2>(0)   (log-log)")
    for v in (0.2, 0.3, 0.5, 1.0):
        p1.ytick(v, f"{v:.1f}")
    for n, rows in data.items():
        r20 = [c["rg2"] for b, c in rows if b == 0][0]
        p1.curve([(b, v / r20) for b, v in positive(rows, "rg2")], COLS[n])
    # b^-0.2 guide through (30, 0.55)
    p1.curve([(b, 0.62 * b ** -0.2) for b in (3, 3000)], GUIDE, r=0)
    cv.text(p1.X(60), p1.Y(0.75), "slope -0.2", GUIDE)

    # --- panel 2: entropy cost (log-log) ---
    p2 = Panel(cv, X0 + PW + GX, Y0, PW, PH, (0.1, 30000), (1e-5, 0.3),
               True, True)
    p2.frame("-dS/n  (nats/cell, log-log)")
    for e in (-4, -3, -2, -1):
        p2.ytick(10 ** e, f"10^{e}")
    for n, rows in data.items():
        ent = entropy(rows, r2ref[n], n)
        p2.curve([(b, v) for b, v in ent if b > 0 and v > 1e-6], COLS[n])

    # --- panel 3: local observables relative to b=0 ---
    p3 = Panel(cv, X0, Y0 + PH + GY, PW, PH, (0.1, 30000), (0.85, 1.45),
               True, False)
    p3.frame("locals relative to b=0")
    for v in (0.9, 1.0, 1.1, 1.2, 1.3, 1.4):
        p3.ytick(v, f"{v:.1f}")
    styles = {"perim": (59, 91, 240), "cyc": (210, 80, 15),
              "mov": (90, 160, 90)}
    for n, rows in data.items():
        base = {k: [c[k] for b, c in rows if b == 0][0]
                for k in ("perim", "cyc", "mov")}
        r = 4 if n == 3000 else 2
        for k, col in styles.items():
            p3.curve([(b, v / base[k]) for b, v in positive(rows, k)], col,
                     r=r)
    cv.text(p3.px + 14, p3.py + 12, "cycles/cell", styles["cyc"])
    cv.text(p3.px + 14, p3.py + 32, "movable frac", styles["mov"])
    cv.text(p3.px + 14, p3.py + 52, "perim/cell", styles["perim"])
    cv.text(p3.px + 14, p3.py + 72, "small dot n=1000, big n=3000", INK)

    # --- panel 4: asphericity ---
    p4 = Panel(cv, X0 + PW + GX, Y0 + PH + GY, PW, PH, (0.1, 30000),
               (0.0, 0.45), True, False)
    p4.frame("asphericity")
    for v in (0.0, 0.1, 0.2, 0.3, 0.4):
        p4.ytick(v, f"{v:.1f}")
    for n, rows in data.items():
        p4.curve(positive(rows, "asph"), COLS[n])

    for i, (n, col) in enumerate(sorted(COLS.items())):
        cv.text(X0 + PW - 150, Y0 + 14 + 20 * i, f"n={n}", col)
    cv.text(X0 + 130, H - 34,
            "b = beta <Rg2>(0)      collapse at b ~ 0.18 n", INK, 2)
    write_png(out, W, H, cv.rows)
    print(out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "../results/betafig.png")
