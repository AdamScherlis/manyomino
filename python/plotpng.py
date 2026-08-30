"""Rasterize the nu-fit log-log plot to PNG (stdlib only): grid, points,
error bars, fit line, and labels in an embedded 5x7 bitmap font."""
import math
import sys

from render import write_png
from svgplot import collect, wls_fit

# 5x7 font: rows as 5-bit ints, MSB = leftmost pixel
FONT = {
    "0": [0b01110,0b10001,0b10011,0b10101,0b11001,0b10001,0b01110],
    "1": [0b00100,0b01100,0b00100,0b00100,0b00100,0b00100,0b01110],
    "2": [0b01110,0b10001,0b00001,0b00010,0b00100,0b01000,0b11111],
    "3": [0b11110,0b00001,0b00001,0b01110,0b00001,0b00001,0b11110],
    "4": [0b00010,0b00110,0b01010,0b10010,0b11111,0b00010,0b00010],
    "5": [0b11111,0b10000,0b11110,0b00001,0b00001,0b10001,0b01110],
    "6": [0b00110,0b01000,0b10000,0b11110,0b10001,0b10001,0b01110],
    "7": [0b11111,0b00001,0b00010,0b00100,0b01000,0b01000,0b01000],
    "8": [0b01110,0b10001,0b10001,0b01110,0b10001,0b10001,0b01110],
    "9": [0b01110,0b10001,0b10001,0b01111,0b00001,0b00010,0b01100],
    ".": [0,0,0,0,0,0b00100,0b00100],
    ",": [0,0,0,0,0,0b00100,0b01000],
    "+": [0,0b00100,0b00100,0b11111,0b00100,0b00100,0],
    "-": [0,0,0,0b11111,0,0,0],
    "=": [0,0,0b11111,0,0b11111,0,0],
    "^": [0b00100,0b01010,0b10001,0,0,0,0],
    "~": [0,0,0b01000,0b10101,0b00010,0,0],
    "(": [0b00010,0b00100,0b01000,0b01000,0b01000,0b00100,0b00010],
    ")": [0b01000,0b00100,0b00010,0b00010,0b00010,0b00100,0b01000],
    " ": [0,0,0,0,0,0,0],
    "R": [0b11110,0b10001,0b10001,0b11110,0b10100,0b10010,0b10001],
    "g": [0,0,0b01111,0b10001,0b01111,0b00001,0b01110],
    "n": [0,0,0b10110,0b11001,0b10001,0b10001,0b10001],
    "u": [0,0,0b10001,0b10001,0b10001,0b10011,0b01101],
    "c": [0,0,0b01110,0b10000,0b10000,0b10001,0b01110],
    "e": [0,0,0b01110,0b10001,0b11111,0b10000,0b01110],
    "l": [0b01100,0b00100,0b00100,0b00100,0b00100,0b00100,0b01110],
    "s": [0,0,0b01111,0b10000,0b01110,0b00001,0b11110],
}


class Canvas:
    def __init__(self, w, h):
        self.w, self.h = w, h
        self.rows = [bytearray(b"\xff" * (3 * w)) for _ in range(h)]

    def px(self, x, y, c):
        if 0 <= x < self.w and 0 <= y < self.h:
            r = self.rows[int(y)]
            o = 3 * int(x)
            r[o], r[o + 1], r[o + 2] = c

    def rect(self, x0, y0, x1, y1, c):
        for y in range(int(y0), int(y1) + 1):
            for x in range(int(x0), int(x1) + 1):
                self.px(x, y, c)

    def line(self, x0, y0, x1, y1, c, wpx=1):
        n = int(max(abs(x1 - x0), abs(y1 - y0))) + 1
        for i in range(n + 1):
            t = i / n
            x, y = x0 + (x1 - x0) * t, y0 + (y1 - y0) * t
            self.rect(x - wpx / 2, y - wpx / 2, x + wpx / 2, y + wpx / 2, c)

    def disc(self, cx, cy, r, c):
        for y in range(int(cy - r), int(cy + r) + 1):
            for x in range(int(cx - r), int(cx + r) + 1):
                if (x - cx) ** 2 + (y - cy) ** 2 <= r * r:
                    self.px(x, y, c)

    def text(self, x, y, s, c, scale=2):
        for ch in s:
            glyph = FONT.get(ch)
            if glyph:
                for gy, bits in enumerate(glyph):
                    for gx in range(5):
                        if bits >> (4 - gx) & 1:
                            self.rect(x + gx * scale, y + gy * scale,
                                      x + gx * scale + scale - 1,
                                      y + gy * scale + scale - 1, c)
            x += 6 * scale


def main(out):
    pts = collect()
    a, nu, snu = wls_fit(pts)
    W, H = 960, 640
    L, R, T, B = 90, 30, 40, 70
    x0, x1 = math.log10(70), math.log10(400000)
    y0, y1 = math.log10(6), math.log10(1100)
    X = lambda n: L + (math.log10(n) - x0) / (x1 - x0) * (W - L - R)
    Y = lambda v: H - B - (math.log10(v) - y0) / (y1 - y0) * (H - T - B)
    cv = Canvas(W, H)
    GRID, INK, ACC, PT = (210, 213, 222), (24, 26, 32), (210, 80, 15), (59, 91, 240)
    for e in range(2, 6):
        cv.line(X(10**e), T, X(10**e), H - B, GRID)
        cv.text(X(10**e) - 18, H - B + 12, f"10^{e}", INK)
    for e in (1, 2, 3):
        cv.line(L, Y(10**e), W - R, Y(10**e), GRID)
        cv.text(L - 52, Y(10**e) - 7, f"10^{e}", INK)
    # frame
    cv.line(L, H - B, W - R, H - B, INK)
    cv.line(L, T, L, H - B, INK)
    # fit line clipped to plot
    nlo, nhi = 80, 220000
    cv.line(X(nlo), Y(math.exp(a) * nlo**nu), X(nhi), Y(math.exp(a) * nhi**nu), ACC, 2)
    for n, m, se in pts:
        cv.line(X(n), Y(m - 2 * se), X(n), Y(m + 2 * se), PT, 2)
        cv.disc(X(n), Y(m), 5, PT)
    cv.text(L + 20, T + 14, f"Rg ~ n^nu,  nu = {nu:.4f} +- {snu:.4f}", INK, 2)
    cv.text((L + W - R) // 2 - 40, H - 30, "n (cells)", INK, 2)
    write_png(out, W, H, cv.rows)
    print(f"{out}: nu = {nu:.4f} +- {snu:.4f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "../results/nufit.png")
