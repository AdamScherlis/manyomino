"""Render a polyomino cell dump ("x y" per line) to PNG.  Stdlib only.

Modes:
  bw    black cells on white, tight bounding box
  dist  cells colored by graph distance from an extremal root cell,
        which makes the branch structure visible

Usage: python3 render.py cells.txt out.png [--mode bw|dist] [--scale S]
Scale defaults to 1 px/cell for n >= 100k bounding-box area, else enough
to make the image ~800 px wide.
"""

import struct
import sys
import zlib
from collections import deque


def load(path):
    cells = []
    with open(path) as f:
        for line in f:
            x, y = line.split()
            cells.append((int(x), int(y)))
    return cells


def write_png(path, width, height, rgb_rows):
    """rgb_rows: list of bytearrays, each 3*width bytes."""
    def chunk(tag, data):
        c = struct.pack(">I", len(data)) + tag + data
        return c + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)

    raw = b"".join(b"\x00" + bytes(row) for row in rgb_rows)
    png = (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(raw, 6))
        + chunk(b"IEND", b"")
    )
    with open(path, "wb") as f:
        f.write(png)


def hsv_to_rgb(h, s, v):
    i = int(h * 6.0) % 6
    f = h * 6.0 - int(h * 6.0)
    p, q, t = v * (1 - s), v * (1 - f * s), v * (1 - (1 - f) * s)
    r, g, b = [(v, t, p), (q, v, p), (p, v, t), (p, q, v), (t, p, v), (v, p, q)][i]
    return int(r * 255), int(g * 255), int(b * 255)


def graph_distances(cells):
    cs = set(cells)
    # root at an extremal cell (min y, then min x) so the gradient sweeps
    root = min(cells, key=lambda c: (c[1], c[0]))
    dist = {root: 0}
    dq = deque([root])
    while dq:
        x, y = dq.popleft()
        d = dist[(x, y)]
        for nb in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
            if nb in cs and nb not in dist:
                dist[nb] = d + 1
                dq.append(nb)
    return dist


def render(cells, out, mode="bw", scale=None):
    minx = min(x for x, y in cells)
    miny = min(y for x, y in cells)
    maxx = max(x for x, y in cells)
    maxy = max(y for x, y in cells)
    bw_, bh = maxx - minx + 1, maxy - miny + 1
    if scale is None:
        scale = max(1, 800 // max(bw_, bh))
    pad = 1
    width, height = (bw_ + 2 * pad) * scale, (bh + 2 * pad) * scale

    rows = [bytearray(b"\xff" * (3 * width)) for _ in range(height)]

    if mode == "dist":
        dist = graph_distances(cells)
        dmax = max(dist.values()) or 1
        colors = {
            c: hsv_to_rgb(0.75 * dist[c] / dmax, 0.85, 0.85) for c in cells
        }
    else:
        colors = {c: (0, 0, 0) for c in cells}

    for (x, y) in cells:
        r, g, b = colors[(x, y)]
        px, py = (x - minx + pad) * scale, (y - miny + pad) * scale
        for dy in range(scale):
            row = rows[py + dy]
            for dx in range(scale):
                o = 3 * (px + dx)
                row[o] = r
                row[o + 1] = g
                row[o + 2] = b

    write_png(out, width, height, rows)
    return width, height


if __name__ == "__main__":
    args = sys.argv[1:]
    src, out = args[0], args[1]
    mode = "bw"
    scale = None
    if "--mode" in args:
        mode = args[args.index("--mode") + 1]
    if "--scale" in args:
        scale = int(args[args.index("--scale") + 1])
    cells = load(src)
    w, h = render(cells, out, mode, scale)
    print(f"{out}: {len(cells)} cells -> {w}x{h} px ({mode})")
