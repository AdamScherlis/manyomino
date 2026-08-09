"""Inflate a snapshot 2x: each cell (x, y) -> the 2x2 block at (2x, 2y).
The result is a valid connected polyomino of 4n cells whose large-scale
geometry matches the source; only sub-inflation-scale (fast local) modes are
out of equilibrium.  Usage: python3 inflate.py in.txt out.txt"""
import sys

cells = [tuple(map(int, l.split())) for l in open(sys.argv[1])]
with open(sys.argv[2], "w") as f:
    for (x, y) in cells:
        for dx in (0, 1):
            for dy in (0, 1):
                f.write(f"{2*x+dx} {2*y+dy}\n")
print(f"{len(cells)} -> {4*len(cells)} cells")
