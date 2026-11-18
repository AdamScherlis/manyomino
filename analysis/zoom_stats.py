"""Scale-by-scale statistics of a shape, for comparing refined shapes (src/zoom.cpp) with
sampled ones.

    local(xy)     quantities per cell that depend only on the neighbourhood of a cell
    boxes(xy)     cells per occupied box of side b, for b = 1, 2, 4, ...  (a local quantity as
                  long as b is well below the size of the shape)
"""
import sys, os, json, glob, math
import numpy as np

BS = [2, 4, 8, 16, 32, 64, 128, 256]


def read(fn):
    a = np.fromfile(fn, dtype="<i4"); n = a[0]; return a[1:1 + 2 * n].reshape(n, 2).astype(np.int64)


def keys(x, y): return (x + (1 << 30)) * (1 << 32) + (y + (1 << 30))


def local(xy):
    x, y = xy[:, 0], xy[:, 1]; k = np.sort(keys(x, y)); n = len(k)
    def has(dx, dy): return np.isin(keys(x + dx, y + dy), k, assume_unique=False)
    E, W, N, S = has(1, 0), has(-1, 0), has(0, 1), has(0, -1)
    deg = E.astype(int) + W + N + S
    sq = (E & N & has(1, 1)).sum()                      # 2x2 blocks
    per = np.unique(np.concatenate([keys(x + dx, y + dy)[~m] for (dx, dy), m in (((1, 0), E), ((-1, 0), W), ((0, 1), N), ((0, -1), S))])).size
    edges = deg.sum() / 2
    out = {"perim": per / n, "squares": sq / n, "loops": (edges - n + 1) / n}
    for d in (1, 2, 3, 4): out[f"deg{d}"] = float((deg == d).mean())
    # straight cells: exactly two neighbours, opposite each other
    out["straight"] = float(((deg == 2) & ((E & W) | (N & S))).mean())
    return {a: float(b) for a, b in out.items()}


def boxes(xy, bs=BS):
    """cells per occupied box, averaged over four positions of the grid of boxes"""
    x, y = xy[:, 0] - xy[:, 0].min(), xy[:, 1] - xy[:, 1].min(); n = len(x); out = {}
    for b in bs:
        c = [np.unique(keys((x + ox) // b, (y + oy) // b)).size for ox, oy in ((0, 0), (b // 2, 0), (0, b // 2), (b // 2, b // 2))]
        out[b] = n / float(np.mean(c))
    return out


def stats(fn):
    xy = read(fn); x, y = xy[:, 0].astype(float), xy[:, 1].astype(float)
    return {"n": len(xy), "rg2": float(x.var() + y.var()), "local": local(xy), "box": boxes(xy)}


if __name__ == "__main__":
    out = {}
    for fn in sys.argv[2:]:
        out[os.path.basename(fn)] = stats(fn); print(os.path.basename(fn), out[os.path.basename(fn)]["n"], flush=True, file=sys.stderr)
    json.dump(out, open(sys.argv[1], "w"))
