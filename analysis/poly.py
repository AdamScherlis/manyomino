"""Shared helpers: load snapshots, graph distances, rendering, compact encoding."""
import struct
import numpy as np
from collections import deque

DX = (1, 0, -1, 0)
DY = (0, 1, 0, -1)


def load_bin(fn):
    with open(fn, "rb") as f:
        n = struct.unpack("<i", f.read(4))[0]
        a = np.frombuffer(f.read(8 * n), dtype="<i4").reshape(n, 2)
    x = a[:, 0] - a[:, 0].min()
    y = a[:, 1] - a[:, 1].min()
    return x.astype(np.int64), y.astype(np.int64)


def neighbours(x, y):
    """nbr[i, d] = index of the cell at (x+DX[d], y+DY[d]) or -1."""
    n = len(x)
    big = int(y.max()) + 3
    key = x * big + y
    order = np.argsort(key)
    skey = key[order]
    nbr = np.full((n, 4), -1, dtype=np.int64)
    for d in range(4):
        k = (x + DX[d]) * big + (y + DY[d])
        pos = np.searchsorted(skey, k)
        pos[pos >= n] = n - 1
        hit = skey[pos] == k
        nbr[hit, d] = order[pos[hit]]
    return nbr


def bfs(nbr, src):
    n = len(nbr)
    dist = np.full(n, -1, dtype=np.int64)
    parent = np.full(n, -1, dtype=np.int64)
    nb = nbr.tolist()
    dl = [-1] * n
    pl = [-1] * n
    dl[src] = 0
    q = deque([src])
    while q:
        v = q.popleft()
        dv = dl[v] + 1
        for w in nb[v]:
            if w >= 0 and dl[w] < 0:
                dl[w] = dv
                pl[w] = v
                q.append(w)
    dist[:] = dl
    parent[:] = pl
    return dist, parent


def centre(nbr):
    """Midpoint of a (double-sweep) longest shortest path; returns (centre, diameter)."""
    d0, _ = bfs(nbr, 0)
    a = int(np.argmax(d0))
    da, pa = bfs(nbr, a)
    b = int(np.argmax(da))
    diam = int(da[b])
    v = b
    for _ in range(diam // 2):
        v = int(pa[v])
    return v, diam


# depth ramp: centre (pale gold) -> orange -> rose -> violet -> blue at the tips
RAMP = np.array([[255, 236, 170], [255, 170, 96], [240, 96, 140], [160, 92, 236], [84, 116, 255]], float)


def ramp(t):
    t = np.clip(np.asarray(t, float), 0, 1) * (len(RAMP) - 1)
    i = np.minimum(t.astype(int), len(RAMP) - 2)
    f = (t - i)[..., None]
    return (RAMP[i] * (1 - f) + RAMP[i + 1] * f).astype(np.uint8)


def render(x, y, depth, px=1, pad=2, bg=(13, 15, 22), max_side=None):
    """RGB image array. px = pixels per cell (>=1), or downsample so the long side <= max_side."""
    from PIL import Image
    w, h = int(x.max()) + 1, int(y.max()) + 1
    col = ramp(depth / max(1, depth.max()))
    if max_side and max(w, h) * px > max_side and px == 1:
        f = int(np.ceil(max(w, h) / max_side))
        W, H = (w + f - 1) // f, (h + f - 1) // f
        # per coarse pixel: coverage and mean colour
        idx = (h - 1 - y) // f * W + x // f
        cnt = np.bincount(idx, minlength=W * H).astype(float)
        rgb = np.stack([np.bincount(idx, weights=col[:, c].astype(float), minlength=W * H) for c in range(3)], 1)
        m = cnt > 0
        rgb[m] /= cnt[m, None]
        # coverage -> opacity, boosted so thin strands stay visible
        a = np.clip(cnt / (f * f) * 3.0, 0, 1) ** 0.45
        out = np.array(bg, float)[None, :] * (1 - a[:, None]) + rgb * a[:, None]
        img = out.reshape(H, W, 3).astype(np.uint8)
        img = np.pad(img, ((pad, pad), (pad, pad), (0, 0)), constant_values=0)
        img[:pad], img[-pad:], img[:, :pad], img[:, -pad:] = bg, bg, bg, bg
        return Image.fromarray(img)
    img = np.empty((h + 2 * pad, w + 2 * pad, 3), np.uint8)
    img[:] = bg
    img[h - 1 - y + pad, x + pad] = col
    im = Image.fromarray(img)
    if px > 1:
        im = im.resize((im.width * px, im.height * px), Image.NEAREST)
    return im


def encode_tree(x, y, nbr, root):
    """Compact encoding: BFS spanning tree from root, one nibble per cell in DFS preorder
    giving which of the 4 directions lead to tree children.  Depth in this tree is the graph
    distance from the root, so the decoder recovers coordinates and depth together."""
    n = len(x)
    dist, parent = bfs(nbr, root)
    mask = np.zeros(n, np.uint8)
    nb = nbr.tolist()
    par = parent.tolist()
    for v in range(n):
        p = par[v]
        if p >= 0:
            mask[p] |= 1 << nb[p].index(v)
    out = np.zeros(n, np.uint8)
    stack = [root]
    i = 0
    ml = mask.tolist()
    while stack:
        v = stack.pop()
        m = ml[v]
        out[i] = m
        i += 1
        for d in (3, 2, 1, 0):
            if m >> d & 1:
                stack.append(nb[v][d])
    assert i == n
    if n % 2:
        out = np.append(out, 0)
    packed = (out[0::2] | (out[1::2] << 4)).astype(np.uint8)
    hdr = struct.pack("<iii", n, int(x[root]), int(y[root]))
    return hdr + packed.tobytes(), dist
