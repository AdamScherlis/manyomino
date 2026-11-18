"""How a doubled shape settles, and how close to a sampled shape it ends up.

Reads the shapes saved by scripts/zoom_relax.sh (zoomtest/relax/TAG_tT.bin: the sample
prod/snap_TAG.bin doubled and relaxed for T sweeps) and, for each T:

  growth    cells of the refined shape over cells of the original
  local     per-cell statistics of the neighbourhood of a cell, relative to sampled shapes
  box       cells per occupied box of side b, relative to sampled shapes of about that size
  changed   share of the original's boxes of side b (in its own cells) whose occupancy differs
            in the refined shape: how far up in scale the picture was redrawn
  clumps    number of loop-connected clumps of each size, relative to equilibrium (64,000 only)

  python3 analysis/zoom_relax.py        # writes analysis/zoom_relax.json
"""
import glob, json, math, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.setrecursionlimit(10000)
from zoom_stats import read, local, boxes, keys
from clump_shapes import clumps
from clump_census import density, SIZES

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
TS = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
BOX = [2, 4, 8, 16, 32, 64]
CH = [1, 2, 4, 8, 16, 32]
NU = 0.6412


def me(v): v = np.asarray(v, float); return float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v)))
def ratio(a, b): return [a[0] / b[0], math.hypot(a[1] / b[0], a[0] * b[1] / b[0] ** 2)]


def changed(parent, fine):
    out = {}
    for b in CH:
        A = np.unique(keys(parent[:, 0] // b, parent[:, 1] // b)); B = np.unique(keys(fine[:, 0] // (2 * b), fine[:, 1] // (2 * b)))
        out[b] = float(np.setxor1d(A, B, assume_unique=True).size / A.size)
    return out


def main():
    out = {"pred": 2 ** (1 / NU), "T": TS, "sets": {}}; eq = density("clumps.json")
    real = {}
    for n in (16000, 64000, 256000):
        fs = sorted(glob.glob(os.path.join(ROOT, "prod", f"snap_?_n{n}_*.bin")))
        real[n] = [dict(local=local(xy), box=boxes(xy, BOX)) for xy in map(read, fs)]; print("sampled", n, len(fs), flush=True)
    loc_names = list(real[16000][0]["local"])
    pooled = {nm: me([r["local"][nm] for n in real for r in real[n]]) for nm in loc_names}
    out["sampled_local"] = pooled
    for n0 in (16000, 64000):
        tags = sorted({re.match(r"(.*)_t1\.bin", os.path.basename(f)).group(1) for f in glob.glob(os.path.join(ROOT, "zoomtest", "relax", f"?_n{n0}_*_t1.bin"))})
        ref = real[4 * n0]; S = {"samples": len(tags), "growth": [], "local": {nm: [] for nm in loc_names}, "box": {b: [] for b in BOX}, "changed": {b: [] for b in CH}, "clumps": {}}
        parents = {t: read(os.path.join(ROOT, "prod", f"snap_{t}.bin")) for t in tags}
        for T in TS:
            rows = []; cl = {s: 0 for s in SIZES}; cells = 0
            for t in tags:
                xy = read(os.path.join(ROOT, "zoomtest", "relax", f"{t}_t{T}.bin"))
                rows.append(dict(m=len(xy) / n0, local=local(xy), box=boxes(xy, BOX), ch=changed(parents[t], xy)))
                if n0 == 64000 and T in (1, 4, 16, 64, 256):
                    sz = np.array([len(c) for c in clumps(xy)]); cells += len(xy)
                    for s in SIZES: cl[s] += int((sz >= s).sum())
            S["growth"].append(me([r["m"] for r in rows]))
            for nm in loc_names: S["local"][nm].append(ratio(me([r["local"][nm] for r in rows]), pooled[nm]))
            for b in BOX: S["box"][b].append(ratio(me([r["box"][b] for r in rows]), me([r["box"][b] for r in ref])))
            for b in CH: S["changed"][b].append(me([r["ch"][b] for r in rows]))
            if cells: S["clumps"][T] = {s: [cl[s], eq[s] * cells] for s in SIZES}
            print(n0, T, "growth %.4f ± %.4f" % tuple(S["growth"][-1]), "squares %.4f ± %.4f" % tuple(S["local"]["squares"][-1]),
                  "changed", " ".join(f"{b}:{S['changed'][b][-1][0]:.3f}" for b in CH), flush=True)
        out["sets"][n0] = S
    json.dump(out, open(os.path.join(ROOT, "analysis", "zoom_relax.json"), "w"))


if __name__ == "__main__":
    main()
