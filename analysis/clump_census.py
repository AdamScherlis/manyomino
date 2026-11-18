"""Do saved shapes hold the equilibrium number of clumps?  Counts the clumps of each size in
snapshots (.bin) and sets them against the densities per cell measured from the long clump
records (analysis/clumps.json, analysis/clumps_a1.json).
    python3 analysis/clump_census.py prod/snap_*_n1024000_*.bin"""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clump_shapes import read, clumps

SIZES = (4, 8, 16, 32, 64, 128, 182, 256)


def density(name):
    """clumps of at least s cells, per cell: weighted mean over the shape sizes of one record"""
    J = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), name))); out = {}
    for s in SIZES:
        v = [(r[1], r[2]) for z in J["sizes"].values() for r in z["density"] if r[0] == s and r[2] > 0]
        if v: w = sum(1 / e ** 2 for _, e in v); out[s] = sum(d / e ** 2 for d, e in v) / w
    return out


if __name__ == "__main__":
    eq = density("clumps.json"); tot = {s: 0 for s in SIZES}; cells = 0; big = []
    for fn in sys.argv[1:]:
        xy = read(fn); sz = np.array([len(c) for c in clumps(xy)]); cells += len(xy); big.append(int(sz.max()))
        for s in SIZES: tot[s] += int((sz >= s).sum())
        print(os.path.basename(fn), " ".join(f"{int((sz >= s).sum()):6d}" for s in SIZES), " largest", big[-1], flush=True)
    print(f"{len(big)} shapes, {cells:,} cells; largest clump {np.mean(big):.0f} ± {np.std(big, ddof=1) / math.sqrt(len(big)):.0f} on average")
    for s in SIZES:
        e = eq[s] * cells
        print(f"  at least {s:3d} cells: {tot[s]:9,d} found, {e:12,.1f} expected, ratio {tot[s] / e:.3f}, {(tot[s] - e) / math.sqrt(e):+.1f} counting errors")
