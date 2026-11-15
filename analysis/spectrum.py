"""The exact relaxation time of the chain on small shapes.

For sizes where every polyomino can be listed, `manyomino -mode trans` counts the transitions
between all pairs of shapes.  The chain is reversible with a uniform stationary law, so the
symmetrised counts give its transition matrix; the second-largest eigenvalue lambda_2 is the
slowest decay rate of ANY property of the shape, not just the ones we thought to measure:
relaxation time = 1 / (1 - lambda_2) moves.

The same matrix also says how much of that slowest mode the observable library in
src/features.h can see: the best linear combination of the observables and their squares has
relaxation time 1 / (1 - mu), with mu its Rayleigh quotient (mu <= lambda_2).

  python3 analysis/spectrum.py            # reads spec/t<n>.{states,trans}, writes analysis/spectrum.json
"""
import glob, json, math, os, re, sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

ROOT = os.path.join(os.path.dirname(__file__), "..")
NAMES = None


def load(prefix):
    global NAMES
    if NAMES is None:
        NAMES = open(glob.glob(os.path.join(ROOT, "conv", "*.feat.names"))[0]).read().split()[1:] if glob.glob(os.path.join(ROOT, "conv", "*.feat.names")) else None
    st = np.fromfile(prefix + ".states")
    tr = np.fromfile(prefix + ".trans").reshape(-1, 3)
    N = int(tr[:, :2].max()) + 1
    st = st.reshape(N, -1)
    return st[:, 0], st[:, 1:], tr


def matrix(visits, tr):
    """symmetric transition matrix from transition counts, using the exact uniform stationary law"""
    N = len(visits); M = visits.sum()
    i = tr[:, 0].astype(np.int64); j = tr[:, 1].astype(np.int64)
    C = sp.coo_matrix((tr[:, 2], (i, j)), shape=(N, N)).tocsr()
    F = (C + C.T) * (N / (2.0 * M))                       # P_ij for i != j
    d = 1.0 - np.asarray(F.sum(1)).ravel()
    return (F + sp.diags(d)).tocsr()


def analyse(n, prefixes):
    parts = [load(p) for p in prefixes]
    # shapes are numbered in order of first visit, differently in each run: match them by their features
    visits, X, tr = parts[0]
    if len(parts) > 1:
        key = lambda row: row.tobytes()
        idx = {key(r): k for k, r in enumerate(np.round(X, 9))}
        for v2, X2, tr2 in parts[1:]:
            m = np.array([idx[key(r)] for r in np.round(X2, 9)])
            visits = visits + np.bincount(m, weights=v2, minlength=len(visits))
            tr = np.vstack([tr, np.column_stack([m[tr2[:, 0].astype(int)], m[tr2[:, 1].astype(int)], tr2[:, 2]])])
    N = len(visits); P = matrix(visits, tr)
    k = min(12, N - 2)
    w, V = eigsh(P, k=k, which="LA") if N > 400 else np.linalg.eigh(P.toarray())
    o = np.argsort(-w); w = w[o]; V = V[:, o]
    lam2 = float(w[1]); trel = 1 / (1 - lam2)
    # how degenerate is the slowest level (lattice symmetry gives multiplets)
    mult = int((w[1] - w[1:8] < 1e-3 * (1 - w[1])).sum())
    # the observable library: features and their squares, centred, constant columns dropped
    Z = X[:, X.std(0) > 0]; Z = (Z - Z.mean(0)) / Z.std(0)
    Phi = np.hstack([Z, Z ** 2]); Phi = Phi - Phi.mean(0)
    u, s, _ = np.linalg.svd(Phi, full_matrices=False); B = u[:, s > 1e-9 * s[0]]       # orthonormal basis of the library's span
    mu = float(np.linalg.eigvalsh(B.T @ (P @ B)).max()); tlib = 1 / (1 - mu)
    # single observables: Rayleigh quotient of each (a lower bound on its slowest component)
    q = np.array([float(z @ (P @ z) / (z @ z)) for z in (Z - Z.mean(0)).T]); tsingle = 1 / (1 - q.max())
    # share of the slowest eigenvectors lying inside the library's span
    r2 = [float(((B.T @ V[:, j]) ** 2).sum()) for j in range(1, 1 + mult)]
    chi = float((visits.std() / visits.mean()) ** 2 * visits.mean())                   # visits should be Poisson-like at best; a sanity number
    out = {"n": n, "shapes": N, "moves": float(visits.sum()), "transitions": int(len(tr)), "lambda2": lam2, "t_rel": trel, "t_rel_per_cell": trel / n,
           "multiplicity": mult, "eigs": [float(x) for x in w[:8]], "t_library": tlib, "t_library_per_cell": tlib / n, "library_dim": int(B.shape[1]),
           "t_single_per_cell": tsingle / n, "r2_slowest": r2, "min_visits": float(visits.min()), "max_visits": float(visits.max())}
    print(f"n={n}: {N} shapes, {visits.sum():.3g} moves, {len(tr)} transitions; relaxation time {trel:.1f} moves = {trel / n:.2f} per cell "
          f"(slowest level {mult}-fold); next levels " + " ".join(f"{1 / (1 - x) / n:.2f}" for x in w[1 + mult:1 + mult + 3]) +
          f"; best library combination {tlib / n:.2f} per cell ({B.shape[1]} dimensions; captures {min(r2):.3f}..{max(r2):.3f} of the slowest eigenvectors); best single observable {tsingle / n:.2f}")
    return out


def main():
    R = []
    for n in range(4, 13):
        pre = sorted(p[:-7] for p in glob.glob(os.path.join(ROOT, "spec", f"t{n}*.states")) if re.fullmatch(rf"t{n}[a-z]?", os.path.basename(p)[:-7]))
        if pre:
            R.append(analyse(n, pre))
            if len(pre) > 1:                              # each run separately: how reproducible is lambda_2?
                R[-1]["t_rel_per_cell_runs"] = [analyse(n, [p])["t_rel_per_cell"] for p in pre]
    json.dump(R, open(os.path.join(ROOT, "analysis", "spectrum.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
