"""Can a classifier tell the zoom out's pictures from real ones?

analysis/e2e_test.py compares a fixed list of numbers.  This is the stricter test: two
classifiers are trained to tell made pictures from real ones and scored on pictures they have
not seen, from samples they have not seen (the 26 samples are split into folds, so that nothing
learnt about one sample's large-scale shape helps).

  summary   a logistic regression on about 150 numbers per picture: cells by distance, counts of
            local patterns in each ring, cells layer by layer across each seam, boxes occupied at
            several scales, and how the picture falls into pieces
  patches   a small neural network that looks at 9x9 windows centred on cells, told which ring
            the cell is in and how far it is from a seam; a picture's score is the mean over
            3,000 of its windows

For each: the area under the ROC curve on held-out pictures (0.5 = cannot tell, 1 = always can)
and the share of pictures called correctly.

  python3 analysis/tell_test.py [TAG=final] [FOLDS=4] [SEED=1]     # reads zoomtest/e2e/TAG_pics.npz; with TAG final writes analysis/tell_test.json
  python3 analysis/tell_test.py TAG FOLDS SEED null                # a control: the real pictures of half the samples against those of the other half
"""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e2e_test import pieces, rad, ROOT, B
from convergence import train_mlp, auc


def summary(g, with_pieces=True):
    R = g.shape[0] // 2; r = rad(R); L = int(round(math.log2(R / B))) - 1; gi = g.astype(np.int32); f = []
    p = np.pad(gi, 1); nb = p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2]; dg = p[2:, 2:] + p[:-2, :-2] + p[2:, :-2] + p[:-2, 2:]
    for k in range(1, L + 1):
        a = B * 2 ** k; m = (r >= a) & (r < 2 * a) & (r < R - 1); c = gi * m; n = max(1, c.sum()); f.append(math.log(n))
        h = np.bincount((nb * 5 + dg)[c > 0], minlength=25)[:25] / n; f += h.tolist()            # by number of side and corner neighbours
        f.append(((nb > 0) * (1 - gi) * m).sum() / n)
        for q in range(4):                                                                        # quarters of the ring
            lo, hi = a * 2 ** (q / 4), a * 2 ** ((q + 1) / 4); f.append(gi[(r >= lo) & (r < hi)].sum() / n)
        for d in range(-4, 4): f.append(gi[r == a + d].sum() / n * a)                             # layer by layer across the seam
        f.append(((r == a - 1) * gi * np.maximum.reduce([np.roll(gi * (r == a), s_, ax) for s_, ax in ((1, 0), (-1, 0), (1, 1), (-1, 1))])).sum() / n * a)
    for s in (2, 4, 8, 16, 32): f.append(math.log(1 + gi.reshape(2 * R // s, s, 2 * R // s, s).max((1, 3)).sum()))
    if with_pieces:
        Rf = R; Rc = 3 * R // 4; gc = g[Rf - Rc:Rf + Rc, Rf - Rc:Rf + Rc]; lab, s = pieces(gc, True); rim = set(np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]])).tolist())
        fl = [m for k, m in enumerate(s) if k + 1 not in rim and lab[Rc, Rc] != k + 1]; s = np.array(s)
        f += [math.log(len(s)), float((s == 1).sum()), float(((s > 1) & (s < 10)).sum()), float(s.max() / s.sum()), float(len(fl)), math.log(1 + sum(fl))]
    return f


def patches(g, n, rng, w=4):
    """n windows of (2w+1)^2 sites centred on cells, with the ring of the cell and its distance from the nearest seam"""
    R = g.shape[0] // 2; r = rad(R); ys, xs = np.nonzero(g & (r >= B + w) & (r < R - w - 1)); k = rng.integers(len(ys), size=n); ys, xs = ys[k], xs[k]
    idx = np.arange(-w, w + 1); P = g[ys[:, None, None] + idx[None, :, None], xs[:, None, None] + idx[None, None, :]].reshape(n, -1).astype(np.float32)
    rr = r[ys, xs].astype(np.float64); lv = np.floor(np.log2(rr / B)); seam = np.minimum(rr - B * 2 ** lv, B * 2 ** (lv + 1) - rr)       # distance to the nearer seam, in cells
    extra = np.stack([lv / 4, np.minimum(seam, 16) / 16, (seam < 2), (seam < 4), (rr - B * 2 ** lv < 4)], 1).astype(np.float32)
    return np.hstack([P * 2 - 1, extra])


def logistic(X, y, l2=1.0, steps=400):
    w = np.zeros(X.shape[1]); b = 0.0; m = np.zeros_like(w); v = np.zeros_like(w); mb = vb = 0.0
    for t in range(1, steps + 1):
        z = X @ w + b; g = 1 / (1 + np.exp(-np.clip(z, -30, 30))) - y; gw = X.T @ g / len(y) + l2 * w / len(y); gb = g.mean()
        m = 0.9 * m + 0.1 * gw; v = 0.999 * v + 0.001 * gw ** 2; mb = 0.9 * mb + 0.1 * gb; vb = 0.999 * vb + 0.001 * gb ** 2
        w -= 0.03 * (m / (1 - 0.9 ** t)) / (np.sqrt(v / (1 - 0.999 ** t)) + 1e-8); b -= 0.03 * (mb / (1 - 0.9 ** t)) / (math.sqrt(vb / (1 - 0.999 ** t)) + 1e-8)
    return lambda Z: Z @ w + b


def rank_auc(s, y):
    """area under the ROC curve, ties counted as half"""
    o = np.argsort(s, kind="stable"); r = np.empty(len(s)); r[o] = np.arange(1, len(s) + 1); ss = s[o]; i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and ss[j + 1] == ss[i]: j += 1
        if j > i: r[o[i:j + 1]] = (i + j + 2) / 2
        i = j + 1
    n1 = y.sum(); n0 = len(y) - n1; return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def report(name, score, y, rng, fold):
    """scores from different folds come from different classifiers, so each fold's are centred on its own median first"""
    score = score.copy()
    for k in np.unique(fold): m = fold == k; score[m] = (score[m] - np.median(score[m])) / (score[m].std() + 1e-12)
    a = rank_auc(score, y); n = len(y); boot = [rank_auc(score[i], y[i]) for i in (rng.integers(n, size=n) for _ in range(300))]
    acc = float(((score > 0) == (y == 1)).mean())
    print(f"  {name:34s} AUC {a:.3f} ± {np.std(boot):.3f}   called correctly {100 * acc:.0f}%"); return {"auc": a, "se": float(np.std(boot)), "correct": acc}


def main():
    null = "null" in sys.argv; args = [a for a in sys.argv if a != "null"]
    tag = args[1] if len(args) > 1 else "final"; K = int(args[2]) if len(args) > 2 else 4; seed = int(args[3]) if len(args) > 3 else 1
    d = np.load(os.path.join(ROOT, "zoomtest", "e2e", tag + "_pics.npz")); rng = np.random.default_rng(seed)
    real = np.unpackbits(d["real"], axis=-1).astype(bool); made = np.unpackbits(d["made"], axis=-1).astype(bool); grp = d["group"]; n = len(real)
    pics = np.concatenate([real, made]); y = np.r_[np.zeros(n), np.ones(n)]; g2 = np.r_[grp, grp]
    folds = rng.permutation(grp.max() + 1) % K; fold = folds[g2]
    if null:                                          # a control: real pictures only, called "made" for half of the samples
        half = rng.permutation(grp.max() + 1) % 2; pics = real; y = half[grp].astype(float); fold = folds[grp]; n = len(real) // 2
    print(f"{tag}{' (control: real against real)' if null else ''}: {n} real and {n} made pictures of {real.shape[1]} cells across, from {grp.max() + 1} samples in {K} folds")
    out = {"pictures": int(n), "size": int(real.shape[1]), "folds": K}
    for name, wp in (("summary", True), ("summary, without the pieces", False)):
        X = np.array([summary(p, wp) for p in pics]); s = np.zeros(len(pics))
        for k in range(K):
            tr = fold != k; mu = X[tr].mean(0); sd = X[tr].std(0) + 1e-9; Z = np.clip((X - mu) / sd, -6, 6); s[~tr] = logistic(Z[tr], y[tr])(Z[~tr])
        out[name] = report(name, s, y, rng, fold)
    NP = 3000; s = np.zeros(len(pics)); P = [patches(p, NP, rng) for p in pics]
    for k in range(K):
        tr = np.nonzero(fold != k)[0]; X0 = np.vstack([P[i] for i in tr if y[i] == 0]); X1 = np.vstack([P[i] for i in tr if y[i] == 1])
        f = train_mlp(X0, X1, epochs=6, hidden=64, seed=seed + k)
        for i in np.nonzero(fold == k)[0]: s[i] = f(P[i]).mean()
    out["patches"] = report("patches", s, y, rng, fold)
    if tag == "final" and not null: json.dump(out, open(os.path.join(ROOT, "analysis", "tell_test.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
