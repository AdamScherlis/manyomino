"""Convergence checks for the two-move chain, from the runs of scripts/convergence.sh.

1. Against samplers with no relaxation time.
   exact   percolation clusters of fixed size grown directly (src/leath.cpp): independent exact
           draws from the ensemble the chain samples with -p.  Compared on every observable in
           src/features.h: means (z-scores), whole distributions (Kolmogorov-Smirnov), and a
           classifier trained to tell the two sets of samples apart.
   growth  uniform polyominoes from a resampled population of growing clusters (src/grow.cpp).
           Correct weights but correlated samples, so only means are compared, with errors
           from the scatter between independent populations.
2. A search for slow modes of the uniform chain.  Each observable's autocorrelation time, and
   the slowest-decaying combinations of observables (and their squares) found by time-lagged
   component analysis on one chain and timed on the other, independent one.  Any observable's
   autocorrelation time is a lower bound on the chain's relaxation time; the search raises it.

  python3 analysis/convergence.py            # writes analysis/convergence.json and .md
"""
import glob, json, math, os, re, sys
import numpy as np
from tau import tau_int, acf_fft

ROOT = os.path.join(os.path.dirname(__file__), "..")
CONV = os.path.join(ROOT, "conv")
rng = np.random.default_rng(1)


def load(fn, lead):
    names = open(fn + ".names").read().split()
    a = np.fromfile(fn)
    a = a[: len(a) // len(names) * len(names)].reshape(-1, len(names))
    return names[lead:], a[:, :lead], a[:, lead:]


def tau_of(x):
    if np.std(x) == 0:
        return 0.5
    return max(0.5, tau_int(x)[0])


def batch_se(x, tau):
    """standard error of the mean of a correlated series, from batches at least 40 tau long"""
    nb = int(min(50, max(8, len(x) / (40 * tau))))
    m = len(x) // nb
    b = x[: m * nb].reshape(nb, m).mean(1)
    return b.std(ddof=1) / math.sqrt(nb)


def ks_p(a, b):
    """two-sample Kolmogorov-Smirnov p-value (asymptotic)"""
    a = np.sort(a); b = np.sort(b); allv = np.concatenate([a, b])
    d = np.abs(np.searchsorted(a, allv, side="right") / len(a) - np.searchsorted(b, allv, side="right") / len(b)).max()
    en = math.sqrt(len(a) * len(b) / (len(a) + len(b)))
    lam = (en + 0.12 + 0.11 / en) * d
    k = np.arange(1, 101)
    return float(min(1.0, max(0.0, 2 * np.sum((-1) ** (k - 1) * np.exp(-2 * k * k * lam * lam)))))


def auc(score, y):
    order = np.argsort(score); r = np.empty(len(score)); r[order] = np.arange(1, len(score) + 1)
    n1 = y.sum(); n0 = len(y) - n1
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def classifier_test(X0, X1, epochs=40, hidden=64, seed=0):
    """Train a small neural network to tell X0 from X1; return held-out AUC and its z-score
    against 0.5 (the value when the two sets come from the same distribution)."""
    r = np.random.default_rng(seed)
    n = min(len(X0), len(X1)); X0 = X0[r.permutation(len(X0))[:n]]; X1 = X1[r.permutation(len(X1))[:n]]
    X = np.vstack([X0, X1]); y = np.r_[np.zeros(n), np.ones(n)]
    sd = X.std(0); keep = sd > 0; X = (X[:, keep] - X[:, keep].mean(0)) / sd[keep]
    X = np.clip(X, -8, 8)
    p = r.permutation(2 * n); X, y = X[p], y[p]; nt = int(0.7 * 2 * n)
    Xt, yt, Xv, yv = X[:nt], y[:nt], X[nt:], y[nt:]
    d = X.shape[1]
    W = [r.normal(0, 1 / math.sqrt(d), (d, hidden)), r.normal(0, 1 / math.sqrt(hidden), (hidden, hidden)), r.normal(0, 1 / math.sqrt(hidden), (hidden, 1))]
    B = [np.zeros(hidden), np.zeros(hidden), np.zeros(1)]
    P = W + B; M = [np.zeros_like(a) for a in P]; V = [np.zeros_like(a) for a in P]; t = 0
    def fwd(x):
        h1 = np.maximum(0, x @ W[0] + B[0]); h2 = np.maximum(0, h1 @ W[1] + B[1]); return h1, h2, (h2 @ W[2] + B[2])[:, 0]
    for ep in range(epochs):
        idx = r.permutation(nt)
        for s in range(0, nt, 256):
            b = idx[s:s + 256]; x = Xt[b]; h1, h2, z = fwd(x)
            g = (1 / (1 + np.exp(-z)) - yt[b])[:, None] / len(b)
            gW2 = h2.T @ g; g2 = (g @ W[2].T) * (h2 > 0); gW1 = h1.T @ g2; g1 = (g2 @ W[1].T) * (h1 > 0); gW0 = x.T @ g1
            G = [gW0 + 1e-4 * W[0], gW1 + 1e-4 * W[1], gW2 + 1e-4 * W[2], g1.sum(0), g2.sum(0), g.sum(0)]
            t += 1
            for i in range(6):
                M[i] = 0.9 * M[i] + 0.1 * G[i]; V[i] = 0.999 * V[i] + 0.001 * G[i] ** 2
                P[i] -= 2e-3 * (M[i] / (1 - 0.9 ** t)) / (np.sqrt(V[i] / (1 - 0.999 ** t)) + 1e-8)
    a = auc(fwd(Xv)[2], yv); n1 = yv.sum(); n0 = len(yv) - n1
    return a, (a - 0.5) / math.sqrt((n1 + n0 + 1) / (12 * n1 * n0)), int(2 * n)


def chain(fn, burn=0.0):
    names, lead, X = load(fn, 1)
    k = int(len(X) * burn)
    return names, X[k:]


def compare_exact(tag, exact_fn, chain_fns, label):
    names, _, E = load(exact_fn, 1)
    chains = [chain(f)[1] for f in chain_fns]
    C = np.vstack(chains)
    rows = []; taus = []
    for j, nm in enumerate(names):
        if E[:, j].std() == 0 and C[:, j].std() == 0:
            continue
        t = float(np.mean([tau_of(c[:, j]) for c in chains])); taus.append(t)
        se_c = math.sqrt(sum(batch_se(c[:, j], t) ** 2 for c in chains)) / len(chains)
        se_e = E[:, j].std(ddof=1) / math.sqrt(len(E))
        mc = float(np.mean([c[:, j].mean() for c in chains])); me = float(E[:, j].mean())
        rows.append({"name": nm, "exact": me, "exact_se": se_e, "chain": mc, "chain_se": se_c, "z": (mc - me) / math.hypot(se_c, se_e),
                     "sd": float(E[:, j].std()), "tau": t, "j": j})
    thin = int(math.ceil(2 * max(taus)))
    T = np.vstack([c[rng.integers(thin)::thin] for c in chains])       # nearly independent chain samples
    for r in rows:
        r["ks_p"] = ks_p(E[:, r["j"]], T[:, r["j"]])
    cols = [r["j"] for r in rows]
    a, z, ntot = classifier_test(E[:, cols], T[:, cols])
    zs = np.array([r["z"] for r in rows]); ps = np.array([r["ks_p"] for r in rows])
    out = {"tag": tag, "label": label, "n_exact": len(E), "n_chain": len(C), "n_chain_thinned": len(T), "thin": thin, "tau_max": max(taus),
           "features": len(rows), "z_rms": float(np.sqrt((zs ** 2).mean())), "z_max": float(np.abs(zs).max()), "z_over2": int((np.abs(zs) > 2).sum()),
           "z_over3": int((np.abs(zs) > 3).sum()), "ks_min_p": float(ps.min()), "ks_below_05": int((ps < 0.05).sum()),
           "auc": a, "auc_z": z, "auc_n": ntot, "rows": [{k: v for k, v in r.items() if k != "j"} for r in rows]}
    worst = sorted(rows, key=lambda r: -abs(r["z"]))[:3]
    print(f"{tag}: exact {len(E)}, chain {len(C)} (thinned {len(T)}, 2 tau_max = {thin} samples); z rms {out['z_rms']:.2f} max {out['z_max']:.2f} "
          f"(|z|>2: {out['z_over2']}/{len(rows)}); KS min p {out['ks_min_p']:.3g} (<0.05: {out['ks_below_05']}); classifier AUC {a:.4f} (z = {z:+.2f}); "
          f"largest: " + ", ".join(f"{r['name']} {r['z']:+.2f}" for r in worst))
    return out


def compare_growth(tag, grow_glob, chain_fns, label, m, burn=0.0):
    files = sorted(f for f in glob.glob(grow_glob) if not f.endswith((".names", ".logz", ".log")))
    reps = []; anc = []; ess = []
    for f in files:
        names, lead, X = load(f, 4)
        sel = lead[:, 0] == m; lead = lead[sel]; X = X[sel]
        if not len(X): continue
        w = np.exp(lead[:, 1] - lead[:, 1].max())
        reps.append((w[:, None] * X).sum(0) / w.sum()); anc.append((len(np.unique(lead[:, 2])), len(np.unique(lead[:, 3])))); ess.append(w.sum() ** 2 / (w ** 2).sum())
    reps = np.array(reps); R = len(reps)
    chains = [chain(f, burn)[1] for f in chain_fns]
    rows = []
    for j, nm in enumerate(names):
        if reps[:, j].std() == 0:
            continue
        t = float(np.mean([tau_of(c[:, j]) for c in chains]))
        se_c = math.sqrt(sum(batch_se(c[:, j], t) ** 2 for c in chains)) / len(chains)
        mc = float(np.mean([c[:, j].mean() for c in chains])); mg = float(reps[:, j].mean()); se_g = float(reps[:, j].std(ddof=1) / math.sqrt(R))
        rows.append({"name": nm, "growth": mg, "growth_se": se_g, "chain": mc, "chain_se": se_c, "t": (mc - mg) / math.hypot(se_c, se_g),
                     "sd": float(np.vstack(chains)[:, j].std())})
    ts = np.array([r["t"] for r in rows])
    out = {"tag": tag, "label": label, "n": m, "replicas": R, "population": int(len(X)), "chain_samples": int(sum(len(c) for c in chains)),
           "median_precision": float(np.median([r["growth_se"] / r["sd"] for r in rows])), "ancestors_m10": [a[0] for a in anc], "ancestors_m2": [a[1] for a in anc], "ess": [float(e) for e in ess],
           "features": len(rows), "t_rms": float(np.sqrt((ts ** 2).mean())), "t_max": float(np.abs(ts).max()), "rows": rows,
           "rel_precision_rg2": next(r["growth_se"] / r["growth"] for r in rows if r["name"] == "rg2")}
    worst = sorted(rows, key=lambda r: -abs(r["t"]))[:3]
    print(f"{tag}: {R} populations, distinct ancestors at n/2: {out['ancestors_m2']}; precision {out['median_precision']:.3f} sd; t rms {out['t_rms']:.2f} max {out['t_max']:.2f} (t with {R - 1} dof); "
          f"largest: " + ", ".join(f"{r['name']} {r['t']:+.2f}" for r in worst))
    return out


def compare_enum(n):
    """the uniform chain against exact averages over every polyomino of n cells (src/enumerate.cpp)"""
    ex = {}
    lines = open(os.path.join(CONV, f"exact_{n}.txt")).read().split("\n")
    count = int(lines[0].split()[1])
    for ln in lines[1:]:
        if ln.strip(): k, m, sd = ln.split(); ex[k] = (float(m), float(sd))
    names, X = chain(os.path.join(CONV, f"e{n}.feat"))
    rows = []
    for j, nm in enumerate(names):
        m, sd = ex[nm]
        if sd == 0: continue
        x = X[:, j]; t = tau_of(x[:400000]); se = batch_se(x, t)
        x2 = (x - m) ** 2; se2 = batch_se(x2, tau_of(x2[:400000]))            # second moment about the exact mean
        rows.append({"name": nm, "exact": m, "exact_sd": sd, "chain": float(x.mean()), "chain_se": se, "z": (x.mean() - m) / se,
                     "z_var": (x2.mean() - sd * sd) / se2, "tau": t})
    z = np.array([r["z"] for r in rows]); zv = np.array([r["z_var"] for r in rows])
    out = {"n": n, "count": count, "samples": len(X), "features": len(rows), "z_rms": float(np.sqrt((z ** 2).mean())), "z_max": float(np.abs(z).max()),
           "z_over2": int((np.abs(z) > 2).sum()), "zvar_rms": float(np.sqrt((zv ** 2).mean())), "zvar_max": float(np.abs(zv).max()),
           "rel_precision": float(np.median([r["chain_se"] / r["exact_sd"] for r in rows])), "rows": rows}
    worst = sorted(rows, key=lambda r: -abs(r["z"]))[:3]
    print(f"enumeration n={n}: {count} shapes, {len(X)} chain samples, {len(rows)} observables; means z rms {out['z_rms']:.2f} max {out['z_max']:.2f} "
          f"(|z|>2: {out['z_over2']}); variances z rms {out['zvar_rms']:.2f} max {out['zvar_max']:.2f}; median precision {out['rel_precision']:.1e} sd; largest: "
          + ", ".join(f"{r['name']} {r['z']:+.2f}" for r in worst))
    return out


def tica(X, lag, eps=1e-8):
    """slowest linear combinations of the columns of X at the given lag: (projection matrix, eigenvalues)"""
    C0 = X.T @ X / len(X)
    w, U = np.linalg.eigh(C0); keep = w > eps * w.max()
    Wh = U[:, keep] / np.sqrt(w[keep])
    Y = X @ Wh
    Ct = Y[:-lag].T @ Y[lag:] / (len(Y) - lag); Ct = (Ct + Ct.T) / 2
    lam, V = np.linalg.eigh(Ct)
    o = np.argsort(-lam)
    return Wh @ V[:, o], lam[o]


def slow_modes(n, fa, fb):
    names, A = chain(fa); _, B = chain(fb)
    T = min(len(A), len(B))
    P = np.vstack([A, B]); mu = P.mean(0); sd = P.std(0); keep = sd > 0
    names = [nm for nm, k in zip(names, keep) if k]
    ZA = (A[:, keep] - mu[keep]) / sd[keep]; ZB = (B[:, keep] - mu[keep]) / sd[keep]
    # single observables (sample spacing is n moves, so times are in moves per cell)
    single = []
    for j, nm in enumerate(names):
        ta, tb = tau_int(ZA[:, j]), tau_int(ZB[:, j])
        single.append({"name": nm, "tau": float((ta[0] + tb[0]) / 2), "err": float(math.hypot(ta[2], tb[2]) / 2)})
    tau_rg = next(s["tau"] for s in single if s["name"] == "rg2")
    # combinations: observables and their squares, fitted on one chain and timed on the other
    def expand(Z):
        Q = Z ** 2; Q = (Q - Q.mean(0)) / (Q.std(0) + 1e-12)
        return np.hstack([Z, Q])
    EA, EB = expand(ZA), expand(ZB)
    mA, sA = EA.mean(0), EA.std(0) + 1e-12
    EA = (EA - mA) / sA; EB = (EB - mA) / sA
    modes = []
    for lagf in (0.5, 2.0):
        lag = max(1, int(round(lagf * tau_rg)))
        for (Ftrain, Ftest, which) in ((EA, EB, "A->B"), (EB, EA, "B->A")):
            Wm, lam = tica(Ftrain, lag)
            for k in range(3):
                y = Ftest @ Wm[:, k]; t = tau_int(y)
                rho = acf_fft(y)
                texp = -lag / math.log(rho[lag]) if 0 < rho[lag] < 1 else float("nan")
                load_ = np.array([abs(np.corrcoef(y, Ftest[:, j])[0, 1]) for j in range(Ftest.shape[1])])
                top = np.argsort(-load_)[:4]
                lab = [(names[j] if j < len(names) else names[j - len(names)] + "^2") + f" {load_[j]:.2f}" for j in top]
                modes.append({"lag": lag, "split": which, "rank": k, "train_eig": float(lam[k]), "tau": float(t[0]), "err": float(t[2]), "tau_exp": float(texp), "loadings": lab})
    best = max(modes, key=lambda m: m["tau"])
    slow_single = sorted(single, key=lambda s: -s["tau"])[:6]
    out = {"n": n, "samples": [len(A), len(B)], "tau_rg2": tau_rg, "run_over_tau_rg2": float((len(A) + len(B)) / tau_rg), "single": single,
           "slowest_single": slow_single, "modes": modes, "slowest_mode": best,
           "ratio_single": slow_single[0]["tau"] / tau_rg, "ratio_mode": best["tau"] / tau_rg}
    print(f"n={n}: {len(A)}+{len(B)} samples; tau(Rg2) = {tau_rg:.1f} moves/cell; slowest single: " +
          ", ".join(f"{s['name']} {s['tau']:.1f}" for s in slow_single[:4]) +
          f"; slowest combination (held out): {best['tau']:.1f} ± {best['err']:.1f} [{'; '.join(best['loadings'][:3])}]")
    return out


def main():
    R = {"enum": [], "exact": [], "growth": [], "slow": [], "ab": []}
    for n in (14, 16):
        if os.path.exists(os.path.join(CONV, f"e{n}.feat")) and os.path.getsize(os.path.join(CONV, f"exact_{n}.txt")) > 100:
            R["enum"].append(compare_enum(n))
    ex = lambda p, n: os.path.join(CONV, f"x{p}_{n}.feat")
    tl = lambda p, n: [os.path.join(CONV, f"t{p}_{n}_{k}.feat") for k in (0, 1)]
    for p, n in (("0.5", 250), ("0.55", 500), ("0.592746", 1000)):
        if os.path.exists(ex(p, n)) and all(os.path.exists(f) for f in tl(p, n)):
            R["exact"].append(compare_exact(f"p={p} n={n}", ex(p, n), tl(p, n), f"{n} cells at density {float(p):.3g}"))
    ctl = os.path.join(CONV, "t0.51_250_0.feat")
    if os.path.exists(ctl) and os.path.exists(ex("0.5", 250)):
        R["control"] = compare_exact("control: chain at p=0.51 against exact p=0.5, n=250", ex("0.5", 250), [ctl], "control")
    uf = lambda n, c: os.path.join(CONV, f"u{n}_{c}.feat")
    for m, g, burn in ((250, "gk1000_?", 0.0), (500, "gk1000_?", 0.05), (1000, "gk1000_?", 0.0), (2000, "gk2000_?", 0.05)):
        if len(glob.glob(os.path.join(CONV, g + ".logz"))) >= 3 and os.path.exists(uf(m, "A")) and os.path.exists(uf(m, "B")):
            R["growth"].append(compare_growth(f"uniform n={m}", os.path.join(CONV, g), [uf(m, "A"), uf(m, "B")], f"{m} cells", m, burn))
    for n in (250, 1000, 4000, 16000, 64000):
        if os.path.exists(uf(n, "A")) and os.path.exists(uf(n, "B")):
            s = slow_modes(n, uf(n, "A"), uf(n, "B")); R["slow"].append(s)
            # two independent chains, thinned to independence: can a classifier tell them apart?
            names, A = chain(uf(n, "A")); _, B = chain(uf(n, "B"))
            thin = int(math.ceil(2 * max(x["tau"] for x in s["single"])))
            a, z, ntot = classifier_test(A[rng.integers(thin)::thin], B[rng.integers(thin)::thin])
            R["ab"].append({"n": n, "auc": a, "z": z, "samples": ntot, "thin": thin})
            print(f"   chain A against chain B, {ntot} thinned samples: classifier AUC {a:.4f} (z = {z:+.2f})")
    json.dump(R, open(os.path.join(ROOT, "analysis", "convergence.json"), "w"), indent=1)
    return R


if __name__ == "__main__":
    main()
