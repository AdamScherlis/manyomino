"""Convergence checks for the two-move chain, from the runs of scripts/convergence.sh.

1. Against samplers with no relaxation time.
   exact   percolation clusters of fixed size grown directly (src/leath.cpp): independent exact
           draws from the ensemble the chain samples with -p.  Compared on every observable in
           src/features.h: means, quartiles of each distribution, and a network trained
           to tell the two sets of samples apart, scored on held-out data.
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


def auc(score, y):
    order = np.argsort(score); r = np.empty(len(score)); r[order] = np.arange(1, len(score) + 1)
    n1 = y.sum(); n0 = len(y) - n1
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def train_mlp(X0, X1, epochs=12, hidden=64, seed=0):
    """Train a small neural network (2 hidden layers) to tell rows of X0 from rows of X1.
    Returns the function row -> score (larger = more like X1)."""
    r = np.random.default_rng(seed)
    n = min(len(X0), len(X1)); X0 = X0[r.permutation(len(X0))[:n]]; X1 = X1[r.permutation(len(X1))[:n]]
    X = np.vstack([X0, X1]); y = np.r_[np.zeros(n), np.ones(n)]
    nt = 2 * n; d = X.shape[1]
    W = [r.normal(0, 1 / math.sqrt(d), (d, hidden)), r.normal(0, 1 / math.sqrt(hidden), (hidden, hidden)), r.normal(0, 1 / math.sqrt(hidden), (hidden, 1))]
    B = [np.zeros(hidden), np.zeros(hidden), np.zeros(1)]
    P = W + B; M = [np.zeros_like(a) for a in P]; V = [np.zeros_like(a) for a in P]; t = 0
    def fwd(x):
        h1 = np.maximum(0, x @ W[0] + B[0]); h2 = np.maximum(0, h1 @ W[1] + B[1]); return h1, h2, (h2 @ W[2] + B[2])[:, 0]
    for ep in range(epochs):
        idx = r.permutation(nt)
        for s in range(0, nt, 256):
            b = idx[s:s + 256]; x = X[b]; h1, h2, z = fwd(x)
            g = (1 / (1 + np.exp(-z)) - y[b])[:, None] / len(b)
            gW2 = h2.T @ g; g2 = (g @ W[2].T) * (h2 > 0); gW1 = h1.T @ g2; g1 = (g2 @ W[1].T) * (h1 > 0); gW0 = x.T @ g1
            G = [gW0 + 1e-4 * W[0], gW1 + 1e-4 * W[1], gW2 + 1e-4 * W[2], g1.sum(0), g2.sum(0), g.sum(0)]
            t += 1
            for i in range(6):
                M[i] = 0.9 * M[i] + 0.1 * G[i]; V[i] = 0.999 * V[i] + 0.001 * G[i] ** 2
                P[i] -= 2e-3 * (M[i] / (1 - 0.9 ** t)) / (np.sqrt(V[i] / (1 - 0.999 ** t)) + 1e-8)
    return lambda x: np.tanh(fwd(x)[2] / 4)          # bounded, so one outlier cannot dominate a mean


def witness_test(ref, chains, ref_is_chain=False, max_train=100000, seed=0):
    """A learned two-sample test that uses every chain sample.

    A network is trained to tell the reference samples from the first half of each chain.  Its
    score is then just one more observable -- the one most likely to differ -- and its mean is
    compared on data the network has not seen: the rest of the reference samples against the
    second half of each chain, with the chain's error from batch means.  Returns z (about N(0,1)
    if the two sources have the same distribution) and the held-out AUC."""
    r = np.random.default_rng(seed)
    refs = ref if ref_is_chain else [ref]
    allref = np.vstack(refs)
    mu = allref.mean(0); sd = allref.std(0); keep = sd > 0
    std = lambda X: np.clip((X[:, keep] - mu[keep]) / sd[keep], -8, 8)
    if ref_is_chain:
        Rtr = np.vstack([x[: len(x) // 2] for x in refs]); Rte = [x[len(x) // 2:] for x in refs]
    else:
        p = r.permutation(len(ref)); Rtr = ref[p[: len(ref) // 2]]; Rte = [ref[p[len(ref) // 2:]]]
    Ctr = np.vstack([c[: len(c) // 2] for c in chains]); Cte = [c[len(c) // 2:] for c in chains]
    n = min(len(Rtr), len(Ctr), max_train)
    f = train_mlp(std(Rtr[r.permutation(len(Rtr))[:n]]), std(Ctr[r.permutation(len(Ctr))[:n]]), seed=seed)
    def mean_se(parts, correlated):
        ms, ses = [], []
        for x in parts:
            s = f(std(x)); ms.append(s.mean())
            ses.append(batch_se(s, tau_of(s)) if correlated else s.std(ddof=1) / math.sqrt(len(s)))
        return float(np.mean(ms)), math.sqrt(sum(e * e for e in ses)) / len(ses)
    mr, sr = mean_se(Rte, ref_is_chain); mc, sc = mean_se(Cte, True)
    k = min(sum(len(x) for x in Rte), sum(len(x) for x in Cte), 50000)
    a0 = np.vstack(Rte); a1 = np.vstack(Cte)
    s0 = f(std(a0[r.permutation(len(a0))[:k]])); s1 = f(std(a1[r.permutation(len(a1))[:k]]))
    return {"z": (mc - mr) / math.hypot(sr, sc), "auc": auc(np.r_[s0, s1], np.r_[np.zeros(k), np.ones(k)]), "train": int(2 * n)}


def chain(fn, burn=0.0, max_rows=None):
    """feature rows of one chain; with max_rows, every k-th row so that at most that many are kept
    (the long chains at 250 cells are sampled far more often than their autocorrelation time)"""
    names = open(fn + ".names").read().split()
    a = np.memmap(fn, dtype=np.float64, mode="r")
    a = a[: len(a) // len(names) * len(names)].reshape(-1, len(names))
    a = a[int(len(a) * burn):]
    k = 1 if not max_rows else max(1, -(-len(a) // max_rows))
    return names[1:], np.array(a[::k, 1:])


def chain_moves(fn):
    """number of moves a chain file covers (first column of its last row)"""
    k = len(open(fn + ".names").read().split())
    a = np.memmap(fn, dtype=np.float64, mode="r")
    return float(a[(len(a) // k - 1) * k])


def between_chains(chains, cols=None):
    """Do independent chains agree within their own error bars?  For every observable, chi-square
    of the chain means about their weighted mean, per degree of freedom (1 if the errors are right)."""
    out = []
    for j in range(chains[0].shape[1]) if cols is None else cols:
        xs = [c[:, j] for c in chains]
        if np.concatenate(xs).std() == 0: continue
        t = float(np.mean([tau_of(x) for x in xs])); ms = np.array([x.mean() for x in xs]); ses = np.array([batch_se(x, t) for x in xs])
        w = 1 / ses ** 2; mu = (w * ms).sum() / w.sum(); out.append(float((((ms - mu) / ses) ** 2).sum() / (len(xs) - 1)))
    return out


def compare_exact(tag, exact_fn, chain_fns, label):
    fns = [exact_fn] if isinstance(exact_fn, str) else exact_fn         # one or several runs of leath
    names = load(fns[0], 1)[0]; E = np.vstack([load(f, 1)[2] for f in fns])
    chains = [chain(f, max_rows=400000)[1] for f in chain_fns]
    wt_ = np.array([chain_moves(f) for f in chain_fns], float); wt_ /= wt_.sum()       # chains weighted by their length in moves
    rows = []; taus = []
    for j, nm in enumerate(names):
        if E[:, j].std() == 0 and C[:, j].std() == 0:
            continue
        t = float(np.mean([tau_of(c[:, j]) for c in chains])); taus.append(t)
        se_c = math.sqrt(sum((w * batch_se(c[:, j], t)) ** 2 for w, c in zip(wt_, chains)))
        se_e = E[:, j].std(ddof=1) / math.sqrt(len(E))
        mc = float(sum(w * c[:, j].mean() for w, c in zip(wt_, chains))); me = float(E[:, j].mean())
        rows.append({"name": nm, "exact": me, "exact_se": se_e, "chain": mc, "chain_se": se_c, "z": (mc - me) / math.hypot(se_c, se_e),
                     "sd": float(E[:, j].std()), "tau": t, "j": j})
    # the shape of each distribution: the fraction of samples below each quartile of the exact
    # sample, with the same batch-means errors (so every chain sample is used, correlated or not)
    for r in rows:
        j = r["j"]; r["zq"] = []
        for q in np.unique(np.quantile(E[:, j], [0.25, 0.5, 0.75])):
            fe = float((E[:, j] <= q).mean())
            if not 0.02 < fe < 0.98: continue
            ind = [(c[:, j] <= q).astype(float) for c in chains]; t = float(np.mean([tau_of(x) for x in ind]))
            fc = float(sum(w * x.mean() for w, x in zip(wt_, ind))); se = math.sqrt(sum((w * batch_se(x, t)) ** 2 for w, x in zip(wt_, ind)))
            r["zq"].append((fc - fe) / math.hypot(se, math.sqrt(fe * (1 - fe) / len(E))))
    zq = np.array([z for r in rows for z in r["zq"]])
    wt = witness_test(E, chains); a, z, ntot = wt["auc"], wt["z"], wt["train"]
    zs = np.array([r["z"] for r in rows]); bc = between_chains(chains, [r["j"] for r in rows])
    out = {"tag": tag, "label": label, "chains": len(chains), "between_chains": float(np.mean(bc)) if len(chains) > 1 else None, "n_exact": len(E), "n_chain": int(sum(len(c) for c in chains)), "tau_max": max(taus),
           "features": len(rows), "z_rms": float(np.sqrt((zs ** 2).mean())), "z_max": float(np.abs(zs).max()), "z_over2": int((np.abs(zs) > 2).sum()),
           "z_over3": int((np.abs(zs) > 3).sum()), "q_n": len(zq), "q_rms": float(np.sqrt((zq ** 2).mean())), "q_max": float(np.abs(zq).max()), "q_over2": int((np.abs(zq) > 2).sum()),
           "auc": a, "auc_z": z, "auc_n": ntot, "rows": [{k: v for k, v in r.items() if k != "j"} for r in rows]}
    worst = sorted(rows, key=lambda r: -abs(r["z"]))[:3]
    print(f"{tag}: exact {len(E)}, chain {sum(len(c) for c in chains)}; means z rms {out['z_rms']:.2f} max {out['z_max']:.2f} (|z|>2: {out['z_over2']}/{len(rows)}); "
          f"quartiles z rms {out['q_rms']:.2f} max {out['q_max']:.2f} (|z|>2: {out['q_over2']}/{len(zq)}); learned test z = {z:+.2f} (AUC {a:.4f}); "
          f"{len(chains)} chains" + (f", between-chain chi2/dof {out['between_chains']:.2f}" if len(chains) > 1 else "") + "; largest: " + ", ".join(f"{r['name']} {r['z']:+.2f}" for r in worst))
    return out


LOG_C, LOG_LAMBDA = math.log(0.316915), math.log(4.0625696)


def log_count(m):
    """log(m A(m)) from the known asymptotic form; the 1/m term is fixed by the exact A(24)."""
    return LOG_C + m * LOG_LAMBDA - 0.70 / m


def compare_growth(tag, grow_glob, chain_fns, label, m, burn=0.0):
    """Grown populations against the uniform chain.

    Each population gives a self-normalised estimate of every average and an unbiased estimate of
    the count m A(m).  Populations are combined weighted by their count estimates (the same as
    pooling all walkers with their true weights), which removes most of the 1/N bias of a single
    population; errors are from a jackknife over populations.  The count itself is known, so it is
    a check of the growth sampler: when the genealogy of a population collapses onto too few
    early ancestors, the count estimate drops below the truth and the averages cannot be trusted."""
    files = sorted(f for f in glob.glob(grow_glob) if not f.endswith((".names", ".logz", ".log")))
    reps = []; anc = []; ess = []; lz = []; npop = 0
    for f in files:
        names, lead, X = load(f, 4)
        sel = lead[:, 0] == m; lead = lead[sel]; X = X[sel]
        if len(X) < 100: continue
        w = np.exp(lead[:, 1] - lead[:, 1].max()); npop = len(X)
        reps.append((w[:, None] * X).sum(0) / w.sum()); anc.append((len(np.unique(lead[:, 2])), len(np.unique(lead[:, 3])))); ess.append(w.sum() ** 2 / (w ** 2).sum())
        lz.append(next(float(l.split()[1]) for l in open(f + ".logz") if int(float(l.split()[0])) == m))
    reps = np.array(reps); R = len(reps); lz = np.array(lz); zw = np.exp(lz - lz.max())
    comb = lambda w, M: (w[:, None] * M).sum(0) / w.sum()
    est = comb(zw, reps)
    jk = np.array([comb(np.delete(zw, i), np.delete(reps, i, 0)) for i in range(R)])
    est_se = np.sqrt((R - 1) / R * ((jk - jk.mean(0)) ** 2).sum(0))
    lmz = lz.max() + math.log(zw.mean())
    jz = np.array([lz.max() + math.log(np.delete(zw, i).mean()) for i in range(R)])
    lmz_se = float(np.sqrt((R - 1) / R * ((jz - jz.mean()) ** 2).sum()))
    chains = [chain(f, burn)[1] for f in chain_fns]
    rows = []
    for j, nm in enumerate(names):
        if reps[:, j].std() == 0:
            continue
        t = float(np.mean([tau_of(c[:, j]) for c in chains]))
        se_c = math.sqrt(sum(batch_se(c[:, j], t) ** 2 for c in chains)) / len(chains)
        mc = float(np.mean([c[:, j].mean() for c in chains])); mg = float(est[j]); se_g = float(est_se[j])
        rows.append({"name": nm, "growth": mg, "growth_se": se_g, "chain": mc, "chain_se": se_c, "t": (mc - mg) / math.hypot(se_c, se_g),
                     "sd": float(np.vstack(chains)[:, j].std())})
    ts = np.array([r["t"] for r in rows])
    out = {"tag": tag, "label": label, "n": m, "replicas": R, "population": int(npop), "chain_samples": int(sum(len(c) for c in chains)),
           "median_precision": float(np.median([r["growth_se"] / r["sd"] for r in rows])), "ancestors_m10": [a[0] for a in anc], "ancestors_m2": [a[1] for a in anc], "ess": [float(e) for e in ess],
           "features": len(rows), "t_rms": float(np.sqrt((ts ** 2).mean())), "t_max": float(np.abs(ts).max()), "t_over2": int((np.abs(ts) > 2).sum()), "rows": rows,
           "rel_precision_rg2": next(r["growth_se"] / r["growth"] for r in rows if r["name"] == "rg2"),
           "count_dev": float(lmz - log_count(m)), "count_dev_se": lmz_se, "logz_sd": float(lz.std(ddof=1))}
    worst = sorted(rows, key=lambda r: -abs(r["t"]))[:3]
    print(f"{tag}: {R} populations of {npop}, distinct ancestors at n/2: {int(np.mean(out['ancestors_m2']))}, at n/10: {int(np.mean(out['ancestors_m10']))}; "
          f"count check {100 * (math.exp(out['count_dev']) - 1):+.1f}% ± {100 * lmz_se:.1f}%; precision {out['median_precision']:.3f} sd; "
          f"t rms {out['t_rms']:.2f} max {out['t_max']:.2f} (|t|>2: {out['t_over2']}/{len(rows)}); largest: " + ", ".join(f"{r['name']} {r['t']:+.2f}" for r in worst))
    return out


ODD = lambda nm: nm in ("q2c", "q2s") or (len(nm) == 3 and nm[0] == "m" and nm[1:].isdigit() and (int(nm[1]) % 2 or int(nm[2]) % 2))


def symmetry_test(n, chain_fns, burn=0.0):
    """Observables that change sign under a reflection or a quarter turn average to exactly zero."""
    chains = [chain(f, burn) for f in chain_fns]; names = chains[0][0]; rows = []
    for j, nm in enumerate(names):
        if not ODD(nm): continue
        xs = [c[1][:, j] for c in chains]
        t = float(np.mean([tau_of(x) for x in xs])); se = math.sqrt(sum(batch_se(x, t) ** 2 for x in xs)) / len(xs)
        rows.append({"name": nm, "z": float(np.mean([x.mean() for x in xs]) / se)})
    z = np.array([r["z"] for r in rows])
    out = {"n": n, "features": len(rows), "z_rms": float(np.sqrt((z ** 2).mean())), "z_max": float(np.abs(z).max()), "rows": rows}
    print(f"symmetry n={n}: {len(rows)} observables that must average to zero: z rms {out['z_rms']:.2f} max {out['z_max']:.2f}")
    return out


def start_test(n, start_fns, ref_fns, burn_moves_per_cell):
    """Chains started from a straight bar and from a square, after burn-in, against the chains
    that came up the size ladder."""
    S = []
    for f in start_fns:
        names, X = chain(f); S.append(X[int(burn_moves_per_cell):])        # one sample per n moves
    refs = [chain(f, 0.05)[1] for f in ref_fns]
    rows = []
    for j, nm in enumerate(names):
        if np.vstack(refs)[:, j].std() == 0: continue
        def ms(parts):
            t = float(np.mean([tau_of(x[:, j]) for x in parts]))
            return float(np.mean([x[:, j].mean() for x in parts])), math.sqrt(sum(batch_se(x[:, j], t) ** 2 for x in parts)) / len(parts)
        a, ea = ms(S); b, eb = ms(refs)
        rows.append({"name": nm, "z": (a - b) / math.hypot(ea, eb)})
    z = np.array([r["z"] for r in rows]); wt = witness_test(refs, S, ref_is_chain=True)
    # how long the start is remembered: Rg^2 of the bar-started chain against its final average
    names, X = chain(start_fns[0]); rg = X[:, names.index("rg2")]; fin = rg[len(rg) // 2:].mean()
    k = np.arange(len(rg)); blocks = np.unique(np.round(np.geomspace(1, len(rg), 60)).astype(int))
    relax = [[float(b0), float(rg[b0 - 1:b1].mean() / fin)] for b0, b1 in zip(blocks[:-1], blocks[1:])]
    out = {"n": n, "burn": burn_moves_per_cell, "samples": [len(x) for x in S], "features": len(rows), "z_rms": float(np.sqrt((z ** 2).mean())),
           "z_max": float(np.abs(z).max()), "z_over2": int((np.abs(z) > 2).sum()), "witness_z": wt["z"], "auc": wt["auc"], "rows": rows, "relax": relax}
    worst = sorted(rows, key=lambda r: -abs(r["z"]))[:3]
    print(f"starts n={n}: bar and square starts after {burn_moves_per_cell} moves per cell against ladder chains: z rms {out['z_rms']:.2f} max {out['z_max']:.2f} "
          f"(|z|>2: {out['z_over2']}/{len(rows)}); learned test z = {wt['z']:+.2f}; largest: " + ", ".join(f"{r['name']} {r['z']:+.2f}" for r in worst))
    return out


def compare_enum(n, p=""):
    """the chain against exact averages over every polyomino of n cells (src/enumerate.cpp);
    with p, both weighted by (1-p)^perimeter"""
    ex = {}; sfx = f"_p{p}" if p else ""
    lines = open(os.path.join(CONV, f"exact_{n}{sfx}.txt")).read().split("\n")
    count = int(lines[0].split()[1])
    for ln in lines[1:]:
        if ln.strip(): k, m, sd = ln.split(); ex[k] = (float(m), float(sd))
    names, X = chain(os.path.join(CONV, f"e{n}{sfx}.feat"))
    rows = []
    for j, nm in enumerate(names):
        m, sd = ex[nm]
        if sd == 0: continue
        x = X[:, j]; t = tau_of(x[:400000]); se = batch_se(x, t)
        x2 = (x - m) ** 2; se2 = batch_se(x2, tau_of(x2[:400000]))            # second moment about the exact mean
        rows.append({"name": nm, "exact": m, "exact_sd": sd, "chain": float(x.mean()), "chain_se": se, "z": (x.mean() - m) / se,
                     "z_var": (x2.mean() - sd * sd) / se2, "tau": t})
    z = np.array([r["z"] for r in rows]); zv = np.array([r["z_var"] for r in rows])
    out = {"n": n, "p": float(p) if p else 0.0, "count": count, "samples": len(X), "features": len(rows), "z_rms": float(np.sqrt((z ** 2).mean())), "z_max": float(np.abs(z).max()),
           "z_over2": int((np.abs(z) > 2).sum()), "zvar_rms": float(np.sqrt((zv ** 2).mean())), "zvar_max": float(np.abs(zv).max()),
           "rel_precision": float(np.median([r["chain_se"] / r["exact_sd"] for r in rows])), "rows": rows}
    worst = sorted(rows, key=lambda r: -abs(r["z"]))[:3]
    print(f"enumeration n={n}{' p=' + p if p else ''}: {count} shapes, {len(X)} chain samples, {len(rows)} observables; means z rms {out['z_rms']:.2f} max {out['z_max']:.2f} "
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


def fit_power(ns, ys, es):
    """weighted least squares of log y on log n: (exponent, its error, prefactor)"""
    x = np.log(np.array(ns, float)); y = np.log(np.array(ys, float)); w = (np.array(ys, float) / np.array(es, float)) ** 2
    xm = (w * x).sum() / w.sum(); ym = (w * y).sum() / w.sum(); sxx = (w * (x - xm) ** 2).sum()
    b = (w * (x - xm) * (y - ym)).sum() / sxx
    chi2 = (w * (y - ym - b * (x - xm)) ** 2).sum(); dof = max(1, len(x) - 2)
    return float(b), float(math.sqrt(1 / sxx) * max(1.0, math.sqrt(chi2 / dof))), float(math.exp(ym - b * xm))


def main():
    R = {"enum": [], "exact": [], "growth": [], "growth_small": [], "slow": [], "ab": [], "symmetry": [], "starts": []}
    for n, p in ((14, ""), (16, ""), (14, "0.5"), (14, "0.592746")):
        sfx = f"_p{p}" if p else ""
        if os.path.exists(os.path.join(CONV, f"e{n}{sfx}.feat")) and os.path.getsize(os.path.join(CONV, f"exact_{n}{sfx}.txt")) > 100:
            R["enum"].append(compare_enum(n, p))
    ex = lambda p, n: sorted(glob.glob(os.path.join(CONV, f"x{p}_{n}.feat")) + glob.glob(os.path.join(CONV, f"x{p}_{n}_s?.feat")))
    tl = lambda p, n: sorted(glob.glob(os.path.join(CONV, f"t{p}_{n}_?.feat")) + glob.glob(os.path.join(CONV, f"r{p}_{n}_?.feat")))
    quick = bool(os.environ.get("CONV_QUICK"))           # a fast partial pass, for testing the pipeline
    for p, n in (("0.55", 500),) if quick else (("0.5", 60), ("0.5", 250), ("0.55", 500), ("0.592746", 1000)):
        if ex(p, n) and len(tl(p, n)) >= 1:
            R["exact"].append(compare_exact(f"p={p} n={n}", ex(p, n), tl(p, n), f"{n} cells at density {float(p):.3g}"))
    ctl = os.path.join(CONV, "t0.51_250_0.feat")
    if os.path.exists(ctl) and ex("0.5", 250):
        R["control"] = compare_exact("control: chain at p=0.51 against exact p=0.5, n=250", ex("0.5", 250), [ctl], "control")
    uf = lambda n, c: os.path.join(CONV, f"u{n}_{c}.feat")
    have = lambda g: len(glob.glob(os.path.join(CONV, g + ".logz"))) >= 3
    for m, g, burn in ((250, "gk1000_?", 0.0), (500, "gk1000_?", 0.05), (1000, "gk1000_?", 0.0), (2000, "gk2000_?", 0.05)):
        if have(g) and os.path.exists(uf(m, "A")) and os.path.exists(uf(m, "B")):
            R["growth"].append(compare_growth(f"uniform n={m}", os.path.join(CONV, g), [uf(m, "A"), uf(m, "B")], f"{m} cells", m, burn))
    # the same comparison with smaller populations: how the growth sampler degrades
    for m, g in ((1000, "gs1000a_*"), (1000, "gs1000b_*"), (2000, "gs2000a_*")):
        if have(g):
            R["growth_small"].append(compare_growth(f"uniform n={m} ({g})", os.path.join(CONV, g), [uf(m, "A"), uf(m, "B")], f"{m} cells", m, 0.05))
    for n in (1000, 4000, 16000) if quick else (250, 1000, 4000, 16000, 64000, 256000):
        if os.path.exists(uf(n, "A")) and os.path.exists(uf(n, "B")):
            s = slow_modes(n, uf(n, "A"), uf(n, "B")); R["slow"].append(s)
            # two independent chains: can a network trained on their first halves tell their second halves apart?
            wt = witness_test([chain(uf(n, "A"))[1]], [chain(uf(n, "B"))[1]], ref_is_chain=True)
            R["ab"].append({"n": n, "auc": wt["auc"], "z": wt["z"]})
            print(f"   chain A against chain B: learned test z = {wt['z']:+.2f} (AUC {wt['auc']:.4f})")
            R["symmetry"].append(symmetry_test(n, [uf(n, "A"), uf(n, "B")]))
    # eight chains from the same equilibrated shape with different random numbers
    mf = sorted(glob.glob(os.path.join(CONV, "m250_?.feat")))
    if len(mf) >= 4:
        chs = [chain(f)[1] for f in mf]; bc = np.array(between_chains(chs))
        R["multi"] = {"n": 250, "chains": len(chs), "moves": float(sum(chain_moves(f) for f in mf)), "features": len(bc), "chi2_mean": float(bc.mean()), "chi2_max": float(bc.max())}
        print(f"{len(chs)} chains at n=250 from one start: between-chain chi2/dof over {len(bc)} observables: mean {bc.mean():.2f}, max {bc.max():.2f} ({len(chs) - 1} dof each)")
    # the largest symmetry violation, measured again on chains that were not used to find it
    if R["symmetry"]:
        zmax, nm, n = max((abs(r["z"]), r["name"], s["n"]) for s in R["symmetry"] for r in s["rows"])
        again = [uf(n, c) for c in ("C", "D") if os.path.exists(uf(n, c))]
        if again:
            def mean_se(fns):
                xs = []
                for f in fns:
                    names, X = chain(f); xs.append(X[:, names.index(nm)])
                sd = np.concatenate(xs).std(); t = float(np.mean([tau_of(x) for x in xs])); w = np.array([len(x) for x in xs], float); w /= w.sum()
                return float(sum(wi * x.mean() for wi, x in zip(w, xs)) / sd), float(math.sqrt(sum((wi * batch_se(x, t)) ** 2 for wi, x in zip(w, xs))) / sd)
            m1, e1 = mean_se([uf(n, "A"), uf(n, "B")]); m2, e2 = mean_se(again)
            R["retest"] = {"n": n, "name": nm, "first": m1, "first_se": e1, "second": m2, "second_se": e2, "z_second": m2 / e2}
            print(f"largest symmetry violation: {nm} at n={n}: {m1:+.3f} ± {e1:.3f} sd ({m1 / e1:+.2f}); on {len(again)} fresh chains {m2:+.3f} ± {e2:.3f} sd ({m2 / e2:+.2f})")
    for n, burn in ((4000, 2000), (16000, 2000)):
        st = [os.path.join(CONV, f"s{n}_{sh}.feat") for sh in ("bar", "rect")]
        if all(os.path.exists(f) for f in st):
            R["starts"].append(start_test(n, st, [uf(n, "A"), uf(n, "B")], burn))
    # how the autocorrelation times grow with size (sizes from 1,000 cells up)
    # (a size whose runs are too short to time its slowest observable to 50% is left out)
    big = [s for s in R["slow"] if s["n"] >= 1000 and s["slowest_single"][0]["err"] < 0.5 * s["slowest_single"][0]["tau"]]
    if len(big) >= 3:
        ns = [s["n"] for s in big]
        R["scaling"] = {}
        for key, get in (("rg2", lambda s: next((x["tau"], x["err"]) for x in s["single"] if x["name"] == "rg2")),
                         ("single", lambda s: (s["slowest_single"][0]["tau"], s["slowest_single"][0]["err"])),
                         ("mode", lambda s: (s["slowest_mode"]["tau"], s["slowest_mode"]["err"]))):
            v = [get(s) for s in big]; b, be, a = fit_power(ns, [x[0] for x in v], [x[1] for x in v])
            R["scaling"][key] = {"exponent": b, "err": be, "at_1M": a * 1024000 ** b}
            print(f"tau({key}) per cell ~ n^{b:+.3f} ± {be:.3f}; extrapolated to 1,024,000 cells: {a * 1024000 ** b:.0f} moves per cell")
    json.dump(R, open(os.path.join(ROOT, "analysis", "convergence.json"), "w"), indent=1)
    return R


if __name__ == "__main__":
    main()
