"""Reduce the production time series to one JSON file for the page."""
import glob, json, re, sys, os
import numpy as np
from tau import tau_int

ROOT = os.path.join(os.path.dirname(__file__), "..")
COLS = "move sec rg2 asph perim cycles noncut leaves w h holes deg1 deg2 deg3 deg4".split()
BURN_N = 300.0          # minimum burn-in of a ladder stage, in units of n moves
BATCH_N = 70.0          # minimum length of an error-bar batch, in units of n moves
LAMBDA = 4.0625696   # Jensen's series estimate of Klarner's constant


def read(fn):
    rows = []
    with open(fn) as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.split()
            if len(p) == len(COLS):
                rows.append(p)
    return np.array(rows, float)


def batches(vals, nb):
    m = len(vals) // nb
    return vals[: m * nb].reshape(nb, m).mean(1)


def batch_err(vals, nb=20):
    """mean and standard error from nb batch means of a stationary series."""
    if len(vals) < nb:
        return float(np.mean(vals)), float("nan")
    b = batches(vals, nb)
    return float(b.mean()), float(b.std(ddof=1) / np.sqrt(nb))


def pooled(bs, ws):
    """weighted mean of batch means from several chains (weight = batch length) and its standard error."""
    b = np.concatenate(bs); w = np.concatenate(ws); w = w / w.sum()
    m = float((w * b).sum())
    return m, float(np.sqrt((w ** 2 * (b - m) ** 2).sum() * len(b) / (len(b) - 1)))


def read_parts(base):
    """Concatenate ts_<base>.part*.txt in order; returns (rows, total seconds)."""
    parts = sorted(glob.glob(base + ".part*.txt"), key=lambda f: int(re.search(r"part(\d+)", f).group(1)))
    ds = [read(f) for f in parts]
    ds = [d for d in ds if len(d) > 1]
    if not ds:
        return np.zeros((0, len(COLS))), 0.0
    return ds, float(sum(d[-1, 1] for d in ds))


def chain_summary(base, burn=0.2, burn_n=0.0):
    ds, secs = read_parts(base)
    n = int(re.search(r"_n(\d+)$", base).group(1))
    total_moves = float(sum(d[-1, 0] for d in ds))
    every = ds[0][1, 0] - ds[0][0, 0]
    # Burn-in: the first 20% of the run, and never less than BURN_N * n moves.  A stage starts from
    # a 2x-inflated sample, whose local structure (perimeter, cut cells) relaxes as exp(-t / ~43n).
    d = np.concatenate(ds)
    nburn = int(max(len(d) * burn, burn_n * n / every))
    d = d[nburn:]
    if len(d) < 30:
        raise ValueError("too short after burn-in")
    out = {"n": n, "file": os.path.basename(base), "samples": len(d), "every": every,
           "moves": total_moves, "sec_per_move": secs / total_moves, "burn_moves": nburn * every}
    c = {k: i for i, k in enumerate(COLS)}
    for k in ("rg2", "asph", "perim", "noncut", "holes"):
        t, w, e = tau_int(d[:, c[k]])
        out["tau_" + k] = float(t * every)
        out["tau_" + k + "_err"] = float(e * every)
    out["neff"] = float(len(d) * every / (2 * out["tau_rg2"]))
    # Error bars from batch means.  Batches must be long compared with the autocorrelation time
    # (about 15n moves for Rg^2), so at most 20 of them and none shorter than BATCH_N * n moves.
    nb = int(min(20, max(4, len(d) * every / (BATCH_N * n))))
    out["nb"] = nb
    out["_b"], out["_w"] = {}, np.full(nb, float(len(d) // nb))
    for k in ("rg2", "asph", "perim", "noncut", "holes", "cycles", "leaves", "deg1", "deg2", "deg3", "deg4"):
        m, e = batch_err(d[:, c[k]], nb)
        out[k], out[k + "_err"] = m, e
        out["_b"][k] = batches(d[:, c[k]], nb)
    m, e = batch_err(np.sqrt(d[:, c["rg2"]]), nb)
    out["rg"], out["rg_err"] = m, e
    out["_b"]["rg"] = batches(np.sqrt(d[:, c["rg2"]]), nb)
    out["rg_sd"] = float(np.sqrt(d[:, c["rg2"]]).std())
    out["rg2_sd"] = float(d[:, c["rg2"]].std())
    # atmospheric estimate of A(n+1)/A(n) = <t>_n / <k>_{n+1}, with <k>_{n+1} = <k>_n + kappa
    out["_tb"], out["_kb"] = out["_b"]["perim"].tolist(), out["_b"]["noncut"].tolist()
    return out


def combine(chains, key):
    """Pool the batches of all chains at one size.  (Inverse-variance weighting of per-chain means
    would favour whichever chain happened to get a small error estimate from its few batches.)"""
    return pooled([c["_b"][key] for c in chains], [c["_w"] for c in chains])


def main():
    prod = os.path.join(ROOT, "prod")
    hyb = {}
    bases = sorted({re.sub(r"\.part\d+\.txt$", "", f) for f in glob.glob(prod + "/ts_[ABC]_n*.part*.txt")})
    for fn in bases:
        try:
            s = chain_summary(fn, burn_n=BURN_N)
        except Exception as ex:
            print("skip", fn, ex, file=sys.stderr); continue
        if s["samples"] < 60 or s["neff"] < 8:
            continue
        hyb.setdefault(s["n"], []).append(s)
    sizes = []
    for n in sorted(hyb):
        ch = hyb[n]
        row = {"n": n, "chains": len(ch), "neff": sum(c["neff"] for c in ch)}
        for key in ("rg2", "rg", "asph", "perim", "noncut", "holes", "cycles", "leaves", "deg1", "deg2", "deg3", "deg4"):
            row[key], row[key + "_err"] = combine(ch, key)
        # autocorrelation times: plain mean over chains (their error bars scale with the estimate
        # itself, so inverse-variance weights would favour whichever chain came out low)
        for key in ("tau_rg2", "tau_perim", "tau_noncut", "tau_asph", "tau_holes"):
            row[key] = float(np.mean([c[key] for c in ch]))
            row[key + "_err"] = float(np.sqrt(sum(c[key + "_err"] ** 2 for c in ch)) / len(ch))
        row["burn_n"] = float(np.mean([c["burn_moves"] for c in ch]) / n)
        row["sec_per_move"] = float(np.mean([c["sec_per_move"] for c in ch]))
        row["rg2_sd"] = float(np.mean([c["rg2_sd"] for c in ch]))
        row["rg_sd"] = float(np.mean([c["rg_sd"] for c in ch]))
        row["moves"] = float(sum(c["moves"] for c in ch))
        row["chain_rg2"] = [[c["rg2"], c["rg2_err"]] for c in ch]
        row["_tb"] = np.concatenate([c["_tb"] for c in ch]).tolist()
        row["_kb"] = np.concatenate([c["_kb"] for c in ch]).tolist()
        row["_bw"] = np.concatenate([c["_w"] for c in ch])
        row["batches"] = int(sum(c["nb"] for c in ch))
        sizes.append(row)
    # kappa = limiting density of non-cut cells: slope of <k> vs n from the two largest sizes
    big = [r for r in sizes if r["n"] >= 4000]
    if len(big) >= 2:
        ns = np.array([r["n"] for r in big], float); ks = np.array([r["noncut"] for r in big])
        kappa = float(np.polyfit(ns, ks, 1)[0])
    else:
        kappa = 0.2942
    for r in sizes:
        t = np.array(r.pop("_tb")); k = np.array(r.pop("_kb"))
        ratio = t / (k + kappa)                       # batchwise A(n+1)/A(n)
        lam = ratio * (1 + 1.0 / r["n"])              # remove the n^-theta factor, theta = 1
        w = r.pop("_bw")
        r["ratio"], r["ratio_err"] = pooled([ratio], [w])
        r["lam"], r["lam_err"] = pooled([lam], [w])
    # a dedicated long run of local observables only (-mode means), folded into its size's estimate
    ded = None
    mf = os.path.join(ROOT, "runs", "means2000_C.txt")
    if os.path.exists(mf):
        lines = [l.split() for l in open(mf) if l.strip()]
        nd, nsamp = int(lines[0][1]), float(lines[0][3])          # "n N samples S"
        m = {l[0]: [float(v) for v in l[1:]] for l in lines[1:]}
        # perimeter and removable-cell errors are anticorrelated at worst: add them linearly
        lam_d = m["perim"][0] / (m["noncut"][0] + kappa) * (1 + 1.0 / nd)
        err_d = lam_d * (m["perim"][1] / m["perim"][0] + m["noncut"][1] / m["noncut"][0])
        ded = {"n": nd, "samples": nsamp, "lam": lam_d, "lam_err": err_d, "rg2": m["rg2"][0], "rg2_err": m["rg2"][1],
               "perim": m["perim"][0], "noncut": m["noncut"][0]}
        for r in sizes:
            if r["n"] == nd:
                r["lam_ladder"], r["lam_ladder_err"] = r["lam"], r["lam_err"]
                w1, w2 = 1 / r["lam_err"] ** 2, 1 / err_d ** 2
                r["lam"], r["lam_err"] = (w1 * r["lam"] + w2 * lam_d) / (w1 + w2), 1 / np.sqrt(w1 + w2)
    reloc = []
    rb = sorted({re.sub(r"\.part\d+\.txt$", "", f) for f in glob.glob(prod + "/ts_reloc_n*.part*.txt")}, key=lambda f: int(re.search(r"_n(\d+)", f).group(1)))
    for fn in rb:
        try:
            s = chain_summary(fn)
        except Exception as ex:
            print("skip", fn, ex, file=sys.stderr); continue
        if s["neff"] < 15:
            continue
        reloc.append({k: s[k] for k in ("n", "tau_rg2", "tau_rg2_err", "tau_perim", "tau_perim_err", "rg2", "rg2_err", "sec_per_move", "neff")})
    out = {"lambda_ref": LAMBDA, "kappa": kappa, "burn_n": BURN_N, "dedicated": ded, "sizes": sizes, "reloc": reloc}
    # the head-to-head from a square at n = 30,000
    p2 = os.path.join(ROOT, "prod2")
    for tag in ("old", "new"):
        ds, _ = read_parts(f"{p2}/ts_{tag}sq_n30000")
        if len(ds):
            off = 0.0; rows = []
            for d in ds:
                rows.append(np.column_stack([d[:, 0] + off, d[:, 2]])); off += d[-1, 0]
            d = np.concatenate(rows)
            idx = np.unique(np.round(np.geomspace(1, len(d), 500)).astype(int) - 1)
            out[f"square_{tag}"] = [[float(d[i, 0]), float(d[i, 1])] for i in idx]
            if tag == "new":
                late = np.sqrt(d[d[:, 0] > 2e7, 1])
                out["square_eq"] = {"rg": float(late.mean()), "rg_sd": float(late.std())}
    json.dump(out, open(os.path.join(ROOT, "analysis", "data.json"), "w"), indent=1)
    for r in sizes:
        print(f"n={r['n']:>8} chains={r['chains']} neff={r['neff']:.0f} rg2={r['rg2']:.5g}±{r['rg2_err']:.2g} rg={r['rg']:.5g} tau={r['tau_rg2']:.3g} ({r['tau_rg2']*r['sec_per_move']:.3g}s) perim/n={r['perim']/r['n']:.5f} noncut/n={r['noncut']/r['n']:.5f} holes/n={r['holes']/r['n']:.5f} lam={r['lam']:.5f}±{r['lam_err']:.5f}  chains: " + " ".join(f"{a:.5g}±{b:.2g}" for a, b in r["chain_rg2"]))
    for r in reloc:
        print(f"reloc n={r['n']} tau={r['tau_rg2']:.3g}±{r['tau_rg2_err']:.2g} neff={r['neff']:.0f}")


    return out


if __name__ == "__main__":
    main()
