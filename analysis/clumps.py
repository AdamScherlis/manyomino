"""How long do clumps last, and how does that depend on their size?

A clump is a maximal set of cells that stays connected when any one adjacency between cells is
cut (the cells "held together by loops").  The largest clump was the one measurement that read
slower at 1,024,000 cells than at smaller sizes (analysis/hero_slow.py).  This script tests the
explanation: clumps are local objects, bigger ones last longer, and a bigger shape holds bigger
ones.

scripts/clumps.sh records, once per move per cell, how many clumps of at least each size a shape
holds, on equilibrium chains of 4,000, 16,000 and 64,000 cells.  From those:

  density      clumps of at least s cells, per cell of the shape (the same at every shape size
               if clumps are local)
  persistence  the autocorrelation time of that count, in moves per cell: how long a clump of
               at least s cells, once seen, goes on being one
  largest      the autocorrelation time of the size of the largest clump, at the sizes run and
               for larger shapes put together from independent stretches of the smaller ones
               (the largest clump of K independent shapes is the largest of their largest)

The relocate move may lift a cell off a loop only if its random connectivity radius K reaches
round the loop, and P(K >= r) = r^-A with A = 2 by default.  scripts/clumps.sh T 1 repeats the
runs with A = 1 (-relA 1), which draws large radii far more often.

  python3 analysis/clumps.py          # reads clump/ts_*_n*.part*.blobs, writes analysis/clumps.json
  python3 analysis/clumps.py clump_a1 # the A = 1 runs, writes analysis/clumps_a1.json
"""
import glob, json, math, os, re, sys
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
NBLOCK = 8


def load(sub="clump"):
    """{n: [one record array per chain]}; a record is (largest clump, counts at each threshold)"""
    out = {}; th = None
    files = sorted(glob.glob(os.path.join(ROOT, sub, "ts_*_n*.part*.blobs")), key=lambda f: (int(re.search(r"_n(\d+)", f).group(1)), os.path.basename(f)[3], int(re.search(r"part(\d+)", f).group(1))))
    for fn in files:
        th = np.array(open(fn + ".sizes").read().split(), int)
        rec = np.dtype([("max", "<u4"), ("c", "<u2", (len(th),))])
        a = np.fromfile(fn, dtype=rec, count=os.path.getsize(fn) // rec.itemsize)
        n = int(re.search(r"_n(\d+)", fn).group(1)); tag = os.path.basename(fn)[3]
        if len(a): out.setdefault(n, {}).setdefault(tag, []).append(a)
    out = {n: [np.concatenate(v) for v in d.values()] for n, d in out.items()}
    return th, {n: [c for c in ch if len(c) >= 64] for n, ch in out.items() if any(len(c) >= 64 for c in ch)}       # too short to time: left out


def acov_sums(x, maxlag):
    """sums of x[i] x[i+t] for t = 0..maxlag, by FFT"""
    n = len(x); m = 1 << int(math.ceil(math.log2(n + maxlag + 1)))
    f = np.fft.rfft(x, m); r = np.fft.irfft(f * np.conj(f), m)[: maxlag + 1]
    return r


def tau_int(series, c=6.0):
    """integrated autocorrelation time, in samples, of several stretches of one stationary
    series about their common mean (Sokal's window)"""
    mean = sum(s.sum() for s in series) / sum(len(s) for s in series)
    maxlag = min(len(s) for s in series) // 2
    S = np.zeros(maxlag + 1); C = np.zeros(maxlag + 1)
    for s in series:
        S += acov_sums(s - mean, maxlag); C += len(s) - np.arange(maxlag + 1)
    if S[0] <= 0: return float("nan"), 0
    rho = (S / C) / (S[0] / C[0])
    cs = 0.5 + np.cumsum(rho[1:]); w = np.arange(1, len(cs) + 1)
    ok = np.nonzero(w >= c * cs)[0]; i = ok[0] if len(ok) else len(cs) - 1
    return float(cs[i]), int(w[i])


def tau_err(chains):
    """time and its error: the chains are cut into blocks, and blocks are left out in turn"""
    blocks = [b for x in chains for b in np.array_split(np.asarray(x, float), NBLOCK)]
    t = tau_int(blocks)[0]
    jk = [tau_int(blocks[:k] + blocks[k + 1:])[0] for k in range(len(blocks))]
    jk = np.array([v for v in jk if v == v]); B = len(jk)
    return t, (math.sqrt((B - 1) / B * ((jk - jk.mean()) ** 2).sum()) if B > 2 else float("nan"))


def mean_err(chains):
    blocks = np.array([b.mean() for x in chains for b in np.array_split(np.asarray(x, float), NBLOCK * 4)])
    return float(np.concatenate([np.asarray(x, float) for x in chains]).mean()), float(blocks.std(ddof=1) / math.sqrt(len(blocks)))


def episodes(chains, k):
    """lengths of the unbroken stretches during which a clump of at least threshold k exists"""
    out = []
    for x in chains:
        on = (x["c"][:, k] > 0).astype(np.int8); d = np.diff(np.concatenate([[0], on, [0]]))
        out += list(np.nonzero(d == -1)[0] - np.nonzero(d == 1)[0])
    return np.array(out)


def synthetic(chains, K, shifts=None):
    """the largest clump of K independent shapes: K stretches of the record, far apart in time,
    laid side by side.  shifts slides each stretch round on itself by a different amount, which
    changes which moments of the K shapes coincide without changing any one of them."""
    x = np.concatenate([c["max"] for c in chains]); P = len(x) // K
    rows = x[: P * K].reshape(K, P)
    if shifts is not None: rows = np.stack([np.roll(r, int(d)) for r, d in zip(rows, shifts)])
    return rows.max(axis=0).astype(float)


NARR = 24          # alignments tried for each put-together shape


def rates(sub):
    """million moves per second of each chain that recorded only clump sizes"""
    out = {}
    for fn in sorted(glob.glob(os.path.join(ROOT, sub, "ts_*_n*.part*.txt"))):
        if os.path.exists(fn.replace(".txt", ".feat")): continue
        last = open(fn).read().strip().split("\n")[-1].split(); n = re.search(r"_n(\d+)", fn).group(1)
        out.setdefault(n, []).append(float(last[0]) / float(last[1]) / 1e6)
    return {n: float(np.mean(v)) for n, v in out.items()}


def feat_compare(sub):
    """autocorrelation time of every library measurement in the chains of this directory that
    recorded them (-feat), next to the same estimate from the default-kernel chains in conv/"""
    def read(fn, burn=0.0):
        names = open(fn + ".names").read().split(); a = np.fromfile(fn, dtype=np.float64)
        a = a[: len(a) // len(names) * len(names)].reshape(-1, len(names)); return names[1:], a[int(len(a) * burn):, 1:]
    out = {}
    for fn in sorted(glob.glob(os.path.join(ROOT, sub, "ts_*_n*.part0.feat"))):
        n = int(re.search(r"_n(\d+)", fn).group(1)); tag = os.path.basename(fn)[3]
        mine = [read(f)[1] for f in sorted(glob.glob(fn.replace("part0", "part*"))) if not f.endswith(".names")]
        names = read(fn)[0]; mine = [np.vstack(mine)]
        base = [read(f, 0.1)[1] for f in sorted(glob.glob(os.path.join(ROOT, "conv", f"u{n}_?.feat")))]
        rows = []
        for j, nm in enumerate(names):
            if any(x[:, j].std() == 0 for x in mine + base): continue
            t, te = tau_err([x[:, j] for x in mine]); b, be = tau_err([x[:, j] for x in base])
            rows.append({"name": nm, "tau": t, "err": te, "base": b, "base_err": be})
        rows.sort(key=lambda q: -q["tau"])
        out[str(n)] = {"chain": tag, "len": int(sum(len(x) for x in mine)), "base_len": int(sum(len(x) for x in base)), "base_chains": len(base), "rows": rows}
    return out


def powerfit(rows, lo, hi):
    """weighted fit of log y = a + b log s over lo <= s <= hi; rows are [s, y, err, ...].
    Returns (b, error of b, y at s = 64 from the fit)"""
    r = np.array([q[:3] for q in rows if lo <= q[0] <= hi and q[1] > 0 and q[2] > 0 and q[2] == q[2]], float)
    if len(r) < 3: return None
    x = np.log(r[:, 0] / 64.0); y = np.log(r[:, 1]); w = (r[:, 1] / r[:, 2]) ** 2
    A = np.array([[w.sum(), (w * x).sum()], [(w * x).sum(), (w * x * x).sum()]]); cov = np.linalg.inv(A)
    a, b = cov @ np.array([(w * y).sum(), (w * x * y).sum()])
    chi2 = float((w * (y - a - b * x) ** 2).sum()); dof = len(r) - 2
    return {"exponent": float(b), "err": float(math.sqrt(cov[1, 1] * max(1.0, chi2 / dof))), "at64": float(math.exp(a)), "chi2": chi2, "dof": dof, "lo": lo, "hi": hi}


def main(sub="clump"):
    th, D = load(sub)
    out = {"dir": sub, "thresholds": th.tolist(), "sizes": {}, "synthetic": []}
    for n in sorted(D):
        ch = D[n]; total = sum(len(c) for c in ch)
        row = {"n": n, "chains": len(ch), "samples": total, "density": [], "tau": [], "episodes": [], "present": []}
        for k, s in enumerate(th):
            cnt = [c["c"][:, k] for c in ch]; m, e = mean_err(cnt)
            ep = episodes(ch, k)
            if m * total < 30 or (m < 0.5 and len(ep) < 20): break         # too few sightings to say anything
            t, te = tau_err(cnt)
            row["present"].append([int(s), float(sum((c > 0).sum() for c in cnt) / total)])
            row["density"].append([int(s), m / n, e / n]); row["tau"].append([int(s), t, te, m * total])
            row["episodes"].append([int(s), len(ep), float(ep.mean()) if len(ep) else 0.0, float((ep.astype(float) ** 2).mean() / (2 * ep.mean())) if len(ep) else 0.0])
        mx = [c["max"].astype(float) for c in ch]; m, e = mean_err(mx); t, te = tau_err(mx)
        row["max"] = {"mean": m, "err": e, "tau": t, "tau_err": te, "median": float(np.median(np.concatenate(mx)))}
        row["tau_fit"] = powerfit(row["tau"], 16, 1e9); row["density_fit"] = powerfit(row["density"], 16, 1e9)
        out["sizes"][str(n)] = row
        print(f"n={n:,}: {len(ch)} chain(s), {total:,} samples; largest clump {m:.1f} ± {e:.1f} cells, time {t:.1f} ± {te:.1f} moves per cell")
        for (s, d, de), (_, t, te, sight), (_, ne, el, er) in zip(row["density"], row["tau"], row["episodes"]):
            print(f"   at least {s:>4} cells: {d:.3e} ± {de:.1e} per cell   lasts {t:6.1f} ± {te:5.1f}   ({sight:,.0f} sightings, {ne:,} episodes of mean length {el:.1f})")
        for nm, key in (("lasting time", "tau_fit"), ("density", "density_fit")):
            f = row[key]
            if f: print(f"   {nm} goes as size^{f['exponent']:.2f} ± {f['err']:.2f} from 16 cells up ({f['at64']:.3g} at 64 cells; chi2 {f['chi2']:.0f} for {f['dof']})")
    # larger shapes put together from independent stretches of smaller ones
    for n in sorted(D):
        K = 4
        rng = np.random.default_rng(1)
        while True:
            P = sum(len(c) for c in D[n]) // K
            if P < 2000: break
            res = []
            for a in range(NARR):
                s = synthetic(D[n], K, None if a == 0 else rng.integers(0, P, K))
                t, te = tau_err([s]); res.append((t, te, s.mean()))
            res = np.array(res); t = float(np.median(res[:, 0])); te = float(max(res[:, 0].std(ddof=1), np.median(res[:, 1])))
            m = float(res[:, 2].mean()); e = mean_err([synthetic(D[n], K)])[1]
            out["synthetic"].append({"from": n, "K": K, "n": n * K, "len": P, "mean": m, "err": e, "tau": t, "tau_err": te,
                                     "tau_lo": float(res[:, 0].min()), "tau_hi": float(res[:, 0].max())})
            print(f"synthetic {n * K:>10,} cells from {K:>5} x {n:,}: largest clump {m:.1f} ± {e:.1f}, time {t:.1f} ± {te:.1f} "
                  f"({res[:, 0].min():.0f} to {res[:, 0].max():.0f} over {NARR} alignments; {P:,} moves per cell)")
            K *= 4
    out["rate"] = rates(sub)
    if sub != "clump":
        out["feat"] = feat_compare(sub)
        for n, r in out["feat"].items():
            print(f"n={int(n):,}: every measurement timed with the default radius law (conv/, {r['base_len']:,} samples) and with this one ({r['len']:,} samples), moves per cell")
            for q in r["rows"][:8] + [q for q in r["rows"][8:] if q["name"] in ("rg2", "trunk01", "blob_max")]:
                print(f"   {q['name']:>10}: {q['base']:6.1f} ± {q['base_err']:4.1f}  ->  {q['tau']:6.1f} ± {q['err']:4.1f}")
    name = "clumps.json" if sub == "clump" else sub.replace("clump", "clumps") + ".json"
    json.dump(out, open(os.path.join(ROOT, "analysis", name), "w"), indent=1)
    return out


if __name__ == "__main__":
    main(*sys.argv[1:2])
