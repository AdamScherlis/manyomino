"""How long do the largest runs remember each measurement?

The slow-mode search in analysis/convergence.py stops at 64,000 cells.  The extended runs at
1,024,000 cells (scripts/extend_hero.sh) record the same measurements every three moves per
cell, on four chains.  That is too few samples to search over combinations, but enough to time
each single measurement: the autocorrelation function is averaged over the chains (about their
common mean), summed with Sokal's window, and the error taken from leaving one chain out.

Runs this short (each chain is some twenty slow relaxation times long) make that estimate read
low.  So the same estimator is also applied to stretches of the 16,000- and 64,000-cell chains
cut to the same length and sampling interval, in moves per cell, where the long-run answer is
known from analysis/convergence.py.  The comparison to make is like with like: the estimate at
1,024,000 cells against the estimates from equally short stretches at the smaller sizes.

  python3 analysis/hero_slow.py        # reads prod/ts_?_n1024000.part*.feat, writes analysis/hero_slow.json
"""
import glob, json, math, os, re
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
N = 1024000


def chains():
    """one array of measurements per chain, its parts joined (a part resumes from a checkpoint
    at most 90 seconds, about one sample, before the end of the previous one)"""
    out = {}
    for fn in sorted(glob.glob(os.path.join(ROOT, "prod", f"ts_?_n{N}.part*.feat")), key=lambda f: (f.split("ts_")[1][0], int(re.search(r"part(\d+)", f).group(1)))):
        names = open(fn + ".names").read().split()
        a = np.fromfile(fn); a = a[: len(a) // len(names) * len(names)].reshape(-1, len(names))
        if len(a) < 2: continue
        tag = os.path.basename(fn)[3]
        out.setdefault(tag, {"rows": [], "step": float(a[1, 0] - a[0, 0]) / N})["rows"].append(a[:, 1:])
    return names[1:], {t: {"X": np.vstack(v["rows"]), "step": v["step"]} for t, v in out.items()}


def acov(x, maxlag, mean):
    """autocovariance sums and pair counts of one series about a given mean"""
    x = x - mean; n = len(x); L = min(maxlag, n - 1)
    return np.array([float(x[: n - t] @ x[t:]) for t in range(L + 1)]), np.arange(n, n - L - 1, -1).astype(float)


def pooled_tau(series, c=6.0):
    """Sokal estimate from the chain-averaged autocorrelation function; in samples"""
    maxlag = min(len(s) for s in series) // 2
    S = np.zeros(maxlag + 1); C = np.zeros(maxlag + 1); mean = np.concatenate(series).mean()
    for s in series:
        a, k = acov(s, maxlag, mean); S[: len(a)] += a; C[: len(k)] += k
    rho = (S / C) / (S[0] / C[0])
    cs = 0.5 + np.cumsum(rho[1:]); w = np.arange(1, len(cs) + 1)
    ok = np.nonzero(w >= c * cs)[0]; i = ok[0] if len(ok) else len(cs) - 1
    return float(cs[i]), int(w[i])


def calibrate(n, seg, step, k, want):
    """the same estimator on stretches of the chains at n cells: sets of k stretches, each seg
    moves per cell long, sampled every step moves per cell.  Returns {name: (mean, sd, sets)}."""
    segs = []
    for fn in sorted(glob.glob(os.path.join(ROOT, "conv", f"u{n}_?.feat"))):
        names = open(fn + ".names").read().split()
        a = np.memmap(fn, dtype=np.float64, mode="r"); a = a[: len(a) // len(names) * len(names)].reshape(-1, len(names))
        L = int(seg)                                            # these chains are sampled every move per cell
        segs += [np.array(a[i:i + L:int(step), 1:]) for i in range(0, len(a) - L + 1, L)]
    names = names[1:]; rng = np.random.default_rng(1); rng.shuffle(segs)
    sets = [segs[i:i + k] for i in range(0, len(segs) - k + 1, k)]
    out = {}
    for nm in want:
        j = names.index(nm); v = [pooled_tau([s[:, j] for s in st])[0] * step for st in sets]
        out[nm] = (float(np.mean(v)), float(np.std(v, ddof=1)), len(v), float(np.mean([s[:, j].mean() for s in segs])))
    return out


def main():
    names, CH = chains()
    tags = sorted(CH); step = CH[tags[0]]["step"]
    total = sum(len(CH[t]["X"]) for t in tags)
    res = []
    for j, nm in enumerate(names):
        ser = [CH[t]["X"][:, j] for t in tags]
        if min(s.std() for s in ser) == 0: continue
        tau, w = pooled_tau(ser)
        # leave one chain out
        jk = [pooled_tau(ser[:k] + ser[k + 1:])[0] for k in range(len(ser))] if len(ser) > 2 else []
        err = math.sqrt((len(jk) - 1) / len(jk) * sum((v - np.mean(jk)) ** 2 for v in jk)) if jk else tau * math.sqrt(2 * (2 * w + 1) / total)
        res.append({"name": nm, "tau": tau * step, "err": err * step, "window": w * step,
                    "own": [pooled_tau([x])[0] * step for x in ser],          # each chain about its own mean (reads low)
                    "means": [float(x.mean()) for x in ser]})
    res.sort(key=lambda r: -r["tau"])
    rg = next(r for r in res if r["name"] == "rg2")
    out = {"n": N, "chains": tags, "samples": {t: len(CH[t]["X"]) for t in tags}, "step": step, "moves_per_cell": total * step,
           "rg2": rg, "slowest": res[:8], "all": res, "calibration": {}}
    # the same estimator on equally short stretches at smaller sizes, and the long-run values there
    want = ["rg2"] + [r["name"] for r in res[:8] if r["name"] != "rg2"]
    try:
        conv = {s["n"]: {x["name"]: x["tau"] for x in s["single"]} for s in json.load(open(os.path.join(ROOT, "analysis", "convergence.json")))["slow"]}
    except Exception:
        conv = {}
    seg = total * step / len(tags)
    for n in (16000, 64000):
        cal = calibrate(n, seg, step, len(tags), want)
        out["calibration"][str(n)] = {nm: {"short": cal[nm][0], "sd": cal[nm][1], "sets": cal[nm][2], "mean": cal[nm][3], "long": conv.get(n, {}).get(nm)} for nm in want}
    # how long the rare largest clumps last: the longest unbroken stretch any chain spends with
    # its largest clump above the level exceeded a tenth of the time
    j = names.index("blob_max"); allb = np.concatenate([CH[t]["X"][:, j] for t in tags]); thr = float(np.percentile(allb, 90)); best = 0
    for t in tags:
        c = 0
        for b in CH[t]["X"][:, j] > thr:
            c = c + 1 if b else 0; best = max(best, c)
    out["blob"] = {"mean": float(allb.mean()), "median": float(np.median(allb)), "p90": thr, "longest_above_p90": best * step}
    json.dump(out, open(os.path.join(ROOT, "analysis", "hero_slow.json"), "w"), indent=1)
    print(f"n={N}: chains {','.join(tags)}, {total} samples every {step:g} moves per cell = {total * step:.0f} moves per cell in all")
    print(f"  squared radius: {rg['tau']:.1f} ± {rg['err']:.1f} moves per cell")
    for r in res[:8]:
        print(f"  {r['name']:>10}: {r['tau']:.1f} ± {r['err']:.1f}   (window {r['window']:.0f})" + "".join(
            f"   at {int(n):,}: short {c[r['name']]['short']:.1f} ± {c[r['name']]['sd']:.1f} (spread over {c[r['name']]['sets']} sets), long {c[r['name']]['long'] or float('nan'):.1f}"
            for n, c in out["calibration"].items() if r["name"] in c))
    for n, c in out["calibration"].items():
        print(f"  squared radius at {int(n):,}: short {c['rg2']['short']:.1f} ± {c['rg2']['sd']:.1f}, long {c['rg2']['long'] or float('nan'):.1f}")
    return out


if __name__ == "__main__":
    main()
