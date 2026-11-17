"""Page numbers and chart data for the clump-lifetime measurements (analysis/clumps.py)."""
import json, math, os

ROOT = os.path.join(os.path.dirname(__file__), "..")


def load(name):
    fn = os.path.join(ROOT, "analysis", name)
    return json.load(open(fn)) if os.path.exists(fn) else None


def pooled(J, lo=8):
    """lasting time at each clump size, combined over the shape sizes run (weights 1/err^2);
    only sizes timed to better than 40%"""
    acc = {}
    for r in J["sizes"].values():
        for s, t, te, _ in r["tau"]:
            if s < lo or not te == te or te <= 0 or te > 0.4 * t: continue
            w = 1 / te ** 2; a = acc.setdefault(s, [0.0, 0.0]); a[0] += w * t; a[1] += w
    return {s: (a[0] / a[1], 1 / math.sqrt(a[1])) for s, a in sorted(acc.items())}


def series(J, lo=8, rel=0.4):
    out = []
    for n, r in sorted(J["sizes"].items(), key=lambda kv: int(kv[0])):
        out.append({"n": int(n),
                    "tau": [[s, round(t, 2), round(te, 2)] for s, t, te, _ in r["tau"] if s >= lo and te == te and 0 < te <= rel * t],
                    "dens": [[s, float(f"{d:.4g}"), float(f"{de:.2g}")] for s, d, de in r["density"] if s >= 4 and 0 < de <= rel * d]})
    return out


def clump_report(hero=None):
    """(chart data, text tokens) for the page; hero is D["hero_slow"] when the largest runs were timed"""
    J = load("clumps.json"); ring = load("clump_ring.json"); A1 = load("clumps_a1.json")
    T = {k: "" for k in ("CL_FIX", "CL_RARE", "CL_GROWTH")}; D = {}
    if not J or not ring: return D, T
    P = pooled(J); sizes = sorted(J["sizes"], key=int); rows = [J["sizes"][n] for n in sizes]
    D["clump"] = {"ring": ring, "base": series(J), "pooled": [[s, round(t, 2), round(e, 2)] for s, (t, e) in P.items()]}
    words = {2: "two", 3: "three", 4: "four", 5: "five"}
    nch = sum(r["chains"] for r in rows); tot = sum(r["samples"] for r in rows)
    T["CL_RING_CELLS"] = str(ring["cells"]); T["CL_NCHAINS"] = words.get(nch, str(nch)); T["CL_READINGS"] = f"{tot / 1e6:.1f} million"
    f = [r["density_fit"]["exponent"] for r in rows if r.get("density_fit")]
    T["CL_DENS_EXP"] = f"{-sum(f) / len(f):.2f}"; T["CL_DENS_FACTOR"] = f"{2 ** (-sum(f) / len(f)):.0f}"
    parts = [f"{r['max']['mean']:.0f} cells on average at {int(n):,}" for n, r in zip(sizes, rows)]
    T["CL_MAX_LIST"] = ", ".join(parts[:-1]) + " and " + parts[-1]
    for s in (16, 64, 128): T[f"CL_T{s}"] = f"{P[s][0]:.0f}"
    top = max(P); T["CL_TOP_S"] = str(top); T["CL_TOP_T"] = f"{round(P[top][0], -1):.0f}"
    g = [r["tau_fit"] for r in rows if r.get("tau_fit")]
    big = [s for s in P if s >= 32]; slope = math.log(P[big[-1]][0] / P[big[0]][0]) / math.log(big[-1] / big[0])
    T["CL_SLOPE"] = f"{slope:.2f}"
    T["CL_GROWTH"] = (f"From {big[0]} cells up, doubling the clump multiplies its lasting time by about {2 ** slope:.1f}, and nothing in the record suggests a ceiling.")
    r4 = rows[0]; pr = dict(map(tuple, r4["present"]))
    T["CL_RARE"] = (f"A {int(sizes[0]):,}-cell shape holds a clump of 128 cells or more {100 * pr[128]:.1f}% of the time, and the error bars already come from "
                    f"the scatter of long blocks, which includes whatever the rings contribute.")
    syn = [x for x in J["synthetic"] if x["from"] == int(sizes[0]) and (hero is None or x["n"] == hero["n"])]
    if syn:
        x = syn[-1]; T["CL_SYN_K"] = str(x["K"]); T["CL_SYN_MEAN"] = f"{x['mean']:.0f}"; T["CL_SYN_ERR"] = f"{x['err']:.0f}"
        T["CL_SYN_TAU"] = f"{x['tau']:.0f}"; T["CL_SYN_TAU_ERR"] = f"{x['tau_err']:.0f}"
    if hero and hero.get("blob"):
        T["CL_HERO_TAU"] = f"{hero['blob'][0]:.0f}"; T["CL_HERO_TAU_ERR"] = f"{hero['blob'][1]:.0f}"; T["CL_HERO_MEAN"] = f"{hero['blob_mean']:.0f}"
    return D, T
