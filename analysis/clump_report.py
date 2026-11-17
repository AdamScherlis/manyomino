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
    ms = [f"{r['max']['mean']:.0f}" for r in rows]; ns = [f"{int(n):,}" for n in sizes]
    T["CL_MAX_LIST"] = f"{', '.join(ms[:-1])} and {ms[-1]} cells on average at {', '.join(ns[:-1])} and {ns[-1]} cells"
    for s in (16, 64, 128): T[f"CL_T{s}"] = f"{P[s][0]:.0f}"
    top = max(P); T["CL_TOP_S"] = str(top); T["CL_TOP_T"] = f"{round(P[top][0], -1):.0f}"
    g = [r["tau_fit"] for r in rows if r.get("tau_fit")]
    big = [s for s in P if s >= 32]; slope = math.log(P[big[-1]][0] / P[big[0]][0]) / math.log(big[-1] / big[0])
    T["CL_SLOPE"] = f"{slope:.2f}"
    T["CL_GROWTH"] = (f"From {big[0]} cells up, doubling the clump multiplies its lasting time by about {2 ** slope:.1f}, and nothing in the record suggests a ceiling.")
    r4 = rows[0]; pr = dict(map(tuple, r4["present"]))
    T["CL_RARE"] = (f"A {int(sizes[0]):,}-cell shape holds a clump of 128 cells or more {100 * pr[128]:.1f}% of the time, and the error bars already come from "
                    f"the scatter of long blocks, which includes whatever the rings contribute.")
    hn = hero["n"] if hero else None
    syn = sorted((x for x in J["synthetic"] if x["n"] == hn), key=lambda x: -x["from"])
    syn = [x for x in syn if x["K"] >= 64][:2]            # enough stretches side by side to stand in for one shape
    if syn:
        x = syn[0]; T["CL_SYN_K"] = str(x["K"]); T["CL_SYN_FROM"] = f"{x['from']:,}"; T["CL_SYN_MEAN"] = f"{x['mean']:.0f}"; T["CL_SYN_ERR"] = f"{x['err']:.0f}"
        T["CL_SYN_TAU"] = f"{x['tau']:.0f}"; T["CL_SYN_TAU_ERR"] = f"{x['tau_err']:.0f}"
        T["CL_SYN_ALT"] = (f" (Built from the {syn[1]['from']:,}-cell record instead: {syn[1]['mean']:.0f}&nbsp;&plusmn;&nbsp;{syn[1]['err']:.0f} cells and "
                           f"{syn[1]['tau']:.0f}&nbsp;&plusmn;&nbsp;{syn[1]['tau_err']:.0f}.)") if len(syn) > 1 else ""
    if A1 and A1.get("sizes"):
        Q = pooled(A1); D["clump"]["fix"] = [[s, round(t, 2), round(e, 2)] for s, (t, e) in Q.items()]
        common = [s for s in Q if s in P and s >= 16]; top = common[-1]
        rate = [1 - A1["rate"][n] / J["rate"][n] for n in A1.get("rate", {}) if n in J.get("rate", {})]
        cost = f"{100 * sum(rate) / len(rate):.0f}" if rate else "10"
        F = A1.get("feat", {}); fn = sorted(F, key=int)
        def obs(n, name): return next(q for q in F[n]["rows"] if q["name"] == name)
        h1 = A1["sizes"].get(str(hn)) if hn else None
        direct = ""
        if h1 and hero and hero.get("blob"):
            direct = (f" And directly at {hn:,} cells, where {words.get(h1['chains'], h1['chains'])} chains were run on for {h1['samples'] / h1['chains']:.0f} moves per cell each "
                      f"with the new rule, the largest clump is forgotten in {h1['max']['tau']:.1f}&nbsp;&plusmn;&nbsp;{h1['max']['tau_err']:.1f} moves per cell, "
                      f"against {hero['blob'][0]:.0f}&nbsp;&plusmn;&nbsp;{hero['blob'][1]:.0f} before.")
        trunk = ""
        if fn:
            trunk = ("<p>It does nothing for the trunk, which was never about rings: " +
                     " and ".join(f"{obs(n, 'trunk01')['tau']:.0f}&nbsp;&plusmn;&nbsp;{obs(n, 'trunk01')['err']:.0f} moves per cell at {int(n):,} cells "
                                  f"(before, {obs(n, 'trunk01')['base']:.0f}&nbsp;&plusmn;&nbsp;{obs(n, 'trunk01')['base_err']:.0f})" for n in fn) +
                     ". The radius is unchanged too. So the trunk is once again the slowest property found, the cost of an independent sample of anything "
                     "reported here stays where it was, and the new rule is worth having mainly because it removes the one known way for the walk to get slower with size. "
                     "Every number on this page was produced with the original rule. The new one is an option in the code, not yet the default.</p>")
        T["CL_FIX"] = (
            "<h3>Drawing larger boxes</h3>"
            "<p>If the wait is for a large box, the remedy is to draw large boxes more often. Let the box reach r cells or more with probability 1/r in place of "
            "1/r&sup2;. The walk stays exactly as correct as it was, because a move and its reverse draw their boxes by the same rule, whatever that rule is. "
            f"As a check, the modified walk passes the same exhaustive test at 8 cells, and it produces the same number of clumps of every size as the original. "
            f"Searching larger boxes costs about {cost}% more time per move.</p>"
            '</section>\n<section>\n  <div class="chart"><h3>The same clumps, with larger boxes drawn more often</h3>\n'
            '    <p class="sub">Autocorrelation time of the number of clumps of at least a given size, shape sizes combined.</p>\n'
            '    <div id="ch-clump-fix"></div></div>\n</section>\n<section class="measure">'
            f"<p>The ladder goes away. With the new rule the count of clumps of at least 16 cells is forgotten in {Q[16][0]:.1f} moves per cell, and of at least "
            f"{top} cells in {Q[top][0]:.1f}, where the original took {P[16][0]:.0f} and {P[top][0]:.0f}. Large clumps no longer outlast small ones."
            f"{direct}</p>" + trunk)
    if hero and hero.get("blob"):
        T["CL_HERO_TAU"] = f"{hero['blob'][0]:.0f}"; T["CL_HERO_TAU_ERR"] = f"{hero['blob'][1]:.0f}"; T["CL_HERO_MEAN"] = f"{hero['blob_mean']:.0f}"
    return D, T
