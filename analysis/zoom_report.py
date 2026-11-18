"""Page data and text for the endless zoom (page/zoom.js, src/zoom.cpp), from
analysis/zoom_relax.json (how a doubled shape settles) and analysis/zoom_test.json (refined
shapes against sampled ones at a million cells).  Run directly, it writes analysis/zoom.md, the
tables for the README."""
import base64, json, math, os
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SWEEPS = 32                   # rounds per level in the browser
NU, NU_ERR = 0.6412, 0.0005


def load(name):
    fn = os.path.join(ROOT, "analysis", name)
    return json.load(open(fn)) if os.path.exists(fn) else None


def first_shape(fn=os.path.join(ROOT, "prod", "snap_A_n16000_p0_0001.bin")):
    a = np.fromfile(fn, dtype="<i4"); n = int(a[0]); xy = a[1:1 + 2 * n].reshape(n, 2)
    x = xy[:, 0] - xy[:, 0].min(); y = xy[:, 1] - xy[:, 1].min(); w = int(x.max()) + 1; h = int(y.max()) + 1
    bits = np.zeros(w * h, np.uint8); bits[y * w + x] = 1
    return {"n": n, "w": w, "h": h, "bits": base64.b64encode(np.packbits(bits, bitorder="little").tobytes()).decode()}


def pooled(pairs):
    w = [1 / e ** 2 for _, e in pairs]; return sum(m * k for (m, _), k in zip(pairs, w)) / sum(w), sum(w) ** -0.5


def growth_by_T(R, X):
    """growth factor at each number of rounds: every measurement of it pooled"""
    out = {}
    for i, T in enumerate(R["T"]): out.setdefault(T, []).extend(tuple(R["sets"][n]["growth"][i]) for n in R["sets"])
    for T, v in (X or {"T": {}})["T"].items(): out.setdefault(int(T), []).extend(tuple(g) for g in v["growth"].values())
    return {T: pooled(v) for T, v in sorted(out.items())}


def zoom_report():
    R, X = load("zoom_relax.json"), load("zoom_test.json")
    pred = 2 ** (1 / NU); pred_err = pred * math.log(2) / NU ** 2 * NU_ERR
    D = {"zoom": dict(first_shape(), sweeps=SWEEPS, growth=pred)}; Z = D["zoom"]
    T = {"ZM_N0": f"{Z['n']:,}", "ZM_T": str(SWEEPS), "ZM_RESULTS": ""}
    if not R: return D, T
    pm = lambda v, d=3: f"{v[0]:.{d}f}&nbsp;&plusmn;&nbsp;{v[1]:.{d}f}"
    G = growth_by_T(R, X); sets = R["sets"]; ti = {t: i for i, t in enumerate(R["T"])}
    Z["pred"] = [pred, pred_err]
    Z["mass"] = [{"name": f"{int(n):,} cells doubled", "pts": [[t, *sets[n]["growth"][i]] for i, t in enumerate(R["T"])]} for n in sorted(sets, key=int)]
    big = sets[max(sets, key=int)]
    loc = lambda nm, t: big["local"][nm][ti[t]]
    worst = lambda t: max(abs(loc(nm, t)[0] - 1) for nm in ("perim", "squares", "deg1", "deg3", "straight"))
    cl = big["clumps"]; ex = lambda t, s: cl[str(t)][str(s)][0] / cl[str(t)][str(s)][1]
    ch = lambda b, t: 100 * big["changed"][str(b)][ti[t]][0]
    gl = lambda t: f"{G[t][0]:.3f}"
    nu_of = lambda m: math.log(2) / math.log(m)
    out = [f"<p>The first is the number of squares. Doubling quadruples it. The walk then takes about a quarter of them away again, and the count comes to rest between 2.94 and 2.95 times "
           f"the original: {gl(32)} after 32 rounds, {gl(64)} after 64 and {gl(256)} after 256, each good to about &plusmn;{max(G[32][1], G[64][1], G[256][1]):.3f}. Scaling predicts a particular number here. "
           f"A shape twice as wide should have <span class=\"nw\">2<sup>1/&nu;</sup>&nbsp;=&nbsp;{pred:.3f}&nbsp;&plusmn;&nbsp;{math.ceil(pred_err * 1000) / 1000:.3f}</span> times the cells, the uncertainty being that of the published "
           f"exponent. Nothing in the procedure knows &nu;. The walk finds for itself how much shape to leave, to within two parts in a thousand, though the count is still creeping "
           f"downward at the longest runs and where it finally stops is not settled to better than that.</p>"]
    out.append('</section>\n<section>\n  <div class="two">\n    <div class="chart"><h3>Cells left after doubling</h3>\n'
               '      <p class="sub">Cells in the refined shape &divide; cells in the original, against rounds of the walk. The band is the prediction from the size exponent.</p>\n'
               '      <div id="ch-zoom-mass"></div></div>\n'
               '    <div class="chart"><h3>Refined four times, against sampled</h3>\n'
               '      <p class="sub">Each dot is one of 60 measurements: how far its average over the refined shapes lies from the sampled average, in standard errors. Dots on the top or bottom edge are further out than 6.</p>\n'
               '      <div id="ch-zoom-z"></div></div>\n  </div>\n</section>\n<section class="measure">')
    out.append(f"<p>The second is the fine texture. One round after doubling, the shape has {loc('squares', 1)[0]:.1f} times too many 2&times;2 blocks and about "
               f"{ex(1, 32):.0f} times too many loops of 32 cells or more. That excess is expected: every loop survives doubling at twice the width, and wide loops ought to be much rarer than "
               f"narrow ones. After 16 rounds the number of loops of 16 cells or more is {ex(16, 16):.2f} of what sampled shapes have, and after {SWEEPS} rounds "
               f"the perimeter, the 2&times;2 blocks, the tips and the branchings per cell are all within {100 * worst(SWEEPS):.1f}% of their sampled values.</p>")
    if X and X["T"]:
        names = X["names"]; L = set(X["local"]); isloc = [nm in L for nm in names]
        Z["test"] = [{"T": int(t), "local": [v["z"][j] for j in range(len(names)) if isloc[j]], "rest": [v["z"][j] for j in range(len(names)) if not isloc[j]]} for t, v in sorted(X["T"].items(), key=lambda kv: int(kv[0]))]
        def verdict(t):
            v = X["T"][str(t)]; z = np.array(v["z"]); lo = np.array(isloc)
            big2 = [(names[j], z[j], 100 * (v["mean"][j] / v["eq"][j] - 1)) for j in np.argsort(-np.abs(z)) if abs(z[j]) > 2.5]
            return v, z, lo, big2
        plain = {"leaves": "tips", "straight": "straight runs", "perim": "the perimeter", "deg2": "cells with two neighbours", "deg3": "branchings", "deg4": "cells with four neighbours",
                 "cycles": "loops", "blocks22": "2&times;2 blocks", "holes": "holes", "noncut": "removable cells", "blob_count": "clumps", "blob_cells": "cells in clumps",
                 "trunk01": "the trunk", "trunk03": "the trunk", "trunk10": "the trunk", "trunk30": "the trunk", "bw": "the width of the outline"}
        Ts = sorted(X["T"], key=int); best = Ts[-1]; v, z, lo, _ = verdict(best); nl = ~lo; nd = int((np.abs(z[nl]) > 2).sum())
        s = (f"<p>The third is everything else. {v['shapes']} samples of 13,650 cells were each refined four times over, to about a million cells, and their {len(names)} measurements "
             f"compared with the million-cell runs. {int(nl.sum())} of them describe the shape at large, from its radius and outline to the length of its trunk. Those agree at every "
             f"number of rounds tried: " + "; ".join(f"with {t} rounds, {int((np.abs(np.array(X['T'][t]['z'])[nl]) > 2).sum())}" for t in Ts)
             + f" lie beyond two standard errors, where chance alone would give about {0.0455 * nl.sum():.0f}. ")
        jb = int(np.argmax(np.where(nl, np.abs(z), 0))); st = X.get("start")
        if st and abs(st["z"][jb]) > 1.5 and st["z"][jb] * z[jb] > 0 and st["percent"][jb] is not None:
            s += (f"The largest, {plain.get(names[jb], names[jb])} at {abs(z[jb]):.1f} with {best} rounds, is not the refinement&rsquo;s doing. The {st['shapes']} starting samples happen to have "
                  f"{'longer trunks' if names[jb].startswith('trunk') else 'a larger value'} than average for their size, by {abs(st['percent'][jb]):.1f}%, and the refined shapes have kept that: "
                  f"{abs(100 * (v['mean'][jb] / v['eq'][jb] - 1)):.1f}%. ")
        def worst(t):
            vv = X["T"][t]; zz = np.array(vv["z"]); c = [(abs(100 * (vv["mean"][k] / vv["eq"][k] - 1)), names[k]) for k in range(len(names)) if lo[k] and abs(zz[k]) > 2]
            return max(c) if c else None
        s += (f"The other {int(lo.sum())} count small features per cell, which a million cells fix to a few parts in ten thousand, so they show how many rounds the fine texture needs. ")
        parts = []
        for t in Ts:
            w = worst(t)
            parts.append((f"with {t} rounds the worst is off by " if not parts else f"with {t} by ") + (f"{w[0]:.1f}%" if w[0] >= 0.095 else f"{w[0]:.2f}%") if w else
                         f"and with {t} none can be told apart (the largest difference is {np.abs(np.array(X['T'][t]['z'])[lo]).max():.1f} standard errors)")
        s += (", ".join(parts) + ". ").replace("w", "W", 1)
        if str(SWEEPS) in X["T"]: s += f"The zoom above uses {SWEEPS} rounds, to keep up with the camera."
        out.append(s + "</p>")
    Dp = load("zoom_deep.json"); deep = ""
    if Dp:
        deep = (f" Run {Dp['levels']} levels deep, a magnification of 10<sup>{Dp['levels'] * math.log10(2):.0f}</sup>, it shows no drift: the perimeter per cell averages "
                f"{Dp['first_half'][0]:.4f} over the first half of the levels and {Dp['second_half'][0]:.4f} over the second, against {R['sampled_local']['perim'][0]:.4f} in sampled shapes, "
                f"and the share of cells that are tips {Dp['first_half'][2]:.4f} and {Dp['second_half'][2]:.4f}, against {R['sampled_local']['deg1'][0]:.4f}.")
    out.append(f"<p>How much of the old picture survives a refinement? Very little at the scale of single squares, and nearly everything at large scales. After {SWEEPS} rounds, "
               f"{ch(4, SWEEPS):.0f}% of the boxes four old squares across have gained or lost their contents, {ch(8, SWEEPS):.0f}% of those eight across, and {ch(32, SWEEPS):.0f}% of those 32 across. "
               f"So the redrawing is local, though there is no sharp scale at which it stops.</p>")
    out.append("<p>The limits should be plain. Nobody has constructed the object that polyominoes of ever more cells converge to, if they do converge, and this does not construct it either. "
               "Each frame is a stand-in for what a sample of that many cells would look like around that point, and the stand-in has been tested as far as four doublings, because "
               "that is where the largest real samples are. Past that depth nothing has been compared with anything. The version in your browser works on a window around the point "
               "being approached, holding a border of doubled squares fixed to stand for the rest of the shape." + deep + "</p>")
    T["ZM_RESULTS"] = "\n  ".join(out)
    return D, T


def zoom_md():
    R, X = load("zoom_relax.json"), load("zoom_test.json"); md = []
    if not R: return ""
    sets = R["sets"]; ns = sorted(sets, key=int); G = growth_by_T(R, X)
    f = lambda v, d=4: f"{v[0]:.{d}f} ± {v[1]:.{d}f}"
    md += ["| sweeps | cells ÷ original, " + " | cells ÷ original, ".join(f"{int(n):,} doubled" for n in ns) + " | all measurements pooled | 2×2 blocks | perimeter | tips | branchings | boxes changed: 2 | 8 | 32 |",
           "|---:|" + "---:|" * (len(ns) + 8)]
    big = sets[ns[-1]]
    for i, t in enumerate(R["T"]):
        md.append(f"| {t} | " + " | ".join(f(sets[n]["growth"][i]) for n in ns) + f" | {f(G[t])} | " + " | ".join(f(big["local"][nm][i], 3) for nm in ("squares", "perim", "deg1", "deg3"))
                  + " | " + " | ".join(f"{100 * big['changed'][b][i][0]:.0f}%" for b in ("2", "8", "32")) + " |")
    md += ["", f"Prediction for the first columns: 2^(1/ν) = {R['pred']:.4f} ± {R['pred'] * math.log(2) / NU ** 2 * NU_ERR:.4f} for ν = {NU} ± {NU_ERR}. The local statistics are per cell, "
           f"relative to sampled shapes, for the {int(ns[-1]):,}-cell samples; boxes are in cells of the original."]
    cl = big["clumps"]
    md += ["", "| sweeps | " + " | ".join(f"clumps ≥ {s}" for s in (4, 8, 16, 32, 64)) + " |", "|---:|" + "---:|" * 5]
    for t in sorted(cl, key=int): md.append(f"| {t} | " + " | ".join(f"{cl[t][str(s)][0] / cl[t][str(s)][1]:.2f}" for s in (4, 8, 16, 32, 64)) + " |")
    md += ["", f"Loop-connected clumps of at least that many cells, relative to the equilibrium number (`analysis/clumps.json`), in {sets[ns[-1]]['samples']} refined samples."]
    if X and X["T"]:
        names = X["names"]; L = set(X["local"]); lo = np.array([nm in L for nm in names])
        md += ["", "| sweeps per level | shapes | cells | growth at levels 1 to 4 | per-cell counts (12): beyond 2 s.e. | largest | the other 48: beyond 2 s.e. | largest | largest differences |", "|---:|---:|---:|---|---:|---:|---:|---:|---|"]
        for t, v in sorted(X["T"].items(), key=lambda kv: int(kv[0])):
            z = np.array(v["z"]); o = np.argsort(-np.abs(z))[:3]
            md.append(f"| {t} | {v['shapes']} | {v['cells'][0]:,.0f} ± {v['cells'][1]:,.0f} | " + ", ".join(f"{v['growth'][k][0]:.4f}" for k in sorted(v["growth"])) +
                      f" | {int((np.abs(z[lo]) > 2).sum())} | {np.abs(z[lo]).max():.1f} | {int((np.abs(z[~lo]) > 2).sum())} | {np.abs(z[~lo]).max():.1f} | " +
                      ", ".join(f"{names[j]} {z[j]:+.1f}" + (f" ({100 * (v['mean'][j] / v['eq'][j] - 1):+.2f}%)" if abs(v["eq"][j]) > 8 * v["eq_err"][j] else "") for j in o) + " |")
        md += ["", f"Samples of 13,650 cells refined four times, against the equilibrium runs at 1,024,000 cells ({X['readings']} readings in {X['blocks']} blocks). "
               "A difference is in standard errors of the mean of that many equilibrium samples."]
    return "\n".join(md) + "\n"


if __name__ == "__main__":
    open(os.path.join(ROOT, "analysis", "zoom.md"), "w").write(zoom_md()); print(zoom_md())
