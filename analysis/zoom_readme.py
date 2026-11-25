"""Put the tables of analysis/zoom_report.py between the markers in README.md."""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zoom_report import zoom_md, load, ROOT
p = os.path.join(ROOT, "README.md"); s = open(p).read(); md = zoom_md()
open(os.path.join(ROOT, "analysis", "zoom.md"), "w").write(md)
s2 = re.sub(r"<!-- zoom tables -->.*?<!-- /zoom tables -->", lambda m: "<!-- zoom tables -->\n" + md + "<!-- /zoom tables -->", s, flags=re.S)
E = load("e2e_test.json"); Pt = load("part_test.json"); Dp = load("out_deep.json")
if E:
    t = [f"{E['roots']} roots in {E['samples']} samples, {E['levels']} levels. *real* and *made* are means over roots (per-cell rows: totals over roots divided), *s.e.* the difference in standard errors "
         "(paired by root; by resampling roots for the per-cell rows). Pieces are counted in the middle three quarters of the picture, so that its edge runs through the middle of the last ring.", "",
         "| | real | made | difference | s.e. |", "|---|---:|---:|---:|---:|"]
    g = lambda v: f"{v:,.0f}" if abs(v) >= 100 else f"{v:.2f}" if abs(v) >= 1 else f"{v:.4f}"
    for k, a_, b_, d, z in E["rows"]:
        if z != z: continue
        t.append(f"| {k} | {g(a_)} | {g(b_)} | {100 * d:+.1f}% | {z:+.1f} |")
    zs = [r[4] for r in E["rows"] if "floating" not in r[0] and r[4] == r[4]]
    t += ["", f"Leaving out the rows for floating pieces (which a real picture cannot have), {sum(1 for z in zs if abs(z) > 2)} of {len(zs)} rows are beyond 2 s.e. and the largest is {max(abs(z) for z in zs):.1f}. "
          f"Per level: {', '.join(f'{x:.1f}' for x in E['seconds'])} s; pieces joined again at the start of a stage (all stages of a level together) {', '.join(f'{x:.1f}' for x in E['loose'])}, with {', '.join(f'{x:.1f}' for x in E['links'])} cells; "
          + (f"of those joins, longer than 8 cells {', '.join(f'{x:.2f}' for x in E['long'])}; " if "long" in E else "") +
          f"parts in the finer world's box {', '.join(f'{x:.2f}' for x in E['fineParts'])}, kept apart in the new world {', '.join(f'{x:.2f}' for x in E['parts'])}."]
    Rf = load("ref_test.json")
    if Rf:
        rr = {r[0]: r for r in Rf["rows"]}; Lr = Rf["levels"]; cell = lambda k: f"{100 * rr[k][3]:+.1f}% ({rr[k][4]:+.1f})" if k in rr else ""
        zr = [r[4] for r in Rf["rows"] if "floating" not in r[0] and r[4] == r[4]]; fl = rr.get("floating pieces")
        t += ["", f"The refining alone (`analysis/e2e_test.py` with `N=0 S=0 ball=1024,250000` and the engine built for larger worlds, {Rf['roots']} roots): the real cells out to 1,024 from a root are "
              f"halved and thinned {Lr + 1} times and the rings refined again, with nothing drawn afresh, so that every world is the real shape at its scale and only the refining is on trial. "
              "Made against real, in per cent, with the difference in standard errors in brackets.", "",
              "| ring | " + " | ".join(str(k) for k in range(1, Lr + 1)) + " |", "|---|" + "---:|" * Lr,
              "| cells in the ring | " + " | ".join(cell(f"L{k} cells") for k in range(1, Lr + 1)) + " |",
              "| cells in the two rows inside its inner seam | " + " | ".join(cell(f"L{k} seam: cells within {1 if k == 1 else 2} inside") for k in range(1, Lr + 1)) + " |",
              "| cells in the two rows outside it | " + " | ".join(cell(f"L{k} seam: cells within {1 if k == 1 else 2} outside") for k in range(1, Lr + 1)) + " |",
              "| joins across it | " + " | ".join(cell(f"L{k} seam: joins across") for k in range(1, Lr + 1)) + " |", "",
              f"(One row either side at the first seam.) {sum(1 for z in zr if abs(z) > 2)} of the {len(zr)} rows of the full list are beyond 2 s.e., the largest {max(abs(z) for z in zr):.1f}"
              + (f"; floating pieces a picture {fl[2]:.3f}" if fl else "") + ". The version before this one, on 104 roots: joins across the five seams -2.5, +2.2, +9.3, -8.8 and -18.7%, "
              "the two rows outside the fifth seam 18.5% short and the fifth ring 5.5% over."]
    St = load("stage_test.json")
    if St:
        Z = St["zones"]
        t += ["", f"Stage by stage (`analysis/stage_test.py`, {St['roots']} roots, again with nothing drawn afresh): cells in the picture of a level at the end of each stage, divided by the true cells at that "
              "resolution, in zones named by layers of the level's own world. The strip is the two layers beyond the ring that are refined with it and not kept; it is thin at its far side, "
              "which is the edge of what is refined.", "",
              "| level | stage | " + " | ".join(Z) + " |", "|---:|---:|" + "---:|" * len(Z)]
        for r in St["rows"]: t.append(f"| {r['level']} | {r['stage']}{' (last)' if r['last'] else ''} | " + " | ".join(f"{r[k][0]:.3f} ± {r[k][1]:.3f}" for k in Z) + " |")
    Ms = load("mass_test.json")
    if Ms:
        t += ["", f"Far out (`analysis/mass_test.py`, {Ms['roots']} roots): real shapes counted in pixels {Ms['pixel']} cells wide, out to 128 cells from a root, and the same shapes halved and thinned "
              "until a cell is a pixel, half a pixel or a quarter of one, counted in the same pixels. Brightness is a pixel's count over the mean count of lit pixels.", "",
              "| cells to a pixel side | lit pixels | of the real number | spread of brightness (s.d. / mean) | brightness against real, pixels lit in both (r) | lit here only | lit in the real only |", "|---|---:|---:|---:|---:|---:|---:|"]
        for r in Ms["stages"]:
            k = r["stages finer than the pixel"]
            t.append(f"| {2 ** k}{' (real cells)' if r['real cells'] else ''} | {r['lit pixels'][0]:.0f} ± {r['lit pixels'][1]:.0f} | {100 * r['lit, of real']:.0f}% | {r['spread of brightness']:.2f} | "
                     + ("" if r["real cells"] or r["correlation with real"] is None or r["correlation with real"] != r["correlation with real"] else f"{r['correlation with real']:.2f}") + f" | {r['lit here only']:.0f} | {r['lit in the real only']:.0f} |")
    Wt = load("world_test.json")
    if Wt:
        L = Wt["layers"]; h = Wt["held to"]; bd = Wt["bands"]; f = lambda v: f"{v[0]:.1f} ± {v[1]:.1f}"; lo = 5
        t += ["", f"The step from one world to the next (`analysis/world_test.py`, {Wt['roots']} roots): cells by layer round the root (layer r is the square ring at distance r) in real shapes at their "
              f"own scale, in the same shapes halved, and in the first world made from them. Halving leaves one cell for {Wt['real cells per halved cell']:.3f} real ones. The sampler holds the layers below {h}.", "",
              "| layer | " + " | ".join(str(k) for k in range(lo, len(L["real"]))) + " |", "|---|" + "---:|" * (len(L["real"]) - lo)]
        for nm, k in (("real", "real"), ("halved", "halved"), ("world", "world")): t.append(f"| {nm} | " + " | ".join(f"{v:.2f}" for v in L[k][lo:]) + " |")
        t += ["", f"In the box: real {f(bd['box']['real'])}, halved {f(bd['box']['halved'])}. In the ring: real {f(bd['ring']['real'])}, halved {f(bd['ring']['halved'])}, world {f(bd['ring']['world'])} "
              f"(world less halved, root by root, {bd['ring']['world - halved'][0]:+.1f} ± {bd['ring']['world - halved'][1]:.1f}). In the first two layers drawn afresh: halved "
              f"{f(bd['first two layers drawn afresh']['halved'])}, world {f(bd['first two layers drawn afresh']['world'])} ({bd['first two layers drawn afresh']['world - halved'][0]:+.2f} ± "
              f"{bd['first two layers drawn afresh']['world - halved'][1]:.2f}). The finer world crosses the edge of the box at {Wt['crossings'][0]:.1f} sites of the halved one."]
    Tl = load("tell_test.json"); Tn = load("tell_null.json")
    if Tl:
        f = lambda k: f"{Tl[k]['auc']:.3f} ± {Tl[k]['se']:.3f} ({100 * Tl[k]['correct']:.0f}% of pictures called correctly)"
        t += ["", f"Classifiers trained to tell the made pictures from the real ones (`analysis/tell_test.py`; {Tl['pictures']} of each, scored on pictures from samples they were not trained on; "
              f"area under the ROC curve, 0.5 = cannot tell): a logistic regression on about 150 summary numbers per picture, {f('summary')}; the same without the numbers that describe "
              f"how the picture falls into pieces, {f('summary, without the pieces')}; a small network looking at 9×9 windows centred on cells, 3,000 per picture, {f('patches')}. "
              + (f"The same classifiers set to tell the real pictures of half the samples from those of the other half, where there is nothing to find, score "
                 f"{Tn['summary']['auc']:.3f} ± {Tn['summary']['se']:.3f}, {Tn['summary, without the pieces']['auc']:.3f} ± {Tn['summary, without the pieces']['se']:.3f} and "
                 f"{Tn['patches']['auc']:.3f} ± {Tn['patches']['se']:.3f}." if Tn else "")]
    if Pt:
        t += ["", f"Where strands rejoin (`analysis/part_test.py`, {Pt['roots']} roots): of the pieces into which the cells within {Pt['radii'][0]} of a root fall, the number still separate within R cells, "
              "in real samples and in the first world made from them (at R/2 of its own cells, and by the world's own rule for what is joined).", "",
              "| R | " + " | ".join(str(r) for r in Pt["radii"]) + " |", "|---|" + "---:|" * len(Pt["radii"]),
              "| real | " + " | ".join(f"{m:.2f} ± {e:.2f}" for m, e in Pt["real"]) + " |", "| made | " + " | ".join(f"{m:.2f} ± {e:.2f}" for m, e in Pt["made"]) + " |"]
    if Dp:
        import numpy as np
        rw = Dp["rows"]; h = len(rw) // 2; mean = lambda rs, k: float(np.mean([r[k] for r in rs]))
        t += ["", f"Run {len(rw)} levels out at the page's resolution (`analysis/out_deep.js`): first half of the levels against second half, cells of the world in its ring "
              f"{mean(rw[:h], 'ringWorld'):.0f} and {mean(rw[h:], 'ringWorld'):.0f}; cells in the finished picture of the ring {mean(rw[:h], 'cells'):,.0f} and {mean(rw[h:], 'cells'):,.0f}; "
              f"empty sites next to the shape per cell {mean(rw[:h], 'empty'):.4f} and {mean(rw[h:], 'empty'):.4f}; cells with one neighbour per cell {mean(rw[:h], 'tips'):.4f} and {mean(rw[h:], 'tips'):.4f}; "
              f"parts kept apart {mean(rw[:h], 'parts'):.2f} and {mean(rw[h:], 'parts'):.2f}. A level takes {mean(rw, 's'):.1f} s, {mean(rw, 'sWorld'):.1f} s of it the sampler. "
              "(The first six levels are left out of the picture numbers: there the picture is no finer than real cells.)" if False else
              f"Run {len(rw)} levels out at the page's resolution (`analysis/out_deep.js`): first half of the levels against second half, cells of the world in its ring "
              f"{mean(rw[:h], 'ringWorld'):.0f} and {mean(rw[h:], 'ringWorld'):.0f}; parts kept apart {mean(rw[:h], 'parts'):.2f} and {mean(rw[h:], 'parts'):.2f}; "
              f"and, from level {Dp['J'] + 1} on, where the picture is as fine as it gets, cells in the finished picture of the ring "
              f"{mean([r for r in rw[:h] if r['level'] > Dp['J']], 'cells'):,.0f} and {mean(rw[h:], 'cells'):,.0f}, empty sites next to the shape per cell "
              f"{mean([r for r in rw[:h] if r['level'] > Dp['J']], 'empty'):.4f} and {mean(rw[h:], 'empty'):.4f}, cells with one neighbour per cell "
              f"{mean([r for r in rw[:h] if r['level'] > Dp['J']], 'tips'):.4f} and {mean(rw[h:], 'tips'):.4f}. A level takes {mean(rw, 's'):.1f} s, {mean(rw, 'sWorld'):.1f} s of it the sampler."]
    Dr = load("out_drift.json")
    if Dr:
        REAL_BOX, REAL_RING = 61.11, 117.90      # cells of a real shape within 8, and from 8 to 16, of one of its cells (1,560 cells of 26 samples)
        w = Dr["what"]; f = lambda k, d=1: f"{w[k]['mean']:.{d}f} ± {w[k]['se']:.{d}f}"; g = lambda k, d=1: f"{w[k]['slope per 100 levels']:+.{d}f} ± {w[k]['slope se']:.{d}f}"
        t += ["", f"Worlds only, {Dr['chains']} runs of {Dr['levels']} levels from cells of {Dr['chains']} different samples (`scripts/out_drift.sh`). Means over all levels: cells of the world in "
              f"its box {f('box')}, in its ring {f('ring')}, parts kept apart {f('parts', 2)}; real shapes have {REAL_BOX:.1f} and {REAL_RING:.1f} cells within 8 and from 8 to 16 cells of one of "
              f"their cells (1,560 cells of 26 samples). A straight line through each run changes by {g('box')} cells per 100 levels in the box and {g('ring')} in the ring. "
              f"Quarter by quarter the ring holds " + ", ".join(f"{v:.0f}" for v in w["ring"]["quarters"]) + "."]
    s2 = re.sub(r"<!-- zoom-out table -->.*?<!-- /zoom-out table -->", lambda m: "<!-- zoom-out table -->\n" + "\n".join(t) + "\n<!-- /zoom-out table -->", s2, flags=re.S)
open(p, "w").write(s2)
