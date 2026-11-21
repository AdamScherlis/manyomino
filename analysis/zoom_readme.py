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
          f"Per level: {', '.join(f'{x:.1f}' for x in E['seconds'])} s; cells added as joins {', '.join(f'{x:.0f}' for x in E['links'])}; "
          f"parts in the finer world's box {', '.join(f'{x:.2f}' for x in E['fineParts'])}, kept apart in the new world {', '.join(f'{x:.2f}' for x in E['parts'])}."]
    if Pt:
        t += ["", f"Where strands rejoin (`analysis/part_test.py`, {Pt['roots']} roots): of the pieces into which the cells within {Pt['radii'][0]} of a root fall, the number still separate within R cells, "
              "in real samples and in the first world made from them (at R/2 of its own cells; in the box, touching cells of different parts are not joined).", "",
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
        q = Dr["quarters"]
        t += ["", f"Worlds only, {Dr['chains']} runs of {Dr['levels']} levels (`scripts/out_drift.sh`), means over the four quarters of the levels: cells of the world in its ring "
              + ", ".join(f"{m:.0f} ± {e:.0f}" for m, e in q["ring"]) + "; in its box " + ", ".join(f"{m:.0f} ± {e:.0f}" for m, e in q["box"]) + "; parts kept apart "
              + ", ".join(f"{m:.1f}" for m, e in q["parts"]) + ". Real samples have " + (f"{[r for r in E['rows'] if r[0] == 'L1 cells'][0][1] / 2.95:.0f}" if E else "about 114")
              + " cells in the ring at any scale (the real count 16 to 32 cells from a root, over 2.95). The first levels start above that (see the table) and the count falls towards it, "
              "by about 8% from the first quarter to the last; in blocks of 15 levels it wanders between 103 and 124, and whether it levels off these runs cannot say."]
    s2 = re.sub(r"<!-- zoom-out table -->.*?<!-- /zoom-out table -->", lambda m: "<!-- zoom-out table -->\n" + "\n".join(t) + "\n<!-- /zoom-out table -->", s2, flags=re.S)
open(p, "w").write(s2)
