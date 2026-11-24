"""Text, numbers and the picture for the zoom-out section of the page.

Reads analysis/e2e_test.json and analysis/e2e_pics.npz (analysis/e2e_test.py ... final),
analysis/part_test.json (analysis/part_test.py) and analysis/out_deep.json (analysis/out_deep.js).
"""
import json, os
import numpy as np
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def load(name):
    p = os.path.join(ROOT, "analysis", name)
    return json.load(open(p)) if os.path.exists(p) else None


def settings():
    """the constants of page/outward.js"""
    import re
    s = open(os.path.join(ROOT, "page", "outward.js")).read(); m = re.search(r"const OUT = \{([^}]*)\}", s)
    return {k.strip(): float(v) for k, v in (kv.split(":") for kv in m.group(1).split(","))}


def pairs_png(path, n=6, crop=192):
    """n panels, half of them real and half made (never both of one cell), in an order fixed by a seed; returns which are real"""
    from PIL import Image
    p = os.path.join(ROOT, "analysis", "e2e_pics.npz")
    if not os.path.exists(p): return None
    d = np.load(p); real = np.unpackbits(d["real"], axis=-1).astype(bool); made = np.unpackbits(d["made"], axis=-1).astype(bool); R = real.shape[1] // 2
    rng = np.random.default_rng(7); which = rng.permutation(n) < n // 2; order = rng.permutation(min(len(real), len(made)))[:n]
    from PIL import ImageDraw
    bg, ink, gap = (13, 15, 22), (233, 235, 242), 8; S = 2 * crop; W = 3 * S + 2 * gap; H = (n // 3) * S + (n // 3 - 1) * gap
    im = np.zeros((H, W, 3), np.uint8); im[...] = (70, 76, 92)
    for k in range(n):
        g = (real if which[k] else made)[order[k]][R - crop:R + crop, R - crop:R + crop]; y0 = (k // 3) * (S + gap); x0 = (k % 3) * (S + gap)
        tile = np.where(g[..., None], np.array(ink, np.uint8), np.array(bg, np.uint8)); im[y0:y0 + S, x0:x0 + S] = tile[::-1]
    img = Image.fromarray(im); dr = ImageDraw.Draw(img)
    for k in range(n):                                                         # a letter in the corner of each
        y0 = (k // 3) * (S + gap); x0 = (k % 3) * (S + gap); dr.rectangle([x0, y0, x0 + 26, y0 + 26], fill=(70, 76, 92))
        t = Image.new("L", (8, 12), 0); ImageDraw.Draw(t).text((1, 0), "abcdefgh"[k], fill=255); t = t.resize((16, 24), Image.NEAREST); img.paste(Image.new("RGB", t.size, ink), (x0 + 6, y0 + 1), t)
    img.save(path, optimize=True)
    return [bool(w) for w in which]


def out_report(assets=None):
    P = settings(); E = load("e2e_test.json"); Pt = load("part_test.json"); Dp = load("out_deep.json")
    B = int(P["B"]); T = {"ZO_N": f"{int(P['N']):,}", "ZO_BOX": str(2 * B), "ZO_REALBOX": str(4 * B), "ZO_REALR": str(2 * B - min(int(P["COLLAR"]), B // 2)), "ZO_TEXT": "", "ZO_PAIRS": ""}
    pm = lambda v, d=2: f"{v[0]:.{d}f}"
    out = []
    out.append(
        f"<p>That gives the surroundings at one scale. To go on without end, the view is handed from one such polyomino to another twice as coarse. Call each one a world: "
        f"a polyomino of {int(P['N']):,} squares that stands for the shape around the starting square, its squares standing for blocks 2, 4, 8&hellip; real squares on a side. "
        f"A new world is made from the last in two steps.</p>")
    out.append(
        f"<p>First it is halved. Every 2&nbsp;&times;&nbsp;2 block that holds anything becomes one square. The result is too fat, so a walk like the one used for zooming in "
        f"thins it, under two rules that keep it lined up with the finer world: no square may sit more than one step from where the finer world has squares, and no place where "
        f"the finer world has squares may be left with none on it or beside it. How hard the walk thins is fixed by one requirement: a polyomino seen at half the scale has "
        f"1/2.95 as many squares, so halving must divide the count by 2.95. Then the middle of the halved world is kept, a box {2 * B} squares wide together with the "
        f"{int(P['HOLD'])} rows around it, and the rest is thrown away. The sampler, with that much held fixed, draws everything else afresh, for {int(P['S'])} moves per square.</p>")
    Wt = load("world_test.json")
    if Wt:
        bd = Wt["bands"]; fr = bd["first two layers drawn afresh"]
        out.append(
            f"<p>Two things can go wrong in that step, and both can be measured ({Wt['roots']} squares of real samples). The halved shape should be an ordinary polyomino at its own "
            f"scale: a real one has {bd['box']['real'][0]:.0f} squares within {B} of one of its squares and {bd['ring']['real'][0]:.0f} between {B} and {2 * B}; the same shapes halved have "
            f"{bd['box']['halved'][0]:.0f} and {bd['ring']['halved'][0]:.0f}. And what the sampler draws should be what was there: in the first two rows it draws afresh the halved real shape has "
            f"{fr['halved'][0]:.1f} squares and the sampled world {fr['world'][0]:.1f} (the difference, shape by shape, is {fr['world - halved'][0]:+.1f}&nbsp;&plusmn;&nbsp;{fr['world - halved'][1]:.1f}).</p>")
    part = ""
    if Pt:
        r, m, R = Pt["real"], Pt["made"], Pt["radii"]; i15 = R.index(3 * B); i2 = R.index(4 * B)
        part = (f" Real shapes do this quickly, and so do the worlds. In real samples the squares within {R[0]} of a given square fall into {pm(r[0])} separate parts on average; "
                f"within {R[i15]} squares {pm(r[i15])} of them are still separate, and within {R[i2]}, {pm(r[i2])}. The first world made from the same squares starts with "
                f"{pm(m[0])} parts (a part of only two or three squares is not kept apart from what it touches) and has {pm(m[i15])} and {pm(m[i2])} at the matching distances "
                f"({Pt['roots']} cells).")
    out.append(
        f"<p>One thing needs care. Two strands that pass within a square of each other without meeting come to touch when the world is halved. If the sampler took them for "
        f"joined it would see no need to join them anywhere else, and the picture would fall into more separate pieces than a real polyomino shows. So each square of the box "
        f"remembers which part of the finer world it stands for, and squares of different parts that touch do not count as joined. The sampler has to join them outside the box, "
        f"as the real shape must. For the same reason a square of the box counts as joined to a square outside it only where the finer world really crosses the edge of the box "
        f"between them.{part}</p>")
    out.append(
        f"<p>The picture takes from each world only a ring, from its box out to a box twice as wide. Inside the ring is the finer world&rsquo;s picture; outside it, the next "
        f"world&rsquo;s. A world&rsquo;s squares are large on the screen, so its ring is drawn finer by the doubling and walk of the last section, several times over, each time "
        f"with the finer picture held in place inside it. The two are joined where the finer picture&rsquo;s strands leave its edge. One square is kept in place there at every "
        f"stage, on each side, and each is always among the four squares its predecessor became, so the strand and the ring meet at every stage with nothing to patch. Until "
        f"the last stage those are the only joins; anything else that touches across the seam merely touches. The walk may also rework a strip up to {int(P['COLLAR'])} squares "
        f"deep just inside the seam, so that the two sides settle against each other: a ring made with nothing outside it comes out thin and short of loose ends along its "
        f"outer edge.</p>")
    if E:
        rows = {r[0]: r for r in E["rows"]}; L = E["levels"]; half = B * 2 ** (L + 1); n = E["roots"]
        def row(key, label, d=1, pc=True):
            if key not in rows: return ""
            _, a, b, diff, z = rows[key]; f = (lambda v: f"{v:,.0f}") if abs(a) >= 100 else (lambda v: f"{v:.{d + 1 if abs(a) < 1 else d}f}" if abs(a) >= 0.1 else f"{v:.4f}")
            return f"<tr><td>{label}</td><td>{f(a)}</td><td>{f(b)}</td><td>{100 * diff:+.1f}%</td><td>{z:+.1f}</td></tr>"
        tb = ['<div class="scroll"><table class="data"><thead><tr><th></th><th>real</th><th>made</th><th>difference</th><th>in standard errors</th></tr></thead><tbody>']
        for k in range(1, L + 1): tb.append(row(f"L{k} cells", f"squares {B * 2 ** k}&ndash;{B * 2 ** (k + 1)} from the middle"))
        tb.append(row(f"L{L} empty sites next to the shape (per cell)", "empty sites beside the shape, per square (outer ring)", 3))
        tb.append(row(f"L{L} with 1 neighbours (per cell)", "squares with one neighbour, per square (outer ring)", 3))
        tb.append(row(f"L{L} full 2x2 (per cell)", "2&nbsp;&times;&nbsp;2 blocks, per square (outer ring)", 3))
        tb.append(row("L1 seam: joins across", f"joins across the first seam ({2 * B} squares out)", 2))
        tb.append(row(f"L{L} seam: joins across", f"joins across the last seam ({B * 2 ** L} squares out)", 1))
        tb.append(row("pieces", "separate pieces in the middle three quarters", 1))
        tb.append(row("share in the largest piece", "share of squares in the largest piece", 3))
        tb.append("</tbody></table></div>")
        zs = [r[4] for r in E["rows"] if "floating" not in r[0] and r[4] == r[4]]; big = sum(1 for z in zs if abs(z) > 2)
        fl = rows.get("cells in floating pieces"); tot = rows.get("cells")
        out.append(
            f"<p>How close is it? For the first few levels it can be set against the real thing. Take a real sample of 1,024,000 cells and one of its cells. Give the procedure "
            f"the real squares near that cell, and let it make {L + 1} levels. At the resolution of single squares the first {L} are a picture {2 * half} squares wide of which only the "
            f"middle {4 * B}&nbsp;&times;&nbsp;{4 * B} was given, {100 * (4 * B) ** 2 / (2 * half) ** 2:.1f}% of its area (and the outer {min(int(P['COLLAR']), B // 2)} rows of that are reworked when the "
            f"first level is joined on). Set it beside the real squares around the same cell, "
            f"for {n} cells in {E.get('samples', 26)} samples.</p>")
        out.append("\n".join(t for t in tb if t))
        seam = [(100 * r[3], r[4]) for r in E["rows"] if "seam" in r[0]]; ring = [abs(100 * rows[f"L{k} cells"][3]) for k in range(1, L + 1)]
        tex = [abs(100 * r[3]) for r in E["rows"] if r[0].endswith("(per cell)") and "floating" not in r[0] and not r[0].startswith("L1 ") and "4 neigh" not in r[0]]
        nfl = round(rows["floating pieces"][2] * n); atseam = sum(1 for d, z in seam if abs(z) > 2)
        near = [(100 * r[3], r[4]) for r in E["rows"] if "seam" in r[0] and r[0][1] in "12"]; far = [(100 * r[3], r[4]) for r in E["rows"] if "seam" in r[0] and r[0][1] not in "12"]
        out.append(
            f"<p>Of {len(zs)} such numbers, {big} differ by more than two standard errors, where chance alone would give about {0.05 * len(zs):.0f}. The seams are where an earlier "
            f"version went wrong: it had 4 to 11% too many squares on either side of every seam, with as many more joins across. Now, at the first two seams the rows of squares "
            f"on either side and the joins across are within {max(abs(d) for d, z in near):.0f}% of the real ones, none of them more than {max(abs(z) for d, z in near):.1f} standard "
            f"errors away. At the deeper seams they fall short, by {min(abs(d) for d, z in far):.0f} to {max(abs(d) for d, z in far):.0f}%, and {sum(1 for d, z in far if abs(z) > 2)} "
            f"of those {len(far)} numbers are beyond two standard errors: a ring that has been drawn finer three or four times meets the picture inside it at too few places. The "
            f"number of squares in each ring is within {max(ring):.0f}% of the real one and the fine texture within {max(tex):.1f}%. In {n} made pictures there "
            f"{'was' if nfl == 1 else 'were'} {nfl if nfl else 'no'} piece{'' if nfl == 1 else 's'} that neither held the starting square nor reached the edge, which a real picture "
            f"cannot have{' (a stray pair of squares each)' if nfl == 2 else ''}.</p>")
    Tl = load("tell_test.json"); Tn = load("tell_null.json")
    if Tl and E:
        pc = lambda k: f"{100 * Tl[k]['correct']:.0f}%"
        out.append(
            f"<p>A list of numbers only finds the differences someone thought to list. So, as for the samples themselves, two classifiers were trained to tell made pictures from "
            f"real ones and then scored on pictures from samples they had not seen. One looks at about 150 summary numbers of a whole picture, the other at 3,000 small windows of "
            f"it, 9 squares wide. A classifier that has found nothing is right half the time. The first calls {pc('summary')} of pictures correctly, the second {pc('patches')}." +
            (f" Set to tell the real pictures of half the samples from those of the other half, where there is nothing to find, they score {100 * Tn['summary']['correct']:.0f}% and "
             f"{100 * Tn['patches']['correct']:.0f}%." if Tn else "") + (" With the earlier version the first classifier was right 59% of the time, and what it had found was the seams. It still does a little better than chance, "
             "and the shortfall at the deeper seams is the likeliest reason.") + "</p>")
    if Dp:
        rw = Dp["rows"]; fine = [r for r in rw if r["level"] > Dp["J"]]; h = len(fine) // 2; a, b = fine[:h], fine[h:]
        ms = lambda rs, k: (float(np.mean([r[k] for r in rs])), float(np.std([r[k] for r in rs], ddof=1) / np.sqrt(len(rs))))
        Dr = load("out_drift.json")
        drift = ""
        if Dr:
            w = Dr["what"]; pm_ = lambda k: f"{w[k]['mean']:.0f}&nbsp;&plusmn;&nbsp;{w[k]['se']:.0f}"
            drift = (f"In {Dr['chains']} runs of {Dr['levels']} levels each, a zoom of 2<sup>{Dr['levels']}</sup>, started from cells of {Dr['chains']} different samples, a world holds "
                     f"{pm_('box')} squares in its box and {pm_('ring')} in its ring on average; a real polyomino has 61 and 118 within the same distances of one of its squares. "
                     f"The amount varies a great deal from level to level, as it does from place to place in a real polyomino, but it does not creep: a straight line fitted through each run "
                     f"changes by {w['ring']['slope per 100 levels']:+.1f}&nbsp;&plusmn;&nbsp;{w['ring']['slope se']:.1f} squares in the ring per hundred levels. ")
        out.append(
            f"<p>Beyond those levels there is nothing real to compare with, since no sample is large enough. What can be checked is that nothing drifts. {drift}"
            f"In the finished pictures of a run of {len(rw)} levels, once they are as fine as they get, there are {ms(a, 'empty')[0]:.3f} empty sites beside the shape per square in the first half and "
            f"{ms(b, 'empty')[0]:.3f} in the second, and {ms(a, 'tips')[0]:.3f} and {ms(b, 'tips')[0]:.3f} squares with a single neighbour; real samples have "
            f"{rows[f'L{L} empty sites next to the shape (per cell)'][1]:.3f} and {rows[f'L{L} with 1 neighbours (per cell)'][1]:.3f}. "
            f"A level takes about {ms(rw, 's')[0]:.0f} seconds to make on one core of a desktop machine, {ms(rw, 'sWorld')[0]:.1f} of them the sampler.</p>" if E else "")
    out.append(
        f"<p>Three limits should be kept in mind. A world has {int(P['N']):,} squares, not infinitely many, which thins its outer parts; within the ring that is drawn the effect "
        f"is under a few per cent. The halving step is an approximation: a shape halved and thinned is not exactly a polyomino drawn at random, only one with the right amount "
        f"of shape in the right places as far as these counts can tell. And each ring is drawn finer by the zoom-in procedure, which is itself only tested, not proved. Two "
        f"numbers in all this are set rather than derived, how hard the halving walk thins and how readily the refining walk adds, both by the one requirement that a "
        f"polyomino at twice the scale has 2.95 times the squares. The sampling in between is exact.</p>")
    T["ZO_TEXT"] = "\n  ".join(out)
    if assets:
        which = pairs_png(os.path.join(assets, "outpairs.png"))
        if which:
            names = "abcdef"; real = ", ".join(names[i] for i, w in enumerate(which) if w)
            T["ZO_PAIRS"] = (
                '<section>\n  <figure>\n    <div class="plate"><img src="explain/outpairs.png" alt="Six pictures of pieces of polyomino, three real and three made by the zoom out" loading="lazy" '
                'style="display:block;width:100%;height:auto"></div>\n'
                f'    <figcaption>Six pictures, each {2 * 192} squares wide, lettered a to f along the rows. Three are cut from real samples of a million cells. In the other three only '
                f'the {4 * B}&nbsp;&times;&nbsp;{4 * B} squares in the middle were given, and the rest was made from them by the procedure above. '
                f'<details style="display:inline"><summary style="display:inline;cursor:pointer;text-decoration:underline">Which are real?</summary> {real}.</details></figcaption>\n  </figure>\n</section>')
    return T


if __name__ == "__main__":
    t = out_report(); print(t["ZO_TEXT"][:3000])
