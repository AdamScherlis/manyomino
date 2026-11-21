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
    B = int(P["B"]); T = {"ZO_N": f"{int(P['N']):,}", "ZO_BOX": str(2 * B), "ZO_REALBOX": str(4 * B), "ZO_REALR": str(2 * B - int(P["COLLAR"])), "ZO_TEXT": "", "ZO_PAIRS": ""}
    pm = lambda v, d=2: f"{v[0]:.{d}f}"
    out = []
    out.append(
        f"<p>That gives the surroundings at one scale. To go on without end, the view is handed from one such polyomino to another twice as coarse. Call each one a world: "
        f"a polyomino of {int(P['N']):,} squares that stands for the shape around the starting square, its squares standing for blocks 2, 4, 8&hellip; real squares on a side. "
        f"A new world is made from the last in two steps.</p>")
    out.append(
        f"<p>First it is halved. Every 2&nbsp;&times;&nbsp;2 block that holds anything becomes one square. The result is too fat, so a walk like the one used for zooming in "
        f"thins it, under two rules that keep it lined up with the finer world: no square may sit more than one step from where the finer world has squares, and no place where "
        f"the finer world has squares may be left with none on it or beside it. Then the middle of the halved world is kept, a box {2 * B} squares wide, and the rest is thrown "
        f"away. The sampler, with the box held fixed, draws everything else afresh, for {int(P['S'])} moves per square.</p>")
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
        f"as the real shape must.{part}</p>")
    out.append(
        f"<p>The picture takes from each world only a ring, from its box out to a box twice as wide. Inside the ring is the finer world&rsquo;s picture; outside it, the next "
        f"world&rsquo;s. A world&rsquo;s squares are large on the screen, so its ring is drawn finer by the doubling and walk of the last section, several times over, each time "
        f"with the finer picture held in place inside it. Where a strand crosses from one to the other the two are kept joined, and the walk may rework a strip "
        f"{int(P['COLLAR'])} squares deep just inside the seam, so that the two sides settle against each other.</p>")
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
            f"the real squares near that cell, and let it make {L} levels. At the resolution of single squares the result is a picture {2 * half} squares wide of which only the "
            f"middle {4 * B}&nbsp;&times;&nbsp;{4 * B} was given, {100 * (4 * B) ** 2 / (2 * half) ** 2:.1f}% of its area (and the outer {int(P['COLLAR'])} rows of that are reworked when the "
            f"first level is joined on). Set it beside the real squares around the same cell, "
            f"for {n} cells in {E.get('samples', 26)} samples.</p>")
        out.append("\n".join(t for t in tb if t))
        seam = [100 * r[3] for r in E["rows"] if "seam" in r[0]]; tex = [abs(100 * r[3]) for r in E["rows"] if r[0].endswith("(per cell)") and "floating" not in r[0] and not r[0].startswith("L1 ") and "4 neigh" not in r[0]]
        out.append(
            f"<p>Of {len(zs)} such numbers, {big} differ by more than two standard errors, where chance alone would give about {0.05 * len(zs):.0f}. So with {n} pictures of each "
            f"kind the two can be told apart, and the differences have a pattern. The made pictures hold {100 * tot[3]:.0f}% more squares overall. The two layers of squares on either "
            f"side of each seam hold {min(seam):.0f} to {max(seam):.0f}% too many, with as many more joins across. The fine texture is off by at most {max(tex):.1f}%. And a real "
            f"picture has no piece that neither holds the starting square nor reaches the edge, where the made ones have {rows['floating pieces'][2]:.1f} on average, small ones, "
            f"holding {100 * fl[2] / tot[2]:.1f}% of the squares: a part near a seam that failed to be joined on. None of it can be seen in a single picture.</p>")
    if Dp:
        rw = Dp["rows"]; fine = [r for r in rw if r["level"] > Dp["J"]]; h = len(fine) // 2; a, b = fine[:h], fine[h:]
        ms = lambda rs, k: (float(np.mean([r[k] for r in rs])), float(np.std([r[k] for r in rs], ddof=1) / np.sqrt(len(rs))))
        Dr = load("out_drift.json")
        drift = ""
        if Dr:
            q = Dr["quarters"]["ring"]; b_ = Dr["quarters"]["box"]
            drift = (f"In {Dr['chains']} runs of {Dr['levels']} levels each, a zoom of 2<sup>{Dr['levels']}</sup>, the ring of a world holds "
                     + ", ".join(f"{m:.0f}&nbsp;&plusmn;&nbsp;{e:.0f}" for m, e in q[:3]) + f" and {q[3][0]:.0f}&nbsp;&plusmn;&nbsp;{q[3][1]:.0f} squares on average over the four quarters of the levels, and its box "
                     + ", ".join(f"{m:.0f}" for m, e in b_[:3]) + f" and {b_[3][0]:.0f}. The amount varies a great deal from level to level, as the amount of shape near a point of a real polyomino does; real samples average "
                     + (f"{rows['L1 cells'][1] / 2.95:.0f}" if E else "114") + f" in the ring at any scale. So the amount falls by about {100 * (1 - q[3][0] / q[0][0]):.0f}% from the first quarter to "
                     "the last, towards the real figure; whether it levels off there these runs cannot say. ")
        out.append(
            f"<p>Beyond those levels there is nothing real to compare with, since no sample is large enough. What can be checked is that nothing drifts. {drift}"
            f"In the finished pictures of a run of {len(rw)} levels, once they are as fine as they get, there are {ms(a, 'empty')[0]:.3f} empty sites beside the shape per square in the first half and "
            f"{ms(b, 'empty')[0]:.3f} in the second, and {ms(a, 'tips')[0]:.3f} and {ms(b, 'tips')[0]:.3f} squares with a single neighbour; real samples have "
            f"{rows[f'L{L} empty sites next to the shape (per cell)'][1]:.3f} and {rows[f'L{L} with 1 neighbours (per cell)'][1]:.3f}. "
            f"A level takes about {ms(rw, 's')[0]:.0f} seconds to make on one core of a desktop machine, most of it the sampler.</p>" if E else "")
    out.append(
        f"<p>Three limits should be kept in mind. A world has {int(P['N']):,} squares, not infinitely many, which thins its outer parts; within the ring that is drawn the effect "
        f"is under a few per cent. The halving step is an approximation, tuned so that what it leaves on the edge of the box matches a typical polyomino. And each ring is drawn "
        f"finer by the zoom-in procedure, which is itself only tested, not proved. The sampling in between is exact.</p>")
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
