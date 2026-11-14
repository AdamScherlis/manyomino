"""Build the explainer page from the production runs.

  python3 analysis/build_page.py SITE_DIR [ARTIFACT_DIR]

Writes SITE_DIR/explain.html plus SITE_DIR/explain/* (figures and the hero sample), and
optionally a body-only copy of the page for publishing as a standalone artifact.
"""
import glob, json, math, os, re, sys, datetime
import numpy as np
import summarize
from poly import load_bin, neighbours, centre, bfs, render_box, encode_tree

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
PAGE = os.path.join(ROOT, "page")
PROD = os.path.join(ROOT, "prod")
P2 = os.path.join(ROOT, "prod2")
LAMBDA, AMP = 4.0625696, 0.316915
# <Rg> published in the earlier gallery (static/manyomino/index.html, Aug 2026)
GALLERY = [[1000, 36.8], [2000, 54.1], [3000, 72.7], [5000, 101.5], [10000, 158.6], [30000, 321.4],
           [100000, 640.5], [400000, 1756], [1600000, 3716]]
EXAMPLE = json.load(open(os.path.join(PAGE, "example.json")))


def wfit(x, y, e):
    """weighted least squares y = a + b x; returns a, b, err_b, chi2/dof"""
    w = 1 / np.asarray(e) ** 2
    x, y = np.asarray(x), np.asarray(y)
    S, Sx, Sy, Sxx, Sxy = w.sum(), (w * x).sum(), (w * y).sum(), (w * x * x).sum(), (w * x * y).sum()
    d = S * Sxx - Sx ** 2
    a, b = (Sxx * Sy - Sx * Sxy) / d, (S * Sxy - Sx * Sy) / d
    chi2 = (w * (y - a - b * x) ** 2).sum() / max(1, len(x) - 2)
    return a, b, math.sqrt(S / d), chi2


def commas(v):
    return f"{int(round(v)):,}"


def sci(v, digits=1):
    e = int(math.floor(math.log10(v)))
    m = v / 10 ** e
    return (f"{m:.{digits}f}&thinsp;&times;&thinsp;" if digits else "") + f"10<sup>{e}</sup>"


def duration(sec):
    if sec < 1.5: return "about a second"
    if sec < 90: return f"about {sec:.0f} seconds"
    if sec < 5400: return f"about {sec / 60:.0f} minutes"
    if sec < 2 * 86400: return f"about {sec / 3600:.0f} hours"
    if sec < 2 * 365 * 86400: return f"about {sec / 86400:.0f} days"
    return f"about {commas(sec / 86400 / 365.25)} years"


def latest_ckpt(n, tag="A"):
    fn = f"{PROD}/ckpt_{tag}_n{n}.bin"
    return fn if os.path.exists(fn) else None


def prep(fn):
    x, y = load_bin(fn)
    nbr = neighbours(x, y)
    c, diam = centre(nbr)
    d, _ = bfs(nbr, c)
    return x, y, nbr, c, d, diam


def main():
    site = sys.argv[1]
    art = sys.argv[2] if len(sys.argv) > 2 else None
    assets = os.path.join(site, "explain")
    os.makedirs(assets, exist_ok=True)
    D = summarize.main()
    sizes, reloc = D["sizes"], D["reloc"]
    T = {}

    # ---- hero: the largest finished ladder size
    hero_n = max(r["n"] for r in sizes if latest_ckpt(r["n"]) and r["neff"] >= 12)
    hx, hy, hnbr, hc, hd, hdiam = prep(latest_ckpt(hero_n))
    enc, _ = encode_tree(hx, hy, hnbr, hc)
    import base64
    open(os.path.join(assets, "hero.txt"), "w").write(base64.b64encode(enc).decode())
    T["NHERO"] = commas(hero_n)
    T["EXP"] = commas(hero_n * math.log10(LAMBDA) + math.log10(AMP) - math.log10(hero_n))
    T["HERO_DIAM"], T["HERO_W"], T["HERO_H"] = commas(hdiam), commas(hx.max() + 1), commas(hy.max() + 1)

    # ---- local composition, from the largest size
    big = max(sizes, key=lambda r: r["n"])
    n = big["n"]
    T["CUT_PCT"] = f"{100 * (1 - big['noncut'] / n):.1f}"
    T["TIP_PCT"] = f"{100 * big['deg1'] / n:.1f}"
    T["PERIM"] = f"{big['perim'] / n:.3f}"
    T["HOLE_PER"] = f"{n / big['holes']:.0f}"

    # ---- scale ladder: four sizes, x16 apart, ending at the hero
    lad = [hero_n // 4096, hero_n // 256, hero_n // 16, hero_n]
    if not latest_ckpt(lad[0]):      # partial data: fall back to x4 spacing
        lad = [hero_n // 64, hero_n // 16, hero_n // 4, hero_n]
    html = ""
    for i, ln in enumerate(lad):
        tag = "B" if i < 3 and latest_ckpt(ln, "B") else "A"
        x, y, nbr, c, d, diam = (hx, hy, hnbr, hc, hd, hdiam) if ln == hero_n else prep(latest_ckpt(ln, tag))
        render_box(x, y, d, size=520).save(os.path.join(assets, f"ladder{i}.png"), optimize=True)
        html += f'<div><div class="plate"><img src="explain/ladder{i}.png" alt="uniform random polyomino with {ln:,} cells" loading="lazy" width="520" height="520"></div><p><b>{ln:,}</b> cells</p></div>\n'
    T["LADDER"] = html

    # ---- mixing-time fits
    rl = [r for r in reloc if r["n"] >= 50]
    a, z, ze, _ = wfit(np.log([r["n"] for r in rl]), np.log([r["tau_rg2"] for r in rl]), [max(r["tau_rg2_err"] / r["tau_rg2"], 0.03) for r in rl])
    D["fit_old"] = {"a": math.exp(a), "z": z, "z_err": ze}
    nw = [r for r in sizes if r["n"] >= 1000]
    a2, z2, ze2, _ = wfit(np.log([r["n"] for r in nw]), np.log([r["tau_rg2"] for r in nw]), [max(r["tau_rg2_err"] / r["tau_rg2"], 0.03) for r in nw])
    D["fit_new"] = {"a": math.exp(a2), "z": z2, "z_err": ze2}
    T["Z_OLD"], T["Z_NEW"] = f"{z:.1f}", f"{z2:.1f}"
    T["OLD_NMIN"], T["OLD_NMAX"] = commas(min(r["n"] for r in rl)), commas(max(r["n"] for r in rl))
    tau_old_hero = math.exp(a) * hero_n ** z
    T["TAU_OLD_HERO"] = sci(tau_old_hero, 0)
    T["OLD_RATE"] = "million"
    T["OLD_YEARS"] = duration(tau_old_hero / 1e6)
    T["TAU_NEW_COEF"] = f"{np.mean([r['tau_rg2'] / r['n'] for r in nw]):.0f}"
    bench = {}
    bf = os.path.join(ROOT, "analysis", "bench.json")
    if os.path.exists(bf):
        bench = {int(k): v for k, v in json.load(open(bf)).items()}
    spm = {r["n"]: (1 / bench[r["n"]] if r["n"] in bench else r["sec_per_move"]) for r in sizes}
    cs = [r for r in sizes if r["n"] >= 4000]
    _, ce, _, _ = wfit(np.log([r["n"] for r in cs]), np.log([spm[r["n"]] for r in cs]), [0.1] * len(cs))
    T["COST_EXP"] = f"{ce:.2f}"
    _, cpe, _, _ = wfit(np.log([r["n"] for r in cs]), np.log([spm[r["n"]] * r["tau_rg2"] for r in cs]), [0.1] * len(cs))
    T["CPU_EXP"] = f"{cpe:.1f}"
    hero_row = next(r for r in sizes if r["n"] == hero_n)
    T["CPU_HERO"] = duration(2 * hero_row["tau_rg2"] * spm[hero_n])
    T["NEFF_HERO"] = commas(hero_row['moves'] / hero_row['tau_rg2'])
    D["cpu"] = [[r["n"], spm[r["n"]], r["tau_rg2"] * spm[r["n"]]] for r in sizes]

    # ---- nu
    NU_NMIN = 4000
    fs = [r for r in sizes if r["n"] >= NU_NMIN]
    a3, nu, nue, chi = wfit(np.log([r["n"] for r in fs]), np.log([r["rg2"] for r in fs]), [r["rg2_err"] / r["rg2"] for r in fs])
    a4, nu_rg, _, _ = wfit(np.log([r["n"] for r in fs]), np.log([r["rg"] for r in fs]), [r["rg_err"] / r["rg"] for r in fs])
    D["fit_nu"] = {"a": math.exp(a4), "nu": nu_rg, "nu2": nu / 2, "err": nue / 2, "chi2": chi}
    T["NU"], T["NU_ERR"], T["NU_NMIN"] = f"{nu / 2:.4f}", f"{nue / 2 * max(1, math.sqrt(chi)):.4f}", commas(NU_NMIN)

    # ---- lambda
    LAM_NMIN = 1000
    ls = [r for r in sizes if r["n"] >= LAM_NMIN]
    w = np.array([1 / r["lam_err"] ** 2 for r in ls])
    lam = float((w * np.array([r["lam"] for r in ls])).sum() / w.sum())
    T["LAM"], T["LAM_ERR"], T["LAM_NMIN"] = f"{lam:.4f}", f"{1 / math.sqrt(w.sum()):.4f}", commas(LAM_NMIN)
    D["lam_pooled"] = [lam, 1 / math.sqrt(w.sum())]

    # ---- agreement between independent ladders
    zs = []
    for r in sizes:
        if len(r["chain_rg2"]) >= 2:
            (m1, e1), (m2, e2) = r["chain_rg2"][:2]
            zs.append(abs(m1 - m2) / math.hypot(e1, e2))
    T["AGREE"] = (f"Two ladders of runs from different seeds were carried through {len(zs)} sizes. Their mean squared radii differ by "
                  f"{np.mean(zs):.1f} standard errors on average and {max(zs):.1f} at worst.")

    # ---- gallery comparison
    gl = []
    for gn, grg in GALLERY:
        pred = math.exp(a4) * gn ** nu_rg
        gl.append((gn, grg, grg / pred))
    D["gallery"] = GALLERY
    D["gallery_ratio"] = gl
    T["GALLERY_VERDICT"] = open(os.path.join(PAGE, "verdict.html")).read().strip() if os.path.exists(os.path.join(PAGE, "verdict.html")) else ""

    # ---- square head-to-head
    sq_moves = None
    for tag in ("old", "new"):
        fn = f"{P2}/snap_{tag}sq_matched.bin"
        if os.path.exists(fn):
            sq_moves = json.load(open(f"{P2}/matched.json"))["moves"]
    strip = ""
    if sq_moves:
        xs = {t: prep(f"{P2}/snap_{t}sq_matched.bin") for t in ("old", "new")}
        span = max(max(v[0].max(), v[1].max()) + 1 for v in xs.values()) / 0.92
        for t, lab in (("old", "One cell at a time"), ("new", "With branch transplants")):
            x, y, nbr, c, d, diam = xs[t]
            render_box(x, y, d, size=640, span=span).save(os.path.join(assets, f"square_{t}.png"), optimize=True)
            rg = math.sqrt(((x - x.mean()) ** 2 + (y - y.mean()) ** 2).mean())
            strip += f'<div><div class="plate"><img src="explain/square_{t}.png" alt="30,000-cell polyomino after {sq_moves:,} moves: {lab.lower()}" loading="lazy" width="640" height="640"></div><p>{lab} &middot; radius {rg:.0f}</p></div>\n'
        T["SQ_MOVES"] = sci(sq_moves, 1)
    T["SQUARE_STRIP"] = strip
    T.setdefault("SQ_MOVES", "")

    demo = next(r for r in sizes if r["n"] == 1000)
    D["demo"] = {"n": 1000, "rg": demo["rg"], "rg_sd": demo["rg_sd"]}
    D["example"] = EXAMPLE
    D["heroUrl"] = "explain/hero.txt"
    D.setdefault("square_old", [[1e5, 5000]]); D.setdefault("square_new", [[1e5, 5000]]); D.setdefault("square_eq", {"rg": 320, "rg_sd": 60})
    T["DATE"] = datetime.date.today().strftime("%B %Y")
    json.dump(D, open(os.path.join(ROOT, "analysis", "data.json"), "w"), indent=1)

    body = open(os.path.join(PAGE, "body.html")).read()
    slim = {k: v for k, v in D.items() if k not in ("gallery_ratio",)}
    for r in slim["sizes"]:
        for k in list(r):
            if isinstance(r[k], float): r[k] = float(f"{r[k]:.7g}")
    js = "const DATA = " + json.dumps(slim, separators=(",", ":")) + ";\n"
    for f in ("charts.js", "viewer.js", "sampler.js", "app.js"):
        js += open(os.path.join(PAGE, f)).read() + "\n"
    cssText = open(os.path.join(PAGE, "page.css")).read()
    title = "Manyomino"

    def fill(b, gallery_href):
        t = dict(T, GALLERY_HREF=gallery_href)
        b = re.sub(r"\{\{(\w+)\}\}", lambda m: str(t[m.group(1)]), b)
        return b

    full = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title>
<meta name="description" content="Uniformly random polyominoes with a million cells, and the Markov chain that makes them.">
<style>
{cssText}</style>
</head>
<body>
{fill(body, "index.html")}
<script>
{js}</script>
</body>
</html>
"""
    open(os.path.join(site, "explain.html"), "w").write(full)
    if art:
        os.makedirs(art, exist_ok=True)
        frag = f"""<title>{title}</title>
<meta name="claude:viewport-fit-sides" content="cover">
<style>
{cssText}</style>
{fill(body, "https://adam.scherl.is/static/manyomino/")}
<script>
{js}</script>
"""
        open(os.path.join(art, "manyomino.html"), "w").write(frag)
    # ---- results table for the README
    md = ["| cells | chains | run length / tau | <Rg^2> | <Rg> | tau / n | perimeter / n | removable / n | holes / n | lambda estimate |",
          "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in sizes:
        md.append(f"| {r['n']:,} | {r['chains']} | {r['moves'] / r['tau_rg2']:,.0f} | {r['rg2']:.5g} ± {r['rg2_err']:.2g} | {r['rg']:.5g} ± {r['rg_err']:.2g} | "
                  f"{r['tau_rg2'] / r['n']:.1f} | {r['perim'] / r['n']:.5f} | {r['noncut'] / r['n']:.5f} | {r['holes'] / r['n']:.5f} | {r['lam']:.4f} ± {r['lam_err']:.4f} |")
    md += ["", f"Relocate-only autocorrelation time (moves), n = " + ", ".join(f"{r['n']}: {r['tau_rg2']:.3g}" for r in rl) + f"; fit tau ∝ n^{z:.2f} (±{ze:.2f}).",
           f"With transplants: tau ∝ n^{z2:.2f} (±{ze2:.2f}), about {T['TAU_NEW_COEF']} moves per cell.",
           f"Size exponent from n ≥ {NU_NMIN:,}: nu = {T['NU']} ± {T['NU_ERR']} (literature 0.6412 ± 0.0005).",
           f"Klarner's constant from n ≥ {LAM_NMIN:,}: lambda = {T['LAM']} ± {T['LAM_ERR']} (accepted 4.0625696)."]
    open(os.path.join(ROOT, "analysis", "results.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(f"{k} = {v}" for k, v in T.items() if k not in ("LADDER", "SQUARE_STRIP")))
    print("gallery / fit:", " ".join(f"{g[0]}:{g[2]:.3f}" for g in gl))
    print("fits", D["fit_old"], D["fit_new"], D["fit_nu"])


if __name__ == "__main__":
    main()
