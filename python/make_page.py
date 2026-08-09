"""Generate the static gallery page (adamscherlis.github.io/static/manyomino)
from the current production runs.

For each n it requires BOTH seeds' chains (bar + rect), measures tau and the
two-seed z on the fly, and only publishes sizes where:
  - the two-seed test passes (|z| < 3, ESS >= 10 per chain), and
  - every published snapshot is at least ~5 tau into its chain.
Snapshots from one chain are spaced >= 2 tau apart.  Everything else is
listed as "still equilibrating".

Usage: python3 make_page.py <site_dir> [--max-per-n 4]
where <site_dir> is the checkout of adamscherlis.github.io.
"""

import glob
import math
import os
import re
import subprocess
import sys
import time

from analyze import analyze
from render import load, render

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROD = os.path.join(ROOT, "results", "prod")
RAW = os.path.join(ROOT, "gallery", "raw")


def chain_stats(n):
    out = {}
    for init in ("bar", "rect"):
        path = os.path.join(PROD, f"n{n}_{init}.csv")
        if not os.path.exists(path) or os.path.getsize(path) < 10000:
            return None
        try:
            out[init] = analyze(path, 0.5)
        except Exception:
            return None
    a, b = out["bar"], out["rect"]
    denom = math.sqrt(a["rg2_stderr"] ** 2 + b["rg2_stderr"] ** 2)
    z = abs(a["rg2_mean"] - b["rg2_mean"]) / denom if denom > 0 else float("inf")
    tau = max(a["tau_steps"], b["tau_steps"])
    wts = [1 / a["rg_stderr"] ** 2, 1 / b["rg_stderr"] ** 2]
    rg = (a["rg_mean"] * wts[0] + b["rg_mean"] * wts[1]) / sum(wts)
    return {
        "bar": a,
        "rect": b,
        "z": z,
        "tau": tau,
        "rg": rg,
        "converged": z < 3.0 and a["ess"] >= 8 and b["ess"] >= 8,
    }


def snapshots(n, tau, max_per_n=4):
    """Pick up to max_per_n snapshot files, >= 5 tau burn-in, >= 2 tau apart,
    alternating chains, newest first."""
    files = []
    for init in ("bar", "rect"):
        cand = []
        for f in glob.glob(os.path.join(RAW, f"n{n}_{init}_*.txt")):
            step = int(re.search(r"_(\d+)\.txt$", f).group(1))
            if step >= 5 * tau:
                cand.append((step, init, f))
        cand.sort(reverse=True)
        kept = []
        for step, init_, f in cand:
            if not kept or kept[-1][0] - step >= 2 * tau:
                kept.append((step, init_, f))
        files.append(kept)
    picked = []
    i = 0
    while len(picked) < max_per_n and (files[0] or files[1]):
        lane = files[i % 2] or files[(i + 1) % 2]
        if lane:
            picked.append(lane.pop(0))
        i += 1
    return picked


CSS = """
:root{
  --ground:#F6F7FA; --panel:#FFFFFF; --ink:#15171D; --mut:#666C7E;
  --rule:#DDE1EA; --acc:#3B5BF0; --accw:#EEF1FE; --ok:#2e7d4f; --pend:#96690f;
}
@media (prefers-color-scheme:dark){
  :root{--ground:#0E1014; --panel:#15181F; --ink:#E9EBF2; --mut:#9AA0B2;
        --rule:#252932; --acc:#7D93FF; --accw:#161B2E; --ok:#6fbf8f; --pend:#d3a34a;}
}
:root[data-theme=dark]{--ground:#0E1014; --panel:#15181F; --ink:#E9EBF2; --mut:#9AA0B2;
        --rule:#252932; --acc:#7D93FF; --accw:#161B2E; --ok:#6fbf8f; --pend:#d3a34a;}
:root[data-theme=light]{--ground:#F6F7FA; --panel:#FFFFFF; --ink:#15171D; --mut:#666C7E;
        --rule:#DDE1EA; --acc:#3B5BF0; --accw:#EEF1FE; --ok:#2e7d4f; --pend:#96690f;}
*{box-sizing:border-box}
body{margin:0;background:var(--ground);color:var(--ink);
  font:17px/1.6 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;}
.wrap{max-width:1060px;margin:0 auto;padding:0 1.25rem 4rem}
h1{font-size:1.7rem;margin:2.2rem 0 .3rem;letter-spacing:-.01em}
.lede{color:var(--mut);max-width:46rem}
code{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:.9em;
  background:color-mix(in srgb,var(--ink) 7%,transparent);padding:.1em .35em;border-radius:4px}
a{color:var(--acc)}
.card{background:var(--panel);border:1px solid var(--rule);border-radius:12px;
  padding:1.1rem 1.25rem;margin:1.4rem 0}
.card h2{margin:.1rem 0 .4rem;font-size:1.25rem}
.sub{color:var(--mut);font-size:.92em;margin:.1rem 0 .6rem}
.badge{display:inline-block;font-size:.72em;letter-spacing:.06em;text-transform:uppercase;
  padding:.15em .55em;border-radius:99px;border:1px solid currentColor;vertical-align:2px;margin-left:.5em}
.ok{color:var(--ok)} .pend{color:var(--pend)}
.stats{display:flex;flex-wrap:wrap;gap:.4rem 1.6rem;margin:.3rem 0 .8rem;
  font-size:.9em;color:var(--mut)}
.stats b{color:var(--ink);font-weight:600}
.panels{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:.9rem}
figure{margin:0;background:#fff;border:1px solid var(--rule);border-radius:8px;padding:.5rem}
figure img{width:100%;height:auto;display:block;image-rendering:pixelated}
figcaption{font-size:.78em;color:var(--mut);margin-top:.35rem}
footer{color:var(--mut);font-size:.85em;border-top:1px solid var(--rule);
  margin-top:2.5rem;padding-top:1rem}
"""


def build(site_dir, max_per_n=4):
    dest = os.path.join(site_dir, "static", "manyomino")
    figdir = os.path.join(dest, "figs")
    os.makedirs(figdir, exist_ok=True)

    sizes = sorted(
        {
            int(m.group(1))
            for f in glob.glob(os.path.join(PROD, "n*_bar.csv"))
            if (m := re.search(r"n(\d+)_bar\.csv$", f))
        },
        reverse=True,
    )

    cards, pending = [], []
    for n in sizes:
        st = chain_stats(n)
        if st is None:
            pending.append((n, "waiting for both chains"))
            continue
        if not st["converged"]:
            pending.append(
                (n, f"two-seed |z| = {st['z']:.1f}, "
                    f"ESS = {st['bar']['ess']:.0f}/{st['rect']['ess']:.0f}")
            )
            continue
        snaps = snapshots(n, st["tau"], max_per_n)
        if not snaps:
            pending.append((n, "converged but no snapshot past 5 tau yet"))
            continue
        figs = []
        for step, init, f in snaps:
            name = f"n{n}_{init}_{step}"
            png = os.path.join(figdir, name + ".png")
            if not os.path.exists(png):
                cells = load(f)
                render(cells, png, mode="dist")
            rg = math.sqrt(analyze_rg2_of_dump(f))
            figs.append(
                f'<figure><img src="figs/{name}.png" alt="random polyomino, '
                f'n={n}"><figcaption>{init}-seeded chain, step {step:.2e}'
                f' &middot; R<sub>g</sub> = {rg:.0f}</figcaption></figure>'
            )
        a, b = st["bar"], st["rect"]
        cards.append(f"""
<section class="card">
<h2>n = {n:,}<span class="badge ok">two-seed converged</span></h2>
<p class="sub">independent chains from a stringy bar seed and a compact
&radic;n&times;&radic;n rectangle seed agree: |z| = {st['z']:.2f} on
&lang;R<sub>g</sub>&sup2;&rang;</p>
<div class="stats">
  <div>&lang;R<sub>g</sub>&rang; <b>{st['rg']:.1f}</b></div>
  <div>&tau; <b>{st['tau']:.1e}</b> moves</div>
  <div>ESS <b>{a['ess']:.0f}</b> (bar) / <b>{b['ess']:.0f}</b> (rect)</div>
  <div>chain length <b>{a['records'] * a['interval'] * 2:.1e}</b> moves each</div>
</div>
<div class="panels">{''.join(figs)}</div>
</section>""")

    pend_html = ""
    if pending:
        items = "".join(
            f"<li><b>n = {n:,}</b> &mdash; {why}</li>" for n, why in pending
        )
        pend_html = (
            '<section class="card"><h2>Still equilibrating'
            '<span class="badge pend">running</span></h2>'
            f"<ul>{items}</ul></section>"
        )

    stamp = time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Uniform random polyominoes</title>
<style>{CSS}</style>
</head>
<body>
<div class="wrap">
<h1>Uniform random polyominoes</h1>
<p class="lede">Each image is a single polyomino drawn <b>uniformly at random</b> among all
fixed polyominoes of its size, sampled by a provably-uniform Markov chain (relocate a random
non-cut cell to a uniform perimeter site; the kernel is symmetric, so the uniform distribution
is stationary).  Uniform polyominoes sit in the branched-polymer universality class: sparse
dendrites with R<sub>g</sub> ~ n<sup>&nu;</sup>, &nu; &asymp; 0.64 &mdash; nothing like compact
blobs.  Cells are colored by graph distance from an extremal cell to show the branch
structure.</p>
<p class="lede" style="font-size:.95em">Sampler, validation suite (exact-count cross-checks,
chi-square uniformity against OEIS <a href="https://oeis.org/A001168">A001168</a>, two-seed
convergence, autocorrelation) and reproduction instructions:
<a href="https://github.com/AdamScherlis/manyomino">AdamScherlis/manyomino</a>.
A size appears here only after two independently-seeded chains agree on
&lang;R<sub>g</sub>&sup2;&rang; and every shown snapshot is &ge;5&tau; into its chain.</p>
{''.join(cards)}
{pend_html}
<footer>Generated {stamp} from the manyomino production runs.</footer>
</div>
</body>
</html>
"""
    with open(os.path.join(dest, "index.html"), "w") as f:
        f.write(html)
    print(f"wrote {dest}/index.html  ({len(cards)} converged sizes, {len(pending)} pending)")
    return len(cards)


def analyze_rg2_of_dump(path):
    cells = load(path)
    n = len(cells)
    sx = sum(x for x, y in cells)
    sy = sum(y for x, y in cells)
    sx2 = sum(x * x for x, y in cells)
    sy2 = sum(y * y for x, y in cells)
    return (sx2 + sy2) / n - (sx / n) ** 2 - (sy / n) ** 2


if __name__ == "__main__":
    site = sys.argv[1]
    mx = 4
    if "--max-per-n" in sys.argv:
        mx = int(sys.argv[sys.argv.index("--max-per-n") + 1])
    build(site, mx)
