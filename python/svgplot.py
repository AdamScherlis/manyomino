"""Hand-rolled log-log SVG plot of <Rg> vs n with the nu fit line.
Writes a theme-friendly standalone SVG (currentColor text, fixed accents)."""
import math
import sys

from analyze import analyze


def collect():
    groups = [
        (100, ["../results/tau_n100.csv"]),
        (200, ["../results/prod/n200_rect.csv"]),
        (300, ["../results/tau_n300.csv"]),
        (500, ["../results/prod/n500_rect.csv"]),
        (1000, ["../results/prod/n1000_bar.csv", "../results/prod/n1000_rect.csv"]),
        (2000, ["../results/prod/n2000_bar.csv", "../results/prod/n2000_rect.csv"]),
        (3000, ["../results/prod/n3000_bar.csv", "../results/prod/n3000_rect.csv"]),
        (5000, ["../results/prod/n5000_bar.csv", "../results/prod/n5000_rect.csv"]),
        (10000, ["../results/prod/n10000_bar.csv", "../results/prod/n10000_rect.csv"]),
        (30000, ["../results/prod/n30000_bar.csv", "../results/prod/n30000_rect.csv"]),
        (100000, sorted(__import__("glob").glob("../results/prod/n100000_*.csv"))),
    ]
    pts = []
    for n, paths in groups:
        ms, ses = [], []
        for p in paths:
            try:
                r = analyze(p, 0.5)
                ms.append(r["rg_mean"])
                ses.append(r["rg_stderr"])
            except Exception:
                pass
        w = [1 / s**2 for s in ses]
        m = sum(mi * wi for mi, wi in zip(ms, w)) / sum(w)
        se = math.sqrt(1 / sum(w))
        pts.append((n, m, se))
    return pts


def wls_fit(pts):
    xs = [math.log(n) for n, m, s in pts]
    ys = [math.log(m) for n, m, s in pts]
    ws = [(m / s) ** 2 for n, m, s in pts]
    sw = sum(ws)
    sx = sum(w * x for w, x in zip(ws, xs))
    sy = sum(w * y for w, y in zip(ws, ys))
    sxx = sum(w * x * x for w, x in zip(ws, xs))
    sxy = sum(w * x * y for w, x, y in zip(ws, xs, ys))
    d = sw * sxx - sx * sx
    b = (sw * sxy - sx * sy) / d
    a = (sxx * sy - sx * sxy) / d
    sb = math.sqrt(sw / d)
    return a, b, sb


def make_svg(pts, out):
    a, nu, snu = wls_fit(pts)
    W, H = 640, 420
    L, R, T, B = 64, 20, 24, 48
    x0, x1 = math.log10(70), math.log10(150000)
    y0, y1 = math.log10(6), math.log10(500)
    def X(n):
        return L + (math.log10(n) - x0) / (x1 - x0) * (W - L - R)
    def Y(v):
        return H - B - (math.log10(v) - y0) / (y1 - y0) * (H - T - B)
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
         f'font-family="system-ui,sans-serif" font-size="13">']
    s.append(f'<rect width="{W}" height="{H}" fill="none"/>')
    # decade grid
    for e in range(2, 6):
        n = 10 ** e
        if x0 <= math.log10(n) <= x1:
            s.append(f'<line x1="{X(n):.1f}" y1="{T}" x2="{X(n):.1f}" y2="{H-B}" '
                     f'stroke="#8888" stroke-width="0.5"/>')
            s.append(f'<text x="{X(n):.1f}" y="{H-B+18}" text-anchor="middle" '
                     f'fill="currentColor">10<tspan dy="-5" font-size="9">{e}</tspan></text>')
    for e in range(1, 3):
        v = 10 ** e
        s.append(f'<line x1="{L}" y1="{Y(v):.1f}" x2="{W-R}" y2="{Y(v):.1f}" '
                 f'stroke="#8888" stroke-width="0.5"/>')
        s.append(f'<text x="{L-8}" y="{Y(v)+4:.1f}" text-anchor="end" '
                 f'fill="currentColor">10<tspan dy="-5" font-size="9">{e}</tspan></text>')
    # fit line
    nlo, nhi = 80, 52000
    s.append(f'<line x1="{X(nlo):.1f}" y1="{Y(math.exp(a)*nlo**nu):.1f}" '
             f'x2="{X(nhi):.1f}" y2="{Y(math.exp(a)*nhi**nu):.1f}" '
             f'stroke="#D2500F" stroke-width="1.8"/>')
    # points
    for n, m, se in pts:
        s.append(f'<line x1="{X(n):.1f}" y1="{Y(m-2*se):.1f}" x2="{X(n):.1f}" '
                 f'y2="{Y(m+2*se):.1f}" stroke="#3B5BF0" stroke-width="1.4"/>')
        s.append(f'<circle cx="{X(n):.1f}" cy="{Y(m):.1f}" r="4" fill="#3B5BF0"/>')
    s.append(f'<text x="{L+10}" y="{T+18}" fill="currentColor" font-size="15">'
             f'&#x27E8;R<tspan dy="3" font-size="10">g</tspan><tspan dy="-3">&#x27E9; '
             f'~ n<tspan dy="-6" font-size="11">&#957;</tspan><tspan dy="6">,  '
             f'&#957; = {nu:.4f} &#177; {snu:.4f}</tspan></text>')
    s.append(f'<text x="{(L+W-R)/2}" y="{H-10}" text-anchor="middle" '
             f'fill="currentColor">n (cells)</text>')
    s.append('</svg>')
    open(out, "w").write("\n".join(s))
    print(f"{out}: nu = {nu:.4f} +- {snu:.4f}")


if __name__ == "__main__":
    pts = collect()
    make_svg(pts, sys.argv[1] if len(sys.argv) > 1 else "../results/nufit.svg")
