"""How long does a piece of the starting shape survive?

Between a starting shape and the shape t moves later, cells keep their identity, so each cell has
a displacement.  A cell that has never been moved has displacement zero; a branch that was
transplanted without being turned keeps one common displacement and its internal arrangement.
The size of the largest group of cells sharing one displacement is therefore the largest piece
of the start that is still intact: a property that no average over the shape would show, and a
direct look at the last thing the walk forgets.

  python3 analysis/fossil.py        # reads surv/ (see scripts/fossil.sh), writes analysis/fossil.json
"""
import glob, json, os, re
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")


def load(f):
    a = np.fromfile(f, dtype=np.int32); n = a[0]
    return a[1:1 + 2 * n].reshape(n, 2)


def largest_group(a, b):
    d = (b - a).astype(np.int64); key = d[:, 0] * (1 << 32) + d[:, 1]
    _, cnt = np.unique(key, return_counts=True)
    return int(cnt.max()), int((key == 0).sum())


def main():
    out = []
    for start in sorted(glob.glob(os.path.join(ROOT, "surv", "s*_start.bin")), key=lambda f: int(re.search(r"s(\d+)([a-z]?)_start", f).group(1))):
        tag = os.path.basename(start)[:-10]; n = int(re.match(r"s(\d+)", tag).group(1))
        step = float(open(os.path.join(ROOT, "surv", tag + ".step")).read())          # moves per cell between snapshots
        a = load(start); rows = []
        for k, f in enumerate(sorted(glob.glob(os.path.join(ROOT, "surv", tag + "_????.bin")))):
            g, z = largest_group(a, load(f)); rows.append([step * (k + 1), g, z])
        # first time the largest surviving piece is below 0.1% of the shape, and at the chance level (a handful of cells)
        t01 = next((r[0] for r in rows if r[1] < 1e-3 * n), None); tend = next((r[0] for r in rows if r[1] <= 12), None)
        tzero = next((r[0] for r in rows if r[2] == 0), None)
        out.append({"tag": tag, "n": n, "rows": rows, "t_below_0.1pct": t01, "t_gone": tend, "t_all_moved": tzero})
        print(f"{tag}: n={n}, {len(rows)} snapshots; all cells moved at least once by {tzero} moves per cell; largest intact piece below 0.1% by {t01}, down to a dozen cells by {tend}")
    json.dump(out, open(os.path.join(ROOT, "analysis", "fossil.json"), "w"))
    return out


if __name__ == "__main__":
    main()
