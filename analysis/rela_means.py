"""Do the two radius laws of the relocate move give the same averages?

Both are reversible with respect to the uniform distribution, so they must.  This sets the
averages of the 60 observables along chains run with -relA 1 (clump_a1/) against chains run
with -relA 2 (conv/), at 16,000 and 64,000 cells, with errors from the scatter of blocks.

  python3 analysis/rela_means.py
"""
import glob, math, os
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
NB = 12        # blocks per chain


def blocks(fn, burn):
    nm = open(fn + ".names").read().split(); a = np.memmap(fn, dtype=np.float64, mode="r")
    a = np.array(a[: len(a) // len(nm) * len(nm)].reshape(-1, len(nm))[:, 1:]); a = a[int(len(a) * burn):]; L = len(a) // NB
    return nm[1:], np.array([a[i * L:(i + 1) * L].mean(0) for i in range(NB)]), len(a)


def side(files, burn):
    B = []; n = 0
    for f in files: names, b, k = blocks(f, burn); B.append(b); n += k
    B = np.vstack(B); return names, B.mean(0), B.std(0, ddof=1) / math.sqrt(len(B)), n, len(files)


for n, pat in ((16000, ["clump_a1/ts_?_n16000.part0.feat"]), (64000, ["clump_a1/ts_?_n64000.part0.feat", "clump_a1/more/m?.feat"])):
    f1 = sorted(f for p in pat for f in glob.glob(os.path.join(ROOT, p))); f2 = sorted(glob.glob(os.path.join(ROOT, "conv", f"u{n}_?.feat")))
    if not f1 or not f2: continue
    names, m1, e1, n1, c1 = side(f1, 0.02); _, m2, e2, n2, c2 = side(f2, 0.1)
    z = (m1 - m2) / np.sqrt(e1 ** 2 + e2 ** 2); o = np.argsort(-np.abs(z))
    print(f"{n:,} cells: -relA 1, {n1:,} readings in {c1} chains; -relA 2, {n2:,} in {c2}. Beyond 2 standard errors: {int((np.abs(z) > 2).sum())} of {len(z)}; "
          f"largest {np.abs(z).max():.1f}; rms {math.sqrt((z ** 2).mean()):.2f}")
    for j in o[:4]: print(f"    {names[j]:>10s} {z[j]:+.1f}   {m1[j]:.6g} ± {e1[j]:.2g}   against {m2[j]:.6g} ± {e2[j]:.2g}")
