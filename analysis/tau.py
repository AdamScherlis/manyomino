"""Integrated autocorrelation time (Sokal automatic windowing) of manyomino time series."""
import sys, numpy as np

def acf_fft(x):
    x = np.asarray(x, float) - np.mean(x)
    n = len(x); m = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(x, m)
    a = np.fft.irfft(f * np.conj(f), m)[:n]
    a /= np.arange(n, 0, -1)
    return a / a[0]

def tau_int(x, c=6.0):
    """returns (tau_int in samples, window, err). tau_int = 1/2 + sum_{t>=1} rho(t)."""
    rho = acf_fft(x)
    cs = 0.5 + np.cumsum(rho[1:])
    w = np.arange(1, len(cs) + 1)
    ok = np.nonzero(w >= c * cs)[0]
    i = ok[0] if len(ok) else len(cs) - 1
    tau = cs[i]
    return tau, w[i], tau * np.sqrt(2.0 * (2 * w[i] + 1) / len(x))

def load(fn, burn_frac=0.1):
    d = np.loadtxt(fn)
    hdr = open(fn).readline()
    k = int(len(d) * burn_frac)
    return d[k:], hdr

if __name__ == "__main__":
    for fn in sys.argv[1:]:
        d, hdr = load(fn)
        every = d[1, 0] - d[0, 0]; sec = (d[-1, 1] - d[0, 1]) / (len(d) - 1)
        out = [fn.split('/')[-1], f"N={len(d)}"]
        for name, col in (("rg2", 2), ("asph", 3), ("perim", 4), ("noncut", 6)):
            t, w, e = tau_int(d[:, col])
            m = d[:, col].mean(); se = d[:, col].std() * np.sqrt(2 * t / len(d))
            out.append(f"{name}: mean={m:.6g}±{se:.2g} tau={t*every:.3g}±{e*every:.2g} moves ({t*sec:.3g}s)")
        print("  ".join(out))
