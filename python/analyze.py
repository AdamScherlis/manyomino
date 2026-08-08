"""Analyze an rg2 time series CSV (step,rg2) from the Rust sampler.

Reports mean R_g^2, mean R_g, integrated autocorrelation time tau (in chain
steps), ESS, and a batch-means stderr cross-check.

Usage: python3 analyze.py ts.csv [burn_frac]
"""

import math
import sys

from stats import autocorr_time, mean_stderr_ess


def load(path):
    steps, vals = [], []
    with open(path) as f:
        next(f)
        for line in f:
            s, v = line.split(",")
            steps.append(int(s))
            vals.append(float(v))
    return steps, vals


def binned(xs, b):
    return [sum(xs[i : i + b]) / b for i in range(0, len(xs) - b + 1, b)]


def tau_in_records(xs):
    """Integrated autocorrelation time in units of the recording interval,
    computed on a binned series (bin chosen so the windowed sum is cheap),
    then unbinned.  Returns tau in records."""
    n = len(xs)
    b = max(1, n // 20000)
    ys = binned(xs, b)
    tau_b, w = autocorr_time(ys)
    tau = tau_b * b
    # if tau is not >> bin, redo with smaller bin for accuracy
    if b > 1 and tau < 5 * b:
        b2 = max(1, int(tau) // 20 + 1)
        if b2 < b:
            ys = binned(xs, b2)
            tau_b, w = autocorr_time(ys)
            tau = tau_b * b2
    return tau


def batch_stderr(xs, nbatch=32):
    n = len(xs)
    b = n // nbatch
    if b < 1:
        return float("nan")
    means = [sum(xs[i * b : (i + 1) * b]) / b for i in range(nbatch)]
    m = sum(means) / nbatch
    var = sum((x - m) ** 2 for x in means) / (nbatch - 1)
    return math.sqrt(var / nbatch)


def analyze(path, burn_frac=0.3):
    steps, vals = load(path)
    interval = steps[1] - steps[0] if len(steps) > 1 else 1
    cut = int(len(vals) * burn_frac)
    xs = vals[cut:]
    tau_rec = tau_in_records(xs)
    tau_steps = tau_rec * interval
    mean2, se2, ess = mean_stderr_ess(xs, tau_rec)
    rg = [math.sqrt(max(v, 0.0)) for v in xs]
    mean1, se1, _ = mean_stderr_ess(rg, tau_rec)
    bse = batch_stderr(xs)
    return {
        "file": path,
        "records": len(xs),
        "interval": interval,
        "tau_steps": tau_steps,
        "ess": ess,
        "rg2_mean": mean2,
        "rg2_stderr": se2,
        "rg2_batch_stderr": bse,
        "rg_mean": mean1,
        "rg_stderr": se1,
    }


if __name__ == "__main__":
    burn = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3
    r = analyze(sys.argv[1], burn)
    print(
        f"{r['file']}: records={r['records']} interval={r['interval']}\n"
        f"  tau = {r['tau_steps']:.3e} steps   ESS = {r['ess']:.0f}\n"
        f"  <Rg^2> = {r['rg2_mean']:.4f} +- {r['rg2_stderr']:.4f} "
        f"(batch stderr {r['rg2_batch_stderr']:.4f})\n"
        f"  <Rg>  = {r['rg_mean']:.4f} +- {r['rg_stderr']:.4f}"
    )
