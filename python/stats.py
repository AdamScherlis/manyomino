"""Small stdlib-only statistics helpers: chi-square p-values, autocorrelation."""

import math


def _gser(a, x, itmax=500, eps=3e-12):
    """Lower regularized incomplete gamma P(a,x) by series (x < a+1)."""
    ap = a
    s = 1.0 / a
    delta = s
    for _ in range(itmax):
        ap += 1.0
        delta *= x / ap
        s += delta
        if abs(delta) < abs(s) * eps:
            break
    return s * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gcf(a, x, itmax=500, eps=3e-12):
    """Upper regularized incomplete gamma Q(a,x) by continued fraction (x >= a+1)."""
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, itmax + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h * math.exp(-x + a * math.log(x) - math.lgamma(a))


def gammaincc(a, x):
    """Upper regularized incomplete gamma Q(a, x)."""
    if x < 0 or a <= 0:
        raise ValueError
    if x == 0:
        return 1.0
    if x < a + 1.0:
        return 1.0 - _gser(a, x)
    return _gcf(a, x)


def chi2_sf(chi2, dof):
    """Survival function (p-value) of the chi-square distribution."""
    return gammaincc(dof / 2.0, chi2 / 2.0)


def chisquare_uniform(counts, total):
    """Chi-square statistic and p-value of `counts` against the uniform
    distribution over len(counts) categories with `total` observations.
    `counts` must include zero entries for unobserved categories."""
    k = len(counts)
    e = total / k
    chi2 = sum((o - e) ** 2 / e for o in counts)
    dof = k - 1
    return chi2, dof, chi2_sf(chi2, dof)


def autocorr_time(xs, c=5.0):
    """Integrated autocorrelation time, tau = 1/2 + sum_{t=1..W} rho(t), with
    Sokal's automatic window: smallest W with W >= c*tau(W).  ESS = N/(2 tau).
    Returns (tau, W). Units: the series spacing."""
    n = len(xs)
    mean = sum(xs) / n
    var = sum((x - mean) ** 2 for x in xs) / n
    if var == 0:
        return 0.5, 0
    tau = 0.5
    for w in range(1, n // 2):
        s = 0.0
        for i in range(n - w):
            s += (xs[i] - mean) * (xs[i + w] - mean)
        rho = s / ((n - w) * var)
        tau += rho
        if w >= c * tau:
            return max(tau, 0.5), w
    return max(tau, 0.5), n // 2


def mean_stderr_ess(xs, tau):
    """Mean and autocorrelation-corrected standard error given tau (in units of
    the series spacing); ESS = N / (2 tau) truncated at N."""
    n = len(xs)
    mean = sum(xs) / n
    var = sum((x - mean) ** 2 for x in xs) / (n - 1) if n > 1 else 0.0
    ess = min(n, n / (2.0 * tau)) if tau > 0 else n
    stderr = math.sqrt(var / ess) if ess > 0 else float("inf")
    return mean, stderr, ess
