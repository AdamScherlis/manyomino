"""Fit a_k ~ C lambda^k / k^theta from the PERM log_ak table.

Klarner's constant lambda ~ 4.0626; theory: theta = 1 for 2D animals.
Fits over k in [kmin, kmax] by least squares on ln a_k."""
import math
import sys

path = sys.argv[1]
kmin = int(sys.argv[2]) if len(sys.argv) > 2 else 30
kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 180

ks, las = [], []
with open(path) as f:
    next(f)
    for line in f:
        k, la = line.split()
        k = int(k)
        if kmin <= k <= kmax:
            ks.append(k)
            las.append(float(la))

# design: ln a_k = c0 + k*ln(lambda) - theta*ln k  -> linear in (1, k, ln k)
import itertools
N = len(ks)
X = [[1.0, k, -math.log(k)] for k in ks]
# normal equations 3x3
A = [[sum(X[i][a] * X[i][b] for i in range(N)) for b in range(3)] for a in range(3)]
b = [sum(X[i][a] * las[i] for i in range(N)) for a in range(3)]
# solve 3x3
import copy
M = [row[:] + [b[i]] for i, row in enumerate(A)]
for col in range(3):
    piv = max(range(col, 3), key=lambda r: abs(M[r][col]))
    M[col], M[piv] = M[piv], M[col]
    for r in range(3):
        if r != col:
            f = M[r][col] / M[col][col]
            for c in range(col, 4):
                M[r][c] -= f * M[col][c]
sol = [M[i][3] / M[i][i] for i in range(3)]
c0, lnlam, theta = sol
resid = [las[i] - (c0 + ks[i] * lnlam - theta * math.log(ks[i])) for i in range(N)]
rms = math.sqrt(sum(r * r for r in resid) / N)
print(f"fit over k = {kmin}..{kmax} ({N} points)")
print(f"lambda = {math.exp(lnlam):.4f}   (Klarner: 4.0626)")
print(f"theta  = {theta:.3f}        (theory: 1)")
print(f"C      = {math.exp(c0):.4f}   rms resid = {rms:.4f}")
