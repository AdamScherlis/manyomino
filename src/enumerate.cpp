// Exact averages of every observable in features.h over ALL fixed polyominoes of N cells
// (Redelmeier's algorithm).  Ground truth for the samplers at sizes where listing is possible.
//   enumerate N [P] > exact_N.txt    lines: name mean sd      (first line: count)
// With P, each shape is weighted by (1-P)^perimeter: the law of a size-N site-percolation cluster.
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include "features.h"
static int N; static const int W = 64;
static uint8_t reached[W * W]; static int32_t X[64], Y[64];
static std::vector<long double> S1, S2; static uint64_t total = 0;
static long double PW[256], SW = 0; static int iperim = 0;
static Features F; static std::vector<double> f;
static int untried[4096];
static const int DX[4] = {1, 0, -1, 0}, DY[4] = {0, 1, 0, -1};
static void rec(int ubegin, int uend, int p) {
  // cells untried[ubegin..uend) may still be added; p cells placed
  while (uend > ubegin) {
    const int c = untried[--uend]; const int x = c % W, y = c / W;
    X[p] = x; Y[p] = y;
    if (p + 1 == N) {
      total++; F.compute(X, Y, N, f);
      const long double w = PW[(int)lround(f[iperim] * N)]; SW += w;
      for (size_t i = 0; i < f.size(); i++) { S1[i] += w * f[i]; S2[i] += w * (long double)f[i] * f[i]; }
      continue;
    }
    int added = 0;
    for (int d = 0; d < 4; d++) {
      const int ax = x + DX[d], ay = y + DY[d], a = ay * W + ax;
      if (ay < 1 || (ay == 1 && ax < W / 2)) continue;          // keep the origin the lowest, then leftmost, cell
      if (reached[a]) continue;
      reached[a] = 1; untried[uend + added++] = a;
    }
    // recurse on a private copy of the list above the current top
    const int nb = uend, ne = uend + added;
    for (int i = ubegin; i < ne; i++) untried[ne + (i - ubegin)] = untried[i];
    rec(ne, ne + (ne - ubegin), p + 1);
    for (int i = 0; i < added; i++) reached[untried[nb + i]] = 0;
  }
}
int main(int argc, char **argv) {
  N = atoi(argv[1]);
  const double P = argc > 2 ? atof(argv[2]) : 0;
  for (int t = 0; t < 256; t++) PW[t] = powl(1 - P, t);
  while (Features::names()[iperim] != "perim") iperim++;
  S1.assign(Features::names().size(), 0); S2.assign(Features::names().size(), 0);
  const int o = 1 * W + W / 2; reached[o] = 1; untried[0] = o;
  rec(0, 1, 0);
  printf("count %" PRIu64 "\n", total);
  for (size_t i = 0; i < S1.size(); i++) {
    long double m = S1[i] / SW, v = S2[i] / SW - m * m;
    printf("%s %.12Lf %.12Lf\n", Features::names()[i].c_str(), m, sqrtl(v > 0 ? v : 0));
  }
  return 0;
}
