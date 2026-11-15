// Exact, independent samples of site-percolation clusters of a given size (Leath growth).
//
//   leath -n N -p P -count C [-tsec S] [-seed K] -out FILE [-snap PREFIX -nsnap M]
//
// Grow the cluster of the origin: every newly reached neighbour site is occupied with
// probability P, otherwise blocked for good.  A finished cluster C of n cells with t blocked
// neighbours (its site perimeter) has probability P^(n-1) (1-P)^t.  Keeping only the clusters
// that finish with exactly N cells therefore gives independent draws of a fixed polyomino with
// probability proportional to (1-P)^perimeter: the ensemble `manyomino -p P` samples by Markov
// chain.  No relaxation time is involved, only a low acceptance rate.
// Output: one row of float64 per sample, [sample index, features...]; names in FILE.names.
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "features.h"
typedef uint64_t u64; typedef uint32_t u32;
struct Rng {
  u64 s[4];
  static u64 sm(u64 &x) { u64 z = (x += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }
  explicit Rng(u64 seed) { for (auto &v : s) v = sm(seed); }
  static u64 rotl(u64 x, int k) { return (x << k) | (x >> (64 - k)); }
  u64 next() { u64 r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17; s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45); return r; }
};
int main(int argc, char **argv) {
  int n = 1000, nsnap = 0; double p = 0.592746, tsec = 0; u64 count = 100, seed = 1; std::string out, snap;
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i]; auto nxt = [&]() { return std::string(argv[++i]); };
    if (a == "-n") n = atoi(nxt().c_str()); else if (a == "-p") p = atof(nxt().c_str());
    else if (a == "-count") count = (u64)atof(nxt().c_str()); else if (a == "-seed") seed = strtoull(nxt().c_str(), 0, 10);
    else if (a == "-tsec") tsec = atof(nxt().c_str()); else if (a == "-out") out = nxt();
    else if (a == "-snap") snap = nxt(); else if (a == "-nsnap") nsnap = atoi(nxt().c_str());
    else { fprintf(stderr, "unknown arg %s\n", a.c_str()); return 1; }
  }
  // sites live in an open-addressing hash table that is emptied by bumping a stamp
  int lg = 4; while ((1 << lg) < 16 * n) lg++;
  const u32 M = (1u << lg) - 1;
  std::vector<u64> key((size_t)M + 1); std::vector<u32> st((size_t)M + 1, 0);
  u32 stamp = 0;
  std::vector<int32_t> X(n + 1), Y(n + 1);
  const int DX[4] = {1, 0, -1, 0}, DY[4] = {0, 1, 0, -1};
  const u64 thr = (u64)(p * 18446744073709551615.0);
  Rng rng(seed); Features F; std::vector<double> f;
  FILE *fo = out.empty() ? nullptr : fopen(out.c_str(), "wb");
  if (fo) { FILE *fn = fopen((out + ".names").c_str(), "w"); fprintf(fn, "index"); for (auto &s : Features::names()) fprintf(fn, " %s", s.c_str()); fprintf(fn, "\n"); fclose(fn); }
  auto now = []() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
  const double t0 = now(); u64 tries = 0, sites = 0, got = 0;
  while (got < count) {
    if ((tries & 4095) == 0 && tsec > 0 && now() - t0 > tsec) break;
    tries++; stamp++;
    int m = 1, head = 0; X[0] = Y[0] = 0;
    { u64 k = 0; u32 h = (u32)((k * 0x9e3779b97f4a7c15ULL) >> (64 - lg)); key[h] = k; st[h] = stamp; }
    bool over = false;
    while (head < m && !over) {
      const int x = X[head], y = Y[head]; head++;
      for (int d = 0; d < 4; d++) {
        const int ax = x + DX[d], ay = y + DY[d];
        const u64 k = ((u64)(u32)ax << 32) | (u32)ay;
        u32 h = (u32)((k * 0x9e3779b97f4a7c15ULL) >> (64 - lg));
        bool seen = false;
        while (st[h] == stamp) { if (key[h] == k) { seen = true; break; } h = (h + 1) & M; }
        if (seen) continue;
        key[h] = k; st[h] = stamp; sites++;
        if (rng.next() < thr) { if (m == n) { over = true; break; } X[m] = ax; Y[m] = ay; m++; }
      }
    }
    if (over || m != n) continue;
    F.compute(X.data(), Y.data(), n, f);
    if (fo) { double idx = (double)got; fwrite(&idx, 8, 1, fo); fwrite(f.data(), 8, f.size(), fo); if ((got & 63) == 0) fflush(fo); }
    if (!snap.empty() && (int)got < nsnap) {
      char fn[512]; snprintf(fn, sizeof fn, "%s_%04d.bin", snap.c_str(), (int)got);
      FILE *fs = fopen(fn, "wb"); int32_t hdr = n; fwrite(&hdr, 4, 1, fs);
      for (int i = 0; i < n; i++) { int32_t xy[2] = {X[i], Y[i]}; fwrite(xy, 4, 2, fs); }
      fclose(fs);
    }
    got++;
  }
  if (fo) fclose(fo);
  const double el = now() - t0;
  fprintf(stderr, "n=%d p=%.6f samples=%" PRIu64 " tries=%" PRIu64 " (accept %.3g) sites/try=%.1f  %.1fs  %.3g s/sample  %.1f Msites/s\n",
          n, p, got, tries, (double)got / tries, (double)sites / tries, el, el / std::max<u64>(1, got), sites / el / 1e6);
  return 0;
}
