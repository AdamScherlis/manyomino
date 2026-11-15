// Uniform polyominoes without a Markov chain: a population of growing clusters with resampling
// ("go with the winners"; a fixed-population relative of the PERM cluster growth that Hsu, Nadler
// and Grassberger used for lattice animals, J. Phys. A 38, 2005).
//
//   grow -nmax N -pop M [-step K] [-sizes 250,1000] -out FILE [-seed S]
//
// A walker is a cluster in growth: its cells in the order they were added, plus the empty sites it
// has ruled out ("blocked").  To add a cell, list the cluster's untried neighbour sites in a fixed
// order (by the cell that touches them, oldest cell first): s_1 .. s_G.  Pick j uniformly, block
// s_1 .. s_(j-1) for good, occupy s_j, and multiply the walker's weight by G.
// Why the weights are right: fix a polyomino C and one of its cells as the origin.  There is
// exactly one history that builds C: at every stage the next cell must be the first listed site
// that belongs to C (sites listed before it are not in C, so blocking them is harmless, and any
// other choice blocks a site of C).  That history is followed with probability 1 / (G_1 G_2 ...),
// so its weight is exactly the reciprocal: every polyomino of m cells is counted m times, the
// average weight at size m estimates m * A(m), and weighted averages are averages over uniformly
// random polyominoes.  When the weights spread out the population is resampled in proportion to
// weight (systematic resampling), which does not bias the counts.
// There is no relaxation time.  The weakness is different: copies share their early history, so
// the walkers are not independent.  The program reports how many distinct ancestors the final
// population has at earlier sizes, and independent runs (seeds) give honest error bars.
// Output: rows of float64 [m, log weight, ancestor at m/10, ancestor at m/2, features...] at the
// requested sizes, names in FILE.names, and FILE.logz with log(m A(m)) estimates.
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <unordered_set>
#include "features.h"
typedef uint64_t u64; typedef uint32_t u32;
struct Rng {
  u64 s[4];
  static u64 sm(u64 &x) { u64 z = (x += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }
  explicit Rng(u64 seed) { for (auto &v : s) v = sm(seed); }
  static u64 rotl(u64 x, int k) { return (x << k) | (x >> (64 - k)); }
  u64 next() { u64 r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17; s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45); return r; }
  double uni() { return (next() >> 11) * (1.0 / 9007199254740992.0); }
};
struct Walker {
  std::vector<u32> cells, blocked;      // grid indices, cells in breadth-first order
  int head = 0;                         // cells before this one have no untried neighbours
  int G = 4;                            // untried neighbour sites of the cluster
  double lw = 0;                        // log weight since the last resampling; -inf if closed off
  std::vector<u32> anc;                 // ancestor label at each checkpoint size
};

int main(int argc, char **argv) {
  int nmax = 1000, pop = 10000, K = 4; double essFrac = 0.6, kappa = 1.0; u64 seed = 1; std::string out, sizesArg = "";
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i]; auto nxt = [&]() { return std::string(argv[++i]); };
    if (a == "-nmax") nmax = atoi(nxt().c_str()); else if (a == "-ess") essFrac = atof(nxt().c_str()); else if (a == "-kappa") kappa = atof(nxt().c_str());
    else if (a == "-pop") pop = (int)atof(nxt().c_str()); else if (a == "-seed") seed = strtoull(nxt().c_str(), 0, 10);
    else if (a == "-step") K = atoi(nxt().c_str()); else if (a == "-out") out = nxt();
    else if (a == "-sizes") sizesArg = nxt();
    else { fprintf(stderr, "unknown arg %s\n", a.c_str()); return 1; }
  }
  std::vector<int> sizes;
  { std::stringstream ss(sizesArg); std::string t; while (std::getline(ss, t, ',')) if (!t.empty()) sizes.push_back(atoi(t.c_str())); }
  // checkpoints at which every walker's ancestor label is reset to itself
  std::vector<int> cps; for (int s : sizes) { cps.push_back(std::max(K, s / 10 / K * K)); cps.push_back(std::max(K, s / 2 / K * K)); }
  std::sort(cps.begin(), cps.end()); cps.erase(std::unique(cps.begin(), cps.end()), cps.end());
  auto cpIndex = [&](int m) { for (size_t i = 0; i < cps.size(); i++) if (cps[i] == m) return (int)i; return -1; };
  const int L = 2 * nmax + 3, C = nmax + 1;
  std::vector<uint8_t> grid((size_t)L * L, 0);
  const long nb[4] = {1, L, -1, -(long)L};
  Rng rng(seed); Features F; std::vector<double> f; std::vector<int32_t> X(nmax), Y(nmax);
  std::vector<Walker> W(pop);
  for (auto &w : W) { w.cells = {(u32)((size_t)C * L + C)}; w.anc.assign(cps.size(), 0); }
  FILE *fo = out.empty() ? nullptr : fopen(out.c_str(), "wb"), *fz = out.empty() ? nullptr : fopen((out + ".logz").c_str(), "w");
  if (fo) { FILE *fn = fopen((out + ".names").c_str(), "w"); fprintf(fn, "m logw anc10 anc2"); for (auto &s : Features::names()) fprintf(fn, " %s", s.c_str()); fprintf(fn, "\n"); fclose(fn); }
  auto now = []() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
  const double t0 = now(); double logZ = 0;                 // log of (sum over polyominoes containing the origin) = log(m A(m))
  std::vector<double> w(pop); std::vector<int> cnt(pop), freeSlots, gen; std::vector<u32> S; std::vector<double> pw; int nwv[4096];
  const double kpow[5] = {1, kappa, kappa * kappa, kappa * kappa * kappa, kappa * kappa * kappa * kappa}, lk = std::log(kappa);
  std::vector<int> special = sizes; special.insert(special.end(), cps.begin(), cps.end());
  double lastNorm = 0;                 // log of the mean weight right after the last resampling (0: all weights 1)
  for (int m = 1; m < nmax;) {
    int target = std::min(nmax, m + K);
    for (int sp : special) if (sp > m && sp < target) target = sp;
    for (auto &wk : W) {
      if (!(wk.lw > -INFINITY)) continue;
      for (u32 g : wk.cells) grid[g] = 1;
      for (u32 g : wk.blocked) grid[g] = 2;
      int mm = (int)wk.cells.size();
      while (mm < target) {
        S.clear(); gen.clear();
        for (int i = wk.head; i < mm; i++) for (int d = 0; d < 4; d++) {
          const size_t g = wk.cells[i] + nb[d];
          if (!grid[g]) { grid[g] = 3; S.push_back((u32)g); gen.push_back(i); }
        }
        const int G = (int)S.size();
        if (!G) { wk.lw = -INFINITY; break; }
        // choose s_j with probability proportional to kappa^(new_j - j), new_j = untried sites s_j would add
        pw.resize(G); double tot = 0, kj = 1;
        for (int i = 0; i < G; i++, kj /= kappa) {
          int nw = 0; for (int d = 0; d < 4; d++) nw += grid[S[i] + nb[d]] == 0;
          nwv[i & 4095] = nw; pw[i] = kj * kpow[nw]; tot += pw[i];
          if (kj < 1e-12 * tot) { for (int r = i + 1; r < G; r++) pw[r] = 0; break; }
        }
        double u = rng.uni() * tot; int j = 0;
        while (j < G - 1 && u >= pw[j]) u -= pw[j++];
        while (pw[j] == 0 && j > 0) j--;
        int nwj = 0; for (int d = 0; d < 4; d++) nwj += grid[S[j] + nb[d]] == 0;
        for (int i = 0; i < j; i++) { grid[S[i]] = 2; wk.blocked.push_back(S[i]); }
        grid[S[j]] = 1; wk.cells.push_back(S[j]); mm++; wk.head = gen[j];
        for (int i = j + 1; i < G; i++) grid[S[i]] = 0;
        wk.lw += std::log(tot / kappa); wk.G = G - j - 1 + nwj;
      }
      for (u32 g : wk.cells) grid[g] = 0;
      for (u32 g : wk.blocked) grid[g] = 0;
    }
    m = target;
    // the population follows counting measure tilted by kappa^G; lc removes the tilt
    double mx = -INFINITY, mxc = -INFINITY; for (auto &wk : W) { mx = std::max(mx, wk.lw); mxc = std::max(mxc, wk.lw - lk * wk.G); }
    double sum = 0, sum2 = 0, sumc = 0, sumc2 = 0; int alive = 0;
    for (int i = 0; i < pop; i++) { w[i] = std::exp(W[i].lw - mx); sum += w[i]; sum2 += w[i] * w[i]; alive += W[i].lw > -INFINITY;
      const double c = std::exp(W[i].lw - lk * W[i].G - mxc); sumc += c; sumc2 += c * c; }
    const double norm = mx + std::log(sum / pop), ess = sum * sum / sum2;     // log mean weight since the last resampling
    const double normc = mxc + std::log(sumc / pop) + lk * 4, essc = sumc * sumc / sumc2;   // same, for the counting measure (kappa^4: the start)
    if (fz) { fprintf(fz, "%d %.10f %d %.1f %.1f\n", m, logZ + normc, alive, ess, essc); fflush(fz); }
    if (std::find(sizes.begin(), sizes.end(), m) != sizes.end() && fo) {
      const int i10 = cpIndex(std::max(K, m / 10 / K * K)), i2 = cpIndex(std::max(K, m / 2 / K * K));
      std::unordered_set<u32> a10, a2;
      for (int i = 0; i < pop; i++) {
        if (!(W[i].lw > -INFINITY)) continue;
        for (int k = 0; k < m; k++) { X[k] = W[i].cells[k] % L; Y[k] = W[i].cells[k] / L; }
        F.compute(X.data(), Y.data(), m, f);
        double hdr[4] = {(double)m, W[i].lw - lk * W[i].G - mxc, (double)W[i].anc[i10], (double)W[i].anc[i2]}; fwrite(hdr, 8, 4, fo); fwrite(f.data(), 8, f.size(), fo);
        a10.insert(W[i].anc[i10]); a2.insert(W[i].anc[i2]);
      }
      fflush(fo);
      fprintf(stderr, "m=%d  alive %d/%d  distinct ancestors at m/10: %zu, at m/2: %zu  log(mA)=%.6f  ESS %.0f  %.0fs\n", m, alive, pop, a10.size(), a2.size(), logZ + normc, essc, now() - t0);
    }
    if (m == nmax) { logZ += normc; break; }
    if (ess > essFrac * pop && cpIndex(m) < 0) continue;
    logZ += norm; (void)lastNorm;
    // systematic resampling to equal weights
    double u = rng.uni(), acc = 0; int made = 0;
    for (int i = 0; i < pop; i++) { acc += w[i] / sum * pop; int c = 0; while (made + u < acc && made < pop) { c++; made++; } cnt[i] = c; }
    for (int i = pop - 1; made < pop; i--) if (w[i] > 0) { cnt[i]++; made++; }
    freeSlots.clear(); for (int i = 0; i < pop; i++) if (cnt[i] == 0) freeSlots.push_back(i);
    for (int i = 0; i < pop; i++) while (cnt[i] > 1) { W[freeSlots.back()] = W[i]; freeSlots.pop_back(); cnt[i]--; }
    for (auto &wk : W) wk.lw = 0;
    const int ci = cpIndex(m);
    if (ci >= 0) for (int i = 0; i < pop; i++) W[i].anc[ci] = i;
  }
  if (fo) fclose(fo);
  if (fz) fclose(fz);
  fprintf(stderr, "nmax=%d pop=%d step=%d  log(n A(n))=%.6f  A(n)~%.6g  %.1fs\n", nmax, pop, K, logZ, std::exp(logZ) / nmax, now() - t0);
  return 0;
}
