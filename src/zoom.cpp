// zoom: refine a polyomino to twice the resolution.
//
// Every cell is replaced by a 2x2 block and the result is relaxed by a walk whose moves are all
// local: a cell is added next to the shape or taken away from it, the number of cells being free
// with weight z per cell.  At z = 1/lambda (lambda = 4.0626, the growth constant of the number of
// polyominoes) this is how a piece of a very large uniformly random polyomino sees the rest of
// it: the conditional law of the cells inside any region, given the cells outside, is
// proportional to z^(cells inside) on the fillings that keep the whole connected.  The walk
// leaves that law invariant, so run for a limited time it rebuilds the fine detail and leaves
// the larger structure where the coarser shape put it.
//
// Move: pick a cell c of the movable region and one of 8 choices.  Four of them name a
// direction d: if s = c + d is empty and movable, continue with probability 1/k (k = movable
// cells next to s) and add s with probability min(1, 4 z n / (n + 1)).  The other four propose
// to remove c, accepted with probability min(1, n / (4 z (n - 1))), provided c has a movable
// neighbour.  Both need a certificate: the cells next to the site must be connected to each
// other, without the site, inside the box of radius K around it, where P(K >= r) = 1/r, capped
// at -kmax (default 256).  The certificate is the same for a move and its reverse, so the walk is reversible
// with respect to z^n on connected shapes, and a removal can never disconnect.  Loops wider
// than the cap are never opened or closed, and doubling widens every loop, so too small a cap
// lets wide loops pile up over repeated refinements: with 64, four refinements to a million
// cells left the largest clump a third too large on average.
//
//   zoom -load F [-inflate 2] [-z Z] [-sweeps T] [-kmax K] [-seed S] [-out F2] [-log E]
//        [-win x0 y0 x1 y1]     only cells with x0 <= x < x1, y0 <= y < y1 (after inflation) move
//        [-crop m]              with -win: write only the cells within m of the window
//        [-snaps PREFIX -snapAt "t1 t2 ..."]   save the shape at those sweeps
//        [-target N [-eps E]]   afterwards steer the cell count to exactly N, with z raised or
//                               lowered by the fraction E (default 0.01) until it gets there
//   zoom -mode exact -z Z -moves M [-kmax K] [-nmax 8]     test: shape frequencies against z^n
//
// A sweep is one proposed move per movable cell.  Shapes are int32 n, then n pairs x y.
#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
typedef uint64_t u64;
typedef uint32_t u32;

struct Rng {
  u64 s[4];
  static u64 splitmix(u64 &x) {
    u64 z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
  explicit Rng(u64 seed) { for (auto &v : s) v = splitmix(seed); }
  static inline u64 rotl(u64 x, int k) { return (x << k) | (x >> (64 - k)); }
  inline u64 next() {
    u64 r = rotl(s[0] + s[3], 23) + s[0], t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
    return r;
  }
  inline double uni() { return (next() >> 11) * (1.0 / 9007199254740992.0); }
  inline u32 below(u32 n) { return (u32)(((next() >> 32) * (u64)n) >> 32); }
};

struct Zoom {
  int W = 0, H = 0, ox = 0, oy = 0;          // grid size; grid (0,0) is lattice (ox, oy)
  std::vector<uint8_t> g;                    // bit 0: occupied; bits 1-3: search label + 1
  std::vector<int> cells;                    // grid index of every movable cell
  int wx0, wy0, wx1, wy1;                    // movable region, grid coordinates
  double z = 1.0 / 4.0625696; int kmax = 256;
  Rng rng{1};
  u64 nAdd = 0, nRem = 0, nCertFail = 0, work = 0;
  std::vector<int> q[4]; int par[4];
  int find(int a) { while (par[a] != a) a = par[a] = par[par[a]]; return a; }

  inline bool occ(int p) const { return g[p] & 1; }
  inline bool movable(int p) const { int x = p % W, y = p / W; return x >= wx0 && x < wx1 && y >= wy0 && y < wy1; }

  // are the occupied neighbours of site p connected to each other, without p, inside the box of
  // radius K around p?
  bool cert(int p, int K) {
    static const int RX[8] = {1,1,0,-1,-1,-1,0,1}, RY[8] = {0,1,1,1,0,-1,-1,-1};
    int px = p % W, py = p / W, m = 0;
    for (int i = 0; i < 8; i++) if (occ(p + RX[i] + RY[i] * W)) m |= 1 << i;
    int deg = ((m >> 0) & 1) + ((m >> 2) & 1) + ((m >> 4) & 1) + ((m >> 6) & 1);
    if (deg <= 1) return true;
    int links = 0;
    for (int i = 0; i < 8; i += 2)
      if ((m >> i & 1) && (m >> ((i + 1) & 7) & 1) && (m >> ((i + 2) & 7) & 1)) links++;
    if (links >= deg - 1) return true;
    if (K <= 1) return false;
    int x0 = std::max(1, px - K), x1 = std::min(W - 2, px + K), y0 = std::max(1, py - K), y1 = std::min(H - 2, py + K);
    bool was = occ(p); g[p] &= ~1;           // searches must not pass through p
    int ng = 0;
    for (int i = 0; i < 8; i += 2) if (m >> i & 1) { int j = p + RX[i] + RY[i] * W; q[ng].clear(); q[ng].push_back(j); g[j] |= (ng + 1) << 1; par[ng] = ng; ng++; }
    int ncls = ng; size_t head[4] = {0, 0, 0, 0}; bool ok = false;
    static const int DX[4] = {1, 0, -1, 0}, DY[4] = {0, 1, 0, -1};
    while (true) {
      bool any = false, dead = false;
      for (int a = 0; a < ng && !ok && !dead; a++) {
        if (head[a] >= q[a].size()) continue;
        int j = q[a][head[a]++]; work++; any = true;
        int x = j % W, y = j / W;
        for (int d = 0; d < 4; d++) {
          int ax = x + DX[d], ay = y + DY[d];
          if (ax < x0 || ax > x1 || ay < y0 || ay > y1) continue;
          int b = ax + ay * W; uint8_t v = g[b];
          if (!(v & 1)) continue;
          int lb = (v >> 1) & 7;
          if (!lb) { g[b] |= (a + 1) << 1; q[a].push_back(b); }
          else { int ra = find(a), rb = find(lb - 1); if (ra != rb) { par[ra] = rb; if (--ncls == 1) ok = true; } }
        }
        // a class whose searches have all run dry can never join the others
        if (!ok) {
          int ra = find(a); bool alive = false;
          for (int c = 0; c < ng; c++) if (find(c) == ra && head[c] < q[c].size()) alive = true;
          if (!alive) dead = true;
        }
      }
      if (ok || dead || !any) break;
    }
    for (int a = 0; a < ng; a++) for (int j : q[a]) g[j] &= 1;
    if (was) g[p] |= 1;
    return ok;
  }

  inline int drawK() { double u = rng.uni(); return u * kmax < 1.0 ? kmax : (int)(1.0 / u); }

  void move() {
    int n = (int)cells.size();
    u32 i = rng.below(n); int c = cells[i]; u32 r = rng.next() >> 61;
    static const int DX[4] = {1, 0, -1, 0}, DY[4] = {0, 1, 0, -1};
    if (r < 4) {
      int s = c + DX[r] + DY[r] * W;
      if (occ(s) || !movable(s)) return;
      int k = 0; for (int d = 0; d < 4; d++) { int t = s + DX[d] + DY[d] * W; if (occ(t) && movable(t)) k++; }
      if (k > 1 && rng.below(k) != 0) return;
      double a = 4.0 * z * n / (n + 1.0);
      if (a < 1.0 && rng.uni() >= a) return;
      if (!cert(s, drawK())) { nCertFail++; return; }
      g[s] |= 1; cells.push_back(s); nAdd++;
    } else {
      if (n <= 1) return;
      int k = 0; for (int d = 0; d < 4; d++) { int t = c + DX[d] + DY[d] * W; if (occ(t) && movable(t)) k++; }
      if (!k) return;
      double a = n / (4.0 * z * (n - 1.0));
      if (a < 1.0 && rng.uni() >= a) return;
      if (!cert(c, drawK())) { nCertFail++; return; }
      g[c] &= ~1; cells[i] = cells.back(); cells.pop_back(); nRem++;
    }
  }
};

static std::vector<int> readShape(const std::string &fn) {
  FILE *f = fopen(fn.c_str(), "rb"); if (!f) { fprintf(stderr, "cannot read %s\n", fn.c_str()); exit(1); }
  int32_t n; if (fread(&n, 4, 1, f) != 1) exit(1);
  std::vector<int> v(2 * (size_t)n); if (fread(v.data(), 4, 2 * (size_t)n, f) != 2 * (size_t)n) exit(1);
  fclose(f); return v;
}

int main(int argc, char **argv) {
  std::string load, out, mode = "zoom", snaps, snapAt; int inflate = 1, nmax = 8, crop = -1; double sweeps = 0, logEvery = 0, moves = 1e8;
  u64 seed = 1; bool haveWin = false; long target = 0; double eps = 0.01; long win[4] = {0, 0, 0, 0}; Zoom Z;
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i]; auto nxt = [&]() { if (i + 1 >= argc) { fprintf(stderr, "missing value for %s\n", a.c_str()); exit(1); } return std::string(argv[++i]); };
    if (a == "-load") load = nxt(); else if (a == "-out") out = nxt(); else if (a == "-inflate") inflate = atoi(nxt().c_str());
    else if (a == "-z") Z.z = atof(nxt().c_str()); else if (a == "-sweeps") sweeps = atof(nxt().c_str());
    else if (a == "-kmax") Z.kmax = atoi(nxt().c_str()); else if (a == "-seed") seed = strtoull(nxt().c_str(), 0, 10);
    else if (a == "-log") logEvery = atof(nxt().c_str()); else if (a == "-mode") mode = nxt(); else if (a == "-moves") moves = atof(nxt().c_str());
    else if (a == "-nmax") nmax = atoi(nxt().c_str()); else if (a == "-crop") crop = atoi(nxt().c_str());
    else if (a == "-target") target = atol(nxt().c_str()); else if (a == "-eps") eps = atof(nxt().c_str());
    else if (a == "-snaps") snaps = nxt(); else if (a == "-snapAt") snapAt = nxt();
    else if (a == "-win") { for (int k = 0; k < 4; k++) win[k] = atol(nxt().c_str()); haveWin = true; }
    else { fprintf(stderr, "unknown option %s\n", a.c_str()); return 1; }
  }
  Z.rng = Rng(seed);

  if (mode == "exact") {
    // a single cell grown at subcritical z: every shape of n cells should turn up in proportion to z^n
    int S = 512; Z.W = Z.H = S; Z.g.assign((size_t)S * S, 0); Z.wx0 = Z.wy0 = 2; Z.wx1 = Z.wy1 = S - 2;
    int c0 = S / 2 + (S / 2) * S; Z.g[c0] = 1; Z.cells.push_back(c0);
    std::map<std::string, u64> hist; std::vector<u64> byN(nmax + 2, 0); u64 every = 200, samples = 0;
    for (u64 m = 0; m < (u64)moves; m++) {
      // keep the shape away from the edge: translate it back to the middle when it drifts
      Z.move();
      if (m % every) continue;
      int n = (int)Z.cells.size(), xmin = S, ymin = S, xmax = 0, ymax = 0;
      for (int p : Z.cells) { xmin = std::min(xmin, p % S); xmax = std::max(xmax, p % S); ymin = std::min(ymin, p / S); ymax = std::max(ymax, p / S); }
      if (xmin < 8 || ymin < 8 || xmax > S - 9 || ymax > S - 9) {
        int dx = S / 2 - (xmin + xmax) / 2, dy = S / 2 - (ymin + ymax) / 2;
        for (int p : Z.cells) Z.g[p] = 0;
        for (int &p : Z.cells) p += dx + dy * S;
        for (int p : Z.cells) Z.g[p] = 1;
        continue;
      }
      samples++;
      if (n > nmax) { byN[nmax + 1]++; continue; }
      byN[n]++;
      std::vector<std::pair<int,int>> c; for (int p : Z.cells) c.push_back({p % S - xmin, p / S - ymin});
      std::sort(c.begin(), c.end());
      std::string key; for (auto &pp : c) { key.push_back((char)pp.first); key.push_back((char)pp.second); }
      hist[key]++;
    }
    // Poisson-style chi-square of the shape counts against A z^n, with A fitted
    std::vector<u64> shapesN(nmax + 1, 0); double num = 0, den = 0;
    for (auto &kv : hist) { int n = (int)kv.first.size() / 2; shapesN[n]++; num += kv.second; den += std::pow(Z.z, n); }
    static const u64 AN[] = {0, 1, 2, 6, 19, 63, 216, 760, 2725, 9910, 36446};
    for (int n = 1; n <= nmax; n++) den += (AN[n] - shapesN[n]) * std::pow(Z.z, n);   // shapes never seen count too
    double A = num / den, chi2 = 0; u64 dof = 0;
    for (auto &kv : hist) { double e = A * std::pow(Z.z, (int)kv.first.size() / 2); chi2 += (kv.second - e) * (kv.second - e) / e; dof++; }
    for (int n = 1; n <= nmax; n++) { chi2 += (AN[n] - shapesN[n]) * A * std::pow(Z.z, n); dof += AN[n] - shapesN[n]; }
    dof -= 1;
    printf("z=%.5f kmax=%d samples=%" PRIu64 " (every %" PRIu64 " moves; correlated, so chi2 is indicative)\n", Z.z, Z.kmax, samples, every);
    for (int n = 1; n <= nmax; n++) printf("  n=%d shapes seen %" PRIu64 " of %" PRIu64 "  share of samples %.5f\n", n, shapesN[n], AN[n], byN[n] / (double)samples);
    printf("chi2=%.1f dof=%" PRIu64 " z=%.2f\n", chi2, dof, (chi2 - dof) / std::sqrt(2.0 * dof));
    // ratios of the totals at successive n, which should be z a_n / a_(n-1)
    for (int n = 2; n <= nmax; n++) printf("  count(%d)/count(%d) = %.5f  expected %.5f\n", n, n - 1, byN[n] / (double)byN[n - 1], Z.z * AN[n] / AN[n - 1]);
    return 0;
  }

  std::vector<int> v = readShape(load); size_t n0 = v.size() / 2;
  long xmin = 1L << 40, ymin = 1L << 40, xmax = -(1L << 40), ymax = -(1L << 40);
  for (size_t i = 0; i < n0; i++) { xmin = std::min<long>(xmin, v[2 * i]); xmax = std::max<long>(xmax, v[2 * i]); ymin = std::min<long>(ymin, v[2 * i + 1]); ymax = std::max<long>(ymax, v[2 * i + 1]); }
  int margin = haveWin ? 4 : 256;
  xmin *= inflate; ymin *= inflate; xmax = xmax * inflate + inflate - 1; ymax = ymax * inflate + inflate - 1;
  Z.ox = (int)xmin - margin; Z.oy = (int)ymin - margin; Z.W = (int)(xmax - xmin + 1 + 2 * margin); Z.H = (int)(ymax - ymin + 1 + 2 * margin);
  if ((double)Z.W * Z.H > 2.0e9) { fprintf(stderr, "grid too large\n"); return 1; }
  Z.g.assign((size_t)Z.W * Z.H, 0);
  if (haveWin) { Z.wx0 = std::max<long>(2, win[0] - Z.ox); Z.wy0 = std::max<long>(2, win[1] - Z.oy); Z.wx1 = std::min<long>(Z.W - 2, win[2] - Z.ox); Z.wy1 = std::min<long>(Z.H - 2, win[3] - Z.oy); }
  else { Z.wx0 = Z.wy0 = 2; Z.wx1 = Z.W - 2; Z.wy1 = Z.H - 2; }
  for (size_t i = 0; i < n0; i++)
    for (int a = 0; a < inflate; a++) for (int b = 0; b < inflate; b++) {
      int p = (v[2 * i] * inflate + a - Z.ox) + (v[2 * i + 1] * inflate + b - Z.oy) * Z.W;
      Z.g[p] = 1; if (Z.movable(p)) Z.cells.push_back(p);
    }
  size_t total0 = n0 * inflate * inflate, fixed = total0 - Z.cells.size();
  auto save = [&](const std::string &fn) {
    std::vector<int32_t> o; o.push_back(0);
    long cx0 = Z.wx0 - crop, cx1 = Z.wx1 + crop, cy0 = Z.wy0 - crop, cy1 = Z.wy1 + crop;
    for (int y = 0; y < Z.H; y++) for (int x = 0; x < Z.W; x++) if (Z.g[x + (size_t)y * Z.W] & 1) {
      if (crop >= 0 && (x < cx0 || x >= cx1 || y < cy0 || y >= cy1)) continue;
      o.push_back(x + Z.ox); o.push_back(y + Z.oy);
    }
    o[0] = (int32_t)(o.size() / 2);
    FILE *f = fopen(fn.c_str(), "wb"); fwrite(o.data(), 4, o.size(), f); fclose(f);
  };
  std::vector<double> snapT; { char *e; const char *s = snapAt.c_str(); while (*s) { double t = strtod(s, &e); if (e == s) break; snapT.push_back(t); s = e; } }
  size_t si = 0; double t = 0, nextLog = 0;
  printf("# start: %zu cells (%zu movable) from %zu, z=%.7f kmax=%d\n# sweeps cells movable added removed cert_fail work_per_move\n", total0, Z.cells.size(), n0, Z.z, Z.kmax);
  u64 done = 0;
  while (true) {
    if (logEvery > 0 && t >= nextLog - 1e-9) {
      printf("%.3f %zu %zu %" PRIu64 " %" PRIu64 " %" PRIu64 " %.2f\n", t, Z.cells.size() + fixed, Z.cells.size(), Z.nAdd, Z.nRem, Z.nCertFail, done ? Z.work / (double)done : 0.0); fflush(stdout);
      nextLog = nextLog < 1 ? 1 : nextLog < 64 ? nextLog * 2 : nextLog + logEvery;
      if (nextLog < logEvery && logEvery < 1) nextLog = t + logEvery;
    }
    while (si < snapT.size() && t >= snapT[si] - 1e-9) { char b[64]; snprintf(b, 64, "_t%g.bin", snapT[si]); save(snaps + b); si++; }
    if (t >= sweeps - 1e-9) break;
    // a quarter of a sweep at a time, so the early thinning is resolved
    u64 batch = std::max<u64>(1, Z.cells.size() / 4);
    for (u64 k = 0; k < batch; k++) Z.move();
    done += batch; t += 0.25;
  }
  double steer = 0;
  if (target > 0) {
    double z0 = Z.z; size_t before = Z.cells.size() + fixed;
    while ((long)(Z.cells.size() + fixed) != target && steer < 1e4) {
      Z.z = z0 * ((long)(Z.cells.size() + fixed) < target ? 1 + eps : 1 - eps);
      u64 batch = std::max<u64>(1, Z.cells.size() / 64);
      for (u64 k = 0; k < batch && (long)(Z.cells.size() + fixed) != target; k++) Z.move();
      steer += 1.0 / 64;
    }
    Z.z = z0; fprintf(stderr, "zoom: steered %zu -> %zu cells in %.2f sweeps at z(1 +- %g)\n", before, Z.cells.size() + fixed, steer, eps);
  }
  if (!out.empty()) save(out);
  fprintf(stderr, "zoom: %zu -> %zu cells (x%.4f of the original %zu) in %.2f sweeps; added %" PRIu64 " removed %" PRIu64 " cert-fail %" PRIu64 "\n",
          total0, Z.cells.size() + fixed, (Z.cells.size() + fixed) / (double)n0, n0, t, Z.nAdd, Z.nRem, Z.nCertFail);
  return 0;
}
