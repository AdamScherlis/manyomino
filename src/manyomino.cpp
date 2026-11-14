// manyomino — uniform sampling of very large fixed polyominoes (site animals on Z^2).
//
// Three symmetric Markov kernels on the set of n-cell polyominoes, mixed at random.
// Each is symmetric (K(P,P') = K(P',P)), so the uniform distribution is stationary for
// any mixture, and RELOCATE alone is irreducible.
//
//  RELOCATE   pick a cell c, a cell u and a direction d uniformly, and a radius K from a
//             fixed heavy-tailed law (P(K >= r) = r^-2).  Let s = u+d and Q = P\c.
//             If s is empty, accept with probability 1/k, k = #neighbours of s in Q.
//             Require both c and s to be "K-locally connected holes" of Q: the
//             neighbours of the site in Q are joined to each other inside the Chebyshev
//             ball of radius K about it.  Move c to s.  The condition is symmetric in
//             (c,s), every valid pair has probability 1/(4n^2) given K, and K-local
//             connectivity at c implies Q is connected.  K exceeds the diameter with
//             positive probability, which gives the unrestricted (irreducible) kernel.
//
//  TRANSPLANT pick a cell c, a direction d, a cell c' and a symmetry g in D4 uniformly.
//             Let B be the component of P\c containing c+d; require 2|B| <= n-1 and
//             c' not in B.  Let T be the rigid motion x -> c' + g(x-c).  Require that
//             T(B) is disjoint from A = P\B and adjacent to A only at c'.  Replace B by
//             T(B).  The reverse move is (c', g(d), c, g^-1), same probability 1/(32n^2).
//
//  PIVOT      the same with c' = c and g != identity.
//
// Implementation note: a TRANSPLANT is decided by breadth-first searches run in lockstep
// from the neighbours of c and from c'.  The images of B's cells are tested as they are
// discovered, and the move is abandoned as soon as it is *certain* to be invalid, so a
// typical rejected proposal costs O(1) rather than O(|B|).  This changes cost only: the
// accept/reject outcome is identical to the naive implementation (ref/manyomino_v1.cpp
// produces bit-identical trajectories from the same seed).
//
// Optional tilt: weight q^(site perimeter) via Metropolis (q=1 is uniform).  With
// q = 1-p this is exactly the law of a size-n site-percolation cluster at density p.
//
// Build:  g++ -O3 -march=native -std=c++17 -o manyomino src/manyomino.cpp

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>
#include <chrono>

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
  inline u32 below(u32 n) { return (u32)(((next() >> 32) * (u64)n) >> 32); }
  inline double uni() { return (next() >> 11) * (1.0 / 9007199254740992.0); }
};

static const int DX[4] = {1, 0, -1, 0}, DY[4] = {0, 1, 0, -1};
// D4 as 2x2 integer matrices (a b; c d)
static const int G[8][4] = {{1,0,0,1},{0,-1,1,0},{-1,0,0,-1},{0,1,-1,0},
                            {1,0,0,-1},{0,1,1,0},{-1,0,0,1},{0,-1,-1,0}};

struct Sampler {
  int n, W, mask, lg;
  // Sparse torus grid: W x W sites in 64x64 blocks allocated on demand, so memory is
  // O(n) while lookups stay cache-friendly.  Entry 0 = empty, else cell index + 1.
  static const int BL = 6, BS = 1 << BL, BM = BS - 1;
  int tl, tm;                       // top-level table is (W/64)^2
  std::vector<u32*> top;
  std::vector<u32> bcount;
  std::vector<u32*> freeBlocks;
  std::vector<int32_t> X, Y;        // unwrapped coordinates of each cell
  std::vector<u32> mark;            // BFS stamp per cell
  std::vector<uint8_t> lab;         // BFS group label per cell
  u32 stamp = 0;
  std::vector<int> q[5];
  size_t qh[5];
  std::vector<u32> pend;            // stamp: cell is a pending collision
  std::vector<int> pendList;
  int relK = 1;                     // 1: random-radius relocate, 0: always exact
  bool refCP = false;               // use the naive reference cut&paste
  Rng rng;
  double logq = 0.0;                // perimeter tilt: weight exp(logq * t)
  // statistics
  u64 nRel = 0, accRel = 0, nCP = 0, accCP = 0, bfsWork = 0, cpMoved = 0;
  u64 cpHist[40] = {0}, cpAccHist[40] = {0};
  int smax;                         // cap on branch size (symmetric rule), default (n-1)/2

  Sampler(int n_, int lgW, u64 seed) : n(n_), W(1 << lgW), mask((1 << lgW) - 1), lg(lgW),
      tl(lgW - BL), tm((1 << (lgW - BL)) - 1), top((size_t)1 << (2 * (lgW - BL)), nullptr),
      bcount((size_t)1 << (2 * (lgW - BL)), 0), X(n_), Y(n_), mark(n_, 0), lab(n_, 0), pend(n_, 0), rng(seed) {
    smax = (n - 1) / 2;
  }
  inline size_t tix(int x, int y) const { return ((size_t)((y >> BL) & tm) << tl) | (size_t)((x >> BL) & tm); }
  inline u32 occ(int x, int y) const {
    const u32 *b = top[tix(x, y)];
    return b ? b[((y & BM) << BL) | (x & BM)] : 0;
  }
  inline void setCell(int x, int y, u32 v) {      // v != 0, site currently empty
    size_t t = tix(x, y);
    u32 *b = top[t];
    if (!b) {
      if (freeBlocks.empty()) b = (u32*)calloc(BS * BS, sizeof(u32));
      else { b = freeBlocks.back(); freeBlocks.pop_back(); }
      top[t] = b;
    }
    b[((y & BM) << BL) | (x & BM)] = v; bcount[t]++;
  }
  inline void clearCell(int x, int y) {           // site currently occupied
    size_t t = tix(x, y);
    u32 *b = top[t];
    b[((y & BM) << BL) | (x & BM)] = 0;
    if (--bcount[t] == 0) { freeBlocks.push_back(b); top[t] = nullptr; }
  }

  void place(const std::vector<std::pair<int,int>> &cells) {
    assert((int)cells.size() == n);
    for (int i = 0; i < n; i++) {
      X[i] = cells[i].first; Y[i] = cells[i].second;
      assert(!occ(X[i], Y[i]));
      setCell(X[i], Y[i], i + 1);
    }
  }

  // --- connectivity race -----------------------------------------------------------
  // BFS simultaneously from the occupied neighbours of cell ic (which is blocked).
  // mode 0: return 1 iff P\c is connected (stop as soon as decided).
  // mode 1: run until at most one class is unexhausted; fills gcls[], exhausted[].
  int par[5], ngrp, gdir[5];
  bool exhausted[5];
  int find(int a) { while (par[a] != a) a = par[a] = par[par[a]]; return a; }

  int race(int ic, int mode) {
    int cx = X[ic], cy = Y[ic];
    stamp++;
    mark[ic] = stamp; lab[ic] = 255;
    ngrp = 0;
    for (int d = 0; d < 4; d++) {
      u32 o = occ(cx + DX[d], cy + DY[d]);
      if (!o) continue;
      int j = o - 1, g = ngrp++;
      q[g].clear(); qh[g] = 0; q[g].push_back(j);
      mark[j] = stamp; lab[j] = g; par[g] = g; gdir[g] = d; exhausted[g] = false;
    }
    int ncls = ngrp, nlive = ngrp;   // classes; classes not yet exhausted
    if (ngrp <= 1) return 1;
    while (true) {
      for (int g = 0; g < ngrp; g++) {
        if (qh[g] >= q[g].size()) continue;
        int j = q[g][qh[g]++];
        bfsWork++;
        int x = X[j], y = Y[j];
        for (int d = 0; d < 4; d++) {
          u32 o = occ(x + DX[d], y + DY[d]);
          if (!o) continue;
          int k = o - 1;
          if (mark[k] != stamp) { mark[k] = stamp; lab[k] = g; q[g].push_back(k); }
          else if (lab[k] != 255) {
            int a = find(lab[k]), b = find(g);
            if (a != b) { par[a] = b; ncls--; nlive--; if (ncls == 1) return 1; if (mode == 1 && nlive <= 1) return 0; }
          }
        }
        if (qh[g] >= q[g].size()) {
          // is the whole class of g exhausted?
          int r = find(g); bool ex = true;
          for (int h = 0; h < ngrp; h++) if (find(h) == r && qh[h] < q[h].size()) { ex = false; break; }
          if (ex) {
            if (mode == 0) return 0;      // a complete component that isn't everything
            for (int h = 0; h < ngrp; h++) if (find(h) == r) exhausted[h] = true;
            nlive--;
            if (nlive <= 1) return 0;
          }
        }
      }
    }
  }

  // local 8-ring test: 1 if occupied 4-neighbours of (cx,cy) are connected through the ring
  inline bool locallyConnected(int cx, int cy, int &deg) const {
    static const int RX[8] = {1,1,0,-1,-1,-1,0,1}, RY[8] = {0,1,1,1,0,-1,-1,-1};
    int m = 0;
    for (int i = 0; i < 8; i++) if (occ(cx + RX[i], cy + RY[i])) m |= 1 << i;
    deg = ((m >> 0) & 1) + ((m >> 2) & 1) + ((m >> 4) & 1) + ((m >> 6) & 1);
    if (deg <= 1) return true;
    // count gaps between consecutive occupied edge-neighbours around the ring:
    // edge neighbours i=0,2,4,6; i and i+2 are locally joined iff corner i+1 is occupied too.
    int links = 0;
    for (int i = 0; i < 8; i += 2)
      if ((m >> i & 1) && (m >> ((i + 1) & 7) & 1) && (m >> ((i + 2) & 7) & 1)) links++;
    return links >= deg - 1;
  }

  // --- perimeter bookkeeping for the tilt -----------------------------------------
  // change in site perimeter if cell at (x,y) is removed (grid must still contain it)
  int dPerimRemove(int x, int y) const {
    int dt = 0, deg = 0;
    for (int d = 0; d < 4; d++) {
      int ax = x + DX[d], ay = y + DY[d];
      if (occ(ax, ay)) { deg++; continue; }
      // empty neighbour: stays perimeter iff it has another occupied neighbour
      int k = 0;
      for (int e = 0; e < 4; e++) if (occ(ax + DX[e], ay + DY[e])) k++;
      if (k == 1) dt--;
    }
    if (deg > 0) dt++;                 // the vacated site becomes perimeter
    return dt;
  }
  int dPerimAdd(int x, int y) const {  // (x,y) empty with >=1 occupied neighbour
    int dt = -1;
    for (int d = 0; d < 4; d++) {
      int ax = x + DX[d], ay = y + DY[d];
      if (occ(ax, ay)) continue;
      int k = 0;
      for (int e = 0; e < 4; e++) if (occ(ax + DX[e], ay + DY[e])) k++;
      if (k == 0) dt++;
    }
    return dt;
  }

  // Are the neighbours (in P minus cell `ex`) of site (cx,cy) joined to each other inside
  // the Chebyshev ball of radius K about it?  (cx,cy) itself is treated as empty.
  bool connAt(int cx, int cy, int ex, int K) {
    static const int RX[8] = {1,1,0,-1,-1,-1,0,1}, RY[8] = {0,1,1,1,0,-1,-1,-1};
    int m = 0;
    for (int i = 0; i < 8; i++) { u32 o = occ(cx + RX[i], cy + RY[i]); if (o && (int)o - 1 != ex) m |= 1 << i; }
    int deg = ((m >> 0) & 1) + ((m >> 2) & 1) + ((m >> 4) & 1) + ((m >> 6) & 1);
    if (deg <= 1) return true;
    int links = 0;
    for (int i = 0; i < 8; i += 2)
      if ((m >> i & 1) && (m >> ((i + 1) & 7) & 1) && (m >> ((i + 2) & 7) & 1)) links++;
    if (links >= deg - 1) return true;
    if (K <= 1) return false;
    // lockstep BFS inside the ball
    stamp++;
    if (ex >= 0) { mark[ex] = stamp; lab[ex] = 255; }
    int ng = 0;
    for (int d = 0; d < 4; d++) {
      u32 o = occ(cx + DX[d], cy + DY[d]);
      if (!o || (int)o - 1 == ex) continue;
      int j = o - 1, g = ng++;
      q[g].clear(); qh[g] = 0; q[g].push_back(j);
      mark[j] = stamp; lab[j] = g; par[g] = g;
    }
    int ncls = ng;
    bool ball = K < (1 << 28);
    while (true) {
      for (int g = 0; g < ng; g++) {
        if (qh[g] >= q[g].size()) continue;
        int j = q[g][qh[g]++];
        bfsWork++;
        int x = X[j], y = Y[j];
        for (int d = 0; d < 4; d++) {
          int ax = x + DX[d], ay = y + DY[d];
          if (ball && (std::abs(ax - cx) > K || std::abs(ay - cy) > K)) continue;
          if (ax == cx && ay == cy) continue;
          u32 o = occ(ax, ay);
          if (!o) continue;
          int k = o - 1;
          if (mark[k] != stamp) { mark[k] = stamp; lab[k] = g; q[g].push_back(k); }
          else if (lab[k] != 255) {
            int a = find(lab[k]), b = find(g);
            if (a != b) { par[a] = b; if (--ncls == 1) return true; }
          }
        }
        if (qh[g] >= q[g].size()) {
          int r = find(g); bool ex2 = true;
          for (int h = 0; h < ng; h++) if (find(h) == r && qh[h] < q[h].size()) { ex2 = false; break; }
          if (ex2) return false;
        }
      }
    }
  }

  // --- RELOCATE ---------------------------------------------------------------------
  void relocate() {
    nRel++;
    int ic = rng.below(n), iu = rng.below(n), d = rng.below(4);
    int K = 1 << 30;
    if (relK) { double u = rng.uni(); K = u < 1e-12 ? (1 << 30) : (int)std::min(1.0e9, 1.0 / std::sqrt(u)); }
    if (iu == ic) return;
    int sx = X[iu] + DX[d], sy = Y[iu] + DY[d];
    int cx = X[ic], cy = Y[ic];
    if (occ(sx, sy)) return;                    // (s == c is occupied -> null move, fine)
    int k = 0;
    for (int e = 0; e < 4; e++) {
      u32 a = occ(sx + DX[e], sy + DY[e]);
      if (a && (int)a - 1 != ic) k++;
    }
    if (k > 1 && rng.below(k) != 0) return;
    if (!connAt(cx, cy, ic, K)) return;
    if (relK && k > 1 && !connAt(sx, sy, ic, K)) return;
    if (logq != 0.0) {
      int dt = dPerimRemove(cx, cy);
      clearCell(cx, cy);
      dt += dPerimAdd(sx, sy);
      if (dt != 0 && rng.uni() >= std::exp(logq * dt)) { setCell(cx, cy, ic + 1); return; }
    } else clearCell(cx, cy);
    X[ic] = sx; Y[ic] = sy; setCell(sx, sy, ic + 1);
    accRel++;
  }

  // --- CUT & PASTE ------------------------------------------------------------------
  std::vector<int> B;
  std::vector<int32_t> nx, ny;
  void cutpasteRef(bool pivot) {
    nCP++;
    int ic = rng.below(n), d = rng.below(4), ia = pivot ? ic : (int)rng.below(n), g = pivot ? 1 + (int)rng.below(7) : (int)rng.below(8);
    int cx = X[ic], cy = Y[ic];
    u32 o = occ(cx + DX[d], cy + DY[d]);
    if (!o) return;
    int deg;
    if (locallyConnected(cx, cy, deg)) return;  // c is not a cut cell: B would be P\c
    if (race(ic, 1)) return;                    // connected after all
    // which group is direction d?
    int gd = -1;
    for (int h = 0; h < ngrp; h++) if (gdir[h] == d) gd = h;
    int r = find(gd);
    size_t bsz = 0;
    if (exhausted[gd]) {
      for (int h = 0; h < ngrp; h++) if (find(h) == r) bsz += q[h].size();
      if ((int)bsz > smax) return;
      if (mark[ia] == stamp && lab[ia] != 255 && find(lab[ia]) == r) return;  // c' in B
      B.clear();
      for (int h = 0; h < ngrp; h++) if (find(h) == r) B.insert(B.end(), q[h].begin(), q[h].end());
    } else {
      size_t ex = 0;
      for (int h = 0; h < ngrp; h++) if (exhausted[h]) ex += q[h].size();
      bsz = (size_t)(n - 1) - ex;
      if ((int)bsz > smax) return;
      // rare: the unexhausted class is itself small enough; finish enumerating it
      bool more = true;
      while (more) {
        more = false;
        for (int h = 0; h < ngrp; h++) {
          if (find(h) != r) continue;
          while (qh[h] < q[h].size()) {
            int j = q[h][qh[h]++]; bfsWork++;
            for (int e = 0; e < 4; e++) {
              u32 a = occ(X[j] + DX[e], Y[j] + DY[e]);
              if (!a) continue;
              int k = a - 1;
              if (mark[k] != stamp) { mark[k] = stamp; lab[k] = h; q[h].push_back(k); more = true; }
            }
          }
        }
      }
      if (mark[ia] == stamp && lab[ia] != 255 && find(lab[ia]) == r) return;
      B.clear();
      for (int h = 0; h < ngrp; h++) if (find(h) == r) B.insert(B.end(), q[h].begin(), q[h].end());
      assert(B.size() == bsz);
    }
    // relabel B cells so membership is a single compare
    for (int j : B) lab[j] = 100;
    int hb = 0; { size_t s = bsz; while (s > 1) { s >>= 1; hb++; } }
    cpHist[hb]++;
    int ax = X[ia], ay = Y[ia];
    const int *M = G[g];
    if (ia == ic && g == 0) return;             // identity
    nx.resize(bsz); ny.resize(bsz);
    // validity: T(B) disjoint from A, and adjacent to A only at c'
    for (size_t i = 0; i < bsz; i++) {
      int j = B[i];
      int dx = X[j] - cx, dy = Y[j] - cy;
      int x = ax + M[0] * dx + M[1] * dy, y = ay + M[2] * dx + M[3] * dy;
      nx[i] = x; ny[i] = y;
      u32 a = occ(x, y);
      if (a) { int k = a - 1; if (!(mark[k] == stamp && lab[k] == 100)) return; }
      for (int e = 0; e < 4; e++) {
        u32 b = occ(x + DX[e], y + DY[e]);
        if (!b) continue;
        int k = b - 1;
        if (k == ia) continue;
        if (mark[k] == stamp && lab[k] == 100) continue;
        return;
      }
    }
    auto apply = [&]() {   // swap B between its current sites and (nx,ny)
      for (int j : B) clearCell(X[j], Y[j]);
      for (size_t i = 0; i < bsz; i++) { int j = B[i]; std::swap(X[j], nx[i]); std::swap(Y[j], ny[i]); }
      for (int j : B) setCell(X[j], Y[j], j + 1);
    };
    if (logq != 0.0) {
      int before = localPerim(B, nx, ny);
      apply();
      int dt = localPerim(B, nx, ny) - before;
      if (dt != 0 && rng.uni() >= std::exp(logq * dt)) { apply(); return; }
    } else apply();
    accCP++; cpMoved += bsz; cpAccHist[hb]++;
  }
  // Lockstep version: identical accept/reject decisions, early certified rejection.
  size_t imgi[5];
  void cutpaste(bool pivot) {
    if (refCP) { cutpasteRef(pivot); return; }
    nCP++;
    int ic = rng.below(n), d = rng.below(4), ia = pivot ? ic : (int)rng.below(n), g = pivot ? 1 + (int)rng.below(7) : (int)rng.below(8);
    int cx = X[ic], cy = Y[ic];
    if (!occ(cx + DX[d], cy + DY[d])) return;
    int deg;
    if (locallyConnected(cx, cy, deg)) return;  // c is not a cut cell
    stamp++;
    mark[ic] = stamp; lab[ic] = 255;
    ngrp = 0; int gd = -1;
    for (int e = 0; e < 4; e++) {
      u32 o = occ(cx + DX[e], cy + DY[e]);
      if (!o) continue;
      int j = o - 1, h = ngrp++;
      q[h].clear(); qh[h] = 0; imgi[h] = 0; q[h].push_back(j);
      mark[j] = stamp; lab[j] = h; par[h] = h; exhausted[h] = false;
      if (e == d) gd = h;
    }
    int ga = -1;
    if (ia != ic) {
      if (mark[ia] == stamp) ga = lab[ia];
      else {
        ga = ngrp++;
        q[ga].clear(); qh[ga] = 0; imgi[ga] = 0; q[ga].push_back(ia);
        mark[ia] = stamp; lab[ia] = ga; par[ga] = ga; exhausted[ga] = false;
      }
      if (ga == gd) return;
    }
    const int ax = X[ia], ay = Y[ia];
    const int *M = G[g];
    pendList.clear();
    bool done = false;
    while (!done) {
      // one BFS step for every group that still has a frontier
      for (int h = 0; h < ngrp; h++) {
        if (qh[h] >= q[h].size()) continue;
        int j = q[h][qh[h]++];
        bfsWork++;
        int x = X[j], y = Y[j];
        for (int e = 0; e < 4; e++) {
          u32 o = occ(x + DX[e], y + DY[e]);
          if (!o) continue;
          int k = o - 1;
          if (mark[k] != stamp) { mark[k] = stamp; lab[k] = h; q[h].push_back(k); }
          else if (lab[k] != 255) { int a = find(lab[k]), b = find(h); if (a != b) par[a] = b; }
        }
      }
      int rB = find(gd), rA = ga >= 0 ? find(ga) : -1;
      if (rA == rB) return;                     // c' is in B
      size_t szB = 0, szEx = 0; bool liveB = false; int liveOther = 0;
      bool live[5] = {false, false, false, false, false};
      for (int h = 0; h < ngrp; h++) if (qh[h] < q[h].size()) live[find(h)] = true;
      for (int h = 0; h < ngrp; h++) {
        int r = find(h);
        if (r == rB) szB += q[h].size(); else if (!live[r]) szEx += q[h].size();
        if (r == h) { if (r == rB) liveB = live[r]; else if (live[r]) liveOther++; }
      }
      if ((int)szB > smax) return;
      if (liveB && liveOther == 0) {            // B is "everything else"
        if ((size_t)(n - 1) - szEx > (size_t)smax) return;
      }
      // images of the B cells discovered so far
      for (int h = 0; h < ngrp; h++) {
        if (find(h) != rB) continue;
        while (imgi[h] < q[h].size()) {
          int j = q[h][imgi[h]++];
          int dx = X[j] - cx, dy = Y[j] - cy;
          int x = ax + M[0] * dx + M[1] * dy, y = ay + M[2] * dx + M[3] * dy;
          for (int e = -1; e < 4; e++) {
            u32 o = e < 0 ? occ(x, y) : occ(x + DX[e], y + DY[e]);
            if (!o) continue;
            int k = o - 1;
            if (k == ia) { if (e < 0) return; continue; }
            if (mark[k] == stamp) {
              if (lab[k] == 255) return;         // hits or touches c
              int r = find(lab[k]);
              if (r == rB) continue;
              if (r == rA) return;               // certainly in A
            }
            if (pend[k] != stamp) { pend[k] = stamp; pendList.push_back(k); }
          }
        }
      }
      if (!liveB) done = true;
    }
    int rB = find(gd);
    for (int k : pendList) if (!(mark[k] == stamp && lab[k] != 255 && find(lab[k]) == rB)) return;
    if (ia == ic && g == 0) return;             // identity
    B.clear();
    for (int h = 0; h < ngrp; h++) if (find(h) == rB) B.insert(B.end(), q[h].begin(), q[h].end());
    size_t bsz = B.size();
    int hb = 0; { size_t t = bsz; while (t > 1) { t >>= 1; hb++; } }
    nx.resize(bsz); ny.resize(bsz);
    for (size_t i = 0; i < bsz; i++) {
      int j = B[i]; int dx = X[j] - cx, dy = Y[j] - cy;
      nx[i] = ax + M[0] * dx + M[1] * dy; ny[i] = ay + M[2] * dx + M[3] * dy;
    }
    auto apply = [&]() {
      for (int j : B) clearCell(X[j], Y[j]);
      for (size_t i = 0; i < bsz; i++) { int j = B[i]; std::swap(X[j], nx[i]); std::swap(Y[j], ny[i]); }
      for (int j : B) setCell(X[j], Y[j], j + 1);
    };
    if (logq != 0.0) {
      int before = localPerim(B, nx, ny);
      apply();
      int dt = localPerim(B, nx, ny) - before;
      if (dt != 0 && rng.uni() >= std::exp(logq * dt)) { apply(); return; }
    } else apply();
    accCP++; cpMoved += bsz; cpAccHist[hb]++;
  }

  // number of perimeter sites among: empty neighbours of current B cells, the sites
  // (ox,oy) and their empty neighbours.  Each site counted once via a small sort.
  std::vector<u64> tmpSites;
  int localPerim(const std::vector<int> &Bc, const std::vector<int32_t> &ox, const std::vector<int32_t> &oy) {
    tmpSites.clear();
    auto push = [&](int x, int y) {
      if (occ(x, y)) return;
      bool p = false;
      for (int e = 0; e < 4; e++) if (occ(x + DX[e], y + DY[e])) { p = true; break; }
      if (p) tmpSites.push_back(((u64)(u32)x << 32) | (u32)y);
    };
    for (size_t i = 0; i < Bc.size(); i++) {
      int j = Bc[i];
      push(X[j], Y[j]);
      for (int e = 0; e < 4; e++) push(X[j] + DX[e], Y[j] + DY[e]);
      push(ox[i], oy[i]);
      for (int e = 0; e < 4; e++) push(ox[i] + DX[e], oy[i] + DY[e]);
    }
    std::sort(tmpSites.begin(), tmpSites.end());
    return (int)(std::unique(tmpSites.begin(), tmpSites.end()) - tmpSites.begin());
  }

  // --- observables -------------------------------------------------------------------
  struct Obs { double rg2, asph; long perim, edges, cycles, noncut, leaves, w, h, holes, deg[5]; };
  std::vector<int> disc, low, stk, pit;
  Obs measure(bool full) {
    Obs o{};
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    int x0 = X[0], y0 = Y[0], xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    for (int i = 0; i < n; i++) {
      double x = X[i] - x0, y = Y[i] - y0;
      sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
      xmin = std::min(xmin, X[i] - x0); xmax = std::max(xmax, X[i] - x0);
      ymin = std::min(ymin, Y[i] - y0); ymax = std::max(ymax, Y[i] - y0);
    }
    o.w = xmax - xmin + 1; o.h = ymax - ymin + 1;
    double mx = sx / n, my = sy / n;
    double cxx = sxx / n - mx * mx, cyy = syy / n - my * my, cxy = sxy / n - mx * my;
    o.rg2 = cxx + cyy;
    double disc2 = std::sqrt((cxx - cyy) * (cxx - cyy) + 4 * cxy * cxy);
    double l1 = 0.5 * (cxx + cyy + disc2), l2 = 0.5 * (cxx + cyy - disc2);
    o.asph = (l1 + l2) > 0 ? (l1 - l2) * (l1 - l2) / ((l1 + l2) * (l1 + l2)) : 0;
    if (!full) return o;
    long perim = 0, edges = 0, leaves = 0, verts = 0;
    for (int i = 0; i < n; i++) {
      int deg = 0;
      for (int d = 0; d < 4; d++) {
        int ax = X[i] + DX[d], ay = Y[i] + DY[d];
        if (occ(ax, ay)) { deg++; continue; }
        // count perimeter site once: from its lexicographically first occupied neighbour dir
        bool first = true;
        for (int e = 0; e < 4; e++) {
          u32 b = occ(ax + DX[e], ay + DY[e]);
          if (b && (int)b - 1 < i) { first = false; break; }
        }
        if (first) perim++;
      }
      edges += deg; if (deg == 1) leaves++;
      o.deg[deg]++;
      // lattice vertices: count each corner once, from the lowest-index cell touching it
      for (int cy2 = 0; cy2 < 2; cy2++) for (int cx2 = 0; cx2 < 2; cx2++) {
        int vx = X[i] + cx2, vy = Y[i] + cy2; bool first = true;
        for (int a = -1; a <= 0 && first; a++) for (int b = -1; b <= 0; b++) {
          u32 c = occ(vx + a, vy + b);
          if (c && (int)c - 1 < i) { first = false; break; }
        }
        if (first) verts++;
      }
    }
    // Euler characteristic of the union of closed squares: V - E + F = 1 - holes
    o.holes = 1 - (verts - (4L * n - edges / 2) + n);
    o.perim = perim; o.edges = edges / 2; o.cycles = o.edges - (n - 1); o.leaves = leaves;
    // articulation points (iterative Tarjan)
    disc.assign(n, 0); low.assign(n, 0); stk.clear(); pit.assign(n, 0);
    std::vector<int> parent(n, -1), childc(n, 0);
    std::vector<uint8_t> isart(n, 0);
    int timer = 0;
    stk.push_back(0); disc[0] = low[0] = ++timer;
    while (!stk.empty()) {
      int v = stk.back();
      if (pit[v] < 4) {
        int d = pit[v]++;
        u32 a = occ(X[v] + DX[d], Y[v] + DY[d]);
        if (!a) continue;
        int w = a - 1;
        if (!disc[w]) { parent[w] = v; childc[v]++; disc[w] = low[w] = ++timer; stk.push_back(w); }
        else if (w != parent[v]) low[v] = std::min(low[v], disc[w]);
      } else {
        stk.pop_back();
        int p = parent[v];
        if (p >= 0) {
          low[p] = std::min(low[p], low[v]);
          if (parent[p] >= 0 && low[v] >= disc[p]) isart[p] = 1;
        }
      }
    }
    if (childc[0] >= 2) isart[0] = 1;
    long art = 0; for (int i = 0; i < n; i++) art += isart[i];
    o.noncut = n - art;
    return o;
  }

  bool selfCheck() {
    // grid/list consistency + connectivity
    for (int i = 0; i < n; i++) if (occ(X[i], Y[i]) != (u32)i + 1) return false;
    std::vector<uint8_t> seen(n, 0); std::vector<int> st{0}; seen[0] = 1; int cnt = 1;
    while (!st.empty()) {
      int v = st.back(); st.pop_back();
      for (int d = 0; d < 4; d++) { u32 a = occ(X[v] + DX[d], Y[v] + DY[d]); if (a && !seen[a - 1]) { seen[a - 1] = 1; cnt++; st.push_back(a - 1); } }
    }
    return cnt == n;
  }
};

static double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

int main(int argc, char **argv) {
  int n = 100, lgW = 0; u64 seed = 1; double pcp = 0.5, tsec = 0; u64 moves = 1000000, every = 0, burn = 0;
  std::string seedShape = "rect", out = "", mode = "run", load = "", snap = ""; double p = 0; u64 snapEvery = 0; int full = 1;
  int smax = -1, relK = 1, refCP = 0, inflate = 1; double ppv = 0, ckptSec = 120; std::string ckpt = "";
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    auto nxt = [&]() { return std::string(argv[++i]); };
    if (a == "-n") n = atoi(nxt().c_str());
    else if (a == "-seed") seed = strtoull(nxt().c_str(), 0, 10);
    else if (a == "-pcp") pcp = atof(nxt().c_str());
    else if (a == "-moves") moves = (u64)atof(nxt().c_str());
    else if (a == "-every") every = (u64)atof(nxt().c_str());
    else if (a == "-burn") burn = (u64)atof(nxt().c_str());
    else if (a == "-shape") seedShape = nxt();
    else if (a == "-out") out = nxt();
    else if (a == "-mode") mode = nxt();
    else if (a == "-lgW") lgW = atoi(nxt().c_str());
    else if (a == "-p") p = atof(nxt().c_str());
    else if (a == "-load") load = nxt();
    else if (a == "-snap") snap = nxt();
    else if (a == "-snapEvery") snapEvery = (u64)atof(nxt().c_str());
    else if (a == "-full") full = atoi(nxt().c_str());
    else if (a == "-smax") smax = atoi(nxt().c_str());
    else if (a == "-ckpt") ckpt = nxt();
    else if (a == "-ckptSec") ckptSec = atof(nxt().c_str());
    else if (a == "-inflate") inflate = atoi(nxt().c_str());
    else if (a == "-relK") relK = atoi(nxt().c_str());
    else if (a == "-ref") refCP = atoi(nxt().c_str());
    else if (a == "-ppv") ppv = atof(nxt().c_str());
    else if (a == "-tsec") tsec = atof(nxt().c_str());
    else { fprintf(stderr, "unknown arg %s\n", a.c_str()); return 1; }
  }
  std::vector<std::pair<int,int>> cells;
  if (!load.empty()) {
    FILE *f = fopen(load.c_str(), "rb"); if (!f) { perror("load"); return 1; }
    int32_t hdr; if (fread(&hdr, 4, 1, f) != 1) return 1; int m = hdr;
    std::vector<int32_t> buf(2 * (size_t)m); if (fread(buf.data(), 4, 2 * (size_t)m, f) != 2 * (size_t)m) return 1; fclose(f);
    // optional k-fold inflation (each cell becomes a k x k block), then pad to n cells by
    // adding random perimeter sites.  Only a starting point: the chain forgets it.
    for (int i = 0; i < m; i++)
      for (int a = 0; a < inflate; a++) for (int b = 0; b < inflate; b++)
        cells.push_back({buf[2 * i] * inflate + a, buf[2 * i + 1] * inflate + b});
    if ((int)cells.size() > n || n == 100) n = (int)cells.size();
    if ((int)cells.size() < n) {
      std::unordered_set<u64> occ;
      auto key = [](int x, int y) { return ((u64)(u32)x << 32) | (u32)y; };
      for (auto &c : cells) occ.insert(key(c.first, c.second));
      Rng r(seed ^ 0xabcdef);
      while ((int)cells.size() < n) {
        auto c = cells[r.below((u32)cells.size())]; int d = r.below(4);
        int x = c.first + DX[d], y = c.second + DY[d];
        if (occ.insert(key(x, y)).second) cells.push_back({x, y});
      }
    }
  } else if (seedShape == "bar") {
    for (int i = 0; i < n; i++) cells.push_back({i, 0});
  } else {
    int w = (int)std::ceil(std::sqrt((double)n));
    for (int i = 0; i < n; i++) cells.push_back({i % w, i / w});
  }
  if (!lgW) {
    // generous: 24 * n^0.64 (and at least the seed extent + margin)
    double ext = 10.0 * std::pow((double)n, 0.64);
    if (seedShape == "bar" && load.empty()) ext = std::max(ext, n + 8.0);
    int xmin = 1 << 30, xmax = -(1 << 30), ymin = 1 << 30, ymax = -(1 << 30);
    for (auto &c : cells) { xmin = std::min(xmin, c.first); xmax = std::max(xmax, c.first); ymin = std::min(ymin, c.second); ymax = std::max(ymax, c.second); }
    ext = std::max(ext, 2.5 * std::max(xmax - xmin, ymax - ymin) + 8.0);
    lgW = 7; while ((1 << lgW) < ext) lgW++;
  }
  Sampler S(n, lgW, seed);
  S.place(cells);
  if (p > 0) S.logq = std::log(1.0 - p);
  if (smax >= 0) S.smax = std::min(smax, (n - 1) / 2);
  S.relK = relK; S.refCP = refCP;
  u64 thr = (u64)(pcp * 18446744073709551615.0), thr2 = (u64)((pcp + ppv) * 18446744073709551615.0);
  auto step = [&]() { u64 r = S.rng.next(); if (r < thr) S.cutpaste(false); else if (r < thr2) S.cutpaste(true); else S.relocate(); };

  if (mode == "exact") {
    if (!every) every = 1;
    // histogram of translation-normalised shapes; chi-square against uniform
    std::map<std::string, u64> hist;
    u64 samples = 0;
    for (u64 m = 0; m < moves; m++) {
      step();
      if (m >= burn && (m % every) == 0) {
        int xmin = 1 << 30, ymin = 1 << 30;
        for (int i = 0; i < n; i++) { xmin = std::min(xmin, S.X[i]); ymin = std::min(ymin, S.Y[i]); }
        std::vector<std::pair<int,int>> c;
        for (int i = 0; i < n; i++) c.push_back({S.X[i] - xmin, S.Y[i] - ymin});
        std::sort(c.begin(), c.end());
        std::string key; for (auto &pp : c) { key.push_back((char)pp.first); key.push_back((char)pp.second); }
        hist[key]++; samples++;
      }
    }
    double expct = 0, chi2 = 0;
    if (p == 0) {
      expct = (double)samples / hist.size();
      for (auto &kv : hist) chi2 += (kv.second - expct) * (kv.second - expct) / expct;
      printf("n=%d shapes=%zu samples=%" PRIu64 " chi2=%.2f dof=%zu z=%.2f\n", n, hist.size(), samples, chi2,
             hist.size() - 1, (chi2 - (hist.size() - 1)) / std::sqrt(2.0 * (hist.size() - 1)));
    } else {
      // tilted: expected weight q^t
      double Z = 0; std::vector<double> w;
      for (auto &kv : hist) {
        // compute perimeter of the shape
        std::vector<std::pair<int,int>> c;
        for (size_t i = 0; i < kv.first.size(); i += 2) c.push_back({kv.first[i], kv.first[i + 1]});
        std::map<std::pair<int,int>, int> oc; for (auto &pp : c) oc[pp] = 1;
        std::map<std::pair<int,int>, int> per;
        for (auto &pp : c) for (int d = 0; d < 4; d++) { std::pair<int,int> s{pp.first + DX[d], pp.second + DY[d]}; if (!oc.count(s)) per[s] = 1; }
        double ww = std::exp(S.logq * per.size()); w.push_back(ww); Z += ww;
      }
      size_t i = 0;
      for (auto &kv : hist) { double e = samples * w[i++] / Z; chi2 += (kv.second - e) * (kv.second - e) / e; }
      printf("n=%d p=%.3f shapes=%zu samples=%" PRIu64 " chi2=%.2f dof=%zu z=%.2f\n", n, p, hist.size(), samples, chi2,
             hist.size() - 1, (chi2 - (hist.size() - 1)) / std::sqrt(2.0 * (hist.size() - 1)));
    }
    printf("acc relocate %.4f  cutpaste %.4f\n", (double)S.accRel / std::max<u64>(1, S.nRel), (double)S.accCP / std::max<u64>(1, S.nCP));
    return 0;
  }

  // run mode: time series
  FILE *fo = out.empty() ? stdout : fopen(out.c_str(), "w");
  if (!every) every = std::max<u64>(1, moves / 1000);
  fprintf(fo, "# n=%d seed=%" PRIu64 " pcp=%.3f shape=%s lgW=%d p=%.4f every=%" PRIu64 "\n", n, seed, pcp, seedShape.c_str(), lgW, p, every);
  fprintf(fo, "# move sec rg2 asph perim cycles noncut leaves w h holes deg1 deg2 deg3 deg4\n");
  double t0 = now(), lastCkpt = t0;
  u64 m = 0; int snapIdx = 0;
  for (; m < moves; m++) {
    step();
    if ((m + 1) % every == 0) {
      auto o = S.measure(full && (m + 1) >= burn);
      if (std::max(o.w, o.h) * 2 > S.W) { fprintf(stderr, "extent %ld exceeds half the torus %d — aborting\n", std::max(o.w, o.h), S.W); return 2; }
      fprintf(fo, "%" PRIu64 " %.2f %.6g %.5f %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld\n", m + 1, now() - t0, o.rg2, o.asph, o.perim, o.cycles, o.noncut, o.leaves, o.w, o.h, o.holes, o.deg[1], o.deg[2], o.deg[3], o.deg[4]);
      fflush(fo);
      if (!ckpt.empty() && now() - lastCkpt > ckptSec) {   // rolling checkpoint, replaced atomically
        std::string tmp = ckpt + ".tmp";
        FILE *f = fopen(tmp.c_str(), "wb"); int32_t hdr = n; fwrite(&hdr, 4, 1, f);
        for (int i = 0; i < n; i++) { int32_t xy[2] = {S.X[i], S.Y[i]}; fwrite(xy, 4, 2, f); }
        fclose(f); rename(tmp.c_str(), ckpt.c_str()); lastCkpt = now();
      }
      if (tsec > 0 && now() - t0 > tsec) { m++; break; }
    }
    if (snapEvery && !snap.empty() && (m + 1) % snapEvery == 0) {
      char fn[512]; snprintf(fn, sizeof fn, "%s_%04d.bin", snap.c_str(), snapIdx++);
      FILE *f = fopen(fn, "wb"); int32_t hdr = n; fwrite(&hdr, 4, 1, f);
      for (int i = 0; i < n; i++) { int32_t xy[2] = {S.X[i], S.Y[i]}; fwrite(xy, 4, 2, f); }
      fclose(f);
    }
  }
  double el = now() - t0;
  if (!S.selfCheck()) { fprintf(stderr, "SELF CHECK FAILED\n"); return 3; }
  if (!snap.empty() && !snapEvery) {
    FILE *f = fopen(snap.c_str(), "wb"); int32_t hdr = n; fwrite(&hdr, 4, 1, f);
    for (int i = 0; i < n; i++) { int32_t xy[2] = {S.X[i], S.Y[i]}; fwrite(xy, 4, 2, f); }
    fclose(f);
  }
  fprintf(stderr, "n=%d moves=%" PRIu64 " %.1fs (%.2f Mmoves/s) accRel=%.4f accCP=%.5f bfs/move=%.1f movedCells/CP=%.2f\n",
          n, m, el, m / el / 1e6, (double)S.accRel / std::max<u64>(1, S.nRel), (double)S.accCP / std::max<u64>(1, S.nCP),
          (double)S.bfsWork / m, (double)S.cpMoved / std::max<u64>(1, S.nCP));
  fprintf(stderr, "cp branch-size histogram (log2 bins: proposed/accepted):");
  for (int i = 0; i < 40; i++) if (S.cpHist[i]) fprintf(stderr, " [%d]%" PRIu64 "/%" PRIu64, i, S.cpHist[i], S.cpAccHist[i]);
  fprintf(stderr, "\n");
  if (fo != stdout) fclose(fo);
  return 0;
}
