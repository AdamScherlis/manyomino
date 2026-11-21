/* outward: the exact step of the zoom out, small enough to run in a browser.
 *
 * A world is a whole polyomino of n cells that stands for the shape around a root at one scale.
 * One step outward:
 *   coarsen()      every 2x2 block holding a cell becomes one cell;
 *   relax(T)       the local walk of src/zoom.cpp (single cells added and removed at weight
 *                  1/lambda each, connectivity certified in a box of random radius) for T sweeps;
 *   pad(N)         cells added at random empty sites next to the shape, outside the held box,
 *                  up to N;
 *   run(M)         M moves of the sampler of src/manyomino.cpp (relocate with a random-radius
 *                  certificate, and branch transplants) with the box held: nothing in the box
 *                  moves and nothing is put there, so the chain samples uniformly from the
 *                  polyominoes of N cells that agree with the world inside the box.
 *
 * Free-standing C (no library), for clang --target=wasm32; `make_outward.sh` builds it.  With
 * -DNATIVE it also builds as an ordinary program for testing.  The lattice is a torus of side
 * W = 2^LG, far wider than a world of the sizes used.
 */
typedef unsigned int u32;
typedef unsigned long long u64;
#ifndef LG
#define LG 11
#endif
#define W (1 << LG)
#define MASK (W - 1)
#ifndef NMAX
#define NMAX 49152
#endif

static int grid[W * W];                       /* index of the cell at a site, plus 1; 0 if empty */
static int X[NMAX], Y[NMAX], n;
static u32 mark[NMAX], pend[NMAX], stamp;
static unsigned char lab[NMAX];
static int Q[5][NMAX], qn[5], qh[5], par[5], imgi[5];
static int Bl[NMAX], nxb[NMAX], nyb[NMAX], pendList[NMAX];
static int frz, fx0, fy0, fx1, fy1;
static u64 rs[4];
static const int DX[4] = {1, 0, -1, 0}, DY[4] = {0, 1, 0, -1};
static const int G8[8][4] = {{1,0,0,1},{0,-1,1,0},{-1,0,0,-1},{0,1,-1,0},{1,0,0,-1},{0,1,1,0},{-1,0,0,1},{0,-1,-1,0}};
static u64 nRel, accRel, nCP, accCP, nGC, accGC;

static inline u64 rotl(u64 x, int k) { return (x << k) | (x >> (64 - k)); }
static inline u64 rnext(void) {
  u64 r = rotl(rs[0] + rs[3], 23) + rs[0], t = rs[1] << 17;
  rs[2] ^= rs[0]; rs[3] ^= rs[1]; rs[1] ^= rs[2]; rs[0] ^= rs[3]; rs[2] ^= t; rs[3] = rotl(rs[3], 45);
  return r;
}
static inline double runi(void) { return (double)(rnext() >> 11) * (1.0 / 9007199254740992.0); }
static inline u32 rbelow(u32 m) { return (u32)(((rnext() >> 32) * (u64)m) >> 32); }
static inline int at(int x, int y) { return ((y & MASK) << LG) | (x & MASK); }
static inline int occ(int x, int y) { return grid[at(x, y)]; }
static inline int held(int x, int y) { return x >= fx0 && x < fx1 && y >= fy0 && y < fy1; }
static inline int find(int a) { while (par[a] != a) a = par[a] = par[par[a]]; return a; }
static inline int iabs(int a) { return a < 0 ? -a : a; }

void seed(u32 a, u32 b) {
  u64 x = ((u64)a << 32) | b;
  for (int i = 0; i < 4; i++) { u64 z = (x += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; rs[i] = z ^ (z >> 31); }
}
int *xs(void) { return X; }
int *ys(void) { return Y; }
int count(void) { return n; }
int capacity(void) { return NMAX; }
int side(void) { return W; }
/* take the m cells written at xs(), ys() as the world */
void place(int m) {
  for (int i = 0; i < W * W; i++) grid[i] = 0;
  n = m; for (int i = 0; i < n; i++) grid[at(X[i], Y[i])] = i + 1;
}
void hold(int x0, int y0, int x1, int y1, int on) { frz = on; fx0 = x0; fy0 = y0; fx1 = x1; fy1 = y1; }

/* Are the neighbours (leaving out cell `ex`) of site (cx, cy) joined to each other inside the box of
   radius K about it?  The site itself counts as empty. */
static int connAt(int cx, int cy, int ex, int K) {
  static const int RX[8] = {1,1,0,-1,-1,-1,0,1}, RY[8] = {0,1,1,1,0,-1,-1,-1};
  int m = 0;
  for (int i = 0; i < 8; i++) { int o = occ(cx + RX[i], cy + RY[i]); if (o && o - 1 != ex) m |= 1 << i; }
  int deg = ((m >> 0) & 1) + ((m >> 2) & 1) + ((m >> 4) & 1) + ((m >> 6) & 1);
  if (deg <= 1) return 1;
  int links = 0;
  for (int i = 0; i < 8; i += 2) if ((m >> i & 1) && (m >> ((i + 1) & 7) & 1) && (m >> ((i + 2) & 7) & 1)) links++;
  if (links >= deg - 1) return 1;
  if (K <= 1) return 0;
  stamp++;
  if (ex >= 0) { mark[ex] = stamp; lab[ex] = 255; }
  int ng = 0;
  for (int d = 0; d < 4; d++) {
    int o = occ(cx + DX[d], cy + DY[d]);
    if (!o || o - 1 == ex) continue;
    int j = o - 1, g = ng++;
    qn[g] = 0; qh[g] = 0; Q[g][qn[g]++] = j; mark[j] = stamp; lab[j] = (unsigned char)g; par[g] = g;
  }
  int ncls = ng;
  for (;;) {
    for (int g = 0; g < ng; g++) {
      if (qh[g] >= qn[g]) continue;
      int j = Q[g][qh[g]++], x = X[j], y = Y[j];
      for (int d = 0; d < 4; d++) {
        int ax = x + DX[d], ay = y + DY[d];
        if (iabs(ax - cx) > K || iabs(ay - cy) > K) continue;
        if (ax == cx && ay == cy) continue;
        int o = occ(ax, ay);
        if (!o) continue;
        int k = o - 1;
        if (mark[k] != stamp) { mark[k] = stamp; lab[k] = (unsigned char)g; Q[g][qn[g]++] = k; }
        else if (lab[k] != 255) { int a = find(lab[k]), b = find(g); if (a != b) { par[a] = b; if (--ncls == 1) return 1; } }
      }
      if (qh[g] >= qn[g]) {
        int r = find(g), dry = 1;
        for (int h = 0; h < ng; h++) if (find(h) == r && qh[h] < qn[h]) { dry = 0; break; }
        if (dry) return 0;
      }
    }
  }
}
static inline int drawK(int cap) { double u = runi(); if (u * cap < 1.0) return cap; return (int)(1.0 / u); }

/* relocate: one cell to an empty site next to the shape */
static void relocate(void) {
  nRel++;
  int ic = rbelow(n), iu = rbelow(n), d = rbelow(4);
  int K = drawK(1 << 28);
  if (iu == ic) return;
  int sx = X[iu] + DX[d], sy = Y[iu] + DY[d], cx = X[ic], cy = Y[ic];
  if (occ(sx, sy)) return;
  if (frz && (held(cx, cy) || held(sx, sy))) return;
  int k = 0;
  for (int e = 0; e < 4; e++) { int a = occ(sx + DX[e], sy + DY[e]); if (a && a - 1 != ic) k++; }
  if (k > 1 && rbelow(k) != 0) return;
  if (!connAt(cx, cy, ic, K)) return;
  if (k > 1 && !connAt(sx, sy, ic, K)) return;
  grid[at(cx, cy)] = 0; X[ic] = sx; Y[ic] = sy; grid[at(sx, sy)] = ic + 1;
  accRel++;
}

static inline int locallyConnected(int cx, int cy) {
  static const int RX[8] = {1,1,0,-1,-1,-1,0,1}, RY[8] = {0,1,1,1,0,-1,-1,-1};
  int m = 0;
  for (int i = 0; i < 8; i++) if (occ(cx + RX[i], cy + RY[i])) m |= 1 << i;
  int deg = ((m >> 0) & 1) + ((m >> 2) & 1) + ((m >> 4) & 1) + ((m >> 6) & 1);
  if (deg <= 1) return 1;
  int links = 0;
  for (int i = 0; i < 8; i += 2) if ((m >> i & 1) && (m >> ((i + 1) & 7) & 1) && (m >> ((i + 2) & 7) & 1)) links++;
  return links >= deg - 1;
}

/* transplant: the branch B hanging off cell c in direction d is moved to hang off cell c' instead,
   turned by one of the eight symmetries.  The searches from the neighbours of c and from c' run in
   step, and the image of each branch cell is tested as it is found. */
static void transplant(void) {
  nCP++;
  int ic = rbelow(n), d = rbelow(4), ia = rbelow(n), g = rbelow(8), smax = (n - 1) / 2;
  int cx = X[ic], cy = Y[ic];
  if (!occ(cx + DX[d], cy + DY[d])) return;
  if (locallyConnected(cx, cy)) return;
  stamp++;
  mark[ic] = stamp; lab[ic] = 255;
  int ngrp = 0, gd = -1;
  for (int e = 0; e < 4; e++) {
    int o = occ(cx + DX[e], cy + DY[e]);
    if (!o) continue;
    int j = o - 1, h = ngrp++;
    qn[h] = 0; qh[h] = 0; imgi[h] = 0; Q[h][qn[h]++] = j; mark[j] = stamp; lab[j] = (unsigned char)h; par[h] = h;
    if (e == d) gd = h;
  }
  int ga = -1;
  if (ia != ic) {
    if (mark[ia] == stamp) ga = lab[ia];
    else { ga = ngrp++; qn[ga] = 0; qh[ga] = 0; imgi[ga] = 0; Q[ga][qn[ga]++] = ia; mark[ia] = stamp; lab[ia] = (unsigned char)ga; par[ga] = ga; }
    if (ga == gd) return;
  }
  const int ax = X[ia], ay = Y[ia]; const int *M = G8[g];
  int npend = 0, done = 0;
  while (!done) {
    for (int h = 0; h < ngrp; h++) {
      if (qh[h] >= qn[h]) continue;
      int j = Q[h][qh[h]++], x = X[j], y = Y[j];
      for (int e = 0; e < 4; e++) {
        int o = occ(x + DX[e], y + DY[e]);
        if (!o) continue;
        int k = o - 1;
        if (mark[k] != stamp) { mark[k] = stamp; lab[k] = (unsigned char)h; Q[h][qn[h]++] = k; }
        else if (lab[k] != 255) { int a = find(lab[k]), b = find(h); if (a != b) par[a] = b; }
      }
    }
    int rB = find(gd), rA = ga >= 0 ? find(ga) : -1;
    if (rA == rB) return;
    int szB = 0, szEx = 0, liveB = 0, liveOther = 0, live[5] = {0, 0, 0, 0, 0};
    for (int h = 0; h < ngrp; h++) if (qh[h] < qn[h]) live[find(h)] = 1;
    for (int h = 0; h < ngrp; h++) {
      int r = find(h);
      if (r == rB) szB += qn[h]; else if (!live[r]) szEx += qn[h];
      if (r == h) { if (r == rB) liveB = live[r]; else if (live[r]) liveOther++; }
    }
    if (szB > smax) return;
    if (liveB && liveOther == 0 && (n - 1) - szEx > smax) return;
    for (int h = 0; h < ngrp; h++) {
      if (find(h) != rB) continue;
      while (imgi[h] < qn[h]) {
        int j = Q[h][imgi[h]++], dx = X[j] - cx, dy = Y[j] - cy;
        int x = ax + M[0] * dx + M[1] * dy, y = ay + M[2] * dx + M[3] * dy;
        for (int e = -1; e < 4; e++) {
          int o = e < 0 ? occ(x, y) : occ(x + DX[e], y + DY[e]);
          if (!o) continue;
          int k = o - 1;
          if (k == ia) { if (e < 0) return; continue; }
          if (mark[k] == stamp) {
            if (lab[k] == 255) return;
            int r = find(lab[k]);
            if (r == rB) continue;
            if (r == rA) return;
          }
          if (pend[k] != stamp) { pend[k] = stamp; pendList[npend++] = k; }
        }
      }
    }
    if (!liveB) done = 1;
  }
  int rB = find(gd);
  for (int i = 0; i < npend; i++) { int k = pendList[i]; if (!(mark[k] == stamp && lab[k] != 255 && find(lab[k]) == rB)) return; }
  if (ia == ic && g == 0) return;
  int bsz = 0;
  for (int h = 0; h < ngrp; h++) if (find(h) == rB) for (int i = 0; i < qn[h]; i++) Bl[bsz++] = Q[h][i];
  for (int i = 0; i < bsz; i++) {
    int j = Bl[i], dx = X[j] - cx, dy = Y[j] - cy;
    nxb[i] = ax + M[0] * dx + M[1] * dy; nyb[i] = ay + M[2] * dx + M[3] * dy;
    if (frz && (held(X[j], Y[j]) || held(nxb[i], nyb[i]))) return;
  }
  for (int i = 0; i < bsz; i++) grid[at(X[Bl[i]], Y[Bl[i]])] = 0;
  for (int i = 0; i < bsz; i++) { int j = Bl[i]; X[j] = nxb[i]; Y[j] = nyb[i]; grid[at(X[j], Y[j])] = j + 1; }
  accCP++;
}

/* m moves of the sampler: transplants 40% of the time */
void run(int m) {
  const u64 thr = (u64)(0.4 * 18446744073709551615.0);
  for (int i = 0; i < m; i++) { if (rnext() < thr) transplant(); else relocate(); }
}

/* half the resolution: a site is occupied if any of the four it stands for is */
int coarsen(void) {
  int m = n;
  for (int i = 0; i < m; i++) grid[at(X[i], Y[i])] = 0;
  n = 0;
  for (int i = 0; i < m; i++) {
    int x = X[i] >> 1, y = Y[i] >> 1;
    if (occ(x, y)) continue;
    X[n] = x; Y[n] = y; grid[at(x, y)] = n + 1; n++;
  }
  return n;
}

/* one move of the local walk: a cell added next to the shape or taken away, the count free */
static void gcMove(void) {
  const double Z = 1.0 / 4.0625696;
  nGC++;
  int i = rbelow(n), r = (int)(rnext() >> 61), cx = X[i], cy = Y[i];
  if (r < 4) {
    int sx = cx + DX[r], sy = cy + DY[r];
    if (occ(sx, sy) || n >= NMAX) return;
    int k = 0; for (int e = 0; e < 4; e++) if (occ(sx + DX[e], sy + DY[e])) k++;
    if (k > 1 && rbelow(k) != 0) return;
    double a = 4.0 * Z * n / (n + 1.0);
    if (a < 1.0 && runi() >= a) return;
    if (!connAt(sx, sy, -1, drawK(256))) return;
    X[n] = sx; Y[n] = sy; grid[at(sx, sy)] = n + 1; n++; accGC++;
  } else {
    if (n <= 1) return;
    int k = 0; for (int e = 0; e < 4; e++) if (occ(cx + DX[e], cy + DY[e])) k++;
    if (!k) return;
    double a = n / (4.0 * Z * (n - 1.0));
    if (a < 1.0 && runi() >= a) return;
    if (!connAt(cx, cy, i, drawK(256))) return;
    grid[at(cx, cy)] = 0; n--;
    if (i != n) { X[i] = X[n]; Y[i] = Y[n]; grid[at(X[i], Y[i])] = i + 1; }
    accGC++;
  }
}
/* the local walk for `quarters` quarter-sweeps (a sweep is one proposal per cell) */
void relax(int quarters) {
  for (int q = 0; q < quarters; q++) { int m = n / 4 > 1 ? n / 4 : 1; for (int k = 0; k < m; k++) gcMove(); }
}

/* cells added at random empty sites next to the shape, outside the held box, until there are `target` */
int pad(int target) {
  if (target > NMAX) target = NMAX;
  u64 tries = 0;
  while (n < target && tries++ < 400ULL * NMAX) {
    int i = rbelow(n), d = rbelow(4), sx = X[i] + DX[d], sy = Y[i] + DY[d];
    if (occ(sx, sy) || (frz && held(sx, sy))) continue;
    X[n] = sx; Y[n] = sy; grid[at(sx, sy)] = n + 1; n++;
  }
  return n;
}
/* acceptance counts, for checking against src/manyomino.cpp: 0 relocations tried, 1 accepted, 2 transplants tried, 3 accepted, 4 local tried, 5 accepted */
double stat(int which) { u64 v[6] = {nRel, accRel, nCP, accCP, nGC, accGC}; return (double)v[which]; }

#ifdef NATIVE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* outward_native exact NCELLS MOVES SEED: every shape of NCELLS cells should turn up equally often
   outward_native held NCELLS MOVES SEED: the same with a 2x2 box held
   outward_native step IN OUT N S TC B SEED: one step outward from the shape IN (root at the origin) */
int main(int argc, char **argv) {
  if (argc > 4 && !strcmp(argv[1], "exact")) {
    int m = atoi(argv[2]); long long moves = atoll(argv[3]); seed(1, atoi(argv[4]));
    for (int i = 0; i < m; i++) { X[i] = i; Y[i] = 0; } place(m);
    /* count shapes by a hash of their cells relative to the lower-left corner */
    enum { H = 1 << 20 }; static unsigned long long keyv[H]; static unsigned cnt[H]; unsigned long long samples = 0; int shapes = 0;
    for (long long t = 0; t < moves; t++) {
      run(1);
      if (t % (10 * m)) continue;
      int x0 = 1 << 30, y0 = 1 << 30; for (int i = 0; i < n; i++) { if (X[i] < x0) x0 = X[i]; if (Y[i] < y0) y0 = Y[i]; }
      unsigned long long h = 0; for (int i = 0; i < n; i++) { unsigned long long c = (unsigned long long)(X[i] - x0) * 64 + (Y[i] - y0); c = (c + 1) * 0x9e3779b97f4a7c15ULL; c ^= c >> 29; c *= 0xbf58476d1ce4e5b9ULL; h += c ^ (c >> 32); }
      unsigned s = (unsigned)(h % H); while (cnt[s] && keyv[s] != h) s = (s + 1) % H;
      if (!cnt[s]) { keyv[s] = h; shapes++; } cnt[s]++; samples++;
    }
    double e = (double)samples / shapes, chi2 = 0; for (int s = 0; s < H; s++) if (cnt[s]) chi2 += (cnt[s] - e) * (cnt[s] - e) / e;
    printf("n=%d shapes=%d samples=%llu chi2=%.1f dof=%d; relocations accepted %.4f, transplants %.4f\n", m, shapes, samples, chi2, shapes - 1, stat(1) / stat(0), stat(3) / stat(2));
    return 0;
  }
  if (argc > 4 && !strcmp(argv[1], "held")) {
    /* the box [0,2)x[0,2) holds cells at (0,0) and (1,1) and two empty sites; the shapes of NCELLS cells
       that agree with it should turn up equally often (count them with analysis/held_count.py) */
    int m = atoi(argv[2]); long long moves = atoll(argv[3]); seed(2, atoi(argv[4]));
    static const int px[7] = {0, -1, -1, -1, 0, 1, 1}, py[7] = {0, 0, 1, 2, 2, 2, 1};
    for (int i = 0; i < 7; i++) { X[i] = px[i]; Y[i] = py[i]; } place(7); hold(0, 0, 2, 2, 1); pad(m);
    enum { H = 1 << 22 }; static unsigned long long keyv[H]; static unsigned cnt[H]; unsigned long long samples = 0; int shapes = 0;
    for (long long t = 0; t < moves; t++) {
      run(1);
      if (t % (10 * m)) continue;
      unsigned long long h = 0; for (int i = 0; i < n; i++) { unsigned long long c = (unsigned long long)(X[i] + 32) * 64 + (Y[i] + 32); c = (c + 1) * 0x9e3779b97f4a7c15ULL; c ^= c >> 29; c *= 0xbf58476d1ce4e5b9ULL; h += c ^ (c >> 32); }
      unsigned s = (unsigned)(h % H); while (cnt[s] && keyv[s] != h) s = (s + 1) % H;
      if (!cnt[s]) { keyv[s] = h; shapes++; } cnt[s]++; samples++;
    }
    double e = (double)samples / shapes, chi2 = 0; for (int s = 0; s < H; s++) if (cnt[s]) chi2 += (cnt[s] - e) * (cnt[s] - e) / e;
    printf("n=%d shapes=%d samples=%llu chi2=%.1f dof=%d; relocations accepted %.4f, transplants %.4f\n", m, shapes, samples, chi2, shapes - 1, stat(1) / stat(0), stat(3) / stat(2));
    return 0;
  }
  if (argc > 8 && !strcmp(argv[1], "step")) {
    FILE *f = fopen(argv[2], "rb"); int m; if (!f || fread(&m, 4, 1, f) != 1 || m > NMAX) return 1;
    static int buf[2 * NMAX]; if (fread(buf, 4, 2 * (size_t)m, f) != 2 * (size_t)m) return 1; fclose(f);
    for (int i = 0; i < m; i++) { X[i] = buf[2 * i]; Y[i] = buf[2 * i + 1]; }
    int N = atoi(argv[4]), S = atoi(argv[5]), TC = atoi(argv[6]), B = atoi(argv[7]); seed(7, atoi(argv[8]));
    place(m); int c = coarsen(); relax(4 * TC); int kept = n;
    hold(-B, -B, B, B, 1); pad(N);
    for (int s = 0; s < S; s++) run(n);
    fprintf(stderr, "outward: %d cells -> %d after halving -> %d after the walk -> %d padded; relocations accepted %.4f, transplants %.4f\n", m, c, kept, n, stat(1) / stat(0), stat(3) / stat(2));
    f = fopen(argv[3], "wb"); fwrite(&n, 4, 1, f); for (int i = 0; i < n; i++) { fwrite(&X[i], 4, 1, f); fwrite(&Y[i], 4, 1, f); } fclose(f);
    return 0;
  }
  fprintf(stderr, "usage: outward_native exact NCELLS MOVES SEED | held NCELLS MOVES SEED | step IN OUT N S TC B SEED\n"); return 1;
}
#endif
