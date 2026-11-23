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
 * The same local walk also draws the pictures.  A world's cells are large on the screen, so its
 * ring is refined several times (each cell becomes four, then the walk runs), every time with the
 * finer level's picture held fixed inside the box and joined to it (pin(), link()).
 *
 * Free-standing C (no library), for clang --target=wasm32; `scripts/make_outward.sh` builds it.
 * With -DNATIVE it also builds as an ordinary program for testing.  The lattice is a torus of
 * side 2^lg (init(lg)), far wider than a world or a picture of the sizes used.
 */
typedef unsigned int u32;
typedef unsigned long long u64;
#ifndef NMAX
#define NMAX 524288                            /* most cells in a picture */
#endif
#ifndef NW
#define NW 65536                               /* most cells in a world */
#endif

static int *grid, LG, W, MASK;                /* index of the cell at a site, plus 1; 0 if empty */
static int X[NMAX], Y[NMAX], n;
static u32 mark[NMAX], pend[NW], stamp;
static unsigned char lab[NMAX];
static int Q[5][NMAX], qn[5], qh[5], par[5], imgi[5];
static int Bl[NW], nxb[NW], nyb[NW], pendList[NW];
static int frz, fx0, fy0, fx1, fy1;
static unsigned char *tubeMap; static int tubeOn, coverOn; static double zmul = 1.0; /* the local walk after halving may be kept to where the finer shape was (bit 1) or next to it (bit 2) */
static int OX[NW], OY[NW], nO;
static int plab[NMAX], cutOn, labmax, lab0, boxOnly;
static int labcnt[NMAX + 2];                   /* for the walk on a picture: cells carrying each of the labels there were at the start */                   /* which part of the finer shape a held cell belongs to; two held cells of different parts that touch are not joined */
static int nfix, regOn, rin, rout;            /* the local walk: cells below nfix stay; cells are added only at rin <= r < rout */
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
static inline int cut(int j, int k) { return cutOn && plab[j] != plab[k] && plab[j] && plab[k] && (cutOn != 4 || (plab[j] <= lab0 && plab[k] <= lab0)); }

void seed(u32 a, u32 b) {
  u64 x = ((u64)a << 32) | b;
  for (int i = 0; i < 4; i++) { u64 z = (x += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; rs[i] = z ^ (z >> 31); }
}
int *xs(void) { return X; }
int *ys(void) { return Y; }
int count(void) { return n; }
int capacity(void) { return NMAX; }
int side(void) { return W; }
#ifdef NATIVE
void *calloc(unsigned long, unsigned long);
static void *room(unsigned long bytes) { return calloc(bytes, 1); }
#else
static void *room(unsigned long bytes) { unsigned long was = __builtin_wasm_memory_grow(0, (bytes + 65535) >> 16); return (void *)(was << 16); }
#endif
/* a lattice of side 2^lg; call once, first */
void init(int lg) { LG = lg; W = 1 << lg; MASK = W - 1; grid = (int *)room((unsigned long)W * W * sizeof(int)); tubeMap = (unsigned char *)room((unsigned long)W * W); n = 0; }
/* the present cells are forgotten: write new ones at xs(), ys() and call place() */
void clear(void) { for (int i = 0; i < n; i++) grid[at(X[i], Y[i])] = 0; n = 0; }
/* take the m cells written at xs(), ys() as the shape (after clear()) */
void place(int m) {
  n = m; for (int i = 0; i < n; i++) grid[at(X[i], Y[i])] = i + 1;
  nfix = 0; regOn = 0; frz = 0; tubeOn = 0; coverOn = 0; zmul = 1.0; cutOn = 0;
  for (int i = 0; i < NMAX; i++) plab[i] = 0;
}
/* Part labels, one per cell (write them after place(); 0 for a cell with none), and how they are heeded:
   1  two cells of the held box that touch are joined only if their labels are the same (the sampler; and the walk, in
      which a new cell in the box takes the label of what it touches there and may not touch two);
   5  as 1 for the walk that thins a halved world, and besides no cell is added in the box unless it touches a labelled
      one, and none is taken away if what it touched in the box would no longer be joined inside the box;
   4  for the walk on a picture: as 1 for the labels there were at the start (two pieces of the inside that only touch
      are not joined), no new cell in the box may touch two of those, the last cell with one of them is never taken
      away, and a new cell that touches none of them gets one new label, the same for all, which is joined to whatever
      it touches. */
int *partlabels(void) { return plab; }
void parts(int on) {
  cutOn = on; labmax = 0; if (on) for (int i = 0; i < n; i++) if (plab[i] > labmax) labmax = plab[i];
  lab0 = labmax < NMAX ? labmax : NMAX;
  if (on == 4) { for (int l = 0; l <= lab0 + 1; l++) labcnt[l] = 0; for (int i = 0; i < n; i++) if (plab[i] > 0 && plab[i] <= lab0) labcnt[plab[i]]++; }
}
/* distance from the root in the sense used throughout: the root's own site is 0, and so are the three sites up and right of it at -1 */
static inline int rad(int x, int y) { int a = x < 0 ? -x - 1 : x, b = y < 0 ? -y - 1 : y; return a > b ? a : b; }
/* the first k cells are fixed for the local walk, which may add cells only at a <= r < b */
void pin(int k, int a, int b) { nfix = k; regOn = 1; rin = a; rout = b; }
void hold(int x0, int y0, int x1, int y1, int on) { frz = on; fx0 = x0; fy0 = y0; fx1 = x1; fy1 = y1; }

/* Are the neighbours (leaving out cell `ex`) of site (cx, cy) joined to each other inside the box of
   radius K about it?  The site itself counts as empty. */
static int connAt(int cx, int cy, int ex, int K) {
  static const int RX[8] = {1,1,0,-1,-1,-1,0,1}, RY[8] = {0,1,1,1,0,-1,-1,-1};
  int m = 0, lb = boxOnly;
  for (int i = 0; i < 8; i++) { int o = occ(cx + RX[i], cy + RY[i]); if (o && o - 1 != ex && !(boxOnly && !held(cx + RX[i], cy + RY[i]))) { m |= 1 << i; if (cutOn && plab[o - 1] && (cutOn != 4 || plab[o - 1] <= lab0)) lb = 1; } }
  int deg = ((m >> 0) & 1) + ((m >> 2) & 1) + ((m >> 4) & 1) + ((m >> 6) & 1);
  if (deg <= 1) return 1;
  int links = 0;
  for (int i = 0; i < 8; i += 2) if ((m >> i & 1) && (m >> ((i + 1) & 7) & 1) && (m >> ((i + 2) & 7) & 1)) links++;
  if (links >= deg - 1 && !lb) return 1;       /* the quick answer takes touching for joined, so not next to labelled cells; */
  if (lb && K < 2) K = 2;                      /* there, search the box of radius 2 at least, which holds everything the quick answer looks at */
  if (K <= 1) return 0;
  stamp++;
  if (ex >= 0) { mark[ex] = stamp; lab[ex] = 255; }
  int ng = 0;
  for (int d = 0; d < 4; d++) {
    int o = occ(cx + DX[d], cy + DY[d]);
    if (!o || o - 1 == ex) continue;
    int j = o - 1, g;
    if (boxOnly && !held(cx + DX[d], cy + DY[d])) continue;
    if (ex >= 0 && cut(ex, j)) continue;       /* it was not joined through the cell that goes */
    g = ng++;
    qn[g] = 0; qh[g] = 0; Q[g][qn[g]++] = j; mark[j] = stamp; lab[j] = (unsigned char)g; par[g] = g;
  }
  if (ng <= 1) return 1;
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
        if (cut(j, k) || (boxOnly && !held(ax, ay))) continue;
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
  if (n > NW) return;
  int ic = rbelow(n), d = rbelow(4), ia = rbelow(n), g = rbelow(8), smax = (n - 1) / 2;
  int cx = X[ic], cy = Y[ic];
  if (!occ(cx + DX[d], cy + DY[d])) return;
  if (!(cutOn & 1)) { if (locallyConnected(cx, cy)) return; }
  else {                                        /* as above: the quick answer only away from labelled cells */
    int lb = plab[ic] != 0;
    for (int a = -1; a <= 1 && !lb; a++) for (int b = -1; b <= 1; b++) { int o = occ(cx + a, cy + b); if (o && plab[o - 1]) { lb = 1; break; } }
    if (!lb && locallyConnected(cx, cy)) return;
    if (cut(ic, occ(cx + DX[d], cy + DY[d]) - 1)) return;
  }
  stamp++;
  mark[ic] = stamp; lab[ic] = 255;
  int ngrp = 0, gd = -1;
  for (int e = 0; e < 4; e++) {
    int o = occ(cx + DX[e], cy + DY[e]);
    if (!o || cut(ic, o - 1)) continue;
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
        if (cut(j, k)) continue;
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
  for (int i = 0; i < nO; i++) for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) tubeMap[at(OX[i] + a, OY[i] + b)] = 0;
  n = 0;
  for (int i = 0; i < m; i++) {
    int x = X[i] >> 1, y = Y[i] >> 1;
    if (occ(x, y)) continue;
    X[n] = x; Y[n] = y; grid[at(x, y)] = n + 1; n++;
  }
  nO = n < NW ? n : 0;
  for (int i = 0; i < nO; i++) { OX[i] = X[i]; OY[i] = Y[i]; for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) tubeMap[at(X[i] + a, Y[i] + b)] |= 2; }
  for (int i = 0; i < nO; i++) tubeMap[at(X[i], Y[i])] |= 1;
  return n;
}
/* after coarsen(): the local walk adds cells only where the finer shape had cells (1), or there and at the
   sites round them (3); 0 lifts the rule */
void tube(int mode) { tubeOn = mode; }
/* with cover on, the walk also never takes a cell away if that would leave a site where the finer shape had cells with
   no cell on it or round it; w multiplies the weight of a cell (to make up for the room the rule takes away) */
void shadow(int cover, double w) { coverOn = cover; zmul = w; }
enum { OB = 64 }; static int orLab[OB * OB];  /* for the sites of the box: the label of the finer shape's part there after halving */
int *orlabels(void) { return orLab; }
static int covered(int qx, int qy, int ex) {
  int want = 0;                                 /* in the box, a part's place must be kept by that part */
  if (cutOn && held(qx, qy) && fx1 - fx0 <= OB) want = orLab[(qy - fy0) * (fx1 - fx0) + (qx - fx0)];
  for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) { int o = occ(qx + a, qy + b); if (o && o - 1 != ex && (!want || plab[o - 1] == want)) return 1; }
  return 0;
}

/* one move of the local walk: a cell added next to the shape or taken away, the count free */
static void gcMove(void) {
  const double Z = zmul / 4.0625696;
  nGC++;
  int i = rbelow(n), r = (int)(rnext() >> 61), cx = X[i], cy = Y[i];
  if (r < 4) {
    int sx = cx + DX[r], sy = cy + DY[r];
    if (occ(sx, sy) || n >= NMAX) return;
    if (regOn) { int r0 = rad(sx, sy); if (r0 < rin || r0 >= rout) return; }
    if (tubeOn && !(tubeMap[at(sx, sy)] & tubeOn)) return;
    int k = 0; for (int e = 0; e < 4; e++) if (occ(sx + DX[e], sy + DY[e])) k++;
    if (k > 1 && rbelow(k) != 0) return;
    double a = 4.0 * Z * n / (n + 1.0);
    if (a < 1.0 && runi() >= a) return;
    int nl = 0;
    if (cutOn) {                                /* a new cell in the box takes the label of what it touches there, and may not touch two */
      if (held(sx, sy)) {
        for (int e = 0; e < 4; e++) { int o = occ(sx + DX[e], sy + DY[e]); if (o && plab[o - 1]) {
          if (cutOn == 4 && plab[o - 1] > lab0) continue;
          if (!nl) nl = plab[o - 1]; else if (nl != plab[o - 1]) return; } }
        if (!nl) { if (cutOn == 5) return; nl = cutOn == 4 ? lab0 + 1 : labmax + 1; }
      }
    }
    if (!connAt(sx, sy, -1, drawK(256))) return;
    if (cutOn) { plab[n] = nl; if (nl > labmax) labmax = nl; if (cutOn == 4 && nl > 0 && nl <= lab0) labcnt[nl]++; }
    X[n] = sx; Y[n] = sy; grid[at(sx, sy)] = n + 1; n++; accGC++;
  } else {
    if (n <= 1 || i < nfix) return;
    if (cutOn == 4 && plab[i] > 0 && plab[i] <= lab0 && labcnt[plab[i]] <= 1) return;      /* the last cell of a piece of the inside stays */
    int k = 0; for (int e = 0; e < 4; e++) if (occ(cx + DX[e], cy + DY[e])) k++;
    if (!k) return;
    double a = n / (4.0 * Z * (n - 1.0));
    if (a < 1.0 && runi() >= a) return;
    if (coverOn) for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) if ((tubeMap[at(cx + a, cy + b)] & 1) && !covered(cx + a, cy + b, i)) return;
    if (!connAt(cx, cy, i, drawK(256))) return;
    if (cutOn == 5 && held(cx, cy)) { boxOnly = 1; int ok = connAt(cx, cy, i, 1 << 20); boxOnly = 0; if (!ok) return; }   /* and a part stays in one piece inside the box */
    grid[at(cx, cy)] = 0; n--;
    if (i != n) { X[i] = X[n]; Y[i] = Y[n]; grid[at(X[i], Y[i])] = i + 1; }
    if (cutOn) { if (cutOn == 4 && plab[i] > 0 && plab[i] <= lab0) labcnt[plab[i]]--; plab[i] = plab[n]; plab[n] = 0; }
    accGC++;
  }
}
/* the local walk for `quarters` quarter-sweeps (a sweep is one proposal per cell) */
void relax(int quarters) {
  for (int q = 0; q < quarters; q++) { int m = n / 4 > 1 ? n / 4 : 1; for (int k = 0; k < m; k++) gcMove(); }
}

/* cells added at random empty sites next to the shape, outside the held box, until there are `target` */
int pad(int target) {
  if (target > NW) target = NW;
  u64 tries = 0;
  while (n < target && tries++ < 400ULL * NW) {
    int i = rbelow(n), d = rbelow(4), sx = X[i] + DX[d], sy = Y[i] + DY[d];
    if (occ(sx, sy) || (frz && held(sx, sy))) continue;
    plab[n] = 0;
    X[n] = sx; Y[n] = sy; grid[at(sx, sy)] = n + 1; n++;
  }
  return n;
}
/* A better start than pad() alone for a world three times the size: turned copies of the shape itself are
   hung on its far tips, a cell at a time outward from the join, until there are `target` cells.  Sites that
   are taken, or in the held box, are skipped (the copy stays joined: it is added in search order and a cell
   goes in only next to one already there). */
int grow(int target) {
  if (target > NW) target = NW;
  int m = n, guard = 0;
  static int sx[NW], sy[NW], order[NW];
  for (int i = 0; i < m; i++) { sx[i] = X[i]; sy[i] = Y[i]; }
  while (n < target && guard++ < 64) {
    /* the cell of the present shape furthest from the root in a random direction, and a source cell far on the other side */
    double th = runi() * 6.283185307179586, ux = 1, uy = 0; { double c = 1, s = 0, a = th; /* cos, sin by series on a reduced angle */
      int q = (int)(a / 1.5707963267948966); a -= q * 1.5707963267948966; double a2 = a * a;
      c = 1 - a2 / 2 + a2 * a2 / 24 - a2 * a2 * a2 / 720 + a2 * a2 * a2 * a2 / 40320; s = a - a * a2 / 6 + a * a2 * a2 / 120 - a * a2 * a2 * a2 / 5040 + a * a2 * a2 * a2 * a2 / 362880;
      for (int k = 0; k < (q & 3); k++) { double t = c; c = -s; s = t; } ux = c; uy = s; }
    int tip = 0; double best = -1e30;
    for (int i = 0; i < n; i++) { double p = X[i] * ux + Y[i] * uy; if (p > best) { best = p; tip = i; } }
    int g = rbelow(8); const int *M = G8[g];
    int src = 0; best = 1e30;
    for (int i = 0; i < m; i++) { double p = (M[0] * sx[i] + M[1] * sy[i]) * ux + (M[2] * sx[i] + M[3] * sy[i]) * uy; if (p < best) { best = p; src = i; } }
    /* the copy's source cell lands on the tip; search the copy outward from there */
    int ox = X[tip] - (M[0] * sx[src] + M[1] * sy[src]), oy = Y[tip] - (M[2] * sx[src] + M[3] * sy[src]);
    /* a temporary index of the source cells by site: use mark[] as visited flags and a small open hash */
    enum { HS = 1 << 17 }; static int hk[HS];
    for (int i = 0; i < HS; i++) hk[i] = -1;
    for (int i = 0; i < m; i++) { u32 h = ((u32)sx[i] * 73856093u ^ (u32)sy[i] * 19349663u) & (HS - 1); while (hk[h] >= 0) h = (h + 1) & (HS - 1); hk[h] = i; }
    stamp++; int qh0 = 0, qn0 = 0; order[qn0++] = src; mark[src] = stamp;
    while (qh0 < qn0 && n < target) {
      int j = order[qh0++];
      int x = ox + M[0] * sx[j] + M[1] * sy[j], y = oy + M[2] * sx[j] + M[3] * sy[j];
      if (!occ(x, y)) {
        if (frz && held(x, y)) continue;
        int nb = 0; for (int e = 0; e < 4; e++) if (occ(x + DX[e], y + DY[e])) nb = 1;
        if (!nb) continue;
        plab[n] = 0; X[n] = x; Y[n] = y; grid[at(x, y)] = n + 1; n++;
      }
      for (int e = 0; e < 4; e++) {
        int qx = sx[j] + DX[e], qy = sy[j] + DY[e]; u32 h = ((u32)qx * 73856093u ^ (u32)qy * 19349663u) & (HS - 1); int k = -1;
        while (hk[h] >= 0) { if (sx[hk[h]] == qx && sy[hk[h]] == qy) { k = hk[h]; break; } h = (h + 1) & (HS - 1); }
        if (k >= 0 && mark[k] != stamp) { mark[k] = stamp; order[qn0++] = k; }
      }
    }
  }
  return n;
}
/* Joining a ring to the finer picture inside it.  The cell (ox, oy) of the coarser stage lay in the ring and
   touched a cell of the inside (r < b in the present, finer lattice) in direction (dx, dy).  Its four cells are now
   in place, and so is the finer picture.  If the four no longer touch it, the shortest path of empty sites inside
   to a cell of it, at most `reach` long, is filled.  With `want`, only a cell of that part (partlabels()) will do.
   Returns the number of cells added, or -1 if there is no such path. */
int link(int ox, int oy, int dx, int dy, int b, int reach, int want) {
  enum { R = 40, S = 2 * R + 1 };
  static int seen[S * S], from[S * S], qx[S * S], qy[S * S]; static int tick;
  if (reach > R - 1) reach = R - 1;
  int cx = dx ? (dx < 0 ? 2 * ox - 1 : 2 * ox + 2) : 2 * ox, cy = dy ? (dy < 0 ? 2 * oy - 1 : 2 * oy + 2) : 2 * oy, h = 0, m = 0; tick++;
  for (int t = 0; t < 2; t++) {            /* the two sites that face the four cells */
    int x = dx ? cx : cx + t, y = dy ? cy : cy + t, o = occ(x, y);
    if (rad(x, y) >= b) continue;
    if (o && (!want || plab[o - 1] == want)) return 0;
    if (o && want && plab[o - 1]) continue;                 /* a cell of another part is no way through: touching it joins nothing */
    int k = (y - cy + R) * S + (x - cx + R); seen[k] = tick; from[k] = -1; qx[m] = x; qy[m] = y; m++;
  }
  int depth = 1, levelEnd = m;
  while (h < m) {                          /* the search may pass through cells of the ring that are there; only empty sites get filled */
    if (h == levelEnd) { depth++; levelEnd = m; if (depth > reach) break; }
    int x = qx[h], y = qy[h], k = (y - cy + R) * S + (x - cx + R); h++;
    for (int e = 0; e < 4; e++) {
      int ax = x + DX[e], ay = y + DY[e], o;
      if (iabs(ax - cx) > R || iabs(ay - cy) > R || rad(ax, ay) >= b) continue;
      o = occ(ax, ay);
      if (o && (!want || plab[o - 1] == want)) {
        int added = 0;
        for (int j = k; j >= 0; j = from[j]) {
          int px = j % S - R + cx, py = j / S - R + cy;
          if (occ(px, py)) continue;
          if (n >= NMAX) return -1;
          plab[n] = want;
          X[n] = px; Y[n] = py; grid[at(px, py)] = n + 1; n++; added++;
        }
        return added;
      }
      if (o && want && plab[o - 1]) continue;
      int j = (ay - cy + R) * S + (ax - cx + R); if (seen[j] == tick) continue;
      seen[j] = tick; from[j] = k; qx[m] = ax; qy[m] = ay; m++;
    }
  }
  return -1;
}
/* one more cell, with a label (for the joins that are worked out outside); 0 if the site is taken */
int add(int x, int y, int lab_) {
  if (occ(x, y) || n >= NMAX) return 0;
  plab[n] = lab_; X[n] = x; Y[n] = y; grid[at(x, y)] = n + 1; n++; return 1;
}
/* acceptance counts, for checking against src/manyomino.cpp: 0 relocations tried, 1 accepted, 2 transplants tried, 3 accepted, 4 local tried, 5 accepted */
double stat(int which) { u64 v[6] = {nRel, accRel, nCP, accCP, nGC, accGC}; return (double)v[which]; }

#ifdef NATIVE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* outward_native exact NCELLS MOVES SEED [THIN]: every shape of NCELLS cells should turn up equally often
     (a shape is noted every THIN moves per cell, 10 unless given; too few and the notes are not independent)
   outward_native held NCELLS MOVES SEED [THIN]: the same with a 2x2 box held
   outward_native cut NCELLS MOVES SEED [THIN]: the same with two held cells that touch without being joined
   outward_native step IN OUT N S TC B SEED: one step outward from the shape IN (root at the origin) */
int main(int argc, char **argv) {
  init(12);
  if (argc > 4 && !strcmp(argv[1], "exact")) {
    int m = atoi(argv[2]), thin = argc > 5 ? atoi(argv[5]) : 10; long long moves = atoll(argv[3]); seed(1, atoi(argv[4]));
    for (int i = 0; i < m; i++) { X[i] = i; Y[i] = 0; } place(m);
    /* count shapes by a hash of their cells relative to the lower-left corner */
    enum { H = 1 << 20 }; static unsigned long long keyv[H]; static unsigned cnt[H]; unsigned long long samples = 0; int shapes = 0;
    for (long long t = 0; t < moves; t++) {
      run(1);
      if (t % ((long long)thin * m)) continue;
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
    int m = atoi(argv[2]), thin = argc > 5 ? atoi(argv[5]) : 10; long long moves = atoll(argv[3]); seed(2, atoi(argv[4]));
    static const int px[7] = {0, -1, -1, -1, 0, 1, 1}, py[7] = {0, 0, 1, 2, 2, 2, 1};
    for (int i = 0; i < 7; i++) { X[i] = px[i]; Y[i] = py[i]; } place(7); hold(0, 0, 2, 2, 1); pad(m);
    enum { H = 1 << 22 }; static unsigned long long keyv[H]; static unsigned cnt[H]; unsigned long long samples = 0; int shapes = 0;
    for (long long t = 0; t < moves; t++) {
      run(1);
      if (t % ((long long)thin * m)) continue;
      unsigned long long h = 0; for (int i = 0; i < n; i++) { unsigned long long c = (unsigned long long)(X[i] + 32) * 64 + (Y[i] + 32); c = (c + 1) * 0x9e3779b97f4a7c15ULL; c ^= c >> 29; c *= 0xbf58476d1ce4e5b9ULL; h += c ^ (c >> 32); }
      unsigned s = (unsigned)(h % H); while (cnt[s] && keyv[s] != h) s = (s + 1) % H;
      if (!cnt[s]) { keyv[s] = h; shapes++; } cnt[s]++; samples++;
    }
    double e = (double)samples / shapes, chi2 = 0; for (int s = 0; s < H; s++) if (cnt[s]) chi2 += (cnt[s] - e) * (cnt[s] - e) / e;
    printf("n=%d shapes=%d samples=%llu chi2=%.1f dof=%d; relocations accepted %.4f, transplants %.4f\n", m, shapes, samples, chi2, shapes - 1, stat(1) / stat(0), stat(3) / stat(2));
    return 0;
  }
  if (argc > 4 && !strcmp(argv[1], "cut")) {
    /* the box [0,2)x[0,1) holds two cells that touch but are of different parts, so the rest must join them:
       the shapes of NCELLS cells that do should turn up equally often (count them with analysis/held_count.py cut) */
    int m = atoi(argv[2]), thin = argc > 5 ? atoi(argv[5]) : 100; long long moves = atoll(argv[3]); seed(3, atoi(argv[4]));
    static const int px[6] = {0, 1, 0, 0, 1, 1}, py[6] = {0, 0, 1, 2, 2, 1};
    for (int i = 0; i < 6; i++) { X[i] = px[i]; Y[i] = py[i]; } place(6); plab[0] = 1; plab[1] = 2; parts(1); hold(0, 0, 2, 1, 1); pad(m);
    enum { H = 1 << 22 }; static unsigned long long keyv[H]; static unsigned cnt[H]; unsigned long long samples = 0; int shapes = 0;
    for (long long t = 0; t < moves; t++) {
      run(1);
      if (t % ((long long)thin * m)) continue;
      unsigned long long h = 0; for (int i = 0; i < n; i++) { unsigned long long c = (unsigned long long)(X[i] + 32) * 64 + (Y[i] + 32); c = (c + 1) * 0x9e3779b97f4a7c15ULL; c ^= c >> 29; c *= 0xbf58476d1ce4e5b9ULL; h += c ^ (c >> 32); }
      unsigned s = (unsigned)(h % H); while (cnt[s] && keyv[s] != h) s = (s + 1) % H;
      if (!cnt[s]) { keyv[s] = h; shapes++; } cnt[s]++; samples++;
    }
    double e = (double)samples / shapes, chi2 = 0; for (int s = 0; s < H; s++) if (cnt[s]) chi2 += (cnt[s] - e) * (cnt[s] - e) / e;
    printf("n=%d shapes=%d samples=%llu chi2=%.1f dof=%d; relocations accepted %.4f, transplants %.4f\n", m, shapes, samples, chi2, shapes - 1, stat(1) / stat(0), stat(3) / stat(2));
    return 0;
  }
  if (argc > 7 && !strcmp(argv[1], "mix")) {
    /* outward_native mix IN N B SEED PAD S1 S2 ...: pad the shape IN (root at the origin) to N cells, hold the box, and
       print the cells in rings round the root after S1, S2, ... moves per cell in all */
    FILE *f = fopen(argv[2], "rb"); int m; if (!f || fread(&m, 4, 1, f) != 1 || m > NMAX) return 1;
    static int buf[2 * NMAX]; if (fread(buf, 4, 2 * (size_t)m, f) != 2 * (size_t)m) return 1; fclose(f);
    for (int i = 0; i < m; i++) { X[i] = buf[2 * i]; Y[i] = buf[2 * i + 1]; }
    int N = atoi(argv[3]), B = atoi(argv[4]), padmode = atoi(argv[6]); seed(9, atoi(argv[5]));
    place(m); hold(-B, -B, B, B, 1); if (padmode == 1) grow(N); pad(N);
    static const double E[8] = {1, 1.5, 2, 2.5, 3, 4, 6, 8}; int done = 0;
    for (int a = 6; a < argc; a++) {
      if (a > 6) { int S = atoi(argv[a]); for (; done < S; done++) run(n); }
      long c[8] = {0}; double r2 = 0, mx = 0, my = 0;
      for (int i = 0; i < n; i++) {
        int rx = X[i] < 0 ? -X[i] - 1 : X[i], ry = Y[i] < 0 ? -Y[i] - 1 : Y[i], r = rx > ry ? rx : ry;
        if (r < B) c[7]++; else for (int k = 0; k < 7; k++) if (r >= E[k] * B && r < E[k + 1] * B) c[k]++;
        mx += X[i]; my += Y[i]; r2 += (double)X[i] * X[i] + (double)Y[i] * Y[i];
      }
      mx /= n; my /= n; r2 = r2 / n - mx * mx - my * my;
      printf("%d %d %ld %ld %ld %ld %ld %ld %ld %ld %.1f\n", done, n, c[7], c[0], c[1], c[2], c[3], c[4], c[5], c[6], r2);
    }
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
  fprintf(stderr, "usage: outward_native exact NCELLS MOVES SEED [THIN] | held NCELLS MOVES SEED [THIN] | cut NCELLS MOVES SEED [THIN] | step IN OUT N S TC B SEED\n"); return 1;
}
#endif
