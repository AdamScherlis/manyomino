// A library of shape observables, computed from nothing but the list of cells.
//
// Used to look for slow modes of the chain and to compare samplers.  Everything here is a
// function of the shape alone (not of cell order).  Three families:
//   geometry  radius of gyration, orientation, all central moments of the mass distribution up
//             to order 6 (in units of Rg), radial mass profile, bounding box;
//   local     perimeter, holes, degrees, 2x2 blocks: quantities that are sums over cells;
//   core      what one-cell cuts do to the shape.  For each cell c, d(c) is the number of cells
//             that come loose from the largest piece when c is removed.  From that: the
//             "centroid" cell (smallest largest piece) and how the mass splits there, the number
//             of cells whose removal detaches at least 1%, 3%, 10%, 30% of the shape (trunk
//             lengths), loop-connected blobs, and graph distances.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

struct Features {
  std::vector<int32_t> grid;
  std::vector<int> disc, low, par, sub, sepSum, maxSep, pit, stk, cstk, dist, que;

  static const std::vector<std::string> &names() {
    static std::vector<std::string> v;
    if (v.empty()) {
      v = {"rg2", "asph", "q2c", "q2s", "bw", "bh"};
      for (int o = 3; o <= 6; o++) for (int p = o; p >= 0; p--) v.push_back("m" + std::to_string(p) + std::to_string(o - p));
      for (const char *s : {"in05", "in10", "out15", "out20", "rmax",
                            "perim", "holes", "cycles", "leaves", "deg2", "deg3", "deg4", "straight", "blocks22", "noncut",
                            "cen_m1", "cen_m2", "cen_m3", "trunk01", "trunk03", "trunk10", "trunk30", "detach",
                            "blob_max", "blob_cells", "blob_count",
                            "depth_mean", "depth_max", "diam", "span", "cen_off", "deep_frac"}) v.push_back(s);
    }
    return v;
  }

  void compute(const int32_t *X, const int32_t *Y, int n, std::vector<double> &f) {
    f.clear();
    // ---- geometry
    double sx = 0, sy = 0; int xmin = X[0], xmax = X[0], ymin = Y[0], ymax = Y[0];
    for (int i = 0; i < n; i++) {
      sx += X[i] - X[0]; sy += Y[i] - Y[0];
      xmin = std::min(xmin, X[i]); xmax = std::max(xmax, X[i]); ymin = std::min(ymin, Y[i]); ymax = std::max(ymax, Y[i]);
    }
    const double mx = sx / n + X[0], my = sy / n + Y[0];
    double mom[7][7] = {{0}};
    for (int i = 0; i < n; i++) {
      double dx = X[i] - mx, dy = Y[i] - my, px[7], py[7];
      px[0] = py[0] = 1; for (int k = 1; k < 7; k++) { px[k] = px[k - 1] * dx; py[k] = py[k - 1] * dy; }
      for (int p = 0; p < 7; p++) for (int q = 0; p + q < 7; q++) mom[p][q] += px[p] * py[q];
    }
    for (int p = 0; p < 7; p++) for (int q = 0; p + q < 7; q++) mom[p][q] /= n;
    const double cxx = mom[2][0], cyy = mom[0][2], cxy = mom[1][1], rg2 = cxx + cyy, rg = std::sqrt(std::max(rg2, 1e-12));
    const double dsc = std::sqrt((cxx - cyy) * (cxx - cyy) + 4 * cxy * cxy);
    f.push_back(rg2); f.push_back(rg2 > 0 ? dsc * dsc / (rg2 * rg2) : 0);
    f.push_back((cxx - cyy) / rg2); f.push_back(2 * cxy / rg2);
    f.push_back((xmax - xmin + 1) / rg); f.push_back((ymax - ymin + 1) / rg);
    for (int o = 3; o <= 6; o++) for (int p = o; p >= 0; p--) f.push_back(mom[p][o - p] / std::pow(rg, o));
    long in05 = 0, in10 = 0, out15 = 0, out20 = 0; double rmax2 = 0;
    for (int i = 0; i < n; i++) {
      double dx = X[i] - mx, dy = Y[i] - my, r2 = (dx * dx + dy * dy) / rg2;
      in05 += r2 < 0.25; in10 += r2 < 1; out15 += r2 > 2.25; out20 += r2 > 4; rmax2 = std::max(rmax2, r2);
    }
    f.push_back((double)in05 / n); f.push_back((double)in10 / n); f.push_back((double)out15 / n); f.push_back((double)out20 / n); f.push_back(std::sqrt(rmax2));

    // ---- local: a dense grid over the bounding box, padded by one site.  0 empty, i+1 cell, -1 counted perimeter
    const int gw = xmax - xmin + 3, gh = ymax - ymin + 3;
    grid.assign((size_t)gw * gh, 0);
    auto at = [&](int i) { return (size_t)(Y[i] - ymin + 1) * gw + (X[i] - xmin + 1); };
    for (int i = 0; i < n; i++) grid[at(i)] = i + 1;
    const long nb[4] = {1, gw, -1, -(long)gw};
    long perim = 0, edges = 0, deg[5] = {0}, verts = 0, straight = 0, blocks = 0;
    for (int i = 0; i < n; i++) {
      const size_t g = at(i); int d = 0, m = 0;
      for (int k = 0; k < 4; k++) {
        int32_t &a = grid[g + nb[k]];
        if (a > 0) { d++; m |= 1 << k; } else if (a == 0) { a = -1; perim++; }
      }
      deg[d]++; edges += d; if (d == 2 && (m == 5 || m == 10)) straight++;
      if (grid[g + 1] > 0 && grid[g + gw] > 0 && grid[g + gw + 1] > 0) blocks++;
      // lattice vertices of the union of closed squares, each counted from the first cell around it
      for (int cy = 0; cy < 2; cy++) for (int cx = 0; cx < 2; cx++) {
        const size_t v = g + cx + (size_t)cy * gw;            // the cell whose lower-left corner is this vertex
        const size_t around[4] = {v - gw - 1, v - gw, v - 1, v};
        bool first = true;
        for (int k = 0; k < 4 && around[k] != g; k++) if (grid[around[k]] > 0) { first = false; break; }
        if (first) verts++;
      }
    }
    edges /= 2;
    const long holes = 1 - (verts - (4L * n - edges) + n);
    // ---- core: one depth-first search gives, for every cell, the pieces its removal leaves
    disc.assign(n, 0); low.assign(n, 0); par.assign(n, -1); sub.assign(n, 1); sepSum.assign(n, 0); maxSep.assign(n, 0); pit.assign(n, 0);
    stk.clear(); cstk.clear();
    int timer = 0; long blobMax = 1, blobCells = 0, blobCount = 0;
    stk.push_back(0); cstk.push_back(0); disc[0] = low[0] = ++timer;
    while (!stk.empty()) {
      const int v = stk.back();
      if (pit[v] < 4) {
        const int32_t a = grid[at(v) + nb[pit[v]++]];
        if (a <= 0) continue;
        const int w = a - 1;
        if (!disc[w]) { par[w] = v; disc[w] = low[w] = ++timer; stk.push_back(w); cstk.push_back(w); }
        else if (w != par[v]) low[v] = std::min(low[v], disc[w]);
      } else {
        stk.pop_back();
        const int p = par[v];
        if (p >= 0) {
          low[p] = std::min(low[p], low[v]); sub[p] += sub[v];
          if (low[v] >= disc[p]) { sepSum[p] += sub[v]; maxSep[p] = std::max(maxSep[p], sub[v]); }
        }
        if (p < 0 || low[v] > disc[p]) {        // the edge to the parent is a bridge: v closes a loop-connected blob
          long c = 0; int u;
          do { u = cstk.back(); cstk.pop_back(); c++; } while (u != v);
          blobMax = std::max(blobMax, c); if (c >= 4) { blobCells += c; blobCount++; }
        }
      }
    }
    long noncut = 0, trunk[4] = {0}; double detach = 0; int cen = 0; long cenL = n;
    const double thr[4] = {0.01 * n, 0.03 * n, 0.1 * n, 0.3 * n};
    for (int v = 0; v < n; v++) {
      const long L = std::max<long>(maxSep[v], n - 1 - sepSum[v]), d = n - 1 - L;
      if (d == 0) noncut++;
      detach += d;
      for (int k = 0; k < 4; k++) if (d >= thr[k]) trunk[k]++;
      // ties go to the lowest, then leftmost, cell: a property of the shape, not of the cell numbering
      if (L < cenL || (L == cenL && (Y[v] < Y[cen] || (Y[v] == Y[cen] && X[v] < X[cen])))) { cenL = L; cen = v; }
    }
    // the pieces at the centroid cell
    long pc[8] = {n - 1 - sepSum[cen], 0, 0, 0, 0, 0, 0, 0}; int npc = 1;
    for (int k = 0; k < 4; k++) {
      const int32_t a = grid[at(cen) + nb[k]];
      if (a > 0 && par[a - 1] == cen && low[a - 1] >= disc[cen]) pc[npc++] = sub[a - 1];
    }
    std::sort(pc, pc + npc, std::greater<long>());
    // ---- graph distances: from the centroid cell, then a double sweep for the diameter
    auto bfs = [&](int s, double *mean) {
      dist.assign(n, -1); que.clear(); que.push_back(s); dist[s] = 0; double tot = 0; size_t h = 0;
      while (h < que.size()) {
        const int v = que[h++]; tot += dist[v];
        for (int k = 0; k < 4; k++) { const int32_t a = grid[at(v) + nb[k]]; if (a > 0 && dist[a - 1] < 0) { dist[a - 1] = dist[v] + 1; que.push_back(a - 1); } }
      }
      if (mean) *mean = tot / n;
      int far = que.back();                 // farthest cell; ties to the lowest, then leftmost
      for (size_t i = que.size(); i-- > 0 && dist[que[i]] == dist[que.back()];) { const int v = que[i]; if (Y[v] < Y[far] || (Y[v] == Y[far] && X[v] < X[far])) far = v; }
      return far;
    };
    double dmean = 0;
    const int a = bfs(cen, &dmean); const int dmax = dist[a];
    long deep = 0; for (int i = 0; i < n; i++) deep += 2 * dist[i] > dmax;
    const int b = bfs(a, nullptr); int diam = dist[b];
    const int c2 = bfs(b, nullptr); diam = std::max(diam, dist[c2]);
    const double span = std::hypot((double)X[a] - X[b], (double)Y[a] - Y[b]) / rg, cenOff = std::hypot(X[cen] - mx, Y[cen] - my) / rg;

    for (double v : {(double)perim / n, (double)holes / n, (double)(edges - (n - 1)) / n, (double)deg[1] / n, (double)deg[2] / n, (double)deg[3] / n,
                     (double)deg[4] / n, (double)straight / n, (double)blocks / n, (double)noncut / n,
                     (double)pc[0] / n, (double)pc[1] / n, (double)pc[2] / n, (double)trunk[0], (double)trunk[1], (double)trunk[2], (double)trunk[3],
                     detach / ((double)n * n), (double)blobMax, (double)blobCells / n, (double)blobCount / n,
                     dmean, (double)dmax, (double)diam, span, cenOff, (double)deep / n}) f.push_back(v);
  }
};
