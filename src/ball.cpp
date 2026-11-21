// ball: cut a connected piece of N cells out of a shape, around a box.
//
// Given a shape, a box (centre and half-width B) and N, writes N cells of the shape that contain
// everything in the box and are connected: first the paths through the shape that join the
// pieces of the box's content to one another, then the cells nearest the box (in steps through
// the shape) up to N.  It is the starting state for sampling the surroundings of the box with
// the box held (manyomino -freeze).  Exits with status 2 if the joining paths alone need more
// than N cells, or the pieces do not join within 60 N cells.
//   ball IN OUT CX CY B N        coordinates in OUT are relative to (CX, CY)
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <deque>
typedef uint64_t u64;
static inline u64 key(int x, int y) { return ((u64)(uint32_t)x << 32) | (uint32_t)y; }
int main(int argc, char **argv) {
  if (argc < 7) { fprintf(stderr, "usage: ball IN OUT CX CY B N\n"); return 1; }
  FILE *f = fopen(argv[1], "rb"); if (!f) return 1;
  int32_t n; if (fread(&n, 4, 1, f) != 1) return 1; std::vector<int32_t> v(2 * (size_t)n); if (fread(v.data(), 4, 2 * (size_t)n, f) != 2 * (size_t)n) return 1; fclose(f);
  int cx = atoi(argv[3]), cy = atoi(argv[4]), B = atoi(argv[5]); size_t N = strtoull(argv[6], 0, 10);
  std::unordered_map<u64, int> idx; idx.reserve(2 * (size_t)n);
  for (int i = 0; i < n; i++) idx[key(v[2 * i], v[2 * i + 1])] = i;
  static const int DX[4] = {1, -1, 0, 0}, DY[4] = {0, 0, 1, -1};
  auto nb = [&](int i, int d) { auto it = idx.find(key(v[2 * i] + DX[d], v[2 * i + 1] + DY[d])); return it == idx.end() ? -1 : it->second; };
  auto inbox = [&](int i) { int x = v[2 * i] - cx, y = v[2 * i + 1] - cy; return x >= -B && x < B && y >= -B && y < B; };
  std::vector<int> src(n, -1), up(n, -1), order; std::vector<char> skel(n, 0);
  // pieces of the box's content
  int npc = 0; std::vector<int> st;
  for (int i = 0; i < n; i++) if (inbox(i) && src[i] < 0) {
    src[i] = npc; st.push_back(i);
    while (!st.empty()) { int a = st.back(); st.pop_back(); order.push_back(a); skel[a] = 1; for (int d = 0; d < 4; d++) { int b = nb(a, d); if (b >= 0 && inbox(b) && src[b] < 0) { src[b] = npc; st.push_back(b); } } }
    npc++;
  }
  if (!npc) { fprintf(stderr, "ball: nothing in the box\n"); return 2; }
  std::vector<int> par(npc); for (int i = 0; i < npc; i++) par[i] = i;
  auto find = [&](int a) { while (par[a] != a) a = par[a] = par[par[a]]; return a; };
  std::deque<int> q(order.begin(), order.end()); int comps = npc; size_t nskel = order.size();
  auto trace = [&](int a) { while (a >= 0 && !skel[a]) { skel[a] = 1; nskel++; a = up[a]; } };
  while (!q.empty() && (comps > 1 || order.size() < N) && order.size() < 60 * N) {
    int a = q.front(); q.pop_front();
    for (int d = 0; d < 4; d++) {
      int b = nb(a, d); if (b < 0) continue;
      if (src[b] < 0) { src[b] = src[a]; up[b] = a; order.push_back(b); q.push_back(b); }
      else if (comps > 1) { int ra = find(src[a]), rb = find(src[b]); if (ra != rb) { par[ra] = rb; comps--; trace(a); trace(b); } }
    }
  }
  if (comps > 1 || nskel > N) return 2;
  std::vector<int32_t> o; o.push_back((int32_t)N); size_t have = nskel;
  for (int i : order) { if (skel[i]) continue; if (have >= N) break; skel[i] = 1; have++; }      // a cell's parent came earlier, so the set stays connected
  if (have != N) return 2;
  for (int i : order) if (skel[i]) { o.push_back(v[2 * i] - cx); o.push_back(v[2 * i + 1] - cy); }
  f = fopen(argv[2], "wb"); fwrite(o.data(), 4, o.size(), f); fclose(f);
  return 0;
}
