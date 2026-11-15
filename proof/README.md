# A quasi-polynomial bound on the relaxation time

**Status.** A complete argument, written out below, using standard tools. Its combinatorial
steps are machine-checked for every pair of shapes with up to 7 cells
(`python3 proof/check_paths.py 7`). It has not been read by anyone else. A literature search
found no earlier bound of this kind for polyominoes or lattice animals, but that search was
short and the claim of novelty is correspondingly weak.

## Statement

Let Ω_n be the set of fixed polyominoes with n cells: finite connected sets of n cells of the
square lattice, two sets being the same polyomino if one is a translate of the other. Write
A(n) = |Ω_n|.

**The relocation chain.** From a polyomino P: choose a cell c of P uniformly, and
independently a cell u of P and one of the four directions d, uniformly. Let s be the
neighbour of u in direction d. If u ≠ c, s is not in P∖{c}, and P∖{c} is connected, then with
probability 1/k, where k is the number of cells of P∖{c} adjacent to s, move to
(P∖{c}) ∪ {s}. Otherwise stay.

For a given pair (c, s), with P∖{c} connected and s an empty site next to P∖{c}, the move has
probability (1/n) · (k/4n) · (1/k) = 1/(4n²). Reversing a move is a move of the same kind, so the
transition matrix is symmetric and the uniform distribution on Ω_n is stationary.

**Theorem.** Let λ₂ be the second-largest eigenvalue of the relocation chain on Ω_n, n ≥ 2. Then

    1 / (1 − λ₂)  ≤  4 n⁵ (L+1)² (2n+2)^(6L),      L = 2 log₂ n + 4.

In particular the relaxation time is n^(O(log n)).

**Consequences.**

1. *No exponentially slow property.* For any function f on Ω_n, the integrated
   autocorrelation time of f along the stationary chain, τ(f) = ½ + Σ_{t≥1} ρ_f(t), is at most
   1/(1 − λ₂). (Write f in the eigenbasis: τ(f) is a weighted average of ½(1+λ)/(1−λ) over the
   eigenvalues λ < 1, and each term is at most 1/(1−λ₂).) So every property of the shape, however
   contrived, is forgotten within n^(O(log n)) moves in this sense.
2. *The sampler in this repository.* Its relocation move makes the same proposal, and accepts
   it when connectivity can be certified within a random radius K with P(K ≥ r) = r⁻²; for
   K ≥ n the certificate is exact. So each transition of the chain above has probability at
   least n⁻² times as large in the sampler's relocation move. The sampler mixes that move with
   transplants, with weight w = 0.6 on relocation. Both kernels are reversible with respect to
   the uniform distribution, so comparing Dirichlet forms gives
   1 − λ₂(sampler) ≥ w n⁻² (1 − λ₂), and the same n^(O(log n)) bound holds for the sampler.
3. *Mixing time.* For the lazy version of the chain, the total-variation mixing time from any
   start is at most (1 − λ₂)⁻¹ (ln A(n) + ln(1/ε)), again n^(O(log n)).

The bound is far from what is observed (about 46 n moves for the sampler, n^2.6 for relocation
alone). What it rules out is a bottleneck: no property can take exponentially long.

## Proof

The proof uses canonical paths (Jerrum and Sinclair). For a reversible chain with stationary
distribution π, choose for every ordered pair x ≠ y a path γ_xy of transitions from x to y.
Then (Sinclair 1992)

    1 / (1 − λ₂)  ≤  ρ̄  :=  max over transitions e = (z, z′) of
                            (1 / (π(z) P(z, z′))) · Σ_{(x,y): γ_xy ∋ e} π(x) π(y) |γ_xy| .

Here π is uniform, P(z, z′) ≥ 1/(4n²), and the paths built below have length at most n, so

    ρ̄  ≤  4 n³ · N / A(n),        N := the largest number of paths through one transition.      (1)

The whole proof is the construction of paths with N ≤ n^(O(log n)) · A(n).

### 1. An order on the cells of a polyomino

Let P be a polyomino and r one of its cells. Fix a spanning tree T of P rooted at r (say the
breadth-first tree, with a fixed rule for ties). For a cell v, its *subtree* is v together with
all its descendants in T.

**Definition.** σ(P, r) is the depth-first preorder of T in which the children of every cell are
visited in increasing order of subtree size (ties by a fixed rule). So the child with the
largest subtree is visited last.

**Lemma 1.** Let σ(P, r) = (σ₁, …, σ_n) and 1 ≤ m < n.
(a) The prefix {σ₁, …, σ_m} is connected and contains r.
(b) The suffix {σ_{m+1}, …, σ_n} is a union of at most L = 2 log₂ n + 4 subtrees of T. Each is a
connected set of cells, and each contains a cell adjacent to the prefix.

*Proof.* (a) A preorder lists every cell after its parent, so the prefix is closed under taking
parents; it is connected through tree edges.

(b) Let v = σ_m and let r = w₀, w₁, …, w_t = v be the tree path to v. At the moment v has just
been visited, the unvisited cells are exactly: the subtrees of the children of v; and, for
each i < t, the subtrees of those children of w_i that come after w_{i+1} in the order. These
are whole subtrees; each is connected through tree edges, and its root is adjacent to its
parent, which is in the prefix.

Count them. The cell v has at most 3 children (4 if v = r, in which case there is no path
above it). A cell w_i has at most 2 children other than w_{i+1} (3 if w_i = r). And w_i has a child
after w_{i+1} only if w_{i+1} is not its last child, in which case some later child has a subtree
at least as large as that of w_{i+1}, so

    size(w_{i+1}) ≤ (size(w_i) − 1) / 2 .

Sizes decrease along the path from n to at least 1, so this halving happens for at most log₂ n
indices i. The total is at most 2 log₂ n + 1 + 3. ∎

### 2. The paths

For a polyomino A let r_A be its rightmost cell, and the topmost of those if there are several.
For a polyomino B let l_B be its leftmost cell, the lowest of those if there are several.

Given A ≠ B, place B so that l_B is the right-hand neighbour of r_A. Every cell of A then has
x-coordinate at most that of r_A, and every cell of B has a larger one, so A and B are disjoint.
Let

    α = σ(A, r_A) = (α₁ = r_A, …, α_n),        β = σ(B, l_B) = (β₁ = l_B, …, β_n),

    Z_k = {α₁, …, α_{n−k}} ∪ {β₁, …, β_k},        k = 0, 1, …, n.

So Z₀ = A and Z_n = B: the path takes A apart in the reverse of the order α while building B
in the order β.

**Lemma 2.** For 1 ≤ k ≤ n, Z_k is obtained from Z_{k−1} by a move of the relocation chain:
remove c = α_{n−k+1}, add s = β_k.

*Proof.* The set Z_{k−1}∖{c} = {α₁..α_{n−k}} ∪ {β₁..β_{k−1}} is connected: each of the two
parts is a prefix, hence connected by Lemma 1(a), and when both are non-empty they are joined
because α₁ and β₁ are neighbours. (It is non-empty since n ≥ 2.) The site s is not in Z_{k−1}
because A and B are disjoint. And s is adjacent to Z_{k−1}∖{c}: for k ≥ 2 its parent in the tree
of B is in {β₁..β_{k−1}}, and for k = 1 it is l_B, adjacent to α₁ = r_A, which is still present
because n − 1 ≥ 1. ∎

If two consecutive sets are translates of each other the step is dropped; the path γ_AB has
at most n transitions.

### 3. Counting the paths through a transition

Fix a transition from Z to Z′, and fix Z as a concrete set of cells. Suppose γ_AB uses it at
step k, so that (after translating A and B together) Z = Z_{k−1}. Let x₀ be the x-coordinate
of r_A. Then

    S_A := {α₁..α_{n−k+1}} = the cells of Z with x ≤ x₀,     S_B := {β₁..β_{k−1}} = the rest of Z,
    A = S_A ∪ R_A,   R_A := {α_{n−k+2}..α_n},        B = S_B ∪ R_B,   R_B := {β_k..β_n}.

So the pair (A, B) is determined by (k, x₀, R_A, R_B), and it remains to count the possible
values of this quadruple. There are at most n values of k and at most n of x₀ (it is the
x-coordinate of a cell of Z).

By Lemma 1(b), R_A is a union of j ≤ L connected sets, each containing a site that is empty in
Z and adjacent to Z (the root of the corresponding subtree; it is adjacent to S_A, is not in S_A,
and is not in S_B because its x-coordinate is at most x₀). The same holds for R_B, with j′ ≤ L
sets. (When k = 1, S_B is empty and R_B = B is a single connected set containing the site to
the right of r_A, which is again an empty site adjacent to Z.)

Two standard facts:

* A set of n cells has at most 2n + 2 empty sites adjacent to it.
* The number of connected sets of m cells that contain a given site is m · A(m): each
  polyomino of m cells has m translates containing the site.

So the number of choices for R_A with given j and given sizes m₁, …, m_j of the pieces is at
most ∏ (2n+2) m_i A(m_i), and likewise for R_B with sizes m′₁, …, m′_j′. There are at most n^j
and n^j′ choices of sizes, and each m_i ≤ n. Therefore

    #(R_A, R_B)  ≤  Σ_{j, j′ ≤ L} ((2n+2) n²)^(j+j′) · max ∏ A(m_i) ∏ A(m′_i) ,

where the sizes in the last product add up to (k−1) + (n−k+1) = n.

**Lemma 3 (Klarner).** A(a) A(b) ≤ A(a+b). Hence ∏ A(m_i) ≤ A(Σ m_i).

*Proof.* Given polyominoes P and Q of a and b cells, place Q so that its leftmost-lowest cell
is the right-hand neighbour of the rightmost-topmost cell of P. The union is a polyomino of
a + b cells, and P is recovered from it as its first a cells in the order "by column, then by
row". So the map (P, Q) ↦ P ∪ Q is injective. ∎

Putting these together,

    N  ≤  n · n · (L+1)² · ((2n+2) n²)^(2L) · A(n)  ≤  n² (L+1)² (2n+2)^(6L) · A(n),

and with (1),

    1 / (1 − λ₂)  ≤  ρ̄  ≤  4 n⁵ (L+1)² (2n+2)^(6L).   ∎

## Remarks

* **Where the logarithm comes from.** If the removed part of A were always a single connected
  piece, the count would be polynomial. That would need an order in which every prefix *and*
  every suffix is connected, which a tree with three long arms does not have. Visiting the
  heaviest subtree last is the classical way to keep the number of pending pieces
  logarithmic, and each pending piece costs a polynomial factor (where it attaches, how large
  it is). For these particular paths the loss looks real: a heuristic count (a long bare spine
  with logarithmically many heavy branches pending along it) gives transitions carrying
  n^(c log n) · A(n) of them. Whether the relaxation time is polynomial is open; it would need
  better paths, or a flow spread over many orders.
* **Other lattices and dimensions.** Nothing is special to the square lattice. On Z^d the same
  proof works with larger constants (a cell has up to 2d − 1 children, a set of n cells at most
  2dn empty neighbours); Lemma 3 holds for lattice animals in every dimension.
* **What is not covered.** The percolation-weighted version of the sampler (weights can differ
  by a factor exponential in n between shapes). And the bound says nothing useful about
  constants: at n = 10⁶ it is astronomically larger than anything a computer will run.

## The machine check

`proof/check_paths.py` enumerates all polyominoes of n cells and all ordered pairs (A, B), builds
the paths exactly as defined above, and asserts: every step is a legal relocation and every
intermediate set is connected (Lemma 2); the removed and missing parts have at most L pieces,
each touching the rest (Lemma 1); the quadruple (k, x₀, R_A, R_B), written relative to Z, is
different for different pairs through the same transition (the injectivity in section 3); the
number of paths through a transition is below the bound; and the exact relaxation time, from
the eigenvalues of the transition matrix, is below the path bound, which is below the theorem's.

| n | shapes | paths | most paths through one transition ÷ A(n) | exact 1/(1−λ₂) | canonical-path bound ρ̄ |
|---:|---:|---:|---:|---:|---:|
| 2 | 2 | 2 | 0.50 | 2.0 | 4 |
| 3 | 6 | 30 | 2.00 | 4.5 | 108 |
| 4 | 19 | 342 | 4.21 | 10.2 | 903 |
| 5 | 63 | 3,906 | 6.95 | 19.8 | 3,429 |
| 6 | 216 | 46,440 | 10.11 | 34.1 | 8,732 |
| 7 | 760 | 576,840 | 14.15 | 53.5 | 19,410 |

This checks the construction, not the asymptotics: at these sizes there are never more than 3
pending pieces.

## References

* A. Sinclair, "Improved bounds for mixing rates of Markov chains and multicommodity flow,"
  *Combinatorics, Probability and Computing* 1 (1992).
* M. Jerrum and A. Sinclair, "Approximating the permanent," *SIAM J. Comput.* 18 (1989).
* D. A. Klarner, "Cell growth problems," *Canadian J. Math.* 19 (1967).
* M. Bressan, "Efficient and near-optimal algorithms for sampling small connected subgraphs,"
  *ACM Trans. Algorithms* (2023): the nearest general result found; its bound for connected
  k-vertex subgraphs carries a factor k^(O(k)), so it does not give this one.
