# Handoff: uniform sampling of very large random polyominoes

## Context

Goal: sample polyominoes (lattice animals on Z², 4-connectivity) **uniformly at random for fixed size n**, at n from ~100 up to ~10,000, and render them. A previous hand-rolled MCMC sampler converged too slowly; this doc specifies a chain that is provably uniform, the engineering needed to make it fast, and the validation that proves it mixed. Target distribution is uniform over **fixed** polyominoes (translation classes), which is what the chain below naturally gives — validate against OEIS **A001168** (fixed), not A000105 (free).

## Deliverables

1. A fast sampler (local-move MCMC as the baseline; see stretch goals for alternatives).
2. A validation suite: exact-count cross-checks, uniformity chi-square, two-seed convergence test, autocorrelation/ESS.
3. A gallery of rendered samples at n ∈ {1000, 3000, 10000}, ≥4 independent samples each.
4. A measured radius-of-gyration exponent ν from the samples.

## What correct output looks like (use as a smell test)

Uniform random polyominoes are in the **branched-polymer universality class**: sparse, wispy, tree-like dendrites with small loops as local decoration. Radius of gyration R_g ~ n^ν with ν ≈ 0.64 in 2D (numerical estimates ≈ 0.6408), i.e. fractal dimension ≈ 1.56. There is no deterministic limit shape; the count grows like λⁿ/n with Klarner's constant λ ≈ 4.0626. **If samples look like compact round blobs, the chain has not mixed or there is a bug.**

## Core algorithm — fixed-n single-cell chain

State: a set A ⊂ Z², |A| = n, 4-connected. One move:

1. Pick a cell c uniformly from the n cells of A.
2. If A − {c} is disconnected: **stay put** (this still counts as a step).
3. Otherwise pick s uniformly from P(A−c), the set of empty sites 4-adjacent to A−c. Note c itself is always in P(A−c), so s = c (a no-op) is a legal outcome and must be allowed.
4. New state: (A − {c}) ∪ {s}.

**Why it's uniform:** for a legal move A → A′ with A′ = (A−c) ∪ {s}, the forward probability is (1/n)·(1/|P(A−c)|). The reverse move removes s and adds c, with probability (1/n)·(1/|P(A′−s)|). Since A−c = A′−s, these are equal: the kernel is symmetric, so the uniform distribution is stationary and **every legal proposal is accepted** — there is no Metropolis ratio to compute. Rejections at step 2 give aperiodicity; irreducibility holds because a non-cut cell always exists (any polyomino can be disassembled into a straight bar by repeatedly relocating leaf cells).

### Correctness details that are easy to get wrong

- **Uniform over perimeter *sites*, not (cell, direction) pairs.** Sampling a random occupied cell and a random direction over-weights sites adjacent to multiple occupied cells and silently biases the chain. Maintain an explicit set of perimeter sites and sample from it uniformly.
- **The perimeter is that of A−c, not A.** Remove c first (updating the perimeter), then sample s, then insert s.
- The chain's centroid drifts (random walk in Z²). Use i64 coordinates and re-center occasionally / at render time. Canonicalize by translating the min corner to the origin only when hashing shapes for validation.
- O(1) uniform sampling from a dynamic set: the standard indexed-set structure (dense array + position hash map, swap-remove on delete). Use it for both the cell set and the perimeter set.

### Fast connectivity test (step 2)

- If c has exactly **one** occupied 4-neighbor, it is never a cut cell → removable, no further work.
- Otherwise do the **3×3 ring test**: if the occupied cells in c's 8-neighborhood form a single connected arc (within the ring) containing all of c's occupied 4-neighbors, removal is safe. This is sufficient but not necessary.
- Otherwise fall back to BFS/DFS in A − {c} starting from one occupied 4-neighbor of c, with early exit once all the other occupied 4-neighbors of c have been reached; c is a cut cell iff the search exhausts without finding them. Worst case O(n), but the cheap tests catch most cells. Dynamic-connectivity structures exist if profiling demands, but almost certainly skip them at these sizes.

## Performance budget

Heuristic relaxation cost for local moves: mass moves diffusively, so full shape decorrelation ≈ n^(1+2ν) ≈ n^2.3 **total moves**. Concretely: n = 3000 → ~10^8 moves; n = 10^4 → ~1.6×10^9 moves. Therefore:

- Write the hot loop in a fast language — **Rust preferred** (or C / Cython / numba). Target ≥10^7 moves/sec, which decorrelates n = 10^4 in minutes. Pure Python (~10^5 moves/sec) is fine only as the reference implementation.
- Build **two implementations**: a slow, obviously-correct Python reference and the fast one. Differential-test them (invariants after every k moves: |A| = n, connectivity, perimeter set equals a from-scratch recompute; plus matching move statistics).

## Validation plan (do this before any big runs)

1. **Exact counts** (fixed polyominoes, A001168), n = 1..10: 1, 2, 6, 19, 63, 216, 760, 2725, 9910, 36446. Write a brute-force enumerator for n ≤ 8–10 as an independent check.
2. **Uniformity chi-square** at n = 5 (63 shapes) and n = 6 (216 shapes): run the chain, record the translation-canonical form every ~10 moves, collect ≥10^6 observations, chi-square against uniform. This catches essentially every proposal-distribution bug, including the perimeter-pair bug above.
3. **Two-seed convergence test** at each production n: chain A seeded from a 1×n bar, chain B from a roughly √n × √n rectangle. Run until the R_g statistics (running mean ± stderr) of the two chains agree. This is the honest mixing test — do not trust single-chain output before it passes.
4. **Autocorrelation**: integrated autocorrelation time τ of R_g² at several n; report effective sample size = N/(2τ). Expect τ to scale roughly like n^2ν in sweeps (n^(1+2ν) in moves).
5. **Exponent check**: fit log⟨R_g⟩ vs log n over n ∈ {100, 200, 500, 1000, 2000, 5000} → slope should come out ≈ 0.64 within a few percent.

## Stretch goals (only if local moves are too slow in practice)

- **Nonlocal cut-and-paste moves** in the style of Janse van Rensburg & Madras's lattice-tree algorithms: detach a branch at an articulation cell, reattach it elsewhere with a random rotation/reflection, with a Hastings ratio counting valid attachment placements in both directions. This is the known fix for slow shape relaxation, but the accept ratio is fiddly to get unbiased — after adding it, re-run validation step 2 before believing anything.
- **PERM** (Rosenbluth growth with pruning/enrichment, Grassberger-style) as an independent non-MCMC sampler. It sidesteps mixing entirely and has different failure modes, so agreement of its R_g statistics with the MCMC is a strong cross-check.

## Rendering

PNG, black cells on white, tight bounding box, ~1 px per cell at large n. Gallery: n ∈ {1000, 3000, 10000}, ≥4 independent samples each (independent = separate chains or ≥ several τ apart). Optional but nice: color cells by graph distance from a root cell to make the branch structure visible.

## Suggested order of work

1. **M1** — Python reference chain + enumerator + validation steps 1–2 pass.
2. **M2** — Fast implementation + differential tests against the reference.
3. **M3** — Diagnostics: steps 3–5 (two-seed, autocorrelation, ν fit).
4. **M4** — Gallery renders.
5. **M5** — Optional: cut-and-paste moves and/or PERM, each followed by re-validation.

Definition of done: validation steps 1–3 pass, ν ≈ 0.64 recovered, gallery produced, README documents how to reproduce each figure.

## Pitfalls checklist

- Sampling perimeter (cell, direction) pairs instead of sites — biases the chain; the n = 5/6 chi-square catches it.
- Computing the perimeter of A instead of A−c during the proposal.
- Forgetting that s = c (no-op) must be a legal outcome — excluding it breaks symmetry.
- Validating against free polyominoes (A000105) instead of fixed (A001168).
- BFS early-exit logic inverted (declaring "cut" when the search merely exited early).
- Incremental perimeter updates drifting out of sync — differential-test against full recomputes.
- Comparing "independent" samples that are actually < τ apart, or reusing RNG streams across chains.
