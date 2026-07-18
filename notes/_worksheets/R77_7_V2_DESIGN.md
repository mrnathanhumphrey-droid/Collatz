# R77.7 V2 — Phase 2: New solver design

**Date:** 2026-05-12. Wilson.

## Disposition of solver choice

**H_CRT_SOLVER_WORKS adopted** (per pre-registration). The Markov kernel is, in fact, **dense** at the level of coprime states — every row has exactly N nonzeros (M = N at q=3, since M = 2·3^(k-1) = N_coprime). Sparsity-exploiting solvers (sparse-LU, Wiedemann, Krylov) give no advantage here; the bottleneck is denominator blowup, not nonzero count.

Therefore the chosen approach is:

> **Dense Gauss elimination over F_p, repeated for multiple primes, then CRT + rational reconstruction.**

This is the project-canonical move when an O(N³)-with-Q solve hits denominator blowup. Each prime gives a numerically-stable O(N³) elimination using only int64 operations (no bigint). After ~25-30 primes the full rational result is recoverable.

## Mathematical setup

We want π_k ∈ Q^N such that

    π_k^T K_k = π_k^T,    π_k > 0,    Σ π_k = 1

with K_k constructed exactly as in `build_markov_rational(k)`. Equivalently, (K_k^T − I) π = 0 with Σπ=1.

Replace the last row of (K_k^T − I) with the all-ones constraint (RHS=1, all others 0):

    A π = e_N            where A = K^T − I with last row replaced

This is the same construction as `stationary_rational` in the original solver.

The matrix entries of K are rationals of the form (1/2^r_v) / Z_v = 2^(M−r_v) / (2^M − 1). They are all over the same denominator 2^M − 1. So we can clear denominators: let

    K = (1 / (2^M − 1)) · K_int

where K_int has integer entries (each entry is a sum of terms 2^(M − r_v) for those r_v whose target hit this column). Then

    (K^T − I) π = 0   ⇔   (K_int^T − (2^M − 1) I) π = 0

This means we can work with an INTEGER matrix in the entire pipeline, no fractions at all in the matrix setup. The stationary equation becomes a linear system over Z with integer coefficients of size at most 2^M. The solution π is still rational of course, but the matrix is integer.

For each prime p:

    A_p = ((K_int)^T − (2^M − 1) I) mod p,   last row replaced with all-ones
    b_p = (0, 0, ..., 0, 1) mod p
    π_p = A_p^{-1} b_p   (in F_p^N)

Then assemble {π_p}_p via CRT to integer π_int, then rational-reconstruct π_int → π_Q.

## Algorithmic components

### 1. Build K_int

Iterate r in coprime states, r_v in 1..M, target = ((3r+1)·inv2^{r_v}) mod 3^k, weight = 2^(M − r_v). Accumulate K_int[r_idx, target_idx] += weight as Python int (since 2^M with M=1458 is a 439-digit bigint, we'll store these accumulated weights as Python ints — but only ONCE at build time; never again during the solve loop). Total accumulation size per entry: bounded by sum of M weights ≤ 2^M.

### 2. Per-prime modular solve

For each prime p (well below 2^63 so int64 arithmetic is safe), build A_p as a numpy `int64` matrix:

    A_p[i, j] = K_int[j, i] mod p       (note transpose)
    A_p[i, i] -= (2^M − 1) mod p
    A_p[N-1, :] = 1                     (last-row replacement)
    b_p[N-1] = 1, others = 0

Then dense Gauss elimination over F_p. The inner ops are:
- Modular inverse of pivot (via `pow(pivot, p-2, p)` — Fermat's little theorem, O(log p) bit ops)
- Row scale: multiply row by inv_pivot mod p, single numpy int64 vector op
- Row subtract: `A_p[other] -= factor * A_p[pivot]; A_p[other] %= p`

Important: to keep numpy int64 from overflowing, we need `p^2 < 2^63 ≈ 9.2·10^18`, so `p < 3·10^9`. We'll choose primes in [10^8, 10^9] giving headroom and giving each prime ~30 bits of denominator coverage.

Per-prime cost: O(N³) int64 ops in numpy. At N=1458, that's ~3·10^9 int64 ops. With numpy's vectorized inner-loop (each row-subtract is a single ~N-element vectorized op), the inner loop becomes ~N² = 2·10^6 vectorized ops, each ~1458 elements. Numpy throughput on int64 add/multiply: ~10^9 ops/sec. So per-prime wall time ≈ 3·10^9 / 10^9 = ~3 sec, plus mod ops. Realistic estimate: 5-15 sec per prime.

### 3. CRT recovery

For each component i of π:
- Have residues (r_1, r_2, ..., r_K) with moduli (p_1, p_2, ..., p_K).
- CRT to integer π_int_i mod P, where P = Π p_k.
- Apply rational reconstruction: find smallest n, d such that π_int_i ≡ n/d (mod P) with |n|, d ≤ √(P/2). The standard half-extended-GCD algorithm.

If P > 2 · max(|num|, den)^2, the reconstruction is unique. For π_k at k=7, num/den sizes are bounded by 2^M = 2^1458 — but π is normalized to sum 1 so denominators are typically MUCH smaller than 2^M. Empirically from cache, ε_6 denominator is ~200 digits ≈ 2^664. So we need P ≈ 2^1330 or larger, i.e., **~50 primes of 30 bits each**, or **~25 primes of 60 bits each**. We'll use ~30 primes of ~30 bits to keep numpy int64 multiply safe (since p^2 needs to be < 2^63).

Cache plan: extend the prime set adaptively. After each new prime, attempt rational reconstruction on a few "sentinel" components. If reconstruction stabilizes (same n/d twice in a row), declare success and verify with a fresh prime. If not, add more primes.

### 4. Verification

After recovering π_k as rational, compute X_k = 3^k · Σ π_k(r)^2 and ε_k = X_k − X_{k-1} − 7/15. For k=2..6, compare against the cache. All must match exactly.

## Why this beats the original

Per-op cost:
- Original: 5–100 μs per Fraction multiply (growing with denominator size). N³ × this = 8.5+ hours.
- New: ~1 ns per int64 op (numpy-vectorized). N³ × this = ~3 sec per prime × 30 primes ≈ 1.5 minutes for the modular solves themselves, plus overhead.

**Total target: under 10 minutes at k=7.** Probably significantly less.

## Risk register

- **R1 (CRT failure on coprime denominators):** if some prime p divides the true denominator of π's components, those components reconstruct wrong mod p. Mitigation: use 30 primes >> minimum; track per-component reconstruction agreement; drop "bad" primes if detected. The probability that ALL 30 primes (drawn from 10^8-10^9 range) coincidentally divide π's denominator is essentially zero for any specific number.

- **R2 (Rational reconstruction failure):** if the true |num|·den exceeds √(P/2), reconstruction is unbounded. Mitigation: use 30 primes ≈ P ≈ 2^900, well above the ε_6-derived denominator-size estimate of 2^664. If a component fails, add more primes.

- **R3 (Modular system rank deficiency):** A_p might be singular over F_p even if A over Q is nonsingular (this happens iff p divides the determinant). Mitigation: detect during elimination (pivot column has no nonzero), skip this prime, log it.

- **R4 (Last-row-replacement breaks rank):** the constraint Σπ=1 substituted for the redundant K^T−I equation. Empirically validated at k=1..6 in the original solver. Should remain valid mod p for almost all p.

- **R5 (Numpy int64 overflow):** if p > √(2^63 − 1) ≈ 3.04·10^9, then `factor * row` can overflow before the mod. We use primes < 2·10^9 to be safe.

- **R6 (K_int magnitude at row-build):** the integer weights 2^(M−r_v) reach up to 2^M ≈ 10^440 digits. We build them as Python ints, then reduce mod p ONCE per prime. The mod cost is O(M/64) per entry × N entries per row × N rows per prime × 30 primes ≈ acceptable.

## Verification protocol (Phase 4)

Per the brief, run k=2..6 first and gate on exact match with cached ε_n. If ANY of those six fail, stop. Only on full pass run k=7. Then write ε_7 to the v2 JSON cache and document.

## Implementation: see `result_77_7_v2.py`
