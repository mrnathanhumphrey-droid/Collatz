# R77.7 V2 — Phase 1: Original solver bottleneck analysis

**Date:** 2026-05-12. Wilson, follow-up to R77.7's 8.5hr kill.

## What the original solver does

`result_77_7_extend_to_k7.py` (preserved as reference) computes the level-k Markov chain stationary distribution π_k over Q at q=3 for k=1..7. The chain lives on the N=2·3^(k−1) coprime states in Z/3^k. The kernel K_k(r → r') is

    K(r, target_of(r, r_v)) += Fraction(1, 2^r_v) / Z_v

with `target = ((3*r + 1) * inv2^r_v) % N`, `M = 2·3^(k-1)`, `Z_v = (2^M − 1)/2^M`. At k=7, this gives a stochastic matrix on N=1458 coprime states with M=1458 nonzero geometric branches per row (one per r_v in 1..M).

Then ε_k = X_k − X_{k−1} − 7/15, where X_k = 3^k · Σ π_k(r)^2.

## How it solves the stationary system

`stationary_rational(K)` builds A = K^T − I, replaces the last row with the all-ones constraint, and runs **textbook Gauss elimination over Q** using `fractions.Fraction`. That is:

1. For each pivot column, find a nonzero pivot row.
2. Divide the pivot row by the pivot.
3. Subtract multiples of the pivot row from every other row to zero out that column.

The inner loop is a single `Fraction` multiply-subtract. Pure Python triple loop, no numpy.

## Why this hits a wall — three compounding factors

### Factor 1: O(N³) operation count

N=1458 at k=7. The pivot loop is N² rows × N columns per pivot column × N pivot columns = ~1.458^3 × 10^9 ≈ 3·10^9 Fraction operations. Even at 1 microsecond per op (pure-Python overhead), that's ~50 minutes; at 5–10 μs each (realistic for non-trivial Fractions), it's hours.

### Factor 2: Denominator blowup (the actual bottleneck)

`Fraction` stores `num` and `den` as arbitrary-precision `int`. After Gauss-elimination over Q, intermediate matrix entries have numerator/denominator sizes that can grow LINEARLY in the elimination depth. By the time you're partway through the 1458-step elimination, each Fraction entry can carry numerators of **hundreds of digits**. The size of the cached k=6 final result confirms this: the ε_6 denominator already has ~200 digits.

Each Fraction multiply at that size is O(d²) bit-operations (schoolbook) or O(d log d) (Karatsuba). The asymptotic operation cost is therefore not O(N³) — it's **O(N³ · d²)** where d grows with N. The 5–10 μs constant becomes 100+ μs by midway through.

This is why the run at 8.5hr was "still grinding at 96% CPU with no signs of stuck" — it wasn't stuck, the per-op cost was just steadily climbing.

### Factor 3: GCD reduction overhead

`Fraction` automatically reduces to lowest terms on every operation via `gcd`. This is O(d log d) per op and consumes a substantial fraction of each multiply/add. Disabling it would speed things up but produces unreduced denominators that bloat the next op even faster.

## Empirical timing at smaller k (extrapolation basis)

From the cache build the small-k times are reproducible. Per the R77.5 logs the times grow roughly:
- k=4 (162 states): O(seconds)
- k=5 (486 states): O(minute)
- k=6 (1458 wait — no, k=6 N = 2·3^5 = 486 states; correction below)
- k=7 (1458 states): 8.5+ hr killed mid-run

Correction on state counts: N_k = 2·3^(k-1). So:
- k=4: 54
- k=5: 162
- k=6: 486
- k=7: 1458

Cube scaling gives factor 27 from k=6 to k=7, multiplied by ~2-4× denominator-size factor. If k=6 takes ~3-5 minutes (matching prior reports), k=7 should land in the 5-10 hour ballpark just from the cube. The denominator growth pushes it further.

## Diagnosis summary

The original solver's bottleneck is **denominator blowup multiplied by O(N³) operation count**. Both are required ingredients:
- Pure O(N³) at constant-cost ops would finish in a small number of hours.
- Denominator blowup alone (in a Krylov method) wouldn't bite because you'd only do O(N·nnz) ops.
- The product is what kills it.

## What the new solver must do differently

To avoid both bottlenecks simultaneously, the new solver must either:

(a) **Eliminate denominator blowup** by working over a finite field F_p (entries are uint64 mod p, fixed cost ops), at the cost of needing CRT recovery to get back to Q. This is the **CRT + modular sparse solve** approach.

(b) **Eliminate the O(N³) factor** by exploiting structure (block-triangular kernel, sparse representation, iterative method) — but this only helps if the iteration converges in O(log N) steps, which is not guaranteed for a stochastic-matrix stationary problem with non-trivial mixing time.

(a) is the cleaner separation: it converts the denominator problem into "do the same O(N³) elimination but over F_p with int64 arithmetic" — each op is now a single machine instruction instead of a Fraction operation with bigint underneath. The CRT layer then assembles the answer. The cost trade is: we do `n_primes × O(N³ over F_p)` instead of `1 × O(N³ over Q with growing d)`. At N=1458 and ~30 primes, this should land 1-2 orders of magnitude faster.

This bottleneck analysis motivates the CRT approach pre-registered in the parent brief.
