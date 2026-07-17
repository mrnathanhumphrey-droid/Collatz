# R77.7 V2 — Phase 5: Wall-time and memory profile

**Date:** 2026-05-12. Wilson. Pre-run estimates; main thread will run and populate.

## Per-prime solve cost model

The dominant per-prime cost is dense Gauss elimination over F_p on an N×N matrix:

- N rows × N elimination columns × O(N) per row-subtract (numpy vectorized)
- Constant: ~3-5 μs per row-subtract in numpy (allocation + int64 op + mod)
- Total: N² · (constant) ≈ N² × 5 μs per prime

| k | N    | est per-prime (sec) | est primes | est total wall |
|---|------|----------------------|------------|----------------|
| 1 | 2    | < 0.01               | 6-10       | < 1s           |
| 2 | 6    | < 0.01               | 6-10       | < 1s           |
| 3 | 18   | ~0.01                | 8-12       | 1s             |
| 4 | 54   | ~0.04                | 12-16      | 1s             |
| 5 | 162  | ~0.4                 | 24-30      | ~10s           |
| 6 | 486  | ~3                   | 50-70      | 2-4 min        |
| 7 | 1458 | **14.4 sec measured**  | 150-200  | **~36 min extrapolated** |

Notes:
- N² × 5 μs is the inner-loop floor; pivot search and pow2_mod build add a few %.
- At k=7, N²=2.13e6 numpy ops × 5 μs = ~10s lower bound per prime; allocations
  push it to ~20-40 sec realistic. **Actual measured: 14.4 sec, near the floor.**

## Measured per-prime wall times (sub-agent run 2026-05-12)

| k | N    | actual per-prime | n_primes (verified) | total wall |
|---|------|------------------|---------------------|------------|
| 1 | 2    | < 1 ms           | 7                   | < 0.01 s   |
| 2 | 6    | < 1 ms           | 7                   | < 0.01 s   |
| 3 | 18   | ~1 ms            | 9                   | 0.01 s     |
| 4 | 54   | ~8 ms            | 13                  | 0.11 s     |
| 5 | 162  | ~80 ms           | 25                  | 2.04 s     |
| 6 | 486  | ~930 ms          | 51                  | 47.35 s    |
| 7 | 1458 | 14.4 s (measured one prime) | ~150 (projected) | ~36 min |

Per-prime cost ratio k=6 → k=7: 14.4 / 0.93 ≈ 15.5×. Theoretical N³ ratio
(1458/486)³ = 27×. Measured ratio is BETTER than theoretical because numpy's
vectorized inner-row-subtract amortizes Python overhead better at larger N.

## Memory profile

K_exp build at k=7:
- N=1458 rows × N=1458 entries per row (dense) × 1 exponent each = 2.1M exponents
- Each Python int (small) ≈ 28 bytes → ~60 MB for K_exp
- After build, this stays resident.

Per-prime numpy A_p:
- N × (N+1) int64 = 1458 × 1459 × 8 bytes ≈ 17 MB
- Freed after each prime, so peak memory ≈ K_exp + 1 A_p ≈ 80 MB.

CRT residue storage:
- 150 primes × N int64 vectors = 150 × 1458 × 8 bytes ≈ 1.8 MB. Negligible.

CRT-combined integer array x_arr:
- N=1458 Python ints, each up to log2(P) ≈ 4500 bits = 560 bytes each.
- Total ≈ 1 MB. Negligible.

Reconstructed Fractions:
- N=1458 Fractions, each ~ 600 bytes for num+den. ~1 MB. Negligible.

**Total peak memory: ~100 MB.** Well within budget.

## CRT prime count

The required CRT modulus product P must satisfy P > 2 · (max num · max den) for
unique rational reconstruction. From cache:
- ε_6 has 200-digit den → π_6 components likely have ~100-150 digit dens → need P > 2^1000
- ε_7 expected to have ~600-800 digit den → π_7 components likely have ~300-500 digit dens → need P > 2^3500

With 30-bit primes (≈ 10^9), this requires ~120 primes for k=7. The script starts
with 150 to provide margin.

## Comparison to original solver

| Metric | Original (Fraction-Gauss) | V2 (CRT-modular) |
|---|---|---|
| k=6 wall time | ~3-5 min | ~2-4 min (similar) |
| k=7 wall time (extrapolated) | 5-15 hr (killed at 8.5hr) | 1-2 hr |
| Speedup at k=7 | 1× | ~5-10× |
| Peak memory | ~2 GB (growing-denom Fractions) | ~100 MB |

The V2 advantage is concentrated at k=7 where:
- The original's per-op cost has blown up (denominator size).
- V2's per-op cost is constant int64.

## Risk: actual wall time could be 2-3× the model

The model assumes pure numpy speed; in practice, Python overhead in pivot
search and `int(M[row, col])` conversions adds 50-100% slowdown for small N
but is amortized at N=1458 (the row-op dominates).

## Verification-phase wall-time budget

k=1..6 should finish in 3-6 minutes total. If it takes > 15 min, investigate
before launching k=7.

## Post-run notes

After main thread runs the script, update this doc with:
- Actual per-prime wall time at each k
- Actual total wall time at each k
- Actual prime count used
- Verification result (all pass / which failed)
- Final ε_7 value and digit count of numerator/denominator
