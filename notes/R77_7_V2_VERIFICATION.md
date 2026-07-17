# R77.7 V2 — Phase 4: Verification protocol & ACTUAL RESULTS

**Date:** 2026-05-12. Wilson.

## Status: VERIFICATION PASSED (k=1..6 against cache)

Sub-agent ran the verification phase end-to-end via short Bash calls (k=7 was
denied due to compute-time policy). Per-k results below. The full k=7 run
must be launched on main thread:

```
python C:/Collatz/result_77_7_v2.py
```

(or with stdout redirected:)

```
python C:/Collatz/result_77_7_v2.py > C:/Collatz/R77_7_V2_RUN_LOG.txt 2>&1
```

## What the script does (verification phase)

For k = 1, 2, 3, 4, 5, 6:

1. Build K_exp at level k (integer matrix as exponent-lists).
2. For an initial batch of primes (sizes from `initial_primes_by_k`):
   - Each prime: build A_p, run Gauss-elim mod p, get π mod p.
3. CRT-combine π mod p_1, ..., π mod p_k.
4. Rational-reconstruct each component to a Fraction.
5. Verify sum(π_k) == 1.
6. Add one extra "witness" prime, solve mod that, verify the reconstructed
   Fractions match the witness's residue.
7. If any step fails: add 10 more primes, retry. Cap at 200 primes.

Then compute ε_k = (3^k · Σ π_k²) − X_{k−1} − 7/15 and compare against the
cached value from `experiments_output/result_77_7_eps_exact_through_k7.json`.

## Pre-registered expected outputs (k=1..6)

The cache contains:

| k | ε_k (decimal)     | ε_k (digits num/den) |
|---|-------------------|----------------------|
| 1 | +2.000000e-01     | 1 / 1                |
| 2 | +9.523810e-03     | 1 / 3                |
| 3 | -5.091987e-03     | 4 / 7                |
| 4 | -2.452108e-03     | 17 / 19              |
| 5 | -1.151758e-03     | 60 / 65              |
| 6 | -4.978049e-04     | 197 / 200            |

The new solver MUST reproduce these exactly as Fractions. The verification
loop (`eps_dict[k] == cached_eps[k]`) checks Fraction equality, not float
approximation — so even a single bit off in the numerator or denominator
counts as FAIL.

## Decision rules from the verification

- **All six match** → proceed to Phase 5 (compute ε_7).
- **Any fail** → halt, save partial cache, write disposition with the specific
  failure mode. Do NOT run k=7 (it would be wasted compute).

## Likely failure modes if any k=2..6 fails

| Failure | Likely cause | Action |
|---|---|---|
| ε_k mismatch but ε_{k-1} OK | Insufficient primes at level k (denom > P/2 reconstruction bound) | Increase initial_primes_by_k[k] |
| ε_k denominator wildly different magnitude | Off-by-one in K_int construction or transpose mistake | Re-inspect build_K_exponents / solve_pi_mod_p |
| sum(π_k) ≠ 1 | Bug in last-row replacement or e_N RHS | Re-inspect solve_pi_mod_p |
| Rational reconstruction returns None | Bound too tight | Increase n_primes_initial |
| Witness mismatch but Fraction equality with cache | Witness prime divides true denominator (extremely unlikely with random primes) | Skip witness, accept reconstruction |

## Reading the run log

The script logs per-prime solve times. From these, the user can extrapolate:
- per-prime cost at k=6 → per-prime cost at k=7 by N³ ratio (× 27).
- total k=7 wall time = (per-prime at k=7) × (n_primes_at_k_7).

If k=6 per-prime is 1.5 sec and uses 50 primes → 75 sec total.
Then k=7 per-prime ≈ 40 sec, × 150 primes = 6000 sec ≈ 100 min ≈ 1.7 hr.
That's within the 2hr budget.

## Wall-time budget reality check

- If k=6 verification finishes in < 5 min total → likely all good, k=7 will be < 2hr.
- If k=6 verification takes 30+ min → k=7 extrapolation > 4hr; halt and reconsider.
- If k=2..5 finishes quickly but k=6 stalls at "solving prime X" forever → likely
  a Python-numpy interaction issue; halt and inspect.

## Decision: proceed to Phase 5?

**YES.** Verification PASSED at k=1..6. The actual run on 2026-05-12 gave:

| k | wall time | primes | match cache |
|---|-----------|--------|-------------|
| 1 | 0.00s     | 7      | True        |
| 2 | 0.00s     | 7      | True        |
| 3 | 0.01s     | 9      | True        |
| 4 | 0.11s     | 13     | True        |
| 5 | 2.04s     | 25     | True        |
| 6 | 47.35s    | 51     | True        |

ε_1..ε_6 reproduced exactly as Fractions from the cache. Decision: proceed
to Phase 5 (run k=7 on main thread).

Single-prime k=7 timing measured separately: **14.4 sec per prime** at
N=1458. With ~150 primes expected for full reconstruction, extrapolated
total k=7 wall ≈ **36 min**. Well under the 2-hour budget; ~14× faster than
the original solver's 8.5+ hr (killed).
