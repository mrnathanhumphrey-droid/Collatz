# R77.7 V2 — DISPOSITION

**Date:** 2026-05-12. Wilson, sub-agent task: design and implement a new exact-rational ε_7 solver to replace the original Fraction-Gauss attempt killed at 8.5 hr.

---

## DISPOSITION: H_CRT_SOLVER_WORKS + DESIGN_COMPLETE_BUT_NOT_RUN (k=1..6 verification PASSED; k=7 awaits main thread)

Sub-agent successfully ran short Bash calls but was denied permission for the
~36-min k=7 run. Phase 4 verification (k=1..6 vs cache) **passed end-to-end**:
all six ε_k reproduced exactly as Fractions from the cache. The solver is
complete and validated; k=7 just needs to be launched on main thread.

| k | wall time | primes used | match cache (Fraction equality) |
|---|-----------|-------------|---------------------------------|
| 1 | 0.00s     | 7           | True                            |
| 2 | 0.00s     | 7           | True                            |
| 3 | 0.01s     | 9           | True                            |
| 4 | 0.11s     | 13          | True                            |
| 5 | 2.04s     | 25          | True                            |
| 6 | 47.35s    | 51          | True                            |

Single-prime k=7 measurement: **14.4 sec/prime** (N=1458). With ~150 primes,
projected k=7 wall ≈ **36 min**. Compared to the original solver's 8.5+ hr
(killed mid-run), V2 is ~14× faster.

**To run:**
```
python C:/Collatz/result_77_7_v2.py
```

Output: `C:/Collatz/experiments_output/result_77_7_eps_exact_through_k7_v2.json`

Optionally tee stdout to `C:/Collatz/R77_7_V2_RUN_LOG.txt`.

---

## Summary of the design

The original solver's bottleneck is dense Gauss-elimination over Q via
`fractions.Fraction`. As elimination proceeds, denominator size grows linearly
in the elimination depth, making per-Fraction-op cost climb from microseconds
to hundreds of microseconds. The O(N³) total operation count multiplied by
this growing per-op cost is the 8.5 hr wall.

The V2 solver eliminates the denominator-growth factor by working over a
sequence of finite fields F_p with `numpy int64`. Each per-prime solve is
fixed-cost O(N³) with vectorized int64 arithmetic. After ~150 primes (for
k=7), the rational result is reconstructed component-wise via CRT + half-
extended Euclidean.

The Markov chain is constructed once over the integers (as
exponent-lists, since each entry is a sum of powers of 2). For each prime, the
entries reduce in O(M) per cell where M=N at q=3. The integer matrix is
A = K_int^T − (2^M − 1) I (last row replaced with all-ones, RHS = e_N).

Adaptive prime-batching: if rational reconstruction fails OR if a held-out
"witness" prime disagrees with the reconstructed Fractions, the solver adds
10 more primes and retries. Max cap: 200 primes per k.

## Pre-run expectations

### Phase 4 verification (k=1..6 vs cache) — RESULT: ALL PASS

All six match the cache exactly as Fractions (sub-agent ran 2026-05-12):
- ε_1 = 1/5 ✓
- ε_2 = 1/105 ✓
- ε_3 = -5191/1019445 ✓
- ε_4 = -11346676448406637/4627031617157687115 ✓
- ε_5 = (60-digit num) / (65-digit den), -1.151747e-3 ✓
- ε_6 = (197-digit num) / (200-digit den), -4.979057e-4 ✓

### Phase 5 expectation (k=7)

Per STATE.md and project's prior numerical (power-iteration) measurement, the
exact ε_7 should agree at the leading ~6 decimal digits with **ε_7 ≈ -1.18e-3**,
and the ratio **|ε_7| · 2^7** should jump to ≈ 0.15 if the prior measurement of
"~4.7× per step" jump holds. (For comparison: |ε_5|·2^5 ≈ 0.0369, |ε_6|·2^6 ≈
0.0319.)

The numerator of ε_7 should be ~600-800 decimal digits and the denominator
similarly sized, based on the cache's 3-4× growth rate per step in ε's
denominator.

### Wall-time expectation

Per `R77_7_V2_PROFILE.md`:
- k=1..6 verification total: 3-6 min
- k=7: 1-2 hr
- Total: under 2.5 hr (under the 4-hr threshold, well under the 8.5-hr original)

## Adversarial checks per pre-registration

- **A1 (project-internal Markov chain fidelity):** The V2 build_K_exponents
  produces the SAME object as the original's `build_markov_rational`, just
  in integer form (with common denominator 2^M − 1). The transformation is
  exact and the same target_of(r, r_v) formula is used.
- **A2 (CRT failure mode if prime divides denom):** Mitigated by using ~150
  primes and the witness check that detects this case and skips the bad
  prime.
- **A3 (Rational reconstruction failure):** Mitigated by adaptive prime
  batching that grows the modulus product P until reconstruction succeeds
  for all components.
- **A4 (Verification scope):** All six of ε_1..ε_6 must match cache before
  k=7 runs. The main loop halts if any fail.

## Deliverables

- `C:/Collatz/R77_7_V2_BOTTLENECK_ANALYSIS.md` — Phase 1
- `C:/Collatz/R77_7_V2_DESIGN.md` — Phase 2
- `C:/Collatz/result_77_7_v2.py` — Phase 3 (runnable from main thread)
- `C:/Collatz/R77_7_V2_VERIFICATION.md` — Phase 4 protocol
- `C:/Collatz/R77_7_V2_PROFILE.md` — Phase 5 estimates
- `C:/Collatz/R77_7_V2_DISPOSITION.md` — this file
- `C:/Collatz/experiments_output/result_77_7_eps_exact_through_k7_v2.json` —
  populated only after main-thread run

## Risk register (open items)

- The wall-time model is an estimate; actual could be 1.5-3× slower if numpy
  overhead dominates more than expected. The internal Python loops in
  `gauss_solve_mod_p` (pivot search, factor extraction) are not vectorized.
- If actual k=6 wall time > 15 min, halt before launching k=7 and tune (use
  scipy `linalg.lu` over F_p via the LU decomposition over Z then mod p?
  Or implement Wiedemann's algorithm? Or use the structure that K = 1/3 ·
  (uniform on cosets) + correction to factor the matrix?).
- If a prime "early" in the sequence divides the determinant, it triggers
  the singular branch; current code skips it and chooses another. Should
  not be a hot path but is logged.

## Hand-off to main thread

The single command to run is:

```
python C:/Collatz/result_77_7_v2.py
```

Run from any directory; the script uses absolute paths. Tee or redirect
stdout if you want a log file:

```
python C:/Collatz/result_77_7_v2.py 2>&1 | tee C:/Collatz/R77_7_V2_RUN_LOG.txt
```

(On Windows PowerShell, use `Tee-Object`:)

```
python C:/Collatz/result_77_7_v2.py 2>&1 | Tee-Object -FilePath C:/Collatz/R77_7_V2_RUN_LOG.txt
```

The script saves the eps cache incrementally after each k completes (so a
mid-run crash retains all completed levels). Expected duration:
- k=1..6 verification: ~50 sec total (dominated by k=6 at ~47 sec)
- k=7: ~36 min projected, ~150 primes at ~14.4 sec each
- Total wall: under 40 min

## What success unlocks

Per Wilson's brief: ε_7 (exact rational) unlocks the [3/3] diagonal Padé
approximant to E(z) = Σ ε_n z^n. With six exact coefficients (ε_2..ε_7),
the m+n=5 Padé approximants become computable, including [3/2] near-diagonal
and [3/3] full diagonal. These sharpen the branch-cut-order discrimination at
z=2 and constrain the subleading singularity location. Per the DELTA and
PADE_EXTENSION dispositions of 2026-05-12, ε_7 is "the gating data" for the
Tauberian framework arc to advance.

The ratio |ε_7|·2^7 also discriminates whether n=7 has entered the
slow-mode regime predicted by the δ_n diagnostic. STATE.md cites a
numerical-only prior estimate of |ε_7|·2^7 with a ~4.7× jump from |ε_6|·2^6;
the exact-rational value will confirm or refute this.
