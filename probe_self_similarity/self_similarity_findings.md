# Self-similarity probe — π_∞ on Z_3 + entropy connection

**Date:** 2026-05-06.  Cached profinite π_k loaded from `probe_profinite/pi_{8,9,10,11}.npz` and `pi_infinity_cylinder_representation.npy` (k=12).  k=5,6,7 computed inline.

## Phase 1.1: T_j domain check

Brief proposed T_j(y) = (2^j y - 1)/3 (inverse Syracuse) as the IFS contraction. In Z_3 metric, that map has Jacobian |2^j|_3 / |3|_3 = 1 / (1/3) = 3 — **expanding**, not contracting. The IFS direction that converges to π_∞ on Z_3 is the FORWARD Syracuse step T_j^fwd(x) = (3x + 1)·2^{-j}, with 3-adic Jacobian 1/3 (uniform contraction over j).

With T_j^fwd, the candidate equation π_∞(A) = Σ_j p_j π_∞((T_j^fwd)^{-1}(A)) reduces to the Markov-chain stationarity equation π_∞ = π_∞ · K_∞ where K_∞ is the framework's transfer kernel (Tao chain). Since π_k is *defined* as the stationary of K_k, the relation holds tautologically at every k.

## Phase 2: residuals (sanity check)

Verifying ||π_k · K_k − π_k||_∞ at each cached k:

| k | n | ||π_k·K_k − π_k||_∞ | Z (= weight normalization) |
|---|---|---|---|
| 5 | 162 | 6.94e-18 | 1.000000000000000 |
| 6 | 486 | 3.47e-18 | 1.000000000000000 |
| 7 | 1458 | 3.47e-18 | 1.000000000000000 |
| 8 | 4374 | 6.07e-18 | 1.000000000000000 |
| 9 | 13122 | 3.04e-18 | 1.000000000000000 |
| 10 | 39366 | 8.67e-19 | 1.000000000000000 |

All residuals at machine precision (≤ 10^-15). The candidate self-similarity is **trivially satisfied at the IFS-as-stationarity level** because π_k is constructed as the stationary measure of K_k.

### Cylinder-set verification

For random cylinders A = {x ≡ r mod 3^{k_a}} at k_a ≤ k_eval−1, compare LHS = π_k(A) to RHS = (π_k · K_k)(A):

| k_eval | n_tests | max |LHS − RHS| | max rel diff |
|---|---|---|---|
| 7 | 32 | 5.55e-17 | 2.27e-16 |
| 8 | 36 | 1.11e-16 | 4.09e-16 |
| 9 | 35 | 1.11e-16 | 9.84e-16 |

## Phase 4: entropy and connection to ρ_slow

Shannon entropy H(π_k) = -Σ π_k(r) log π_k(r) for k = 5..12, computed from cached π_k arrays. Per-level increment ΔH_k := H(π_{k+1}) − H(π_k) approaches log 3 = 1.0986 from below. Deficit Δ_k := log 3 − ΔH_k.

| k | H(π_k) | n_k | H_uniform | D_KL = H_uni - H | ΔH_k | Δ_k = log3 - ΔH_k |
|---|---|---|---|---|---|---|
| 5 | 4.607200 | 162 | 5.087596 | 0.480396 | 1.049423 | 4.918976e-02 |
| 6 | 5.656623 | 486 | 6.186209 | 0.529586 | 1.058474 | 4.013817e-02 |
| 7 | 6.715097 | 1458 | 7.284821 | 0.569724 | 1.065094 | 3.351783e-02 |
| 8 | 7.780191 | 4374 | 8.383433 | 0.603242 | 1.070221 | 2.839140e-02 |
| 9 | 8.850412 | 13122 | 9.482045 | 0.631633 | 1.074259 | 2.435307e-02 |
| 10 | 9.924671 | 39366 | 10.580658 | 0.655986 | 1.077562 | 2.104983e-02 |
| 11 | 11.002234 | 118098 | 11.679270 | 0.677036 | 1.080293 | 1.831918e-02 |
| 12 | 12.082527 | 354294 | 12.777882 | 0.695355 | — | — |

### Geometric decay fit

OLS fit log Δ_k = a + b·k over k = 5..11:

- intercept a = -2.195373
- slope b = log ρ = -0.163339
- **decay rate ρ_Δ = exp(b) = 0.849303**
- R² = 0.995597

Comparison to ρ_slow ≈ 0.826934 (order-3 recurrence root from `result_renormalization_recurrence_fits.csv`, fitted on ε_2..ε_10 in 2026-05-05; **note: this fit is now known to be window-unstable** — see `probe_epsilon_13/epsilon_13_findings.md`. The L¹/TV inverse-limit rate ρ ≈ 0.834 from `probe_profinite` is the more reliable identification of ρ_slow at the same scale):

- |ρ_Δ − ρ_slow| = 0.022369  (vs recurrence-fit value 0.826934)
- ρ_Δ / ρ_slow = 1.027050

## Verdict

**Phase 2 outcome:** the candidate self-similarity holds at machine precision, but this is **tautological** — the equation π = Σ_j p_j (T_j)_* π is exactly the stationarity equation π = π·K_k, and π_k is defined as the stationary. So π_∞ IS a self-similar measure under the IFS (T_1^fwd, T_2^fwd, ...) with weights (1/2, 1/4, 1/8, ...) — but this is a structural fact about how the framework was set up, not new empirical evidence for the Eberhard-Varjú class.

**The non-trivial Eberhard-Varjú-class question** is whether π_∞ has additional regularity properties (absolute continuity, dimension, etc.) implied by the entropy/contraction theory of self-similar measures. That requires moment computations, dimension estimation, or scaling exponent measurements that are out of scope for this probe.

**Phase 4 outcome — partial entropy connection (Outcome B/C, ambiguous):** the entropy
deficit Δ_k = log 3 − [H(π_{k+1}) − H(π_k)] decays geometrically with R² = 0.996
at empirical rate **ρ_Δ ≈ 0.8493**. The order-3-recurrence ρ_slow ≈ 0.8269
sits 0.022 below it (2.7% relative gap). Same scale (~0.83), not exact match.

Two diagnostic concerns about treating this as identifying the same mode:

1. **The per-step ratio Δ_{k+1}/Δ_k drifts monotonically upward** across the
   6 transitions: 0.816, 0.835, 0.847, 0.858, 0.864, 0.870. The geometric-fit
   rate ρ_Δ ≈ 0.85 is the average of these; the *terminal* rate at k=11→12
   is 0.870, and the trajectory suggests it may continue rising rather than
   approach 0.827. If true, ρ_Δ → some asymptote distinct from ρ_slow.

2. **Δ_k is monotone-positive (entropy deficit always positive and decreasing)
   while ε_k oscillates** with sign flips (most recently at k=9→10). Two
   sequences sharing an *approximate* decay rate at k=5..11 but with
   structurally different time profiles (monotone vs oscillatory) suggests
   they may track different underlying modes that happen to have similar
   magnitudes at the tested k.

Compatible interpretations, in order of decreasing optimism:

- **(B) Same mode, different observables.** ρ_Δ converges to ρ_slow as k→∞;
  the rising-ratio trend is a finite-k correction. The order-3 recurrence
  captures the oscillatory aspect of the same rate that drives the monotone
  entropy-deficit decay. If true, this is a substantive structural connection.
- **(C) Different modes, similar magnitudes.** ρ_Δ → some asymptote ≠ ρ_slow
  but in the same ballpark. They coincide approximately by accident of
  finite-k scale.
- **(D) Spurious match.** ρ_Δ from a 7-point geometric fit; ρ_slow from a
  10-point order-3 recurrence with R² = 0.797. Both estimators are noisy
  enough that a 2.7% gap doesn't distinguish "same mode" from "different
  modes that happen to be close in magnitude".

**Higher-k extension (k=13, 14) decides the question — Outcome C.** Built
sparse K_k truncated at v_max = 60 (tail contribution below 2⁻⁶⁰ ≈ 10⁻¹⁸,
below double precision), power-iterated to ‖δ‖₁ < 1e-15:

| k | H(π_k) | ΔH_k = H(π_{k+1}) − H(π_k) | Δ_k = log 3 − ΔH_k | ratio Δ_{k+1}/Δ_k |
|---|---|---|---|---|
| 12 | 12.082527 | 1.0825877 | 1.6025e-2 | 0.8747 |
| 13 | 13.165115 | 1.0845271 | 1.4085e-2 | 0.8790 |
| 14 | 14.249642 | — | — | — |

The full ratio sequence k=5→6 through k=12→13 (eight transitions):
**0.8160, 0.8351, 0.8470, 0.8578, 0.8643, 0.8703, 0.8747, 0.8790** —
strictly monotone rising. The trend continues; no saturation in the tested
range.

Refit results:
- Full fit (k=5..13, 9 points): ρ_Δ = **0.8569**, R² = 0.9947
- Late fit (k=9..13, 5 points): ρ_Δ = **0.8722**, R² = 0.9994

The late fit shows a *cleaner* geometric (R² = 0.9994 vs 0.9947) at a *larger*
rate (0.872 vs 0.857). The gap to ρ_slow ≈ 0.827 is **widening with more data**:

| fit window | ρ_Δ | gap from ρ_slow | R² |
|---|---|---|---|
| k=5..11 (original 7-point) | 0.8493 | +0.022 (2.7%) | 0.9956 |
| k=5..13 (full 9-point) | 0.8569 | +0.030 (3.6%) | 0.9947 |
| k=9..13 (late 5-point) | 0.8722 | +0.045 (5.4%) | 0.9994 |

**This is the signature of two different modes, not one mode tracked
imprecisely.** If ρ_Δ were the same as ρ_slow, the gap would close (or at
least stay flat) as more data was collected — the OLS estimator would
converge to the true rate. Instead, the gap *widens* and the late-fit
estimator is *better* than the full-fit. ρ_Δ has its own asymptotic limit,
distinct from ρ_slow, somewhere ≥ 0.88.

**Verdict: Outcome (C) confirmed.** The entropy-deficit decay rate and the
order-3-recurrence rate are **different structural quantities** that
happened to be in the same order of magnitude at small k. The 0.83 vs 0.85
proximity at k≤11 was a finite-k coincidence; with k=12, 13, 14 in hand
the rates demonstrably separate.

What ρ_Δ ≈ 0.87 *could* be is now an open question — candidates:
- 7/8 = 0.875 (suggestive but no derivation)
- A different eigenvalue/singular value of the framework's operators
- An asymptotic limit not yet reached at k=13

Resolving this is a separate open piece, not load-bearing for the original
ρ_slow ≈ 0.83 question.

What the probe **does** establish:
- π_∞ satisfies the Markov-IFS stationarity equation by construction (Phase 2,
  not new content — tautological).
- Per-level entropy of π_k approaches log 3 from below with a clean monotone
  geometric decay at rate ~0.85 over k=5..11 (R² = 0.996).
- This rate sits within 3% of ρ_slow ≈ 0.83, but the connection isn't pinned
  down at the precision of these probes.

**Other entropy candidates from the brief, ruled out:**
- exp(−H(π_∞)) ≈ ρ_slow → predicts H_per_step ≈ 0.186, way smaller than the
  empirical ΔH_k ≈ 1.08. **No match.**
- ρ_slow = ΔH/log 3 → predicts ρ_slow ≈ 0.984 (matches the now-walked-back
  two-mode rate at k=2..10, not the current ρ_slow ≈ 0.83). **No match for
  current ρ_slow.**
- ρ_slow = 1 − H_per_step / log 3 → predicts ρ_slow ≈ 0.017. **No match.**

The only candidate that lands in the right neighborhood is the entropy-deficit
decay rate itself ≈ 0.85 vs ρ_slow ≈ 0.83 — discussed above as ambiguous.

## Files

- `test_results_geometric_V.csv` — cylinder-set self-similarity tests
- `entropy_computation.csv` — per-k entropy, KL deficit, ΔH_k, Δ_k
- `self_similarity_findings.md` — this writeup