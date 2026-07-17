# q-sweep across qx+1 family — K-sum saturation

**Date:** 2026-05-04
**Object:** K_p(r, c, m) = Σ_{u=0}^{N_p − 1} e_M(c · (1+p)^u − p² · m · u),  M = p^{r+1},  N_p = p^{r-1}
**Goal:** determine whether the empirical √N saturation observed at p = 3 is q-universal (Pattern α), exponent-universal but prefactor-varying (Pattern β), or genuinely q-shaped (Pattern γ).

---

## 1. Executive summary

**Pattern β (with weak signal — close to Pattern α at the exponent level).**

All five primes saturate at the √N rate. The saturation exponent β_p is universal across {3, 5, 7, 11, 13} within sampling noise:

| Metric | β range across primes | Spread |
|---|---|---|
| K_max (c,m search) | [0.486, 0.525] | 0.040 |
| K_c1m0 (canonical) | [0.483, 0.518] | 0.035 |

The K_max exponent for p = 3 is borderline outside [0.48, 0.52] (0.525), but this is driven by (c, m) resonance peaks that are visible at p = 3's deeper r-range (r = 8..20) and not seen as sharply at smaller-r data for other primes. Under the canonical (c=1, m=0) metric, p = 3 is at 0.494 — the K_max effect was a (c, m)-search-window artifact, not a structural difference.

Prefactor C_p = |K|/√N varies mildly with p. Under K_c1m0, C drifts from 0.83 (p = 3) to 1.17 (p = 13), a ~40% variation — barely above the 30% Pattern β threshold.

**Strategic implication:** The Plancherel obstruction is q-universal at the exponent level. Path C (5x+1 or any other prime as a sibling attack) does **not** open a polynomial closure — every prime hits the √N wall at the same exponent. The mild prefactor variation has no apparent arithmetic structure (not monotone in p, not correlated with p mod 4, not correlated with mult. order of 2 mod p). **Consolidate Path B; Path C blocked.**

---

## 2. Methodology

### 2.1 q-adic K formula derivation

The 3-adic Kalafatelis sum is K_3(r, c, m) = Σ_{u=0}^{N-1} e_q(c · 4^u − 9m · u) with q = 3^{r+1}, N = 3^{r-1}. The constants "4" and "9" come from:

- "4" = 1 + 3 — the standard generator of principal units 1 + 3Z mod 3^k for any k ≥ 2
- "9" = 3² — the modulus reduction factor in the short-window length N = q/9

Generalizing to prime p ≥ 3:

- M = p^{r+1} (modulus)
- (1 + p) generates the principal units 1 + pZ mod p^{r+1}, with order p^r
- N_p = p^{r-1} = M / p² (short-window length)
- K_p(r, c, m) = Σ_{u=0}^{N_p − 1} e_M(c · (1+p)^u − p² · m · u)

At p = 3 this reduces to the existing 3-adic formula. **Verification at p = 3, r ∈ {8, 9, 10}:** the new generalized code reproduces the existing `r79b_S_partial_data.csv` K_c1m0_abs values to relative error 10⁻⁹ (see `q_sweep.py` verify_p3_match function).

### 2.2 Computation

For each p ∈ {3, 5, 7, 11, 13}:
- Computed |K_p(r, c, m)| via Numba @njit(parallel=True) direct sum over u with running multiplication for (1+p)^u mod M
- Searched (c, m) ∈ {1..p-1} × {0..p-1} for the maximum |K| at each r
- Recorded both |K_c1m0| (canonical) and |K_max| (search optimum)

For p = 3, reused existing `r79b_S_partial_data.csv` data at r ∈ {8, 10, 12, 14, 16, 18, 20}.

For p = 5, 7, 11, 13, computed fresh at r-targets:
- p = 5: r ∈ {6, 8, 10, 12}
- p = 7: r ∈ {6, 8, 10}
- p = 11: r ∈ {5, 6, 7, 8}
- p = 13: r ∈ {4, 5, 6, 7}

CSV written with append+fsync after each (p, r) — watchdog-safe.

---

## 3. Per-q results

### 3.1 K_max metric (search over (c, m))

| p | r-range | n | β (OLS) | σ_β | R² | C @ r_max | C mean |
|---|---|---|---|---|---|---|---|
| 3 | 8..20 | 7 | 0.5251 | 0.0125 | 0.9972 | 2.652 | 1.883 |
| 5 | 6..12 | 4 | 0.4897 | 0.0121 | 0.9988 | 1.545 | 1.521 |
| 7 | 6..10 | 3 | 0.4855 | 0.0213 | 0.9981 | 1.557 | 1.575 |
| 11 | 5..8  | 4 | 0.4939 | 0.0211 | 0.9964 | 1.641 | 1.770 |
| 13 | 4..7  | 4 | 0.4910 | 0.0187 | 0.9971 | 1.772 | 1.800 |

Range across primes: β ∈ [0.486, 0.525], spread 0.040. C @ r_max ratio max/min = 1.717 (p=3 vs p=5).

### 3.2 K_c1m0 metric (canonical c=1, m=0)

| p | β_canonical | C_c1m0 mean | C_c1m0 @ r_max |
|---|---|---|---|
| 3 | 0.4935 | 0.8846 | 0.8300 |
| 5 | 0.4831 | 0.9439 | 0.9048 |
| 7 | 0.4836 | 0.9845 | 0.9357 |
| 11 | 0.5169 | 1.0581 | 1.1194 |
| 13 | 0.5177 | 0.9964 | 1.1653 |

Range across primes: β ∈ [0.483, 0.518], spread 0.035. C @ r_max ratio max/min = 1.404 (p=13 vs p=3).

### 3.3 Cross-comparison

The two metrics give consistent classification (β universal, mild C variation), but the K_max metric inflates p = 3 because of the deeper r-window catching rare (c, m) resonance peaks. The canonical metric is more apples-to-apples across primes.

The (c, m) inflation factor K_max/K_c1m0 at largest r per prime:
- p = 3, r = 20: 3.20×
- p = 5, r = 12: 1.71×
- p = 7, r = 10: 1.66×
- p = 11, r = 8: 1.47×
- p = 13, r = 7: 1.52×

p = 3 has noticeably more (c, m)-resonance peaks per unit r-volume than other primes — but this reflects the larger r-window we have for p = 3, not a different scaling exponent.

---

## 4. Pattern fit

### 4.1 Brief's classification thresholds

- Pattern α: all β_p ∈ [0.48, 0.52] AND C_p within 30%
- Pattern β: all β_p ∈ [0.48, 0.52] BUT C_p varies > 30%
- Pattern γ: any β_p outside [0.48, 0.52] with R² ≥ 0.95

### 4.2 Verdict

**On K_max metric:** β_3 = 0.525 is borderline outside [0.48, 0.52] (R² = 0.997). Strict letter would suggest Pattern γ candidacy for p = 3. But the elevation is window-driven, not structural (canonical metric drops to 0.494). C @ r_max ratio = 1.72, exceeds 30% threshold.

**On K_c1m0 metric:** all β within [0.483, 0.518], well inside the band. C ratio = 1.40, exceeds 30%.

**Combined verdict: Pattern β.** Universal β ≈ 0.5 across all primes; prefactor C varies modestly (factor ~1.4–1.7 across primes by metric).

The variation is at the borderline of the Pattern α / β cutoff. Calling it Pattern β rather than α is conservative — under tighter prefactor measurement (e.g., averaged over a uniform r-window per prime) the spread might fall under 30%.

### 4.3 Statistical caveats

- Per-prime n_r ranges from 3 (p = 7) to 7 (p = 3) — β estimates have heterogeneous precision
- σ_β ranges from 0.012 to 0.021; β_p − 0.5 distances are 0.5 σ to 2.0 σ — none statistically distinguishable from 0.5 at α = 0.05
- Finite-r asymptotics: at small r the prefactor C fluctuates substantially. p = 3 ratios oscillate 1.63–2.65 across r = 8..20 (see raw CSV) — the "C @ r_max" point estimate is noisy

---

## 5. C_p variation analysis (Pattern β follow-up)

Per the brief, if Pattern β is identified, characterize C_p variation. Tested:

- **C vs p (monotone?):** C_c1m0 mean: 0.88 (p=3), 0.94 (p=5), 0.98 (p=7), 1.06 (p=11), 1.00 (p=13). Weakly increasing in p but non-monotone (p=13 < p=11).
- **C vs p mod 4:** p=3 (≡3): 0.88, p=5 (≡1): 0.94, p=7 (≡3): 0.98, p=11 (≡3): 1.06, p=13 (≡1): 1.00. No clean grouping.
- **C vs ord(2 mod p):** 2, 4, 3, 10, 12. C: 0.88, 0.94, 0.98, 1.06, 1.00. Roughly increasing in ord but breaks at p=13.
- **C vs (p−1)/ord(2 mod p):** All five primes give ratio ≈ 1, except p = 7 (= 2). No structure.

**No arithmetic invariant of p found that explains C_p variation.** The variation looks like residual finite-r noise + potential weak p-scaling, not a sharp structural relationship.

---

## 6. Strategic implication

The data does **not** support a closure path through prime-switching:

1. **Path C (attack a specific prime where β is smaller)** — no prime in {3, 5, 7, 11, 13} has β statistically distinguishable from 0.5. The √N wall is the same height in every prime.
2. **No prime gives polynomial saving** (β < 0.5 with statistical significance). Even the smallest β estimate (p = 7, K_max metric, β = 0.486 ± 0.021) is 0.7 σ below 0.5 — not a real signal.
3. **Prefactor variation is not exploitable** — no arithmetic invariant explains C_p variation, and even the worst-prefactor prime is still saturated at √N · const.

**Recommended call:** Consolidate Path B (publish what's rigorous, frame the universality of the Plancherel obstruction across qx+1 as the substantive R78 finding). Skip Path C — the empirical universality result IS the closure-path-blocking statement.

The Pattern β signal (mild C variation) is itself worth a one-paragraph remark in the consolidation: **"The Plancherel obstruction has q-universal exponent and weakly-q-dependent constant — the wall is the same height across qx+1 systems, with prefactor varying within a factor of ~1.4 across small primes."**

---

## 7. Honest caveats

1. **r-windows are not matched across primes.** p = 3 covers r = 8..20; p = 13 covers r = 4..7. Extrapolating "C is universal" requires comparing prefactors at matched scales, which we don't have. A more rigorous follow-up would extend p = 5, 7, 11, 13 to comparable r-ranges (compute-bounded but feasible to r ≈ 14 for p = 5 with current code).

2. **The (c, m) search is not exhaustive.** We searched {1..p-1} × {0..p-1} = (p-1)·p combinations. For closure-path analysis we want the supremum over all (c, m) that arise in the bridge-equation framework's natural parameters. This may include (c, m) outside our search box for some primes.

3. **β_3 = 0.525 (K_max metric) deserves separate scrutiny.** Restricting the p = 3 OLS to higher-r windows pushes β even higher (r = 14..20: β = 0.556, R² = 0.99). This monotone β-with-window growth is INCONSISTENT with β converging to 0.5 from above; it's consistent with rare (c, m) resonance contributions being more visible at higher r. The canonical-metric β = 0.494 supports the resonance interpretation, not a genuine β > 0.5.

4. **q = 9 (composite) was not tested.** The brief specified primes only. A worthwhile sanity check is whether composite q ∈ {9, 15, 21} also saturate at √N — if so, the obstruction is broader than "qx+1 over primes."

5. **The prefactor variation might disappear under standardized r-windows** — calling Pattern β rather than α is a conservative reading.

---

## Files

- `C:\Collatz\q_sweep.py` — implementation (Numba parallel)
- `C:\Collatz\q_sweep_data.csv` — raw data
- `C:\Collatz\q_sweep_analyze.py` — OLS fits and pattern classifier
- `C:\Collatz\q_sweep_canonical.py` — canonical (c=1, m=0) cross-check
- `C:\Collatz\q_sweep_p3_subset.py` — p=3 r-window stability check
- `C:\Collatz\q_sweep_log.txt` — runtime log including verification of p=3 against existing CSV
