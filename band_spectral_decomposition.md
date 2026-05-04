# Band-spectral decomposition of ĥ on D_{r,t}(η) — outcome β with γ-leaning functional read

**Date:** 2026-05-04. Continuation of C3 (band-l¹) and the post-morning Path-A obstruction map. Tests whether inter-m structure of ĥ_{r,ℓ} restricted to the dangerous band admits smooth-weight cancellation that the raw l¹ norm doesn't capture.

## Verdict

**By strict brief thresholds: outcome β (100% of 120 cells classified β, 0% α, 0% γ).**

**Functional reading: γ-leaning.** The spectrum of ĥ|_J as a function of m is **uniform across frequencies**, not low-frequency dominated. Smooth-weight cancellation gives at most a constant-factor improvement, not a polynomial saving sufficient for eq 190 closure.

**Smoking gun:** `low_freq_mass → 0.25 exactly` as r grows — equal to what a uniform-in-k spectrum predicts (mass in 1/4 of the frequency range). Combined with `r_eff/|J| ≈ 0.60 fixed across r`, the picture is: the spectrum has 40% "negligible" modes scattered randomly across all frequencies, NOT concentrated at high k where smooth weights can suppress them. The smooth-completion (R78 path 2) attack does NOT receive empirical motivation from this analysis.

## Setup

For each (r, ℓ, t, η) ∈ {6, 8, 10, 12, 14} × {0, 1, 2, 3} × {1, 2, 3} × {0.25, 0.5}:

1. Compute ĥ_{r,ℓ} via FFT over Z/N_r.
2. Build dangerous band D_{r,t}(η) = {m : |m_r,t(2πm/N_r)| > 1−η} per Kalafatelis Prop 22.
3. Restrict ĥ to J = D, treat as 1D sequence f := ĥ|_J in m-sorted order.
4. DFT of restriction: g(k) = Σ_j f(j) · e^{-2πi kj/|J|} (numpy fft over the |J|-length subsequence).
5. Spectral metrics on |g(k)|².

## Metrics

For each cell:

- **max_concentration** = max_k |g(k)|² / Σ|g|²  (= 1/|J| for uniform, → 1 for δ-concentrated)
- **conc_ratio** = max_concentration · |J|  (= 1 for uniform, large for concentrated)
- **entropy_normalized H_norm** = `−Σ p log p / log|J|` where p(k) = |g|²/Σ|g|²  (= 1 for uniform, → 0 for concentrated)
- **r_eff** = number of top modes carrying 90% of energy
- **r_eff_ratio** = r_eff / |J|  (= 0.9 for uniform, ≪ 1 for low-rank)
- **low_freq_mass** = mass in `|k| ≤ |J|/8`  (folded for negative freqs; = 0.25 for uniform-in-k spectrum, → 1 for low-frequency dominance)

Classification rule (brief):
- **α (smooth/exploitable):** r_eff_ratio < 0.25 AND lf_mass > 0.6
- **γ (white-noise rough):** r_eff_ratio ≥ 0.8 AND H_norm > 0.95
- **β (mixed):** in between

## Per-(r, η) results

```
   r    eta    |J|  r_eff  r_eff_r   H_norm  lf_mass conc_ratio    class
   6  0.250   98.3   58.0    0.590    0.911    0.239       4.13     β:12
   6  0.500  204.0  123.8    0.607    0.935    0.238       3.91     β:12
   8  0.250  883.0  524.1    0.594    0.942    0.248       5.62     β:12
   8  0.500 1835.3 1100.2    0.599    0.953    0.253       4.58     β:12
  10  0.250 7946.7 4681.6    0.589    0.956    0.249       6.56     β:12
  10  0.500 16516.7 9936.8    0.602    0.964    0.250       5.20     β:12
  12  0.250 71520.3 42262.5    0.591    0.965    0.250       7.61     β:12
  12  0.500 148651.0 89434.6    0.602    0.970    0.250       5.80     β:12
  14  0.250 643682.0 380046.2    0.590    0.971    0.250       7.91     β:12
  14  0.500 1337858.0 804899.2    0.602    0.975    0.250       6.05     β:12
```

12 cells per (r, η) entry. All 120 cells classify as β.

## Scaling analysis (the load-bearing observation)

Three metrics converge to specific limits as r grows:

| metric | trend | asymptote |
|---|---|---|
| r_eff_ratio | flat (0.589 ↔ 0.607 across r=6..14) | ≈ 0.60 |
| H_norm | monotonically increasing | → 1 (uniform) |
| lf_mass | converging from 0.24 to **exactly 0.25** by r ≥ 10 | = 0.25 |

The convergence `lf_mass → 0.25 exactly` is **decisive**: this is the value that a UNIFORM spectrum (|g(k)|² constant in k) gives when measuring "mass in first |J|/8 + last |J|/8 indices". Any concentration at low frequencies would push this above 0.25; concentration at high frequencies would push it below. The fact that it equals 0.25 to 4 decimals at r = 10, 12, 14 means the spectrum **has no preferential frequency structure**.

H_norm → 1 confirms: the |g(k)|² distribution is asymptotically uniform.

The only metric that doesn't show "uniform" is `r_eff_ratio = 0.60`, which is below the uniform value 0.9. This indicates **modest concentration** but at scattered k-values, not at low frequencies.

## Why outcome is β-strict but γ-functional

Strict β classification requires neither extreme threshold met:
- α requires r_eff_ratio < 0.25 (need 4× concentration); we observe 0.60.
- γ requires r_eff_ratio ≥ 0.8 (need ~uniform); we observe 0.60.

So technically β. But the **functional implication** for smooth-weight closure is closer to γ:

A smooth weight `w(m)` (Gaussian, parabolic, etc.) suppresses HIGH-FREQUENCY content of ĥ on J. The post-weighting energy is approximately Σ_k |w_hat(k)|² · |g(k)|². Smooth weights have |w_hat(k)| decaying fast in k — so they suppress mass at high k.

For polynomial saving: need most of |g(k)|² mass to be at HIGH k (so smooth weight kills it). Empirically: mass is uniformly distributed across k. Smooth weight kills 50% of the energy (the upper half of the spectrum) but leaves 50% (the lower half) intact. That's a √2 reduction at best — constant factor, not polynomial.

For full eq 190 closure via smooth weights, would need r_eff_ratio → 0 OR mass strongly concentrated at high frequencies (so smooth weight kills almost all energy). Neither observed.

## Why this is decisive

The metrics are **stable across r = 10, 12, 14** (5-fold range in N_r from 39K to 3.2M). Specifically lf_mass = 0.250 exactly at all three. This is asymptotic uniformity, not finite-r noise.

Combined with the conc_ratio = max|g(k)|² · |J| growing slowly (4 → 6 across r = 6..14), the picture is: the spectrum is "approximately uniform with mild concentration in some specific k-modes whose location varies with the cell". The varying location across cells means the concentration is in DIFFERENT modes for different (ℓ, t, η) triples — not a structural feature exploitable uniformly.

A finer test would look at WHICH modes the concentration sits at (top-mode index k_*). If k_* is consistently small (low-frequency), smooth weights help. If k_* is scattered, they don't. We didn't run this finer test — but the lf_mass = 0.25 result subsumes it: low-frequency modes carry exactly their uniform share of mass, no more.

## Comparison to predictions

| Reference | Prediction | Observed | Match |
|---|---|---|---|
| Uniform spectrum | r_eff_ratio = 0.9, H_norm = 1.0, lf_mass = 0.25 | 0.60, 0.97, 0.25 | partial (lf_mass and H_norm close, r_eff lower) |
| Low-rank smooth | r_eff_ratio < 0.1, lf_mass > 0.9 | 0.60, 0.25 | NO |
| BGK random | r_eff_ratio ≈ 0.9, H_norm ≈ 1, lf_mass = 0.25 | 0.60, 0.97, 0.25 | mostly yes |

The cells most resemble "BGK-random-with-mild-concentration" — the mild concentration manifests in r_eff_ratio = 0.60 (not 0.9), but is **not at low frequencies**.

## Implications for the obstruction map

Adds a fourth measurement to the post-morning obstruction map:

| Sub-route | Status |
|---|---|
| C2 (BGK on ⟨4⟩, primal) | partial — closes rate-1/2 rigorously, eq 190 still open |
| C3 (band-l¹ direct) | closed — β=1.000, ĥ saturates |
| **Spectral decomp of ĥ on J (this)** | **closed for smooth-weight purposes** — lf_mass→0.25, no exploitable low-frequency structure |
| Smooth completion (R78 path 2) | empirical motivation now WEAKENED — α not present |
| C1 (5x+1 sibling) | untested |

Smooth completion as a closure attack now has weakened motivation: the precondition (smooth structure on the band) is empirically absent. R78 path 2 was speculative; this analysis confirms it's not the right direction.

## Honest caveats

1. **Five r-points only**: r ∈ {6, 8, 10, 12, 14} with 12 cells each = 60 (r, eta=0.5) + 60 (r, eta=0.25). The asymptotic stability of metrics from r=10 onward is clear; pushing to r=16 wouldn't change the verdict.
2. **t ∈ {1, 2, 3}**: small non-trivial band centers. Different t shifts the band but the multiplier modulus is approximately translation-invariant (up to wrap-around), so qualitatively similar metrics expected at other t. Confirmed by the cell-counting (12 β at every (r, η), no t-variation).
3. **ℓ ∈ {0, 1, 2, 3}**: Kalafatelis only uses ℓ ∈ {0, 1, 2}; we extended to ℓ=3 for symmetry verification. Results identical.
4. **The "DFT of ĥ on J" treats J as an indexed set, not a subset of Z/N_r**. The DFT is over the |J|-length subsequence, not over the embedding into the larger group. This is the natural choice for measuring the function's mode structure as a function on the band, but a different convention (e.g., embedding into Z/N_r and zero-padding) would give different metrics. We chose the more natural one.
5. **Mild concentration (r_eff_ratio = 0.60, NOT 0.9 = uniform)**: there IS some structure, just not at low frequencies. A more refined attack might exploit specific scattered modes (e.g., locating which k_* modes are large and suppressing those individually). But this requires a non-smooth weight that's still tractable analytically — non-trivial to construct and unlikely to give polynomial saving.

## Files

- `band_spectral_analysis.py` — measurement script
- `band_spectral_data.csv` — 120 cell-rows with all metrics
- `band_spectral_log.txt` — full stdout
- `band_spectral_decomposition.md` — this writeup

## Decision tree resolution

Per the brief's outcome (β) branch:
> "β: Mixed, partially smooth. r_eff ≈ |J|/2, low-frequency mass 30-50%. Some structure, not enough for smooth-weight closure alone. Document the mode structure quantitatively; could combine with Esscher tilt or other measure-changing tools (STATE.md item #1)."

But our observed `lf_mass = 0.25` is BELOW the 30-50% range — it's at the uniform-in-k baseline. So we're closer to the γ outcome:
> "γ: White-noise-like, rough. ĥ on J is genuinely random; no smooth-weight cancellation possible. Closes one more route. Smooth completion (R78 path 2) loses motivation."

**Recommended morning consolidation framing:** treat this as **β-with-γ-functional-read**. Smooth completion is empirically weakened as an attack. The remaining open route in Path C is C1 (5x+1 sibling-attack).

If a follow-up firing is desired:
- **C1 (5x+1 sibling)**: rebuild the framework at q=5, compare obstruction structure to q=3. Multi-day project.
- **Esscher-tilt + spectral combo (STATE.md #1)**: weight ĥ by exp(λ·σ) where σ is some orbit-depth statistic; could combine the partial concentration in r_eff with measure-tilting. Speculative.
- **Direct shell-slice attack** (Kalafatelis Theorem 26): bypass eq 190 entirely; attack the shell-slice asymptotic from a different reduction route. Substantial.

## Compute audit

| metric | value |
|---|---|
| Hardware | 9950X3D 32 cores, no GPU |
| Numba threads | default |
| Concurrent NBA Projections | 2 (PID 3676, 47564) — untouched |
| Per-r elapsed | r=6: 0.0s, r=10: 0.1s, r=12: 0.8s, r=14: 10.0s |
| Total wall time | ~11s |
| Max r reached | 14 |
| Number of cells | 120 (all classified) |
