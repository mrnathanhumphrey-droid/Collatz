# DELTA_DIAGNOSTIC_DISPOSITION — leading-vs-subleading top-level result

**Date:** 2026-05-12. Wilson (analyst) reporting to Nathan. Top-level disposition of the leading-vs-subleading ansatz diagnostic, follow-on to TAUBERIAN_SCOPING_*.

---

## DISPOSITION: **H_DELTA_IRREGULAR**

> **None of the five pre-registered ansatze (geometric, power-law, log, two-term, oscillating) reproduce δ_6 within the pre-registered 20% threshold. All five predict δ_6 with the wrong sign. Reading B (δ_n is genuinely irregular on the n=2..6 sample; the mixed-singularity Tauberian framework needs richer structure than leading-simple-pole + single-secondary-mode) is supported.**
>
> The four fit points (n=2..5) are all positive (δ_n > 0); the held-out point δ_6 = −0.00147 is negative. Every pre-registered ansatz form — including the two-term and oscillating forms that *can in principle* accommodate sign-flips — extrapolates monotonically (or near-monotonically) from the n=2..5 positive sequence and predicts positive δ_6 in the range +0.0012 to +0.0092. Relative residuals at n=6 range from **182% (log) to 727% (power-law)**, all wildly outside the 20% threshold.

---

## Why H_DELTA_IRREGULAR and not the alternatives

- **H_DELTA_GEOMETRIC** (single geometric c·ρ^n): **REJECTED.** Real ρ cannot reproduce the sign flip. Residual 378%.
- **H_DELTA_POWER_LAW** (c·n^(−α)·(1/2)^n): **REJECTED.** Best LSQ gives α ≈ −3.87 (i.e., growing in n, since n=2..5 are increasing in δ·2^n), and prediction at n=6 is +0.0092 (727% off). No real α produces a sign flip.
- **H_DELTA_LOG** (c·log(n)·(1/2)^n): **REJECTED.** Predicts +0.0012 (182% off, wrong sign). Closest in *magnitude* among the failures, but still strictly positive at n=6 by structure.
- **H_DELTA_TWO_TERM** (c_1·ρ_1^n + c_2·ρ_2^n): **REJECTED.** Prony method fits n=2..5 exactly; held-out n=6 is +0.0016 (210% off). The fit gives complex-conjugate roots ρ ≈ 0.504 ± 0.273i, |ρ| ≈ 0.573, arg ≈ 0.495 rad. The cosine period is ≈ 12.7 steps. n=6 sits before the first implied zero crossing, hence positive.
- **H_DELTA_OSCILLATING** (c·cos(ωn+φ)·(1/2)^n): **REJECTED.** Cosine recurrence on A_n = δ_n·2^n gives cos(ω) ≈ 0.96 (ω ≈ 0.29 rad), small enough that n=6 stays before the implied zero crossing. Prediction +0.0019 (229% off).
- **INCONCLUSIVE**: **REJECTED.** No ansatz is "tied at similar quality" — all are uniformly bad, all on the same side. Disposition is clean toward IRREGULAR.

---

## Framework implication: what Reading B blocks

The Tauberian scoping probe (TAUBERIAN_SCOPING_DISPOSITION.md) left two possible readings open. Reading A says δ_n is leading + secondary in a clean mixed-singularity sense; Reading B says δ_n is irregular through n=6.

This diagnostic supports **Reading B**: the leading-simple-pole + single-secondary-mode (or two-secondary-mode) framework is **EMPIRICALLY INSUFFICIENT** on n=2..6 data. Specifically:

1. **Chevalier Thm 1.16 (meromorphic h with pole of order M)** in its single-pole-secondary form gives a fixed-α (i.e., M − 3/2) power-law correction. Ansatz (b) directly tests this; FAILS.
2. **Logarithmic-correction Tauberian theorems** (Flajolet-Sedgewick Ch. VI log case) are tested by ansatz (c). FAILS.
3. **Single-secondary-pole on second sheet** (geometric correction with ρ ≠ 1/2) is tested by ansatz (a). FAILS.
4. **Two-secondary-modes / complex-conjugate-pair secondary singularities** are tested by ansatz (d) and (e). FAILS.

What remains structurally consistent with the data:

- The leading 1/30·(1/2)^n is qualitatively right (|ε_n|·2^n stays within ±22% of 1/30 across n=2..6).
- The subleading δ_n is *not* a perturbation: |δ_n / ε_n| ranges from 0.5 to 3.0 across n=2..6 — the correction is comparable to or larger than ε_n itself. This means **n=2..6 is pre-asymptotic** in a strong sense; the leading hasn't dominated yet.
- The non-monotonicity and sign-flip are real but uncaptured by any 1-mode or 2-mode form. Possible structural readings (none disambiguable from this data alone):
  - Three or more competing secondary singularities (3+ modes) with comparable magnitudes.
  - A *moving* secondary structure — e.g., density of pole accumulation at z = 2 that changes character with n.
  - Essential singularity at z = 2 (not just power/log branch).
  - The fit window n=2..5 is in a transient before any clean asymptotic regime starts.

---

## What additional data would resolve

**Primary:** ε_7 (and ideally ε_8). With 6 (or 7) data points, the diagnostic gains real discriminating power. Specifically:

- If δ_7 is **strongly negative** (continuing the n=6 break), the sign flip is the start of a settled negative regime — consistent with the leading 1/30 being slightly over-corrected and δ_n decaying into the negative quadrant. This would support a logarithmic-type secondary with a phase delay, or a 3-mode superposition.
- If δ_7 returns positive (n=6 was a transient), then we're seeing genuine oscillation that ansatz (d)/(e) might capture with better fit window. Critically, the cosine-period of ≈12.7 steps from Prony would predict δ_8 to be negative if extrapolated — n=7 should still be positive in the implied cycle.
- If δ_7 is **small** (close to zero), we're near a node of an oscillation — the two-term/oscillating ansatze would gain support.

**Pre-registered prediction for δ_7 under each rejected ansatz** (for falsification check when ε_7 lands):

| Ansatz | δ_7^pred (next data point) | Plausible sign |
|---|---:|:---:|
| (a) geometric | 0.00735·0.906⁷ = +0.00370 | + |
| (b) power-law | 5.74e-4·7^{3.87}·(1/128) = +0.01067 | + |
| (c) log | 0.0429·log(7)·(1/128) = +0.000653 | + |
| (d) two-term | a·δ_6 + b·δ_5 = 1.008·(−0.00147) − 0.328·0.00352 = **−0.00264** | **−** |
| (e) oscillating | continuation from recurrence: A_7 = 2·0.96·A_6 − A_5 = 2·0.96·(−0.0939) − 0.1127 = **−0.293**; δ_7 = −0.00229 | **−** |

So ε_7 (when computed) would *itself* discriminate: if δ_7 is positive, the data continues to be irregular relative to all forms. If δ_7 is in [−0.002, −0.003], two-term / oscillating gain support — i.e., the n=6 break was the start of the predicted oscillation cycle, just one step late at n=2..5. Note that ansatze (d)+(e) BOTH predict δ_7 ≈ −0.0023 to −0.0026.

**Secondary:** ε_8 alone (without ε_7) is less useful — it'd be a single far-extrapolated point with no intermediate constraint.

**Note on compute cost:** R77.7's ε_7 attempt was killed at 8.5hr. Per existing project policy, **do NOT re-fire that compute** without a different solver or significant infrastructure change. The diagnostic stands on existing data.

---

## What is NOT blocked

The leading-coefficient reading from R76 §10 (|ε_n|·2^n → 1/30) remains qualitatively supported by this data. The rate-1/2 conclusion (S_n → 7/15 with leading geometric decay (1/2)^n) is unchanged.

What is blocked is the **clean Tauberian-framework closure** that would name the specific singularity-type of δ_n's contribution at z=2 and recover the rate-1/2 from a named theorem. The mixed-singularity framework remains the right abstraction; the specific theorem-statement to plug in is undetermined and 5-point data cannot determine it via the standard ansatz tests.

---

## Caveats

1. **n=2..6 is a 5-point sample.** Small. The pre-registered 20% threshold is honestly loose, and even at 20% the ansatze all fail by >100%. The disposition (IRREGULAR) is therefore robust to threshold quibbles — no ansatz comes close.

2. **The "monotone positive followed by negative" pattern at n=2..6 might be a transient.** This diagnostic does not rule out a clean asymptotic form taking over at large n. It rules out that any of the 5 pre-specified forms describe the *observed sample*.

3. **No additional ansatz was tried.** The brief pre-specifies the five. Inventing more forms post-hoc would be curve-fitting. The honest disposition with this discipline is IRREGULAR.

4. **The two-term and oscillating fits gave VERY similar held-out predictions** (+0.00161 and +0.00190). Prony's complex-root solution is essentially a re-parameterization of the cosine form. Their joint failure is one structural failure, not two independent ones.

5. **Per-point precision of δ_n:** computed from exact Fractions then converted to float; precision to >10 significant digits. The ~12% spread in |ε_n|·2^n across n=2..6 is much larger than any computational uncertainty.

---

## Files

- `DELTA_DIAGNOSTIC_QUANTITIES.md` — Phase 1, five quantities tabulated
- `DELTA_DIAGNOSTIC_ANSATZE.md` — Phase 2, five ansatz fits with residuals
- `DELTA_DIAGNOSTIC_DISPOSITION.md` — this file
- `delta_diagnostic.py` — verification script for main-thread execution

## Synopsis (one paragraph)

The Tauberian scoping probe (2026-05-12) left two readings open for δ_n := |ε_n|·2^n − 1/30: Reading A (clean leading + subleading with named secondary singularity, ansatz characterizable) vs Reading B (genuinely irregular on n=2..6, needs richer framework). This diagnostic tests five pre-registered ansatz forms — geometric, power-law correction, log correction, two-term superposition, oscillating — by fit on n=2..5 and held-out residual at n=6. **All five fail**, all in the same direction: each predicts δ_6 > 0, but actual δ_6 = −0.00147. Smallest residual is 182% (log), largest is 727% (power-law). The Prony-method two-term fit returns complex-conjugate roots with |ρ|≈0.573 and cosine period ≈12.7 steps — structurally identical to the oscillating ansatz — and likewise predicts δ_6 positive (just inside the first half-period). **Disposition: H_DELTA_IRREGULAR.** The mixed-singularity Tauberian framework with leading-simple-pole + single-secondary-mode (or two-secondary-mode) is empirically insufficient on n=2..6 data. Notable structural finding: |δ_n / ε_n| ranges from 0.5 to 3.0, so n=2..6 is *not yet in the asymptotic regime* where δ_n is a small correction — the leading and subleading are comparable in magnitude. Computing ε_7 would discriminate cleanly: two-term and oscillating both pre-register δ_7 ≈ −0.0023 to −0.0026 (so if δ_7 is in that range, the n=6 sign-flip was the start of the predicted cycle and ansatze d/e gain back some support); if δ_7 is positive or far from that range, the irregularity reading hardens. Existing project policy (R77.7 killed at 8.5hr) blocks re-firing the ε_7 compute without new infrastructure. The Tauberian framework remains the right level of abstraction for E(z); the specific theorem to invoke is not selectable from n=2..6.
