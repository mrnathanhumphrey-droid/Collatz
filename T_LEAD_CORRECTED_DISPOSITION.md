# T_LEAD_CORRECTED_DISPOSITION — top-level disposition

**Date:** 2026-05-12. Wilson, T_lead Nisoli-bypass re-evaluation at the corrected asymptotic rate. Follow-up to:
- seven-probe spectral trajectory (closed at H_M_RECURSION_UNDERSPECIFIED)
- PADE_NUMERICAL_DISPOSITION (z=2 REFUTED as leading singularity, slow-mode at z≈1.016 candidate)
- R77.7 v2 ε_7 exact-rational compute (39 min, witness-verified)
- Framework-reopening map (SESSION_DISPOSITIONS_2026_05_12 + STATE.md)

---

## DISPOSITION: **H_T_LEAD_CARRIES_DIFFERENT_RATE** (with structural-meaning content)

> **Headline.** T_lead = T_diag + Off_lin, computed exactly over Q from cross-freq machinery's W_+(g) = 2^{−g+1}/15 weights, has spectrum **{43/45, 0}** on (P_+, P_-) with the leading eigenvalue 43/45 ≈ **0.9556** sitting on the (1, 4) eigenvector. The eigenvalue is a clean rational over Q, derived without numerical fits, and structurally meaningful as **1 − Σ_g W_+(g) = 1 − 2/45** — the sum of off-diagonal cross-frequency weights contracting T_diag's λ=1 mode.
>
> **This eigenvalue is NEITHER 1/2 (the R77.3-falsified rate) NOR ρ ≈ 0.984 (the slow-mode rate from STATE.md's prior-session two-mode fit).** It's a distinct value with explicit algebraic origin. The 2.9% gap to 0.984 and 0.5% gap to its inverse-radius z = 45/43 ≈ 1.046 are within the noise of where PADE_NUMERICAL_DISPOSITION's Hadamard trend (n=13 → 1.57, inward-trending) is heading at n→∞, but the empirical data through n=13 is too transient to confirm 43/45 over 0.984.
>
> **Framework reopening status:** The R77.3 falsification at rate-1/2 becomes a SCOPE finding ("rate-1/2 wrong target; rate ≈ 0.956 is the within-level cross-freq answer"). Nisoli bypass at corrected rate is **possible in principle but significantly harder** than at rate-1/2 (M_3'' ≈ 24 vs ~2-3 at rate-1/2). Requires the structural assumption that X̄_n(g≥2) projects onto (1,4) of (P_+, P_-) — an open assumption per T_V_DISPOSITION's findings.

---

## Pre-registered hypotheses, decided

| Hypothesis | Status |
|---|---|
| H_T_LEAD_HITS_CORRECTED_RATE (λ = ρ_empirical ≈ 0.984 within ~5%) | **PARTIAL** — T_lead's eigenvalue is 43/45 ≈ 0.9556, which is 2.9% from 0.984 (within 5% threshold). However, n=13 empirical data is too transient to confirm. |
| **H_T_LEAD_CARRIES_DIFFERENT_RATE** | **CHOSEN** — λ = 43/45 is a clean rational over Q, distinct from both 1/2 (falsified) and 0.984 (prior-session two-mode fit). Structural meaning: sum of off-diagonal cross-freq weights. |
| H_T_LEAD_NO_CORRECTED_EIGENVALUE | NO — clean Q-spectrum derived. |
| H_OFFLIN_INSUFFICIENT_FOR_MATRIX | NO — cross-freq §7 gives the (1,4)-image of Off_lin's rank-1 structure (sufficient to determine T_lead's spectrum). |
| H_NEED_NEW_OPERATOR | NO — T_lead exists, has clean spectrum. The OPEN QUESTION is whether 43/45 is the asymptotic rate. |

H_T_LEAD_HITS_CORRECTED_RATE was the pre-registered favored hypothesis (specifically the complex-conjugate-pair version). The actual outcome is a REAL eigenvalue at 43/45 (T_lead is rank-1, so only one non-trivial eigenvalue; the (1,4) direction yields a REAL value not a complex pair). The 2.9% gap to 0.984 is small enough that H_T_LEAD_HITS_CORRECTED_RATE is a defensible interpretation — but the cleanest reading is "T_lead carries 43/45, which is in the same ballpark as 0.984 but is a structurally-distinct value."

The complex-conjugate-pair prediction (period-9 oscillation) is NOT consistent with T_lead's real spectrum; the period-9 phenomenon must live in a DIFFERENT operator (inter-level R̃ perhaps, or in the cascade of moments outside V_M per T_V findings).

---

## Phase summary

### Phase 1: Off_lin sum from cross-freq weights, T_lead's action on (1, 4)

  W_+(g) = 2^{−g+1}/15 for g ∈ {2, 4, 6, ...}
  Σ_g W_+(g) = (1/30) · 1/(1 − 1/4) = (1/30) · (4/3) = **2/45**
  Σ_g W_-(g) = 4 · (2/45) = **8/45**

  Off_lin's (1, 4)-eigenvalue = −2/45 (sign forced by empirical S_n → 7/15)
  T_lead's (1, 4)-eigenvalue = 1 − 2/45 = **43/45**

### Phase 2: T_lead's 2D spectrum

  T_lead = [[7/45, 9/45], [28/45, 36/45]] (rank-1)
  
  Characteristic polynomial: λ² − (43/45) λ = 0
  
  **Spectrum {43/45, 0}** with eigenvectors (1, 4) and (9, −7).

### Phase 4: empirical cross-check

|Prediction | T_lead (43/45) | 0.984 (STATE prior) | 1/2 (R77.3 falsified) |
|---|---|---|---|
|Late |ε_n/ε_{n-1}| geomean k=11..13 = 1.6 | gap 67% | gap 63% | gap 220% |
|Hadamard radius at n=13 = 1.57 | gap 50% | gap 54% | gap 22% (but trending past) |
|Hadamard prediction n→∞ | 1.046 | 1.016 | 2.0 (refuted by n=10) |

Through n=13, the empirical data is in a transient with Hadamard inward-trending from 2.06 (n=10) to 1.57 (n=13). Both 43/45 and 0.984 are consistent with the inward direction; the data does not distinguish them.

### Phase 5: M_3'' resolvent norm

  ||(I − T_lead)^{−1}|| ≈ **24.4** (operator ℓ² norm)

Significantly larger than the rate-1/2 case (M_3'' ~ 2-3). Nisoli closure inequality |K|·K^{−A}·M_3'' < 1 is possible in principle but requires Nisoli's A > 1 and r-iteration to dominate.

---

## Trajectory placement (eight probes now)

| Probe | Object | Disposition |
|---|---|---|
| T_3 (R77.3) | 3x3 companion matrix at rate-1/2 | FALSIFIED |
| R_k (R77.4 erratum §1) | Inter-level residual operator | H_R_K_INTRACTABLE |
| Candidate A | W_k φ_n bilinear-pair-form | H_CANDIDATE_A_FALSIFIES_F2 |
| R76 §11 2D | T_diag + Off conjectural at rate-1/2 | INCONCLUSIVE |
| T_N construction | T_diag + Off_lin as 2x2 at rate-1/2 | H_OFF_LIN_UNDERSPECIFIED |
| Cross-freq closure | Closure space for Off_lin at rate-1/2 | H_CROSS_CLOSES_ON_ENLARGED_SPAN |
| T_V spectrum | T_V on V_M^{(g_max)} at rate-1/2 | H_M_RECURSION_UNDERSPECIFIED |
| **T_lead corrected (this)** | T_lead at corrected rate | **H_T_LEAD_CARRIES_DIFFERENT_RATE (43/45 over Q)** |

The eighth probe is **the first positive structural finding at the corrected rate**. It produces a clean Q-rational eigenvalue 43/45 with explicit algebraic origin (1 − Σ_g W_+(g)). This is meaningfully different from the prior seven probes (which all hit obstructions at the rate-1/2 target). The shift from rate-1/2 to rate-43/45 is the SCOPE correction that the framework-reopening map predicted.

---

## Reconciliation with prior findings

### With PADE_NUMERICAL_DISPOSITION's branch-cut/multi-singularity reading

PADE_NUMERICAL_DISPOSITION reads "two singularities visible at n=13; leading at |z| ≈ 1.5..1.7, sub-leading at z=2; long-term asymptotic at z ≈ 1.016". T_lead's eigenvalue 43/45 corresponds to singularity at |z| = 45/43 ≈ **1.0465**.

If the long-term asymptotic singularity is at z = 1.016 (radius 1/0.984), then T_lead's 1.046 is ABOUT 3% OFF. This could mean:
- (a) T_lead's 43/45 captures the WITHIN-LEVEL contribution, with a small inter-level correction shifting the final asymptote from 1.046 to 1.016 (Δ_inter ≈ -0.030 / 1.046 ≈ -2.9%).
- (b) The "0.984" was a finite-window fit artifact; the true asymptote is 43/45 ≈ 0.956 (and the slow-mode fit at k=7..10 in the prior session was contaminated by the transient).

Both readings are consistent with the available data. Resolving (a) vs (b) requires either k≥18 ε measurements (out of reach for now) or a separate cross-check via the inter-level R̃ operator.

### With T_V_DISPOSITION's H_M_RECURSION_UNDERSPECIFIED

T_V_DISPOSITION found that V_M doesn't close under iteration (phase + parity obstructions). T_lead is a WITHIN-LEVEL 2x2 operator that bypasses this issue (no iteration; just one-step recursion). So T_lead's 43/45 is well-defined within-level, but its prediction for the asymptotic rate INVOKES an implicit assumption that the within-level eigenvalue carries over to the n→∞ limit — which T_V's obstructions cast doubt on.

The CLEANEST reading: T_lead at 43/45 is the within-level rate-carrier prediction, and the asymptotic rate may differ from it due to the cascade of moments outside V_M (the inter-level renormalization piece). The slow-mode rate 0.984 might encode this inter-level correction.

### With R77.6's branch-cut interpretation

R77.6 reads the rate-1/2 phenomenon as branch-cut at z=2 (continuous spectrum endpoint, not discrete eigenvalue). PADE_NUMERICAL_DISPOSITION refutes z=2 as the leading singularity at n=13. T_lead's corrected eigenvalue 43/45 (singularity at z=1.046) is consistent with a DIFFERENT discrete eigenvalue closer to z=1. Whether this is genuinely discrete (per Nisoli framework) or a branch-cut endpoint at z=1.046 is unresolved.

---

## Routing recommendations (surfaced for Nathan, not chosen)

### Route A: Push Nisoli closure at λ = 43/45

Compute Nisoli's A explicitly (from project's R77 framework), check whether |K|·K^{-A}·M_3'' < 1 holds for r = 3, 4, 5 with M_3'' = 24.4 and known |K| ≤ 2√N at r ≤ 3. If the closure inequality is satisfiable at some r, the corrected-rate Nisoli closure WORKS and c = 7/45 has its first rigorous spectral closure (at rate 43/45, not 1/2).

If the closure inequality fails at all attainable r, then either:
- T_lead's 43/45 is NOT the asymptotic rate (needs to be smaller to give better M_3''), or
- The Nisoli framework structural-fails here too (consistent with the prior probes).

Estimated effort: 1-2 sessions.

### Route B: Compute the inter-level R̃ operator at the corrected rate

result_R_operator_spectrum.md tested R̃_k for ρ ≈ 0.984 and found |λ_2| at 10^-5 to 10^-3 (algebraically because R̃ as defined reduces to K_k). At the corrected rate 43/45, the R̃ operator might also have a CC pair near 43/45 that the prior probe didn't look for. Specifically, the period-9 oscillation phenomenon (sign pattern + + - - - - - - - + + + +) might live in R̃ with arg ≈ 2π/9 ≈ 0.683 rad and modulus near 43/45 OR 0.984.

But R̃ in its current form (L · K^m · P) algebraically reduces to K_k, so the slow mode is unlikely to live there. A DIFFERENT inter-level operator (the "non-square deviation operator R: level-k → level-{k+1}" alternative noted in result_R_operator_spectrum §"Three places the actual slow mode might live") may be the right object.

Estimated effort: 2-3 sessions for the alternative R operator construction.

### Route C: Acknowledge the partial structural finding and document for paper

T_lead's eigenvalue 43/45 is the first POSITIVE algebraic result from the spectral-probe trajectory. Even if the asymptotic rate isn't exactly 43/45, the closed-form 1 − Σ_g W_+(g) = 1 − 2/45 is a structural milestone: it identifies the within-level off-diagonal cross-frequency contribution as the contraction mechanism, and computes it exactly.

For paper purposes, this is a publishable structural anchor:

> "The within-level T_lead = T_diag + Off_lin (constructed exactly from R77 sketch §5 + cross-frequency closure machinery) has spectrum {43/45, 0} over Q on the (P_+, P_-) class-resolved moment space, with the (1, 4) eigenvector preserving the squared-class-mass ratio (1/3)² : (2/3)² = 1 : 4 from R64.B. The eigenvalue 43/45 = 1 − Σ_g W_+(g) is the explicit closed-form contraction from cross-frequency bilinear weights."

Combined with the bilinear bound (user's memory `project_collatz_r78_bilinear_cracked`) and the seven prior probes' negative spectral findings, this adds a NINTH spectral anchor (positive, this time) to the c = 7/45 structural ledger.

Estimated effort: 0.5 session (writeup) if not pushing further; 1-2 sessions if combining with Route A.

---

## Adversarial check outcomes (probe-level)

**(A1) Cross-freq fidelity.** Re-derived W_+(g) = 2^{-g+1}/15 directly from cross_freq §6 verbatim. Geometric sum Σ_g W_+(g) = 2/45 via standard geometric series. T_lead eigenvalue 43/45 follows by adding T_diag's (1,4)-eigenvalue 1 and Off_lin's (-2/45) contraction. All steps traced. ✓

**(A2) Period-9 vs single real eigenvalue.** T_lead has REAL spectrum {43/45, 0} on (P_+, P_-). The empirical period-9 oscillation cannot live in T_lead. It must live in a different operator (inter-level R̃ or similar). The seven-probe trajectory's negative R̃ finding at ρ = 0.984 (now Tier 5, walked back) leaves this unresolved. **The corrected-rate T_lead does NOT explain period-9 oscillation.** ✓ — flagged honestly.

**(A3) Exact rationals throughout.** All Off_lin entries, T_lead entries, eigenvalues, and resolvent norm computation use Fractions where applicable. Numerical operator norm computation uses floats; spectral radius 22.5 is exact. ✓

**(A4) Stability across n.** T_lead's matrix entries are absolute constants (cross_freq §6 weights are n-independent). So eigenvalue 43/45 is the same at every n ≥ 2. ✓ — but this is a tautology of the construction; whether 43/45 is the n→∞ asymptotic rate is the open question.

**(A5) Empirical fit at n=7..13.** Honest read: data is transient; doesn't decisively confirm 43/45 over 0.984 or alternatives. Both 43/45 and 0.984 are consistent with the inward Hadamard trend through n=13. ✓ — flagged.

**(A6) Reconciliation with R77.6 + cross_freq + T_V.** Documented in §"Reconciliation with prior findings". T_lead's discrete eigenvalue 43/45 may be a within-level approximation of an asymptotic operator whose true spectrum has the slow-mode at 0.984 + period-9 oscillation. The 2.9% gap may encode the inter-level correction. ✓

---

## Deliverables

In C:/Collatz/:

- T_LEAD_CORRECTED_PHASE1.md — cross-freq weight sums + (1,4)-action
- T_LEAD_CORRECTED_SPECTRUM.md — T_lead 2D spectrum {43/45, 0}
- T_LEAD_CORRECTED_EMPIRICAL.md — empirical cross-check at ε_1..ε_13
- T_LEAD_CORRECTED_CLOSURE.md — Phase 5 resolvent norm M_3'' ≈ 24.4
- T_LEAD_CORRECTED_DISPOSITION.md (this file) — top-level
- t_lead_corrected.py — main-thread verification script

---

## Synopsis (one paragraph)

The eighth spectral probe re-evaluated T_lead = T_diag + Off_lin at the corrected rate (not 1/2, per R77.3 falsification + Padé extension's z=2 refutation). The cross-freq machinery's W_+(g) = 2^{-g+1}/15 weights sum to Σ_g W_+(g) = 2/45 over g ∈ {2, 4, 6, ...}. T_lead's eigenvalue on (1, 4) is exactly **43/45 ≈ 0.9556 over Q** (rank-1 operator, spectrum {43/45, 0}). This is distinct from both 1/2 (R77.3 falsified) and 0.984 (STATE.md prior-session two-mode fit), structurally meaningful as "1 minus the sum of off-diagonal cross-freq weights". The empirical data through n=13 is in a Padé-numerical-confirmed transient and does NOT decisively distinguish 43/45 from 0.984; both are consistent with the Hadamard-inward-trending direction at n=13 (radius 1.57 trending toward [1.02, 1.05]). T_lead at 43/45 is the first POSITIVE algebraic result from the spectral-probe trajectory — a clean Q-rational with explicit closed-form origin. Nisoli resolvent norm M_3'' = ||(I − T_lead)^{-1}|| ≈ 24.4, significantly larger than at rate-1/2 (~2-3) but still admitting closure in principle if Nisoli's A > 1. The honest read: this opens the corrected-rate framework reopening with a concrete eigenvalue but doesn't yet close c = 7/45 rigorously — the (1, 4)-projection of X̄_n(g≥2) remains a structural assumption (cross_freq §7), and the period-9 empirical oscillation isn't captured by T_lead's real spectrum (must live elsewhere, likely in an inter-level operator).
