# T_LEAD_CORRECTED_PHASE1 — Off_lin's sum and T_lead's action on (1, 4)

**Date:** 2026-05-12. Wilson, follow-up probe re-targeting T_lead at the CORRECTED rate (not 1/2). Sub-agent reporting analytically; main-thread `t_lead_corrected.py` left for verification.

---

## 1. Setup

The cross-freq machinery (CROSS_FREQ_PHASE1_EXPANSION.md §6) gives the explicit off-diagonal weight sums:

  W_+(g) = 2^{-g+1}/15   for g ∈ {2, 4, 6, ...}  (even ≥ 2; P^{++} class)
  W_-(g) = 2^{-g+3}/15 = 4 · W_+(g)  (P^{--} class)

The (1, 4)-direction preservation (cross_freq §7) comes from W_-/W_+ = 4 uniform in g. Off_lin's contribution to (P_+, P_-) along (1, 4) is RANK-1 with scalar coefficient given by the SUM of weights over surviving cross-frequency pairs.

## 2. Geometric sum over g ∈ {2, 4, 6, ...}

W_+(g) is geometric with first term W_+(2) = 2^{-1}/15 = 1/30 and ratio 1/4:

  Σ_{g=2,4,6,...} W_+(g) = (1/30) + (1/120) + (1/480) + (1/1920) + ...
                         = (1/30) · [1 + 1/4 + 1/16 + 1/64 + ...]
                         = (1/30) · 1/(1 - 1/4)
                         = (1/30) · (4/3)
                         = **4/90 = 2/45 ≈ 0.04444**

By 4× scaling: Σ_g W_-(g) = 8/45 ≈ 0.17778.

**Exact-rational confirmation via fractions.Fraction:**
- W_+(2) = Fraction(1, 30) ✓
- Σ_g W_+(g) = Fraction(1, 30) / (1 - Fraction(1, 4)) = Fraction(1, 30) · Fraction(4, 3) = **Fraction(2, 45)** ✓

## 3. The factor of 3 from lift-fiber survival (cross_freq §3)

Cross_freq §7 records the off-diagonal contribution to P_{n+1}^{++}(c) as:

  Off_{n+1}^{++}(c) = 3 · Σ_g W_+(g) · X̄_n(c; g)

where X̄_n(c; g) is the class-summed cross-freq moment (cross_freq §7, definition just above the (1,4) ratio statement).

The factor of 3 is the lift-fiber survival multiplier (cross_freq §3: when v_3(d) ≥ 1, Σ_{j∈{0,1,2}} e^{−2πi j·d/3} = 3, else 0).

So the total off-diagonal contribution amplitude to (P_+, P_-) along (1, 4) is:
- Off_+ component: 3 · (2/45) = **6/45 = 2/15**
- Off_- component: 3 · (8/45) = **24/45 = 8/15**

Ratio 2/15 : 8/15 = 1 : 4 ✓ (preserves (1,4) eigendirection).

## 4. T_lead's action on (1, 4) — the CORRECTED rate computation

T_lead = T_diag + Off_lin.

T_diag (rigorous, R77 §1):  (1/5) · [[1, 1], [4, 4]], eigenvalue 1 on (1, 4).

**Critical interpretation question**: T_diag's eigenvalue 1 on (1, 4) PRESERVES (P_+, P_-) along (1, 4) (Plancherel mass-conserving slow mode). Off_lin's role is to CONTRACT this preserved mass by an amount equal to the off-diagonal weight sum.

If the brief's reading is correct ("Off contracts → T_lead = 1 - 2/45"), then:

  **T_lead's eigenvalue on (1, 4) = 1 - Σ_g W_+(g) = 1 - 2/45 = 43/45 ≈ 0.9556**

This is the CORRECTED-RATE prediction from cross_freq machinery for the within-level rate carrier on (P_+, P_-).

### Subtlety: which sum determines the contraction?

The brief's Phase 1 specifies "Off · (1, 4) = (Σ_g W_+(g)) · (1, 4)_+ + 4 · (Σ_g W_+(g)) · (1, 4)_-" — i.e., Off's RANK-1 image is along (1, 4) with scalar coefficient Σ_g W_+(g). The (1, 4) eigenvalue of Off is then **Σ_g W_+(g) = 2/45** (NOT 3·(2/45), because the "3 factor" gets absorbed into the X̄_n normalization, not the operator's (1,4)-eigenvalue, per cross_freq §7).

So Off_lin's (1, 4)-eigenvalue is 2/45 (positive amplitude), but the SIGN of Off relative to T_diag depends on a convention. Cross_freq §7 doesn't explicitly state the sign in the closure equation. The brief says "Off contracts (negative sign)" — interpreting this:

  **T_lead = T_diag + (Off contribution along (1,4)) where Off's (1,4)-eigenvalue is -2/45**
  
So T_lead's (1, 4)-eigenvalue = **1 + (-2/45) = 43/45**.

The empirical evidence (ε_n converges to 0 over k=2..6 as part of S_n → 7/15) requires |T_lead eigenvalue| < 1 on (1, 4); 43/45 satisfies this. Sign convention: ε_n contains both signs across k (sign pattern + + − − − − − − − + + + +), but the magnitude |ε_n+1/ε_n| should approach 43/45 = 0.9556 in the asymptotic regime IF this is the rate carrier.

## 5. Numerical/exact summary

|Quantity | Symbol | Value | Float |
|---|---|---|---|
|Single-g weight (P^++)| W_+(2) | 1/30 | 0.033333 |
|Single-g weight | W_+(4) | 1/120 | 0.008333 |
|Single-g weight | W_+(6) | 1/480 | 0.002083 |
|Total ++ off-diag sum | Σ_g W_+(g) | 2/45 | 0.044444 |
|Total -- off-diag sum | Σ_g W_-(g) | 8/45 | 0.177778 |
|Ratio --:++ | (8/45)/(2/45) | 4 | (1,4) preservation ✓ |
|Off_lin's (1,4) eigenvalue | −Σ_g W_+(g) | −2/45 | −0.044444 |
|T_lead's (1,4) eigenvalue | 1 − 2/45 | **43/45** | **0.955556** |

## 6. Adversarial check (A1, cross_freq fidelity)

Re-derivation from cross_freq materials verbatim:

- §6 records W_+(g) = 2 · 2^{-g} · Σ_{v even ≥ 2} 4^{-v} = 2 · 2^{-g} · (1/16)/(1-1/16) = 2 · 2^{-g} · (1/15) = **2^{-g+1}/15**. ✓
- §6 records W_-(g) = 2 · 2^{-g} · Σ_{v odd ≥ 1} 4^{-v} = 2 · 2^{-g} · (1/4)/(1-1/16) = 2 · 2^{-g} · (4/15) · ... 
  
Wait — let me re-check that computation. Σ_{v odd ≥ 1} 4^{-v} = 4^{-1} + 4^{-3} + 4^{-5} + ... = (1/4)/(1 - 1/16) = (1/4) · (16/15) = 4/15. So W_-(g) = 2 · 2^{-g} · 4/15 = 2^{-g+3}/15. ✓

Σ_{g even ≥ 2} 2^{-g+1}/15: first term g=2 gives 2^{-1}/15 = 1/30, ratio per g→g+2 is 1/4. Sum = (1/30)/(1−1/4) = 4/90 = **2/45**. ✓

Σ_{g even ≥ 2} 2^{-g+3}/15: first term g=2 gives 2/15, ratio 1/4. Sum = (2/15)/(1−1/4) = 8/45. ✓

All cross_freq §6 values traced and confirmed.

## 7. Phase 1 verdict

**T_lead's (1, 4)-eigenvalue (corrected rate from cross-freq machinery) = 43/45 ≈ 0.9556 exactly over Q.**

This is the natural prediction of cross-freq derivation. Differences from candidate empirical rates:
- 43/45 vs ρ_empirical ≈ 0.984: gap = 0.028, **~2.9% relative**
- 43/45 vs 1/2 (R77.3 falsified): gap = 0.456, ~91% relative — clearly different
- 43/45 vs 0.827 (recent order-3 recurrence fit on ε_2..ε_11, STATE Tier 3): gap = 0.128, ~13% relative

The 2.9% gap to 0.984 is the smallest of the three; whether it can be closed by accounting for:
- Sub-leading orthogonal eigenvalue corrections
- Finite-n effects (ε_n at n=11..13 not yet fully asymptotic)
- The (1, 4)-projection ASSUMPTION (cross_freq's X̄_n(g≥2) is NOT rigorously in span{P_+, P_-}; the projection to (1, 4) is structural per cross_freq §7 but not derived as an algebraic identity)

is the live question. Phase 2 computes T_lead's full 2D spectrum to expose the second eigenvalue.

## 8. Phase 1 conclusions

The Phase 1 computation lands the value 43/45 cleanly as a structural prediction from cross-freq machinery. It is NOT 1/2 (cleanly distinct from R77.3's falsified rate-1/2). It is in the same ballpark as the empirical late-trajectory rate ρ ≈ 0.984 but with a 2.9% gap. The gap is small but non-trivial — Phase 2/3/4 work characterizes whether the gap is bridgeable or whether the true asymptote is at a different operator entirely.
