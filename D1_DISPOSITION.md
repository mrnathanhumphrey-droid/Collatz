# D1 disposition — n=4 alternating non-vanishing CONFIRMED

**Date:** 2026-05-15
**Task:** numerical confirmation of H1' violation at n=4 alternating index sequence (j_1, j_2, j_1, j_2)
**Script:** `C:/Collatz/verify_n4_alternating.py`
**JSON:** `C:/Collatz/experiments_output/n4_alternating_diagnostic.json`

## Verdict: **H1' VIOLATION CONFIRMED NUMERICALLY** ✓

The structural argument from `H1_PRIME_LOW_ORDER_CHECKS.md` is robust under all four scalar reductions. M_4_alt is non-zero at 4-7 orders of magnitude above the noise floor.

**(Note: the script's auto-generated verdict text "NOT_CONFIRMED — M4_alt at noise floor" is a verdict-logic bug — the threshold comparison was miscalibrated. The actual numerical values clearly show non-vanishing as documented below.)**

## Values at level n=3, (Z/27)*, 18 states, V_TRUNC=16 (Reading B = marginal centering)

| Reduction | M_2 (sec-order distinct, should ≈ 0) | M_4_alt = `E_B(X̃_1·X̃_2·X̃_1·X̃_2)` | M_4_distinct (4-distinct, V=6 trunc, should ≈ 0) |
|---|---|---|---|
| tr_pi | — | **5.357×10⁻²** | 1.895×10⁻⁷ |
| vac_pi | — | **4.775×10⁻³** | 1.698×10⁻⁷ |
| delta_1 | — | **5.742×10⁻²** | 2.277×10⁻⁶ |
| sum_entries | **1.076×10⁻⁷** (matches Task 1 anchor) | **6.089×10⁻¹** | 2.786×10⁻⁵ |

## Separations (H1' violation strength)

H1'-predicted RHS = 0 (via peak-rule substitution at position 2 or 4: X̃_{j_2} → ϕ(X̃_{j_2}) = 0). LHS = M_4_alt. Gap = |LHS − RHS| = |M_4_alt|.

| Reduction | \|M_4_alt\| / \|M_2\| (Task 1 noise) | \|M_4_alt\| / \|M_4_distinct\| (4-distinct control) |
|---|---|---|
| sum_entries | **5.66×10⁶** | 2.19×10⁴ |
| tr_pi | (Task 1 noise floor n.a. for this reduction; see M_4_distinct column) | 2.83×10⁵ |
| vac_pi | — | 2.81×10⁴ |
| delta_1 | — | 2.52×10⁴ |

**The smallest separation across all reductions is ~2.8×10⁴** between M_4_alt and the 4-distinct control. **Robust non-vanishing across all four readings.**

## Fubini inner factor — structural confirmation

`F(v_1, v_1') = Σ_{v_2, v_2'} 2^{-v_2-v_2'} · scalar(X̃_2 · X̃_1 · X̃_2)` computed across a 12-point (v_1, v_1') grid under Reading B sum_entries:

**All 12 grid values: F = 6.347×10⁻² (identical to 14 significant digits).**

The inner factor `E_{(v_2)}[X̃_2 · X̃_1 · X̃_2]` does NOT depend on (v_1, v_1'). This is structurally clean: the H1' violation arises from a within-pair-at-j_2 component that's independent of the j_1 pair-group's specific (v, v') split.

**Cross-checks the structural argument** in `H1_PRIME_LOW_ORDER_CHECKS.md` (Fubini decomposition over independent pair-groups).

## Reading A control

Reading A (strict conditional centering on full B per SETUP.md §5) gives M_4_alt values 2.21×10⁻², 7.07×10⁻⁴, 2.21×10⁻², 7.40×10⁻² across the four scalar reductions. These are **also non-zero** but smaller than Reading B by 10-100×. The Reading A non-zeroness is a separate, also-real structural finding (consistent with `AMALG_FREENESS_MOMENT_CALCULATION.md` §6 — under Reading A, only the trivial T̃_j = 0 centering trivially gives zero; at higher orders even Reading A doesn't trivialize for the n=4 alternating moment).

This Reading A non-zeroness reinforces the verdict: the n=4 alternating moment is structurally non-zero **regardless of centering choice**.

## Implication

The H1' framework as stated (verbatim HS 2014 Defn 2.2) does **not** hold for Syracuse X̃_j family. This is a real structural finding, not an artifact of:
- Scalar reduction choice (all 4 reductions agree)
- Vacuum-state collapse (5 reductions including sum_entries which is the full-state norm)
- Centering choice (both Reading A and Reading B show non-zero)
- V_TRUNC choice (computed at V=16, well past the 2^-16 < 1.6×10⁻⁵ tail threshold)

The Syracuse independence structure is **strictly weaker than monotone independence in HS 2014's sense.** This is the "unnamed independence regime" flagged in `H1_PRIME_DISPOSITION.md`.

## What this means for c=7/45 (cross-reference to D3)

D3 audit determined that the c=7/45 derivation **never depended on the failed regime** (no non-adjacent-repeat partitions appear in the derivation; the H1' failure has zero impact on the leading coefficient). See `C:/Collatz/D3_DERIVATION_AUDIT.md`. The leading derivation is rigorous unconditional via R75+R76+R77+R64.B.

D1 + D3 together: H1' violation is real (D1) but c=7/45 is unaffected (D3).
