# HENSEL Phase 3 — Numerical verification

**Date:** 2026-05-11. Analyst: Wilson. Reports to: Nathan.

## Status: SCRIPT_WRITTEN_PYTHON_DENIED (one hand check passes)

Python execution denied in this session. Verification script `hensel_approach_a_verify.py` is written, self-contained, and ready for main-thread execution.

One hand-computation cell was verified: **(p=3, r=4, a=1)** — the closed-form polynomial prediction `P_a(s*(r=4)) ≡ -p²·s*²/2 + p³·s*³/6 - p^4·s*^4/12 mod p^5` matched the direct computation `P_a(s*(r=4)) = p·s* − C_a·L_3(1+p·s*) mod p^5` to integer precision (both gave 81 mod 243).

## Hand check details (p=3, r=4, a=1)

Computed via fractions and modular arithmetic:

| Quantity | Value |
|---|---|
| J_p (truncation level) | 5 |
| L_3(4) | 717/20 |
| L̃_3 = L_3(4)/3 mod 81 | 16 |
| L̃_3^{-1} mod 81 | 76 |
| C_a = a·L̃_3^{-1} mod 81 (a=1) | 76 |
| s*(r=4) = (C_a−1)/3 mod 27 | 25 |
| Base-3 digits of s* (s_0, s_1, s_2) | (1, 2, 2) |
| L_3(1+3·25) = L_3(76) mod 243 | 147 |
| **P_a(s*(r=4)) = 3·s* − 76·L_3(76) mod 243** | **81** |
| Closed form: −9·s*²/2 + 27·s*³/6 − 81·s*^4/12 mod 243 | **81** ✓ |

**Match:** the closed-form polynomial expansion produces exactly the same integer mod p^{r+1} as the direct evaluation. ✓

This confirms the Phase 1 articulation's `(1+y) log(1+y)` series identity reduction at one cell.

## Predicted vs actual G_p(a=1) at p=3, r=4

Predicted: `G_p(a=1) = p^{(r+1)/2} · η_p(r=4) · e_q(P_a(s*(r=4)))`

- `p^{(r+1)/2} = 3^{5/2} = 9·√3 ≈ 15.588`
- `η_3(r=4) = (1/√3)·Σ_{h=0}^{2} e_3(2·h²) = -i` (computed by hand)
- `e_{243}(81) = e^{2πi·81/243} = e^{2πi/3} = -1/2 + i·√3/2`

`Predicted G_p(a=1) = 9√3 · (-i) · (-1/2 + i√3/2)`
`= 9√3 · (i/2 + √3/2)`     [since -i·(-1/2) = i/2 and -i·(i√3/2) = -i²·√3/2 = √3/2]
`= (9√3 · √3)/2 + (9√3/2)·i`
`= 27/2 + (9√3/2)·i`
`≈ 13.5 + 7.794·i`
`|Predicted G_p| = √(13.5² + 7.794²) = √(243) = 15.588` ✓ magnitude saturation

**Actual G_p(a=1) at p=3, r=4** requires summing 81 terms; this is the numerical-verification target. Predicted value 13.5 + 7.794·i is the testable claim.

## Verification cells (script `hensel_approach_a_verify.py`)

| (p, r) | q = p^{r+1} | period = p^r | |supp| | J_p | Estimated runtime |
|---|---|---|---|---|---|
| (3, 4) | 243 | 81 | 27 | 5 | < 1 sec |
| (3, 5) | 729 | 243 | 81 | 6 | < 5 sec |
| (3, 6) | 2187 | 729 | 243 | 7 | < 30 sec |
| (5, 4) | 3125 | 625 | 125 | 5 | < 10 sec |
| (5, 5) | 15625 | 3125 | 625 | 5 | < 1 min |
| (7, 4) | 16807 | 2401 | 343 | 5 | < 1 min |
| (7, 5) | 117649 | 16807 | 2401 | 6 | < 5 min |
| (7, 6) | 823543 | 117649 | 16807 | 7 | < 30 min |
| (11, 4) | 161051 | 14641 | 1331 | 5 | < 5 min |
| (11, 5) | 1771561 | 161051 | 14641 | 6 | < 30 min |

**Verification target:** max over a ∈ supp of `|G_p(a)_actual − predicted_via_closed_form| / |G_p(a)_actual|` < `1e-12` (machine precision with mpmath dps=50).

## Pass/fail criteria

- **PASS at cell (p, r):** max rel dev < 1e-12. Confirms Approach A's closed form is exact at this cell.
- **PASS but with eta_p(r=4) absorbed into phase:** if predicted vs actual differ only by a constant a-independent factor of unit modulus, the closed form is correct modulo the η_p convention. Recompute η_p carefully.
- **FAIL at small p (p=3, p=5) only:** small-prime caveat from Phase 1 articulation; the closed form has p-specific corrections. Document and re-derive at small-p.
- **FAIL at large p (p ≥ 7):** Approach A derivation has an error; redo by hand at the failing cell.

## What this would establish

If the script PASSES at p ≥ 7 (clean range) for r ∈ {4, 5, 6}:

> **The family-level Hensel-lifted closed form is empirically verified at the (clean) cells. Combined with the structural Phase 1 derivation, this is strong evidence for H_HENSEL_CLOSES at family level, conditional on the open bilinear closure step (the nested inner-Plancherel chain).**

If the script PASSES at p=3 r=4..6 ALSO (small-prime caveats handled):

> **Family-level + small-prime closure. Strongest possible structural outcome.**

If the script FAILS at p=3 only:

> **Family-level closure for p ≥ 5; small-prime caveats at p=3 (the j=3 stratum-merge identified in Phase 1) require separate treatment, consistent with R79b's open-at-p=3 framing.**

## Action item for main thread

> **Run `python C:/Collatz/hensel_approach_a_verify.py` and inspect output `HENSEL_APPROACH_A_VERIFICATION.csv`.** Estimated total runtime: 1-2 hours (dominated by p=7, r=6 cell at ~30 min).
>
> If short on compute, run only `(p, r) ∈ {(3, 4), (5, 4), (7, 4), (11, 4)}` first (under 5 minutes total) to confirm the r=4 family-level closure.

## Caveats

1. **Hand check covered only the polynomial P_a(s*(r=4)) match, not the full G_p(a) match.** The Gauss-sum saturation magnitude is established structurally (Theorem 78.3); the η_p factor and the residual quadratic Gauss sum are derived by the digit-chain argument in Approach A. Numerical verification covers the COMBINED prediction.

2. **At p=3, the cubic-term stratum-merge** (Phase 1 articulation §"Small-prime caveats") may produce a different closed-form polynomial structure at p=3 than at p ≥ 5. The script's p=3 cells will reveal this if so.

3. **At p=5, the quintic-term j=5 case (1/5 has v_5 = -1)** lifts the j=5 contribution into the same stratum as j=4. Similar to p=3 cubic merge. Watch for p=5, r=5 cells.

4. **At p=7, r=7 the j=7 case kicks in** — same pattern. The script doesn't cover p=7, r=7 (would take too long); the structural pattern continues.

## Outcome categorization (pending script run)

> **PROVISIONAL APPROACH_A_EXACT.** Structural derivation produces the closed form; one hand-cell match at p=3, r=4 confirms the polynomial part. Full empirical verification awaits Python run.
