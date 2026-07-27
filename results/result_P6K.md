# RESULT — P6K: MICROCOSM (λ-deformation) — λ=½ is critical; Collatz sits on the no-turnover side of a sign-flip (2026-07-26)

**MODEL-FAMILY RESULT — at λ≠½ the map is NOT Syracuse.** Every statement here is about the model family's behaviour at
and near the physical point λ=½, not about Collatz directly.

**Probe:** `probe_p6k.py`. Deform the valuation law to `P(v) ∝ λᵛ` (λ=½ = Syracuse `2⁻ᵛ`), keep the `÷2ᵛ·(3x+1)` map.
Drift `E[3·2⁻ᵛ] = 3(1−λ)/(2−λ)` = 1 at λ=½ (critical), <1 for λ>½ (converge), >1 for λ<½ (diverge). Object:
`T_i(λ) = 3ⁱΣ_{k≥1}4⁻ᵏ⟨ρ_λ,shiftₖρ_λ⟩` (= ½S at λ=½). `stationary_trunc_lam` = gated refactor (weight λᵛ).

## M-A — construction valid at the physical point
`T_i(0.5)` reproduces the certified S-ladder (P6H) to **1e-9** (power-iteration + truncation floor). The deformation
passes through the real object.

## M-B / M-D — λ=½ is the critical point (criticality identification CONFIRMED)
| λ | drift 3(1−λ)/(2−λ) | behaviour | conv. rate |
|---|---|---|---|
| 0.40 | 1.125 | diverges | >1 |
| 0.45 | 1.064 | diverges | ~1 |
| 0.48 | 1.026 | diverges (slow) | ~1 |
| **0.50** | **1.000** | **marginal (critical)** | ~1 (no conv) |
| 0.52 | 0.973 | converges | 0.30 |
| 0.55 | 0.931 | converges | 0.54 |
| 0.60 | 0.857 | converges | 0.50 |
**M-D sanity PASSES:** convergence speeds up sharply off-critical on the λ>½ side (rate 0.3–0.5 vs marginal at ½); λ<½
diverges exactly as `drift>1` predicts. The gap opens away from ½ — the criticality identification is correct.

## M-C — the three readings
**(1) Λ_i>0 is NOT universal.** All-positive on the divergent side (λ=0.40, 0.45, 0.48); mixed-sign on the convergent
side (λ≥0.50). So the positivity we see at Collatz is not a structural feature of the whole family.

**(2) The approach direction flips — and Collatz is on the no-turnover side.** `sign(lim Λ)` = the direction T approaches
its limit: **+ (from below) for λ≤0.52, − (from above) for λ≥0.55** — boundary at **λ≈0.53–0.54**. **λ=½ (Collatz) sits
~0.03–0.04 below the flip, on the Λ>0 side.** In the family's structure the physical point approaches its limit from
below with Λ>0 — i.e. **no turnover at λ=½**; a Λ<0 (turnover) regime exists only at λ≳0.54, which is not Collatz. This is
the family-level echo of P6H/P6I's "T_i rising, no turnover."

**(3) No dramatic kink; the critical value is (within uncertainty) the right-limit.** Converged S(λ)=2·limT: λ=0.505→0.470,
0.510→0.468, 0.540→0.446, 0.550→0.445 (min), 0.600→0.452 — non-monotone with a minimum near λ≈0.55. As λ→½⁺ the two
nearest points (0.505, 0.510) head to ~0.47, consistent with the critical `S_∞≈0.477` (P6I) **continuously** — no
evidence of a jump (the earlier apparent kink was coarse λ-spacing + unconverged levels). Caveat: the near-critical zone
(0.51–0.53) is numerically delicate (rate→1, slow convergence), so the right-limit is only pinned to ~0.47±0.01. The left
side (λ<½) diverges, so there is no left-limit to compare.

## Verdict — the family places Collatz just below a turnover boundary it does not reach
λ=½ is confirmed critical (M-D), and the sign-flip of the approach direction sits at λ≈0.53–0.54, with Collatz on the
Λ>0 / no-turnover side ~0.03 away (M-C2). No kink at the critical point (M-C3). This corroborates the depleted-side
reading (P6I/P6J): at the physical point the residual stays positive; the turnover that 7/15 needs lives at λ>0.54, off
the physical point. Model family only — not a statement about Collatz at λ≠½. Not at stake: P6D–P6J, S_{i+1}=2T_i, R1–R30.
