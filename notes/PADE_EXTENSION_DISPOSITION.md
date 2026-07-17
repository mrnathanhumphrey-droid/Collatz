# PADE_EXTENSION_DISPOSITION — Phase 4 synthesis

**Date:** 2026-05-12. Wilson, follow-up probe to R77.6.

---

## DISPOSITION: **H_AMBIGUOUS / INCONCLUSIVE** (with two sub-findings)

> **No new Padé approximants are computable from ε_n through k=6** (the m+n ≤ 4 budget is already exhausted by R77.6's eight tabulated approximants). The probe re-examines those for full pole structure and adds the ratio diagnostic |ε_n|/|ε_{n-1}|. Both diagnostics CONVERGE on the same reading and CLOSE TWO of the pre-registered hypotheses while leaving the residual branch-cut order open:
>
> **SUB-FINDING 1 — H_COMPLEX_SECONDARY is REJECTED.** Every Padé primary pole (closest-to-z=2 across (1,1)..(3,1)) is on the positive real axis. No complex-conjugate primary structure appears across approximants. The only off-axis poles are in [0/4]'s all-pole approximant (two cc-pairs on |z|≈1.10 and |z|≈1.22), which are classical artifacts of forcing the numerator to a constant for a function with branch structure. The ratio diagnostic (Phase 3) shows NO oscillation pattern — instead a clean sign change from + at n=3 to − at n=4,5,6 with monotone-descending deviation. A complex secondary singularity at z = re^{±iθ} would imprint as oscillating ratios; nothing of the kind is observed.
>
> **SUB-FINDING 2 — H_PURE_SIMPLE_POLE is REJECTED.** Padé diagonal still shows drift (2.0764 → 2.0513), not stability at z=2. Ratio diagnostic shows ratios LEAVING 0.5 in absolute terms (departure 0.018 → 0.030 → 0.068 across n=4..6), accelerating downward. A pure simple pole would give ratios converging to 0.5 with geometric speed; we see acceleration AWAY from 0.5. Both diagnostics rule out pure rate-1/2 simple-pole structure.
>
> **RESIDUAL: H_BRANCH_CUT_ORDER_DETERMINED remains undetermined.** The diagonal Padé sequence has only two points ([1/1], [2/2]) — not enough to fit a branch-cut convergence law (needs ≥3 diagonal points). The ratio diagnostic shows a monotone-descending pattern at n=4..6 but the function form (1/n vs 1/n^α vs log) is not separable from 3 data points. **H_SECONDARY_SINGULARITY_LOCATED is NOT supported** — no specific secondary location can be named from this data; the secondary is consistent with the branch cut at z=2 itself rather than a distinct discrete singularity at z = 1/ρ for some ρ < 1/2.

---

## What the two probes say jointly

The Padé diagnostic and the ratio diagnostic measure DIFFERENT aspects of E(z):
- **Padé** localizes the geometric position of the dominant singularity in the complex plane.
- **Ratio** measures the asymptotic rate of |ε_n| as n grows; sensitive to ALL singularities at distance ρ from origin via their combined contribution at each n.

For the readings to be **consistent**, both must agree on:
1. Position of dominant singularity (z=2 — both agree).
2. Whether dominant singularity is a simple pole (both reject) or branch-cut (both consistent with).
3. Whether a complex secondary singularity exists (both reject).

These are all consistent.

The DIVERGENCE (or what looks like one):
- Padé says diagonals are monotone-descending toward z=2 (a "well-behaved" branch-cut convergence).
- Ratio diagnostic at n=4..6 says ratios are monotone-descending AWAY from 0.5 (NOT well-behaved if interpreted as approach to simple pole).

These are NOT contradictory once you read them correctly: the Padé monotonicity is about the diagonal sequence as m+n grows; the ratio descent is about a single asymptotic relation at increasing n. Both are saying the singularity is NOT a clean simple pole, but they don't yet name what it IS.

---

## Hypothesis-by-hypothesis disposition

- **H_SECONDARY_SINGULARITY_LOCATED**: REJECTED. No consistent secondary location identified. The [2/2] secondary at z=0.69 is artifact; [0/4]'s cc-pairs are artifacts. No Padé approximant clusters a secondary pole at any consistent location near or off z=2.

- **H_BRANCH_CUT_ORDER_DETERMINED**: NOT RESOLVED. The branch-cut nature is established (per R77.6); the order α is not determinable from N=5 data. The ratio diagnostic's monotone-descending |r_n − 0.5| AT n=4..6 is consistent with a power-law branch but does not pin α.

- **H_COMPLEX_SECONDARY**: REJECTED. No complex-conjugate primary structure; no oscillation in ratios.

- **H_PURE_SIMPLE_POLE**: REJECTED. R77.6's reading is reinforced — diagonals still drift, ratios are leaving 0.5 not approaching it geometrically.

- **H_IRREGULAR_STRUCTURE**: NOT NEEDED. The pattern is stable enough to rule out two hypotheses cleanly. Not labeled irregular.

- **INCONCLUSIVE (default)**: CHOSEN. The data closes two of the six pre-registered options and leaves the residual branch-cut-order question open. Additional ε_n needed (Route A in R77.6 / TAUBERIAN_SCOPING_DISPOSITION.md).

---

## Connection to Approach 1's δ_n diagnostic

A parallel agent is running Approach 1 (DELTA_DIAGNOSTIC_*), which directly probes the subleading correction δ_n = |ε_n|·2^n − 1/30 for multi-term fit structure (per TAUBERIAN_SCOPING_VERIFICATION.md Phase 4 finding: δ_n is non-monotone and sign-flips between n=5 and n=6).

This probe's ratio diagnostic and Approach 1's δ_n diagnostic are mathematically related but not identical:
- δ_n probes the **value** of the subleading correction at each n.
- |ε_n|/|ε_{n-1}| probes the **rate** (logarithmic derivative) of |ε_n|.

Both should converge on the same singularity structure. Predicted alignment:
- δ_n sign flip at n=5→6 (TAUBERIAN_SCOPING_VERIFICATION.md): corresponds in ratio space to the inflection at n=4 (where r_n crosses 0.5 from above to below) and the subsequent acceleration of departure downward through n=5, n=6.
- Both readings point to: leading rate-1/2 with NEGATIVE subleading correction that grows in magnitude with n — exactly the signature of a branch-cut at z=2 with finite-N convergence to the asymptotic rate.

If Approach 1 finds a clean multi-term ansatz (e.g., A·(1/2)^n + B·(1/2)^n·n^{-α} with specific A, B, α), the ratio diagnostic predicts: ratios r_n at large n should approach 0.5 from below with a power-law correction matching α. The current n=4..6 trend (accelerating departure from 0.5) is consistent with α ∈ (0, 1) — i.e., a relatively slow branch-cut correction. **The two probes should converge on this.**

If Approach 1 finds a complex/oscillatory ansatz, the ratio diagnostic would DISCONFIRM (the sign pattern + - - - rules out oscillation). This is a **non-trivial cross-check**.

If Approach 1 also lands at INCONCLUSIVE / multi-term ambiguous, both probes agree that N=5 is the limiting factor. **Route A (compute ε_7) is required to advance.**

---

## Files
- [PADE_EXTENSION_TABLE.md](PADE_EXTENSION_TABLE.md) — Phase 1 (extended pole table with full lists)
- [PADE_EXTENSION_PATTERN.md](PADE_EXTENSION_PATTERN.md) — Phase 2 (pattern across approximants)
- [PADE_EXTENSION_RATIOS.md](PADE_EXTENSION_RATIOS.md) — Phase 3 (ratio diagnostic)
- [PADE_EXTENSION_DISPOSITION.md](PADE_EXTENSION_DISPOSITION.md) — this file (Phase 4 synthesis)
- [pade_extension.py](pade_extension.py) — main-thread verification script

## Connection to existing scoping documents
- [TAUBERIAN_SCOPING_DISPOSITION.md](TAUBERIAN_SCOPING_DISPOSITION.md) — baseline disposition (H_AMBIGUOUS, recommends Route A)
- [TAUBERIAN_SCOPING_VERIFICATION.md](TAUBERIAN_SCOPING_VERIFICATION.md) — δ_n non-monotone, sign-flip at n=5→6
- [result_77_6_generating_function.md](result_77_6_generating_function.md) — R77.6 (G-branch-cut, type indeterminate)
