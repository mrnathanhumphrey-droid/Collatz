# Probe R27 — the recurrence fit — **jackpot dead; |λ₂|≈0.5 (real) reinforces gap-survival; ≥3 modes, period-9 demoted**

**Date:** 2026-07-22  Exact-rational (A/B/D) + float subcritical (C). Probe `probes/probe_recurrence_R27.py`.
Wilson's amendment: fit the recurrence on **Λ_r** (→0, carries only subdominant modes — the difference kills the
eigenvalue-1 constant that gave R26 its spurious 1.12), on the **exact rationals** Λ₁…Λ₇ (no drift possible).
Λ_r = (ε_{r+1}−ε_r)/2, signs −,−,+,+,+,−,+.

## R27-A — TWO-TERM RECURRENCE, EXACT: **no exact recurrence; clean solves give |λ₂|≈0.5 real**
`Λ_{r+1} = a·Λ_r + b·Λ_{r−1}`, four exact solves from consecutive pairs:

| solve | a | b | root \|λ₂\| | note |
|---|---|---|---|---|
| {2,3} | +0.2919 | −0.0363 | 0.190 (complex, **period 9.010**) | transient-contaminated (Λ₁,Λ₂ large) |
| {3,4} | +0.5001 | +0.0013 | **0.503 (real)** | clean window (r=3–6) |
| {4,5} | −76.2 | +37.8 | 0.493 | **near-singular** (Λ₄²≈Λ₃Λ₅, det≈0) — garbage a,b |
| {5,6} | −0.766 | −0.136 | ≈0.49 | clean-ish (r=4–7) |

**The four solves do NOT agree** (a spans −76→+0.5) — so there is **no exact finite 2-term recurrence**, and the
jackpot branch (rational `L(z)=ΣΛ_r z^r`, closed-form ΣΛ) is **dead**. But the well-conditioned solves converge on a
**dominant subdominant `|λ₂| ≈ 0.5`, real** ({3,4}=0.503, {5,6}≈0.49) — matching **R18-A's exact ratio 0.493/0.503**
exactly. A *real* ≈0.5 mode is far from the unit circle: the "complex pair near the circle shuts the gap" picture is
not the right one. **Period-9 appears** (solve {2,3}: 9.010) but as a *secondary* complex mode (|λ|≈0.19 < 0.5) from
a transient-contaminated fit — so period-9 is **demoted, not settled** (against the pre-registered expectation of a
clean miss, it survives, but only as a minor mode).

## R27-B — ORDER TEST: **≥3 modes** (not an exact 3-term recurrence either)
3-term solve from r=3,4,5 predicts Λ₇ with residual **−7.3e−4** (≈3× Λ₇ itself) — nonzero, and large. The 2-term
{2,3} residuals predicting Λ₅,₆,₇ are +1.9e−4, −4.1e−4, +3.3e−4 (same order as the Λ's). So the Λ sequence over
r=1…7 is **not** a clean 2- or 3-mode signal — ≥3 significant modes, and the short transient-heavy window defeats
clean modal extraction. (This is why R27-A doesn't converge.)

## R27-D — SIGN-PATTERN CROSS-CHECK: 2-mode model **fails** ⟹ third mode non-negligible
Propagating each fitted 2-term recurrence from Λ₁,Λ₂ gives signs `--+++++` (solve{2,3}) / `---+-+-` (solve{4,5}) vs
the banked `--+++-+` — **mismatch both times**. A two-mode model does not reproduce the observed sign pattern, so
the third mode matters (R27-B is the operative reading), consistent with R27-A's non-convergence.

## R27-C — κ = (|λ₂|/ρ)/(2λ²): **not cleanly monotone; subcritical extraction is method-dependent**
2-term fit on subcritical Λ_r (float, deep r): ρ(fit) matches ρ_pred well (leading mode correct), but κ =
1.261, 1.054, 1.025, 0.925 at ε = 0.05, 0.075, 0.08, 0.1 — **overshoots 1** as ε→0 rather than approaching it
cleanly, and the extracted |λ₂|/ρ (0.763 at ε=0.05) **disagrees with R26's Δp-rate (0.588)**. So the subdominant
extraction is method-dependent and noisy subcritically; the 2λ² law's finite-ε correction (κ) is not pinned. (The
disagreement is the ≥3-mode contamination again — different estimators weight the modes differently.)

## Status
**R27: the gatekeeper stays favorable, the jackpot dies, the fine structure is unsettled.** **A** — no exact
recurrence (four solves disagree), so no closed-form ΣΛ; but the well-conditioned solves give **|λ₂| ≈ 0.5, real**
(={3,4}=0.503, {5,6}≈0.49, = R18-A 0.493/0.503) — a third independent line on |λ₂|<1, and *real*, so the dominant
subdominant is **not** a near-circle complex pair (gap-survival reinforced, the shutting-gap picture refuted more
strongly). **B/D** — ≥3 modes (3-term residual nonzero, 2-mode sign test fails), so the short exact window can't
resolve the full modal structure. **C** — subcritical |λ₂|/ρ is method-dependent (0.59 vs 0.76), κ overshoots 1;
not pinned. **Period-9** appears as a *secondary* mode (|λ|≈0.19, period 9.010) — demoted, neither cleanly confirmed
nor killed.

**Consequence for the crux (owed to the pen):** the gatekeeper conclusion **|λ₂| ≈ 0.5 < 1 (gap survives, route
closes)** is *reinforced* — now anchored on the exact-rational critical Λ recurrence (clean solves) + R18-A + R26,
and with the dominant mode identified as **real ≈0.5** (not a near-circle pair). What R27 does *not* deliver: a clean
recurrence / closed form (jackpot dead), a settled period (period-9 demoted to a minor mode), or a pinned finite-ε
κ. The Λ sequence has ≥3 modes over the reachable exact window, so the *fine* structure needs Wilson's analytic
gap-operator (R27-E) — the exact κ(λ), the mode count, and whether κ(½)=1 — which the 7-point data cannot settle. No
fitting of the leading rate (ρ exact); exact-rational solves reported with their spread and singularities flagged;
the non-convergence reported as non-convergence, not smoothed into a false clean answer.
