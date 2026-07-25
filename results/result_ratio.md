# Probe RATIO — the successive-ratio sequence of the dominant mode Re δ̂_r(1) — **the ratio ρ_r = Re δ̂_{r+1}(1)/Re δ̂_r(1) is NOT monotone: it turns DOWN at r=6 and r=10 (gap 4), sign pattern of d(ρ) = `+ + + − + + + − +`, four flips. By Wilson's pre-registered rule this is OSCILLATION ⟹ a complex pair in the n=1 mode is NOT excluded ⟹ the monotone-rise conditional for `Re δ̂(1)>0` is VOID and the 7/15 branch is live again. This is TSW's real-vs-complex question asked of the one object we can compute exactly, and it answers "not provably real."**

**Date:** 2026-07-25. Probe `probes/probe_ratio.py`. Reuses the validated `delta_from_nu` (build_nu→dlog→|FFT|², r=2..12) and reads `Re δ̂_r(1) = FFT(δ_r)[1].real`. Gates Wilson's pre-registered conditional: *if* `Re δ̂_r(1)=Σ_i c_i ρ_i^r` with finitely many **real** ρ_i>0, and the positive sequence has successive ratios **rising monotonically** to a limit, then the largest-ρ component carries c>0 (a negative dominant coefficient forces the sequence through zero, at which point ratios **fall** to zero and go negative rather than rising) ⟹ no approaching sign change of the dominant mode. **A complex pair shows up as OSCILLATION in ρ_r, not monotone rise** — the pre-registered falsifier.

## Re δ̂_r(1) — the dominant fluctuation mode (all eleven levels)
Positive for every r=2..12 (the single-signed dominant mode from MODES holds):
`+5.714e−2, +1.983e−2, +8.643e−3, +8.039e−3, +7.742e−3, +5.787e−3, +4.666e−3, +3.895e−3, +3.601e−3, +3.211e−3, +2.964e−3` (r=2..12).

## READ 1&2 — ρ_r and its successive change (the decisive table)
| r | ρ_r = d̂_{r+1}(1)/d̂_r(1) | d(ρ) | |
|---|---|---|---|
| 2 | 0.3470 | — | |
| 3 | 0.4359 | +0.0889 | RISE |
| 4 | 0.9301 | +0.4942 | RISE |
| 5 | 0.9631 | +0.0330 | RISE |
| **6** | **0.7474** | **−0.2157** | **DOWNTURN** |
| 7 | 0.8063 | +0.0589 | RISE |
| 8 | 0.8347 | +0.0284 | RISE |
| 9 | 0.9245 | +0.0898 | RISE |
| **10** | **0.8918** | **−0.0327** | **DOWNTURN** |
| 11 | 0.9229 | +0.0311 | RISE |

- **Sign pattern of d(ρ): `+ + + − + + + − +` — four direction flips.**
- **Downturns at r=6 and r=10 (gap 4).** The monitor does NOT stay rising.
- Per Wilson's rule: oscillation ⟹ **complex pair not excluded ⟹ conditional VOID ⟹ 7/15 live.**

**The downturns sit on the lengthening-crossing cadence.** r=6 is exactly the n=3 g_r crossing (from MODES); the next ρ-dip is at r=10, gap 4. The oscillation is **intrinsic to the n=1 coefficient's own ratio** — it is NOT the n=2/n=3 excursions bleeding in (those are separate Fourier coefficients). So the real-spectrum assumption Wilson's proof rests on **fails the test on its own object**: Re δ̂(1) does not decay by a sum of pure real geometric modes; there is a ~period-4 oscillation in its ratio.

## READ 3 — the limit ρ_r curves toward
- Raw ρ tail (r=2..11): `0.3470 0.4359 0.9301 0.9631 0.7474 0.8063 0.8347 0.9245 0.8918 0.9229`.
- Aitken Δ² acceleration: `0.3275 0.9654 0.9345 0.7937 0.8611 0.7932 0.9005 0.9077` — **wobbles (does not settle to a clean fixed point), last value ≈0.908.**
- **Curves toward ~0.90–0.91, NOT the banked ~0.984** within the observed window. Tail-sum factor `1/(1−ρ)` is ≈10.75 at 0.907 (vs 62.5 at 0.984) — so the 0.477 extrapolation is **not** recomputed at a much slower rate here; the value story does not unify to 0.984 by r=12.

## What this settles and what it does not
- **The clean proof Wilson hoped for does NOT certify.** The conditional required a monotone-rising ρ_r; ρ_r is not monotone (downturns at r=6, r=10). So "rising ratios ⟹ no approaching sign change" cannot be invoked — **the antecedent is false in the data.**
- **A complex pair (period ~4, φ~π/2) is not excluded** — indeed the downturn cadence (gap 4) is exactly its signature. If real, the next ρ-dip is due ~r=14, unobserved past the compute wall.
- **#30 caution, both ways.** Two downturns = one gap. A genuine period-4 complex pair is indistinguishable here from two transient dips of the same lengthening family we already see in n=2/n=3. Either reading leaves the sign of the deep tail undetermined: complex pair ⟹ g_r oscillates ⟹ 7/15 route alive; two transients ⟹ eventual monotone rise still possible ⟹ 0.477.
- **Net for the arc:** MODES' `Re δ̂(1)>0` for r=2..12 is unshaken (all eleven values positive). But the RATIO test — the exact-computable form of TSW's real-vs-complex question, asked of a scalar sequence with no operator to construct — comes back **oscillatory, not monotone.** So `Re δ̂(1)>0` is **observed to r=12, NOT proved to persist**; the monotone-rise route to proving it is closed. The honest surviving statement is unchanged from MODES but now with the proof-route ruled out: **leans 0.477 (dominant mode positive, envelope of the n=2/n=3 spikes decreasing), does NOT establish it; 7/15 requires a spike — or a genuine n=1 complex-pair dip — that clears threshold past r≈14–16, unobserved.**

## Status
**RATIO:** `Re δ̂_r(1)>0` all r=2..12 (single-signed dominant mode confirmed at all eleven levels). But `ρ_r = Re δ̂_{r+1}(1)/Re δ̂_r(1)` is **NOT monotone** — downturns at r=6 and r=10 (gap 4), d(ρ) sign pattern `+ + + − + + + − +`, four flips = **OSCILLATION** by Wilson's pre-registered rule ⟹ a **complex pair in the n=1 mode is not excluded** ⟹ the monotone-rise conditional proof of `Re δ̂(1)>0` is **VOID**, 7/15 branch live. Aitken limit ~0.90–0.91 (not the banked 0.984). The ratio oscillation is intrinsic to the n=1 coefficient (not n=2/n=3 bleed) and rides the lengthening-crossing cadence (dip at the r=6 crossing, next at r=10). #30 both ways: 2 downturns = 1 gap, complex-pair-period-4 indistinguishable from two transient dips ⟹ deep-tail sign undetermined. **This is TSW's real-vs-complex question asked of the one exactly-computable object, and it answers "not provably real."** Leans 0.477 (unchanged), does NOT establish it; the monotone-rise proof route is now closed. Not at stake: R1–R30, R80–R82, the g_r ladder to r=12, Λ^unif closed form, MODES' `Re δ̂(1)>0` r=2..12. Corrects nothing prior — sharpens MODES: the surviving proof burden (an analytic lower bound on Re δ̂_r(1) from the ν_r construction) cannot be discharged by the ratio-monotonicity shortcut.
