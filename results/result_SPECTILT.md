# Probe SPECTILT — placement of the non-3-divisible power spectrum (Wilson's spectral identity) — **the sign lives in a nearly-N/4-symmetric nonnegative spectrum: `Re δ̂_r(1) = ⟨cos(2πm/N)⟩` over the probability measure `Ĉ(m)/S` on 3∤m. Item 2: the mass placement is FROZEN (deciles identical in |m|/N across r=12..16, median ≈ N/4). Item 3: the +margin is a residue of two LARGE opposing cos-bands (±0.319) at 0.61% and TIGHTENING (0.93%→0.61%, r12→16) — the mass split →0.5/0.5, a 0.26% low-freq excess. This is the residue-of-opposing-bands case (the harder proof), the obstruction relocated a THIRD time: Ĉ≥0 bounds mass, not placement. Item 1: my margin-rate note (~0.90) was WRONG — the margin decays ~0.848 late; the identity NumRate=d1Rate·DenRate closes per-level (0.900·0.942=0.848). Wilson's flat-self-spike mechanism is FALSE: a flat Ĉ gives Num=0 exactly (Σ_{3∤m}cos=0).**

**Date:** 2026-07-25. Probe `probes/probe_spectilt.py`, log `logs/spectilt_run.log`. No build_nu: Item 1 from the gated AC-B lags (r=2..16); Items 2/3 from the saved spectra `scratchpad/rho_12..16.npy`. Wilson's identity (replica-peak variables): `Re δ̂_r(1) = [Σ_{3∤m} Ĉ(m)cos(2πm/N)] / [Σ_{3∤m} Ĉ(m)]`, `Ĉ(m)=|ρ̂(m)|² = prof[m] ≥ 0` (power spectrum, even) — **the cos-weighted mean of the nonnegative fluctuation spectrum on 3∤m, a probability measure.** Denominator `= 2[C(0)−C(N/3)] =` total fluctuation spectral mass; its vanishing = forced-uniformity (R14), not cancellation.

## Item 1 — rates and the identity (the discrepancy was the margin rate)
`Den=C(0)−C(N/3)`, `Num=2C(1)−C(N/3−1)−C(N/3+1)`, `d1=Num/(2·Den)` (all /C0). Geometric rates:
| window | Den | Num (margin) | d1 | identity: d1·Den |
|---|---|---|---|---|
| full-span r2→16 | 0.904 | 0.710 | 0.785 | 0.710 ✓ |
| **late r12→16** | **0.942** | **0.848** | **0.900** | 0.900·0.942 = **0.848** ✓ |

All three rates are level-dependent and **slowing** (Den 0.80→0.947 across levels). **`Re δ̂(1) = Num/(2·Den)` is an identity, so `rate(Num)=rate(d1)·rate(Den)`, checked per-level.** The AC-LAGS note's "margin decays ~0.90" was **wrong** — the late margin rate is **0.848**; Wilson's 0.81 estimate used the full-span Den 0.904 rather than the late 0.942. Resolved: the margin rate is the object that was mis-stated.

## Item 2 — mass placement is FROZEN, not migrating
Deciles of `Σ_{3∤m}Ĉ(m)` in `|m|/N` (=min(m,N−m)/N; N/4 → 0.25), r=12..16 **identical**:
`10%→0.050, 25%→0.124, 50%→0.249, 75%→0.374, 90%→0.450`.
Self-similar in rescaled frequency; **median ≈ N/4** (the cos zero-crossing). The measure's shape is stable across the levels — nothing migrates outward. What shrinks is the asymmetry, not the placement.

## Item 3 — the +margin is a residue of two large opposing bands (the hard case, tightening)
Numerator by cos-sign band (contribution to ⟨cos⟩ = d1):
| r | pos (\|m\|<N/4) | neg (\|m\|>N/4) | residue = d1 | \|residue\|/\|pos\| | mass(\|m\|<N/4) |
|---|---|---|---|---|---|
| 12 | +0.31972 | −0.31676 | +2.964e−3 | 0.0093 | 0.5022 |
| 14 | +0.31951 | −0.31707 | +2.441e−3 | 0.0076 | 0.5016 |
| 16 | +0.31931 | −0.31737 | +1.939e−3 | **0.0061** | 0.5013 |

The two cos-weighted half-contributions (±0.319) are each **~165× the residue**, and the balance is **tightening** (residue 0.93%→0.61% of a band, r12→16). The mass split (0.5022→0.5013) → 0.5/0.5. **So the measure is nearly symmetric about N/4, and d1 is a ~0.26% low-frequency mass excess.** This is the **residue-of-opposing-bands** case (Wilson's two strategies: this is the harder one — proving a whisker-thin excess of two O(1) quantities), NOT a genuine tilt of a spread measure.

## The flat-self-spike mechanism is FALSE
A flat component `Ĉ(m)=c` (all m) gives `Num = c·Σ_{3∤m}cos(2πm/N) = c·(−1 − (−1)) = 0` exactly — **the k=0 self-spike is annihilated by the 3∤m projector**, not positively tilted. Wilson's conjectured mechanism (flat → positive tilt on 3∤m) does not hold. The tilt lives entirely in the **non-flat** structure of Ĉ on 3∤m — the 0.26% low-frequency mass excess.

## The obstruction, third time, now fully concrete
The pen must prove **`⟨cos(2πm/N)⟩ > 0`** over the probability measure `Ĉ(m)/S` on 3∤m — a measure that is nearly symmetric about N/4 and, under `Ĉ≥0` alone, **free to be exactly symmetric** (⟨cos⟩=0). Positivity bounds the total mass; it does not forbid the symmetric placement. The spectral chart did not dissolve the near-cancellation — it **relocated** it from the divergent G−K split into the two cos-bands and made it a single measured distribution, whose low-side excess (0.26% mass, 0.61%-of-band cos residue, both shrinking) is the entire sign question.

## Status
**SPECTILT (measurement; sharpens, does not decide):** Wilson's identity `Re δ̂_r(1)=⟨cos(2πm/N)⟩_{Ĉ/S, 3∤m}` confirmed (d1_spec matches RATIO-2 d1 to 5 digits at r=12..16). **Item 1:** identity closes per-level; the AC-LAGS margin-rate note (~0.90) was wrong — late margin rate 0.848 = d1(0.900)·Den(0.942); Den rate slows 0.80→0.947. **Item 2:** mass placement FROZEN (deciles identical in |m|/N, median ≈ N/4) — stable, not migrating. **Item 3:** the +margin is a **residue of two large opposing cos-bands** (±0.319, ~165:1), TIGHTENING (0.93%→0.61%, r12→16); measure nearly symmetric about N/4 (mass →0.5/0.5), d1 = 0.26% low-freq excess — the **harder** proof strategy. **Flat-self-spike mechanism FALSE** (flat Ĉ → Num=0). **Obstruction recurs a THIRD time: Ĉ≥0 constrains mass, not placement** — a symmetric-about-N/4 spectrum is admissible and gives ⟨cos⟩=0. Concrete target for the pen: ⟨cos(2πm/N)⟩>0 on the 3∤m spectrum, i.e. a strict low-freq mass excess. AC-E "A_r ~ RATIO-2 oscillation" left as a named conjecture, not a finding (Wilson's flag). Not at stake: R1–R30, R80–R82, d1 ladder to r=16, MODES' Re δ̂(1)>0, RATIO-2's subdominance, AC-LAGS' inequality (now in spectral form). rho_12..16 + shape grids remain in scratchpad.
