# RESULT — SOLSTICE: the drift is transient convergence to a single mode ρ≈0.908 — turnover NOT corroborated (2026-07-26)

**Probe:** `probe_solstice.py`. Does the deparitied-rate drift persist (→ turnover) or flatten? Two-mode model
`Λ_i = A ρ₁ⁱ − B ρ₂ⁱ` (ρ₂>ρ₁) ⟹ rate falls monotonically to a crossing. Certified base-2/ν_e machinery only.

## S-A — gate passes; extension walls at i=16
`Λ₁₂..₁₅ = 0.33677, 0.31971, 0.28672, 0.26193 ×10⁻³` reproduced to **1e-9** (bit-for-bit); `S₁₆ = 0.471352`; `Λ₁₆ =
+0.00023426`. **Precision floor ≈ 3¹⁵·1e-16 ≈ 1.4e-9**, far below the ~3e-3 rate signal — precision is not the limit.
**Depth is:** `stationary_trunc` needs **41 GiB at n=17**, so i=16 is the exact-method wall. i=17,18 require more RAM (Lambda).

## S-B — decider: drift present, even-i clears 3σ, classes agree in sign
| class | rates (i) | drift/level | SE | \|slope/SE\| |
|---|---|---|---|---|
| even-i | 0.9285(12), 0.9227(14), 0.9039(16) | **−0.00616** | 0.00188 | **3.28** |
| odd-i | 0.9097(13), 0.9051(15) | −0.00226 | — (2 pts) | — |
Row-1 (drift persists, \|slope/SE\|>3) fires for even-i; odd-i weaker but same sign. **Not** row-3 (no sign disagreement).

## S-C — curvature: even-i rate accelerating down (−0.013), consistent with two-mode… but see S-E
## S-D — mirror FAILS: enriched Λ⁺ converged flat (\|Λ⁺\|<1e-6 by i≈13), depleted still falling
A genuine measure-global slow mode must show in both hemispheres. It does not — **the slow rate (0.91) is depleted-only.**
Per the pre-registration this **falsifies the clean two-mode / global-slow-mode reading.**

## S-E — the gate (held-out validated): SINGLE mode, no crossing
4-parameter free fit `Λ = A ρ₁ⁱ − B ρ₂ⁱ` on i=10..16 **collapses to one mode: ρ₁ = 0.908, ρ₂ → 0** (the second term is a
decaying spike, not a growing mode). **Held-out** (fit i=10..14, predict Λ₁₅,Λ₁₆ before the date): predicts +2.60e-4,
+2.36e-4 vs actual +2.62e-4, +2.34e-4 — **0.7% / 0.8% miss, PASSES.** A single positive geometric ρ=0.908 **never crosses
zero**: no turnover, no date. (Contrast #45: that fit mispredicted the held-out; this one nails it — for a single mode.)

## Verdict — which S-B outcome fired, honestly: the drift is real but it is CONVERGENCE, not a fall to a crossing
The raw even-i drift clears 3σ (S-B row-1), but the two gating tests built to validate the turnover interpretation both
reject it: **S-E** finds one geometric mode ρ≈0.908 with no second/growing mode and a clean held-out prediction (no
crossing), and **S-D** finds no mirror in the enriched hemisphere. The measured "falling drift" is the deparitied rate
**settling onto ρ₁≈0.908 from above** (0.9285→0.9227→0.9039, → 0.908), i.e. transient convergence — which will *flatten*
at 0.908, not continue to zero. That gives `S_∞ ≈ 2·(T_∞) ≈ 0.476` with **no turnover** — 7/15's only positive indicator
(the raw drift) does not survive the model-free fit. **The i=17,18 extension directly discriminates:** rate flattening at
≈0.908 confirms single-mode/no-turnover; rate falling below 0.908 (against the S-E fit) would resurrect the two-mode date
(~i≈30). Blocked locally by the 41 GiB wall → **running on Lambda.** Value unaffected: `S_∞ ∈ [0.4714, 0.478]` stands.
Not at stake: P6D–P6K identities, S_{i+1}=2T_i, R1–R30.
