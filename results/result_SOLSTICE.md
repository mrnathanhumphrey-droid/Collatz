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

## EXTENSION i=17,18 (Lambda, matrix-free GPU) — the local single-mode reading is FALSIFIED
`T_17 = 0.23611673, T_18 = 0.23629630` (matrix-free power iteration on A100, torch; `build_base2`/`stationary_trunc`
refactored matrix-free, gated locally vs certified T to 3e-9; on-instance gate `T_15,T_16` diff 1.7e-9, 2.3e-9 PASS).
`Λ_17 = +2.0665e-4, Λ_18 = +1.7957e-4` — **still positive, still decreasing.**
- **Deparitied rate fell BELOW the local ρ=0.908 and the drift STEEPENED.** even-i: 0.9285(12), 0.9227(14), 0.9039(16),
  **0.8755(18)** — successive diffs −0.0058, −0.0188, −0.0284 (**accelerating**); drift −0.0118/level, **|slope/SE| = 8.53**
  (was 3.28 at i≤16); curvature −0.0096. odd-i: 0.9051(15), **0.8882(17)**.
- **S-E single-mode fit FALSIFIED out-of-sample:** the ρ=0.908 fit on i=10..16 predicts `Λ_18 = +1.947e-4`, actual
  `+1.796e-4` — **8.4% miss.** The "settling at 0.908" reading of the local probe is wrong; the true dominant decay is
  lower.
- **Best two-mode fit (i=12..18, held-out validated):** `Λ_i = A·0.867ⁱ − B·0.628ⁱ` (A=2.44e-3, B=2.72e-2, resid
  1.9e-6); ρ₂<ρ₁ ⟹ **no crossing** (two *decaying* modes, dominant ρ₁≈0.867), held-out (fit 12..16 → predict 17,18)
  within 1.4%. BUT the data's late acceleration is marginally steeper than this no-crossing model, so a crossing
  (ρ₂>ρ₁) is **not excluded** — only not preferred.

## REVISED VERDICT — no turnover observed; value pulled down to ≈0.475; sign still unresolved at the knife's edge
The extension **kills the clean local reading** (rate ≠ 0.908; it is ~0.87 and still falling). But `Λ` remains **positive**
through i=18 — **no turnover has occurred**, and 7/15 needs `Σ_{i≥19}Λ_i = −0.00296` (a sustained negative run) against a
current `Λ_18 = +1.8e-4`. **Value revised down:** exact floor `2·T_18 = 0.47259`; geometric tail at rate ≈0.867–0.88 gives
`S_∞ ≈ 0.4749–0.4752` (down from ≈0.476 — the faster decay shrinks the tail), landing just above the `<0.475` headline.
**Sign:** the accelerating rate-fall is best explained by a second *decaying* mode (ρ₂≈0.63), not a growing one — the LS
optimum has no crossing — but with 7 points and a late steepening, a crossing is not ruled out. **i=19,20 (Lambda) would
separate asymptote-at-≈0.87 (no turnover) from continued fall (turnover).** Which S-B outcome fired: **row-1 (drift
persists, |slope/SE|=8.53), reversing the i≤16 single-mode reading** — but the two-mode date is *not* robustly
determined (S-E prefers no-crossing). Not at stake: P6D–P6K identities, S_{i+1}=2T_i, the value floor 2·T_18=0.4726, R1–R30.
