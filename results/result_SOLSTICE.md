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

## EXTENSION i=19 (Lambda, tiled matrix-free A100 40GB) — the ambiguous middle; acceleration STALLED; intermittency refuted
`T_19 = 0.23645299645774767` (int32-index squeeze + tiled `index_add_` scatter to fit 40 GB; on-instance gate
`T_15,T_16,T_18` **bit-identical** to bank, diff 1.7e-9/2.3e-9/3.0e-9 PASS — int32 index path proven precision-neutral).
`Λ_19 = +1.56700e-4` — **positive; no turnover.**

**Three PRE-REGISTERED bins (pinned on-disk before the run) for rate(19)=(Λ_19/Λ_17)^½:** >0.888 = intermittency turn-up
(power-law tail i^{−β}, continuum at 1 — the only shape that saves 7/15); 0.871–0.888 = geometric asymptote at ρ₁ (no
crossing); ≤0.871 = crossing lean. **Landed: rate(19) = 0.8708** → the **LINEAR/ambiguous** boundary (predicted linear
continuation 0.8713; miss 0.0005).
- **Intermittency REFUTED one level further.** A power-law tail forces the deparitied rate to *rise* toward 1 (fitted β=1.55
  ⟹ rate rises 0.868→0.913 over i=12..18); measured rate *falls* 0.929→0.871, opposite sign. No turn-up at i=19.
- **The odd-class acceleration STALLED.** Per-level drift: −0.00228 (i15), −0.00845 (i17), **−0.00872 (i19)** — the fall
  stopped accelerating (Δcurv = −0.00026 ≈ 0). Wilson's "accelerating fall ⟹ crossing" indicator **did not fire**; the
  −0.0087/level drift sits on his "fit stands / no-crossing" value (−0.0074), far from "crossing live" (−0.014).
- **1-term tail gate:** Λ_19 / geometric-continuation(0.8882) = 0.15670/0.16304 = **0.961** — slightly below constant-rate
  ("sustained," fall continued), not above ("decelerated"). Consistent with the linear read.
- **Value floor rises, model-free:** `2·T_19 = 0.472906` (supersedes 2·T_18=0.47259). Model-dependent tail (rate ≈0.87,
  UNGATED per Wilson) → `S_∞ ≈ 0.475–0.476`.

**READ (honest):** the fall *continues* (intermittency dead) but its acceleration is *dead too* — landing exactly on linear
continuation. This **weakly favors asymptote/no-crossing** (drift on the fit line, acceleration stalled, crossing-drift
−0.014 not seen), but it is genuinely the ambiguous bin: rate(19) sits only ~0.004 above ρ₁≈0.867, and the **even class was
still accelerating at i=18** (−0.0058,−0.0188,−0.0284) — the two classes now **disagree on curvature**, and no like-for-like
compare exists without Λ_20. **i=21 (odd) separates flatten-at-0.867 (asymptote) from continue-below (crossing); i=20 (even)
tests whether even also stalls.** Both are behind the base-4 dlog-table wall (q=3²¹≈84 GB to materialize; int32 index dies
at n=20 since 3²⁰>2³¹) — a separate CPU-big-RAM or walk-without-DL build, not a bigger-card rerun. Not at stake: P6D–P6K
identities, S_{i+1}=2T_i, the value floor 2·T_19=0.472906, R1–R30.

## EXTENSION i=20 (Lambda, CPU big-RAM parallel) — THE CROSSING IS DEAD; even class decelerates; no turnover
`T_20 = 0.2365885345278137` (parallel matrix-free CPU on a100_sxm4 200 GB RAM: numba `parallel=True` transfer with 6
private accumulators + closed-form coprime-to-q rank `d=idx−idx//q−1` + **uint64** index arithmetic — int64 overflows at
n=20, `base·p→1.2e19 > 9.2e18` — + circular autocorr. Kernel gated locally to **0.00e+00** vs banked T_15; on-instance
gate `T_15,T_16,T_18` diff ≤3e-9 PASS. n=20 ncp=2.32e9, ~20 power-iters × ~920 s; run decoupled from the launcher cap and
monitored to completion, $12.58/6.3 h.) `Λ_20 = +1.35538e-4` — **positive; no turnover.**

**The decider (even class, the one live crossing indicator).** Even 2-step rates 0.9285, 0.9227, 0.9039, 0.8755, **0.8688**;
successive falls **−0.0058, −0.0188, −0.0284, −0.0067** — the acceleration that was steepening toward a crossing **REVERSED
to a 4× deceleration.** Both parity classes now decelerate (odd stalled at i=19, even collapses at i=20).
- **Pre-registered rate(20) = 0.8688** → the **decelerate / no-crossing** band (below the 0.8755 turn-up line, well above the
  0.847 crossing line); settling onto ρ₁≈0.867 from above.
- **Wilson falsifier:** even per-level drift(20) = **−0.00337**, *below* the "fit-stands" −0.0074 and far from "crossing-live"
  −0.014 ⟹ **decisive no-crossing.**
- **Model-free tail gate (gates the tail arithmetic):** geometric-from-Λ_18 predicted Λ_20 = 0.13765e-3, actual 0.13554e-3,
  **ratio 0.985** — the tail arithmetic (and the 47–150× exclusion built on it) HOLDS.
- **Phase-matched (oscillation-free) 8-level trend `(Λ_20/Λ_12)^(1/8)` = 0.8925.**

**VERDICT — no turnover; the crossing is dead; S∞ ≈ 0.475.** New model-free floor `2·T_20 = 0.473177` (supersedes 0.472906);
geometric tail at ≈0.87 (still falling) gives **S∞ ≈ 0.4750**. 7/15 = 0.4667 now requires an unobserved *late sign flip*
(`Σ_{i≥19}Λ = −0.00296`) that **both** parity classes' deceleration actively argues against — the accelerating even-class
fall was 7/15's only surviving indicator, and i=20 killed it. Consistent with the "no crossing, value decelerates to a
non-rational limit" reading (irrational-rotation / quasi-periodic-tail hypothesis). Not at stake: P6D–P6K identities,
S_{i+1}=2T_i, the value floor 2·T_20=0.473177, R1–R30.
