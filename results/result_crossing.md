# Probe CROSSING — direct Λ_r past the ε-data (r=12..16 via ν) — **two-monotone-difference model KILLED; the bore is a single LONG positive lobe, period ≳ 26; and this explains why every period fit failed (the data hasn't completed even a half-cycle)**

**Date:** 2026-07-23. Probe `probes/probe_crossing.py`. Key capability: `Λ_r = Σ_m 4^{−m} A_r(m)` with
`A_r(m) = γ_r(τ_m) − γ_{r−1}(τ_m)` is computable **directly from the ν measure (build_nu) at any r**, NOT limited by
the ε-data (which stopped at ε₁₂ ⟹ Λ₁₁). So this reaches **Λ₁₂…Λ₁₆ — five levels past the ε ladder.** Tests
Wilson's two-monotone-mode model `Λ_r ≈ aρ₁^r − bρ₂^r` (a difference of decaying exponentials crosses zero exactly
once and stays negative) against the long-period model.

**Validation:** ν-computed Λ matches the ε-telescoping Λ to **0.00%** at r=3..11 — the r≥12 values are trustworthy.

## The crossing — Λ_r stays POSITIVE through r=16
| r | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|---|
| Λ_r (×10⁻⁴) | +3.69 | +3.64 | +3.91 | +3.86 | **+3.37** | +3.20 | +2.87 | +2.62 | +2.34 |
| source | ε+ν | ε+ν | ε+ν | ε+ν | **ν only** | ν | ν | ν | ν |

**Λ is positive and monotonically declining from its peak (~r=10) all the way to r=16 — no sign change.**

## Two-monotone-difference model: **KILLED**
The one-number test: `Λ₁₂ = +3.37×10⁻⁴`, **not** the predicted `T·(1−ρ) = −1.8×10⁻⁵`. There is **no crossing at
r=12**, and none through r=16. A difference of two decaying exponentials would have crossed once and stayed
negative; instead Λ is a slowly-declining positive lobe. The **MONO/OSC split** shows why: `MONO = 0.25·A(1) +
(1/64)·A(3)` is dominated by the **positive-monotone m=1 channel**, and the negative pieces (m=3, and OSC =
`0.0625·A(2)+(1/256)·A(4)`) are far too small to overtake it — `|OSC/MONO| ≤ 0.25` for r≥7, → 0.01 by r=12. MONO
never crosses. The model is dead.

## What the direct extension actually shows: a single LONG positive lobe
- Λ signs: `+(3,4,5)  −(6)  +(7…16)`. The `−` at r=6 is a **single isolated term** — the transient/bore crossover
  (Probe D's r=6 crossover), not a bore half-period.
- The bore is a **single positive lobe spanning r≈7 through past r=16**, still declining and not yet crossed. A lobe
  ≥10 levels wide, un-crossed at r=16, means the half-period is **≳ 13–18**, so the **period is ≳ 26** — a *lower
  bound* that is **longer than every prior point estimate** (9, 17.6, 22) and consistent with R26's ">16" and the
  retired float-era ~37.
- `|Λ_r|/0.984^r` declines 0.46 → 0.30 over r=10..16 — i.e. Λ decays *faster* than 0.984^r locally, the signature
  of a cosine heading toward its zero-crossing (somewhere r≈22–28), not a change in the envelope rate.

## Why every period fit failed — resolved
Reconstructed ε keeps **rising**: `ε₁₂ = 2.27e−3 → ε₁₆ ≈ 4.7e−3` (via `ε_{r+1}=ε_r+2Λ_r`). The measured S_k is still
moving *away* from 7/15 at k=16 — **we are still inside the rising half of the first bore lobe, pre-peak.** There has
been **no completed bore oscillation anywhere in the data** (the first bore crossing hasn't happened by r=16). That
is the real reason Probe D, SAT, and MOON all failed to pin a period: **you cannot fit a frequency to data that
hasn't completed even a half-cycle.** The wall is not just "window < 1 period" — it is "window < ½ period."

## Is 7/15 threatened? No.
For `S_∞ = 7/15`, ε must turn over (Λ crosses negative ~r≈22–28) and return to 0 — a large, slow excursion (peak
ε ~ +5×10⁻³ near k≈21). The amplitude **is** decaying (`|Λ_r|/0.984^r` declining, local rate ~0.91 < 1), so the
excursion will turn over — **consistent with 7/15**, just via a big slow oscillation whose rising half is all we've
seen. No alarm on the constant; the convergence is simply much slower and longer-period than the low-r transient
suggested.

## Status
**CROSSING: two-monotone-difference model killed; bore is a single long (period ≳ 26) positive lobe; the
period-fit failures are explained as "less than ½ a bore period in the data."** Wilson's crossing prediction
(`Λ₁₂ ≈ −1.8e−5`, negative) is **refuted** — `Λ₁₂ = +3.37e−4`, positive, and Λ stays positive through r=16 (ν
validated to 0.00% vs ε at r≤11). The MONO channel (m=1 dominant, positive) never crosses; OSC is negligible. The
r=6 `−` is the transient crossover, not a bore zero. **New results:** (1) Λ is directly computable past the ε-data
via ν, reaching r=16; (2) the bore is a single un-crossed positive lobe ⟹ **period ≳ 26**, longer than all prior
estimates, consistent with the retired ~37; (3) ε is still rising (pre-peak) at k=16, so **no completed bore
oscillation exists in the data** — the true reason every period probe (D/SAT/MOON) railed. `S_∞ = 7/15` is not
threatened (amplitude decays, local rate ~0.91 < 1). No fitting; ν-vs-ε validation at 0.00%; the crossing extension
reported as direct measurement, the period as a bounded inference (≳26), the 7/15 non-threat argued from the
decaying envelope. Not at stake: R1–R30, R80–R82, Probe D, SAT, MOON — all consistent, now explained.
