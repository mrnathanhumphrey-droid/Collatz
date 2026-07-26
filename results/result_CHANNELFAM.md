# RESULT — PROBE CHANNELFAM: the channel-constant family and the tail-free determination of S_∞ (2026-07-26)

**Probes:** `probe_channelfam.py` (+ `probe_alpha.py` closing the sign route). Wilson's channel-family reframe:
`γ_r(k)=3^r p_r(k)` (`p_r(k)`=lag-k autocorr of the dlog profile), `A_r(k)=γ_r(k)−γ_{r-1}(k)`, `Λ_r=Σ_{k≥1}4^{-k}A_r(k)`.
Self-check derived last turn: **`Σ_{k≥1} 4^{-k} γ_∞(k) = S_∞/2`**, giving S_∞ from the channel constants — never touching the ε tail.

## STEP 0 — α closes the sign route (nondescript phase, "done looking")
`probe_alpha.py`: the reduction reconstructs `q−⅓` exactly (s=−1 MATCH, +4.1789e-4). But `arg α` is **nondescript** —
dominant cell c=1 (|α|=0.39) at −158.6°, cells scatter −129°…+31°, nowhere near 4π/9 or any small-denominator angle.
By Wilson's own criterion: a number, not a structural constant. Confirms irreducibility; validates the pivot to constants.

## STEP 1 — the normalization gate PASSES (Wilson's flagged step is sound)
The relation reduces (two telescopes of `Λ_r`) to `Σ 4^{-k}γ_∞(k) = 1/3 + (S_∞−S_1)/2`, so **`= S_∞/2` iff `S_1=2/3`** —
sidestepping the Cauchy-transform factor Wilson was unsure of. Measured:
- **(a) EXACT:** `p_1(0)=‖ρ_1‖²=5/9`, `γ_1(0)=5/3`, `S_1=A_1(0)=2/3` ✓. **The relation is `S_∞/2` exactly.**
  Bonus clean seed: `γ_1(0..3) = [5/3, 2/3, 2/3, 5/3]`.
- **(b) bookkeeping** `Σ_{r≤R}Λ_r == Σ_k 4^{-k}(γ_R(k)−1)`: worst |diff| = **0.00e+00**.
- **(c) finite-R** `2·Σ4^{-k}γ_R(k)` vs `S_R`: R=16 → 0.47182 vs 0.47135 (rel 1e-3), converging together.

Convention confirmed (STEP 2): `Re δ̂_2(1)=A(1)/A(0)=2/35=+0.057143` (banked). `γ_r(1)→0.730` (matches CHANNEL_ID).

## STEP 3 — the channel constants (r=16 + deparitied per-channel rate)

| k | 3∣k | γ_16(k) | A_16(k) | rate ρ_k | γ_∞(k) |
|---|-----|---------|---------|----------|--------|
| 1 |     | 0.730013 | +9.14e-4 | 0.90 (transitional↓) | ≥0.730, ~0.733–0.744 |
| 2 |     | 0.473408 | +7.4e-5  | ~1.4 (oscillating)   | ~0.47–0.55 (loose) |
| 3 | ✓   | 1.237070 | +4.1e-5  | 0.73 | **1.2372** (tight) |
| 4 |     | 0.860520 | +1.2e-4  | 0.66 | **0.8608** (tight) |
| 5 |     | 0.764975 | −1.7e-5  | 0.47 | **0.76496** (tight) |
| 6 | ✓   | 1.371654 | +3.7e-5  | 0.51 | **1.3717** (tight) |

**k≥3 converge fast (rates <0.73) ⟹ tightly determined** — new constants nobody had computed. k=1 (75% of the weight)
carries the transitional rate; k=2 oscillates. So the uncertainty is concentrated exactly in the two slow channels.

## STEP 4 — HEADLINE: the tail-free determination SEPARATES from 7/15

`S_∞ = 2·Σ_{k≥1}4^{-k}γ_∞(k)`. Two robust readings and one elegant equivalence:

1. **Finite-R anchor, no extrapolation, no rate assumption:** `2·Σ4^{-k}γ_16(k) = 0.47182` — already **above 7/15=0.46667**
   by 0.0052 (and ≈ S_16, as the relation demands).
2. **Clean lower bound excluding 7/15.** Split `S_∞ = ½·γ_∞(1) + [2·Σ_{k≥2}4^{-k}γ_∞(k)]`, tight part ≈ **0.107**.
   Since `γ_r(1)` is increasing (`A_r(1)>0` = d1>0, established to r=16), `γ_∞(1) ≥ γ_16(1)=0.730`, so
   **`S_∞ ≥ ½(0.730)+0.107 = 0.472 > 7/15`.** 7/15 would require `γ_∞(1)=0.719 < γ_16(1)` — i.e. the k=1 channel
   reversing, which contradicts d1>0. **7/15 is excluded unless the m=1 channel turns over.**
3. **The elegant equivalence (Wilson's two numbers, now forced to agree):** the relation makes
   **`γ_∞(1)=11/15 ⟺ S_∞≈0.4737`.** 11/15 and 0.473 are no longer independent coincidences — they are the *same*
   hypothesis under the exact relation, and it needs the k=1 rate ≈0.80 (the steepened deparity rate), not 0.90.

Full extrapolated band (biased high by the transitional k=1 rate 0.90 and the loose k=2): `S_∞ = 0.483 [0.474, 0.493]`.
Corrected for the physical k=1 rate (~0.80 ⟹ γ_∞(1)≈11/15): `S_∞ ≈ 0.474`. Every reading lands on the **0.473–0.477
side, above 7/15** — a **third independent method** (after S2 splitting and TSW/deparity), and the only one that never
touches the ε tail.

## Net
- **Normalization verified** (`S_1=2/3` exact ⟹ relation is `S_∞/2` exactly; bookkeeping exact). Wilson's concern closed.
- **α closed the sign-phase route** (nondescript, done looking) — the constants are where the signal is.
- **New tight constants** `γ_∞(3)=1.2372, γ_∞(4)=0.8608, γ_∞(5)=0.76496, γ_∞(6)=1.3717` — STEP 5 (PSLQ) targets.
- **7/15 excluded** given the established d1>0 (no k=1 turnover); channel route leans **0.473–0.477**, third independent method.
- The residual uncertainty is the SAME transitional-rate/turnover question, now living in γ_∞(1); but 25% of S_∞ (the
  k≥2 part) is now pinned tail-free, and `11/15 ⟺ 0.473` welds Wilson's candidate to the value.

**Precision caveat:** γ_∞(k) are FLOAT (build_nu rho, ~1e-6) + modest extrapolation → ~4-5 digits; exact rationals only
reach r≤5 here (γ_1=[5/3,2/3,2/3,5/3], γ_2(1)=0.69388…). Real PSLQ (STEP 5) needs higher-precision γ_r(k) than the
vmax-truncated exact builder gives — flag before leaning on any closed form for the tight constants.

**Not at stake:** CHANNEL_ID/CARRYLEMMA identities, R1–R30, R80–R82. Cheap (cached rho + build_nu(11), 3.7s).
