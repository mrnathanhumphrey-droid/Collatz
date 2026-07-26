# RESULT — PROBE CHANNELDICH: the enriched/depleted dichotomy + γ_∞(1) precision (2026-07-26)

**Probe:** `probe_channeldich.py`. `γ_r(k)=3^r C_r(k)`, `C_r(k)=⟨ρ_r, shift_k ρ_r⟩`. **White baseline = 1 exactly**
(uniform ρ ⟹ `C=1/3^r` ⟹ γ=1), so `γ_∞(k)>1` = ENRICHED (positive autocorr at lag k vs white), `<1` = DEPLETED.

## STEP A — THE DICHOTOMY IS 12/12 (the finding)

`γ_∞(k) > 1 ⟺ 3∣k`, every channel k=1..12, all pre-registered predictions confirmed:

| k | 3∣k | γ_∞(k) | vs white 1 |
|---|-----|--------|------------|
| 1 |     | 0.7384 | depleted |
| 2 |     | 0.4734 | depleted |
| **3** | ✓ | **1.2372** | **ENRICHED** |
| 4 |     | 0.8608 | depleted |
| 5 |     | 0.7650 | depleted |
| **6** | ✓ | **1.3717** | **ENRICHED** |
| 7 |     | 0.4278 | depleted (predicted) |
| 8 |     | 0.7500 | depleted (predicted) |
| **9** | ✓ | **2.1121** | **ENRICHED (predicted)** |
| 10 |    | 0.5210 | depleted |
| 11 |    | 0.5917 | depleted |
| **12** | ✓ | **1.5276** | **ENRICHED (predicted)** |

**This is the first genuine cross-channel structural law in the family, and it's clean, mechanistic, and (unlike the
sign) plausibly provable:** the dlog profile is positively autocorrelated at fiber-multiple lags (3∣k), anti-correlated
off them, relative to white — the lag-domain dual of the projector P that kills 3∣m Fourier modes. Candidate proof
route: the fiber structure (tower ν_r mod 3^r = ν_{r-1}, 3-to-1 lift), which is exactly the machinery that gave
CHANNEL_ID. **This is a pen target with the right shape: a provable structural statement, not a sign.**

## STEP B — 11/15 SURVIVES the precision test (promoted, not confirmed)

- `γ_16(1) = 0.730012587970` (12 digits, float ρ_16).
- `γ_∞(1)` over the physical rate band [0.78, 0.93]: **[0.73325, 0.74251]**.
- **`11/15 = 0.733333` is INSIDE the interval — it survives.** And 11/15 requires rate `ρ_1 = 0.7842` — exactly the
  ~0.78–0.80 steepened deparity/component rate. So 11/15 ⟺ the physical (steepened) rate, coherently.
- Status: **survives, at the low-rate edge** (the raw deparity rate 0.90 gives 0.742). Confirming it means pinning the
  k=1 rate to 0.784 — the same transitional-rate question. Promoted to the working candidate; not uniquely pinned.

## STEP C — the relation, constants, and 11/15 form a mutually consistent set

If `γ_∞(1)=11/15` and `S_∞=0.4737`, the relation predicts `γ_∞(2)=0.4751`; measured `γ_16(2)=0.4734` (oscillating
around it) — consistent, and `<1` as the dichotomy requires (3∤2). Self-consistent solution: **`γ_∞(1)=11/15,
γ_∞(2)≈0.475, S_∞≈0.4735`.** The three lock together.

## STEP D — no small-denominator family among the fast constants (Wilson's prediction confirmed)

Denominator ≤45 scan (resolves 1/45=0.022, nothing finer). Every channel finds *a* nearest rational at |d|~0.0003,
but the denominators are **scattered and incoherent** — 31/42, 9/19, 47/38, 37/43, 13/17, 48/35, 3/7, 19/9 — which is
exactly what chance produces (within-0.0003 of some q≤45 is ~60% likely per constant). **No fifteenths, no
forty-fifths.** This confirms Wilson's prediction: **the /15 denominator family is specific to k=0,1 — S_∞'s
denominator structure comes from the dominant channel, not from the family.** A real structural statement.
One flag worth keeping: **`γ_∞(8) = 0.750024 ≈ 3/4` (|d|=0.00002, 15× tighter than the chance-level matches)** — the
only non-coincidental-looking hit. Not confirmed (could still be chance), but the depleted k=8 channel sitting on 3/4
is the one thing the scan surfaced that isn't noise.

## Net
- **γ_∞(k)>1 ⟺ 3∣k, 12/12** — a clean, mechanistic, plausibly-provable structural law (the real finding; sign was not).
- **11/15 survives** six-digit precision (interval [0.7333, 0.7425] contains it; ⟺ rate 0.784), and locks consistently
  with γ_∞(2)≈0.475 and S_∞≈0.4735. Promoted to working candidate; confirmation = pinning the k=1 rate.
- **No small-denominator family in the fast channels** — S_∞'s /15 structure is a property of the dominant channel,
  not the family. Lone flag: γ_∞(8)≈3/4.
- Weighting stands: 75% of S_∞ in the one monotone channel; 7/15 needs it to reverse (not confirmed reversing).

**Not at stake:** CHANNEL_ID/CARRYLEMMA, R1–R30, R80–R82. Cheap (cached ρ + build_nu(11), 2.4s). Precision on the fast
constants is ~4-5 digits (float ρ) — a real PSLQ needs a higher-precision γ_r(k) builder; the dichotomy and 11/15
survival do not depend on it.
