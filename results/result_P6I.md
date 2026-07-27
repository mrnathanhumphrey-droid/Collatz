# RESULT — P6I: the deparitied residual rate is ~0.91, not 0.977 — S_∞ ≈ 0.476–0.478 (raw-rate scare dispelled) (2026-07-26)

**Probe:** `probe_p6i.py`. The one open number after P6H's telescope collapse is the decay rate of the residual
`Λ_i = T_i − T_{i−1}`, `T_i = S_{i+1}/2`. Raw endpoint gave 0.977/level (→ S_∞≈0.50); **HARD RULE: deparity before
quoting**. Vectorized `build_base2` (numpy bincount, same certified SINGLEREC one-step) to push i to 15.

## Wilson's factor-4 reconciliation (structural — it's Re w's numerator)
With `P(a)=2⁻ᵃ`, `ν_e = Σ_{b≥1}4⁻ᵇ(shift₋ᵦ)ρ` (a=2b), so `R̂_e(φ) = Ĉ_ρ(φ)/(17−8cosφ)` and `Re w = (4cosφ−1)/(17−8cosφ)`
— denominators cancel: `4R_e(1) − R_e(0) = ⟨C_ρ, Re w⟩ = Σ_k4⁻ᵏC_ρ(k)`. The three-point kernel `(4,−1)` **is Re w's
numerator**; ν_e absorbs its denominator (the even-branch smoothing). No discrepancy — the hand-conversion was the slip.
(Consequence, dropped as a target: `4·slope(R_e(2))=slope(R_e(0))` asymptotically is just `Λ_i→0`, forced by S converging.)

## GATE — vectorized == certified, bit-for-bit
`build_base2_fast == build_base2` to **1e-17** (i=2..7). Refactor, not reconstruction. Pushed to i=15 (i=15 build 99s).

## The residual settles — the raw 0.977 was a transient wobble
`Λ_i` (×10⁻³): i=8..15 = 0.369, 0.364, 0.391, 0.386, **0.337, 0.320, 0.287, 0.262**. The i=8–11 values wobble
(2-cycle); from **i=12 the decay is clean and monotone**. Deparitied two-step rate `(Λ_i/Λ_{i−2})^{1/2}`:

| i | 12 | 13 | 14 | 15 |
|---|----|----|----|----|
| deparitied rate | 0.929 | 0.910 | 0.923 | 0.905 |

Both subsequences agree: even-i (8,10,12,14) → 0.923, odd-i (9,11,13,15) → 0.905. **Deparitied rate ρ ≈ 0.91**, not the
raw 0.977 — the raw endpoint straddled the wobble, exactly the raw-rate trap the deparity rule guards against.

## Implied S_∞ ≈ 0.476–0.478 — the 0.50 scare is a raw-rate artifact
Geometric tail from i=15 at the deparitied rate:
| ρ | tail Σ | S_∞ |
|---|--------|-----|
| 0.90 | +0.00236 | **0.47607** |
| 0.93 | +0.00348 | 0.47831 |
| 0.977 (raw) | +0.01113 | 0.49360 |

At the deparitied `ρ≈0.90–0.93`, **S_∞ ≈ 0.476–0.478** — landing on the standing headline (`≈0.4766`, 7/15 superseded),
**not** 0.50 (which required the raw 0.977). `T_i` is still rising (0.2357 at i=15, +0.0023 above 7/30) with no turnover,
so **7/15 = 0.4667 remains conditional on a deep turnover not seen at i≤15**; the deparitied trajectory points above it.

## Verdict — the crux number is the deparitied rate ≈0.91 ⟹ S_∞≈0.476–0.478
The deparity rule earned its keep: raw 0.977 → S_∞≈0.50 was a wobble artifact; the deparitied rate ≈0.91 gives
S_∞≈0.476–0.478, consistent across even/odd subsequences and with the corpus. **Caveat (Wilson's "burned three times"):**
still only 4 deparitied points (i=12..15); the rate drifts slightly down (0.929→0.905), so S_∞ could sit at the lower
end. i=16,17 pending to firm it. Not at stake: P6D–P6H identities (all exact), the S_{i+1}=2T_i collapse, R1–R30.
