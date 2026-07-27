# RESULT — P6D: the cross-term COLLAPSE is EXACT — ν_o = ½(×2⁻¹)ν_e + 2⁻¹(m₁)ν (2026-07-26)

**Probe:** `probe_p6d.py`. Test Wilson's collapse `ν_o = ½(×2⁻¹)_*ν_e` in the **certified base-2 coordinate** he
authorized: 2 is a primitive root mod 3^{n+1}, so `dlog₂: units → ℤ/(2·3^n)` is a bijection with `×2 = +1`, and CRT
`ℤ/(2·3^n) ≅ ℤ/2 (⟨4⟩-coset = branch parity = x mod 3) × ℤ/3^n (base-4 dlog)`. In this coordinate `×2⁻¹` is a plain
shift by −1. P6B's roll-scan failed ONLY because it lived in the folded base-4 coordinate, where `×2⁻¹` flips the
coset (not a roll). Here: no fold.

## The identity is EXACT, not an approximate roll — and Wilson's boundary term is forced
`×2⁻¹·push_a = push_{a+1}` exactly (one more division by 2). With `ν_o = Σ_{a odd≥1}2⁻ᵃpush_a`,
`ν_e = Σ_{a even≥2}2⁻ᵃpush_a` (a=0 excluded, `stationary_trunc` uses a≥1):
```
½(×2⁻¹)ν_e = ½ Σ_{a even≥2} 2⁻ᵃ push_{a+1} = Σ_{b odd≥3} 2⁻ᵇ push_b = ν_o − 2⁻¹ push_1(ν)
⟹  ν_o = ½(×2⁻¹)ν_e + 2⁻¹(m₁)ν      [Wilson's boundary term ½(m₁)ν, EXACT, carried]
```

## GATE + TEST — machine precision, every level n=2..6
Built as a certified SINGLEREC one-step in the residue domain, indexed by base-2 dlog (unfolded):

| n | GATE-1 fixed-point `\|reduce(R)−ν\|` | mass (a-even, a-odd) | COLLAPSE `\|R_o − ½shift(R_e) − B\|` | drop B → residual |
|---|---|---|---|---|
| 2 | 5.6e-17 | 0.33333, 0.66667 | **2.3e-16** | 0.983 (dropped == B: 2e-16) |
| 3 | 2.8e-17 | 0.33333, 0.66667 | **7.7e-17** | 0.990 |
| 4 | 6.9e-18 | 0.33333, 0.66667 | **3.1e-16** | 0.993 |
| 5 | 1.4e-17 | 0.33333, 0.66667 | **1.5e-16** | 0.996 |
| 6 | 6.9e-18 | 0.33333, 0.66667 | **3.1e-16** | 0.998 |

GATE-1: ν is the exact fixed point of the one-step (⟹ it IS the certified Syracuse step, reindexed). Mass split is
**exactly (1/3, 2/3) = the seed** `P(a even)=Σ4⁻ᵃ=1/3`. COLLAPSE holds to **1e-16**; without B the residual is ~1
and the dropped term is **exactly B** — the boundary is load-bearing, and it is precisely `2⁻¹(m₁)ν`.

## Verdict — Wilson's collapse CONFIRMED exactly; P6B "inconclusive" RESOLVED (it was the coordinate)
`ν_o` carries **no information** beyond `ν_e` plus the a=1 boundary. The odd sub-measure is a shifted half-copy of the
even one. **Consequence:** `d₁ = A_j(1) = ⟨even×odd cross-parity⟩` reduces to a functional of the **single** sub-measure
ν_e (its lag autocorrelation) plus a boundary correction — the cascade target `Σ_{i≥2}Λ_i = −1/210` is now expressible
through one sub-measure's autocorrelation, back in the covariance family (the m=0 proof's family). P6B's roll-scan was
not a refutation — it was the folded coordinate; the base-2 coordinate closes it. **Pen next (Wilson):** write d₁ /
A_j(1) explicitly as `½·(autocorr ν_e at the shifted lag) + boundary`, and see whether the `−1/210` falls out.
Not at stake: P6/P6B findings, P6C, P1LVL, BRIDGE2, CHANNEL_ID, MEAN1, dichotomy, R1–R30. Cheap (0s).
