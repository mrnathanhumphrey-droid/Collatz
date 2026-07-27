# RESULT — P6E: collapse→autocorrelation EXACT + D̃ identity; but the base-4↔base-2 lag map is a real parity gap (2026-07-26)

**Probe:** `probe_p6e.py` (+ inline check). Gate Wilson's explicit collapsed cross-parity autocorrelation form, read
off the reindex constants, and connect to Λ. Two clean gates, one honest catch.

## SOLID — (A) the collapse propagates exactly to the base-2 cross autocorrelation
With `ν_o = ½ν_e(·+1) + ½β`, β=(m₁)ν, R_e := autocorr(ν_e), the cross-parity autocorrelation
`X = ν_e⋆ν̃_o + ν_o⋆ν̃_e` (base-2 lag m) satisfies **X(m) = ½[R_e(m+1)+R_e(m−1)] + boundary = R_e(m) + ½Δ²R_e(m) +
boundary** to **1e-17..1e-19, EXACT, j=2..6**. Wilson's two phrasings (flanking-average and second-difference) are the
same and both hold. (`boundary = ½[(ν_e⋆β̃)+(β⋆ν̃_e)]`.)

## SOLID — (C) the D̃ odd-part identity
`Σ_{n odd} 2⁻ⁿ zⁿ = 2z/(4−z²) = D̃(z)−D̃(−z)`, `D̃(u)=Σ_{v≥0}2⁻⁽ᵛ⁺¹⁾uᵛ = 1/(2−u)` — coefficient-exact (Fractions) to
z¹¹. The odd-lag weight IS twice the odd part of the certified branch generating function.

## CONFIRMED (arithmetic) — (B) Wilson's reindex constants
Reindexing `Σ_{k≥1}4⁻ᵏ·½[R_e(2k+1)+R_e(2k−1)]` gives coeff of R_e(n): **c_n = (5/4)2⁻ⁿ (odd n≥3), 1/8 (n=1)**, exactly
Wilson's numbers. Structure: `c_n = 2⁻ⁿ + 2⁻⁽ⁿ⁺²⁾` = the two flanks; 5/4 = 1+¼ (flank spacing); 1/8 = the n=1 defect
(no k=0 lower flank). So the kernel is **odd-part-of-D̃ acted on by the `(I+½Δ²)` flanking kernel** — the 5/4 and 1/8 he
was unsure of are exactly this second-difference smearing, not an error. NOT bare 2⁻ⁿ.

## THE CATCH — (D) the base-4↔base-2 conversion is a LAG-PARITY GAP, and the channel is a different measure
Wilson flagged "the base-4↔base-2 conversion needs care." It's more than care:
- **ν_e lives entirely on EVEN base-2 positions (coset-1), ν_o entirely on ODD (coset-2)** (mass 1/3 even, 2/3 odd —
  the CRT ℤ/2 factor, exact). ⟹ **cross autocorr X is nonzero ONLY on ODD base-2 lags**; same-parity R_ee/R_oo only on
  EVEN lags.
- base-4 lag k → base-2 lag 2k (×4ᵏ=×2²ᵏ) is **EVEN**, where **X ≡ 0**. So **`A(k)=X(2k)` is degenerate** — the
  premise under (B)'s reindex. Wilson's own instinct ("reindex puts weight on ODD n") is the tell: the cross channel is
  an odd-base-2-lag object, so the map is not k→2k; pinning it (2k±1? how the boundary distributes) is UNRESOLVED.
- **Second seam:** the channel `A_j(k)` is the autocorrelation of the **3x+1 NUMERATOR profile ρ** (bridge domain,
  `γ_r(k)−1=Σ|ρ̂|²e(ak/N)`), whereas ν_e/ν_o are the **DIVIDED forward measure**. Different measures, joined only by the
  bridge (with its |·|² and fold). So the transport from the base-2 collapse to Λ is the bridge, **not a lag relabel**.

## Verdict — the collapse→R_e reduction is exact on the base-2 (divided) side; the transport to Λ is NOT yet a relabel
The base-2 collapse (P6D) and its autocorrelation form (A) are exact and clean, and Wilson's constants (B) and D̃
identity (C) are confirmed. **But `d₁ / Λ = clean R_e functional` is NOT established:** connecting `Σ4⁻ᵏA(k)` to R_e has
two unpinned seams — the base-4-channel↔base-2 lag-parity map (cross = odd base-2 lag, so k↔2k±1 not 2k), and the
numerator-profile-ρ vs divided-ν seam (the bridge). Per the coordinate-map guardrail I hand both back to the pen rather
than guess. **−1/210 still does not drop out** (agreed — it becomes ⟨R_e, odd-part-of-D̃⟩ once the map is pinned).
Not at stake: P6D collapse, P6/P6B/P6C, P1LVL, BRIDGE2, CHANNEL_ID, dichotomy, R1–R30.
