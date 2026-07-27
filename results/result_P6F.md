# RESULT — P6F: the chain CLOSES — C_ρ (hence the whole cascade) is a functional of the single sub-measure ν_e (2026-07-26)

**Probe:** `probe_p6f.py`. Wilson pinned both seams from P6E: the numerator-vs-divided seam is a **multiplier**, and the
base-4↔base-2 lag map is **derived** from the branch-difference parity. This probe gates all of it and closes the chain.

## The two seams, both closed
Since `W = 2⁻ᵃY`, `ν = Σ_a P(a)(×2⁻ᵃ)_*ρ` (P(a)=2⁻ᵃ, a≥1) — so the numerator profile ρ and the divided measure ν differ
by a **pointwise multiplier**, not a transform:
```
C_ν(m) = Σ_d q(d) C_ρ(m+d),   q(d) = Σ_a P(a)P(a+d) = 2⁻|d|/3   (branch-difference law)
ν̂ = D̃·ρ̂,   |D̃|² = q̂ = 1/(5−4cosθ) ∈ [1/9, 1]   (never vanishes ⟹ deconvolution well-defined + stable)
```

## GATES — machine precision, every level
| gate | result |
|------|--------|
| (1) `C_ν = q⋆C_ρ`, j=2..6 | **1e-17..1e-19**; `\|mult\|²==1/(5−4cosθ)` 1e-15; `q̂==1/(5−4cosθ)` 1e-15 |
| (2) deconv round-trip `ρ̂=ν̂/mult` → certified ρ | **2e-17** (n=3,5); `\|mult\|² ∈ [0.1111, 1.0] = [1/9, 1]` |
| (3) parity (derived) | `C_ρ`: even-lag L1=1.0, odd-lag **1e-16**; `C_ν^cross`: even-lag **1e-17**, odd-lag L1=0.444 |
| (4) **CHAIN CLOSURE** | reconstruct certified `C_ρ` from `ν_e`+boundary → **1e-17, every n=2..6, CHAIN CLOSED** |

Gate (3) is Wilson's derivation made numerical: cross-parity ⟺ branch-difference `d` odd, and ρ on `⟨4⟩` = even
base-2 positions ⟹ `C_ρ` on even lags ⟹ `C_ν^same` needs m even, `C_ν^cross` needs m odd. The lag map is
`k=(m±d)/2` over odd d with weight `q(d)=2⁻|d|/3`; dominant `d=±1` gives `k=(m±1)/2` — the "k↔2k±1", now derived with
weights attached, not guessed.

## The chain, every link explicit and certified
```
ν_e ──[collapse P6D]──▶ ν_o = ½·shift(ν_e) + ½β ──[sum]──▶ ν ──[autocorr]──▶ C_ν
    ──[deconv q  = ×(5−4cosθ)]──▶ C_ρ ──[fold→base-4, k-channel]──▶ γ_n(k) ──[Σ_k 4⁻ᵏ, tower]──▶ Λ_i ──[Σ_{i≥2}]──▶ target
```
Every arrow is a fixed linear map, invertible where the path needs it (the deconv, via `|D̃|²∈[1/9,1]`). The only inputs
are `ν_e` and the fixed boundary sub-measure `β = ½(m₁)ν`. **Gate (4) verifies the composition end-to-end:** the
certified `C_ρ` is reproduced from `ν_e`+β to machine precision at every level.

## Verdict — first time the target sits over one sub-measure with no free maps in the path
`Σ_{i≥2}Λ_i = F(ν_e)` for an explicit composed operator F = (Σ_{i≥2})∘(Σ_k4⁻ᵏ, tower)∘(fold)∘(×(5−4cosθ))∘(autocorr)∘
(collapse). Equivalently a quadratic functional `⟨R_e, K⟩ + boundary` whose kernel K is the composition of the flanking
`(I+½Δ²)`, the deconvolution `(5−4cosθ)`, the fold, and the `4⁻ᵏ` channel weighting — Wilson's `⟨R_e, odd-part-of-D̃⟩`
shape with the full kernel now assembled and every map verified. **The closed-form K, and whether `−1/210` falls out of
it, is the pen's** — the path no longer has a coordinate or transport uncertainty in it. `−1/210` does not "drop out"
mechanically; it is now a definite value of one autocorrelation paired against one explicit kernel.
Not at stake: P6D/P6E collapse+identity, P6/P6B/P6C, P1LVL, BRIDGE2, CHANNEL_ID, dichotomy, R1–R30.
