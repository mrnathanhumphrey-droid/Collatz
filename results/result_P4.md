# RESULT — P4: γ_∞(1) is carried by the COARSE 3-adic frequencies, dominated by a=N/3 (= the class mean) (2026-07-26)

**Probe:** `probe_p4.py`. Wilson's "change the pitch of part of them" made trivial: since
`γ_r(k)−1 = Σ_{a≠0}|ρ̂(a)|²e(ak/3^r)` is linear (Wiener–Khinchin on the certified `γ_r(k)=3^r⟨ρ,shift_k ρ⟩`), the
per-frequency influence on channel k IS the summand. Sort the pair-contributions `2|ρ̂(a)|²cos(2πa/3^r)` and read off
which frequencies carry `γ_∞(1)=0.733`. ρ = cached certified dlog profile (r=12,14,16). π̂ (additive) does NOT enter.

## GATE — the summand identity reproduces γ (machine precision, all r)
`1 + Σ_{a≠0}|ρ̂(a)|²e(+2πia/N)` = `3^r⟨ρ,roll(ρ,−1)⟩` to <1e-9 at r=12,14,16 (γ = 0.725655, 0.728069, 0.730013).
Sign convention fixed (numpy `ifft` gives `e^{+2πiak/N}`). So the decomposition is exact.

## THE ANSWER — one coarse frequency dominates; it IS the class mean
The dominant carrier of `γ_r(1)−1` at every r is the **single frequency `a* = 3^{r−1} = N/3`** (a/N=1/3, v₃=r−1, the
MAXIMAL 3-adic valuation — the coarsest frequency):
- **`|ρ̂(N/3)|² = 1/3` EXACTLY, r-invariant.** This is the frozen mod-3 marginal `ν₁=(0,⅓,⅔)`: analytically
  `ρ̂(N/3)=Σρ[j]e^{−2πij/3}=(⅓)ζ+(⅔)ζ²`, `|·|² = 5/9+(4/9)cos(4π/3) = 5/9−2/9 = 1/3`.
- `cos(2π/3) = −½`, so its contribution is **exactly `2·(1/3)·(−½) = −1/3`** — 123% of `γ−1 = −0.270` (the rest is a
  cascade that partially cancels back).
- **`1 − 1/3 = 2/3 = M₋`**, the depleted class mean (3∤k). So the leading Fourier term reproduces the class mean
  exactly, and `γ_∞(1) = M₋ + (high-v₃ cascade) = 2/3 + 0.067 = 0.733`. This ties P4 directly to MEAN1.

## The cascade — coarse (high-v₃) frequencies, r-stable, alternating
Ranks 2–10 (identical structure at r=12,14,16, only a rescales by 3^{Δr}):
| rank | a/N | v₃(a) | cos | \|ρ̂\|² | pair contribution |
|------|-----|-------|-----|--------|-------------------|
| 1 | 1/3 | max | −0.500 | 0.3333 | −0.3333 |
| 2 | 4/9 | max−1 | −0.940 | 0.0494 | −0.0928 |
| 3 | 13/27 | max−2 | −0.993 | 0.0448 | −0.0890 |
| 4 | 1/9 | max−1 | +0.766 | 0.0459 | +0.0704 |
| 5 | 2/27 | max−2 | +0.894 | 0.0344 | +0.0615 |
| 6 | 2/9 | max−1 | +0.174 | 0.1428 | +0.0496 |
| 7 | 4/27 | max−2 | +0.597 | 0.0357 | +0.0426 |

All carriers are **3-adic-rational frequencies `a = c·3^{r−j}` (high v₃, low 3-adic denominator)** — the coarse end.
The cumulative fraction oscillates (1.23 → 1.58 → 1.91 → 1.65 → 1.42 → …) converging to 1: a slowly-convergent
alternating series over the coarse frequencies. **top-1% carriers are 100% v₃≥2**, median a/N ≈ 0.22–0.24.

## NOT diffuse, NOT the π̂-sup — a structured coarse set
- The heavy tail of |π̂| (HOMOG-D) was diffuse; the **carriers of the channel are the opposite — a sharp, structured,
  r-stable set of coarse (high-v₃) 3-adic frequencies**, dominated by a=N/3.
- The **π̂-sup** (fine 2^m, v₃=0) contributes ~nothing to γ — confirmed. (Note: the ρ̂-sup IS a=N/3, the dominant
  carrier — the additive-sup and multiplicative-sup are different objects, and only the latter carries the channel.)

## Net
First clean read of **which frequencies carry a channel**: the coarse 3-adic ones, led by a=N/3 whose weight is the
exact mod-3 marginal (|ρ̂|²=1/3), giving the leading term `1−⅓ = 2/3 = M₋` (the class mean), with a structured high-v₃
cascade converging to `γ_∞(1)=0.733`. This unifies the channel value with MEAN1 (class mean) and the v₃-HIERARCHY
(carriers are high-v₃), and it says the channel is an **amplitude-on-coarse-frequencies** object — which P1 (phase
scramble) will now test directly. Not at stake: CHANNEL_ID, v₃ HIERARCHY, MEAN1, R1–R30. Cheap (4.3s).
