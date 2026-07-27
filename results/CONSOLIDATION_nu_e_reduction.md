# Consolidation — the ν_e reduction of S∞ (P6D–P6K, SOLSTICE), for the pen

All numerics are machine-verified at the stated precision; "gate" = agreement with a certified prior object.

## 0. Objects and conventions

Level `n`. Units of `ℤ/3^{n+1}`. `2` is a primitive root mod `3^{n+1}`, so
`dlog₂ : (ℤ/3^{n+1})* → ℤ/(2·3ⁿ)` is a bijection with `×2 = +1`, and by CRT
`ℤ/(2·3ⁿ) ≅ ℤ/2 × ℤ/3ⁿ`, where the `ℤ/2` factor is the `⟨4⟩`-coset = branch parity = `x mod 3`, and the `ℤ/3ⁿ`
factor is the base-4 dlog.

- `ν` = stationary forward Syracuse measure of `x ↦ (3x+1)/2ᵛ`, `v = v₂(3x+1)`, `P(v) = 2⁻ᵛ` (v ≥ 1).
- `push_a` = pushforward by one step conditioned on `v = a`. `ν = Σ_{a≥1} P(a)·(×2⁻ᵃ)∘push_a·ν` (fixed point).
- `ν_e, ν_o` = branch-even / branch-odd sub-measures (`a` even / odd). Masses `Σ_{a even}2⁻ᵃ = 1/3`, `Σ_{a odd} = 2/3`.
- `β := (m₁)_*ν` = the `a=1` pushforward (the boundary sub-measure).
- `ρ` = numerator profile = pushforward of `ν` by `x ↦ 3x+1` in the base-4 dlog (coset-1 = `⟨4⟩`). The channel object.
- Channel `γ_r(k) = 3ʳ⟨ρ_r, shift_k ρ_r⟩` (base-4 lag `k`). Dichotomy `γ_∞(k) → M₊ = 5/3` for `3|k`, `M₋ = 2/3` for `3∤k`.
- `R_e(m) = Σ_t ν_e(t)ν_e(t−m)` = autocorrelation of `ν_e` (base-2 lags). Base-4 lag `k` = base-2 lag `2k`.
- S-ladder: `S₁ = 2/3`, `S₂ = 10/21`; `S^{(j)} = 2/3 + 2Σ_{i≤j}Λ_i = S_{j+1}`; `Λ_i = Σ_{k≥1}4⁻ᵏA_i(k)`,
  `A_i(k) = γ_i(k) − γ_{i−1}(k)` (primitive shell, telescoping). Target: `S∞ = 7/15 ⟺ Σ_{i≥2}Λ_i = −1/210`, `Λ₁ = −2/21`.

## 1. Collapse (P6D)

**`ν_o = ½(×2⁻¹)_*ν_e + ½·(m₁)_*ν`**, exact. In base-2 dlog: `ν_o[t] = ½ν_e[t+1] + ½β[t]`. Forced by
`(×2⁻¹)∘push_a = push_{a+1}` (one more division by 2), together with `ν_e = Σ_{a even≥2}2⁻ᵃpush_a`,
`ν_o = Σ_{a odd≥1}2⁻ᵃpush_a`. Gates (n=2..6): identity residual `1e-16`; fixed-point `|reduce(ν_e+ν_o) − ν| = 1e-17`;
masses `(1/3, 2/3)` exact. ⟹ `ν_o` carries no information beyond `ν_e` and the fixed boundary `β`.

## 2. Autocorrelation propagation (P6E)

Cross-parity autocorrelation `X(m) := (ν_e⋆ν̃_o + ν_o⋆ν̃_e)(m)` satisfies
**`X(m) = ½[R_e(m+1)+R_e(m−1)] + boundary = R_e(m) + ½Δ²R_e(m) + boundary`** (gate `1e-17`, all m, n=2..6),
`boundary(m) = ½[(ν_e⋆β̃)+(β⋆ν̃_e)](m)`. Branch generating function identity:
`Σ_{n odd}2⁻ⁿzⁿ = D̃(z) − D̃(−z)`, `D̃(u) = Σ_{v≥0}2⁻⁽ᵛ⁺¹⁾uᵛ = 1/(2−u)` (coefficient-exact).
`ν_e` lives on even base-2 positions, `ν_o` on odd ⟹ `X` supported on **odd** base-2 lags, same-parity autocorr on **even**.

## 3. The two seams close (P6F)

(i) **Numerator↔divided is a multiplier, not a transform.** `ν = Σ_a P(a)(×2⁻ᵃ)_*ρ` ⟹ **`ν̂ = D̃·ρ̂`** pointwise;
equivalently `C_ν = q ⋆ C_ρ` with `q(d) = Σ_a P(a)P(a+d) = 2⁻|d|/3` and **`q̂ = |D̃|² = 1/(5−4cosθ) ∈ [1/9,1]`** —
never vanishes, so the deconvolution `Ĉ_ρ = (5−4cosθ)·Ĉ_ν` is well-defined and stable.
(ii) **Lag map derived:** cross-parity ⟺ branch-difference `d` odd; the base-4 channel `k` draws `C_ρ` at lags
`k = (m±d)/2` with weight `q(d)`, dominant `d=±1 ⟹ k=(m±1)/2`.
Gates (n=2..6): `C_ν = q⋆C_ρ` `1e-17`; deconv round-trip `ρ̂ = ν̂/D̃ → ρ` `2e-17`; **chain closure** — certified `C_ρ`
reconstructed from `ν_e + β` alone (`ν_o` via §1) `1e-17`.

## 4. The three-point kernel (P6G, with factor-4 reconcile from P6I)

In Fourier, `c := cosθ`: flanking `= c`; `Ĉ_ν = (5/4 + c)R̂_e + B̂`; deconvolution `(5−4c)`; channel weight
`Re ŵ = (4c−1)/(17−8c)`. Composition: `R̂_e = Ĉ_ρ/(17−8c)`, so the `(17−8c)` denominators cancel and the pairing is
**`⟨C_ρ, Re w⟩ = 4R_e(2) − R_e(0)`** (base-2 lags 2, 0). Kernel `K = {±2: 2, 0: −1}`, `K̂ = 4cos2θ − 1` = the
**numerator** of `Re w`; `ν_e` absorbs the denominator. Telescoping:
**`Λ_i = 3ⁱ[4R_e⁽ⁱ⁾(2) − R_e⁽ⁱ⁾(0)] − 3ⁱ⁻¹[4R_e⁽ⁱ⁻¹⁾(2) − R_e⁽ⁱ⁻¹⁾(0)]`**.
Gate vs certified `Λ_i` (shellA) `1e-16` (i=2..6); `Λ₁ = −2/21` exact; `3ⁱ` normalization validated
(`3ⁱ⟨ρ,shift₁⟩ → 0.733 = M₋ + cascade`); **boundary term `= 1e-18` (identically zero)**.
**Positivity:** `⟨C_ρ,Re w⟩ > 0 ⟺ R_e(2) > ¼R_e(0)` — Chebyshev / covariance family (the `m=0` proof's family;
first time the `m=1` side is in-family). Holds every level i=1..6.

## 5. The telescope (P6H)

`T_i := 3ⁱ[4R_e⁽ⁱ⁾(2) − R_e⁽ⁱ⁾(0)] = Σ_{k≥1}4⁻ᵏγ_i(k) = ½S^{(i)} = ½S_{i+1}`. Hence
**`S_{i+1} = 2·3ⁱ[4R_e⁽ⁱ⁾(2) − R_e⁽ⁱ⁾(0)]`**, `S∞ = 2·lim_i T_i`, and **`7/15 ⟺ lim T_i = 7/30`**.
Anchors: `T₀ = 1/3 = ½S₁`, `T₁ = 5/21 = ½S₂`. Gate: `S_{i+1} = 2T_i` vs certified S-ladder `1e-15` (i≤6).
`3ⁱR_e(0)` and `4·3ⁱR_e(2)` each diverge ~linearly (slopes `0.03121`, `0.03158`, like `X_i = 3ⁱΣν²`); **`T_i` is the
bounded residue of their cancellation** (Wilson). The asymptotic slope equality `4·slope(R_e(2)) = slope(R_e(0))` is a
consequence of `Λ_i → 0`, not an independent condition.

## 6. Rate and the value/sign separation (P6I, SOLSTICE)

Matrix-free vectorization of the level build, gated bit-for-bit `1e-17`, reaches i=16 locally / i=18 on GPU.
`Λ_i` (×10⁻³): `i=12..18 = 0.33677, 0.31971, 0.28672, 0.26193, 0.23426, 0.20665, 0.17957` — all positive, decreasing.
**Value/sign separation:** `S∞ = S₁₆ + 2Σ_{i≥16}Λ_i` with **`S₁₆ = 2T₁₅ = 0.471352` exact**. Any `Λ_i ≥ 0` past 16 forces
`S∞ ≥ 0.4714 > 7/15`; the rate moves the value only *within* the bracket. **Only `Λ < 0` (a sign change) reaches 7/15**,
requiring `Σ_{i≥19}Λ_i = −0.00296` (≈ a sustained negative run) against `Λ₁₈ = +1.8e-4`.

## 7. Hemispheres (P6J)

`T_i = T⁺_i + T⁻_i`, `T± = Σ_{3|k / 3∤k}4⁻ᵏγ_i(k)`. Gates: `(1/3)M₊ + (2/3)M₋ = 1` exact all i (average of `γ_i(k)`
over all k is exactly 1); `γ_∞(3,6,9) → 1.2370, 1.3716, 2.1115` (certified 1.2372, 1.3717, 2.112).
**`T⁺` (enriched, `3|k`) converges** — `Λ⁺ → 0` by i≈13, `lim T⁺ = 0.019672`. **`T⁻` (depleted, `3∤k`) carries the entire
slow residual**, rate ≈ the whole-system rate. `S∞/2 = 0.01967 (enriched, done) + 0.2186 (depleted)`. The open sign
question is confined to the depleted hemisphere.

## 8. λ-family (P6K) — MODEL FAMILY, not Collatz at λ≠½

Deform `P(v) ∝ λᵛ`; drift `E[3·2⁻ᵛ] = 3(1−λ)/(2−λ) = 1` at λ=½ (critical). Gate M-A: `T_i(½)` = certified S-ladder
`1e-9`. `λ>½` converges (rate `3(1−λ)/(2−λ) < 1`), `λ<½` diverges. **The approach direction `sign(lim Λ)` flips at
`λ ≈ 0.53–0.54`; λ=½ sits ~0.03 below on the `Λ>0` (no-turnover) side.** No kink at the critical point
(`S(λ→½⁺) → ~0.47 ≈` critical, continuous within near-critical numerical uncertainty).

## 9. Sign decider (SOLSTICE, i=17,18 on Lambda A100, matrix-free torch, gated on-instance `T₁₅,T₁₆` `2e-9`)

`T₁₇ = 0.23611673`, `T₁₈ = 0.23629630`; `Λ₁₇ = +2.0665e-4`, `Λ₁₈ = +1.7957e-4`.
- **Deparitied two-step rate `(Λ_i/Λ_{i−2})^{1/2}`** — even-i: `0.9285(12), 0.9227(14), 0.9039(16), 0.8755(18)`;
  odd-i: `0.9051(15), 0.8882(17)`. Even-i drift `−0.0118/level`, **`|slope/SE| = 8.53`** (was 3.28 at i≤16),
  curvature `−0.0096` (accelerating).
- **Single-mode `ρ = 0.908` (the i≤16 reading) is falsified out-of-sample:** it predicts `Λ₁₈ = +1.947e-4` vs actual
  `+1.796e-4` (**8.4% miss**).
- **Best two-mode LS fit i=12..18:** `Λ_i = A·0.867ⁱ − B·0.628ⁱ` (A=2.44e-3, B=2.72e-2, resid `1.9e-6`); `ρ₂ < ρ₁`
  ⟹ two *decaying* modes, **no crossing**, dominant `ρ₁ ≈ 0.867`; held-out (fit 12..16 → predict 17,18) within 1.4%.
  The data's late acceleration is marginally steeper than this no-crossing model, so a crossing (`ρ₂ > ρ₁`) is **not
  excluded, only not preferred.**

**Status.** `Λ_i > 0` through i=18 — **no turnover has occurred**. The clean single-mode/no-turnover reading of the
i≤16 probe is falsified (true rate ≈0.867, still falling), which **revises the value down**: exact floor `2T₁₈ =
0.47259`; geometric tail at rate `0.867–0.88` gives **`S∞ ≈ 0.4749–0.4752`** (from ≈0.476 — faster decay, smaller tail).
The sign is unresolved at the knife's edge: the LS optimum has no crossing, but with 7 points and a persistent
accelerating drift a crossing is not ruled out. `i=19,20` separate "asymptote at ≈0.87" (no turnover, `S∞ ≈ 0.475`) from
"continued fall" (turnover toward 7/15).

## 10. Net for the pen

`S∞` is a single limit of one explicit sequence in two autocorrelation values of the mass-`1/3` even-branch
sub-measure: `S_{i+1} = 2·3ⁱ[4R_e⁽ⁱ⁾(2) − R_e⁽ⁱ⁾(0)]`, kernel three points `{±2:2, 0:−1}`, no free maps, no fitted
constant, no boundary. The **value** is bracketed `S∞ ∈ [0.4726, 0.4752]` (floor exact; upper from rate). The **sign** —
whether `Λ` ever goes negative — is where 7/15 lives; through i=18 it has not, the LS fit prefers two decaying modes
(no crossing), and 7/15 requires a sustained negative run of `Σ_{i≥19}Λ = −0.00296`. Open exact statement:
`7/15 ⟺ lim_i 3ⁱ[4R_e⁽ⁱ⁾(2) − R_e⁽ⁱ⁾(0)] = 7/30`; equivalently, whether the depleted-hemisphere autocorrelation
`R_e(2) − ¼R_e(0)` telescoped with weight `3ⁱ` turns its sign.
