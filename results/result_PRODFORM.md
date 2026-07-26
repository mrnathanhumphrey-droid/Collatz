# Probe PRODFORM — pinning Wilson's product form against R16-A — **FALSIFIED: `Ĉ(m) ≠ Π_j w(2⁻ʲm/N)`. The log-correlation of the measured spectrum with the lacunary product is ~0 at EVERY number of factors J=1..2r, BOTH dilation directions (2⁻ʲ and 2⁺ʲ), and single-w at any k (|corr|≤0.012 vs a sanity self-corr of 1.0). So the "exact-resonance vs partial-resonance" interpretation of A(K) is VOID (Wilson's caveat (a) fires; R16-A's non-diagonality was decisive). BUT the core of the direction survives on a firmer footing: A(K)=C(K)=autocorr(ρ)≥0 (in fact ≥0.045·C0, strictly positive) via ρ≥0 alone — NO product form needed — so `⟨cos⟩>0 ⟺ A(1)+A(−1) beats the partial-resonance shifts` is a genuine nonneg-vs-nonneg race (= the AC-LAGS inequality `C(1)>½[C(N/3−1)+C(N/3+1)]`). What's lost is the MECHANISM (why A(1) beats shifts); the product form was to be that mechanism, and it's false. Aggregate check: the raw 3∤m mass S is CONSTANT (~0.47, slightly increasing), not exponentially decaying — the S/C0 decay is C0 growing (normalization artifact), so no Tao-beating aggregate result.**

**Date:** 2026-07-25. Probe `probes/probe_prodform.py`, log `logs/prodform_run.log`. No build_nu (scratchpad rho_12..16). Pins the one line Wilson owed: is the accumulated modulation the lacunary product `Π_j w(2⁻ʲm/N)`, `w(u)=1/(5−4cos2πu)`?

## (1) The product form is FALSIFIED
`corr(log Ĉ, log Π_{j<J} w(2⁻ʲm/N))` over 3∤m, both r=12 and r=16, all J:
| J | 1 | 2 | 4 | r/2 | r−1 | r | r+2 | 2r |
|---|---|---|---|---|---|---|---|---|
| r=12 | −.006 | +.003 | +.004 | +.003 | +.002 | +.002 | +.001 | −.000 |
| r=16 | −.003 | .000 | +.001 | .000 | .000 | .000 | .000 | .000 |

**No climb toward 1 at any J.** Due-diligence variants (r=12): 2⁺ʲ product best |corr|=0.006; single `w(2⁻ᵏm)` best |corr|=0.012 at any k; sanity `corr(log Ĉ, log Ĉ)=1.000`. **`|ρ̂(m)|²` simply does not factor as a product of w's along the ×2 orbit** — exactly what R16-A predicted (the transport `ν̂_r(ξ)=e(·)E_v[ν̂_{r−1}(ξ2⁻ᵛ)]` is a coherent SUM over the whole Geom ladder, non-diagonal, not a per-level single-frequency factor). **Consequence: the resonance interpretation `A(K)=Σ_{Σn_j2⁻ʲ≡K}Π ŵ(n_j)` is VOID** — Wilson's caveat (a) fires. The "exact resonance beats partial resonance" *mechanism* is not how Ĉ is built.

## (2) But A(K) ≥ 0 survives — and needs only ρ ≥ 0
`A(K)` are the Fourier coefficients of `Ĉ(m)=|ρ̂(m)|²`, i.e. `A(K)=invFFT(Ĉ)(K)=C(K)=Σ_s ρ(s)ρ(s+K)` — the **autocorrelation of ρ**. Since ρ ≥ 0 (a density), `A(K)=C(K) ≥ 0` **trivially**, no product form required. Measured `min_K C(K)/C(0)`: **+0.057 (r=12), +0.045 (r=16), zero negatives** — in fact strictly positive, bounded below by ~4.5% of the peak. So Wilson's decomposition is valid:
$$\mathrm{Num}\propto \tfrac{2N}{3}[A(1)+A(-1)] - \tfrac{N}{3}\!\!\sum_{\substack{K\equiv\pm1\,(N/3)\\ K\not\equiv\pm1\,(N)}}\!\! A(K),\qquad A(K)\ge0,$$
and `⟨cos⟩>0 ⟺ A(1)+A(−1) beats the partial-resonance shifts` — a genuine nonneg-vs-nonneg race. This is exactly the AC-LAGS inequality `C(1) > ½[C(N/3−1)+C(N/3+1)]` (A(±1)=C(1), partials = C(N/3±1)), now with the positivity of every term **manifest** and sourced from ρ≥0. **"The positivity we have is the positivity we need" is TRUE** — but it is autocorrelation positivity (ρ≥0), not ŵ>0 through the (false) product.

## What is lost, stated plainly
The product form was to supply the **mechanism**: "exact resonances concentrate mass, partial resonances smear, and the exact wins." That mechanism is **void** — Ĉ is not that product. The inequality `C(1)>½[C(N/3±1)]` is still TRUE (AC-LAGS: holds by 0.19% at r=16) but its REASON is not the resonance structure. Worse for the strategy: Wilson's aim was to make the inequality **about the transport kernel, not ν_r** (kernel statements are provable). With the product form false, the inequality **remains a property of ρ** (the Syracuse dlog profile) — a comparison of autocorrelation values of the density itself. The positivity reframing is real and clean; it does not, by itself, escape the "property of ν_r" character.

## (3) The aggregate contribution does NOT pan out (normalization artifact)
| r | S = Σ_{3∤m}Ĉ | C0 = Σ_m Ĉ | S/C0 | S/C0 rate |
|---|---|---|---|---|
| 12 | 0.46894 | 6.803 | 0.06893 | — |
| 14 | 0.47025 | 7.743 | 0.06073 | 0.9405 |
| 16 | 0.47135 | 8.685 | 0.05427 | 0.9468 |

**The raw fluctuation mass S is essentially CONSTANT (~0.469→0.471, slightly INCREASING), not decaying.** C0 grows (~1.06/level). So `S/C0` decays ~0.94 **only because C0 grows** — a normalization artifact, not an exponential decay of the aggregate. **No result beating Tao's pointwise superpolynomial in the aggregate** — the constancy of S is just the whiteness/Parseval bookkeeping (R14 forced-uniformity), not a new bound. Wilson's hoped contribution doesn't materialize.

## Tao's theorem — the obstruction in final form (Wilson's, upheld)
Tao Prop 1.14/1.17: `|E e^{−2πiξ·Syrac/3ⁿ}| ≪_A n^{−A}` uniform in n AND ξ (3∤ξ) — our exact Ĉ shell. **Useless here for a reason sharper than "magnitude-only": ⟨cos⟩ is a normalized RATIO, and a ξ-uniform bound appears identically in numerator and denominator and cancels exactly.** Scale-uniformity is invisible to a scale-invariant functional. Banked as the final form of the obstruction — the reason not to go shopping the decay shelf again. (Note also: the +1/4 sits at n=±2, not n=0 — the n=0/white mode contributes exactly 0; confirmed, Wilson's audit-correction stands.)

## Status
**PRODFORM:** ⚠️**Product form FALSIFIED** — `corr(log Ĉ, log Π_j w(2⁻ʲm))≈0` at all J=1..2r, both dilation directions, single-w any k (|corr|≤0.012 vs self-corr 1.0). `|ρ̂|²` does not factor as a ×2-orbit product of w's (R16-A non-diagonality decisive). **Resonance interpretation of A(K) VOID** (Wilson caveat (a)). ⭐**A(K)=C(K)=autocorr(ρ)≥0 SURVIVES via ρ≥0 alone** (min ≥0.045·C0, no product needed) ⟹ `⟨cos⟩>0 ⟺ A(1)+A(−1) beats partial shifts` = AC-LAGS `C(1)>½[C(N/3±1)]` as a nonneg-vs-nonneg race, positivity manifest. **"Positivity we have = positivity we need" TRUE (from ρ≥0), but the MECHANISM (why A(1) beats shifts) is void** and the inequality remains a property of ρ, not the kernel. ⚠️**Aggregate mass S CONSTANT (~0.47), not decaying — S/C0 decay is C0 growth = normalization artifact; no Tao-beat.** ⭐**Tao Prop 1.14/1.17 = our exact object; cancels in the normalized ratio ⟨cos⟩ (scale-uniformity invisible to scale-invariant functional) = obstruction's final form.** Not at stake: R1–R30, R80–R82, MODES/RATIO-2/AC-LAGS/SPECTILT/POINCARE/MODREG/LITHUNT. Owed to pen: the inequality `C(1)>½[C(N/3±1)]` (autocorr of ρ, all terms ≥0, margin 0.19%) still needs a mechanism — the product/resonance route is closed. rho_12..16 in scratchpad. commit pending.
