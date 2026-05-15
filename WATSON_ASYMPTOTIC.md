# WATSON_ASYMPTOTIC — Phase 2 results: saddle-point applied to R78/R79 bilinear + Darboux applied to ε_k

**Date:** 2026-05-14. Phase 2 of WATSON probe. Pre-registration locked in WATSON_PRE_REGISTRATION.md.

## Summary of executed analyses

Two parallel threads per WATSON_R78_SUM_FORM.md:

- **Thread α**: p-adic saddle-point on T_p(r) at p=3, r=2..6 (direct computation)
- **Thread β**: Darboux/multi-saddle on f(z) = Σ ε_k z^k using ε_k k=1..13 data

## Thread α — Saddle-point on T_p(r) bilinear sum, p=3, r=2..6

### Direct computation of T_3(r)

Using R78.6's closed form `ψ(a) = F̂(3a) / (3√q)`, we computed:

```
T_3(r) = Σ_{a ≡ 1 mod 3 in Z/3^r} 1̂(3·a) · ψ(a)
```

| r | N = 3^{r−1} | |T_3(r)| | κ = log\|T\| / log N |
|---|---|---|---|
| 2 | 3   | 1.5231  | 0.3830 |
| 3 | 9   | 8.1200  | 0.9532 |
| 4 | 27  | 13.4584 | 0.7888 |
| 5 | 81  | 78.0446 | 0.9915 |
| 6 | 243 | 301.5609 | 1.0393 |

**Log-linear fit r=2..6:** `log|T_3(r)| = -0.816 + 1.169 · log(N)`, so `κ_fit = 1.169`.

### Comparison to R79b empirical r=8..20

R79b reported `|K(r)| ∝ N^{0.522 ± 0.008}` at r=8..20 (R² = 0.9976), where K is essentially T_3.

**The TRANSIENT regime (r=2..6) has κ ≈ 1.17 (super-trivial!).** The asymptotic regime (r=8..20)
has κ ≈ 0.522 (square-root cancellation). The transition occurs at r ≈ 6..8.

### Saddle-exactness check at r=3, 4, 5

The R78.6 saddle prediction ψ_lead(a) = e_q(P_a(s*(C_a))):

- **r=3**: ψ_lead matches ψ_true to machine precision (R78.6 RIGOROUS).
- **r=4**: phase deviation arg(ψ_true / ψ_lead) up to 88°, magnitude deviation up to 1.4.
- **r=5**: phase deviation up to 160°, magnitude deviation up to 2.0.

This confirms PATH2_BILINEAR's diagnosis: **the saddle is exact at r=3 but Hensel-corrections at
r ≥ 4 are not captured by leading-order saddle**. R79b documents this same gap as `|Σ 1̂·ψ_lead| /
|Σ 1̂·ψ_true| ≈ 0.4-0.6` (factor-2 gap).

### Saddle-point yields the wrong asymptotic exponent for r=2..6

If saddle-point gave the CORRECT exponent, T_3(r) should follow R79b's κ ≈ 0.522 even at small r.
Instead, small-r data has κ ≈ 1.17 (essentially linear in N) — i.e., NO cancellation in the
saddle-point representation at small r.

**The saddle-point captures the LEADING phase but the cancellation requires the SUM OVER a to be
nontrivial — and the saddle representation collapses the cancellation by treating each a's
contribution as a phase factor of unit modulus.** The square-root cancellation κ=0.522 must come
from phase cancellation BETWEEN a's, not from any single a's phase magnitude.

**Conclusion (Thread α):** Saddle-point on T_p(r) ALONE does not produce the asymptotic rate
κ=0.522. The Hensel-lifted saddle at r ≥ 4 plus inter-a phase cancellation is required, and
neither classical saddle-point nor Watson lemma supplies the inter-a structure.

This matches what PATH2_BILINEAR_FROM_CLOSED_FORM.md derived: the cleanest rigorous bound is
|T_p(r)| ≤ 2N (i.e., κ ≤ 1) at r=3 using Inner-Plancherel; the empirical κ=0.522 is finer and is
NOT captured by saddle + Cauchy-Schwarz + Plancherel within the family-level machinery.

## Thread β — Darboux/multi-saddle on f(z) = Σ ε_k z^k

### Hadamard radius trajectory

| k | \|ε_k\|^{1/k} | ρ_k = 1/\|ε_k\|^{1/k} |
|---|---|---|
| 10 | 0.4850 | 2.062 |
| 11 | 0.5538 | 1.806 |
| 12 | 0.6022 | 1.661 |
| 13 | 0.6388 | 1.565 |

**Trend: ρ_k is decreasing**. PADE_NUMERICAL_DISPOSITION reports ρ at n=13 = 1.57, consistent.
Linear extrapolation of ρ_k vs 1/k from k=10..13 gives ρ_∞ ≈ -0.11 (negative, meaningless), which
means the trend HASN'T STABILIZED at k=13.

### Ratio test |ε_{k+1}/ε_k|

| k → k+1 | ratio |
|---|---|
| 10 → 11 | 2.084 |
| 11 → 12 | 1.515 |
| 12 → 13 | 1.296 |

Ratios are FALLING from ≫ 1 toward 1 from above. f(z) at k=10..13 is in the TRANSIENT growth
regime. The asymptotic singularity ρ_∞ (slow-mode, predicted by PADE as ρ=1.016) has not yet
been reached at k=13.

### Multi-saddle / Darboux fits

**Constrained fit holding PADE prediction (ρ_1=1.016, ρ_2=1.57, θ_2=2π/9.2=0.683 rad):**

```
ε_k ≈ A · (1.016)^{-k} + R · (1.57)^{-k} · cos(k · 0.683 + φ)
```

Best fit: A = 1.06e-3, R = -3.95e-2, φ = 2.91 rad, RSS = 4.14e-5.

Per-k fit comparison:

| k | ε_k actual | Constrained fit | A · 1.016^{-k} | pair · 1.57^{-k} |
|---|---|---|---|---|
| 1 | +2.0e-1 | +2.4e-2 | +1.0e-3 | +2.3e-2 |
| 2 | +9.5e-3 | +7.8e-3 | +1.0e-3 | +6.8e-3 |
| 3 | −5.1e-3 | −1.5e-3 | +1.0e-3 | −2.5e-3 |
| 4 | −2.5e-3 | −4.2e-3 | +9.9e-4 | −5.2e-3 |
| 5 | −1.2e-3 | −3.2e-3 | +9.8e-4 | −4.1e-3 |
| 6 | −5.0e-4 | −1.0e-3 | +9.6e-4 | −2.0e-3 |
| 7 | −1.2e-3 | +6.7e-4 | +9.5e-4 | −2.8e-4 |
| 8 | −7.5e-4 | +1.5e-3 | +9.3e-4 | +5.3e-4 |
| 9 | −7.5e-6 | +1.6e-3 | +9.2e-4 | +6.4e-4 |
| 10 | +7.2e-4 | +1.3e-3 | +9.0e-4 | +4.1e-4 |
| 11 | +1.5e-3 | +1.0e-3 | +8.9e-4 | +1.5e-4 |
| 12 | +2.3e-3 | +8.6e-4 | +8.8e-4 | −1.9e-5 |
| 13 | +2.9e-3 | +7.8e-4 | +8.6e-4 | −8.0e-5 |

**The PADE-prediction model fits k=1..6 ROUGHLY but fails k=7..13 by a factor of 2-100.**
Specifically: the fit predicts k=13 at +7.8e-4 but actual is +2.9e-3 — off by 4×.

### Unconstrained 5-parameter fit (ρ_2 held at 1.57)

Best fit finds ρ_1 = 0.54 (NOT 1.016 — pulls a GROWING mode), θ_2 = 1.60 rad period 3.9 (NOT
period 9.2). RSS = 1.64e-5.

This is the multi-saddle's BEST fit to k=2..13 data when we hold ρ_2=1.57. The asymptotic real
singularity is at ρ_1 = 0.54 (modulus < 1!) which means **f(z) has a singularity INSIDE the unit
disk** at k=2..13. This contradicts PADE's claim of ρ_∞ ≈ 1.016 (just OUTSIDE the unit disk).

### Pure complex-pair fit (k=4..13)

Best fit: ρ = 0.937, θ = 7e-4 rad (= 0 essentially), R = 0.43, RSS = 1.35e-6.

This is a SLOWLY-VARYING decaying mode with ρ ≈ 1 (almost no oscillation). It fits the
k=4..13 data well (RSS = 1.35e-6 ≪ 4.1e-5 for PADE-constrained).

**This is the strongest fit but it doesn't match the PADE prediction.**

### Why PADE prediction fails to match the fit

Multiple competing fits all have RSS ~ 1e-5..1e-6 — meaning the data at k=2..13 (12 points,
6 free parameters in two-mode model) is **underdetermined** for multi-mode fitting. The
near-zero ε_9 (which forces sign change near n=9-10) is a STRONG constraint, but with only 12
data points and the sign change happening at k=9, ANY model with a single zero crossing fits.

The PADE prediction (period 9.2) requires a SHARP oscillation, which the data doesn't have
yet (only 1 sign change). At k=20+, if data continues, multiple oscillations would constrain
the period — but at k=13 we can't distinguish period 9.2 from period 30+.

## Phase 2 numerical verification at k=8..13

Using the BEST two-component fit (free 5-param: ρ_1=0.54, ρ_2=1.57, θ_2=1.60):

| k | ε_k actual | model | rel err |
|---|---|---|---|
| 8  | −7.46e-4 | varies (see below) | — |
| 9  | −7.52e-6 | varies | — |
| 10 | +7.21e-4 | +5.75e-4 | −20% |
| 11 | +1.50e-3 | +1.03e-3 | −32% |
| 12 | +2.27e-3 | +1.83e-3 | −20% |
| 13 | +2.95e-3 | +3.25e-3 | +10% |

Even the BEST 5-param fit has 20-30% errors at k=10..12.

**No multi-saddle Darboux model fits the available ε_k data within tight tolerance.**

## Synthesis

### What saddle-point + Darboux DO match

- **Hadamard radius trend** ρ_k → 1.57 at k=13 is CONSISTENT with the PADE leading-singularity
  picture (ρ_2 ≈ 1.57 at transient).
- **Sign pattern + + − − − − − − − + + + +** with zero-crossing near k=9 is CONSISTENT with
  cos(k θ + φ) modulation (one half-cycle in k=4..13).
- **Saddle-point at r=3 is EXACT** (R78.6) — confirms the saddle technique applies at the
  level-r side.

### What saddle-point + Darboux DO NOT match

- **Asymptotic exponent κ for T_p(r):** small-r computation gives κ≈1.17 (super-trivial), large-r
  empirical R79b gives κ≈0.522. Transition at r≈6..8. Saddle alone doesn't predict κ; it
  needs the additional inter-a phase cancellation that R78/R79 identify as OPEN.
- **PADE period 9.2** does NOT appear as the best-fit period to ε_k data (best fit prefers
  period 3.9 in 5-param free fit, or period ∞ in pair-only fit).
- **Slow-mode ρ_1 = 1.016**: not visible in ε_k data k=2..13; the data is in TRANSIENT regime
  before slow-mode asymptotic kicks in. Free 5-param fit prefers ρ_1 = 0.54 (a growing mode).
- **Faure prediction ρ = √3 ≈ 1.732**: held-ρ fit with √3 gives RSS = 1.75e-5, comparable to
  but worse than held-1.57 (RSS = 1.71e-5). Both are POOR fits compared to free fits.

### The closure-conversion question

Does the Darboux asymptotic convert to a polynomial-in-A bound on |μ̂_n(ξ)| via R75/R76/R77?

The conversion mechanism: f(z) = Σ ε_k z^k. The singularity at |z| = ρ_∞ corresponds to a
spectral radius of the transfer operator T_M of 1/ρ_∞. If ρ_∞ > 1 (slow-mode at 1.016), then
T_M has spectral radius < 1, hence ε_k → 0 exponentially and |μ̂_n(ξ)| → 0.

**For this conversion to be RIGOROUS**, we need:
- An analytic argument that f(z) IS analytic in |z| < ρ_∞.
- A bound on f(z) growth as |z| → ρ_∞.
- A characterization of the singularity type (algebraic order, branch cuts).

**None of these are provided by saddle-point + Darboux.** They are properties of f(z) that
require transfer-operator analytic theory — exactly the profinite semiclassical infrastructure
that FAURE_DISPOSITION identifies as missing.

## Phase 2 conclusion (preview for disposition)

**Saddle-point applied to T_p(r) bilinear sum (Thread α):**
- Recovers the R78.6 saddle-exact closed form at r=3 (already RIGOROUS)
- Does NOT extend to r ≥ 4 without Hensel-lifted correction (OPEN per R79b/R78.6)
- Does NOT produce the empirical rate κ = 0.522 at r=8..20 (which requires inter-a cancellation
  beyond saddle's reach)

**Darboux multi-saddle applied to f(z) (Thread β):**
- Hadamard tail trajectory and sign pattern are CONSISTENT with PADE picture (and with Faure
  spectral-gap framework qualitatively)
- Multi-saddle fits to ε_k k=2..13 data are UNDERDETERMINED — period 9.2 is not strongly
  preferred by data alone
- Predicted asymptotic ε_k → A · ρ_∞^{-k} requires ρ_∞ KNOWN (rigorously) — saddle doesn't
  provide this
- Conversion to |μ̂_n(ξ)| bound requires transfer-operator analyticity that's NOT supplied by
  saddle-point/Darboux machinery

**Neither thread alone closes c=7/45.** The closure target (polynomial-in-A bound on |μ̂_n|)
requires:
1. Analytic structure of f(z) (transfer operator theory — Faure-style smooth analog needed)
2. Rigorous identification of ρ_∞ (which Padé estimates from finite data; rigorous proof open)
3. Inter-a phase cancellation in T_p(r) at r ≥ 4 (saddle alone misses this; Plancherel/Cauchy
   give only κ ≤ 1, the trivial bound)

Saddle-point + Darboux supply the EXPECTED FUNCTIONAL FORM (complex-pair singularity → cos(kθ+φ)
modulation, real singularity → exponential decay) but not the RIGOROUS INPUT (location, type,
existence) that the form requires.

This is the **PARTIAL** disposition, sharpened to a specific gap.
