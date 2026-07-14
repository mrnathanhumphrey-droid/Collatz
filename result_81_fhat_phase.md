# R81 disposition — F-hat phase profile on the R78 support

**Date:** 2026-07-13. **Verdict: H_GROWING_DEGREE (refutes H_QUAD; not fixed-degree, not pseudo-random).**

Probe `result_81_fhat_phase_profile.py`; data `result_81_fhat_phase_data.csv`; log `result_81_log.txt`.

Smoke (Th 78.3, |F̂|=3√q constant; support a clean coset mod 3, size 3^(r-1)): **PASS**. All 6 c_{ℓ,ε}=2^ε·(1+3^r)^ℓ tested (ε=0 → support a≡1 mod 3, ε=1 → a≡2 mod 3).

## Exact method

- **Phase certified exactly (no √q):** F̂(3a)²/(9q) is a root of unity ζ_q^s in Z[ζ_q], certified by integer cyclotomic reduction mod Φ_q=x^(2·3^r)+x^(3^r)+1. Sign σ∈{±1}: σ=−1 for r even is the √−3 quadratic-Gauss twist. Exact phase index **J₄∈Z/4q**, arg F̂ = 2π·J₄/4q.
- **Congruence fit (branch-free):** support is consecutive in b (a=a0+3b), so polynomiality in a ⇔ polynomiality in b of equal degree, tested EXACTLY by finite differences (Δ^(g+1)≡0 mod D ⇔ degree ≤ g). This sidesteps the singular-mod-3 Vandermonde entirely. r=2 excluded (|A₂|=3). No magnitude filter.

## Phase degree in b (finite-difference, mod 4q) — the result

| r | degree(s) across 6 c-combos | Δ^deg (const) | v₃(Δ^deg) |
|---|---|---|---|
| 2 *(excluded)* | — | | |
| 3 | [3] | 108 | 3 |
| 4 | [4] | 324 | 4 |
| 5 | [4] | 324 | 4 |
| 6 | [5] | 5832 | 6 |

Degree pattern (r=3,4,5,6): **3, 4, 4, 5** — grows ≈ ⌊r/2⌋+2, unbounded. Leading finite differences are 3-adically deep (v₃ = 3,4,4,6), the fingerprint of a **3-adic-analytic** phase, not a fixed-degree polynomial.

## Decision (§3′ rule)

- **H_QUAD (degree 2 in a at every r≥3): REFUTED.** Degree is ≥3 at every r≥3 (mod 4q AND mod q agree). No denominator D rescues degree 2.
- **H_LIN: refuted** (degree ≥3).
- **H_POLY_HIGHER (fixed degree 3–4): does NOT hold** — the degree is not fixed; it grows with r (3→4→4→5). A uniform bounded-degree Weyl/van-der-Corput route therefore does not exist.
- **H_PSEUDO (equidistributed random): refuted** — the phase is exactly polynomial at each r (finite differences close), i.e. fully deterministic and structured, not equidistribution-random. The obstruction is *growing degree*, not randomness.

## Mechanism (derived + verified)

Collapsing the u-sum (4 has order d=3^r mod q; 3ad≡0 mod q) gives the **exact identity** (self-tested):

&nbsp;&nbsp;&nbsp;&nbsp;**F̂(3a) = 3·Σ_{j=0}^{d−1} e_q(c·4^j)·e_d(−aj) = 3·ĝ(a)**,

so F̂(3a)/3 is the d-point DFT of the **exponential chirp** g(j)=e_q(c·4^j). A *quadratic* chirp e(αj²) has a flat DFT with quadratic phase (a Gauss sum); the exponential chirp 4^j gives a flat magnitude (Th 78.3) but a 3-adic-analytic phase of r-growing polynomial degree. That is precisely why the Gauss-sum/H_QUAD picture fails.

## Routing (which of the three paper routes this opens/closes)

**Closes the smooth-completion / stationary-phase route as a *uniform* square-root mechanism.** R78 §‘Crucial observation’ needed the saving in Σ 1̂(3a)·F̂(3a) to come from phase cancellation in the product; H_QUAD was the hope that arg F̂ is a fixed quadratic Gauss sum enabling completing-the-square. It is not: arg F̂ is a polynomial whose degree grows ≈ ⌊r/2⌋+2, so any Weyl/van-der-Corput completion needs Θ(r) differencing steps and yields no uniform √-saving. This is a **certifying negative** for the fixed-degree smooth-completion route (complementary to band-ℓ1 CLOSED and BGK random-like): the residual bilinear bound genuinely needs Burgess-strength input; the Burgess wall is real for this route. Theorems 78.1–78.3 are unaffected (this only concerns the phase's *degree*, not its magnitude). The one door left ajar: the phase being 3-adic-analytic (an explicit exponential chirp) is more structure than 'random' — a p-adic stationary-phase / oscillatory-integral treatment of ĝ(a) is a distinct, non-Weyl avenue, but it is outside R78's current route list and outside this probe's scope.

_Reporting discipline: outcome reported as fired, including the negative. r=2 carries no evidential weight. No within-support magnitude filter applied. A refutation of H_QUAD is stated as a refutation, not a partial._
