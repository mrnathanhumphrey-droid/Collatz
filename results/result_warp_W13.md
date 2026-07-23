# Warp W1/W3 — the subordination "costume check" — **W3 confirmed (affine shift = genuine deformation, four lines safe); W1 value refuted (a moment-confusion); ½ is second-moment**

**Date:** 2026-07-22. Probe `probes/probe_warp_W13.py` (fixed-level first-moment model on ν̂). Wilson's warp
(pen): our random-scale renewal is the fixed-λ Bernoulli operator run for a random *geometric* time — a
resolvent `T_ours = (½B)(I−½B)⁻¹` of the single-step operator B. The certified R28 operator splits as
`ν̂_r(ξ) = e(ξ/3^{r+1})·E_v[ν̂_{r−1}(ξ2^{−v})]` = **affine dressing `A_*` = diag(e(ξ/3^{r+1}))** ∘ **linear core**
`E_v[·] = Σ_v 2^{−v}(halving)^v`. Pre-registrations: **W1** core subdominant → 1/3; **W3** is `A_*` a
spectrum-moving deformation (½ real) or a spectrum-preserving conjugation (½ an artifact, would reverse four lines).

## W1 — the unshifted core has **no spectral gap** (not 1/3)
`CORE = Σ_{v≥1} 2^{−v}(ξ ↦ ξ2^{−v} mod 3^r)`, the exact resolvent of the halving permutation:

| r | core \|λ₂\| |
|---|---|
| 3 | 0.898 |
| 4 | 0.987 |
| 5 | 0.9985 |
| 6 | 0.99983 |

**|λ₂| → 1: no gap.** The eigenvalues on the units are exactly `1/(2χ(2)−1)` over characters χ. The
smallest-order-3 character gives `1/(2e^{2πi/3}−1) = 1/√7 ≈ 0.378`, but the **generator** character
(χ(2)=e^{2πi/φ}, closest to 1) gives an eigenvalue → 1. So the halving-permutation resolvent has **no gap and
no 1/3.** The pre-registered `f(½)=1/3` does **not** describe this operator: `cos(π/3)=½` is the *first-moment
Bernoulli-convolution* factor, but the actual unshifted `T_ours` is the *permutation* resolvent, whose spectrum
is on the unit circle, not at ½. **The two "B"s are different operators — 1/3 is a first-moment number applied
to the wrong first-moment operator (the moment-confusion Wilson flagged, made concrete).**

## W3 — the affine phase is a **genuine deformation** (decisive; four lines safe)
`SHIFTED = diag(e^{2πiξ/3^{r+1}})·CORE`:

| r | shifted \|λ₂\| |
|---|---|
| 3 | 0.97008 |
| 4 | 0.97008 |
| 5 | 0.97008 |
| 6 | 0.97008 |

**The phase creates a gap** — from |λ₂| → 1 (core, no gap) to a **stable |λ₂| = 0.970**. So `A_*` is a genuine
**spectrum-moving deformation, not a conjugation.** This is the decisive W3, and it lands on the non-scary side:
the shift is **load-bearing**, the gap is **real and shift-created**, and **no numerical line gets reversed** —
remove the tailor (`A_*`), the gap vanishes. The "+1 is the tailor" reading is confirmed at the operator level.

## But the value here is 0.970, not ½ — because ½ is **second-moment**
Neither 1/3 nor ½ appears. The shift-created first-moment gap sits at **0.970**, the neighborhood of the *old
envelope / period-9 family* (the core's units-only subdominant at r=3 is literally `0.898 @ 38°` → period 9.5 =
the R77 period-9.2 mode), **not** the |λ₂|=½ that governs S_r. Honest resolution of the moment-confusion:
- **|λ₂| = ½ is a second-moment (Plancherel `Σ|ν̂|²`) rate.** It cannot appear in a first-moment operator on ν̂,
  and it doesn't.
- The first-moment operator's slow mode (~0.97, shift-created) is the **envelope/period-9** object — a different,
  slower mode.

## Status
**Warp structure: right. Warp value: misplaced (moment-confusion). Four numerical lines: safe.** **W3 GATE** —
`A_*` is a genuine deformation: the unshifted core has no gap (|λ₂|→1), the shift creates a stable gap (0.970),
so ½ is not a conjugation artifact and the four |λ₂|≈½ lines are not reversed. **W1** — the pre-registered 1/3
is **not** in the first-moment core (which has no gap); `f(½)=1/3` conflates the Bernoulli-convolution factor
(first moment) with the halving-permutation resolvent (the actual unshifted operator). **The ½ lives strictly in
the SECOND moment**, which the first-moment costume check provably cannot carry.

**Consequence for the crux (owed to the pen):** the warp correctly identifies the affine dressing as the
gap-source (deformation, not conjugation — that much is now gated), but the "1/3 corrected to ½ by `A_*`"
arithmetic can only be tested on the **second-moment operator** — exactly the R29 growing-spaces object. So **W2
(the tailor: connect via the rigorous 76.2 −2 pairing) must be built on the second moment**; the first-moment
version provably can't reach the ½. No fitting; exact analytic core-eigenvalue identity `1/(2χ(2)−1)` and stable
shifted |λ₂|=0.970 reported plainly; the W1 refutation reported as an honest negative, the W3 confirmation as the
positive. Fixed-level model caveat stated (proxy for the true growing operator; qualitative conclusions —
core-no-gap, shift-creates-gap — are robust; the 0.970 value is proxy-level).
