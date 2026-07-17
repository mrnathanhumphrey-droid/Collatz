# Phase 3: Canonical-relation to μ_n — the transfer gate

The surviving candidates from Phase 2 are **C5** (class-symmetry-broken two-branch IFS Φ' = {φ_+, φ_-, φ_2}) and, marginally, **C4** (parity-density Gaussian). The question this phase answers is:

**Does the self-conformal measure ν of Φ' have a canonical relation to the Syracuse stationary measure μ_n on Z_3 (or its limit on Z_3)?**

This is the gate where Probe 2 (SL_2 embedding) and Probe 3 (Cocycle Dolgopyat candidate b) both fell. The pre-registration favors H_DELTA_EXISTS_TRANSFER_BROKEN for that reason.

---

## (a) The natural lift: is ν → μ_n via projection?

Map κ: ℤ_3 → [0,1] via base-3 expansion: an element x = Σ a_i 3^i ∈ Z_3 maps to Σ a_i 3^{-i-1} ∈ [0,1].

For the **untwisted** base-3 IFS Φ_3 = {x ↦ (x+a)/3}, the self-conformal measure with Bernoulli weights p = (p_0, p_1, p_2) on {0,1,2} is the κ-pushforward of the corresponding Bernoulli measure on Z_3:

> ν_Φ_3,p = κ_∗ (Bernoulli(p)^⊗ℕ on Z_3)

But μ_n on Z_3 is **not Bernoulli** in the digit representation — it's the stationary measure of the Syracuse Markov chain on Z_3, which has correlations between digits introduced by the carry structure of the (3x+1)/2^{v_2} map. So even for the untwisted base-3 IFS, ν_Φ_3,p with any Bernoulli p ≠ μ_n's marginals would NOT equal κ_∗ μ_n.

**For C5's twisted IFS Φ' = {φ_+, φ_-, φ_2}, the self-conformal measure ν_Φ' lives on the attractor K_{Φ'} ⊆ [0,1].** The attractor K_{Φ'} ≠ [0,1] generically (it's a Cantor set if the cylinder unions are not space-filling), and the natural correspondence to Z_3 via base-3 expansion is **broken**: a point in K_{Φ'} is encoded by codings (a_1, a_2, ...) ∈ {+, −, 2}^ℕ, but the encoding is via the twisted maps φ_+, φ_-, φ_2, not via affine base-3 digits.

**Specifically: the encoding map Π: {+, −, 2}^ℕ → K_{Φ'} is**

> Π(a_1, a_2, ...) = lim_n φ_{a_1} ◦ ⋯ ◦ φ_{a_n}(0)

For the affine base-3 IFS this is exactly κ. For C5's twisted IFS, Π differs from κ by a (small) C² conjugacy h: [0,1] → [0,1] satisfying h ◦ φ_a = (untwisted-affine f_a) ◦ h for a ∈ {0, 1, 2} (where we identify {+, −, 2} with {0, 1, 2}).

**The conjugacy h is precisely the conjugacy ARHW Theorem 1.1's hypothesis FORBIDS.** If h exists C² and conjugates Φ' to Φ_3, then Φ' IS C² conjugate to a linear IFS — and ARHW Thm 1.1 does NOT apply to Φ'.

**Therefore:**
- If h exists ⟹ Φ' is conjugate to linear ⟹ Thm 1.1 does NOT fire (UNI fails) ⟹ no Fourier decay conclusion.
- If h does NOT exist ⟹ ν_Φ' is on a topologically-distinct attractor K_{Φ'}, and there is no canonical map ν_Φ' → μ_n.

This is the **structural pincer for canonical relation**. The non-linearity hypothesis (which is the hypothesis we need for ARHW) is exactly the obstruction to a canonical-conjugacy-based transfer.

---

## (b) Fourier-transform incompatibility — even granting a projection

Suppose, generously, we accept a non-canonical measurable map π: K_{Φ'} → Z_3 such that π_∗ ν_Φ' = μ_n (i.e. we choose to define one). Does real-Fourier decay of ν_Φ' transfer to 3-adic Fourier decay of μ_n?

ARHW Theorem 1.1 conclusion (verbatim, p. 2): |F_q(gν)| = O(1/|q|^α), where F_q is the **real Fourier transform** on ℝ:

> F_q(gν) = ∫ e^{-2πiqx} g(x) dν(x), q ∈ ℝ.

The 3-adic Fourier transform on Z_3 (or Z/3^n Z) uses characters:

> μ̂_n(ξ) = ∫_{Z/3^n} e^{2πi(ξ·y)/3^n} dμ_n(y), ξ ∈ Z/3^n Z.

The characters of Z_3 are **the dual group Q_3 / Z_3 ≅ Z[1/3] / Z**, with characters χ_q(x) = e^{2πi{qx}_3} where {·}_3 denotes 3-adic fractional part.

**These two transforms are incompatible.** Under the embedding κ: Z_3 ↪ [0,1] via base-3 expansion, a 3-adic character χ_q does NOT restrict to a real character e^{-2πi q'x}. Concretely:

- κ identifies x = Σ a_i 3^i with Σ a_i 3^{-i-1}.
- The 3-adic character χ_{1/3^k}(x) = e^{2πi a_{k-1}/3} depends on the (k-1)-th digit.
- The real character e^{-2πi q·κ(x)} = exp(-2πi q · Σ a_i 3^{-i-1}) is a **convolution over all digits** — does not factor through a single digit.

So even with a measurable bijection π, |F_q(ν_Φ')| → 0 for large real q tells us nothing about |μ̂_n(ξ)| for ξ in the dual group of Z/3^n Z.

**This is the exact same T1-transfer failure mode that closed Probe 2 (SL_2 embedding) and Probe 3 (cocycle Candidate b).** The framework computes decay on the wrong Fourier object.

The mechanism is structurally non-fixable by perturbation: real-Fourier and 3-adic Fourier are **different harmonic analyses on different topological groups**, and there is no general theorem converting decay of one to decay of the other (the only honest converter would be via uniform distribution / equidistribution mod 1, which doesn't preserve polynomial-rate quantitative bounds).

---

## (c) Base-3 expansion correspondence — does it bridge them?

The honest question to consider: is there a **partial** correspondence, e.g. via finite truncations Z/3^n Z embedded in [0, 1]?

Consider the truncation: μ_n on Z/3^n Z lifts to a measure μ_n^[0,1] on [0,1] supported on the 3^n equally-spaced points {j/3^n}. Its real Fourier transform at frequency q:

> F_q(μ_n^[0,1]) = Σ_j μ_n(j) · e^{-2πi q · j/3^n}

For q = ξ (integer) this becomes a discrete Fourier sum. **For q = ξ·3^n** this matches the 3-adic Fourier transform up to normalization. **For q = ξ·3^k with k < n** it captures a partial 3-adic frequency.

So real-Fourier decay of μ_n^[0,1] at SCALES q = ξ·3^k DOES correspond to 3-adic Fourier decay of μ_n at depth k. But:
1. ARHW gives decay on ν, **not on μ_n^[0,1]**. To transfer requires ν → μ_n^[0,1] canonically — and that's the broken canonical-relation gate from (a).
2. Even granting it, the polynomial rate |q|^{-α} at q = ξ·3^k gives a bound on 3-adic Fourier with exponent α — but the polynomial-in-A bound needed for Tao Prop 1.17 requires α as a function of the cyclic-structure size A = 3^n, which ARHW doesn't extract (per `COCYCLE_DOLGOPYAT_TRANSFER.md` (c)).

So the partial correspondence exists in principle, but the canonical-relation gate (a) still blocks it.

---

## (d) Summary: where each candidate falls

| Candidate | Has C² + UC + UNI? | Canonical ν → μ_n? | Verdict |
|---|---|---|---|
| C4: parity-Gaussian | yes (marginal) | NO — Gaussian convolution decorrelates from Syracuse structure | TRANSFER BROKEN |
| C5: class-broken two-branch | yes | NO — see pincer in (a) | TRANSFER BROKEN |

**Both surviving candidates fall at the canonical-relation gate.** This is exactly H_DELTA_EXISTS_TRANSFER_BROKEN.

The pincer in (a) is structural, not technical: the non-linearity hypothesis ARHW needs (h does NOT conjugate Φ' to Φ_3) is the same hypothesis that breaks the natural μ_n-projection. They are dual conditions: ν canonically relates to μ_n iff the IFS is conjugate to base-3 affine iff Thm 1.1 does NOT fire.

---

## (e) Adversarial check on Phase 3

(A3) **Transfer-mechanism honesty.** Confirmed broken via two independent routes: (a) canonical-conjugacy pincer, (b) real-vs-3-adic Fourier incompatibility. Both routes reproduce the Probe 2 / Probe 3 failure mode. We are not papering over.

(A4) **No terminology-overlap mistake.** ν as a self-conformal measure on [0,1] and μ_n as a Syracuse stationary measure on Z/3^n Z are different categories of object. Real Fourier ≠ 3-adic Fourier. This is verified, not assumed.

(A5) **Not repeating Probe 3 Candidate b.** C5 was honestly Syracuse-derived (T_lead's (1, 4) eigenvector) — but the derivation-fidelity passing does not save it from the transfer-gate failure. That confirms the strategic point: even a structurally Syracuse-encoding δ_a does not bridge the real-Fourier / 3-adic-Fourier gap, because the gap is **about the object the framework computes decay on**, not about the perturbation that gets us past the entry gate.
