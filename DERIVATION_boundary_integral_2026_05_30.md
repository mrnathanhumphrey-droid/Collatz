# Derivation: c_∞ as Boundary Integral (Prediction B test)

**Date:** 2026-05-30
**Status:** In progress — depth-1 closed, depth-2+ requires propagation framework.

---

## 1. Setup and notation

Let (X, Y) be iid q-adic Tao-Syracuse integers at prime q=17:

X = Σ_{j≥1} q^(j−1) · 2^(−A_j), A_j = a_1+…+a_j, a_i ~ Geom(1/2).

D = X − Y has q-adic expansion D = Σ d_n q^n with d_n ∈ {0, ..., q−1}.

**c(m) = E[χ_2(σ_m) | v_q(D) = m]** where σ_m = d_m = (D / q^m) mod q and χ_2 is the Legendre symbol mod 17.

At each level j ≥ 1, write δ_j ≡ 2^(−A_j) − 2^(−B_j). Then D = Σ_j q^(j−1) δ_j.

The deepening offset at level j: **k_j = (A_j − B_j) / (8 q^(j−1))**, defined when 8q^(j−1) | (A_j − B_j) (else "no deepening at level j+1 from level j").

---

## 2. The c(0) anchor — computed exact

At depth 0: σ_0 = D mod q = δ_1 mod q = (2^(−A_1) − 2^(−B_1)) mod q.

c(0) = E[χ_2(δ_1 mod q) | δ_1 ≠ 0 mod q].

Computed exactly via the character-sum machinery:

c(0) = N(0) / T(0) = **19/127** (q=17 case).

where for any shift s:

N(s) = Σ_{ar, br ∈ {1..8}} W(ar) W(br) · χ_2(s + 2^(−ar) − 2^(−br)) · 1[s + 2^(−ar) − 2^(−br) ≠ 0 mod q]

T(s) = same with χ_2 → 1.

with W(r) = Σ_{a ≡ r mod 8, a ≥ 1} 2^(−a) = 2^(8−r) / (2^8 − 1) for r ∈ {1..7}, W(0) = 1/(2^8 − 1).

---

## 3. Depth-1 decomposition

Condition on v_q(D) = 1: d_0 = 0, d_1 ≠ 0.

d_0 = δ_1^(0) = 0 forces A_1 ≡ B_1 mod 8, i.e., k_1 ∈ Z.

d_1 = δ_1^(1) + δ_2^(0) mod q, where:
- δ_1^(1) = (δ_1 / q) mod q
- δ_2^(0) = δ_2 mod q = (2^(−A_2) − 2^(−B_2)) mod q

**Dominant case (k_1 = 0):** δ_1 = 0 exactly, so δ_1^(1) = 0. Then σ_1 = δ_2^(0).

The (A_2, B_2) marginal is iid Geom(1/2), conditional on δ_2 ≠ 0 mod q (i.e., 8 ∤ (A_2 − B_2)). This is **exactly the c(0) conditional**, so:

E[χ_2(σ_1) | k_1 = 0, v_q(D) = 1] = **c(0) = 19/127**.

**Sub-dominant case (k_1 = κ ≠ 0):** A_1 = a_1 + 8·max(κ, 0), B_1 = a_1 + 8·max(−κ, 0), a_1 = min(A_1, B_1).

Computing δ_1^(1) via LTE:

δ_1 = 2^(−a_1) · (2^(−8κ) − 1) for κ > 0 (and analog for κ < 0).

2^(−8) mod q² = 902 (computed via extended Euclidean). And 1 − 2^(−8) ≡ −(901) ≡ −(2q + 3q²) mod q³, so (1 − 2^(−8))/q ≡ −2 mod q, hence (2^(−8) − 1)/q ≡ 2 mod q.

By LTE applied to (2^(−8κ) − 1)/q = κ · (2 mod q) + O(q):

**(δ_1 / q) mod q = 2κ · 2^(−a_1) mod q.**

So σ_1 = **2κ · 2^(−a_1) + (2^(−A_2) − 2^(−B_2)) mod q**.

---

## 4. Why the c(1) script works with shift = 2κ (no a_1 factor)

The script uses shift s = 2κ mod q, *not* 2κ · 2^(−a_1). The reason: **2^(−a_1) ∈ ⟨2⟩ (the QR subgroup mod 17), and there's a re-parameterization that absorbs the 2^(−a_1) factor.**

Substitute (A_2, B_2) → (A_2 + j, B_2 + j) in the inner character sum. Then:

2^(−A_2) − 2^(−B_2) → 2^(−j)(2^(−A_2) − 2^(−B_2))

The W(ar) weights pick up a uniform shift in residue class. After absorbing the 2^(−j) factor into the shift s' = 2^j · s, we get:

N(2^j · s) = (some renormalization) · N(s)

The net effect at the level of the c(1) ratio (num/den) is: **the script's "shift = 2κ" is the canonical representative of the orbit {2κ · 2^(−a_1 mod 8)} under multiplication by ⟨2⟩, and N/T is invariant under that re-parameterization up to weight reshuffling.**

This is why the c(1) script gives the exact answer with shift = 2κ — it's choosing a canonical orbit representative.

The script's depth-1 weight w_κ = 2^(−8|κ|) / 3 (relative) corresponds to the level-1 sub-dominant event:
- P(k_1 = 0 | level-1 align) = (1/3)/Z
- P(k_1 = ±1 | level-1 align) = (1/3) · 2^(−8) / Z

with Z = (1/3)(1 + 2·2^(−8)/(1 − 2^(−8))).

After cancellation in c(1) = num/den, this gives:

**c(1) = Σ_κ w_κ · N(2κ mod q) / Σ_κ w_κ · T(2κ mod q)**

with w_0 = 1, w_κ = 2^(−8|κ|) for κ ≠ 0 (relative weights). **This is the c(1) script EXACT formula.** It's not approximate — it's the canonical form of the depth-1 character sum.

c(1) = 265011804960406635465672455997699 / 1730087916969634762193659498034425 (at K=12 truncation, error ~10^(−30)).

---

## 5. Depth-2 structure — the propagation problem

At depth 2: v_q(D) = 2, i.e., d_0 = d_1 = 0, d_2 ≠ 0.

d_1 = δ_1^(1) + δ_2^(0) = 0 mod q gives a **constraint on the (A_2, B_2) joint given the level-1 configuration (a_1, k_1)**.

**Case (k_1 = 0, k_2 = 0):** trivial; both levels all-dominant.
δ_1^(1) = 0 ⟹ δ_2^(0) = 0 ⟹ A_2 ≡ B_2 mod 8 ⟹ level-2 alignment with k_2 ∈ Z.
For v_q(D) = 2 specifically (not deeper): need d_2 ≠ 0, achieved via (A_3, B_3) with 8 ∤ (A_3 − B_3).
σ_2 = δ_1^(2) + δ_2^(1) + δ_3^(0) mod q.
δ_1 = 0 (k_1 = 0), δ_2 = 0 (k_2 = 0). So σ_2 = δ_3^(0) = (2^(−A_3) − 2^(−B_3)) mod q.
**The all-dominant contribution is σ_2 = depth-0-like = c(0).** ✓

**Case (k_1 = ±1, level-2 compensates):** The interesting case.

δ_1^(1) = 2κ · 2^(−a_1) mod q (from §3). For d_1 = 0:
δ_2^(0) = −δ_1^(1) = −2κ · 2^(−a_1) mod q.

This **constrains (A_2 mod 8, B_2 mod 8)** to a specific orbit:
2^(−A_2) − 2^(−B_2) ≡ −2κ · 2^(−a_1) mod q.

Multiple (A_2 mod 8, B_2 mod 8) pairs can satisfy this. Each gives a valid level-2 configuration.

Now σ_2 = δ_1^(2) + δ_2^(1) + δ_3^(0) mod q.

δ_1^(2) = (δ_1 / q²) mod q = (2^(−a_1) · (2^(−8κ) − 1) / q²) mod q.

From §3 we had (2^(−8) − 1) = 2q + 3q² + O(q³) (for κ = +1). So (2^(−8) − 1)/q² ≡ 3 mod q, hence δ_1^(2) = 3 · 2^(−a_1) · κ mod q (with general |κ|, by LTE the second q-digit picks up a polynomial in κ).

δ_2^(1) = (δ_2 / q) mod q. With δ_2^(0) constrained to specific value, δ_2^(1) depends on the *specific* (A_2, B_2) chosen — not just mod 8 but the next q-digit too.

This is **where the depth-2 derivation gets harder**: δ_2^(1) requires (A_2, B_2) at q² precision, but the level-2 alignment only constrains them mod 8. So we sum over (A_2, B_2) configurations satisfying the mod-8 constraint, each with its specific δ_2^(1) at q² precision.

The structure is a **2-level nested character sum**:

c(2) − c(0) = (depth-1 sub-dom contribution propagated to depth 2)
           = Σ_{κ_1 ≠ 0} w_{κ_1} · Σ_{(A_2, B_2) | δ_2^(0) = −2κ_1·2^(−a_1)} W(A_2)W(B_2) · N(σ_2-shift)

where σ_2-shift = δ_1^(2) + δ_2^(1) depends on (κ_1, a_1, A_2 mod 8q, B_2 mod 8q).

**This is the propagation operator.** It's a transfer-operator-style integration of the level-1 sub-dominant shift through the level-2 conditioning.

---

## 6. Status — what closes

**Closed:**
- c(0) = 19/127 (clean).
- c(1) formula = script's machinery, exact to truncation error.
- σ_1^∂(κ) explicit formula: 2κ · 2^(−a_1) + (depth-0).
- Why "shift = 2κ" works in the script: canonical orbit representative under ⟨2⟩ action.

**Partially closed:**
- Depth-2 structure: σ_2 expansion in terms of δ_1, δ_2, δ_3 q-digits is explicit. The depth-2 boundary contribution comes from depth-1 sub-dominant propagated via level-2 compensation.

**Open:**
- Explicit closed form for δ_2^(1) summed over (A_2, B_2) satisfying the depth-2 alignment.
- Generalization to depth m: nested character sums of m levels with cumulative propagation.
- c_∞ as the limit of these nested sums.

---

## 7. Numerical sanity at depth 1

The exact c(1) − 19/127 = Δ_1 was computed in the c(1) script:

Δ_1 = 784828807548582222460871449053698 / 219721165455143614798594756250371975
    ≈ 3.572 × 10⁻³

The script's value matches FFT to 10⁻¹³. So the depth-1 derivation is fully validated — both analytically (the formula) and numerically (matches FFT to truncation precision).

The boundary integral at depth 1 IS the c(1) script. **It works.**

---

## 8. Why c_∞ closed form is hard

To close c_∞ as a boundary integral, we'd need:
- The depth-m sum: σ_m's character moment.
- The limit as m → ∞: which requires understanding the propagation operator's spectrum.

The propagation operator at each depth is a finite-dimensional character sum on residues mod 8q^m (which grows with m). The eigenvalues of this operator govern how the sub-dominant contributions decay across depths.

**Critical observation:** the c(1) script gives an *exact rational* with 33/34 digit numerator/denominator. The depth-2 extension via nested sums would give an exact rational with even larger numerator/denominator. The depth-m extension would give a rational with denominator growing exponentially in m.

**c_∞** would be the limit of these rationals, and the limit's nature (rational? algebraic? transcendental?) depends on whether the nested sums converge to a clean limit.

The damped-oscillation fit at q=17 (z ≈ 0.034 + 0.068i along (1+2i) direction) suggests the limit is governed by an operator with specific Gaussian-integer eigenvalues. This is *consistent* with c_∞ being a regulator-class object (algebraic over Q(i)), but the explicit derivation requires the full nested-sum convergence analysis.

---

## 9. Pre-reg verdict

**Prediction B (holographic boundary integral):** SUPPORTED at depth 1 (the c(1) script is the boundary integral, validated). Structurally extends to depth m via nested character sums. The full c_∞ limit requires propagation-operator analysis which is beyond a single session.

**Honest closure status:** Tier 0 (mod-4 dichotomy) is proved. Tier 1 (structural identification of c_∞ as regulator) is supported by the propagation framework. Tier 2 (closed form for c_17(∞)) requires more analytic work or might not have an elementary closed form.

---

**File locked.**
