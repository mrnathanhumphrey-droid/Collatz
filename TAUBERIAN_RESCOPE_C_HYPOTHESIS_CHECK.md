# TAUBERIAN_RESCOPE_C_HYPOTHESIS_CHECK (Korevaar 2002 Wiener-Ikehara + Newman × inputs)

**Date:** 2026-05-13.

Two theorem readings: 4.2 (W-I full form) and 6.1/8.1 (Newman-Korevaar Dirichlet form).

---

## h × I matrix — Theorem 4.2 (Wiener-Ikehara, Laplace-Stieltjes form)

Constructions needed: a nonnegative, nondecreasing S(t) and an A and an analytic-continuation condition.

| Hypothesis | (1) ε_k | (2) C1 | (3) C2 | (4) BT |
|---|---|---|---|---|
| h_1: S(t) = 0 for t < 0 | Construct S(t) := Σ_{k ≤ e^t} |ε_k| · 2^k (cumulative normalized magnitudes); S(t) = 0 for t < 0. SATISFIED. | N/A | N/A | N/A |
| h_2: S(t) ≥ 0 for t ≥ 0 | SATISFIED by |ε_k|·2^k ≥ 0. | | | |
| h_3: S(t) nondecreasing | SATISFIED — cumulative sum of nonnegative terms is monotone nondecreasing. | | | |
| h_4: Laplace-Stieltjes transform f(z) = ∫ e^{-zt} dS(t) exists for Re z > 1 | UNVERIFIABLE — dS(t) is supported on logarithmic spacing t = log k, so f(z) = Σ_k |ε_k|·2^k · k^{−z}. Convergence for Re z > 1 requires |ε_k|·2^k = O(k^{1−ε}). Empirically |ε_k|·2^k is O(0.2) for k ≤ 8; the *full* sequence behavior is unknown — but even if bounded, this only gives convergence for Re z > 1 (consistent). h_4 plausibly SATISFIED if |ε_k|·2^k is bounded. | (2) C1 renewal-walk gives no direct argument about the Dirichlet series. | (3) BMP F_1 diffraction is on ℝ as Fourier transform, different object. | (4) BT archimedean place: bears on whether the Dirichlet series sees the closure. |
| h_5: g(z) = f(z) − A/(z−1) has continuous extension to Re z ≥ 1 | UNVERIFIABLE — this is the **load-bearing analytic-continuation condition** of Wiener-Ikehara. It is essentially the closure target dressed as a Dirichlet-series analytic property. For Syracuse μ_n / ε_k there is no proven Dirichlet-series functional equation. **Mode H circular: the condition we'd derive (decay of partial sums) requires the analytic continuation of the Dirichlet series we're trying to use.** | | | |

Output for SELECTED case would be: e^{-t} S(t) → A, i.e. Σ_{k ≤ e^t} |ε_k|·2^k · e^{-t} → A. This is a 1-term asymptotic for the *partial sum*, not for ε_k itself.

**Theorem 4.2 disposition: BLOCKER (h_5 Mode H circular) PARTIAL (everything else satisfiable).**

---

## h × I matrix — Theorem 6.1 / 8.1 (Newman-Korevaar Dirichlet)

Reading f(z) = Σ a_n n^{-z} for a_n = ε_n (or |ε_n|·2^n).

| Hypothesis | (1) ε_k | (2) C1 | (3) C2 | (4) BT |
|---|---|---|---|---|
| h_1': |a_n| ≤ C bounded | UNVERIFIABLE for full sequence; satisfied for k ≤ 8 (|ε_n| ≤ 0.2; |ε_n|·2^n ≤ 0.4). | N/A | N/A | N/A |
| h_2': f(z) = Σ a_n n^{-z} converges for Re z > 1 | UNVERIFIABLE — requires |a_n| bounded *and* the Dirichlet series convergence at Re z > 1. For a_n = ε_n the latter requires |ε_n| ≤ C n^{1−ε}, weaker than boundedness — plausibly SATISFIED. For a_n = |ε_n|·2^n, requires |ε_n|·2^n ≤ C n^{1-ε}, plausible from observed data. | | | |
| h_3': g(z) = f(z) − A/(z-1) has holomorphic extension to Re z ≥ 1 | UNVERIFIABLE — Mode H circular. **Same load-bearing condition as Theorem 4.2.** | | | |

Output for SELECTED case would be: Σ a_n / n converges to f(1) − A. This is even *weaker* than 4.2's S(t)/e^t → A — it just says the partial-sum-at-z=1 series converges, with no rate.

**Theorem 6.1 disposition: BLOCKER (h_3' Mode H circular).**

---

## Mode H summary for C

The Wiener-Ikehara and Newman-Korevaar theorems both require knowledge that the relevant Dirichlet (or Laplace-Stieltjes) transform admits **continuous/holomorphic extension to Re z ≥ 1** — equivalently, that the generating function has a controlled pole-only-at-z=1 singular structure on the boundary line.

For Syracuse μ_n / ε_k there is NO proven Dirichlet-series functional equation, no proven analytic continuation past Re z = σ_a, and the "no zeros on Re z = 1" property analogous to ζ(s) does not have any known analog for the Syracuse-derived Dirichlet series.

This is a clean Mode H circularity: the analytic property used as Tauberian *hypothesis* is essentially the target object's polynomial-in-A regularity.

**Aggregate disposition for C: BLOCKER (Mode H circular).**
