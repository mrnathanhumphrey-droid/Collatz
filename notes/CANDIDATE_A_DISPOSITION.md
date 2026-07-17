# CANDIDATE_A_DISPOSITION — Reading A scoping probe minimum-viable test outcome

**Date:** 2026-05-12. Wilson (analyst) reporting to Nathan. Top-level disposition.

---

## DISPOSITION: **H_CANDIDATE_A_FALSIFIES_F2**

> **The W_k multi-resolution filtration on L²(Ẑ_3^×) does not carry rate-1/2 via the bilinear-pair-form moment functional φ_n.** Exact-Q computation of c_{n,k} := ⟨φ_n, lift_n(R_k)⟩ for all (n, k) with 0 ≤ k < n ≤ 6 shows c_{n, k} = 0/1 exactly for every k < n − 1. The entire moment functional φ_n lives in the single finest-scale subspace W_{n−1}, with c_{n, n−1} = S_n converging to 7/15 (not decaying to zero).
>
> Pre-registered patterns A1 (dominant-k rate-0.5) and A2 (signed-sum cancellation rate-0.5) are both **decisively rejected**. Pattern F2 (single-level k = n − 1 dominance) holds in the strongest possible form — as an exact rational identity, not an approximate one.
>
> Per the scoping recommendation, **route to Candidate B (Kozyrev p-adic wavelets) as the next probe**. The L²(Ẑ_3^×) framing survives — the W_k basis is the wrong decomposition for this moment functional, but a different basis on the same Hilbert space (Kozyrev's wavelet basis, which mixes dilation with translation) may carry the rate.

---

## Decision rule applied

Per the pre-registration decision tree:

- **H_CANDIDATE_A_FALSIFIES_F2:** Single-level k = n − 1 dominance. The pre-registration explicitly noted: "Conflicts with R77.4; either R77.4 needs re-examination or the framework genuinely lives in single-level structure that R77.4 missed."

In this case, **no conflict with R77.4** is incurred. R77.4 ruled out single-level rate-1/2 in the *spectrum* of K_n (within-level Markov operator). This probe's F2 says φ_n lives in W_{n−1} (a property of the moment functional's W_k decomposition, not of K_n's spectrum). The two findings are independent and both stand.

---

## Load-bearing data (the c_{n,k} table)

Exact rationals (full table in `candidate_a_c_nk_table.csv`):

| n | k = 0 | k = 1 | k = 2 | k = 3 | k = 4 | k = 5 |
|---|-------|-------|-------|-------|-------|-------|
| 1 | 1/6   |       |       |       |       |       |
| 2 | **0/1**   | 10/21 |       |       |       |       |
| 3 | **0/1**   | **0/1**   | 31370/67963 |       |       |       |
| 4 | **0/1**   | **0/1**   | **0/1**   | S_4 (large) |       |       |
| 5 | **0/1**   | **0/1**   | **0/1**   | **0/1**   | S_5 (60+ digits) |       |
| 6 | **0/1**   | **0/1**   | **0/1**   | **0/1**   | **0/1**   | S_6 (217+ digits) |

15 of 21 entries are exactly 0/1 over Q. The 6 nonzero entries are all on the diagonal k = n − 1, equal to S_n exactly (since ⟨φ_n, π_∞^{(n)}⟩ = 0 for n ≥ 2).

### Ratios across n = 2..6 — A1 / A2 rejection

|c_{n, k*(n)}| / |c_{n−1, k*(n−1)}| ratios:

| n | ratio    | target if A1/A2 |
|---|----------|------------------|
| 2 | 2.857143 | 0.5              |
| 3 | 0.969307 | 0.5              |
| 4 | 1.005719 | 0.5              |
| 5 | 1.002802 | 0.5              |
| 6 | 1.001405 | 0.5              |

Ratios converge to **1.0**, not 0.5. The dominant-k value c_{n, n−1} = S_n is converging to a *constant* (7/15), not decaying. The deviation S_n − 7/15 = ε_n is what decays at rate ≈ 0.5 — but that's the project's known rate-1/2 phenomenon (R77.6), not a new W_k-filtration finding.

---

## Why F2 holds — structural diagnosis

The bilinear-pair-form moment functional has the explicit form

  φ_n(r) = Σ_s π_n(s) K_n(r − s)

with kernel K_n(d) := Σ_{ξ ∈ Z/3^n, 3∤ξ} e^{−2πi d ξ / 3^n}, which evaluates to:

  K_n(0) = 2 · 3^{n−1};  K_n(±3^{n−1}) = −3^{n−1};  K_n elsewhere = 0.

So K_n is supported on the cosets of 3^{n−1} Z mod 3^n. Hence

  φ_n(r) = 3^n · (π_n(r) − π̄_n(r))

where π̄_n(r) is the 3-fiber-average of π_n at scale 3^{n−1}. **φ_n has zero 3-fiber-mean at scale 3^{n−1}, so φ_n ∈ W_{n−1} by definition.**

Inner product against any T-lift from a coarser level k < n − 1 is exactly zero by the structural orthogonality W_{n−1} ⊥ T^{j}(W_k) for k ≠ n − 1 (R77.5 §3).

The exact-Q zero is automatic from this structural fact; it's not a contingent feature of the Markov dynamics.

---

## Recommendation: route to Candidate B (Kozyrev wavelets)

The Reading A scoping recommendation (READING_A_SCOPING_RECOMMENDATION.md) pre-specified the F2/F3 outcome as the trigger for Candidate B. The case for Candidate B:

1. **Kozyrev wavelets diagonalize joint translation + dilation** on L²(Q_p), unlike the W_k filtration which captures only dilation. φ_n, while in W_{n−1}, has non-trivial *translation* structure within W_{n−1}. The Kozyrev basis might localize this translation structure to specific wavelet indices, and rate-1/2 might appear as the decay of |⟨φ_n, ψ_{j*, n_*, ε_*}⟩| as n grows.

2. **R77.6's branch-cut framing at z = 2** suggests rate-1/2 is a spectral-density feature, not an eigenvalue. Kozyrev wavelets diagonalize the Vladimirov fractional-derivative operator D^α, whose spectrum is continuous — making it a natural diagonalizer for spectral-density-type phenomena.

3. **The L²(Ẑ_3^×) framing survives.** The probe didn't fail at C1 (decomposition sanity): Σ_k c_{n,k} = ⟨φ_n, π_n − π_∞^{(n)}⟩ exactly. The W_k filtration is a valid orthogonal decomposition; it's just that φ_n is concentrated entirely in the highest-frequency band of that decomposition. Switching to Kozyrev gives a different decomposition where φ_n may spread across multiple wavelet indices.

### Honest caveat

It is possible that Kozyrev wavelets *also* concentrate φ_n at the finest scale (since φ_n's frequency support is exactly the coprime-to-3 frequencies at level 3^n, which is precisely the finest-scale wavelet index in Kozyrev's basis). If so, the W_k = F2 finding generalizes to "any natural orthonormal basis on L²(Ẑ_3^×) concentrates φ_n at the finest scale," and the rate-1/2 phenomenon is fundamentally **not** a wave-localization feature.

In that case, the next probe would route to **Candidate C (transfer-operator on Syracuse coherent extension)**, but with the prior step-(a) gap to bridge (specifying the Syracuse extension to Ẑ_3 — see READING_A_SCOPING_CANDIDATES.md §C.a).

The scoping cost for Candidate B is articulated as ~5-6 weeks of construction work after a confirming probe; the *probe itself* is roughly 1 focused session of the same kind that completed Candidate A's probe, with the additional setup cost of Kozyrev basis arithmetic over Q (requires roots-of-unity arithmetic, ~1-2 sessions of basis construction first).

### Alternative recommendation: pause to reconsider φ_n articulation

Before committing to Candidate B, one could re-examine whether the linearized φ_n used here is *the right φ_n* for the spectral-completion question. R76 introduces φ_n abstractly as "the bilinear pair-form moment functional"; the probe used φ_n = K_n * π_n (the natural linearization, depending on π_n itself). An alternative articulation: define φ_n as the **gradient of S_n with respect to π_n**, which is 2 · K_n * π_n (factor of 2 from the quadratic form). This rescales c_{n,k} by 2 but does not change the F2 finding.

A more substantive alternative: use **R76's class-resolved decomposition** (R76 §11) where π_n splits into π_{+, n} + π_{−, n} on the two mod-3 classes, and the (1, 4)-eigendirection is the rate-1/2 carrier. That framework is finite-dimensional (2D at each n) and might be a better testing ground for "where does rate-1/2 live spectrally" than the L²(Ẑ_3^×) framework. The R76 §11 framework was the original "operator at finite truncation" framing of R77.2; R77.4 ruled out the K_n form of that operator, but the (P_+, P_−) 2D recursion remains an active anchor.

If the next move is to revisit R76 §11's class-resolved (P_+, P_−) 2D recursion as the spectral object, the scoping work shifts away from Reading A entirely. This is a viable alternative to the Reading-A → Candidate-B trajectory.

---

## Summary for Nathan (orientation only — not load-bearing)

The probe was a clean run, computed exactly over Q in one session (~10 min total: 7 min for π_6 stationary plus ~30 s for everything else). The data is unambiguous:

- **Phase 1 unblocked.** φ_n construction worked, ε_n reproduced exactly at n = 1..6.
- **Phase 2 unblocked.** Lift orthogonality verified exactly over Q at all 34 pairs.
- **Phase 3 unambiguous.** 15 of 21 c_{n,k} are exactly 0/1 over Q.
- **Phase 4 unambiguous.** F2 holds; A1, A2, F1 rejected.

The disposition has no internal ambiguity: this is not a "noisy data, ambiguous" situation. The c = 0 over Q is *structural* — it falls out of the support of K_n's kernel function on 3^{n−1}-cosets.

The scoping probe trajectory now becomes:
- **M_3** (R77.3 — falsified) → **R_K** (R_K probe — intractable) → **W_k filtration via φ_n** (this probe — F2 falsified). Three negative findings, each carving the "where rate-1/2 lives" question more precisely.
- Next probe: **Candidate B Kozyrev wavelets**, with the honest caveat that the finest-scale concentration of φ_n may generalize to any natural orthonormal basis. If it does, route to Candidate C with the step-(a) extension specification.

**The framework remains rate-1/2 = branch-cut singularity at z = 2 (R77.6).** Each probe rules out a candidate operator/decomposition; none has yet found a positive carrier.

---

## File index

- `CANDIDATE_A_PHI_CONSTRUCTION.md` — Phase 1 explicit construction
- `CANDIDATE_A_LIFT_CONSTRUCTION.md` — Phase 2 lift apparatus
- `CANDIDATE_A_C_NK.md` — Phase 3 c_{n,k} table
- `CANDIDATE_A_PATTERN_MATCH.md` — Phase 4 pre-registered pattern outcomes
- `CANDIDATE_A_DISPOSITION.md` — this file
- `candidate_a_compute.py` — reproducible script (Q-arithmetic only, no external packages)
- `candidate_a_c_nk_table.csv` — exact-rational table
- `candidate_a_phi_n_verify.csv` — verification table
- `candidate_a_lift_orthogonality.csv` — orthogonality table
