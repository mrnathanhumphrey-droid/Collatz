# CANDIDATE_A_PATTERN_MATCH — pre-registered pattern matching

**Date:** 2026-05-12. Wilson (analyst) reporting to Nathan. Phase 4 of the Reading A scoping probe.

---

## Per-pattern disposition

### A1 — Dominant-k rate-0.5 decay

**Disposition: NO.**

Pre-reg statement: "For each n, there exists a dominant k = k*(n) such that |c_{n, k*(n)}| dominates the sum Σ_k |c_{n, k}|, AND the dominant value scales like (1/2)^n as n grows from n = 1 to n = 6."

- First clause **holds in the strongest possible form**: k*(n) = n − 1 carries 100% of Σ_k |c| (all other c are exactly zero over Q).
- Second clause **decisively rejected**: |c_{n, k*(n)}| = S_n which converges to 7/15 ≈ 0.4667, NOT decaying. Ratios across consecutive n converge to 1.0, not 0.5.

The dominant-k structure exists, but its value does **not** scale like (1/2)^n — it scales like (constant + O((1/2)^n)).

### A2 — Phase-cancellation signed-sum rate-0.5

**Disposition: NO.**

Pre-reg statement: "The signed sum S_n := Σ_k c_{n,k} (without absolute values) shows rate-0.5 decay even though individual |c_{n,k}| terms may not."

- All nonzero c_{n, k} are positive (sign pattern: "0 0 ... 0 +"). No cross-k sign opposition. No cancellation mechanism.
- Σ_k c_{n,k} = S_n → 7/15. No decay to zero. Ratios converge to 1.0, not 0.5.

A2 is decisively rejected. The pre-registration's anticipated mechanism (telescoping signed contributions across W_k levels, motivated by R77.6's branch-cut at z = 2) **does not appear in the c_{n,k} structure**.

### F1 — All c_{n,k} approximately zero

**Disposition: NO.**

c_{n, n−1} is large and nonzero at every n ≥ 2. Specifically, c_{n, n−1} ≈ 0.46 ≈ 7/15, comparable to S_n itself. So F1 is rejected.

### F2 — c_{n,k} dominated by k = n − 1 (single-level effect)

**Disposition: YES — in the strongest possible form.**

c_{n, k} = 0/1 *exactly* for every k with 0 ≤ k < n − 1, at every n ∈ {2, 3, 4, 5, 6}. The single-level k = n − 1 carries 100% of the contribution; no other W_k subspace contributes anything.

The dominance is not merely approximate or asymptotic — it's an **exact rational identity** at every n in the test range.

### F3 — c_{n,k} scales without rate-0.5 signature

**Disposition: PARTIALLY YES (subsumed by F2).**

The single nonzero c_{n, n−1} scales at **rate 1.0 toward the limit 7/15**, not at rate 0.5. So in the sense of "the pattern exists but at a rate other than 1/2," F3 is also a match. But this is structurally a restatement of F2: there's only one nonzero c per level, it converges to a constant, and the deviation from that constant carries rate ≈ 0.5 (which is ε_n itself, by construction).

---

## Overall disposition: H_CANDIDATE_A_FALSIFIES_F2

**The W_k multiresolution decomposition does NOT distribute φ_n's overlap with π_n − π_∞ across multiple scales. The entire bilinear-pair-form moment lives in the single finest-scale subspace W_{n−1} at every level.**

The pre-registration's hypothesized "rate-1/2 emerges from cross-W_k structure" picture is rejected: there is no cross-W_k structure to emerge from. φ_n is entirely in W_{n−1}.

### Why F2 holds — structural diagnosis

This is **not** a quirk of the Markov chain dynamics. It's a structural property of the bilinear-pair-form moment functional φ_n itself:

  K_n(d) is supported on d ≡ 0 mod 3^{n−1}

  ⇒ φ_n(r) = Σ_s π_n(s) K_n(r − s) has zero 3-fiber-mean at scale 3^{n−1}

  ⇒ φ_n ∈ W_{n−1} (by definition of W_{n−1} as functions with mean-zero 3-fibers at the finest scale)

  ⇒ ⟨φ_n, lift_n(R_k)⟩ = 0 for every R_k that lies in a coarser scale (T^{n−k−1}(W_k) for k < n−1 sits in scales 0, 1, ..., n−2, all orthogonal to W_{n−1})

  ⇒ Only k = n−1 contributes.

This is a consequence of the specific form of S_n as a bilinear moment with frequency support on {ξ : 3∤ξ in Z/3^n} — the frequency restriction concentrates the moment functional at the finest scale.

### Reconciling with R77.4 (single-level was supposedly ruled out)

R77.4 ruled out that **K_n itself** (the level-n Markov operator) has rate-1/2 in its spectrum. F2 here says something different: it says **φ_n** (the moment functional) lives entirely in W_{n−1}, not that K_n has a rate-1/2 eigenvalue at level n−1. These are consistent claims:

- K_n's spectrum: no eigenvalue at 1/2 (R77.4, ruled out).
- φ_n's W_k decomposition: entirely in W_{n−1} (this probe).
- ε_n's evolution: rate ≈ 0.5 (project anchor, R77.6 supports branch-cut singular structure).

The three are compatible because rate-1/2 lives in the n-evolution of *⟨φ_n, π_n − π_∞⟩* itself, which is now seen to equal ⟨φ_n, P_{W_{n−1}}(π_n − π_∞)⟩ — i.e., the projection of π_n − π_∞ onto the moving finest-scale subspace W_{n−1} at each level. As n grows, W_{n−1} moves to finer scales, and π_n's overlap with that moving subspace decays at rate ≈ 0.5.

This is consistent with R77.6's branch-cut framing: there's no fixed Hilbert subspace of L²(Ẑ_3^×) carrying the rate; the rate lives in the *flow* of the moment functional through the multiresolution filtration as n → ∞.

---

## Adversarial check outcomes

**(A1) φ_n construction fidelity.** ⟨φ_n, π_n⟩ = S_n verified exactly at n = 1, 2, 3 against R76 values, and ε_n = S_n − 7/15 reproduced for n = 1..6 against R77.6 values. ✓

**(A2) Lift orthogonality verification.** ⟨lift_n(R_{k1}), lift_n(R_{k2})⟩ = 0 exactly over Q for all 34 cross-pairs (n, k1, k2). ✓

**(A3) Exact rationals throughout.** `fractions.Fraction` end-to-end. The c_{n,k} = 0/1 results are exact rational equalities, not numerical roundoff. ✓

**(A4) A1 vs A2 honest disambiguation.** Both rejected. No motivated reasoning to declare a confirm. ✓

**(A5) Range honesty.** Pattern holds exactly through n = 6 (exact-Q c = 0 for k < n−1 at n = 2, 3, 4, 5, 6 — five independent levels). The probe does NOT extrapolate beyond n = 6, but the pattern's exact-rational nature (c = 0/1, not c ≈ 0) means it's a structural identity, not a finite-budget effect. ✓

---

## Disposition handoff

H_CANDIDATE_A_FALSIFIES_F2 routes the project to:

1. **Candidate B (p-adic wavelets, Kozyrev)** per the scoping recommendation. Kozyrev wavelets provide a *different* multiresolution decomposition on L²(Ẑ_3^×) where translation and dilation act jointly; the W_k filtration here is purely the dilation axis. If φ_n has a richer decomposition in the Kozyrev basis (where dilation is mixed with translation modes), the rate-1/2 might localize there.

2. **Reconcile with R77.6's branch-cut framing.** The F2 result is *consistent* with rate-1/2 being a branch-cut feature of the generating function rather than an eigenvalue-localization phenomenon. The W_k filtration is not "the wrong basis for an eigenvalue" — it's "any basis, applied to a non-eigenvalue feature." Candidate B may or may not help.

3. **Note for R77.4 erratum.** R77.4 ruled out the *spectrum* of K_n having an eigenvalue near 1/2. This probe doesn't contradict R77.4 — it shows the moment functional structure, not the operator spectrum. Both findings stand.

---

## Files

- `candidate_a_compute.py`
- `candidate_a_c_nk_table.csv`
- `candidate_a_phi_n_verify.csv`
- `candidate_a_lift_orthogonality.csv`
