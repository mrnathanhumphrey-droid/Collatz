# READING_A_SCOPING_MIN_VIABLE_TEST — Phase 2 pre-registration-ready test for Candidate A

**Date:** 2026-05-11. Wilson (analyst) reporting to Nathan. Phase-2 deliverable of the Reading A scoping probe.

Candidate selected: **(A) Hilbert spaces of locally constant functions on Ẑ_3^×** (highest tractability per Phase 1).

---

## Goal of the minimum-viable test

Determine, empirically and over Q, whether the multi-resolution decomposition

```
L²(Ẑ_3^×, μ) = W_0 ⊕ W_1 ⊕ W_2 ⊕ ...    (orthogonal direct sum in 3-adic Haar L²)
```

carries the rate-1/2 envelope of ε_n through the moment-functional decomposition

```
ε_n = ⟨φ_n, π_n − π_∞⟩  =  Σ_{k ≥ 0} ⟨φ_n, lift_n(P_{W_k}(π_n − π_∞))⟩.
```

The test does NOT attempt to derive Φ_∞'s spectrum or rate-1/2 from first principles. It tests whether the **per-W_k contributions** to ε_n exhibit the structural pattern that the Candidate-A framework would predict.

---

## What's computed

### Definitions (all over Q via fractions.Fraction)

Let π_n ∈ V_n be the level-n stationary distribution from R77.5's existing machinery (`build_markov_rational`, `stationary_rational`).

Let R_k ∈ W_k ⊂ V_{k+1} be the lift residual from R77.5 (`pi_dict`, `lift_pi`, R_k = π_{k+1} − T(π_k)).

For n ≥ k+1, let `lift_n(R_k)` ∈ V_n be the iterated lift of R_k to level n (applying T_{j→j+1} successively from level k+1 to level n).

Let φ_n be the moment functional from R76 such that ε_n = ⟨φ_n, π_n⟩ − 7/15 = ⟨φ_n, π_n − π_∞⟩ (the bilinear pair-form moment). It's a function `φ_n : coprime states in Z/3^n → Q`, already exact over Q in R76's existing code.

Define for each (n, k) with 0 ≤ k < n:

```
c_{n,k} := ⟨φ_n, lift_n(R_k)⟩ ∈ Q.
```

(Inner product over the level-n coprime states with uniform counting measure, exact rational.)

### Quantities tabulated

For n = 2, 3, 4, 5, 6 and each valid k = 0, 1, ..., n−1:

1. **Per-(n,k) raw contribution:** c_{n,k} (exact Q, float).
2. **Decomposition check:** Σ_{k=0}^{n−1} c_{n,k} should equal ε_n + (rest term from level k=0 boundary handling). This is a sanity check that the W_k decomposition is exhaustive at level n.
3. **Per-k geometric ratio:** for each fixed k, compute r_k(n) := c_{n,k} / c_{n−1,k} as n grows (telling us how the level-k contribution scales as we lift further).
4. **Cross-n ratio at fixed offset:** ρ(n) := ε_n / ε_{n−1} for n = 3, 4, 5, 6 — already known to approach ≈ 0.522 empirically (per R79b). Reproduce as anchor.
5. **Per-k contribution magnitude:** |c_{n,k}| as a function of k at fixed n. Locate the dominant k = k*(n) — the level whose W_k subspace contributes most to ε_n.

### Compute cost

Each π_n at n ≤ 6 takes ~minutes to hours over Q (n=6 alone was 446 s per R77.6). All R_k for k=1..5 are already computed and stored exactly (R77.5 anchors). φ_n at n ≤ 6 already exists in R76 / R77.x.

The new compute is: assembling lift_n(R_k) for each pair (n, k) and computing ⟨φ_n, lift_n(R_k)⟩. Each is O(3^n) rational operations. Total tabulation ≤ 6 × 5 × O(3^6) = ~22,000 rational dot products. Wall-clock estimate: **~1 hour for the full (n=2..6, k=0..n−1) table**, dominated by the n=6 inner products.

Estimated session compute budget: **~2-3 hours total**, including verification and CSV output. Well within one focused session.

---

## What's verified (signatures consistent with Candidate-A framework)

The Candidate-A framework predicts: rate-1/2 of ε_n is carried by *the projection of φ_n onto Σ_k W_k*, NOT by the per-W_k norms (which contract at rate 1/√3 per R77.5 Stage 1).

Specifically, **at least one of the following two patterns** must appear if Candidate A is the right entry point:

**Pattern A1 — Per-k geometric decay.** For some specific k* (possibly k*(n) growing with n, possibly fixed), the dominant contribution c_{n,k*} satisfies

```
|c_{n,k*}|  ~  C · 2^{−n}    as n grows
```

while subdominant contributions decay strictly faster. This would localize rate-1/2 to a specific W_k subspace.

**Pattern A2 — Cross-k phase-cancellation envelope.** The per-k contributions individually decay at rates ≠ 1/2, but their *signed sum* exhibits cancellation that produces rate-1/2 in ε_n. Specifically:

```
Σ_k c_{n,k} = ε_n ~ C · 2^{−n}    while    max_k |c_{n,k}| / 2^{−n} → ∞ as n grows.
```

This would correspond to rate-1/2 being a *spectral-density* feature of Φ_∞ rather than an isolated eigenvalue — consistent with R77.6's branch-cut finding at z=2.

### Confirming signatures

The Candidate-A framework is **confirmed as the right entry point** if:

(C1) The exact-Q decomposition `Σ_{k=0}^{n−1} c_{n,k} = ε_n` (up to W_0 boundary handling) holds at n=2..6. This validates that W_k filtration is a complete carrier of ε_n's structure.

(C2) Either Pattern A1 or Pattern A2 appears, with the rate ratio matching ≈ 0.522 (empirical) or 1/2 (theory) within the level-budget noise.

(C3) The per-k contributions can be interpreted as inner products on a *single* Hilbert space L²(Ẑ_3^×, μ) via the natural lift-isometry — i.e., no rectangular-operator pathology of the kind that killed the R_K probe.

---

## What's falsified (signatures inconsistent with Candidate-A framework)

The Candidate-A framework is **falsified as a useful entry point** if:

(F1) The decomposition `Σ_{k=0}^{n−1} c_{n,k} = ε_n` fails at exact-Q level. (Should not happen if the W_k construction is correct — but a strong sanity check.)

(F2) **All** per-k contributions decay at rate 1/√3 ≈ 0.577 (= the L² norm contraction rate from R77.5 Stage 1), reflecting only the trivial cardinality scaling. No Pattern A1, no Pattern A2. ε_n's rate-1/2 in this case would NOT be encoded in the W_k filtration structure — it would be living somewhere outside the L²(Ẑ_3^×) framing entirely.

(F3) The dominant k*(n) wanders erratically with n, with no consistent assignment of rate-1/2 to any subspace family. Suggests the W_k decomposition is "the wrong basis" — the right framework might still be L²(Ẑ_3^×) but in a different decomposition (which would push us toward Candidate B's wavelet basis as the next probe).

### Falsifying signature

If (F2) or (F3) hold, Candidate A is *the right Hilbert space but the wrong basis*. The probe's recommendation would be to next test Candidate B (p-adic wavelets) on the same L²(Ẑ_3^×), with the Kozyrev basis as the alternative decomposition.

If (F1) holds, the W_k construction itself is broken — would invalidate R77.5 Stage 2 (c_k = 0 over Q). Extremely unlikely; (F1) is a sanity check.

---

## Test specification (pre-registration-ready)

### Hypotheses

- **H_A_CONFIRMED:** (C1) ∧ ((C2 pattern A1) ∨ (C2 pattern A2)) ∧ (C3). Candidate-A framework is the right entry point; full Reading A construction should proceed.
- **H_A_FALSIFIED_WRONG_BASIS:** (C1) holds but (F2) or (F3). L²(Ẑ_3^×) is right, W_k filtration is wrong basis; next probe is Candidate B.
- **H_A_FALSIFIED_WRONG_HILBERT_SPACE:** (F1) — unexpected, would require deeper structural review.
- **H_A_INCONCLUSIVE:** levels n=2..6 give a mixed signal — e.g., partial pattern A1 at k ≤ 2 but n-budget too short to discriminate from cardinality scaling. Recommendation in this case: extend to n=7 (substantial: ~hours of Q linear algebra) before deciding.

### Pre-registration favored outcome

`H_A_CONFIRMED` or `H_A_INCONCLUSIVE`. The R77.6 generating-function evidence (branch cut at z=2, not simple pole) lightly favors Pattern A2 (phase-cancellation envelope) over Pattern A1 (isolated eigenvalue). Single-W_k localization of rate-1/2 would be inconsistent with branch-cut singularity — that finding makes Pattern A2 the expected confirming pattern.

`H_A_FALSIFIED_WRONG_BASIS` is a coherent possibility — would route the project to Candidate B.

`H_A_FALSIFIED_WRONG_HILBERT_SPACE` is essentially ruled out by R77.5's structural identities; (F1) would be a bug.

### Procedure

1. **Phase A:** Extend `result_77_5_compute_R_k.py` to compute and store `lift_n(R_k)` for all (n, k) with 1 ≤ k+1 ≤ n ≤ 6, in fractions.Fraction format.
2. **Phase B:** Load existing φ_n at n=2..6 from R76 / R77.x storage.
3. **Phase C:** Compute c_{n,k} = ⟨φ_n, lift_n(R_k)⟩ over Q for all (n, k).
4. **Phase D:** Verify (C1) — decomposition sum equals ε_n at exact-Q. If fails: (F1) triggered.
5. **Phase E:** Tabulate per-k ratios r_k(n) and per-n dominant k*(n). Classify outcome per hypotheses above.

### Adversarial checks for the test itself

- **Don't curve-fit.** Decide rate by ratio of consecutive ε_n / ε_{n−1}, not by least-squares slope. R77.4-trap precedent.
- **Don't conflate confirmation patterns.** Pattern A1 (single-k localization) and Pattern A2 (cross-k cancellation) are different framework directions. Report which one (if either) the data supports.
- **Don't reify k*(n) without n-stability.** If k*(n) varies erratically with n, that's evidence of (F3), not "we found the rate-carrying scale."
- **Honest INCONCLUSIVE call.** n=2..6 is a 5-level budget; if patterns require n=7+ to discriminate, say so. Don't force a confirmation under a sparse signal.

---

## Compute cost and session feasibility

| Phase | Estimate | Notes |
|---|---|---|
| A: lift_n(R_k) for all valid (n,k), n ≤ 6 | ~20 min | All Q, extension of existing R77.5 code |
| B: load φ_n at n=2..6 | ~5 min | Already exists in R76 |
| C: c_{n,k} computation | ~30-60 min | Dominated by n=6 inner products over Q |
| D: decomposition sanity check | ~5 min | Trivial Q sum |
| E: pattern classification, write-up | ~30 min | Manual / inspection |

**Total: ~2 hours of compute + ~1 hour of analysis = one focused session.** Well within the scoping-probe budget.

Note: this is ONE TEST. It is NOT the full Reading A construction. A confirming outcome (H_A_CONFIRMED) tells us *which framework to spend the multi-week construction on*; it does NOT close c=7/45. The full Reading A construction (deriving Φ_∞ on L²(Ẑ_3^×), proving its spectrum has rate-1/2 in the appropriate sense, feeding into a Nisoli-style closure or its analog) is a separate project per R_K_DISPOSITION recommendation (a).

---

## What this test does NOT do

- Does not derive Φ_∞ or its spectrum.
- Does not establish that rate-1/2 is a spectral feature of any operator on L²(Ẑ_3^×) — only tests whether the W_k filtration is *consistent* with carrying that rate.
- Does not provide a rigorous bound; just empirical pattern.
- Does not produce a publishable result on its own. It produces a *decision* — pursue full Candidate-A construction (multi-week), or pivot to Candidate B as the next scoping probe.

This is intentional. The probe is "where to start," not "construct it."
