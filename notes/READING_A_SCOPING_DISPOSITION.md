# READING_A_SCOPING_DISPOSITION — top-level disposition of R77.5 §7 scoping probe

**Date:** 2026-05-11. Wilson (analyst) reporting to Nathan. Top-level disposition of the Reading A function-space scoping probe for c = 7/45 spectral-completion question.

---

## DISPOSITION: **H_CANDIDATE_A_TRACTABLE**

> **Candidate A — Hilbert spaces of locally constant functions on Ẑ_3^× — is the most-viable entry point for the Reading A construction.** Construction is essentially implicit in R77.5's existing V_k / W_k apparatus; explicit basis at small k is a direct extension of `result_77_5_compute_R_k.py`'s exact-Q tooling; φ_n projection onto Σ_k W_k is a finite linear-algebra problem at each (n, k); external reading required for the full construction is *light, standard, citable* (Folland Ch. 3-6, Tate §2, Vladimirov-Volovich-Zelenov Ch. 6).
>
> A pre-registration-ready minimum-viable test is articulated in READING_A_SCOPING_MIN_VIABLE_TEST.md. **Estimated session cost: ~2-3 hours of focused computation.** Outcome of the test gates whether to invest 3-4 weeks in the full Candidate-A construction, or pivot to Candidate B (p-adic wavelets) as the next probe.

Pre-registered favored hypotheses were `H_CANDIDATE_A_TRACTABLE` or `H_NONE_TRACTABLE_AT_THIS_SCOPE`. Outcome matches the more optimistic of the two — Candidate A is tractable at the scoping level.

---

## Three-obstruction landscape — where this probe sits

| Obstruction | Status | Reference |
|---|---|---|
| (1) Tao Prop 1.17 effective C_A — Route 1 | INFEASIBLE | BOOKKEEPING_PHASE1_DISPOSITION |
| (2) Bilinear \|K\| | DELIVERED (strict 2√N at r≤3, polylog-free 2√p·√N at r≥4) | PATH2_DISPOSITION + HENSEL_DISPOSITION |
| (3) Spectral M_3 — within-level T_3 | FALSIFIED (R77.3) | M3_DISPOSITION |
| (3) Spectral M_3 — inter-level R_k (rectangular operator) | INTRACTABLE (R77.5 / R_K probe) | R_K_DISPOSITION |
| **(3) Spectral M_3 — function-space Φ_∞ on L²(Ẑ_3^×)** | **scoping probe positive — Candidate A entry point identified, full construction is a 3-4 week subproject** | **this disposition** |

This is the third candidate-operator probe after **M_3** (literal T_3) and **R_K** (inter-level residual). The pattern across the three probes is consistent: the framework's spectral content lives at a level of mathematical abstraction the project hasn't yet developed. The Reading A scoping probe identifies the first concrete entry point for that development — the L²(Ẑ_3^×) function-space framework that R77.5 §7 named but didn't construct.

---

## Tractability assessment summary (Phase 1 outputs)

| Criterion | (A) Locally constant fns | (B) p-adic wavelets | (C) Transfer operator |
|---|---|---|---|
| Construction tractability | HIGH (essentially done in R77.5) | MEDIUM (Kozyrev 2002 + Khrennikov 2009) | LOW (no fixed Syracuse-Ẑ_3 extension) |
| Basis tractability k=2..5 | HIGH (Q linalg, exact rationals) | MEDIUM (roots-of-unity arithmetic) | N/A |
| φ_n articulation | MEDIUM (finite linalg) | LOW-MEDIUM (advantage conditional) | POTENTIALLY HIGH (gated by Syracuse-Ẑ_3 extension) |
| Min-viable test | YES (one session) | YES, heavier | NOT YET ARTICULABLE |
| External reading | LOW (standard, ~3 named books / papers) | MEDIUM (2 specific named papers + standards) | HIGH and scattered |
| R77.5 §7 fidelity | direct | direct | direct (but specification gap) |
| Single Hilbert space framing | ✓ | ✓ | ✓ if defined |
| Falsifiable | ✓ Pattern A1 vs A2 vs (F2 / F3 / F1) | ✓ similar shape | not yet |

**Verdict:** Candidate A wins on every assessed dimension. Candidate B is the natural fallback if Candidate A's test falsifies the W_k basis (signatures F2 / F3). Candidate C is conceptually clean but blocked at specification level.

---

## Minimum-viable test specification (Phase 2)

**For Candidate A (the top-ranked entry point):**

### What's computed

For each (n, k) with 0 ≤ k < n and n ∈ {2, 3, 4, 5, 6}, compute

```
c_{n,k} := ⟨φ_n, lift_n(R_k)⟩  ∈ Q (exact rational).
```

All quantities live in existing R77.5 / R76 machinery — `pi_dict`, `lift_pi`, R_k = π_{k+1} − T(π_k), φ_n from R76's bilinear pair-form moment definition. Compute is `fractions.Fraction` linear algebra at level n, dominated by n=6 (~hour).

### What's verified (confirms Candidate A)

- **(C1) Decomposition sanity:** `Σ_{k} c_{n,k} = ε_n` over Q. (Sanity check on the W_k filtration's completeness.)
- **(C2) Rate-1/2 signature** in one of two forms:
  - **Pattern A1:** dominant single-k contribution decays at rate ~0.5 per level; others decay strictly faster.
  - **Pattern A2:** per-k contributions decay slower than 0.5 (likely at the trivial 1/√3 cardinality rate), but their *signed sum* exhibits phase-cancellation that produces rate-0.5 in ε_n.
- **(C3) Single-Hilbert-space framing preserved** — no rectangular-operator pathology.

### What's falsified (rules out Candidate A or reroutes)

- **(F1) Sum ≠ ε_n** — would indicate R77.5 structural identity has a bug (very unlikely).
- **(F2) All c_{n,k} decay at trivial 1/√3 rate, no cross-k cancellation** — would mean L²(Ẑ_3^×) framing is right but W_k basis is wrong, route to Candidate B.
- **(F3) Erratic k*(n)** — same conclusion as (F2), route to Candidate B.

### Compute cost

**~2-3 hours session compute**, all over Q, no external packages beyond `fractions`. Well within scoping-probe budget.

---

## Full-construction project scope (Phase 3, Rank #1 — Candidate A)

If H_A_CONFIRMED from the minimum-viable test, the full Reading A construction is:

- **Phase I — Standard-apparatus port (2-3 days):** establish L²(Ẑ_3^×, μ) with 3-adic Haar measure, identify R77.5's V_k / W_k as the standard scale-filtration of L²(profinite abelian group).
- **Phase II — Φ_∞ definition (3-5 days):** define the projective-limit operator on L²(Ẑ_3^×) that restricts to K_k on each V_k; verify bounded-operator extension.
- **Phase III — Spectral characterization (1-2 weeks):** compute spectrum at small-k truncations, identify rate-1/2 as (a) isolated eigenvalue (unlikely per R77.6 branch-cut evidence), (b) branch-cut feature, or (c) spectral density.
- **Phase IV — Closure attempt (1-2 weeks):** either Nisoli-style resolvent bound on Φ_∞, or direct Tauberian extraction from generating function (R77.6 branch-cut framing).

**Total: ~3-4 weeks of focused construction.** This is a separate research subproject, not a session task. The scoping probe identifies it; future work executes it.

### External reading required for full construction

- **Folland G.B., _A Course in Abstract Harmonic Analysis_, 2nd ed. (2016)**, Ch. 3-6 — Haar measure, profinite groups, Plancherel.
- **Tate J. (1950 thesis, in Cassels-Frohlich), §2** — local L²-theory on Ẑ_p^×.
- **Vladimirov V.S., Volovich I.V., Zelenov E.I., _P-adic Analysis and Mathematical Physics_** (World Scientific 1994), Ch. 6 — integral transforms on Ẑ_p.

Standard, named, citable. Reading scope ~2-3 days before construction begins.

---

## Adversarial check outcomes

**(A1) R77.5 §7 fidelity.** Three candidates assessed = three R77.5 §7 candidates. No fourth proposed. ✓

**(A2) Avoid the operator-on-fixed-Hilbert-space-versus-rectangular trap.** All three candidates respect L²(Ẑ_3^×) as a single Hilbert space (Φ_∞ on the same space, not rectangular maps W_{k−1} → W_k). The R_K probe's failure mode is structurally avoided. ✓

**(A3) External-machinery honesty.** Candidate A external reading is light and named (Folland, Tate, Vladimirov-Volovich-Zelenov). Candidate B names Kozyrev (2002) and Khrennikov-Shelkovich-Skopina (2009) specifically. Candidate C honestly flagged as having a specification gap — the "Syracuse coherent extension to Ẑ_3" is not in named literature in the form needed. No "p-adic wavelets are well-known" handwaving. ✓

**(A4) Min-viable-test falsifiability.** Candidate A's test distinguishes Pattern A1 vs Pattern A2 vs (F1, F2, F3). Each has a distinct empirical signature in c_{n,k}'s scaling. No "compatible with anything" pathology. ✓

---

## Recommendation

**Run Candidate A's minimum-viable test in a single focused session.** It's ~2-3 hours of `fractions`-arithmetic compute. The test gates whether to invest 3-4 weeks in the full Reading A construction.

If the test **confirms** Candidate A (Pattern A1 or A2), the spectral M_3 obstruction (#3 in the three-obstruction landscape) becomes a well-defined construction project rather than an open structural gap. This would unlock a *potential* path to c=7/45 closure (combined with delivered bilinear bound from PATH2 / HENSEL), or alternatively a standalone-publishable result on the spectral theory of profinite Markov chains.

If the test **falsifies** Candidate A in the F2 / F3 sense (W_k basis is wrong), Candidate B (p-adic wavelets) becomes the next scoping probe. The L²(Ẑ_3^×) framing survives; only the basis choice changes.

If the test is **inconclusive** at n=2..6 budget, extension to n=7 (~hours of Q linear algebra) is the natural next step before deciding.

**Do not commit to the full 3-4 week construction without first running the minimum-viable test.**

---

## Pre-registration mapping

Pre-registered hypotheses → outcome:

- **H_CANDIDATE_A_TRACTABLE — ACCEPTED.** Construction implicit in R77.5; basis tractable at small k via existing Q machinery; external reading light and standard.
- **H_CANDIDATE_B_TRACTABLE — REJECTED (as entry point).** Becomes preferred only if Candidate A's test falsifies the W_k basis.
- **H_CANDIDATE_C_TRACTABLE — REJECTED.** Construction step (a) — choosing the Syracuse coherent extension to Ẑ_3 — is unresolved in R77.5 and in the project. Blocks scoping articulation.
- **H_NONE_TRACTABLE_AT_THIS_SCOPE — REJECTED.** Candidate A's minimum-viable test is articulable and project-internal at the scoping level.
- **INCONCLUSIVE — REJECTED.** R77.5 §7 / §10 articulation suffices to assess the three candidates; the probe produced definite ranking.

The pre-registration favored `H_CANDIDATE_A_TRACTABLE` or `H_NONE_TRACTABLE_AT_THIS_SCOPE`. Outcome matches the more optimistic case.

---

## Disposition file references

- `READING_A_SCOPING_CANDIDATES.md` — Phase 1, candidate-by-candidate tractability assessment.
- `READING_A_SCOPING_MIN_VIABLE_TEST.md` — Phase 2, pre-registration-ready test specification for Candidate A.
- `READING_A_SCOPING_RECOMMENDATION.md` — Phase 3, ranked ordering and rationale.
- `READING_A_SCOPING_DISPOSITION.md` — this file.

### Anchors

- `result_77_5_inter_level_residual.md` — R77.5 §7 source for the three candidates; §10 names the open construction steps.
- `result_77_5_compute_R_k.py` — Q-arithmetic tooling that the minimum-viable test extends.
- `result_77_6_generating_function.md` — branch-cut evidence at z=2 (informs Pattern A2 expectation; rules out Pattern A1 if R77.6 holds at higher orders).
- `M3_DISPOSITION.md`, `R_K_DISPOSITION.md` — precedent probes establishing the structural shape of the spectral M_3 obstruction.

---

## What this probe does and doesn't do

**Does:**
- Assesses tractability of R77.5 §7's three candidates as entry points.
- Identifies Candidate A as most tractable, with specific rationale per criterion.
- Articulates a falsifiable minimum-viable test for Candidate A.
- Names external reading for the full Reading A construction.
- Estimates full-construction scope (3-4 weeks) honestly.

**Does NOT:**
- Attempt the full Reading A construction.
- Construct any function space, basis, or operator.
- Claim that Candidate A *will* deliver rate-1/2 in the W_k filtration — that's the test's job to answer.
- Propose any fourth framework not in R77.5 §7.
- Read external papers; only names them.

The probe is "where to start," not "construct it."
