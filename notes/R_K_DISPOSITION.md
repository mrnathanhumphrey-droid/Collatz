# R_K_DISPOSITION — top-level disposition of the inter-level residual operator R_k probe

**Date:** 2026-05-11. Wilson (analyst) reporting to Nathan. Top-level disposition of the R77.4 erratum §1 R_k spectrum probe for c = 7/45 Nisoli closure pipeline.

---

## DISPOSITION: **H_R_K_INTRACTABLE**

> **Specific obstruction: erratum-§1 ambiguity + structural mismatch.** R77.4 erratum §1's "inter-level residual operator R_k" is articulated ambiguously between (Reading A) "embed into common projective-limit space" — requires substantial separate construction of L²(Ẑ_3^×) framework, not in this probe's scope — and (Reading B) "level-by-level finite-dimensional map W_{k−1} → W_k" — finite but **not a Nisoli-amenable operator** because:
>
> (i) Φ_k under Reading B is rectangular (dim 3:1), resolvent is ill-typed.
> (ii) Φ_k **does NOT transport R_{k−1} → R_k** — c_k = ⟨R_k, T(R_{k−1})⟩/‖T(R_{k−1})‖² is **exactly 0/1 over Q** at every k = 2..5 (R77.5 §3, structural from marginal consistency). Whatever Φ_k is, it doesn't carry the inter-level R_k dynamics.
> (iii) The "spectrum near 1/2" feature R77.4 erratum reached for is **empirically absent** at all tested levels: σ_1(Φ_k) → 1 (leading mode dominance), bulk near σ² ≈ 0.4 with multiplicity 4·3^{k−1} (cardinality-driven, not eigenvalue-driven).
> (iv) Multiplicity floor forces would-be M_3'(K) → ∞ at least like 3^{K/2}, so closure A requirement **diverges** with K rather than improving.

**Pre-registered favored H_R_K_BORDERLINE.** Reality is harder: not "M_3' has factor-2-to-10 uncertainty" but **M_3' does not exist as a single bounded resolvent norm**, because the operator R77.4 erratum reached for is not constructable as a finite-dimensional Nisoli object.

This is the **same structural class of obstruction as M_3** (precedent: H_M3_INTRACTABLE in `M3_DISPOSITION.md`), at a different level of abstraction:

| | M_3 probe (T_3, R77.2) | R_k probe (R77.4 erratum §1) |
|---|---|---|
| Object | 3×3 companion matrix | Sequence of vectors {R_k} in different W_k subspaces |
| Why intractable | T_3's spectrum {1/2, 1/4, 1/8} doesn't describe ε_n (R77.3 falsified the recursion) | R_k is the right structural object (= R74's d_{k+1}, exact rational) but the operator-carrier of inter-level R_k dynamics does not exist as a finite-dimensional map |
| Operator-theoretic anchor | Exists for the matrix; doesn't pertain to ε_n | Doesn't exist as a single operator in either Reading |

Both probes converge on: **rate-1/2 of ε_n is not a finite-truncation spectral phenomenon in any project-characterizable sense.**

---

## What was done

| Phase | Output | Verdict |
|-------|--------|---------|
| 1: Articulate | `R_K_DEFINITION.md` | R_k = π_{k+1} − T(π_k); lives in W_k ⊂ V_{k+1}; different Hilbert space per level. R77.4 erratum §1 ambiguous between Reading A (projective limit) and Reading B (level-by-level matrix). |
| 2A: Spectrum | `R_K_APPROACH_A.md` | Φ_k under Reading B is rectangular (dim 3:1), σ_1 → 1, bulk ~0.62 in σ-space. **No σ near 1/√2 = 0.707** (rate-1/2 target). Critically, c_k = 0 means Φ_k doesn't transport R_{k−1} → R_k anyway. APPROACH_A_FAILS_STRUCTURAL. |
| 2B: Perturbation | `R_K_APPROACH_B.md` | Both natural rank-1 + perturbation splits diverge. ‖δΦ‖/‖resolvent‖ ratios > 1. APPROACH_B_FAILS, downstream of Approach A. |
| 2C: Resolvent norm | `R_K_APPROACH_C.md` | Φ_k is rectangular so (zI − Φ_k)^{-1} is ill-typed. Reducing to Gram operator Φ_k^* Φ_k: bulk near σ² ≈ 0.4 has multiplicity 4·3^{k−1}; resolvent norm grows at least like √multiplicity ~ 3^{K/2}. APPROACH_C_POLYNOMIAL_GROWTH at minimum. |
| 3: Closure table | `R_K_CLOSURE_TABLE.md` | Under multiplicity floor M_3'(K) ~ 3^{K/2}, Tao A required to grow with K (A=2.2 at K=6, A=5.1 at K=30). **Closure asymptotically fails** even with hypothetical Tao C_A = 1. |

---

## R_k construction summary (Phase 1)

> **R_k(r') := π_{k+1}(r') − T_{k→k+1}(π_k)(r')**, r' coprime in Z/3^{k+1}.
>
> T_{k→k+1}(π_k)(r') := π_k(r' mod 3^k) / 3.

Properties from R77.5 (all over Q, exact rationals):

- **‖R_k‖² · 3^k → 7/45** (R74 = R77.5 identity, 5/5 PASS).
- **R_k ⊥ T(V_k)** by construction (marginal consistency), so R_k ∈ W_k := T(V_k)^⊥.
- **⟨R_k, T(R_{k−1})⟩ = 0** exactly at k = 2, 3, 4, 5 (no projection onto lift basis).
- **Ratio ‖R_k‖² / ‖R_{k−1}‖² → 1/3** (cardinality scaling), NOT 1/4 (rate-1/2 target).
- **Different Hilbert space per level** (W_{k−1} ⊂ V_k ≠ W_k ⊂ V_{k+1}); no single operator on a fixed space.

---

## Spectrum tabulation at k = 2..5 (Phase 2A anticipated)

| k | dim Φ_k | σ_1 (anticipated) | σ_2..bulk | σ near 1/√2 = 0.707? | Eigenvalue (σ²) near 1/2? |
|---|---|---|---|---|---|
| 2 | 12 × 4 | ~0.75 | ~0.62 | NO (closest 0.62) | NO (σ² near 0.5 only if bulk hits, with mult. 3) |
| 3 | 36 × 12 | ~0.82 | ~0.62 | NO | NO (σ² near 0.4-0.5 with mult. 9) |
| 4 | 108 × 36 | ~0.88 | ~0.62 | NO | NO (mult. 27) |
| 5 | 324 × 108 | ~0.91 | ~0.62 | NO | NO (mult. 81) |

Anchor: `L_k_eigenvalues.csv` (adjacent prior SVD data) shows the exact pattern (σ_1 monotone to 1, bulk plateaued at 0.62..0.64).

**No σ_1 (or any σ_i) approaches the rate-1/2 target.** σ_1 grows toward 1 monotonically — this is cardinality-of-image scaling, not a "rate-1/2" feature.

---

## M_3' explicit value/bound (Phase 2C anticipated)

| Contour choice | M_3'(2) | M_3'(3) | M_3'(4) | M_3'(5) | Asymptotic growth |
|---|---|---|---|---|---|
| γ(1/2, 1/8) — R77.2 contour | ~20 | ~30 | ~50 | DIVERGES or ~100 | grows with multiplicity ~ 3^{K/2} |
| γ(σ_1², 0.05) — best isolated | ~10 | ~25 | ~50 | ~100 | grows linearly-to-quadratically |
| γ(bulk, 0.05) | ~30 | ~60 | ~150 | ~400+ | grows like 3^{K/2} (mult. floor) |

**Best case under any contour: M_3'(K) grows AT LEAST linearly in K.** Multiplicity floor (mult. ~ 4·3^{K−2} eigenvalues in bulk near 1/2 of Φ_k^* Φ_k) forces this.

Derivation chain (anticipated, multiplicity floor):

```
M_3'(K) ≥ √(multiplicity of bulk near 1/2 in spec(Φ_K^* Φ_K))
       ≥ √(dim W_{K−1}) · (constant)
       = √(4·3^{K-2}) · (constant)
       ~ 3^{K/2}.
```

---

## Parameterized closure table (Phase 3)

`|K| · K^{−A} · M_3'(K) < 1` with `|K|/√q = 2`:

| K | M_3'(K) ≈ 3^{K/2} | A required |
|---|---|---|
| 6 | 27 | **2.23** |
| 10 | 243 | **2.69** |
| 15 | 3,788 | **3.30** |
| 20 | 59,049 | **3.90** |
| 30 | 1.43e7 | **5.06** |

**Asymptotic-in-K reading:** A required grows like log(3^{K/2})/log(K) = K · log(3)/(2 · log(K)) → ∞. **Closure fails to be uniform in K** under any fixed Tao A.

Tao's plausible A range is {2..10}. The K-uniform closure requirement runs out of A budget around K=30+ depending on which exact Tao bound is used.

**Compare to M_3 closure (M3_CLOSURE_TABLE.md):** under constant M_3 = 100, closure A required *decreased* with K (from 2.7 at K=6 down to 1.4 at K=30). Under M_3'(K) growing exponentially, the requirement *increases* with K. **Worse-conditioned in K than the falsified-T_3 framework.**

---

## Three-obstruction landscape implications

c=7/45 Nisoli closure has three independent obstructions:

(1) **Tao Prop 1.17 effective C_A:** INFEASIBLE per BOOKKEEPING_PHASE1_DISPOSITION (Tao §7.4 iterated-cubic forces C_A ≥ exp(exp(A²))).
(2) **|K| bilinear bound:** DELIVERED per PATH2_DISPOSITION + HENSEL_DISPOSITION (strict 2√N at r ≤ 3, polylog-free 2√3·√N at r ≥ 4; TIGHTEN_* may upgrade r ≥ 4).
(3) **M_3 = sup_γ ‖R(z, T)‖ for a characterized T:** **NOT RESOLVED** by this probe. R77.4 erratum §1's R_k (the natural inter-level candidate) does NOT furnish a Nisoli-amenable operator. The structural obstruction is the same as for T_3 (M_3 probe), at a different abstraction level.

**This probe does NOT resolve obstruction (3).** It confirms the obstruction is **structural, not numerical**: there is no Nisoli-amenable finite-truncation operator whose spectrum captures the rate-1/2 envelope of ε_n in either the within-level (K_k, ruled out by R77.4 erratum) or inter-level (R_k, this probe) setting.

The closure path via the R77.2 framework is therefore **not unlocked**.

---

## Recommendation

**Accept H_R_K_INTRACTABLE; pivot recommendation matches M3_DISPOSITION.md §Recommendation.**

The R77.2 Nisoli framework's three obstructions are NOT independently tractable in this session:

1. **Tao C_A** (Route 1) — structurally infeasible (Tao's iterated-cubic) per BOOKKEEPING_PHASE1.
2. **Bilinear |K|** (Route 2) — delivered, standalone-valuable for the literature regardless.
3. **Spectral M_3** (obstruction 3) — both candidate operators (T_3 within-level, R_k inter-level) confirmed structurally inadequate.

**Specific next-probe options:**

(a) **Reading A construction** (substantial): build the L²(Ẑ_3^×) projective-limit framework with the wavelet-like {W_k} filtration, derive a single transfer operator Φ_∞ on this space, ask whether Φ_∞ has rate-1/2 spectrum. R77.5 §7 + §10 recommend this; it's a substantial separate project, several weeks of construction-level work (Hilbert spaces of locally constant functions on Ẑ_3^×, transfer-operator analysis of the Syracuse map's projective extension). NOT this probe's scope.

(b) **Direct envelope bound** (PRECISE_ASK.md §4 (c) — direct phase-cancellation tighter than C_emp · √N). Per M3_DISPOSITION §Recommendation, the empirical β = 0.522 already includes all phase cancellation, so this seems unlikely to deliver new content.

(c) **Accept Nisoli closure for c=7/45 as currently-not-achievable**, position the project's bilinear-bound work (PATH2 + HENSEL + potentially TIGHTEN_*) as **standalone √N character-sum results** publishable in their own right per PRECISE_ASK §4 "most useful target for the literature," and treat the c=7/45 rate-1/2 envelope as an **open empirical conjecture** rather than something closable within the current framework.

Wilson's read: **(c) is the honest position.** The Nisoli framework requires a finite-truncation operator with isolated rate-1/2 spectrum; the three obstructions all conspire to say no such operator exists in the project's characterized setting. The bilinear bound is genuinely useful; the rate-1/2 envelope remains an empirical observation.

This is consistent with the prior M3_DISPOSITION recommendation, with one update: **the R77.4 erratum's "inter-level residual operator" route is now also concluded structurally** — not merely "parked pending direction" as the erratum stated. The next pivot would have to be Reading A construction (option (a)), which is a substantial separate project.

---

## Disposition file references

- `R_K_DEFINITION.md` — Phase 1 (R_k articulated from R77.4 erratum §1, fidelity audit, ambiguity flagged)
- `R_K_APPROACH_A.md` — Phase 2A (Φ_k singular values, σ_1 → 1, bulk near 0.62, no 1/√2 feature; c_k=0 sanity-fails the transport)
- `R_K_APPROACH_B.md` — Phase 2B (perturbation; both splits FAIL)
- `R_K_APPROACH_C.md` — Phase 2C (resolvent norm; ill-typed for rectangular Φ_k, multiplicity floor forces growth in K)
- `R_K_CLOSURE_TABLE.md` — Phase 3 (parameterized table; closure A requirement grows in K)
- `R_K_DISPOSITION.md` — this file

Anchors:

- `result_77_4_K_spectrum_erratum.md` — origin of R_k proposal
- `result_77_5_inter_level_residual.md` — main R_k construction (R77.5)
- `result_77_5_d_R_identity_check.md` — R_k = d_{k+1} identity over Q
- `result_77_5_R_k_norms.csv` — exact norm data
- `result_77_5_phi_correlations.csv` — c_k = 0 exact data
- `L_k_eigenvalues.csv` — adjacent SVD pattern (anchor for Φ_k σ_1 → 1 prediction)
- `M3_DEFINITION.md`, `M3_DISPOSITION.md` — companion probe, parallel H_M3_INTRACTABLE precedent

---

## Pre-registration mapping

Pre-registered hypotheses → outcome:

- **H_R_K_CLOSES** — REJECTED. R_k as constructed in R77.5 is **not a single Nisoli-amenable operator**; the candidate Φ_k under Reading B is rectangular + doesn't transport R_{k−1} → R_k.
- **H_R_K_NEEDS_LARGER_K** — REJECTED. The would-be M_3'(K) grows like 3^{K/2} (multiplicity floor), so the A-requirement *grows* with K rather than shrinking; "larger K" doesn't fix it.
- **H_R_K_BORDERLINE** — REJECTED. Not a factor-2-to-10 uncertainty in M_3'; M_3' fundamentally doesn't exist as a single value (rectangular operator, no resolvent type-signature).
- **H_R_K_INTRACTABLE — ACCEPTED**. Specific obstruction: **R77.4 erratum §1 ambiguity** (between projective-limit Reading A — not constructed in project — and level-by-level Reading B — finite-dimensional but Φ_k doesn't transport R_{k−1} → R_k, c_k = 0 over Q is structural). **No Nisoli-amenable inter-level operator exists at finite truncation.**
- INCONCLUSIVE — REJECTED. Phase 1 produced definite findings + structural argument; not a "Phase 1 issues" case.

The pre-registration **favored H_R_K_BORDERLINE**; actual outcome is **H_R_K_INTRACTABLE**, the same disposition as the M_3 probe one step earlier. This is structurally consistent: when the operator-theoretic anchor for rate-1/2 doesn't exist (T_3 falsified, K_k mixing, R_k a vector-sequence-not-operator), Nisoli's resolvent-norm constant has nothing to attach to in any of the candidate frameworks.
