# M3_DISPOSITION — top-level disposition of M_3 explicit-constants probe

**Date:** 2026-05-11. Wilson (analyst) reporting to Nathan. Top-level disposition of the M_3 extraction probe for c = 7/45 Nisoli closure pipeline.

---

## DISPOSITION: **H_M3_INTRACTABLE**

> The M_3 numerical bound is computable as a fact about the literal R77.2 matrix T_3 (`M_3 ∈ [50, 200]` anticipated numerical, `≤ 944` rigorous upper bound from κ(V) × 8, vs R77.2's loose 800–1000). **But this M_3 doesn't pertain to any operator describing ε_n's actual dynamics**, because R77.3 falsified the 3-mode {1/2, 1/4, 1/8} spectrum that T_3 was built to encode.

**Specific obstruction:** the load-bearing spectrum that M_3's contour `γ = |z − 1/2| = 1/8` was designed to isolate **does not exist on any project-characterized operator**. R77.3 falsifies the candidate T_3 (3-mode model); R77.4 erratum rules out the natural within-level Markov K_k (no eigenvalue near 1/2 at any k ∈ {3..6}); R77.4 erratum §"What this DOES change" lists alternative candidate operators (inter-level residual, generating-function singularity) as **parked pending direction** — none currently characterized.

The probe's pre-registered favorite was H_M3_PARTIAL (BORDERLINE or NEEDS_LARGER_K). Reality is harder: M_3 is computable for the *wrong* T, and there is no *right* T currently in the project.

---

## What was done

| Phase | Output | Verdict |
|-------|--------|---------|
| 1: Articulate | `M3_DEFINITION.md` | M_3 = sup_γ ‖R(z, T_3)‖_op. T_3 = R77.2 companion matrix, spec **falsified** by R77.3 as descriptor of ε_n. |
| 2A: Spectral | `M3_APPROACH_A.md` | Spectral-radius lower bound = 8; κ(V) × 8 upper bound = 944. Bracket: M_3 ∈ [8, 944]. T_3 non-normal, so spectral radius is strict lower bound. |
| 2B: Perturbation | `M3_APPROACH_B.md` | Diverges. Both natural splits (diagonal, shift-block) give `‖T_1‖ · ‖R(z, T_0)‖ ≫ 1`, no refinement of A. |
| 2C: Numerical | `M3_APPROACH_C.md` | Specified (python denied this task). Anticipated `M_3 ∈ [50, 200]`. k-extension to higher modulus also falsified (R77.3 §7). |
| 3: Closure table | `M3_CLOSURE_TABLE.md` | Parameterized over (K, A) and (M_3 = 100, 1000); closure feasible at A ≥ 3 (M=100) or A ≥ 5 (M=1000) for K ≥ 6, but **conditional** on Tao C_A = O(1) and conditional on T_3 actually describing ε_n. |

---

## M_3 explicit value / bound

**For the literal R77.2 T_3** (3×3 companion of `(7/8, −7/32, 1/64)`, spectrum {1/2, 1/4, 1/8} as a fact about the matrix):

- **Rigorous lower bound (spectral radius on γ):** M_3 ≥ 8.
- **Rigorous upper bound (Approach A, κ(V) method):** M_3 ≤ 944.
- **Best numerical estimate (Approach C anticipated):** M_3 ∈ [50, 200].
- **R77.2 §3.3 quoted:** "800–1000" (loose) up to "11320" (crude).

Derivation chain (Approach A):

```
‖R(z, T_3)‖_op  ≤  κ(V) · max_λ 1/|z−λ|
                ≤  (‖V‖₂ · ‖V⁻¹‖₂) · 8        on γ
                ≤  1.843 · 64 · 8
                ≈  944.
```

Lower bound: spectral radius equals `max_λ 1/|z−λ| = 8` on γ. T_3 non-normal so this is strict lower bound.

**Caveat (load-bearing):** the matrix this M_3 describes was conjectured to encode ε_n's recursion. R77.3 §3–§5 falsifies the recursion at n = 1 (predicted ε_4 = −222733/65244480 ≈ −3.41×10⁻³; actual ≈ −2.45×10⁻³, residual ~28%). The matrix is mathematically real; the operator-theoretic claim about ε_n is not.

---

## Parameterized closure table (abbreviated)

`|K| · K^{−A} · M_3 < 1`. With `|K|/√q = 2` (from polylog-free 2√p·√N at r ≥ 4):

| K | A | Closes at M_3=100? | Closes at M_3=1000? |
|---|---|---------------------|----------------------|
| 6 | 2 | NO (5.56) | NO (55.6) |
| 6 | 3 | YES (0.93) | NO (9.3) |
| 6 | 5 | YES (0.026) | YES (0.26) |
| 10 | 3 | YES (0.20) | NO (2.0) |
| 10 | 5 | YES (0.002) | YES (0.020) |
| 15 | 2 | YES (0.89) | NO (8.9) |
| 15 | 3 | YES (0.06) | YES (0.59) |
| 20 | 2 | YES (0.50) | NO (5.0) |
| 30 | 2 | YES (0.22) | NO (2.2) |
| 30 | 3 | YES (0.007) | YES (0.074) |

**Reading:** at M_3 = 100 and Tao's plausible A = 3, closure holds **everywhere from K = 6 up**. At M_3 = 1000 and A = 2, closure **fails everywhere**.

**Tao A delivery:** Tao Prop 1.17's effective C_A is **INFEASIBLE this session** per R77.2 §3.4 and the user task statement. The parameterization is honest "if Tao gives us A = X" framing; the table does **not** claim closure.

---

## Why H_M3_INTRACTABLE (not BORDERLINE)

The honest read is that **the operator T behind M_3 is undefined**:

1. R77.2's T_3 (conjectured via 3-mode model) is **falsified** by R77.3 — recursion fails at n=1 with 28% relative residual in Q.
2. The natural within-level transition K_k has **no eigenvalue near 1/2** at any k ∈ {3..6} (R77.4 erratum §"Empirical evidence").
3. The 4-mode extension is **also falsified** (R77.3 §7, 5–10% residuals at n=5,6).
4. Alternative candidate operators (inter-level residual R_k, generating-function singularity) are **parked pending direction** (R77.4 erratum §recommended next moves).

So while one can compute M_3 ≈ 100 for the falsified T_3, this number doesn't feed any actually-correct Nisoli application. **It is operationally moot.**

This is qualitatively different from BORDERLINE (where M_3 exists, has factor-of-K uncertainty, and tightening would help) or NEEDS_LARGER_K (where M_3 exists, closure works at larger K). In both BORDERLINE and NEEDS_LARGER_K, the operator-theoretic anchor exists. Here it does not.

---

## Recommendation

**Pivot: accept M_3 as an additional open question on the c=7/45 closure roadmap.**

The c = 7/45 Nisoli closure path has at least **three** independent open obstructions:

1. **Tao Prop 1.17 effective C_A** (R77.2 §3.4 INFEASIBLE this session; standalone re-derivation of Tao §7.2–7.3 with bookkeeping).
2. **|K| bilinear bound** (HENSEL_DISPOSITION polylog-free 2√p·√N delivered; TIGHTEN_* may upgrade to strict 2√N).
3. **M_3 = ‖R(z, T)‖_op for a characterized T** (this probe: T currently unspecified after R77.3 falsification).

Item 1 (Tao C_A) and item 3 (M_3 anchor) are both pre-conditions for Nisoli closure even with item 2 fully delivered. Item 3 is the more structural obstruction: even effective Tao C_A doesn't help if the operator the C_A applies to isn't characterized.

**What this means for the existing toolkit:**

- The bilinear-bound work (PATH2, HENSEL, in-flight TIGHTEN_*) continues to be useful as a **standalone √N character-sum result** — publishable in its own right per PRECISE_ASK §4 "most useful target for the literature."
- The Nisoli pipeline as a whole is **not currently closable via R77.2's framework**; it would require either (a) a new candidate operator T whose spectrum is characterized (R77.4 erratum §1: inter-level residual operator is the natural next probe), or (b) abandoning the spectral-gap framing in favor of a different rigor route (e.g., direct envelope bound from PRECISE_ASK §4 (c)).

**Recommended sequencing:**

1. **Accept H_M3_INTRACTABLE for now.** Don't burn cycles re-deriving M_3 against more candidate T's until one is operator-theoretically characterized.
2. **Continue TIGHTEN_* and CHAIN_RECON_* in parallel.** Bilinear-bound work has standalone value regardless of Nisoli status.
3. **If pursuing Nisoli closure remains the goal:** the next natural probe is R77.4 erratum §1 (inter-level residual operator R_k spectrum). That's a discrete, finite-compute task; if R_k has a clean eigenvalue near 1/2 with explicit eigenvector basis, the M_3 calculation regenerates with operator-theoretic backing.
4. **If c=7/45 rigor can be achieved by other means** (e.g., PRECISE_ASK §4 (c) phase-cancellation tighter than C_emp · √N — but the empirical β = 0.522 already includes all empirical phase cancellation, so this seems unlikely), Nisoli is unnecessary anyway.

---

## Disposition file references

- `M3_DEFINITION.md` — Phase 1 (M_3 articulated, R77.3 falsification documented)
- `M3_APPROACH_A.md` — Phase 2A (spectral; M_3 ∈ [8, 944])
- `M3_APPROACH_B.md` — Phase 2B (perturbation; diverges, no refinement)
- `M3_APPROACH_C.md` — Phase 2C (numerical; specified, M_3 ∈ [50, 200] anticipated)
- `M3_CLOSURE_TABLE.md` — Phase 3 (parameterized table, with full caveats)
- `M3_DISPOSITION.md` — this file

## Pre-registration mapping

Pre-registered hypotheses → outcome:

- H_M3_EXTRACTABLE_AND_CLOSES — REJECTED (M_3 extractable for wrong T; closure inequality at currently-feasible K requires A=3 with M=100 or A=5 with M=1000; conditional on Tao C_A = O(1) which is INFEASIBLE)
- H_M3_EXTRACTABLE_NEEDS_LARGER_K — REJECTED (closure works at K=6 already if Tao A=5; the obstruction is not K size, it's T characterization)
- H_M3_BORDERLINE — REJECTED (M_3 has explicit estimate but operator-theoretic anchor is missing)
- **H_M3_INTRACTABLE — ACCEPTED** (R77.3 falsified the spectrum; no alternative T characterized; load-bearing obstruction is the operator definition, not the resolvent calculation)
- INCONCLUSIVE — REJECTED (Phase 1 produced definite findings)

The pre-registration favored H_M3_PARTIAL. Actual outcome is harder than that — the issue isn't M_3 magnitude, it's that the M_3 calculation has nothing to attach to.
