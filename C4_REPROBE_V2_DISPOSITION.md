# C4_REPROBE_V2_DISPOSITION

**Date:** 2026-05-14. Focused re-probe v2 of Component 4 (= Ch 11 renewal-Egorov composition formula) against:
- **Primary:** Cébron 2013 "Free Convolution Operators and Free Hall Transform" (arXiv:1304.1713v3, 55 pp)
- **Secondary:** Goldsheid-Margulis 1989 "Lyapunov indices of a product of random matrices" (Russian Math Surveys, 61 pp)

**Sources:** `C:/Users/Nate/OneDrive/Documents/closure hunt/{1304.1713v3.pdf, rm1893_eng.pdf}`. Verbatim extractions in `C:/Collatz/_cebron_pages/` (55 files) and `C:/Collatz/_goldsheid_pages/` (61 files).

**Background:** Voiculescu-Dykema-Nica 1992 "Free Random Variables" (CRM Monograph Vol. 1) is **scanned PDF** (no extractable text via pypdf); all references to VDN are via Cébron's citations — marked [NON_VERBATIM_VDN] in companion files.

Cross-refs (not re-extracted): `C4_REPROBE_TAO_RMT_DISPOSITION.md` (v1 finding), `PROFINITE_TRANSFER_OPERATOR_BLUEPRINT.md`, `C1_TAO_RECURSION_FORM.md`.

---

## Status: **C4_PARTIAL_OUTSIDE**

The corpus **delivers a sharper version of the v1 finding** but does **not close C4**. The needed next step — operator-valued / amalgamated free probability for non-free random variables — is **outside this corpus**. The Ch 11 effort estimate **slightly tightens** (lower bound shifts from "stays open" to "stays open with one specific paper now identified as the next probe target"), but the load-bearing chapter remains research-scope.

---

## Closest-fit theorem (the cleanest match in the corpus)

**Cébron 2013, Theorem 2.13 (§2.7, p. 24):**

> "Let I be an arbitrary index set. Let A = (A_i)_{i ∈ I} ∈ A^I be such that τ(A_i) ≠ 0 for all i ∈ I. For all P ∈ C{X_i : i ∈ I}, and all B = (B_i)_{i ∈ I} ∈ A^I **free from (A_i)_{i ∈ I}** and such that τ(B_i) ≠ 0 for all i ∈ I, we have
>
> τ(P(AB)|B) = (e^{D_A} P)(B)."

This is the EXACT multiplicative composition formula that Tao RMT §2.5 Remark 2.5.24 punted to Speicher. Cébron supplies it as a named theorem in a 55-page paper with full proof.

**Supplemental:** Cébron Theorem 4.6 (§4.4, p. 50) gives the **large-N (matrix size) convergence** of the classical Brownian motion on GL_N(C) to the free circular multiplicative BM, with rate O(1/N²) — driven by smoothness of exp on End(C^d{X_i}).

**Supplemental:** Goldsheid-Margulis Theorem 1.2 (§1, p. 19) gives the **scalar Lyapunov spectrum** for iid products on GL(m, R), under ln^+‖A‖ ∈ L_1. This closes the **leading-order scalar** question for Syracuse (already known: λ_Syracuse = log 3 − 2 log 2 ≈ −0.288).

---

## Why C4 does NOT close

Two obstructions, both inherited from v1:

### Obstruction 1 (freeness, sharpened)

Cébron Thm 2.13 **explicitly requires A and B to be FREE in (A, τ)**. The Syracuse iterated Tao recursion (Tao 1909.03562 §7.1, eq 7.5) has step operators T_j and T_{j+1} that **share the 2-adic exponent b_{[1, j]}** — they are NOT free, they are arithmetically coupled. This is **the same R77 cross-frequency v ≠ v' bilinear coupling** identified as the non-freeness signature in the v1 Tao-RMT disposition.

**v2 sharpening:** the obstruction is now **explicit** (Cébron's Thm 2.13 names the freeness hypothesis and uses it crucially in the proof via Proposition 2.1 / NC(2n) factorization, pp. 24-25), where v1 left it implicit. **What's needed is amalgamated-free probability** (Voiculescu 1995, Speicher 1998), which is **outside this corpus**.

### Obstruction 2 (profinite portability, sharpened via §C-5/C-6 analysis)

Cébron's construction of the free circular multiplicative BM (G_t) uses:
- The Laplacian Δ_{GL_N(C)} = Σ_a Z̃_a² as a sum of squares of right-invariant **smooth vector fields** on the Lie group (Cébron p. 45).
- Itô's formula in the smooth-manifold sense (Cébron eq 3.15, p. 35).
- The differentiability of the **exponential map on End(C^d)** to drop O(1/N²) corrections (Thm 4.6 proof, p. 50-51).

Profinite groups (Z/3^n)* / Z_3* have **no smooth structure**. The closest profinite analog (Vladimirov operator / p-adic Brownian motion, Kozyrev 2007) is **conceptually distant** from Cébron's free-circular BM. The "free profinite multiplicative Brownian motion" object is **research-scope to construct** (per §C-5, C-6 marked TRANSFER_RESEARCH in `C4_REPROBE_V2_PROFINITE_PORT.md`).

### Single most load-bearing portability question

> Can a "free profinite multiplicative Brownian motion" be constructed on Z_3* such that the analog of Cébron Thm 4.6 (large-cyclic-level limit) reproduces the PADE-observed spectral radius 1/√3 — AND can the freeness hypothesis of Cébron Thm 2.13 be relaxed to amalgamation over a subalgebra capturing the 2-adic/3-adic coupling?

This is **three pieces of original mathematics**, the third depending on the first two. **Research-monograph scope** — exactly the Ch 11 ceiling.

---

## Updated Ch 11 effort estimate

Per `PROFINITE_TRANSFER_OPERATOR_BLUEPRINT.md` line 222: Ch 11 (renewal Egorov, LOAD-BEARING) = 3-6 months at user pace / 9-18 months at typical pace.

**Updated:** **unchanged in magnitude** (still 3-6 months / 9-18 months), but **the next probe target is now identified specifically** — Voiculescu's 1995 operator-valued free probability paper (or Speicher 1998 Memoirs AMS), addressing the freeness obstruction directly.

If the operator-valued generalization closes the freeness obstruction without introducing new smoothness requirements, **the Ch 11 lower bound could plausibly tighten to 2-4 months at user pace** (subject to a v3 probe of Voiculescu 1995). But this is **a forward-looking estimate, not a confirmed reduction**.

---

## Specific next-probe target

**Status of C4 = C4_PARTIAL_OUTSIDE** because the missing piece is identifiable and external to this corpus:

### Primary target

**Voiculescu, D. — "Operations on certain non-commutative operator-valued random variables." Astérisque 232 (1995), pp. 243-275.** [SUSPECTED REFERENCE — citation should be verified]

This paper introduces the operator-valued / B-amalgamated free probability framework that generalizes Cébron's scalar τ to a conditional-expectation-valued E_B : A → B for B ⊂ A a subalgebra. Under amalgamation, **the freeness condition is replaced by "free with amalgamation over B"** — appropriate for Syracuse's R77 coupling (where the subalgebra B would capture the shared 2-adic exponent structure).

### Secondary target (canonical reference)

**Speicher, R. — "Combinatorial theory of the free product with amalgamation and operator-valued free probability theory." Memoirs of the American Mathematical Society, no. 627, 1998, 88 pp.**

The combinatorial / Hopf-algebra formulation of operator-valued free probability. This is the natural extension of Cébron's free log-cumulants L_κ on the Hopf algebra Y(k) to **B-valued cumulants**.

### Tertiary target (recent application)

**Das 2025 "Free Probabilistic Framework for Denoising Diffusion Models"** (in `closure hunt/`, 2510.22778v2.pdf) is a recent application of operator-valued free probability + free entropy to diffusion models. The framework is structurally relevant (it uses operator-valued free probability extensively); specific theorems are downstream. Cite as **framework-confirmation only**, NOT a load-bearing reference.

### Negative results in corpus

- **Goldsheid-Sodin 2022** (2012.03017v2.pdf, "Lower bounds on Anderson-localised eigenfunctions on a strip") — spectral-gap / Lyapunov-decay context. **Does not address composition formula.** Confirmed not load-bearing for C4.

- **Armentano-Chinta-Sahi-Shub 2024** (random-and-mean-lyapunov-...pdf, ETDS) — mean vs random Lyapunov exponents on GL_n(R) via spherical polynomials. **Relevant to GM extension, but downstream of Ch 11.** Not load-bearing for the composition formula itself.

- **Voiculescu-Dykema-Nica 1992** (scanned, unreadable) — the foundational free-probability monograph. Cited by Cébron repeatedly. The free-probability operator-algebra setup is **standard** by Cébron's time; nothing in this scanned book changes the C4 disposition. **[NON_VERBATIM_VDN]: per standard reference knowledge, VDN 1992 sets up (A, τ) and free independence but does NOT develop operator-valued / amalgamated free probability — that is Voiculescu's later 1995 work.**

---

## What changed from v1 (Tao RMT) to v2 (Cébron)

| Aspect | v1 (Tao RMT, 2026-05-14) | v2 (Cébron, 2026-05-14) |
|---|---|---|
| Multiplicative composition formula | Punted to Speicher, no statement | **Explicit: Thm 2.13** |
| Required hypothesis | Implicit (Tao RMT mentions "free random variables") | **Explicit: A and B FREE in (A, τ)** |
| Status of freeness for Syracuse | Implicit non-freeness via R77 | **Explicit: arithmetic coupling via shared 2-adic exponent** |
| Profinite portability | Not addressed (Tao RMT is single-base) | **TRANSFER_RESEARCH** at construction level (C-5, C-6) |
| Next probe target | "Speicher's survey" (loose) | **Voiculescu 1995 + Speicher 1998 (Memoirs AMS)** |
| Ch 11 lower bound | 3-6 months at user pace | 3-6 months at user pace (unchanged), potentially 2-4 months pending v3 Voiculescu probe |

**Net:** **the v1 → v2 transition tightens the picture but does not close C4.** Cébron 2013 is the *right* multiplicative-free-probability paper to read at the scalar level; the freeness obstruction is now sharp and the next step (operator-valued) is identifiable.

---

## File index

- `C:/Collatz/C4_REPROBE_V2_CEBRON_HYPOTHESES.md` — verbatim Cébron extracts (Thm 2.13, Thm 4.6, Prop 3.6, Lemma 4.3)
- `C:/Collatz/C4_REPROBE_V2_GOLDSHEID_MARGULIS_HYPOTHESES.md` — verbatim Goldsheid-Margulis extracts (Thm 1.2 MET, Thm 5.4, Thm 6.6, Thm 6.11 Kotani)
- `C:/Collatz/C4_REPROBE_V2_HYPOTHESIS_CHECK.md` — match table against Syracuse iterated recursion
- `C:/Collatz/C4_REPROBE_V2_PROFINITE_PORT.md` — component-by-component profinite portability marks
- `C:/Collatz/C4_REPROBE_V2_DISPOSITION.md` — this file
- `C:/Collatz/_cebron_pages/` — 55 page-text files from arXiv:1304.1713v3
- `C:/Collatz/_goldsheid_pages/` — 61 page-text files from Russian Math Surveys 44:5

No git operations performed.

---

## Bottom line

**C4 = Ch 11 remains the load-bearing original chapter.** Cébron 2013 is the **structurally correct paper** for the scalar-valued multiplicative free convolution; its Theorem 2.13 IS the composition formula Tao RMT punted to Speicher. But Cébron's freeness hypothesis **fails for Syracuse** (R77 cross-frequency coupling = explicit non-freeness via shared 2-adic exponents), and Cébron's construction uses **smooth-Lie-group machinery** with **no obvious profinite analog**.

**Next probe target identified:** Voiculescu 1995 "Operations on Certain Non-Commutative Operator-Valued Random Variables" (Astérisque 232) OR Speicher 1998 "Combinatorial theory of the free product with amalgamation" (Memoirs AMS 627) — the operator-valued / amalgamated free probability framework, which is precisely what relaxes Cébron's freeness hypothesis.

**The disposition is C4_PARTIAL_OUTSIDE.**
