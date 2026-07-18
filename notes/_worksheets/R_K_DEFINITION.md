# R_K_DEFINITION — explicit articulation of R_k from R77.4 erratum §1

**Date:** 2026-05-11. Phase 1 of `R_K_*` probe. Articulates the inter-level residual operator R_k proposed in R77.4 erratum §"What this DOES change — reframing" as a candidate replacement for the falsified T_3. Companion files: `R_K_APPROACH_A.md`, `R_K_APPROACH_B.md`, `R_K_APPROACH_C.md`, `R_K_CLOSURE_TABLE.md`, `R_K_DISPOSITION.md`.

## 0. Where R_k comes from

R77.4 erratum (`result_77_4_K_spectrum_erratum.md`) discovered that the natural within-level Markov transition K_k has **no spectral feature near 1/2 at any k ∈ {3..6}** (closest eigenvalue at every level is λ_2 ≈ 0). Erratum §"What this DOES change — reframing" therefore reframed the rate-1/2 question:

> "The right operator-theoretic question is no longer 'what spectrum does K_k have near 1/2?' but rather: **What is the operator governing the inter-level refinement π_k → π_{k+1}, and what does its spectrum look like near 1/2?**"
>
> "The state space at level k+1 is a 3-fold refinement of level k's: each coprime residue r mod 3^k corresponds to {r, r+3^k, r+2·3^k} mod 3^{k+1}. The 'lift' operator L_{k→k+1}: R^{N_k} → R^{N_{k+1}} sends π_k to the uniform extension over the 3 fibers; the actual π_{k+1} differs from L_{k→k+1}·π_k by a residual whose magnitude empirically scales like (1/2)^n."
>
> "The spectrum of this **residual operator** (or its companion on the projective limit) is what R77.x was actually trying to characterize."
>
> [R77.4 erratum §"What this DOES change — reframing"; companion §"Recommended next moves" item 1: "Build the inter-level residual operator R_k: π_k ↦ π_{k+1} − L_{k→k+1}·π_k (after suitable embedding to a common space) and compute its spectrum across k."]

## 1. Explicit definition (R77.5 articulation)

R77.5 (`result_77_5_inter_level_residual.md`) carried out R77.4 erratum's recommendation. The construction proceeded in two parts:

### 1.1 The natural lift T_{k→k+1}

> T_{k→k+1}: V_k → V_{k+1},   T(π_k)(r') := π_k[r' mod 3^k] / 3   for r' coprime in Z/3^{k+1}.
>
> where V_k := R^{N_k}, N_k := 2·3^{k−1} (number of coprime residues mod 3^k), and r' ranges over coprime residues mod 3^{k+1}.

This is the cardinality-uniform extension: each coprime r ∈ Z/3^k has exactly 3 coprime preimages r' ∈ Z/3^{k+1} (namely r, r+3^k, r+2·3^k), each receiving mass π_k(r)/3.

T_{k→k+1} is an **isometric embedding up to a √3 factor**:

> ⟨T(u), T(v)⟩_{V_{k+1}} = ⟨u, v⟩_{V_k} / 3.

### 1.2 The lift residual

> **R_k(r') := π_{k+1}(r') − T(π_k)(r')**       (r' coprime in Z/3^{k+1})

This is a **vector**, not an operator-on-a-fixed-Hilbert-space. Each R_k lives in V_{k+1}, depending on k. Properties (R77.5 §§2-3, all verified over Q via fractions.Fraction):

- **Σ_{r'} R_k(r') = 0** (signed total).
- **R_k ∈ W_k := T(V_k)^⊥** by construction (vanishes on each 3-fiber {r, r+3^k, r+2·3^k}, by marginal consistency of the projective Markov system).
- **‖R_k‖² = 0.155 · 3^{−k} + O(3^{-k} · 2^{-k})** for k = 1..5; ratio ‖R_k‖² / ‖R_{k−1}‖² → 0.334... ≈ **1/3** (not 1/4).
- ‖R_k‖² = ‖d_{k+1}‖²_R74 as an **exact rational identity** (R77.5 follow-up `result_77_5_d_R_identity_check.md` 5/5 PASS). R_k is the R74 deviation vector in clean wavelet-like notation; no new mathematical content beyond R74 in the L²-norm channel.

### 1.3 The "operator" Φ — domain/codomain mismatch

If we want to read R_k as the action of an operator, the candidate is

> **Φ : W_{k−1} → W_k,   R_{k−1} ↦ R_k**

But:

- W_{k−1} ⊂ V_k is **different from** W_k ⊂ V_{k+1}. They are different Hilbert spaces, with dim(W_k) = 4·3^{k−1} growing by factor 3 per level.
- There is **no single Hilbert space on which Φ acts at finite truncation**. Φ is a **sequence** of inter-level maps, not a single operator with a spectrum in the operator-theoretic sense.
- The natural setting where Φ becomes a single operator is the **projective limit** (inverse limit) on Ẑ_3^×, where {V_k} → L²(Ẑ_3^×, π_∞) and {W_k} → orthogonal-complement filtration. The level-by-level decomposition is wavelet-like (Haar on Ẑ_3^× adapted to the {W_k} filtration).
- **No finite-dimensional matrix realization of Φ exists** — only level-by-level R_{k−1} → R_k. Computing Φ's "matrix" requires a basis for {W_k} and an explicit description of the map R_{k−1} ↦ R_k in that basis, which is a substantial multi-resolution analysis (R77.5 §4 flagged this as "substantial, separate work").

This is the critical structural fact: **R_k is not a self-map of one space.**

## 2. Relationship to T_3 (the falsified within-level operator)

| | T_3 (R77.2, falsified) | R_k (R77.5, this probe) |
|---|---|---|
| Type | 3×3 matrix on C^3 | Vector R^{N_{k+1}} (sequence indexed by k) |
| Domain/codomain | C^3 → C^3 (self-map) | W_{k−1} ⊂ V_k → W_k ⊂ V_{k+1} (level-dependent) |
| Spectrum | {1/2, 1/4, 1/8} (algebraic fact of matrix) | Not directly defined; only via Φ sequence (no spectrum without single Hilbert space) |
| Operator norm | ≤ 944 (M_3 probe, against falsified spectrum) | Per-level: ‖R_k‖_{L²} ∼ √0.155 · 3^{−k/2}; "Φ" has no global operator norm |
| Rate-1/2 status | Recursion falsified by R77.3 (28% residual at n=1) | Norm ratio stabilizes at **1/3, NOT 1/4** (R77.5 §2) |
| Connection to ε_n | Was supposed to encode ε_n's recursion (falsified) | R_k = d_{k+1} (R74 deviation vector, exact rational identity); ε_n lives in moment-functional projection, NOT in ‖R_k‖ itself |
| Reduces to other in limit? | No — different mathematical type (matrix vs vector sequence) | No |

R_k does **NOT** reduce to T_3 in any limit. They are categorically different objects. R_k is the right object structurally (it sits on the inter-level lift residual, where R77.4 erratum suggested looking) but its mathematical type is "sequence of vectors in different orthogonal-complement subspaces", not "operator with spectrum on a fixed Hilbert space."

## 3. Relationship to ε_n via Nisoli framework

R77.5 §5 established the decomposition

> π_n − π_∞ = Σ_{k ≥ 1} (lift to level n of R_k component)
>
> ε_n = ⟨φ_n, π_n − π_∞⟩ = Σ_{k} ⟨φ_n, lift_n(R_k)⟩

The W_k subspaces are mutually orthogonal in V_n (after lift). Per-level mass: ‖lift_n(R_k)‖² = ‖R_k‖² · 3^{n−k−1} = 0.155 · 3^{n−2k−1}. **Sum over k diverges** (3^{n−1} per level for k=0 floor), so ε_n's rate-1/2 cannot be carried by R_k's norm alone.

**The rate-1/2 lives in the projection structure ⟨φ_n, lift_n(R_k)⟩**, which is a property of the moment functional φ_n (bilinear pair-form from R76), not of R_k itself. R77.5 §5 explicitly stated:

> "Rate-1/2 of ε_n cannot live in ‖R_k‖ — it must live in the projection of the bilinear pair-form moment functional φ_n onto Σ_k W_k. R77.2's 'find a 1/2-eigenvalue at finite truncation' framing is structurally displaced; the right framework is multi-resolution / transfer-operator analysis on the projective limit."

This is the structural obstruction. **The Nisoli framework requires a single operator T with a spectrum; R_k provides a multi-resolution filtration, not a single operator.**

## 4. Structural properties of R_k

(a) **Type:** Sequence of vectors {R_k}_{k≥1}, R_k ∈ V_{k+1}, where V_j = R^{N_j} with N_j = 2·3^{j−1}. Live in mutually distinct ambient spaces.

(b) **Norms:** Exact rationals at k=1..5; ‖R_k‖² = 10/189 at k=1, 31370/1835001 at k=2, etc. Pattern: ‖R_k‖² · 3^k → 7/45 (rigorous identity inherited from R74's S_∞ = 7/15 conjecture).

(c) **Orthogonality structure:** R_k ∈ W_k := T(V_k)^⊥. The decomposition V_{k+1} = T(V_k) ⊕ W_k is a clean orthogonal direct sum (R77.5 §3.1).

(d) **Inter-level regression:** ⟨R_k, T(R_{k−1})⟩ = 0 exactly over Q at every k = 2, 3, 4, 5 (R77.5 §3, c_k = 0/1). This is **structural**, from marginal consistency Σ_{r' lifts of r} π_{k+1}(r') = π_k(r), not a learned dynamical fact. The "operator" Φ has no preferred lift-basis projection — its image in T(W_{k−1}) embedded into V_{k+1} is zero by construction.

(e) **Self-adjointness/normality:** N/A — Φ is not a self-map of a single space, so "self-adjoint" and "normal" don't apply.

(f) **Boundedness:** ‖R_k‖_{L²} → 0 like 3^{−k/2}; this is "cardinality contraction" (mass spreads into more coordinates), not eigenvalue contraction.

(g) **Resolvent (z·I − R_k)^{−1}:** undefined for the same reason — R_k is a vector, not an operator, so no resolvent. For Φ as a sequence W_{k−1} → W_k, "resolvent on a contour around λ = 1/2" doesn't exist because the domain and codomain are different Hilbert spaces.

## 5. R77.4 erratum §1 verbatim fidelity check (Phase 4 A1)

Erratum §"What this DOES change — reframing" verbatim:

> "What is the operator governing the inter-level refinement π_k → π_{k+1}, and what does its spectrum look like near 1/2?"

Erratum §"Recommended next moves" item 1 verbatim:

> "Build the inter-level residual operator R_k: π_k ↦ π_{k+1} − L_{k→k+1}·π_k (after suitable embedding to a common space) and compute its spectrum across k. Spectral feature near 1/2 here would be the actual empirical signature R77.x was reaching for."

**Fidelity audit:**

- "**π_k ↦ π_{k+1} − L_{k→k+1}·π_k**" — matches R77.5's R_k := π_{k+1} − T(π_k) exactly (T_{k→k+1} = L_{k→k+1}).
- "**after suitable embedding to a common space**" — this is where the erratum's articulation is **inherently ambiguous**. The two natural readings are:
  - **Reading A:** Embed each R_k (which lives in V_{k+1}) into a common ambient space (e.g., the projective limit's L²(Ẑ_3^×, π_∞)), then view {R_k} as an orthogonal filtration in a single space, and ask for the spectrum of "the operator that maps level-k components to level-(k+1) components." Under this reading, Φ : ⊕_k W_k → ⊕_k W_k shifts level-k to level-(k+1); spectrum on a single space exists in principle but requires construction of the L²(Ẑ_3^×) framework (not yet done in the project).
  - **Reading B:** Treat each R_k as a finite-dimensional vector and ask for the spectrum of Φ_k := R_k → R_{k+1} as a level-by-level map (a matrix W_{k−1} → W_k for each k). Under this reading, Φ is not a single operator; it's a sequence of inter-level transfer matrices.
- "**compute its spectrum across k**" — this phrasing leans more toward Reading B (spectrum-as-it-varies-with-k), but spectrum of a non-self-map across mismatched Hilbert spaces requires choosing a basis and asking about singular values, not eigenvalues.

**Verdict:** The erratum's articulation is ambiguous between "embed into common projective-limit space" (Reading A — requires substantial separate construction; not done in the project; cf. R77.5 §7) and "level-by-level matrix in some basis" (Reading B — finite-dimensional but spectrum-meaning unclear, see Phase 2). R77.5 chose to compute ‖R_k‖² and the lift regression c_k as the most natural finite-level proxies. **Phase 2 attempts to compute R_k's spectrum under Reading B** (level-by-level matrix in a canonical basis), since Reading A is not constructable inside this probe's scope.

The ambiguity is flagged as a primary obstruction (A1).

## 6. What's known about R_k that R77.5 already established

The R77.5 work (already done, on disk) gives:

- **‖R_k‖² values at k=1..5** as exact rationals (`result_77_5_R_k_norms.csv`).
- **Ratio stabilizes at 1/3** (NOT 1/4). **Disqualifies R_k's L²-norm sequence as carrier of rate-1/2.**
- **c_k = 0 exactly** — Φ has no projection onto T(R_{k−1}), the natural "lift basis" direction.
- **R_k = d_{k+1} (R74)** as exact rational identity, so this is geometric reframing of R74, not new content.

This already strongly suggests **H_R_K_INTRACTABLE** at the structural level — the operator R77.4 erratum proposed doesn't have the type-signature that Nisoli closure needs.

## 7. Going into Phase 2

Phase 2 will compute (under Reading B) the singular values of the inter-level transfer map Φ_k : W_{k−1} → W_k for k = 2..5. If a singular value near 1/√2 emerges (so σ² ≈ 1/2, the "spectrum near 1/2" target), Reading B may yet rescue the framework. If singular values cluster near 1/√3 (cardinality scaling) with no σ² ≈ 1/2 feature, Reading B is also falsified and the disposition is H_R_K_INTRACTABLE on structural grounds.

This is **Approach A** (direct numerical "spectrum"). Approach B (Neumann series) and Approach C (resolvent norm) are subsequent only if Approach A surfaces something to bound.

## 8. Files

- `R_K_DEFINITION.md` (this file) — Phase 1 articulation
- `R_K_APPROACH_A.md` — Phase 2A direct numerical spectrum
- `R_K_APPROACH_B.md` — Phase 2B perturbation (if A insufficient)
- `R_K_APPROACH_C.md` — Phase 2C resolvent norm (if A, B insufficient)
- `R_K_CLOSURE_TABLE.md` — Phase 3 parameterized closure
- `R_K_DISPOSITION.md` — top-level disposition

Cross-references:

- `result_77_4_K_spectrum_erratum.md` — origin of the R_k proposal
- `result_77_5_inter_level_residual.md` — main R_k construction and analysis
- `result_77_5_d_R_identity_check.md` — R_k = d_{k+1} identity proof
- `result_77_5_R_k_norms.csv` — exact norm data k=1..5
- `result_77_5_phi_correlations.csv` — c_k = 0 exact data k=2..5
- `M3_DEFINITION.md` — companion probe (M_3 against T_3, falsified)
- `M3_DISPOSITION.md` — H_M3_INTRACTABLE precedent
