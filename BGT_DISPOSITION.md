# BGT_DISPOSITION

**Date:** 2026-05-13.
**Probe:** BGT (Bingham-Goldie-Teugels) regular-variation single-theorem selection on the ε_k sequence.
**Mode:** E (verbatim theorem hypotheses from PDF, no inheritance from prior project files).
**Pre-registration:** `C:/Collatz/BGT_PRE_REGISTRATION.md`.

---

## Headline

**No SELECTED.** The BGT corpus closes **NO_FIT-dominant with one PARTIAL (E, sequential RV)** and one BLOCKER (I, multi-singularity, UNVERIFIABLE-PHASE-0).

The k=7 jump in |ε_k|·2^k is the load-bearing obstruction across every candidate: single-regime BGT theorems (Karamata representation, Karamata theorem, 2RV, Π-variation, Potter bounds, monotone density, hazard-rate) categorically fail because their hypothesis (slow variation / first-order RV / Π-class) is violated by L(6→7) = 4.72 — emphatically not in [1−δ, 1+δ] for any small δ. Sequence-native candidate E (Kendall / BO sequential RV) PARTIAL because the convergence hypothesis a_n·f(x_n·λ) → g(λ) cannot be verified at N=8 across the jump; ε_9..ε_K for K ≥ 15-20 would resolve.

The BGT re-categorization to ε_k sequence object DID side-step the chain-side category obstructions (no group structure, reversibility, proximality required, per FG_DISPOSITION's pre-clearance rationale) — but ran into the **same multi-regime obstruction** that closed TAUBERIAN_RESCOPE. The k=7 jump is structurally incompatible with single-regime asymptotic frameworks across two independent literatures (Tauberian + BGT regular variation).

This **completes the seventh category-of-object barrier** in the systematic obstruction map for c=7/45 closure.

---

## Summary table

| Code | Theorem | Phase 0 | Phase 1 | Phase 2 | Phase 3 | Disposition | File |
|---|---|---|---|---|---|---|---|
| A | Karamata representation theorem (Kevei Thm 6) | Extracted | FAILED at h_2 (jump) | SHAPE_MISMATCH | STRUCTURALLY_BLOCKED | **NO_FIT** | BGT_A_HYPOTHESES.md |
| B | Karamata Tauberian theorem (Kevei Thm 18) | Extracted | FAILED h_1 + MODE_H | SHAPE_MISMATCH (bare-asymptotic, no rate) | STRUCTURALLY_BLOCKED (for remainder version) | **MODE_H_CIRCULAR + NO_FIT** | BGT_B_HYPOTHESES.md |
| C | De Haan Π-variation (Hawkes Def 2.4) | Extracted | FAILED at h_2 (additive jump) | SHAPE_PARTIAL (log correction, not poly-A) | STRUCTURALLY_BLOCKED | **NO_FIT** | BGT_C_HYPOTHESES.md |
| D | Second-order RV (Hawkes Def 2.2) | Extracted | FAILED h_1 (|ε_k| not RV) | SHAPE_MISMATCH (prefactor refinement) | STRUCTURALLY_BLOCKED | **NO_FIT** | BGT_D_HYPOTHESES.md |
| E | Sequential RV (BO Thm K2, Kendall) | Extracted | FAILED h_2 across jump (plateau-only SATISFIED) | SHAPE_PARTIAL | STRUCTURALLY_BLOCKED at N=8 (would unblock with K ≥ 15-20) | **PARTIAL** | BGT_E_HYPOTHESES.md |
| F | Karamata monotone density (Kevei Thm 14) | Extracted | FAILED h_3 (non-monotone) | N/A | N/A | **NO_FIT (categorical)** | BGT_F_HYPOTHESES.md |
| G | Karamata integration theorem (Kevei Thm 12) | Extracted | FAILED h_1 (same as D) | SHAPE_MISMATCH | STRUCTURALLY_BLOCKED | **NO_FIT** | BGT_G_HYPOTHESES.md |
| H | Goldie subexp / RV tail (Jessen-Mikosch) | Extracted | FAILED h_1 + h_2 (categorical) | N/A | N/A | **NO_FIT (categorical)** | BGT_H_HYPOTHESES.md |
| I | Multi-singularity / two-singularity ext | UNVERIFIABLE-PHASE-0 | — | — | — | **BLOCKER** | BGT_I_HYPOTHESES.md |
| J | Hazard-rate framework (arxiv 2504.11655) | Extracted | FAILED h_1 (categorical) + von Mises at jump | N/A | N/A | **NO_FIT (categorical)** | BGT_J_HYPOTHESES.md |
| K | 2RV Karamata representation (Hawkes Thm 3.2) | Extracted | FAILED h_1 (inherits D) | SHAPE_MISMATCH | STRUCTURALLY_BLOCKED | **NO_FIT** | BGT_K_HYPOTHESES.md |
| L | Extended 2RV Karamata (Hawkes Thm 3.3) | Extracted | FAILED h_1 + h_2 | N/A | N/A | **NO_FIT** | BGT_L_HYPOTHESES.md |
| M | Potter bounds (Kevei Thm 11) | Extracted | FAILED at jump (same as A) | SHAPE_MISMATCH | STRUCTURALLY_BLOCKED | **NO_FIT** | BGT_M_HYPOTHESES.md |

Total candidates with extractable Phase-0 statements: **12 of 13** (only I UNVERIFIABLE-PHASE-0).

Net dispositions:
- 1 PARTIAL: E (Sequential RV / Kendall)
- 1 MODE_H_CIRCULAR + NO_FIT: B
- 9 NO_FIT: A, C, D, F, G, H, J, K, L, M (10 if M counted; corrected: A, C, D, F, G, H, J, K, L, M = 10 NO_FIT)
- 1 BLOCKER (UNVERIFIABLE-PHASE-0): I

Re-count: 1 PARTIAL + 1 MODE_H_CIRCULAR + 10 NO_FIT + 1 BLOCKER = 13. ✓

---

## Final disposition: **PARTIAL** (candidate E only; all others NO_FIT or BLOCKER)

The single PARTIAL candidate is **Bingham-Ostaszewski sequential regular variation (Kendall's Theorem K2)**:

- Phase 0: extracted verbatim.
- Phase 1 (plateau-restricted): would fire cleanly within k=2..6 plateau with L(k) slowly varying at index ρ = 0 (constant ~0.04). a_n·L(n·λ) → g(λ) achievable on N ≤ 6 with a_n ≡ 1 and g ≡ 0.04 ✓.
- Phase 1 (full sequence k=1..8): fails at h_2 due to k=7 jump.
- Phase 2: SHAPE_PARTIAL — would deliver sequence-native slow-variation identification, but conversion to polynomial-in-A bound on |μ̂_n(ξ)| via R75/R76/R77 yields prefactor refinement only, not a new rate.
- Phase 3: STRUCTURALLY_BLOCKED at N = 8 data; would unblock with ε_9..ε_K for K ≥ 15-20 IF post-jump regime stabilizes to a new plateau. If post-jump regime continues to escalate (L(k) growing for k ≥ 9), then PARTIAL closes to NO_FIT.

**Realizable additional input that would convert PARTIAL → SELECTED-eligible: ε_9..ε_K via running k=9, k=10 Markov chains.** Per R77.6, each ε_k computation scales O(N^3) on Gauss elimination over Q with N = 3^{k-1}; k=7 took hours, k=8 was the input here, k=9 would scale 27× further (~tens of hours), k=10 would be ~hundreds of hours. Realistically deliverable for k=9 and possibly k=10; k ≥ 11 unrealistic without algorithmic improvements.

**Even with K ≥ 15-20, the SELECTED branch is contingent on the post-jump regime stabilizing.** No theorem in the BGT corpus accommodates a true two-regime asymptotic (plateau-then-new-plateau or plateau-then-escalation) — that would require either (a) a multi-singularity Tauberian theorem (candidate I, currently BLOCKER), or (b) a "piecewise sequential RV" extension not in BO 2020.

---

## What category of theorem is missing

Same description as TAUBERIAN_RESCOPE Section "What category of theorem is missing":

The closure target requires:
1. A theorem operating on a sequence object (BGT sequential RV achieves this — good).
2. That handles **multi-regime asymptotics** (the k=7 jump in ε_k) — NO BGT theorem in corpus does this.
3. That delivers a **polynomial-in-A rate refinement** on the (1/2)^n envelope — the natural BGT output is slowly-varying-prefactor refinement, not rate.
4. Whose hypotheses are checkable from a finite sequence (ε_k available at k ≤ 8 currently).

The BGT corpus achieves (1) and (4) cleanly via candidate E (Kendall / BO), satisfying the FG_DISPOSITION rationale that BGT side-steps the chain-side category obstructions. The corpus FAILS at (2) — every BGT theorem assumes single-regime — and FAILS at (3) — BGT delivers prefactor refinement, not rate enhancement.

**Missing-category description:** a theorem in the BGT family that admits **multi-regime sequence asymptotics with explicit second-order rate corrections** (not just second-order *prefactor* corrections). The 2RV literature (Hawkes 2RV paper et al.) refines prefactors of an already-known first-order rate; it does not produce the rate.

---

## Surprises in the inputs

### Surprise 1: BGT slow-variation hypothesis hits the SAME k=7 obstruction as Tauberian Chevalier 1.16

The TAUBERIAN_RESCOPE finding that "no integer M ≥ 1 fits the k=7 jump in Chevalier 1.16's b_n ~ D·n^{M−3/2}" is structurally parallel to the BGT finding that "L(k) is not slowly varying across the k=7 jump." Both literatures encode multi-regime asymptotics via a single-singularity / single-regime hypothesis, which the empirical ε_k data violates.

The k=7 jump is now **independently confirmed as the load-bearing obstruction across two separate framework families** (complex-Tauberian + BGT regular variation), which strongly suggests it is a *structural* feature of the Syracuse μ_n stationary measure, not a finite-N artifact.

### Surprise 2: Within plateau k=2..6, BGT slow variation holds well

L(k) values for k=2..6: 0.038, 0.041, 0.039, 0.037, 0.032. L-ratios L(k+1)/L(k) drift from 1.07 down to 0.86 with a clean downward trend. Within-plateau, this is consistent with slow variation with a *gentle secular drift* — i.e., possibly 2RV with auxiliary A(t) decaying slowly. The 5-point fit is too short to confirm 2RV, but the trend is plausible.

If the post-jump regime k ≥ 7 turns out to be a *second* plateau at higher L, then BGT sequential RV would fire on the union of two plateaus via a hypothetical multi-regime extension (not in literature). If post-jump escalates further, then no BGT extension fires.

### Surprise 3: The k=7 jump is exactly at the threshold where R77 conjectured a third spectral mode

R77.5/R77.6 documented a branch-cut singularity at z=2 in the generating function f̃(z), with type indeterminate at N=5 coefficients. The k=7 transition is suspiciously coincident with the level where R77.5's renormalization step would predict a *third mode* surfacing beyond the (1,4) and (1,-1) modes of T_diag. This is structural, not coincidental:

- R77 T_diag two-mode eigenstructure {0, 1} on (1,−1) and (1,4) is a 2×2 reduction.
- The off-diagonal rate-½ correction is a *third mode* at level n+1 to level n; R77.6 calls this branch-cut at z = 2 = 1/(1/2).
- The k=7 jump might be the level at which the third mode's contribution surfaces above the noise of the two-mode envelope — consistent with both candidate D (2RV at second order) and candidate E (sequential RV with a regime change).

**Interpretation:** the empirical k=7 jump is plausibly the signature of R77.6's branch-cut at z=2 becoming visible in ε_k partial-data. If this is correct, K ≥ 15-20 coefficients would reveal a *single new plateau* at the higher L level — making the BGT sequential RV PARTIAL → SELECTED-eligible (modulo multi-regime extension being constructible).

### Surprise 4: All natural BGT first-order RV candidates fail at the |ε_k| level categorically

The reason is straightforward but worth flagging: |ε_k| ≤ const · 2^{-k} is *geometric* decay, which in BGT terminology is *rapid* variation (de Haan), NOT regular variation with finite index. First-order RV requires power-law decay (|ε_k| ~ C·k^{-α}·ℓ(k)). The 2^k normalization removes the geometric factor and exposes the slowly-varying / 2RV structure in L = |ε_k|·2^k — this re-categorization is correct, and BGT operates on L not on |ε_k| directly.

This means the candidate-D and candidate-G "h_1 categorically fails" verdicts are real but cleanable: the right BGT input is L, not |ε_k|. With L: A, C, E, M apply (and fail at the k=7 jump for the same reason).

---

## SECONDARY ROUTING

Per pre-registration, the routes flagged in advance:

1. **Igusa local zeta / functional equation** — operates on the (1+3)^u 3-adic structure of R78. Categorically distinct from BGT (operates on a generating-function p-adic integral). **PRIORITY: HIGH.** The (1+3)^u algebraic root in R78's D=0 disambiguation is concrete and unexplored. Igusa Z(s; f) = ∫_{Q_3} |f(x)|^s dx gives meromorphic continuation + functional equation; if ε_k's generating function = Igusa local zeta of a specific f, the closure rate emerges from meromorphic-continuation poles. Requires explicit identification of f, which is plausible given R78's polynomial identification.

2. **Adelic Mellin construction** — restore archimedean place per BT_DISPOSITION. **PRIORITY: MODERATE-HIGH.** Technically demanding; payoff is the only route explicitly addressing archimedean-visibility. May connect with c=7/45 being structurally archimedean.

3. **Watson lemma / saddle-point on R78/R79 bilinear** — operates closer to chain-side. **PRIORITY: MODERATE.** Complementary to BGT; gives saddle-point asymptotic for the bilinear off-diagonal sum, which is where the k=7 third-mode contribution lives. Doesn't directly close c=7/45 but tightens the structural picture.

4. **Faure 2009 semiclassical spectral gap** — downstream R78 work, NOT closure-route. **PRIORITY: LOW for closure.**

5. **Heat-kernel on profinite tree** — already pre-classified as STRUCTURALLY_BLOCKED in FG_DISPOSITION. **PRIORITY: LOW.**

6. **Differential / D-Padé approximants on f̃(z)** — refinement of R77.6's branch-cut diagnosis. Would discriminate (G-power) from (G-log) at z=2. **PRIORITY: CONTINGENT (only if BGT PARTIAL becomes the active line) — would tighten what additional ε_k data is needed.**

7. **NEW candidate surfaced during BGT probe:** A **two-plateau / regime-change Tauberian** theorem — explicitly multi-regime sequence asymptotic. No verbatim statement found in 12-PDF BGT corpus or 20-PDF Tauberian corpus. Would need to scan Flajolet-Sedgewick Ch. VI (multi-singularity), BGT book Ch. 4 (Abelian-Tauberian extensions), or Korevaar's Tauberian Theory (specialized chapters). **PRIORITY: MODERATE-LOW** (corpus would need extension).

### Recommended top-priority secondary route: **Igusa local zeta**.

Rationale:
- Categorically distinct from BGT (which operates on the sequence-level slow-variation structure) — addresses a *different aspect* of ε_k (its 3-adic / algebraic structure from R78's polynomial identification).
- The (1+3)^u algebraic root is concrete and unexplored: Igusa Z(s; (1+3)^u − c) for c the algebraic root would directly probe the analytic continuation.
- Operates on a *single complex variable s* (no chain-side category mismatch, no proximality, no flag variety).
- Bypasses the multi-regime obstruction entirely: Igusa local zeta delivers a meromorphic function with poles encoding the rate, regardless of whether ε_k has plateau-then-jump structure.

---

## Files produced

- `C:/Collatz/BGT_PRE_REGISTRATION.md` (pre-reg, locked before selection)
- `C:/Collatz/BGT_A_HYPOTHESES.md` (Karamata representation)
- `C:/Collatz/BGT_B_HYPOTHESES.md` (Karamata Tauberian)
- `C:/Collatz/BGT_C_HYPOTHESES.md` (de Haan Π-variation)
- `C:/Collatz/BGT_D_HYPOTHESES.md` (2RV)
- `C:/Collatz/BGT_E_HYPOTHESES.md` (Sequential RV — the PARTIAL)
- `C:/Collatz/BGT_F_HYPOTHESES.md` (Monotone density)
- `C:/Collatz/BGT_G_HYPOTHESES.md` (Karamata integration)
- `C:/Collatz/BGT_H_HYPOTHESES.md` (Goldie subexp)
- `C:/Collatz/BGT_I_HYPOTHESES.md` (Multi-singularity — UNVERIFIABLE)
- `C:/Collatz/BGT_J_HYPOTHESES.md` (Hazard rate)
- `C:/Collatz/BGT_K_HYPOTHESES.md` (2RV Karamata representation, Hawkes Thm 3.2)
- `C:/Collatz/BGT_L_HYPOTHESES.md` (Extended 2RV Karamata, Hawkes Thm 3.3)
- `C:/Collatz/BGT_M_HYPOTHESES.md` (Potter bounds)
- `C:/Collatz/BGT_DISPOSITION.md` (this file)

PDF extractions (UTF-8 from pypdf): `C:/tmp/bgt/*.txt` (12 PDFs, all extracted; 2 had null bytes cleaned).

No git operations performed (per discipline).

---

End disposition.
