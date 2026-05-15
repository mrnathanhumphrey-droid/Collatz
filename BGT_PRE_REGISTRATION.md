# BGT_PRE_REGISTRATION

**Date:** 2026-05-13.
**Probe:** BGT (Bingham-Goldie-Teugels) regular-variation single-theorem selection on the ε_k sequence.
**Mode:** E (verbatim theorem hypotheses from PDF; no inheritance from prior project files).
**Status:** Pre-registration locked BEFORE candidate selection.

---

## ε_k input verification (verbatim from result_77_7_eps_exact_through_k8_v2_vec_pool.json)

Exact-rational signs and normalized magnitudes (computed from the JSON):

| k | sign | |ε_k|·2^k |
|---|---|---|
| 1 | + | 0.400000 |
| 2 | + | 0.038095 |
| 3 | − | 0.040736 |
| 4 | − | 0.039236 |
| 5 | − | 0.036856 |
| 6 | − | 0.031866 |
| 7 | − | 0.150430 |
| 8 | − | 0.190860 |

Sign pattern: (+, +, −, −, −, −, −, −). Matches input spec exactly.
Plateau k=2..6 at ~0.04; jump to 0.150 at k=7; k=8 0.191 confirms post-jump regime.

## Diagnostic ratios

Ratios |ε_{k+1}| / |ε_k|:
- k=1→2: 0.0476
- k=2→3: 0.5347
- k=3→4: 0.4816
- k=4→5: 0.4697
- k=5→6: 0.4323
- k=6→7: 2.3604  ← jump
- k=7→8: 0.6344

Slow-variation test, L(k) := |ε_k|·2^k, L-ratios L(k+1)/L(k):
- k=2→3: 1.0693
- k=3→4: 0.9632
- k=4→5: 0.9393
- k=5→6: 0.8646
- k=6→7: 4.7207  ← jump
- k=7→8: 1.2688

Partial sums S_n = Σ_{k≤n} ε_k:
- n=1: +0.200000
- n=2: +0.209524
- n=3: +0.204432
- n=4: +0.201980
- n=5: +0.200828
- n=6: +0.200330
- n=7: +0.199155
- n=8: +0.198409

The partial-sum sequence converges (numerically near 0.198) but oscillates / changes direction at the k=7 jump.

---

## Candidate list (locked)

A. **Karamata's representation theorem** (Kevei Thm 6 = BGT 1.3.1)
B. **Karamata's Tauberian theorem** (Kevei Thm 18; BGT 1.7.1; Feller XIII §5)
C. **De Haan Π-variation / Π-class** (2RV-Hawkes Def 2.4; BGT Ch. 3)
D. **Second-order regular variation (2RV)** (2RV-Hawkes Def 2.2; BGT Ch. 3; de Haan-Stadtmüller)
E. **Bingham-Ostaszewski sequential regular variation** (BO 2020 Thm 1, Thm K2 Kendall)
F. **Karamata's monotone density theorem** (Kevei Thm 14)
G. **Karamata's integration theorem (Karamata's lemma)** (Kevei Thm 12, Prop 4)
H. **Goldie subexponential / regularly varying tail** (Jessen-Mikosch §4; arXiv 2001.05420)
I. **Multi-singularity extension** (looking for two-singularity-aware BGT theorems; specifically z=1 AND z=−1)
J. **Hazard-rate framework** (arXiv 2504.11655 Thm 2.4, 2.7, 2.11, 2.15)

Additional surfaced during Phase 0 reading:
K. **2RV-Hawkes Theorem 3.2 (second-order Karamata representation)** — closer in spirit to the k=7 jump signature than first-order Karamata.
L. **2RV-Hawkes Theorem 3.3 (extended second-order Karamata)** — covers integral transforms.
M. **Potter bounds** (Kevei Thm 11) — derivable bound, sub-route for conclusion-shape conversion.

---

## Pre-registered priors per candidate

| code | theorem | SELECTED | PARTIAL | NO_FIT | BLOCKER |
|---|---|---|---|---|---|
| A | Karamata representation | 25 | 35 | 30 | 10 |
| B | Karamata Tauberian | 10 | 25 | 35 | 30 |
| C | De Haan Π-variation | 15 | 30 | 40 | 15 |
| D | 2RV (second-order RV) | 30 | 35 | 25 | 10 |
| E | Sequential RV (Kendall + BO) | 30 | 30 | 30 | 10 |
| F | Karamata monotone density | 5 | 15 | 75 | 5 |
| G | Karamata integration | 10 | 30 | 55 | 5 |
| H | Goldie subexponential | 5 | 10 | 75 | 10 |
| I | Multi-singularity / FS VI.5-style | 10 | 25 | 40 | 25 |
| J | Hazard-rate framework | 5 | 15 | 65 | 15 |
| K | 2RV Karamata representation (Hawkes Thm 3.2) | 25 | 35 | 30 | 10 |
| L | Extended 2RV Karamata (Hawkes Thm 3.3) | 15 | 30 | 45 | 10 |
| M | Potter bounds | 5 | 25 | 65 | 5 |

**Rationale highlights:**
- **D, E, K (highest SELECTED priors)** — 2RV and sequential RV are the categorically right tools: 2RV admits second-order corrections to a first-order RV envelope (the k=7 jump might fit as a regime change in the auxiliary function A(t)), and sequential RV is sequence-native (ε_k IS a sequence, not a continuous function). Hawkes Thm 3.2 has the cleanest representation for second-order corrections.
- **F, H (lowest SELECTED, highest NO_FIT)** — F requires positive monotone density (ε_k is sign-mixed and non-monotone). H is subexponential / heavy-tail, wrong category (ε_k is light-tail-bounded by 2^{-k}).
- **B (high BLOCKER)** — Tauberian theorem inherits the Mode H trap from the Tauberian re-scope probe: requires Laplace-side analytic-continuation hypothesis equivalent to the closure target.
- **I (high BLOCKER)** — Multi-singularity theorems in scanned corpus (Flajolet-Sedgewick VI.5 paywalled, BGT 4.x book not pullable, Bingham 2007 CDAM very-slowly-varying focused on a different question) may not have verbatim statements in 12-PDF corpus.

---

## Decision rules

**SELECTED** — ε_k input data satisfies all extractable hypotheses (or fails NEEDS_PROOF gracefully), AND the conclusion produces:
  (i) a polynomial-in-A Fourier-decay bound on |μ̂_n(ξ)| via the R75/R76/R77 conversion chain, OR
  (ii) a structural asymptotic identification of ε_k (e.g., ε_k ~ C·r^k·L(k) for explicit r, slowly-varying L) convertible to (i) via Plancherel + R76 leading-mode identity + R77 T_diag + R78/R79 chain-side bilinear analysis.

**PARTIAL** — Hypotheses SATISFIED but conclusion falls short of (i)/(ii). Report what additional input (ε_9..ε_K) would close the gap.

**NO_FIT** — At least one hypothesis FAILED for every candidate (categorical mismatch).

**BLOCKER** — Theorem statement UNVERIFIABLE in the 12-PDF corpus, OR a load-bearing hypothesis is itself a target-object (Mode H circular — the closure target dressed as a hypothesis).

**MODE_H_CIRCULAR** — Specific sub-type of BLOCKER. Theorem applies in principle, but its load-bearing hypothesis (e.g., "f(z) extends analytically past z=1", "Laplace transform has known asymptotic at s=0") is itself the closure target.

## Locked thresholds for decision categories

- SELECTED triggers only if at least one candidate produces (i) directly OR (ii) with a chain spelled out using R74-R77 statements that are already proved.
- PARTIAL requires explicit identification of K-value (number of additional ε_k coefficients needed) and explicit prediction of which candidate would fire post-extension.
- NO_FIT requires at least one explicit failed-hypothesis citation per candidate.
- BLOCKER requires explicit identification of which theorem statement is non-locatable or which hypothesis is target-object.
- Multi-regime test (Phase 3) is gating: if a candidate's Phase 3 verdict is STRUCTURALLY_BLOCKED, the candidate disposition cannot exceed PARTIAL regardless of Phase 1-2.

---

## Pre-registered secondary routing

If overall disposition is NO_FIT or PARTIAL:

1. **Igusa local zeta / functional equation** — (1+3)^u algebraic structure in R78 D=0 disambiguation. Operates on the generating function of ε_k as a p-adic integral. Categorically distinct from BGT.

2. **Heat-kernel / Brownian motion on profinite tree** — re-pre-classified as VSC-style, pre-registered LOW priority because of compact-vs-noncompact volume growth gap (per FG_H disposition).

3. **Adelic Mellin construction** — restore archimedean place. Operates on a global object (Mellin transform on Q*_A / Q*). MODERATE priority; technically demanding.

4. **Faure 2009 semiclassical spectral gap** — for rigorous spectral-gap on R77 T_diag's abelianized FG transfer operator (downstream R78 work, not closure-route).

5. **New candidates surfaced during Phase 0:**
   a. Differential / D-Padé approximants on f̃(z) (R77.6 mentioned but didn't run; could discriminate (G-power) from (G-log) in branch-cut signature at z=2).
   b. Local-singularity asymptotic via Mellin-Barnes contour on the Padé branch-cut structure already identified in R77.6 — operates on the same object (ε_k generating function) but uses different machinery.
   c. Watson lemma / saddle-point on the bilinear off-diagonal sum (R77.3-R77.7) — operates closer to the chain-side, complementary to BGT.

**Priority ordering** (if BGT closes NO_FIT or PARTIAL):
- HIGHEST: Igusa local zeta (categorically distinct from prior 6 probes; the (1+3)^u structure is fresh).
- MODERATE: Adelic Mellin (archimedean closure target).
- MODERATE: Watson lemma / saddle-point on R78/R79 bilinear (closer to existing chain-side machinery).
- LOW: heat-kernel, Faure 2009 (covered by prior probes).
- CONTINGENT: D-Padé / Mellin-Barnes (only fires if BGT closure is PARTIAL on "need higher K", since this is just diagnostic refinement of the branch-cut at z=2, not a new closure path).

---

## Mode H awareness specific to BGT

BGT theorems typically have hypotheses of the form "f is in RV_ρ" or "f is in Π(A)" — these are OPERATIONALLY checkable from sequence/function data (via ratio limits, finite differences). These are the candidates where Phase 1 can resolve SATISFIED/FAILED directly.

The Mode H trap candidates are:
- **B (Karamata Tauberian)**: Hypothesis is "Laplace transform Û(s) has known asymptotic as s→0" — equivalent to the closure target after conversion.
- **Any theorem whose conclusion-shape conversion to |μ̂_n(ξ)| requires an analytic-continuation step not derivable from R74-R77.**

Filter strategy:
- For each candidate, identify the load-bearing hypothesis.
- If load-bearing hypothesis is "sequence/function property" — operational check, Phase 1 resolves.
- If load-bearing hypothesis is "transform-side asymptotic" — Mode H check, document as BLOCKER if transform's asymptotic is itself the closure target.

---

## Files to produce per candidate

For each K ∈ {A, B, C, D, E, F, G, H, I, J, K_2RV_rep, L_2RV_int, M_Potter}:
- `BGT_{K}_HYPOTHESES.md` — verbatim theorem statement + hypothesis types (Phase 0)
- `BGT_{K}_HYPOTHESIS_CHECK.md` — hypothesis × ε_k matrix with SATISFIED / FAILED / UNVERIFIABLE / N/A / NEEDS_PROOF (Phase 1)
- `BGT_{K}_CONCLUSION_SHAPE.md` — conclusion vs |μ̂_n(ξ)| target (Phase 2)
- `BGT_{K}_MULTIREGIME_CHECK.md` — k=7 jump accommodation (Phase 3)

Final disposition: `BGT_DISPOSITION.md`.

---

End of pre-registration. Locked.
