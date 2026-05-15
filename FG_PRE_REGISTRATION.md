# FG_PRE_REGISTRATION

**Date:** 2026-05-13.
**Probe:** Furstenberg-Guivarc'h random-walks-on-locally-compact-groups expressibility test for Syracuse mu_infinity Fourier decay.
**Mode:** E (verbatim theorem hypotheses from PDF, no inheritance from prior project files — re-read everything).
**Commit policy:** Nathan commits manually; this file is written but NOT pushed.

---

## Context (carried in once, then dropped)

Six prior probe arcs closed on category-of-object grounds (5-probe modern Fourier-decay, Cluster C1, Cluster C2, Bruhat-Tits adelic, Tauberian re-scope). The category-correct candidate framework — random walks on locally compact / profinite groups — has not yet been probed.

Syracuse mu_n(xi) lives on the chain (Z/3^n)* under a 2-adic Geom(2) step distribution; the natural ambient is the profinite group lim_n (Z/3^n)* (a Z_3-module of rank ≤ 1, isomorphic to Z_3 x Z/2). Furstenberg-Guivarc'h theory operates on this category natively for several sub-cases; whether any specific theorem accepts Syracuse mu_infinity as input and yields polynomial-in-A decay is the probe's question.

---

## Inputs (verbatim references, read at session start)

- (1) Syracuse mu_n form: `C:/Collatz/C1_TAO_RECURSION_FORM.md` — Tao 1909.03562 Section 7.1, mu_n(xi) = E chi(2-adic exponential of Geom(2)^n tuple-space), product over 2D renewal walk after pair-grouping; **additive character only**, summation domain is Geom(2)^n tuple-space NOT (Z/3^n)*, phase function 2-adic exponential in tuple variables.
- (2) C1 renewal-walk structural decomposition: `C:/Collatz/C1_DISPOSITION.md` — Tao's bound is product of conditional expectations along a 2D walk indexed by partial sums b_{[1,j]}; **NOT a standard i.i.d. renewal walk**; the "renewal" walk is 2-D with phase factor evaluated at 3^{2j-2} * 2^{-b_{[1,j]}}.
- (3) R75/R76/R77 chain-side rigorous data:
  - R75 (`c_seven_forty_fifth.md`): Plancherel decomposition S_n = Sum_{xi : 3 nmid xi} |mu_n(xi)|^2.
  - R76 (`result_76_conservation_law.md`): conservation Sum_{j=0,1,2} M_{n+1}(eta_0 + j*3^n) = 0; leading-mode identity S_{n+1} = -2 * M_{n+1}(1+3^n) = -2 * M_{n+1}(1+2*3^n).
  - R77 (`result_77_T_lead_spectrum.md`): T_diag = (1/5)*[[1,1],[4,4]] with eigenvalues {0,1}, eigenvectors (1,-1) and (1,4); off-diagonal rate 1/2 empirically verified k=2..6.
- (4) eps_k exact-rational sequence k=1..8: `experiments_output/result_77_7_eps_exact_through_k8_v2_vec_pool.json`. Normalized |eps_k|*2^k sequence with values approximately:
  - k=1: 0.400, k=2: 0.019, k=3: 0.041, k=4: 0.039, k=5: 0.037, k=6: 0.032, k=7: **0.150** (sharp jump), k=8: **0.191** (confirmed).
  - Sign pattern: (+, +, -, -, -, -, -, -).
  - **Multi-regime asymptotic structure** (single stationary measure with rate-1/2 envelope would predict monotone behavior; the k=7 jump breaks monotone single-pole pattern).

---

## Candidate list (locked before selection)

Candidates A-H from brief, plus any added during Phase 0 reading. Each candidate is named by primary theorem family; specific theorem statements extracted in Phase 0.

| Code | Candidate (theorem family) | Primary corpus PDF |
|---|---|---|
| A | Furstenberg-Kesten / Furstenberg formula (Lyapunov on random products) | Furstenberg_1963_Noncommuting_Random_Products.pdf |
| B | Guivarch-Le Page renewal theory (asymptotic moments stationary measure) | Guivarch_1980_Loi_Grands_Nombres.pdf + Le_Page_1982 + Guivarch_Raugi_1985 |
| C | BFLM 2007/2011 (quantitative equidistribution / Fourier decay stationary measure on p-adic / torus) | BFLM_2007 + BFLM_2011 |
| D | Benoist-Quint (stationary measures on homogeneous spaces, quantitative) | Benoist_Quint_Random_Walks_Reductive_Groups + Benoist_Quint_2016 |
| E | Bourgain discretized sum-product Fourier decay on R^n | Discretized_SumProduct_Fourier_Decay_Rn.pdf |
| F | Le Page CLT / local CLT for products of group elements | Le_Page_1982 + Li_Fourier_Decay_Renewal_Spectral_Gaps |
| G | Saloff-Coste / Diaconis-Saloff-Coste (mixing of finite groups, Nash / log-Sobolev) | Saloff_Coste_2004_Random_Walks_Finite_Groups |
| H | Varopoulos-Saloff-Coste-Coulhon heat-kernel | **UNVERIFIABLE-PHASE-0** (book paywalled, not in corpus) |

Additional candidates that surface during Phase 0 reading get appended.

---

## Pre-registered priors per candidate

(Phase 1 disposition probabilities: SELECTED / PARTIAL / NO_FIT / BLOCKER. Each row sums to 100.)

| Code | SELECTED | PARTIAL | NO_FIT | BLOCKER | Rationale |
|---|---|---|---|---|---|
| A | 3 | 12 | 60 | 25 | Chain has no obvious matrix-product representation; (Z/3^n)* is abelian; Lyapunov = trivial structure for abelian. Almost certainly NO_FIT on h_walk (matrix products) and h_conclusion (positive exponent). |
| B | 8 | 35 | 35 | 22 | Renewal structure of inputs (1)+(2) is non-standard (2D walk over partial sums); G-L-P targets stationary measures of random walks on R^d / homogeneous spaces. Profinite extension uncertain. |
| C | 12 | 30 | 28 | 30 | **Pre-cleared categorically** (BFLM targets p-adic stationary measures with quantitative Fourier decay) — main risk is irreducibility/stiffness and proximality hypotheses; Syracuse chain has no obvious proximality structure. |
| D | 5 | 25 | 45 | 25 | Homogeneous-space framing requires group action on quotient G/H; Syracuse's (Z/3^n)* with self-multiplication is degenerate as homogeneous space (G acts on itself = trivial). Stationarity yes; quantitative rate is the open Q. |
| E | 15 | 30 | 30 | 25 | **Pre-cleared categorically** in BT_DISPOSITION as live secondary route. Bourgain sum-product on multiplicative subgroups (2^v mod 3^n) is the closest in spirit. Main risk: hypotheses require multiplicative-group orbit, not tuple-space; Syracuse phase is over Geom(2)^n tuple. |
| F | 10 | 35 | 35 | 20 | Local CLT for products gives Fourier decay of n-step product distribution. Profinite extension and additive-vs-multiplicative bookkeeping are open. |
| G | 2 | 10 | 70 | 18 | Mixing-of-finite-groups is the wrong asymptotic (transient vs stationary); same problem as Probe R2 drift. Saloff-Coste 2004 specifically is finite-group, not profinite-limit. Almost certainly NO_FIT. |
| H | — | — | — | 100 | UNVERIFIABLE-PHASE-0 (book paywalled). Mark BLOCKER unless secondary-citation reconstruction in Phase 0 yields a workable statement. |

Aggregate prior over SELECTED: ~12% (matches honest prior in brief, ~10-15%).
Aggregate prior over PARTIAL: ~28% (matches "most likely outcome" framing).
Aggregate prior over NO_FIT: ~38%.
Aggregate prior over BLOCKER: ~23%.

---

## Decision rules (locked)

For each candidate K:

- **SELECTED.** Phase 1 hypothesis x input matrix shows ALL hypotheses SATISFIED (or SATISFIED-with-routine-structural-argument flagged NEEDS_PROOF that the inputs make tractable), AND Phase 2 conclusion-shape delivers polynomial-in-A bound (or a bound straightforwardly convertible to polynomial-in-A by the chain-side data in (3)), AND Phase 3 extension to profinite is STRAIGHTFORWARD or candidate is profinite-native.
- **NO_FIT.** At least one hypothesis FAILED by verbatim reading of input (1), (2), (3), or (4). Multiple-hypothesis failure increases confidence but doesn't change the category.
- **PARTIAL.** Phase 1 passes but Phase 2 conclusion is qualitative-only / non-polynomial / mismatched output. Report what additional structure or additional theorem would close the gap.
- **BLOCKER.** Theorem statement UNVERIFIABLE in available corpus (Phase 0 incomplete). Specify what literature would unblock.
- **MODE_H_CIRCULAR.** A hypothesis turns out to require polynomial-in-A Fourier decay (or an obvious equivalent of the closure target) as INPUT. Flag immediately and move to next candidate. Do NOT spend Phase 1/2/3 effort.

### Locked thresholds for disposition labels

- "Polynomial-in-A" means: a bound of form |mu_n(xi)| <= C(A) * f(n, xi)^{-A} for any A > 0, with C(A) depending only on A (and possibly on chain parameters fixed by R75-R77), or equivalent. The "Tao bound" target is S_chi(n) <<_A n^{-A}.
- "Convertible to polynomial-in-A" means: the theorem delivers a quantitative bound where the rate parameter is in principle improvable by sharpening one of: (a) an irreducibility/spectral-gap input, (b) a moment hypothesis on the step distribution, (c) a density-of-states / regularity input on the stationary measure. The chain-side data (3) (R76 conservation, R77 T_diag) is allowed to supply (a), (b), or (c).
- "Mode H circular" means: a hypothesis would be discharged only by knowing |mu_n(xi)| <= n^{-A} (or its Sobolev / Besov / regularity equivalent) — i.e., the conclusion is hidden in the hypothesis.

---

## Pre-registered secondary routing (if NO_FIT or PARTIAL dominant)

If FG closes negative or partial, these are the routing options flagged in advance:

1. **BGT regular variation (Bingham-Goldie-Teugels).** Operates on the k=7 jump signature: if eps_k has regular variation with index alpha (not in 8 coefficients yet), de Haan / Karamata machinery gives asymptotic recovery. The two-regime k<7 vs k>=7 signature suggests slowly-varying correction. Suitability: ~moderate.

2. **Igusa local zeta / functional equation.** The (1+3)^u algebraic structure (R78 disambiguation, "4=1+3 in Z_3") points to Igusa local zeta function of a polynomial in Q_3 — Z(s; f) with f(u) = (1+3)^u or related. Functional equation and meromorphic continuation give a generating-series transform that may evade Mode H. Suitability: ~moderate, requires identifying f explicitly.

3. **Heat-kernel / Brownian motion on profinite tree.** Saloff-Coste-style heat-kernel bounds, but on the 3-adic tree T_3 (= the Bruhat-Tits tree for SL_2(Q_3)). Different from cluster BT (closed-negative) because we're not embedding chain into SL_2 dynamics; we're using the tree as a metric profinite group on which standard heat-kernel exists. Suitability: ~moderate, depends on whether tree heat-kernel respects (Z/3^n)* fiber.

4. **Adelic Mellin construction.** Restore the global / archimedean place via adelic packaging: the chain on (Z/3^n)* lifts to a chain on the finite-idele class group A^*_Q,f / Q^*. Mellin transform across all places gives a multiplicative analog of Fourier; Tate-style local-global factorization may apply. Suitability: ~low-to-moderate, technically demanding.

5. **NEW candidate surfaced during Phase 0.** Specifically reserved slot for any theorem family encountered during PDF reading that's not in the A-H list.

Selection of secondary route is made AT FG_DISPOSITION time based on which one best matches the gap surfaced.

---

## Procedure discipline locks

- Mode E: theorem hypotheses VERBATIM from PDF text extracted via pypdf to UTF-8 file. Do NOT inherit assumptions from prior project files (memory, STATE.md, dispositions). Re-read.
- Mode H awareness: filter on first contact. If hypothesis = "stationary measure has Fourier decay at polynomial rate" or any equivalent, flag MODE_H_CIRCULAR and skip Phase 1/2/3.
- For each candidate: write FG_K_HYPOTHESES.md (verbatim), FG_K_HYPOTHESIS_CHECK.md (input matrix), FG_K_CONCLUSION_SHAPE.md (if Phase 1 passes), FG_K_EXTENSION_CHECK.md (if non-profinite-native).
- Final: FG_DISPOSITION.md.

---

## Honest prior outcome distribution (locked)

- SELECTED: ~12%
- PARTIAL: ~28%
- NO_FIT: ~38%
- BLOCKER: ~22%

Matches brief's honest prior. Most-likely outcome: PARTIAL with a specific theorem identified whose profinite extension requires moderate technical work.

---

End pre-registration. Proceed to Phase 0.
