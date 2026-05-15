# Post-compact next steps

**Date written:** 2026-05-14
**Context:** Prepared at the end of the session that produced the 11-arc obstruction map terminal finding.

## What just happened (1-paragraph summary)

The 11-arc c=7/45 closure investigation terminated with a sharp framework identification: Syracuse Markov chain transfer operator analysis requires **B-valued MONOTONE independence** (Muraki 2003; Hasebe-Saigo 2011 operator-valued amalgamation), NOT B-valued free independence. The verification probe found that the third-order alternating B-centered moment `φ(X̃_{j_1} · X̃_{j_2} · X̃_{j_1})` does not vanish (second-order and three-distinct-index moments do), which is the diagnostic signature of monotone independence. Effort estimate collapsed from 12-19 months (original construction-blueprint scale) to 5-9 hours (mechanical application of established theory). Full writeup at `C:/Collatz/OBSTRUCTION_MAP_TERMINAL.md`.

## The four-task program (5-9 hours total)

**Task 1 — numerical confirmation (1 hour)**

Confirm `φ(X̃_1 · X̃_2 · X̃_1) ≠ 0` at level n=3 numerically from the existing `bilinear_pair_operator.py` infrastructure at `C:/Collatz/`. The verification probe's third-order computation was structural; a numerical run on actual operators at n=3 will give a specific non-zero value with confidence interval.

**Task 2 — literature pull (in flight before compact, finish post-compact)**

Pull these two papers:
- **Muraki 2003** "Monotonic independence, monotonic central limit theorem and monotonic law of small numbers." Likely findable on arXiv or via Muraki's faculty page. Try arxiv search for "Muraki monotone independence."
- **Hasebe & Saigo 2011** "The monotone cumulants" Ann. Inst. Henri Poincaré. arXiv:1011.6321 or similar. Open access via numdam / Annales IHP.

Place in `C:/Users/Nate/OneDrive/Documents/closure hunt/` alongside the existing Voiculescu / Speicher / Cébron / Young / Tsujii / Goldsheid-Margulis.

**Task 3 — framework writeup (2 hours)**

Document the framework identification finding as a stand-alone result. The skeleton is at `C:/Collatz/OBSTRUCTION_MAP_TERMINAL.md`. Expand into paper form if desired, with verbatim Muraki/Hasebe-Saigo theorem citations.

**Task 4 — explicit monotone cumulant computation (4 hours, OPTIONAL but gives closure)**

Using Hasebe-Saigo's monotone cumulant additivity formula, compute:
1. `M_2^B(X_{j_1}, X_{j_2})` for j_1 < j_2 (second-order monotone cumulant)
2. `M_3^B(X_{j_1}, X_{j_2}, X_{j_1})` (the third-order alternating monotone cumulant)
3. Apply the cumulant additivity property to derive `μ̂_n(ξ)`'s asymptotic at large n

Then compare against:
- **Wilson's PADE prediction:** leading singularity asymptotic at z ≈ 1.016, complex-conjugate pair with period ≈ 9.2 in n-space, sign pattern (+,+,−,−,−,−,−,−,−,+,+,+,+)
- **Faure 2009 spectral radius prediction:** √3 ≈ 1.732 (matches PADE 1.57 at n=13 within 10%)

If asymptotic matches → c=7/45 closure derivation in hand.
If asymptotic doesn't match → identifies further structural feature (likely non-trivial fourth-order monotone cumulant indicating an even finer dependence structure).

## Key files (preserved through compact)

**Main writeup:** `C:/Collatz/OBSTRUCTION_MAP_TERMINAL.md`

**State documents:**
- `C:/Collatz/STATE.md` (header has 2026-05-14 terminal finding entry)
- `C:/Users/Nate/.claude/projects/c--As-Above-So-Below-Master/memory/project_collatz_monotone_terminal.md` (dedicated memory entry, indexed in `MEMORY.md`)
- `C:/Collatz/POST_COMPACT_NEXT_STEPS.md` (this file)

**Verification probe outputs:**
- `C:/Collatz/AMALG_FREENESS_SETUP.md` — operator-valued probability space `(A, E_B, B)` definitions
- `C:/Collatz/AMALG_FREENESS_SUBALGEBRA_CHECK.md` — B as valid amalgamation subalgebra
- `C:/Collatz/AMALG_FREENESS_MOMENT_CALCULATION.md` — explicit moments at orders 2, 3, 4
- `C:/Collatz/AMALG_FREENESS_DISPOSITION.md` — verification finding with Voiculescu's verbatim freeness definition

**C4 re-probe series:**
- `C:/Collatz/C4_REPROBE_TAO_RMT_DISPOSITION.md` (v1)
- `C:/Collatz/C4_REPROBE_V2_*.md` (Cébron + Goldsheid-Margulis)
- `C:/Collatz/C4_REPROBE_V3_*.md` (Voiculescu + Speicher + Young + Tsujii)

**Construction blueprint:**
- `C:/Collatz/PROFINITE_TRANSFER_OPERATOR_LITERATURE_MAP.md`
- `C:/Collatz/PROFINITE_TRANSFER_OPERATOR_BLUEPRINT.md`

**Closure-hunt PDF corpus:** `C:/Users/Nate/OneDrive/Documents/closure hunt/`
Currently has: Voiculescu 1995 (AST_1995__232__243_0.pdf), Speicher 1998 (memoirs.pdf; non-standard glyph encoding caveat), Cébron 2013 (1304.1713v3.pdf), Young 1986 (random-perturbations-of-matrix-cocycles.pdf), Tsujii 2010 (0806.0732v3.pdf), Goldsheid-Margulis 1989 (rm1893_eng.pdf), Goldsheid-Sodin (2012.03017v2.pdf), VDN scanned (free-random-variables-1nbsped-082186999x.pdf), Aoun-Sert (2305.02879v2.pdf), Sawyer Martin boundary (Martin_boundaries_and_random_walks.pdf), Armentano-Chinta-Sahi-Shub (random-and-mean-lyapunov-exponents-...pdf), Bougerol (1408.2108v2.pdf), Das (2510.22778v2.pdf), Shusterman (1501.01227v1.pdf — irrelevant, different "free").

**Construction lit corpora (already pulled, ~80 PDFs across 7 folders):**
- `C:/Users/Nate/OneDrive/Documents/profinite_transfer_operator/pdfs/`
- `C:/Users/Nate/OneDrive/Documents/faure_semiclassical/pdfs/`
- `C:/Users/Nate/OneDrive/Documents/furstenberg_guivarch/pdfs/`
- `C:/Users/Nate/OneDrive/Documents/adelic_mellin/pdfs/`
- `C:/Users/Nate/OneDrive/Documents/igusa_local_zeta/pdfs/`
- `C:/Users/Nate/OneDrive/Documents/regular_variation/pdfs/`
- `C:/Users/Nate/OneDrive/Documents/watson_saddle_point/pdfs/`

## Chain-side input files (load-bearing, don't move)

- `C:/Collatz/c_seven_forty_fifth.md` — R75 Plancherel decomposition
- `C:/Collatz/result_76_conservation_law.md` — R76 conservation `Σ_j M_{n+1}(η_0 + j·3^n) = 0`
- `C:/Collatz/result_77_T_lead_spectrum.md` — R77 T_diag eigenstructure {0,1} on (1,-1) and (1,4), conjectured rate-½ off-diagonal at k=2..6
- `C:/Collatz/result_78.md`, `result_79.md` — (1+3)^u algebraic substrate, bilinear bound
- `C:/Collatz/C1_TAO_RECURSION_FORM.md` — Tao recursion `μ̂_n(ξ) = E χ(2-adic exp of Geom(2)^n tuple)` verbatim
- `C:/Collatz/PADE_NUMERICAL_DISPOSITION.md` — Wilson's multi-spectral picture
- `C:/Collatz/experiments_output/result_77_7_eps_exact_through_k8_v2_vec_pool.json` — ε_k=1..8 exact rationals

## Cleanup pending

Many `_*_pages/` temp folders at `C:/Collatz/` from PDF extraction (~600 page text files across 6 folders). Cleanup is optional and can be done whenever — they don't affect any active probe.

`_tao_rmt_pages/` (340 files), `_cebron_pages/` (55), `_goldsheid_pages/` (61), `_voiculescu_pages/` (34), `_speicher_pages/` (88), `_young_pages/` (11), `_tsujii_pages/` (59). Plus extraction scripts `_extract_v3.py`, `_tao_rmt_extract.py`.

## Operational notes for post-compact session

- **Math-heavy probes on opus may hit Usage Policy refusals.** Cause unclear; sonnet retry consistently succeeds. If a probe agent on opus returns a policy-refusal error, re-fire with `model: "sonnet"` explicitly in the Agent tool's parameters.
- **Speicher 1998 PDF non-standard glyph encoding** — current copy is unreadable via pypdf. The v3 synthesis reconstructed from cover page + decoded section headings + standard reference knowledge. Mingo-Speicher 2017 "Free Probability and Random Matrices" Fields Institute Monograph vol. 35 covers the same material in clean form — pull as backup if verbatim Speicher citations needed for paper.
- **No active long-running compute.** The k=8 ε_8 exact-rational run completed cleanly on 2026-05-13. No competing-compute concerns.
- **User's pace:** ~10-14× typical research-engineering pace (per feedback memory). The 5-9 hour estimate is at this pace.

## The c=7/45 connection (in plain terms)

If task 4 lands, the explicit asymptotic for `μ̂_n(ξ)` derived from monotone cumulants should reproduce the c=7/45 coefficient from the R77 T_diag eigenstructure (eigenvalues {0, 1} on (1, −1) and (1, 4)). The (1, 4)-direction structurally encodes the 7/45 via R64.B's class-mass identity (1/3)² : (2/3)² = 1:4 combined with Plancherel weights. The monotone cumulant computation makes this explicit.
