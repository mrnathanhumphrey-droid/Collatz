# Bookkeeping Phase 1 — Pre-Registration

**Date locked:** 2026-05-11
**Project:** C:/Collatz/ — Tao 2022 §7.2–7.4 effective C_A bookkeeping feasibility
**Continues:** R77.2 (result_77_2_nisoli_certification.md, Stage 2 outcome (δ))
**Status:** locked BEFORE any Phase 1a/b/c work; commit precedes results commit

---

## 1. Motivation

R77.2 stalled at Stage 2 outcome (δ): the Nisoli Lemma 2.9 chain closes for the order-3 companion T_3 with explicit M = sup_γ ‖R(z, T_3)‖ ≈ 800–1000 (γ = circle of radius 1/8 around λ_2 = 1/2), but the required ε_K input — translated from Tao 2022 Prop 1.17, |E e^{-2πi ξ Syrac(Z/3^n Z)/3^n}| ≪_A n^{-A} — is given in Tao 2022 only qualitatively (footnote, p. 12: "uniform in n and ξ, though as indicated we permit it to depend on A").

Making C_A effective means redoing the proof of Prop 1.17 (Tao §7.2–§7.4) with line-by-line constant bookkeeping. This pre-registration locks the rules under which we will decide whether that bookkeeping is **tractable for this project's resources before committing to Phase 2 (actual bookkeeping).**

## 2. Hypotheses (pre-registered)

**H_TRACTABLE:** Tao §7.2–7.4 bookkeeping is tractable with project resources. The constants can be tracked through to an effective C_A as a function of A, producing a usable number for Nisoli ε_K input.

**H_NULL:** Bookkeeping requires arithmetic-combinatorics / analytic-number-theory expertise the project does not have, OR the proof structure absorbs constants in ways that resist line-by-line tracking, OR the resulting C_A would be too loose to satisfy η = ε_K · M_3 < 1 at the verified K range.

**Prior:** Pre-registered favoring NULL. Override to FEASIBLE requires positive feasibility on three independent checks (constant map clean, tractability classification mostly TRIVIAL/MODERATE, looseness projection satisfies η < 1 at K ≤ 10).

## 3. Procedure — Phase 1 (feasibility only, NOT bookkeeping itself)

### 3.1 Phase 1a: Proof structure map
Produce `TAO_PROOF_CONSTANT_MAP.md` enumerating every constant introduced, absorbed, or made implicit in Tao §7.2–§7.4. For each entry record:
- (i) the estimate it appears in (specific lemma/equation #s from arXiv 1909.03562 / Forum of Math Pi version);
- (ii) explicit (numeric), named-but-unspecified ("absolute constant c > 0"), or absorbed in Vinogradov ≪;
- (iii) parameter dependence (A, α, ε, m, n, P, …).

### 3.2 Phase 1b: Tractability assessment
For each constant in 1a, classify in `TAO_BOOKKEEPING_TRACTABILITY.md` as:
- **TRIVIAL** — explicit or one-step from cited reference;
- **MODERATE** — careful reading, no new math;
- **HARD** — requires deeper analytic technique or auxiliary estimate;
- **BLOCKED** — needs arithmetic-combinatorics / analytic-number-theory expertise the project lacks.

### 3.3 Phase 1c: Looseness projection
Produce `TAO_CA_LOOSENESS_PROJECTION.md` projecting a conservative upper bound on C_A under pessimistic absorption (each ≪ constant at the largest plausible value implied by its derivation). Compare against the Nisoli requirement `ε_K · M_3 < 1` with M_3 ≈ 800–1000 (R77.2 §3.3) and `ε_K ≈ C_A · K^{-A}` (or the Plancherel-corrected `≈ C_A^2 · K^{-2A}` depending on how the norm bound on T − T_K propagates) at K = 6 and K = 10.

### 3.4 Decision rule (pick ONE, no hedging)

- **FEASIBLE:** complete map; most constants TRIVIAL/MODERATE; HARDs have identifiable extraction paths; projected C_A satisfies Nisoli η < 1 at K ≤ 10 even under looseness pessimism. Phase 2 justified.
- **TRACTABLE_BUT_LOOSE:** 1a/1b clean; but projected C_A wouldn't satisfy η < 1 at K ≤ 10. Phase 2 useful only with adjusted expectations.
- **BLOCKED_BY_EXPERTISE:** enough HARD/BLOCKED constants that project can't complete bookkeeping. Phase 2 not justified.
- **INFEASIBLE:** proof structure resists line-by-line tracking, OR projected C_A too loose at any K. Bookkeeping route closes.

## 4. Adversarial safeguards (binding)

- **A1 — OCR fidelity:** local plaintext `C:/Collatz/tao2022.txt` is OCR'd. Greek letters and many in-line symbols render approximately. For any load-bearing estimate, cross-reference arXiv 1909.03562 via WebFetch and flag OCR discrepancies in 1a.
- **A2 — Absorbed-constant honesty:** every Vinogradov ≪ has specific content. Do not treat ≪ as benign. Either extract the implicit constant or record the obstacle. Tao's convention (p. 13) that c can "vary from line to line" is itself a load-bearing source of looseness — log every line where this happens.
- **A3 — Expertise honesty:** if a step requires arithmetic-combinatorics / analytic-number-theory technique the project lacks, mark BLOCKED. Speculative bookkeeping is worse than reporting a blocker.
- **A4 — Deviation logging:** pre-reg expects NULL-favored. If Phase 1 lands positive (FEASIBLE), the disposition file must document why each adversarial check passes — not just assert FEASIBLE.

## 5. Constraints (binding)

- Phase 1 is **feasibility assessment only**, NOT bookkeeping itself. We do NOT extract specific C_A values in Phase 1.
- arXiv 1909.03562 is canonical reference. Citations use Tao 2022 (Forum of Math Pi 10, e12, 2022 — published version contains §7.4) lemma/equation #s.
- Disposition is ONE of the four categories listed in §3.4. No hedging, no "between FEASIBLE and TRACTABLE_BUT_LOOSE", no extra branches.

## 6. Deliverables

| # | File | Phase |
|---|---|---|
| 1 | `BOOKKEEPING_PHASE1_PRE_REGISTRATION.md` (this file) | pre-reg — committed FIRST |
| 2 | `TAO_PROOF_CONSTANT_MAP.md` | 1a |
| 3 | `TAO_BOOKKEEPING_TRACTABILITY.md` | 1b |
| 4 | `TAO_CA_LOOSENESS_PROJECTION.md` | 1c |
| 5 | `BOOKKEEPING_PHASE1_DISPOSITION.md` | disposition + summary |

Commit structure: pre-reg first (separate hash), then 1a/1b/1c/disposition together (second hash). Mirrors recent project pattern (PROP_TEST_PRE_REGISTRATION → PROP_TEST_RESULTS).

## 7. Closure

After commits land, return a report under 400 words containing: disposition category, one-paragraph rationale, both commit hashes, highest-difficulty constants encountered with specific lemma/eq citations, any A1 OCR discrepancies found, Phase 2 justified / conditionally justified / closed.

— End of pre-registration. Locked.
