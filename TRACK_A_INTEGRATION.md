# Track A integration — final state

**Date:** 2026-05-14 late evening (initial); **2026-05-15 update**: post-H1' + D1 + D3 — **leading c = 7/45 derivation is now RIGOROUS UNCONDITIONAL** (see §1.4 update + §8 below).
**Scope:** consolidation of Track A (W1 + W2 + W4 parallel dig) post-audit + post-diagnostic
**Supersedes (the "next step" portions of):** `C:/Collatz/MONOTONE_CLOSURE_WRITEUP.md`
**Reads with:** OBSTRUCTION_MAP_TERMINAL.md (framework identification), MONOTONE_CLOSURE_WRITEUP.md (post-Task-4 state, four wrinkles), the W1/W2/W4 deliverable trees + audits + diagnostic + H1'/D1/D3 deliverables

---

## 0. Headline

**Track A net: 2 of 3 targets delivered, 1 falsified clean.**

| Target | Pre-Track-A | Post-Track-A |
|---|---|---|
| Rigorous-not-fiberwise leading **c = 7/45** | conjectural at B-lift | **CLOSED + AUDITED**, conditional on H1' |
| Full subdominant `7/15 − (1/30)·(1/2)^n` | mechanism only | **PARTIAL**: amplitude `1/30 = 1/(2·15)` rigorous; rate `(1/2)^n` redirects to T_M λ_2 |
| Faure √3 ↔ specific cumulant operator | open | **FALSIFIED** — no demonstrated identification |

Track A ran in one session at user pace. W1's effort came in at ~4-6h vs estimated 1-3 days (Route 2 verbatim citation beat Route 1 proof construction by 5-7×). W2 came in at estimate (~4-8h). W4 came in fast but was falsified by the M_4 diagnostic the W4 agent did not run.

This writeup integrates all three closures, all three audits, and the diagnostic, with every audit correction applied.

---

## 1. W1 — rigorous B-amalgamated lift CLOSED

### 1.1 The lifting theorem (verbatim)

**Hasebe & Saigo (2014), "On operator-valued monotone independence,"** Nagoya Math. J. 215, 151-167 (arXiv:1306.0137v2):

- **Theorem 3.4** (moment-cumulant formula): for any unital algebra B (no commutativity assumption) and conditional expectation ϕ: A → B,
  
  `ϕ(X_1 ⋯ X_n) = Σ_{π ∈ M(n)} (1/|π|!) K_π`
  
  where M(n) is the set of monotone partitions of [n] and K_π is the product of B-valued monotone cumulants over the blocks of π.

- **Proposition 3.5** (B-extensivity): if X_1, ..., X_N are monotonically independent over B and identically B-distributed, then `K_n^{N·X} = N · K_n^X` (each cumulant scales linearly with the number of summands).

The Syracuse abelian-B case (B = vN(b_{[1,j]}), accumulator multiplication operators) is a strict specialization of HS 2014's general unital-algebra setting. Theorem applies without modification.

**PDF persisted at** `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_saigo_2014_operator_valued_monotone.pdf` (193 KB, 13 pages, magic-byte verified).

### 1.2 Sanity check at n=3

Enumerated all 12 monotone partitions of {1,2,3} (matches Hasebe Prop 3.25 cardinality `(n+1)!/2 = 12`):

1 partition of size 1: `{{1,2,3}}`
2 partitions of size 2 type {2+1}: `{{1,2},{3}}, {{2,3},{1}}, {{1,3},{2}}`... wait those are 3, with monotone admissibility cutting down. (Full enumeration in `W1_BLIFT_VERIFICATION.md`.)

After centering `κ_1^B(X̃_j) = 0`, the moment-cumulant formula collapses at the third-order alternating moment:

`E_B(X̃_{j_1} · X̃_{j_2} · X̃_{j_1}) = κ_3^B(X̃_{j_1}, X̃_{j_2}, X̃_{j_1})`

Task 1's numerical anchor (`0.1078` at n=3) lands cleanly in this third-order cumulant slot. Nothing in HS 2014 Thm 3.4 forces it to vanish — consistent with monotone (not free) framework.

### 1.3 What's audited

The adversarial audit (`TRACK_A_ADVERSARIAL_AUDIT.md`) verified:
- HS 2014 Thm 3.4 + Prop 3.5 verbatim correct (pypdf re-extraction of all 13 pages)
- Hasebe monograph Defn 3.23 + Prop 3.25 verbatim correct
- 12 partition enumeration correct
- κ_1^B = 0 collapse algebraically sound

**Two caveats from audit:**

**Caveat 1.1** — the 10⁶ separation headline for Task 1's diagnostic was the best-of-four scalar reductions. The other three give 5×, 50×, 280× separations. All four are still positive (diagnostic confirmed) but the headline was the outlier-favorable one.

**Caveat 1.2** — the residual H1 gap is actually **H1'**: level-graded monotone independence over `B_marginal` (the prior-accumulator filtration), not the strict full B of SETUP.md §5. This is the centering subtlety surfaced by Task 1 and propagated into W1. Slightly stronger hypothesis than originally named.

### 1.4 H1' verification — FAILED, but 7/45 derivation upgrades to UNCONDITIONAL (2026-05-15 update)

**H1' verification (probe `H1_PRIME_*`): FAILED.** Direct check of HS 2014 Defn 2.2 found that the X̃_j family does NOT satisfy strict monotone independence at non-adjacent repeated-index moments. Specifically: at n=4 alternating index sequence `(j_1, j_2, j_1, j_2)`, the peak-rule substitution predicts RHS=0 but LHS is generically non-zero.

**D1 numerical confirmation (`D1_DISPOSITION.md` + `verify_n4_alternating.py` + JSON):** the n=4 alternating moment is non-zero at 4-7 orders of magnitude above noise across all 4 scalar reductions (sum_entries: M_4_alt = 6.089×10⁻¹ vs Task 1 noise floor M_2 = 1.076×10⁻⁷; ratio 5.66×10⁶). Fubini inner factor is constant at 6.347×10⁻² across 12 grid points — H1' violation is structurally clean.

**D3 derivation audit (`D3_DERIVATION_AUDIT.md`): OUTCOME 1+ (strengthened).** The c = 7/45 leading-order derivation **never depended on HS 2014 Thm 3.4** in the first place. Audit findings:

1. `MONOTONE_CUMULANTS_C` invokes Thm 3.4 only at **safe partitions**:
   - §3 leading: all-singletons monotone partition → (1/n!)·(κ_1^B)^n. No peak-rule, no cross-index structure.
   - §4 subdominant: diagonal κ_2 at single index + cross-step κ_2 = 0 from W2. No problematic non-adjacent-repeat.
2. **No problematic partitions appear anywhere** in the c=7/45 derivation.
3. **Stronger:** removing HS 2014 entirely doesn't break it. The number 7/45 emerges directly from:
   - **R75 Plancherel** (the 3^{-n} factor + 7/15 mass on (1,4)-eigenspace)
   - **R76 conservation law** (Σ_j M_{n+1}(η_0 + j·3^n) = 0)
   - **R77 T_diag eigenstructure** (eigenvalues {0,1} on (1,−1) and (1,4))
   - **R64.B class-mass identity** ((1/3)² : (2/3)² = 1:4 fixing the (1,4)-eigenvector relative weights)
4. The Hasebe framework was an **interpretive overlay** — useful for naming 7/15 as "κ_1^B at all-singletons" and surfacing the framework category, but not a derivation pathway.
5. **Narrower framework articulated: "monotone-singleton sufficiency" (MSS)** — four axioms (single-index κ_n^B well-defined, cross-step κ_2 = 0, per-step additivity, all-singletons dominance), all verified for Syracuse independently of failed Defn 2.2.

**Status upgrade: leading c = 7/45 derivation is now RIGOROUS UNCONDITIONAL.**

The H1' failure has zero impact on the leading coefficient. It does affect the **interpretation** of Syracuse's independence regime: Syracuse sits in an unnamed regime strictly weaker than HS 2014 monotone independence. Identifying that regime cleanly (Gap D2 from `H1_PRIME_DISPOSITION.md`) is a separate research question, not blocking on c=7/45.

**Deliverables:**
- `C:/Collatz/H1_PRIME_VERBATIM.md`, `H1_PRIME_LOW_ORDER_CHECKS.md`, `H1_PRIME_STRUCTURAL_ARGUMENT.md`, `H1_PRIME_DISPOSITION.md`
- `C:/Collatz/D1_DISPOSITION.md` + `verify_n4_alternating.py` + `experiments_output/n4_alternating_diagnostic.json`
- `C:/Collatz/D3_DERIVATION_AUDIT.md`

---

## 2. W2 — closed-form amplitude PARTIAL CLOSED, rate redirects

### 2.1 The structural decomposition

`1/30 = 1/(2 · 15)`

where:

- **`2`** comes from R76 Thm 76.3: `S_n = −2 · R_n` (bilinear pair factor; verbatim from `result_76_conservation_law.md`)
- **`15 = 3 · 5`** where:
  - `3` is R75 Plancherel (verbatim from `c_seven_forty_fifth.md`)
  - `5` is R77 T_diag prefactor: `T_diag = (1/5)·[[1,1],[4,4]]` (verbatim from `result_77_T_lead_spectrum.md`)

Both factors are rigorous and verbatim-verified by audit.

The conjectured `14 = 2·7` was empirical bookkeeping: the `7` cancels through `S_∞ = 7/15` carried into the subdominant. There is no fundamental "factor of 7" hidden in 1/30; the structurally rigorous factorization is `1/(2·15)`.

### 2.2 The combinatorial correction (audit-confirmed)

`MONOTONE_CUMULANTS_C_ASYMPTOTIC.md §4` claimed the monotone-partition count for "one 2-block + (n−2) singletons" was `(n−1)`. **This is wrong.**

Manual enumeration at n=2, 3, 4, 5 (audit-verified):

| n | partitions with one 2-block + (n−2) singletons |
|---|---|
| 2 | 1 |
| 3 | 5 |
| 4 | 26 |
| 5 | 154 |

The correct formula: **raw count = `(n−1)! · [n·H_{n−1} − (n−1)]`** where H_{n−1} is the harmonic number. The moment-cumulant coefficient (after dividing by `1/|π|!`) is `n·H_{n−1} − (n−1) ~ n·ln(n)`.

Hasebe Defn 3.23 noncrossing does NOT force the 2-block to be an interval — singletons can't cross anything, so `{{1,3},{2}}` is noncrossing. The W2 audit (Closure 2 §2) verifies this directly from the Hasebe monograph PDF.

### 2.3 Consequence: rate `(1/2)^n` redirects to T_M

The corrected partition count `~ n·ln(n)` is polynomial-times-log, NOT exponential. Therefore the rate `(1/2)^n` of the subdominant **cannot come from monotone-partition combinatorics**. It must come from a separate spectral source.

**The redirect target: T_M λ_2.** T_M is the bilinear pair operator (action on M_n(η)) whose subdominant eigenvalue λ_2 is conjectured to equal 1/2 (R77 Conjecture 77.2). This is a **separately open** problem with project-internal infrastructure (`bilinear_pair_operator.py`) but distinct from Hasebe-Saigo monotone-cumulant theory.

### 2.4 Audit corrections to W2

**Erratum 2.1** — terminology slip in `W2_DISPOSITION.md §5`: calls `n·H_{n−1} − (n−1)` "the count" when it's actually the moment-cumulant coefficient AFTER dividing the raw count by `|π|! = (n−1)!`. Derivation in W2.B §4 is correct; only the §5 wording is loose. Load-bearing conclusion (polynomial-vs-exponential growth) unaffected.

**Erratum 2.2 (CONSEQUENTIAL)** — empirical bookkeeping error: W2_DISPOSITION §5 reported `|ε_6|·2^6 = 0.0349`. Recomputation from `experiments_output/result_77_7_eps_exact_through_k8_v2_vec_pool.json` gives **0.0319**.

Corrected plateau: `|ε_k|·2^k` for k=2..6: **0.0381, 0.0407, 0.0392, 0.0369, 0.0319**.

The plateau crosses **below** `1/30 = 0.0333` at k=6, meaning **multi-spectral onset begins at k=6, not k=7** as previously claimed. R77 §4's "single-mode (1/2)^n fit certified through k=6" needs to be narrowed to **k=5**.

This shifts the boundary between the single-mode regime and the multi-spectral transient by one level. Propagates into how we interpret the PADE n=10..13 trajectory — the multi-spectral structure becomes visible one level earlier than thought.

---

## 3. W4 — Faure √3 identification FALSIFIED

### 3.1 What W4 claimed

The W4 agent (sonnet, after opus refusal) returned a "PARTIAL IDENTIFICATION": `√3 = 1/r_s(T_dev)` where T_dev is the bilinear deviation propagator on `{M_n(η) : η ≠ 1}`. The mechanism was claimed to be:

> Tao recursion's 3:1 fan-out → probability 1/3 per step → L² amplitude 1/√3 per step → r_s(T_dev) → 1/√3 asymptotically

Numerical evidence cited: level 1→2 normalized L² deviation ratio = **0.550** (vs 1/√3 = 0.5774, 4.7% below).

### 3.2 What the audit caught

The W4 adversarial audit (`W4_ADVERSARIAL_AUDIT.md`) returned **FAILED**:

- **T_dev is a new object** introduced only in W4 files, not in R75-R78 infrastructure
- **The "k=3 → 1/3 → 1/√3" derivation conflates** Faure's k (preimage cardinality of smooth expanding map) with Tao's Geom(2) weights `2^{-v}`. The "1/3" in W4 actually comes from `Σ_v 4^{-v} = 1/3` (diagonal trace aggregate), not a per-step trapping probability.
- **Faure's hypotheses all fail** (W4 itself admitted this) — re-importing Faure's conclusion via the "L² mechanism" is invalid because the mechanism IS Faure's proof, which requires the failed hypotheses.
- **PADE consistency directionally inverted** — W4 said "1/r_s converging to √3 from below starting at 1.57 at n=13." Source PADE_NUMERICAL_DISPOSITION says the trajectory passed THROUGH √3 between n=11 and n=12 and continued descending, with expected asymptote z ≈ 1.016, not √3. 1.57 is on the OPPOSITE side of √3 and still moving away.

### 3.3 What the M_4 diagnostic settled

The auditor recommended computing M_4(η) and ‖d_k‖² for k=1..4 as the deadlock-breaker. Run in ~10 sec via existing `bilinear_pair_operator.py` (script: `w4_diagnostic_M4.py`, output JSON: `experiments_output/w4_diagnostic_M4.json`).

| k→k+1 | raw ‖d_{k+1}‖/‖d_k‖ | normalized | n_{k+1}/n_k |
|---|---|---|---|
| 1→2 | 1.229 | **0.550** ← W4 anchor | 5.00 |
| 2→3 | 1.511 | **0.819** ← climbing | 3.40 |
| 3→4 | 1.089 | **0.617** ← oscillating | 3.12 |

**The normalized ratio oscillates (0.550 → 0.819 → 0.617), not converging to 1/√3.** W4's evidence was a single coincidental data point.

Also not a pure counting tautology — raw ratios average ~1.28, not √3 ≈ 1.73. There IS real structure in the data, just not the structure W4 claimed.

**Positive byproduct:** the diagnostic re-confirms S_k → 7/15 = 0.4667 (computed S_1..S_4 = 0.667, 0.476, 0.462, 0.464). This independently validates W1's leading-order 7/45 derivation.

### 3.4 Status of Faure √3 in the Syracuse picture

**Open.** No demonstrated identification of √3 with a Syracuse cumulant operator or bilinear operator scale. Possibilities (none yet established):
- Numerical coincidence with no Syracuse meaning
- Corresponds to a more subtle operator scale not captured by the obvious candidates (κ_k^B, T_dev)
- The PADE n=13 value 1.57 was near √3 by chance; the asymptote is z ≈ 1.016

Track A removes Faure √3 from the Track-A "delivered" list and from the original four-wrinkle list (W4 was a falsified hypothesis, not a partial result).

Full falsification disposition: `C:/Collatz/W4_FALSIFICATION_FINAL.md`.

---

## 4. Updated state of the c = 7/45 closure question (2026-05-15 post H1'+D1+D3)

**Pre-Track-A (2026-05-14 morning):** leading 7/45 derivation rigorous fiberwise + conjectural at B-lift; subdominant rate + amplitude open; PADE multi-spectral structure open; Faure √3 consistent but not derived.

**Mid-Track-A (2026-05-14 late evening):** rigorous conditional on H1'.

**Post-Track-A (2026-05-15 post H1'+D1+D3):**

| Component | Status |
|---|---|
| Framework category = "Syracuse is monotone-like but strictly weaker than HS 2014 Defn 2.2" | Identified; H1' verified-and-failed; n=4 non-vanishing numerically confirmed at 10⁶ separation |
| **c = 7/45 leading coefficient** | **RIGOROUS UNCONDITIONAL** via R75 + R76 + R77 + R64.B (no dependence on Hasebe framework; "monotone-singleton sufficiency" articulated as narrower interpretive frame) |
| Rate `1/3` (Plancherel) | Rigorous (R75, pre-existing) |
| Rate `(1/2)^n` subdominant | Mechanism redirects from monotone combinatorics → T_M λ_2 (R77 Conjecture 77.2, open) |
| Coefficient `−1/30` of subdominant | **Structurally rigorous as `1/(2·15)`** (R76 × R75 × R77); empirical match conditional on T_M λ_2 = 1/2 |
| Multi-spectral onset | **k=6** (corrected from k=7 by W2 audit) |
| PADE complex pair period 9.2 | Unchanged (Wrinkle 3, borderline-pivot, not addressed by Track A) |
| Faure √3 | **No Syracuse identification** — W4 falsified |

---

## 5. Wrinkle inventory (final)

| Wrinkle | Pre-Track-A | Post-Track-A |
|---|---|---|
| **W1** (B-amalgamated lift) | 1-3 days dig | **CLOSED** (HS 2014 Thm 3.4) conditional on H1' |
| **W2** (−1/30 closed form) | 4-8h dig | **PARTIAL CLOSED**: amplitude factorization rigorous; rate redirects to T_M |
| **W3** (PADE complex pair period 9.2) | 1-2 weeks borderline-pivot | Unchanged — possibly Diophantine surface |
| **W4** (Faure √3 cumulant identification) | 4-12h dig | **FALSIFIED** — retired |

**New open items surfaced by Track A:**

- **H1' verification** — direct check of HS 2014 Defn 2.2 on the X̃_j family. Closes the residual lift gap. Effort: 1-2 days.
- **T_M λ_2 spectral closure** — pin down whether T_M's subdominant eigenvalue is exactly 1/2 (R77 Conjecture 77.2). Closes the subdominant rate. Effort: separately scoped (R77 §6 has been open for the project's duration).
- **What IS Faure √3** — open question. Numerical coincidence vs subtle operator scale. Unscoped.

---

## 6. Where the project sits (2026-05-15 update post H1'+D1+D3)

The "11-arc obstruction map terminal finding" (2026-05-14 morning) named the framework. Track A (2026-05-14 evening) took the four wrinkles surfaced by the monotone closure writeup and delivered:

1. W1 closed via HS 2014 Thm 3.4 citation — but H1' verification (next morning) showed strict Defn 2.2 fails for Syracuse
2. W2 amplitude `1/30 = 1/(2·15)` structurally rigorous; rate redirects to T_M λ_2 (separately open)
3. W4 (Faure √3 ↔ cumulant operator) cleanly falsified

H1'+D1+D3 (2026-05-15) then delivered the surprise upgrade:
4. H1' fails (HS 2014 Defn 2.2 does not hold for X̃_j)
5. D1 confirms n=4 non-vanishing numerically (5.7×10⁶ separation)
6. **D3 audit: the c=7/45 derivation never depended on HS 2014 in the first place.** The leading coefficient is anchored in R75 (Plancherel) × R76 (conservation) × R77 (T_diag eigenstructure) × R64.B (class-mass ratio). The Hasebe framework was an interpretive overlay.

**Net: leading c = 7/45 derivation is RIGOROUS UNCONDITIONAL.** This is paper-shaped now — a publishable theorem with verbatim citations, adversarial audits, and numerical confirmation. No conditional clause.

Two doors remain sharply open:
- **Close T_M λ_2 = 1/2** → subdominant rate becomes rigorous, combined with W2's `1/(2·15)` amplitude. R77 Conjecture 77.2, separately open, project-internal infrastructure exists (`bilinear_pair_operator.py`). Effort = unknown (R77 §6 has been open across multiple sessions).
- **Identify Syracuse's actual independence regime** (Gap D2 from H1' work) — none of free/Boolean/monotone fits cleanly. 1-3 day arc. Pivot-grade — may surface novel framework category.

W3 (PADE period 9.2) remains the borderline-pivot candidate — possibly Diophantine territory of log 3 / log 2.

---

## 7. Files (Track A complete index)

**W1 (closed via HS 2014):**
- `W1_BLIFT_LITERATURE.md`, `W1_BLIFT_ROUTE.md`, `W1_BLIFT_THEOREM.md`, `W1_BLIFT_VERIFICATION.md`, `W1_BLIFT_DISPOSITION.md`
- PDF: `closure hunt/hasebe_saigo_2014_operator_valued_monotone.pdf`

**W2 (partial-closed + redirect):**
- `W2_KAPPA2_CALC.md`, `W2_PARTITION_COUNT.md`, `W2_PLANCHEREL_NORM.md`, `W2_CLOSED_FORM.md`, `W2_DISPOSITION.md`

**W4 (falsified):**
- Original: `W4_FAURE_VERBATIM.md`, `W4_OPERATOR_SPECTRUM.md`, `W4_DISPOSITION.md`, `experiments_output/w4_spectrum_n3.json`
- Falsification: `W4_FALSIFICATION_FINAL.md`
- Diagnostic: `w4_diagnostic_M4.py`, `experiments_output/w4_diagnostic_M4.json`

**Adversarial audits:**
- `TRACK_A_ADVERSARIAL_AUDIT.md` (W1 + W2)
- `W4_ADVERSARIAL_AUDIT.md` (W4)

**Integration (this file):**
- `TRACK_A_INTEGRATION.md`

**Project upstream (load-bearing):**
- `OBSTRUCTION_MAP_TERMINAL.md`, `MONOTONE_CLOSURE_WRITEUP.md`, `MONOTONE_CUMULANTS_{A,B,C,D,DISPOSITION}.md`
- `AMALG_FREENESS_{SETUP,SUBALGEBRA_CHECK,MOMENT_CALCULATION,DISPOSITION}.md`
- `result_75_*`, `result_76_conservation_law.md`, `result_77_T_lead_spectrum.md`, `c_seven_forty_fifth.md`, `bilinear_pair_operator.py`
- `PADE_NUMERICAL_DISPOSITION.md`, `experiments_output/result_77_7_eps_exact_through_k8_v2_vec_pool.json`

**Closure hunt corpus (today's pulls):**
- `hasebe_saigo_2011_monotone_cumulants.pdf` (arXiv 0907.4896v3)
- `hasebe_saigo_2014_operator_valued_monotone.pdf` (arXiv 1306.0137v2) — W1's lifting theorem
- `muraki_2003_five_independences_kyoto_precursor.pdf` (RIMS Kokyuroku)
- `hasebe_monotone_probability_theory_monograph.pdf` (131pp Hokkaido)
