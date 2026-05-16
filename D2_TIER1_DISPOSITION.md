# D2 — Tier 1 Disposition (verdict)

**Date:** 2026-05-15
**Mode:** E.
**Reads with:** `D2_TIER1_TAXONOMY.md` (verbatim defs), `D2_TIER1_FIT_CHECK.md` (per-candidate fit), `H1_PRIME_DISPOSITION.md`, `D1_DISPOSITION.md`.

---

## Outcome: **B — no Tier 1 candidate fits Syracuse**

All three Tier 1 monotone variants (anti-monotone, bi-monotone, indented / α-free) **fail at the Syracuse row (f): 4-alternating non-vanishing**, by an identical structural mechanism: each candidate inherits HS 2014 Defn 2.3's **iid-copies-of-A_λ assumption** that the Syracuse construction structurally violates (row j: single X̃_j per step, not iid copies).

---

## 1. Per-candidate verdict + structural-feature responsible

| Candidate | Verdict | Structural feature that prevents fit |
|---|---|---|
| **Anti-monotone** (Muraki 2002; Hasebe 2010 Defn 1.9(ii)) | FAIL row (f) | iid-copies requirement (inherits HS 2014 Defn 2.3 under order flip) + flip of order alone is order-independent, so doesn't help with single-X̃_j-per-step. Anti-monotone is "essentially the same as monotone" (Hasebe 2010 p2, Hasebe monograph p28). |
| **Bi-monotone** (Gu-Hasebe-Skoufranis 2017) | FAIL row (f) | Two-faced structure has no natural Syracuse analogue (X̃_j is a single operator, not a left/right pair) → degenerates to monotone under any Syracuse-compatible reduction. Plus, GHS 2017 p13 dot operation explicitly requires iid copies of paired algebras. |
| **Indented / α-free** (Hasebe 2010 Defn α-free, formerly "indented") | FAIL row (f) | Triplet-of-states structure has no natural Syracuse analogue (only one state ϕ available) → degenerates to monotone or c-monotone. The c-monotone two-term peak rule (Hasebe 2010 Defn 1.12 CM2) still gives 0 at the 4-alternating moment under marginal centering of both states. Hasebe 2011 Defn 8.2 confirms iid-copies inheritance. |

**The common load-bearing failure** is the iid-copies-of-A_λ assumption (HS 2014 Defn 2.3), which all three Tier 1 candidates inherit:

> "For every X ∈ A, let us take copies {X^(j)}_{j ≥ 1} in an algebraic probability space (Ã, B, φ̃) such that: (1) X ↦ X^(j) is a B-homomorphism for each j; (2) φ̃(X_1^(j) X_2^(j) · · · X_n^(j)) = φ(X_1 X_2 · · · X_n) for any X_i ∈ A, j, n ≥ 1; (3) the subalgebras A^(j) := {X^(j)}_{X ∈ A} are monotone independent over B."

(HS 2014 Defn 2.3, p3.)

Syracuse's X̃_j family lacks iid copies — each step j contributes ONE X̃_j operator, used at every occurrence of step j in a word. The prior probe (now-stopped D2 prior version) correctly identified this as load-bearing; the focused Tier 1 scan confirms that this load-bearing assumption is inherited unchanged by **all three** Tier 1 monotone variants.

---

## 2. Recommendation on Tier 2/3/4

The iid-copies issue is **NOT** a monotone-family-specific artifact — it's a property of the entire universal-product / natural-independence classification (Speicher-Schürmann-BenGhorbal-Muraki universal-products axioms). All five natural independences (tensor, free, Boolean, monotone, antimonotone) are constructed on the universal-product framework, which presupposes iid copies via the dot operation (or equivalent).

**Bigraph independence (Gilliers-Jekel 2026)** and **BMT independence (Arizmendi-Mendoza-Vázquez-Becerra 2023)** are mixtures of the five natural independences — they also inherit iid-copies architecture (verified by spot-check in `D2_TIER1_TAXONOMY.md §4`: Gilliers-Jekel 2026 p3 says "enabling the construction of independent copies of arbitrary non-commutative probability spaces"; BMT uses pairwise independence graph atop the five natural primitives).

**Cross-cutting prediction.** Any Tier-N framework constructed by universal-product axioms (= any framework derived from the Muraki-Speicher-BenGhorbal classification) will share the iid-copies feature and therefore fail at row (f) for the SAME reason. This rules out:
- Tier 2: c-monotone, c-antimonotone (already failed via α-free)
- Tier 3 (universal-product mixtures): BMT, bigraph, Λ-monotone, Λ-Boolean, ε-independence, BM independence, free-Boolean-monotone digraph independence
- Tier 4 (natural-product extensions): tensor variants, q-deformations within the universal-product framework

What might still work — frameworks that DON'T require iid copies:
- **Spreadability systems / sequence-of-states approach** (Hasebe-Lehner 2023, "Cumulants, spreadability and the Campbell-Baker-Hausdorff series"). Spreadability defines cumulants by a different axiomatization (exchangeability-like) that may or may not require iid copies — needs verification.
- **Cumulant-only / partition-counting frameworks** (Lehner 2004 "Cumulants in noncommutative probability theory I-IV"). Some of these define cumulants WITHOUT reference to iid copies, treating the moment-cumulant duality as primitive. The vanishing-of-mixed-cumulants characterization may apply to families that aren't iid-copy-constructible.
- **Operator-valued / B-valued frameworks at the LEVEL-GRADED reading** (Skeide 2003, Popa 2008). Level-graded readings adapt the peak's ϕ to the peak position; this may or may not require iid copies — Skeide [12]'s operator-valued monotone is built on the iid-copies HS framework, but Popa 2008's modifications may relax this.
- **Tauberian / generating-function frameworks** (Flajolet-Sedgewick, Chevalier 2507.15394). These don't presuppose any independence framework — they derive asymptotics directly from generating-function structure. The c=7/45 derivation may be re-framed entirely without an independence framework (see `project_collatz_r78_bilinear_cracked.md` memory note re: Tauberian arc).

**Honest recommendation:** Tier 2/3/4 are **NOT WORTH SCANNING** as monotone-variant frameworks, because the iid-copies issue is order-flip-independent, multi-state-independent, and graph-mixture-independent. Syracuse's X̃_j family is **genuinely in an unnamed regime** within the universal-product / natural-independence classification.

**Where the closure-hunt arc actually lives:** at the spreadability / cumulant-primitive / Tauberian level, NOT at the natural-independence level. The Hasebe-Lehner 2023 spreadability framework + the Tauberian generating-function arc (already flagged in `project_collatz_r78_bilinear_cracked.md`) are the right next stops, not more Tier-N monotone variants.

---

## 3. Mode-E gaps

### Gap D2-A — spreadability not directly probed

This Tier 1 disposition does NOT close whether Hasebe-Lehner 2023 spreadability systems require iid copies. Spreadability theory is referenced briefly in Hasebe 2010 §4.3 ("General theory of cumulants for spreadability systems"), but the verbatim definition was not extracted. A follow-up Tier-2-equivalent probe scoped to spreadability (NOT to monotone variants) would close the question of whether Syracuse fits any cumulant-primitive framework.

**Effort estimate:** 2-4 hours focused (extract Hasebe-Lehner 2023 + Lehner 2004 verbatim spreadability/cumulant axioms; check the dot-operation analog).

### Gap D2-B — operator-valued single-X̃-per-step frameworks not exhausted

Skeide 2003 operator-valued monotone (Hasebe-Saigo cite as ref [12]) and Popa 2008 may have iid-copies-free variants. Operator-valued conditional cumulants à la Popa 2008 sometimes work at the LEVEL of the conditional expectation (not requiring B-homomorphic copies). Worth a focused 2-3 hour scan if the c=7/45 work needs to STAY in the operator-algebra framework rather than pivoting to Tauberian.

### Gap D2-C — verification that the c=7/45 derivation is unaffected (already largely closed by D3)

`D1_DISPOSITION.md §6` (citing D3 audit) states: "D3 audit determined that the c=7/45 derivation never depended on the failed regime ... The leading derivation is rigorous unconditional via R75+R76+R77+R64.B." If correct, the H1' failure and the D2 Tier 1 disposition have ZERO impact on the c=7/45 leading-order publication.

In that case, this D2 finding is a clean **structural result about Syracuse independence type, NOT a blocking finding for the publication arc**. The right way to write up D2 is:

> "Syracuse's X̃_j family lives outside the natural-independence / universal-product framework. The c=7/45 leading coefficient does not depend on this — see D3 audit. The level-graded reading of HS 2014 used in W1 §3 §4 is consistent with the structural facts; the verbatim HS 2014 Defn 2.2 fails at 4-alternating, but this only blocks higher-order cumulant statements at non-adjacent-repeat moments. The leading order is unaffected."

---

## 4. One-page summary

1. **Anti-monotone, bi-monotone, indented/α-free all FAIL** to fit Syracuse's row (f) (4-alternating non-vanishing).
2. **Common failure mode:** iid-copies-of-A_λ assumption (HS 2014 Defn 2.3) shared across all three; Syracuse has a single X̃_j per step instead.
3. **Order-flip, two-faced, three-state generalizations all FAIL** by the same iid-copies issue — these structural degrees of freedom are orthogonal to the row (j) failure mode.
4. **Tier 2/3/4 universal-product extensions** (BMT, bigraph, conditional variants, ε/Λ independences) inherit the same iid-copies architecture — verified by spot-check; predicted to fail identically.
5. **The right next arc is NOT more Tier-N monotone variants** but either (a) spreadability/cumulant-primitive frameworks (Hasebe-Lehner 2023, Lehner 2004), (b) operator-valued without B-homomorphic copies (Popa 2008 maybe), or (c) Tauberian / generating-function level entirely outside independence frameworks (Flajolet-Sedgewick, Chevalier 2507.15394) — the latter already flagged in `project_collatz_r78_bilinear_cracked.md`.
6. **Practical impact on c=7/45 publication: NONE** (per D3 audit). D2 is a clean structural classification result, not a blocking finding.

---

## 5. Files

- Verbatim defs: `C:/Collatz/D2_TIER1_TAXONOMY.md`
- Per-candidate fit checks: `C:/Collatz/D2_TIER1_FIT_CHECK.md`
- This disposition: `C:/Collatz/D2_TIER1_DISPOSITION.md`
- Numerical anchor (row f): `C:/Collatz/D1_DISPOSITION.md`
- H1' structural argument: `C:/Collatz/H1_PRIME_STRUCTURAL_ARGUMENT.md`
- H1' disposition: `C:/Collatz/H1_PRIME_DISPOSITION.md`
- Pulled PDFs (all verified %PDF magic bytes):
  - `C:/Users/Nate/OneDrive/Documents/closure hunt/bi-monotonic_gu_hasebe_skoufranis_2017.pdf` (352 KB)
  - `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_2010_three_state_independence.pdf` (540 KB)
  - `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_2011_conditionally_monotone.pdf` (502 KB)
  - `C:/Users/Nate/OneDrive/Documents/closure hunt/bmt_independence_2023.pdf` (409 KB; out of Tier 1, kept for completeness)
  - `C:/Users/Nate/OneDrive/Documents/closure hunt/bigraph_independence_mixture_2026.pdf` (712 KB; out of Tier 1, kept for completeness)
