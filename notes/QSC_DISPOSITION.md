# QSC_DISPOSITION — Hudson-Parthasarathy / Attal-Pautrat QSC vs Syracuse

**Date:** 2026-05-15
**Mode:** E. Self-adversarial. Honest verdict per brief's three-outcome menu.
**Companion:** `QSC_VERBATIM.md`, `QSC_SYRACUSE_IDENTIFICATION.md`, `QSC_MOMENT_PREDICTIONS.md`.

---

## 0. One-sentence verdict

**Outcome C — DOES NOT FIT.** Quantum stochastic calculus (Hudson-Parthasarathy 1984, Attal-Pautrat 2006) does NOT capture Syracuse's independence structure. The mismatch is not "near fit with a small extension" (Outcome B). It is structural: QSC's time-local Itô table + vacuum-centered increments + non-commutative filtration is geometrically opposite to Syracuse's classical-filtration + cross-time-coupled transfer operators. Syracuse genuinely needs a new framework.

The brief's prior expectation ("Attal-Pautrat 2006 is the obvious fit") is not borne out.

---

## 1. Verdict in detail

### Key candidates evaluated

(a) **HP 1984 — continuous-time QSC on Boson Fock space.** Three differentials `dA, dA†, dΛ` + time `dt`, with a 4×4 Itô table. Adapted processes via `Γ = Γ_{(t]} ⊗ Γ_{(t}` tensor splitting.

(b) **Attal-Pautrat 2006 — discrete-time toy QSC on atom-chain `⊗_k ℂ^{n+1}`.** Four matrix-unit basis operators `a_k^{ij}` per site, satisfying a discrete Itô table. Filtration `F_k = vN({a_j^{αβ} : j < k})`.

(c) **Köstler-Speicher 2008 — free with amalgamation over tail algebra.** Eliminated upstream (Syracuse fails free axiom at third cumulant).

### Why HP and AP both fail

1. **Filtration mismatch.** Syracuse's filtration `B_j = vN({M_{b_{[1, k]}} : k ≤ j})` is ABELIAN (built from real-valued running 2-adic-valuation sums on a classical Geom(2) sample space). HP's `Γ_{(t]}` and AP's `F_k` are NON-COMMUTATIVE (Fock or qubit-chain operator algebras). Syracuse's "past" is a classical probability space; QSC's "past" is a quantum operator algebra. The two filtration geometries are categorically different.

2. **Increment mismatch.** Syracuse's X̃_j is a single integral operator on `L²((Z/3^n)*)` derived from the (3-adic, 2-adic) arithmetic structure with Σ_{v ≠ v'} 2^{-v-v'}-weighted Geom(2) sum, χ-phase, and shift content. It has no natural decomposition into the HP triple `(dA, dA†, dΛ)` or the AP matrix-unit basis. The phase factor χ_j has `b_{[1, j-1]}` coupled multiplicatively to (2^{-v} − 2^{-v'}) INSIDE the exponential, so X̃_j does NOT factor as `(B-measurable scalar) · (pure quantum increment)`.

3. **Itô-table prediction mismatch.** Under every plausible identification (X̃_j ↔ dA + dA†, X̃_j ↔ dΛ, X̃_j ↔ a_j^{ij}, X̃_j ↔ F_j · dQ_j adapted integrand), the QSC Itô table predicts `ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1}) = 0` for j_1 ≠ j_2 (because the middle factor is vacuum-centered at a different time slice from the surrounding factors, and time-local independence kills the moment).

   Syracuse measures `ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1}) = 0.108`. The mismatch is 4–7 orders of magnitude above noise, robust across 4 scalar reductions, and structurally fundamental.

4. **The constant Fubini factor F(v_1, v_1') = 6.347×10⁻²** (constant across all 12 grid points) tightly constrains any framework: the inner factor `E_{(v_2)}[X̃_{j_2} X̃_{j_1} X̃_{j_2}]` is INDEPENDENT of the j_1-fiber data. QSC's adapted-process predictions would have the inner factor depend on the j_1-fiber via X̃_{j_1} content. Syracuse's constancy is **tighter** than QSC delivers, not looser.

5. **Categorical structure of non-commutativity.** QSC localizes non-commutativity AT each time-slice (the Fock-space site has non-commuting `dA, dA†`). Syracuse localizes non-commutativity in the per-step TRANSFER OPERATORS T_j on `H = L²((Z/3^n)*)`, with the filtration B remaining classical. The two geometries are **dual** in the wrong direction.

### Why no extension closes the gap

A "Syracuse-friendly QSC" would need:
- Abelian filtration B_j (matching Syracuse) — but then it's not "quantum" in HP/AP's sense.
- Non-time-local Itô products that couple distinct time-slices via B-measurable phases — but then it's not "Itô" in HP/AP's sense.
- Increments that mix B-measurable content INSIDE the noise integral (not as multiplicative integrand) — but then it's not "adapted process" in HP/AP's sense.

What's left after removing those three load-bearing pieces of HP/AP is **not QSC anymore**. It's a different framework that uses the words "filtration" and "increment" but has none of QSC's defining structure.

---

## 2. What Syracuse DOES need (positive characterization)

From the D1/D2/H1' findings, Syracuse's "right home" framework needs:

(P1) **Abelian filtration** indexed by ℕ.
(P2) **Single fixed operator per step** (no iid copies).
(P3) **Cross-step second moment = 0** (no κ_2 coupling distinct slices).
(P4) **Distinct-index higher moments = 0** (3 distinct indices ⇒ moment vanishes).
(P5) **Non-adjacent repeat moments NON-ZERO**, driven by level-graded phase twist Δ_{j_2}(b_{[1, j_1]}).
(P6) **B-measurable phase content mixed INSIDE the noise integral** (not factored as integrand × differential).
(P7) **Constant Fubini inner factor** — the inner-pair integration produces an effective scalar that's independent of the outer fiber.

(P1)–(P4) overlap with QSC structurally. (P5)–(P7) are where Syracuse departs.

### Closest analogs in the literature

| Property | Closest QSC analog | Syracuse needs |
|---|---|---|
| Abelian past filtration | None in HP/AP (their pasts are non-commutative) | Yes |
| Adapted increment | HP `F_t dA_t` with F_t ∈ Γ_{(t]} | Syracuse has B_{j-1}-coefficient INSIDE the increment, not multiplying it |
| Cross-time coupling at non-adjacent repeats | Forbidden by HP/AP Itô table (only time-local products) | Required by Syracuse measurement |

The picture: **Syracuse is closer to a CLASSICAL Markov / adapted-process setting on the abelian B-filtration**, but with non-commutative T_j-content carried separately on the H-side. The cross-time coupling P5/P6 is a **classical conditional-expectation phenomenon**, not a quantum Itô one — but the T_j's are genuinely non-commutative, so it's also not purely classical.

This suggests the right framework may be:

(i) A **classical-quantum hybrid**: classical (abelian) filtration B carrying the phase data, with quantum (non-commutative) operators T_j coupled via a conditional-expectation tower.

(ii) Or a **"quantum Markov chain" on a commutative state space** — these exist in the literature (e.g. Accardi-Frigerio-Lewis 1982 "Quantum Markov processes"), but typically have the same time-locality property that HP/AP have. The Accardi framework is worth checking but is unlikely to capture P5 cleanly.

(iii) Or genuinely a NEW framework that hasn't been published — the project's `OBSTRUCTION_MAP_TERMINAL.md` per the MEMORY entry points to "B-valued MONOTONE independence (Muraki 2003 / Hasebe-Saigo 2011)" as the load-bearing candidate; H1' verification showed Defn 2.2 verbatim doesn't hold, so even monotone-B isn't the literal fit. A **level-graded refinement of monotone-B** with the operator-valued conditional-expectation tower from §1.2 of `H1_PRIME_STRUCTURAL_ARGUMENT.md` may be the right home.

---

## 3. Implication for the c = 7/45 derivation

**No change.** Per `D3_DERIVATION_AUDIT.md` and `THEOREM_C_745.md`, the leading c = 7/45 coefficient is rigorous UNCONDITIONAL via R75 + R76 + R77 + R64.B + HR74. No QSC framework is needed for c = 7/45.

The framework question is about higher-order corrections, the subdominant rate, and the SHAPE of the asymptotic expansion — none of which depend on c = 7/45 alone.

---

## 4. Mode-E gaps remaining

| Gap | Description | Effort |
|-----|------|----|
| QSC-G1 | Verbatim equation-numbered quotes from HP 1984 / AP 2006 PDFs pending user download (sandbox blocked). Verbatim quotes refine the §1 description in QSC_VERBATIM.md but do NOT change the disposition. | 30 min user fetch |
| QSC-G2 | A direct moment computation for the AP framework with site-Hilbert-space identification `n + 1 = 3^n` (forcing Syracuse Z/3^n → AP site basis) to see if row (d) becomes non-zero under that bespoke choice. Predicted to remain 0 by general time-local argument, but a clean Mode-E confirmation would close the loop. | 4-6 hours scripting |
| QSC-G3 | Accardi-Frigerio-Lewis 1982 "Quantum Markov chains" framework check — does that have the right "classical past + quantum increments" structure? Not done in this session. | 1-2 day lit pull |
| QSC-G4 | The "level-graded refinement of monotone-B" hypothesis — does Muraki 2003 or Hasebe-Saigo 2011 admit a refinement with an operator-valued conditional-expectation tower E_{B_j} that satisfies a peak rule at the prior-accumulator level (vs the verbatim B_∞ peak rule which fails)? This is what `H1_PRIME_STRUCTURAL_ARGUMENT.md §1.4-1.7` argued for but `H1_PRIME_DISPOSITION.md` ultimately rejected via D1's confirmation. A refinement may exist with weaker axioms. | 1-3 days |

---

## 5. Recommended next step (not done in this session)

Switch the closure hunt to:

(a) **Accardi-Frigerio-Lewis 1982 "Quantum Markov Processes"** (and follow-ups by Accardi, Kümmerer, Maassen) — these provide quantum Markov chains over commutative state spaces with the structural geometry Syracuse needs (classical past, quantum increments).

(b) **Sauvageot 1986 / Bhat-Parthasarathy 1995 — Markov semigroups on operator algebras** — the operator-algebra Markov property may match P1-P7.

(c) **Belavkin 1989 / Belavkin-Staszewski quantum filtering** — the brief mentioned this; specifically, the quantum filtering equations have an abelian observation filtration coupled to non-commutative system operators, which matches Syracuse's P1 + non-commutative T_j separation. **This may be the closest published structural fit.**

The Belavkin filtering candidate is the most promising next pull. It's a NATURAL extension of HP 1984 that explicitly handles the case where one wants a CLASSICAL observation filtration (abelian) feeding back into non-commutative system dynamics — exactly Syracuse's P1+P5+P6 combination.

---

## 6. Files

- This disposition: `C:/Collatz/QSC_DISPOSITION.md`
- Verbatim defns: `C:/Collatz/QSC_VERBATIM.md`
- Identification attempt: `C:/Collatz/QSC_SYRACUSE_IDENTIFICATION.md`
- Moment predictions: `C:/Collatz/QSC_MOMENT_PREDICTIONS.md`
- Project context: `H1_PRIME_DISPOSITION.md`, `D1_DISPOSITION.md`, `D2_BMT_BIGRAPH_DISPOSITION.md`, `THEOREM_C_745.md`, `C1_TAO_RECURSION_FORM.md`, `AMALG_FREENESS_SETUP.md`
- Open PDF URLs (user-fetch follow-up): see QSC_VERBATIM.md §0
