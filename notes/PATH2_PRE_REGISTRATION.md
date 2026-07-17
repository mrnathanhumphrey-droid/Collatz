# PATH2 Pre-Registration — Direct Construction at Family Level

**Timestamp:** 2026-05-11 (UTC date locked at session start).
**Author:** subagent under user-direct-construction brief.
**Status:** LOCKED before any computation. No retroactive edits.

## Strategic context

Five lit extraction attempts ruled out the "translate an existing theorem" path for c=7/45 closure (GARCIA_YOUNG_DISPOSITION.md, PASCADI_DISPOSITION.md). Remaining structural path is **direct construction** via project's R78.4-78.6 q=3 machinery extended to family prime p ≥ 3.

## Pre-registered hypotheses

- **H_DIRECT_WORKS** — Family-level extension is mechanical (Cochrane log generalizes p-blindly + saddle structurally same + bijection extends); substitution into bilinear sum reveals phase-cancellation mechanism giving rigorous |Σ 1̂·F̂| ≤ C·N·√q with explicit C matching empirical (≈ 2.0 up to polylog).
- **H_PARTIAL** — Extension works; bilinear bound has explicit log factor / suboptimal constant / restricted r range.
- **H_NEEDS_NEW_MATH** — Extension works; bilinear substitution doesn't admit elementary bound. Requires new analytical input (specified explicitly).
- **H_CLOSES** — Closed-form approach reveals fundamental obstruction.

**Priors:** H_PARTIAL or H_NEEDS_NEW_MATH favored. H_DIRECT_WORKS unlikely given R79b's β=0.522 empirical regime with no sub-Weyl saving. H_CLOSES informative.

## Locked procedure

### Phase 1 — Family extension of R78.4-78.6

Replace `3 → p` throughout R78.4-78.6 statements. Identify q=3 specifics. Define:
- `L_p(1+ps) = Σ_{j=1}^{J_p} (-1)^{j-1}/j · (ps)^j` (Cochrane truncated p-adic log)
- `J_p` = max j with j - v_p(j) < r+1
- `L̃_p = L_p(1+p)/p` (unit after stripping single factor of p)
- `C_a = a · L̃_p^{-1} mod p^r`
- `P_a(s) = ps − C_a · L_p(1+ps)` in coefficient-form mod q = p^{r+1}
- Saddle: `s*(C_a) = (C_a − 1)/p mod p`
- `F̂_p(p·a) = p · e_q(c) · G_p(a)` for a ≡ c mod p

Deliverable: PATH2_FAMILY_EXTENSION.md.

### Phase 2 — Empirical verification of extension

Cells: `(p, r) ∈ {(3, 2), (3, 3), (5, 2), (5, 3), (7, 2), (7, 3), (11, 2), (11, 3)}`.

Per cell, verify (all to ≤ 1e-9 absolute):
- C1: bijection a ↔ C_a holds (computational exhaustive check; collision count = 0).
- C2: |G_p(a)| = p^{(r+1)/2} for all a in support (within float ~1e-12).
- C3: At r=2 with saddle prediction G_p(a) = √q · e_q(P_a(s*(C_a))) — magnitude equality `|RHS| = √q` exact; phase agreement at r=2 might pick up Gaussian-integration factor (analogous to q=3 r=2 picking up e^{iπ/6}); document precisely.
- C4: At r=3, saddle phase prediction should match exactly (per q=3 r=3 exact match in R78.6). Tabulate per a.

Deliverable: PATH2_FAMILY_EXTENSION_VERIFICATION.md + CSV.

**Stop rule:** If C1/C2 fails at any cell, disposition lands H_CLOSES with the cell-specific obstruction documented. Skip Phases 3-5.

### Phase 3 — Substitute into bilinear sum

Define
- S_p := Σ_{a ≡ c mod p in Z/p^r} 1̂(p·a) · F̂_p(p·a)
- T_p := Σ_{a ≡ c mod p in Z/p^r} 1̂(p·a) · e_q(P_a(s*(C_a)))

with `1̂(p·a) = Σ_{u=0}^{N-1} e_q(p·a·u)`, N = p^{r-1}.

Goal: bound |T_p| from explicit form. Attempts to try:
- A: substitute C_a = a · L̃_p^{-1}, get phase explicit in a.
- B: partition by s*-class (j=0,1,2,...,p-1), use R79b's piecewise-linear-in-a structure.
- C: Cauchy-Schwarz on |T_p|², separate diagonal/off-diagonal.
- D: Poisson on the a-sum if structure suggests.
- E: empirical computation of T_p at small cells, compare to empirical |K|.

Deliverable: PATH2_BILINEAR_FROM_CLOSED_FORM.md.

### Phase 4 — Adversarial checks

- A1: numerical T_p vs empirical |K|/√N saturation. Magnitude mismatch = red flag.
- A2: Hensel safety at r ≥ 4. Does the analysis depend on r ≤ 3 exact closed form?
- A3: Cubic-character speculation walkback (R79b). Phase 3 must NOT reproduce the falsified "ψ_lead is cubic in a" framing.
- A4: For H_NEEDS_NEW_MATH, name the specific new input (literature analog if any).

### Phase 5 — Disposition

Deliverable: PATH2_DISPOSITION.md with disposition at top + 1-paragraph rationale + arguments / scope / new input / obstruction.

## Decision rules (locked)

- **H_DIRECT_WORKS:** Phase 3 produces rigorous |T_p| ≤ C·N (hence |S_p| ≤ C·N·√q) with explicit C **AND** A1 magnitude match within factor 2 **AND** A2 survives Hensel at r ≥ 4 **AND** A3 not relying on walked-back claim.
- **H_PARTIAL:** Bound holds with explicit log factor OR at restricted r range OR with explicit asymptotic-only constant; scope documented precisely.
- **H_NEEDS_NEW_MATH:** Phase 3 attempts hit specific elementary-bound failure; SPECIFIC analytical input named (lemma type, estimate framework). Not vague "more work needed".
- **H_CLOSES:** Phase 1/2 reveals q=3 specific in the extension OR Phase 3 reveals fundamental structural obstruction (phase doesn't cancel / Hensel destroys structure / R79b's anomalous j=0 class breaks all bilinear attempts).

## Adherence rules

- Procedure followed as locked. Method deviations logged explicitly with reasoning.
- No retroactive hypothesis changes.
- No threshold relaxation.
- Honest scope: H_NEEDS_NEW_MATH requires specifying WHAT input; can't be a hedge.
- If Phase 1 or 2 fails, STOP at H_CLOSES — don't waste Phase 3 effort.

## Files committed to writing

1. PATH2_PRE_REGISTRATION.md (this file — WRITTEN FIRST)
2. PATH2_FAMILY_EXTENSION.md
3. PATH2_FAMILY_EXTENSION_VERIFICATION.md (+ CSV)
4. PATH2_BILINEAR_FROM_CLOSED_FORM.md
5. PATH2_DISPOSITION.md

No git commits. User handles commits.
