# Result — swap-involution J REFUTED (walk-back #14); the extracted constraint; board state at close of Phase-2b turn.

**Date:** 2026-07-16. Nathan's turn-close: a second candidate symmetry (the swap involution J) refuted for the same reason as the rotation S, one genuine constraint extracted from the double failure, and the Phase-2b board consolidated. Banked as reported (the refutations and the L-A/L-B/skeleton items are Nathan's session results; provenance tagged below).

## J (swap involution) REFUTED — walk-back #14
`[M, J] = 0.26` at L=2 — **numerically identical to S's commutator** (both `(0.444, 0.258)`). The naive swap fails for the same reason the rotation did, and the matching norms say it is the **same breakage**: the carry's floor arithmetic `(g+T)//q` on representatives in `[0, q^L)` **does not commute with modular negation** any more than with modular rotation.
- ⇒ the **odd-k selection rule is NOT earned**. No fishing for twisted variants until one sticks — that is the decoration anti-pattern, and the record is unambiguous where it leads.
- The **odd-k zero-amplitude pattern** (from H1) joins the **R32 near-λ₁ tower cluster** as the **second banked open mystery** (selection-rule explanation refuted as posed).

## The genuine hint (a constraint, not a theory — filed as such)
> **Any candidate invariance of this operator must act on the carry as an INTEGER map, not a modular one.** The representative arithmetic is load-bearing; both refutations (S rotation, J swap) died at exactly that point.

For whoever hunts the true symmetry later: modular actions on γ are ruled out; the invariance (if any) lives in the integer/representative structure of the carry.

## Board state — close of Phase-2b turn (consolidated)

**PROVEN (Claude-verified + Nathan-written):**
- **THEOREM D1 complete** — first spectral gap of the program: `r(λ) = (1−λ²)/(1+λ²)`, maximality via nilpotence of the e=−1 block (`result_phase2b_Dmax.md`, `result_phase2b_F.md`).
- **Five substrate lemmas** (Phase 1): FORGET / ONE-STEP / INTERTWINE / REFINE / PYTHAGORAS.

**MACHINE-VERIFIED, awaiting proof:**
- **Full circulant family {c_k} exact at every level** — 2/6/**18** complete, LU-confirmed (`result_phase2b_H.md`). *[Claude-verified.]*
- **L-A no-return** (C ↛ Δ; one-step M[Δ,C]=0 + BFS at L=2,3). *[**Claude-verified 2026-07-16**, `result_phase2b_LALB.md`.]*
- **L-B gauge factorization** (k=0 co-invariant ℓ₀ factors through (ρ,γ) exactly, dev 9e-16; 54 orbits = ≤54). *[**Claude-verified 2026-07-16**, `result_phase2b_LALB.md`.]*

**REAL-T1 SKELETON stands** *[Nathan-reported]* — three lemmas, keystone verified, object reduced to **≤ 54 gauge classes**, D1's digit-automaton technology ready to redeploy. (Posed in `result_phase2b_H.md`: construct ℓ_k as the c_k-discounted Δ-return functional, verify Mᵀℓ_k = c_k ℓ_k by the gate algebra.)

**OPEN MYSTERIES (honestly held):**
1. **odd-k zero-amplitude** (H1) — selection-rule explanation refuted as posed (S and J both die on carry-integer-arithmetic).
2. **R32's near-λ₁ tower cluster** — NOT the circulant family (H1 showed the zero-amplitude family members are mid-spectrum).

**STANDING PROBE:**
- **G** — the L=4 partner (braid point three, rate-law point four). Flagged **hard**; needs the **direct-solver discipline** (no iterative eigensolvers).

**INSTRUMENT LAW (thrice-earned):** near the exceptional point use **direct/LU only** — the boundary degrades every iterative tool aimed at it (ESPRIT on the Jordan / Q6-B; eig+inv at G0c; ARPACK shift-invert at H2). Reconfirmed here indirectly (the symmetry commutators are dense-computed).

## Next session — the invariant hunt (Nathan, by hand)
The **21 closed classes** of the reduced **(ρ, γ)-chain**, looking for the conserved quantity that proves **no path from carried-diagonal ever returns to Δ**. That lemma closes **L-A**; L-A closes the **skeleton**; the skeleton delivers the **entire kinematic half** of the q=3 boundary in closed form — at which point the entrance exam is half over and the only thing left standing in Phase 2 is the **dynamical partner** itself.

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, THEOREM D1, Request F, the H-gates. No `r_q` value changes. The dynamical partner (non-family coalescing member) remains the D3 unknown.

_Reporting discipline: J is banked as a refutation (walk-back #14), diagnosed by the matching-norm identity with S (same carry-arithmetic breakage), not spun. The constraint extracted is filed as a constraint, not a theory. Provenance is tagged: D1 / circulant-family-completeness / H-gates are Claude-verified; L-A / L-B / the ≤54-gauge-class reduction / the S,J refutations are Nathan-reported session results (no probe artifact in-tree). The odd-k and R32 items are held open, not closed by narrative._
