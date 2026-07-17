# FG_DISPOSITION

**Date:** 2026-05-13.
**Probe:** Furstenberg-Guivarc'h random-walks-on-locally-compact-groups expressibility test for Syracuse mu_infinity Fourier decay.
**Mode:** E (verbatim theorem hypotheses from PDF, no inheritance from prior project files).
**Pre-registration:** `C:/Collatz/FG_PRE_REGISTRATION.md`.

---

## Headline

**No SELECTED.** The Furstenberg-Guivarc'h corpus — even with the most pre-cleared candidates (BFLM, Li-Bourgain, Li 2018 Theorem 1.7) — closes **NO_FIT-dominant** with two MODE_H_CIRCULAR sub-cases at Phase 3 extension and one BLOCKER (UNVERIFIABLE-PHASE-0).

The category-of-object barrier is the same one identified in the prior 5 probe arcs: **Syracuse on (Z/3^n)* is an abelian profinite chain**; FG-school theorems target **non-abelian semisimple Lie or p-adic Lie groups with proximal / non-arithmetic structure**. The abelian Syracuse setting either:
1. **Categorically fails** hypotheses (no proximality, no flag variety, no matrix-product, no reversibility, no Lyapunov exponent), or
2. **At Phase 3 extension**, reduces to the closure target itself (the imaginary-line spectral gap of the transfer operator IS the Fourier decay of the chain stationary measure).

This is the **sixth category-of-object barrier** in the systematic obstruction map for c=7/45 closure, parallel to and reinforcing the prior five.

---

## Summary table

| Code | Theorem(s) | Phase 0 | Phase 1 | Phase 2 | Phase 3 | Disposition | File |
|---|---|---|---|---|---|---|---|
| A | Furstenberg 1963 Thm 8.5/8.6 (Lyapunov, SL(m,R)) | Extracted | NO_FIT (group/walk) | SHAPE_MISMATCH (Lyapunov ≠ Fourier decay) | STRUCTURALLY_BLOCKED (abelian → trivial Lyapunov) | **NO_FIT** | FG_A_HYPOTHESES.md |
| B | Le Page 1982 Thm 1, 6 + Guivarc'h 1980 (renewal/contraction/escape rate) | Extracted | NO_FIT (SL(d,R) group) | SHAPE_MISMATCH | REDUCES_TO_TARGET (profinite extension = R75/R76/R77 Plancherel) | **MODE_H_CIRCULAR** | FG_B_HYPOTHESES.md |
| C | BFLM 2007/2011 Thm A + Cor C (effective stationary measure on T^d) | Extracted | NO_FIT (SL_d(R) on T^d, fixed-d, proximality) | STRONG MATCH | STRUCTURALLY_BLOCKED (proximality + fixed-d) | **NO_FIT** | FG_C_HYPOTHESES.md |
| D | Benoist-Quint 2016 Thm 1.1, 1.3 (qualitative stationary = Haar) | Extracted | NO_FIT (semisimple Lie, Zariski-density vacuous in abelian) | SHAPE_MISMATCH (qualitative not quantitative) | STRUCTURALLY_BLOCKED (abelian → semisimple gap) | **NO_FIT** | FG_D_HYPOTHESES.md |
| E | Li 2018 / Bourgain 2010 (discretized sum-product Fourier decay R^n) | Extracted | NEEDS_PROOF / NO_FIT (R^n vs (Z/3^n)* ambient + iteration mapping) | PARTIAL (polynomial decay with fixed exponent, not any-A) | STRUCTURALLY_BLOCKED (k-fold mult. conv. doesn't reconstruct Syracuse stationary) | **NO_FIT** | FG_E_HYPOTHESES.md |
| F | Li 2018 Thm 1.7 + Le Page CLT (Fourier decay of Furstenberg measure on flag variety) | Extracted | NO_FIT (R-split reductive R-group, no flag variety on abelian) | STRONG MATCH | REDUCES_TO_TARGET (imaginary-line spectral gap = closure target) | **MODE_H_CIRCULAR** | FG_F_HYPOTHESES.md |
| G | Saloff-Coste 2004 (mixing time, spectral gap, log-Sobolev, finite group RW) | Extracted | NO_FIT (reversibility / detailed balance) | SHAPE_MISMATCH (mixing time, not Fourier decay) | STRUCTURALLY_BLOCKED (non-reversible chain) | **NO_FIT** | FG_G_HYPOTHESES.md |
| H | Varopoulos-Saloff-Coste-Coulhon heat-kernel | UNVERIFIABLE | — | — | predicted STRUCTURALLY_BLOCKED (compact-vs-noncompact volume growth) | **BLOCKER** | FG_H_HYPOTHESES.md |

(Phase 0 = verbatim extraction status. Phase 1 = hypothesis × input matrix. Phase 2 = conclusion-shape vs closure target. Phase 3 = profinite extension feasibility.)

Total: **7 candidates with extractable statements** in Phase 0 (A through G). 1 candidate (H) UNVERIFIABLE-PHASE-0 / BLOCKER.

Net dispositions: 5 NO_FIT (A, C, D, E, G), 2 MODE_H_CIRCULAR (B, F), 1 BLOCKER (H).

---

## Final disposition: **NO_FIT** (with secondary MODE_H_CIRCULAR finding).

No theorem in the Furstenberg-Guivarc'h corpus accepts Syracuse mu_n on (Z/3^n)* as a hypothesis-satisfying instance under Mode E discipline. The closest candidates (BFLM, Li 2018 Thm 1.7) deliver the right conclusion shape (polynomial-in-A Fourier decay of stationary measure), but their hypotheses categorically fail on the abelian profinite setting; the natural Phase 3 extension reduces to the closure target itself (Mode H circular).

This **completes the sixth category-of-object barrier** in the obstruction map. The unifying gap surfaced across all six probe arcs is now:

> **No theorem in any scanned corpus operates on the correct category-of-object for Syracuse mu_infinity: an abelian profinite multiplicative-group Markov chain with non-reversible 2-adic Geom step distribution and a multi-regime asymptotic deviation eps_k.**

---

## Secondary routing recommendation

Per pre-registration, the routes flagged in advance:

1. **BGT regular variation (Bingham-Goldie-Teugels)** — operates on the k=7 jump signature in eps_k. The two-regime k<7 vs k≥7 break in |eps_k|*2^k suggests slowly-varying correction to the rate-½ envelope. **PRIORITY: HIGH.** BGT machinery is *categorically correct* for sequences with regularly-varying tails; the input is **directly the eps_k sequence**, no chain-side category issues. The Karamata / de Haan analysis on the *sequence* eps_k (a real-number sequence) is the cleanest re-categorization.

2. **Igusa local zeta / functional equation** — the (1+3)^u algebraic root in R78 D=0 disambiguation. Igusa Z(s; f) = ∫_{Q_3} |f(x)|^s dx for f a polynomial in Q_3[u]; gives meromorphic continuation + functional equation. If the *generating function* of eps_k can be identified with an Igusa local zeta (via the (1+3^n) leading-mode-identity bookkeeping), the closure rate emerges from meromorphic-continuation poles. **PRIORITY: MODERATE-TO-HIGH.** Requires explicit identification of f.

3. **Heat-kernel on profinite tree** — pre-classified as VSC-style, structurally blocked at the compact-vs-noncompact gap. **PRIORITY: LOW** (covered by Candidate H disposition).

4. **Adelic Mellin construction** — restore the archimedean place. Connects with the BT_DISPOSITION finding that the c=7/45 constant is archimedean. The Mellin transform on Q*_A / Q* with idele-class-character structure may give a global statement. **PRIORITY: MODERATE.** Technically demanding; payoff is the only route that explicitly addresses the archimedean-visibility issue.

5. **NEW candidate surfaced during Phase 0:** The **transfer operator on (Z/3^n)*** as a direct re-instantiation of the Li-Dolgopyat spectral-gap framework. This is what FG_F Phase 3 surfaces as Mode H circular — but it's also the *most natural* path to a fully rigorous proof of the rate-½ off-diagonal eigenvalue conjectured in R77 (Conjecture 77.2). **PRIORITY for downstream R78 work, not for closure routing**: this is what Nisoli Theorem 2.15 application (R77 § 6) was supposed to do; the FG arc confirms the framework is right but the rigorous spectral analysis on the abelian transfer operator IS the open analytical step, not the closure route.

### Recommended top-priority secondary route: **BGT regular variation on eps_k**.

Rationale: BGT operates **directly on the eps_k sequence as a real-number object**. No chain-side category mismatch, no group structure required, no reversibility, no proximality. The 8-coefficient data + k=7 jump signature is exactly the input BGT/de Haan machinery is designed for. Whether the slow-variation analysis matches is testable from existing eps_k data alone (no new chain computation needed).

Igusa local zeta is the second-priority route — requires more explicit algebraic identification of f, but operates on a *different* object (the generating function), categorically distinct from chain-on-(Z/3^n)*.

---

## Surprises in the inputs

### Surprise 1: the k=7 jump is a multi-regime signature that NO single FG theorem handles

All FG theorems deliver **single-rate** Fourier-decay bounds (Li 1.7: ξ^{-ǫ_1}; BFLM: e^{-c_2 n/M}). The Syracuse eps_k sequence has:
- k=2..6 plateau at |eps_k|*2^k ≈ 0.04 (rate-½ envelope).
- k=7 jump to 0.150, k=8 to 0.191.
- Sign pattern (+, +, −, −, −, −, −, −).

This **multi-regime structure** is incompatible with single-spectral-gap-eigenvalue Fourier-decay theorems. Even if BFLM/Li 1.7 hypotheses were satisfied (they aren't), the conclusion would predict monotone single-power decay — which the eps_k data violates at k=7.

The R77 conjecture of "off-diagonal mode at rate ½ overlaid on T_diag rank-1" naturally accommodates a multi-regime structure: T_diag's eigenvalue 0 mode (the (1,-1) direction) is exactly killed; the rate-½ mode lives on (1,4); a *third* mode (off-diagonal cross-frequency, currently the off-diagonal correction in R77) is a candidate explanation for the k=7 jump.

**Interpretation:** the FG framework is the right *category* but with the wrong *number of modes*. A two-mode or three-mode analog of Li 1.7 — Fourier decay of a stationary measure under a *non-simple* spectral structure (multiple Lyapunov-like exponents) — would be the right shape. This is not in the FG corpus; it's a research direction.

### Surprise 2: R77 T_diag spectrum already implements the FG transfer operator structure on (Z/3^n)*

The R77 derivation T_diag = (1/5)·[[1,1],[4,4]] with eigenvalues {0,1} on (1,-1) and (1,4) eigenvectors **IS** the abelianized version of the Furstenberg-school transfer-operator spectrum, restricted to the 2-dim "class-resolved" subspace. The eigenvalue 1 on (1,4) preserves Plancherel mass (= mass of trivial / unique-stationary mode); the eigenvalue 0 on (1,-1) is "instantly killed" mode (= mass-conservation null space).

This means **R77 already does what an FG-on-(Z/3^n)* analog would do** at the diagonal-transfer-operator level. The off-diagonal corrections (cross-frequency bilinear coupling at rate ½) are exactly the **Dolgopyat-oscillation** content of Li 1.7's Phase 2 (the spectral-gap-on-imaginary-line proof) — but at level of *characters* on (Z/3^n)* rather than non-arithmeticity on cocycles in SL(d,R).

**Interpretation:** the FG framework's machinery, when transferred to abelian profinite, IS the R77 framework. The remaining open step in R77 (Conjecture 77.2: rigorous rate-½ off-diagonal eigenvalue) is the structural analog of "spectral gap of P_{ib} for b ≠ 0" in Li 2018. The FG corpus doesn't supply this for the abelian profinite case, but the *form* of the analytical step needed is now clear.

### Surprise 3: BGT (regular variation) is uniquely categorically clean

Of the four secondary routes, BGT is the only one that **operates on the eps_k sequence as a real-number object**. It has no group-side hypotheses to fail, no chain-side iteration mapping to verify. The eps_k sequence is what's actually computed (input 4), the k=7 jump is the structural surprise — BGT's slowly-varying-correction-to-power-law is *exactly* the test for this kind of pattern.

This wasn't surprising to the pre-registration ("priority: high"), but it's worth flagging that **BGT is the only secondary route that side-steps the entire category-of-object obstruction map** by re-categorizing the problem to the eps_k sequence object.

---

## Files produced

- `C:/Collatz/FG_PRE_REGISTRATION.md` (pre-reg, locked before selection)
- `C:/Collatz/FG_A_HYPOTHESES.md` (Furstenberg 1963 Lyapunov)
- `C:/Collatz/FG_B_HYPOTHESES.md` (Le Page + Guivarc'h)
- `C:/Collatz/FG_C_HYPOTHESES.md` (BFLM)
- `C:/Collatz/FG_D_HYPOTHESES.md` (Benoist-Quint)
- `C:/Collatz/FG_E_HYPOTHESES.md` (Li-Bourgain sum-product)
- `C:/Collatz/FG_F_HYPOTHESES.md` (Li 1.7 + Le Page CLT — most-pre-cleared)
- `C:/Collatz/FG_G_HYPOTHESES.md` (Saloff-Coste)
- `C:/Collatz/FG_H_HYPOTHESES.md` (VSC heat-kernel — UNVERIFIABLE)
- `C:/Collatz/FG_DISPOSITION.md` (this file)

PDF extractions (UTF-8 from pypdf 6.10.2):
- `C:/tmp/fg/furstenberg_1963.txt`, `guivarch_1980.txt`, `le_page_1982.txt`, `guivarch_raugi_1985.txt`, `benoist_quint_book.txt`, `benoist_quint_2016.txt`, `bflm_2007.txt`, `bflm_2011.txt`, `bourgain_sumproduct_rn.txt`, `li_renewal_spectral.txt`, `saloff_coste_2004.txt`, `aldous_fill.txt`, `lindenstrauss_mohammadi_2022.txt`.

No git operations performed (per discipline).

---

End disposition.
