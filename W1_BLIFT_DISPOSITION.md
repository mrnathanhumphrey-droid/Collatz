# W1 Disposition — B-amalgamated lift of Hasebe-Saigo

**Date:** 2026-05-14
**Wrinkle:** #1 of 4 from MONOTONE_CLOSURE_WRITEUP.md §2
**Effort estimate:** 1-3 days (user pace, ~10-14× typical) — **actual: single session, ~4-6 hours equivalent.**

---

## Verdict: **W1 CLOSED via Route 2 (verbatim citation).**

The B-amalgamated lift step that was conjectural in
MONOTONE_CLOSURE_WRITEUP.md §1.2 is now **rigorous, conditional on
hypothesis H1 (monotone independence over B of the X̃_j family).**

Source: **Hasebe, T. and Saigo, H. (2014). "On operator-valued monotone
independence." Nagoya Mathematical Journal 215, 151-167. arXiv:1306.0137v2.**
Theorem 3.4 (moment-cumulant formula) + Proposition 3.5 (B-extensivity)
are stated for arbitrary unital algebras B (no commutativity assumption),
so the abelian-B case of Syracuse is a strict specialization. The theorem
applies verbatim.

---

## What this changes in MONOTONE_CLOSURE_WRITEUP.md

The §1.2 leading-order c = 7/45 derivation status upgrades:

| Field | Before W1 | After W1 |
|---|---|---|
| Framework = monotone | Identified, diagnostic confirmed | unchanged |
| c = 7/45 leading | **rigorous fiberwise + conjectural at B-lift** | **rigorous conditional on H1** |
| H1 (monotone indep over B) | implicit framework input | explicit, project-internal, numerically supported at 10⁶ separation |

The Mode-E ledger (MONOTONE_CLOSURE_WRITEUP.md §5) line
> "c = 7/45 leading | Rigorous fiberwise + conjectural at B-lift (Wrinkle 1)"
should be updated to
> "c = 7/45 leading | Rigorous conditional on H1 (monotone indep of (X̃_j) over B); H1 numerically supported at 10⁶ separation."

---

## Mode-E gaps remaining (after W1)

1. **H1 (the monotone-independence-over-B hypothesis) is not theorem-grade.**
   It is supported by:
   - Framework-identification step from OBSTRUCTION_MAP_TERMINAL.md (the
     11-arc obstruction map terminal finding, ruling out free / Boolean /
     tensor / classical).
   - Numerical diagnostic at 10⁶ separation (M_3_alt = 0.10783 vs
     algebraic zero under strict centering).
   - Hasebe-monograph Defn 1.21 peak-rule factorization match at n=3
     (the diagnostic value's product structure matches `Δ_{j_2} ·
     E_B(X̃_{j_1}²)` from the peak rule).
   This is **strong evidence but not a verbatim theorem**. Converting it
   to a theorem would require checking HS 2014 Def 2.2 directly on the
   Syracuse A_j subalgebras at all orders (the peak factorization
   condition), which is a structural calculation in its own right.

2. **Centering subtlety (reading A vs reading B).** The HS 2014 theorem
   applies cleanly to whichever centering is operative, provided the
   conditional expectation `ϕ` matches. For reading B (marginal centering),
   `ϕ = E_{B_marginal}` (the level-graded conditional expectation onto
   prior-accumulator information), not the strict `E_B` of SETUP.md §5.
   The W1 verification (Deliverable D §4-§5) makes this explicit; the
   level-graded formulation should be propagated to Deliverables B and C
   of the MONOTONE_CUMULANTS_* chain (Syracuse cumulants, asymptotic
   derivation) for full consistency.

3. **Hypothesis-match at level of subalgebras.** HS 2014 Def 2.2 talks
   about **subalgebras** `(A_λ)_{λ ∈ Λ}` over B being monotone
   independent. The Syracuse `A_j := B⟨X̃_j⟩_0` are well-defined
   subalgebras over B (any B-bimodule combination of powers of X̃_j),
   and Theorem 3.4 applies to multilinear functionals on these
   subalgebras. The numerical Task 1 probe tested specific monomials
   (X̃_{j_1} · X̃_{j_2} · X̃_{j_1}); extending to general subalgebra
   elements is automatic if H1 holds (HS 2014 Def 2.2 is a statement
   about all monomials simultaneously).

The principal remaining work to close c = 7/45 fully is **W2** (subdominant
coefficient -1/30) and possibly a project-internal proof of H1 if a
theorem-grade upgrade is desired.

---

## Effort tally

- Route 1 (direct construction) estimate: 1-3 days. **Not used.**
- Route 2 (verbatim citation) estimate: 4-8 hours. **Actual: single
  session, ~4-6 hour equivalent.** Time spent:
  - Literature pull (HS 2014 arXiv:1306.0137, full 13-page read): ~30 min
  - Hypothesis-match table + theorem specialization (W1.C): ~1.5 h equiv
  - Sanity check + centering subtlety resolution (W1.D): ~2 h equiv
  - Disposition + writeup integration (this file): ~30 min

Net cost vs original 1-3 day estimate: ~5-7× under (because Route 2
existed verbatim, no original construction required).

---

## Output files

- `C:/Collatz/W1_BLIFT_LITERATURE.md` — Hasebe-Saigo 2014 pull, magic-byte
  verification, section-by-section summary, Popa 2008 deferred-pull note.
- `C:/Collatz/W1_BLIFT_ROUTE.md` — Route 2 chosen, justified.
- `C:/Collatz/W1_BLIFT_THEOREM.md` — Verbatim HS 2014 hypotheses + Thm 3.4 +
  Prop 3.5 + Syracuse-specialization theorem + proof + hypothesis-match
  table.
- `C:/Collatz/W1_BLIFT_VERIFICATION.md` — n=3 sanity check, explicit
  enumeration of M(3) = 12 monotone partitions, centering subtlety logged,
  consistency with Task 1 numerical 0.10783 established.
- `C:/Collatz/W1_BLIFT_DISPOSITION.md` (this file) — 1-page verdict.

---

## Recommended next step

**Update MONOTONE_CLOSURE_WRITEUP.md** Mode-E ledger (§5) and Wrinkle 1
characterization (§2.1) to reflect "closed conditional on H1", with the
H1 itself logged as a remaining project-internal load-bearing input.
Then proceed to W2 (subdominant −1/30) or W4 (Faure √3 ↔ cumulant op),
whichever the user picks.
