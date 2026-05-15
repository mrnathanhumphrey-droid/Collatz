# Monotone Cumulants — Final Disposition (1-page)

**Date:** 2026-05-14
**Task:** Apply Hasebe-Saigo 2011 monotone cumulants to the Tao Syracuse
recursion; derive an asymptotic for `μ̂_n(ξ)`; compare against PADE / Faure
numerics. Task 4 of POST_COMPACT_NEXT_STEPS.md.

---

## Headline

**Leading-order closure derivation: IN HAND.** The coefficient `c = 7/45` is
reproduced by the monotone-cumulant framework as `κ_1^B(Off_j)` projected onto
the (1, 4)-eigenvector of R77's T_diag, with the (1/3)²:(2/3)² = 1:4 class-mass
identity and R75 Plancherel normalization. The derivation cites Hasebe-Saigo
2011 (Thm 4.5 cumulant definition, Thm 4.8 moment formula) and Hasebe monograph
Thm 3.26 (moment-cumulant formula via monotone partitions) verbatim.

**Full asymptotic closure: NOT IN HAND.** The rate-1/2 subdominant exponent
is *mechanistically* explained (per-step cumulant additivity + B-measurable
phase-twist factor decays at 1/2) but not *quantitatively* derived in the
monotone framework — its exponent rests on R77 §3 outside the cumulant
machinery. The coefficient −1/30 of the subdominant has no closed-form
derivation; it appears as `S_∞/14 = 7/(15·14)` with `14 = 2·7` conjecturally
from Plancherel bilinear normalization (R77 §6 open). PADE multi-spectral
structure (z ≈ 1.5..1.7 complex pair, period 9.2, Faure √3, k=7,8 ε deviation)
is *consistent* with the framework but not *derived* from it.

---

## Match against numerical anchors

| Predicted | Derived from monotone cumulants? | Numerical / spectral value |
|---|---|---|
| `c = 7/45` | **YES (rigorous fiberwise, conjectural at B-lift)** | 7/45 |
| Rate `1/3` | YES (R75 Plancherel, pre-existing) | 1/3 |
| Rate `1/2` subdominant | Mechanism only | (1/2)^n at k=2..6 |
| `−1/30` coefficient | Mechanism only | −1/30 (numerical fit) |
| PADE complex pair period 9.2 | Consistent, not derived | 9.2 |
| Faure √3 ≈ 1.732 | Consistent (intermediate cumulant scale) | √3 |

---

## Closure verdict

**c = 7/45 closure derivation IS in hand** at the leading-coefficient level,
modulo the conjectural B-valued lift of HS Thm 3.26. The lift is operational
and consistent with the 2026-05-14 numerical probe (`M_3_alt = 0.1078`,
`M_2 ≈ 10⁻⁷`, `M_3_distinct ≈ 10⁻⁵`), and the abelian-B Syracuse setting
makes the lift especially tractable (reduces to fiberwise scalar HS plus
integration over accumulator history).

The framework also explains the **pattern of failure of B-freeness** identified
on 2026-05-14: the alternating-with-repeated-index pattern (j_1, j_2, j_1) gives
a non-zero monotone cumulant signature via the peak-rule factorization
`E_B(X̃_{j_1} X̃_{j_2} X̃_{j_1}) = E_B(X̃_{j_2}) · E_B(X̃_{j_1}²)`, where under
marginal centering `E_B(X̃_{j_2})` retains a B-measurable phase-twist `Δ_{j_2}`
through the accumulator coupling. This is **exactly** the diagnostic predicted
by Hasebe monograph Defn 1.21.

---

## Caveats / Mode-E gaps

1. **B-amalgamated lift of HS Thm.** Hasebe-Saigo 2011 + Hasebe 2024 monograph
   develop monotone cumulants in the scalar case (state φ: A → C). The
   operator-valued / B-amalgamated extension used in Deliverables B-C is
   *conjectural*. No verbatim B-valued theorem in the closure-hunt corpus.

2. **Mixed monotone cumulants.** HS 2011 §6 references unpublished "in
   preparation" work for mixed cumulants. The mixed cumulant `M_3^B(X_{j_1},
   X_{j_2}, X_{j_1})` of Task 4 deliverable B is not a strict HS cumulant —
   it is a mixed moment that the framework computes via Hasebe Defn 1.21
   peak-rule factorization, not via a published mixed-cumulant formula.

3. **Closed-form 1/30.** The combinatorial factor 14 = 2 · 7 in
   1/30 = S_∞/14 is conjectured from Plancherel structure; not derived in
   monotone partition counting.

4. **PADE complex pair / Faure √3.** Framework-consistent but framework
   does not derive these values. They emerge from higher-cumulant /
   multi-spectral operator structure that the framework supports but has
   not been used to compute.

---

## Files

- `C:/Collatz/MONOTONE_CUMULANTS_A_VERBATIM.md` — Deliverable A: verbatim HS / Muraki / Hasebe-monograph definitions
- `C:/Collatz/MONOTONE_CUMULANTS_B_SYRACUSE.md` — Deliverable B: per-step Syracuse monotone cumulants
- `C:/Collatz/MONOTONE_CUMULANTS_C_ASYMPTOTIC.md` — Deliverable C: asymptotic derivation
- `C:/Collatz/MONOTONE_CUMULANTS_D_COMPARISON.md` — Deliverable D: comparison to PADE + Faure
- `C:/Collatz/MONOTONE_CUMULANTS_DISPOSITION.md` — this file
- `C:/Collatz/verify_monotone_diagnostic.py` + `experiments_output/monotone_diagnostic_n3.json` — n=3 numerical anchor (M_3_alt = 0.1078)
- `C:/Collatz/AMALG_FREENESS_MOMENT_CALCULATION.md` §8 — marginal-centering reading (load-bearing for B-lift)
- `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_saigo_2011_monotone_cumulants.pdf`
- `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_monotone_probability_theory_monograph.pdf`
- `C:/Users/Nate/OneDrive/Documents/closure hunt/muraki_2003_five_independences_kyoto_precursor.pdf`
