# Probe W4 — the EP witnesses at L=4 (overlap, defectiveness) + same-instrument re-base

**Date:** 2026-07-18  Lambda A10 (terminated). FULL operator. Goal: overlap₄ and g₄ as EP coalescence
witnesses, with L=2,3 re-based on one instrument. **Outcome: W3 (re-base) is a clean win; the L=4 witnesses
are NOT reliably obtainable in double precision — the near-EP makes the c₀ eigenvector ill-conditioned, a
validated finding, not a number.**

## W3 — same-instrument re-base (dense eig, L=2,3): VALIDATED
The banked 0.998/0.99999 (fork-b) and 0.0505/0.0188 (F2-4) now sit on one instrument (dense full-operator eig):

| L | partner | c₀-mode | **overlap** \|⟨r₀,rₚ⟩\| | **g** = B[tow,kin] | B[kin,tow] | defectiveness g/\|Δ_c₀\| |
|---|---|---|---|---|---|---|
| 2 | 0.3468267 | 0.3439153 | **0.99834443** | **+0.050531** | −0 (≤5e-9) | 0.050531/2.911e-3 = **17.4** |
| 3 | 0.3332363 | 0.3333359 | **0.99998600** | **+0.018820** | −0 | 0.018820/9.958e-5 = **189.0** |

Reproduces the banked series exactly, nails the convention (**g = B[tow,kin]**, the tower→kinematic coupling;
the 2×2 is lower-triangular, B[kin,tow]=0), and confirms the right eigenvectors are coalescing (overlap → 1).
Cross-method caveat removed.

## L=4 — the split, and why the witnesses are not numbers
**Method:** block-2 orthogonal power iteration on the full operator (Wilson's spec) → **failed**: the partner
and c₀ right eigenvectors are numerically parallel (overlap→1), so the orthogonal frame's second vector
collapses onto the complex pair, not c₀. Switched to **oblique deflation** (partner right rₚ + left ℓₚ from
the tower, deflate, power-iterate for c₀).

**The eigenvalues split cleanly (gate PASS):** partner = **0.33349990132** (= G4), c₀ = **0.3333333333** (=
closed form, machine precision). Their gap = **1.6657e-4 = |Δ₄| (detuning-vs-c₀)** — an *independent third
confirmation* of the partner↔c₀ gap that the braid sequence carries.

**But the overlap/g witnesses are artifacts, proven by cross-validation at L=3:** running the *same deflation
pipeline* at L=3 recovers the c₀ **eigenvalue** exactly (0.3333358765) but its **eigenvector** has only
**1.25% overlap** with the true (dense) c₀ eigenvector (`⟨r0_defl, r0_dense⟩ = 0.0125`), giving overlap 0.018
where the true value is **0.99999**. The tower-embedded partner is exact (`⟨rₚ_tower, rₚ_dense⟩ = 1.0`, dense
partner has 1e-24 mass on γ=0), so the failure is entirely in the **c₀ eigenvector**: at overlap→1 its
condition number ~1/(1−overlap²) ≈ 5·10⁴ (L3) and larger at L4, so double-precision power iteration returns a
vector with the right Rayleigh quotient but the wrong direction.

⟹ The reported L=4 numbers overlap₄ = 0.034 and g₄ = 0.0046 (defect 27.7) are **method artifacts, discarded** —
NOT a drop in the coalescence. Both orthogonal block-2 and oblique deflation fail the eigenvector for the same
reason: the witness quantity is ~1 and therefore ill-conditioned.

## The finding (report loudly)
The overlap witness **cannot be evaluated at L=4 in double precision** — and *that failure is itself the
signature*: c₀ and the partner right eigenvectors have merged to the point where the separate c₀ eigenvector is
numerically non-recoverable (a defective-EP approach). What survives cleanly is the **eigenvalue** structure:
partner and c₀ split by exactly 1.666e-4, the same gap the braid tracks. The defectiveness meter's third point
(17 → 189 → ?) is **inaccessible numerically**; to obtain it one needs **extended precision** (mpmath) — but
dense eig at 236,196 states is infeasible even there, so the L=4 overlap/g require an **analytic/perturbative**
route (the effective-model unfolding), not a numerical one. The eigenvalue gap is the accessible L=4 witness,
and it is confirmed.

## Deliverables
W3 series (dense, one instrument) is the solid product. L=4: eigenvalue split confirmed (gate pass, gap =
1.666e-4); overlap₄/g₄ withheld as numerically ill-conditioned (validated at L=3). W4 c₀-γ-profile likewise
rests on the corrupted eigenvector and is withheld.

Probes `probes/w4_ep_witness.py` (re-base + block-2), `probes/w4b_deflate.py` (deflation),
`probes/w4c_validate.py` (the L=3 validation that caught the artifact); logs in `logs/`.
