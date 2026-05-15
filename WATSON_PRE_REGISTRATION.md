# WATSON_PRE_REGISTRATION — saddle-point / Watson lemma probe on R78/R79 bilinear sum

**Date:** 2026-05-14. Locked BEFORE Phase 2 execution. Probe slot: secondary routing from FAURE_DISPOSITION priority-1.

## Probe target

Apply Watson lemma / saddle-point / steepest descent / stationary phase to the R78/R79 bilinear
off-diagonal sum and derive a leading-order asymptotic for ε_k (or for the generating function
f(z) = Σ_k ε_k z^k).

Compare against:
1. **PADE_NUMERICAL_DISPOSITION** prediction: leading singularity at |z|≈1.57 (n=13 transient),
   asymptotic slow-mode at z≈1.016, complex-conjugate pair period 9.2 (θ≈0.68 rad), sign pattern
   `+ + − − − − − − − + + + +` with zero-crossing at n=9→10.
2. **FAURE_DISPOSITION** prediction: spectral radius ≤ 1/√E_min in Faure 2009 partially-expanding
   theorem; for Syracuse with 3-adic expansion the prediction is 1/√3 → singularity at |z|=√3≈1.732.
3. **R75/R76/R77 conversion**: does the asymptotic on ε_k convert to a polynomial-in-A bound on
   |μ̂_n(ξ)| via the bilinear pair operator T_M?

## Object being analyzed (verbatim from project)

**R78/R79 bilinear off-diagonal sum** (at q=3, the load-bearing case for c=7/45):

```
T_p := Σ_{a ∈ supp(F̂_p)} 1̂(p·a) · ψ_true(a)
```

where, per `result_78_extended.md` Theorems 78.4–78.6 (at r=3, p=3, J=3):

- `q = p^{r+1}, period = p^r, N = p^{r-1}, supp = {a ≡ 1 mod p in Z/p^r}, |supp| = p^{r-1}`
- `F̂_p(p·a) = p · e_q(1) · G(a)` with `|G(a)| = p^{(r+1)/2}` (FHAT_THEOREM_VERIFIED across 33 cells)
- `G(a) = Σ_{s=0}^{p^r − 1} e_q(P_a(s))` with `P_a(s) = ps − C_a · L_p(1 + ps)`
- `L_p(1 + ps) = Σ_{j=1}^{J} (-1)^{j-1}/j · (ps)^j` (Cochrane truncated p-adic log)
- `C_a = a · L̃_p^{-1} mod p^r`, `L̃_p = L_p(1 + p) / p`
- `s*(C_a) = (C_a − 1)/p mod p` (Cochrane Prop 4 saddle)
- `ψ_lead(a) = e_q(P_a(s*(C_a)))` (saddle-point leading order)
- `1̂(p·a) = Σ_{u=0}^{N-1} e_q(p·a·u)` (Dirichlet kernel)

The connection to ε_k goes through the Tao-recursion pair operator T_M and bilinear pair moments
`M_n(η) = Σ_ξ μ̂_n(ξ) μ̂_n(ξ·η)*` (per bilinear_pair_operator.py). The rate at which |ε_n|·2^n
varies as n grows is governed by the spectral radius of T_M on a deviation subspace.

At p=3, the sum `T_3 = Σ_a 1̂(3a) · ψ(a)` is EMPIRICALLY |T_3(r)| ∝ N^{0.522 ± 0.008}, R²=0.9976,
at r=8..20 (R79b). The rigorous bound from PATH2_BILINEAR is |T_3| ≤ 2N + O(p log p / N) using
inner Plancherel on c_2 + Dirichlet-kernel summation identity.

## Candidate techniques (Mode E to extract, then apply)

| Code | Technique | Pre-reg prior | Rationale |
|---|---|---|---|
| A | Watson's lemma (1918, real Laplace ∫_0^∞ f(t) e^{-λt} dt) | LOW | The bilinear sum is over a finite set Z/p^r, no Laplace transform structure. |
| B | Laplace's method (∫ φ exp(N f)) | LOW | Same as A — discrete sum, not exponentially-weighted integral. |
| C | Saddle-point / steepest descent (∫_C f(z) e^{x φ(z)} dz with saddle φ'(z*)=0) | **HIGH** | Cochrane Prop 4 IS a saddle-point method; already used at q=3 in R78.6. Family-level extension is the natural target. |
| D | Stationary phase (∫ φ e^{i N g} dt) | MODERATE | The bilinear sum DOES have oscillatory character phase e_q(P_a(s)); however the "asymptotic parameter" is r (level), and characters have unit modulus — saddle is degenerate. |
| E | Mellin-Barnes contour (∫ Γ(s)/(2πi) ds) | **HIGH** | Direct mechanism to convert coefficient asymptotic ε_k ↔ generating-function singularity. Required if we want ε_k asymptotic from f(z) Padé picture. |
| F | Multi-saddle-point with complex pair | **HIGH** | PADE finds complex pair at θ≈0.68 rad. Multi-saddle theory (two real saddles or complex-conjugate-pair saddles) is the exact machinery for cos(nθ+φ) modulation. |
| G | p-adic saddle-point (Cochrane Prop 4) | **HIGH** | Already used at r=3 in R78.6; family-level extension is the load-bearing technical question. |
| H | Watson + Borel resummation | LOW | Borel resummation handles divergent series; our ε_k sequence is bounded, no divergence. |
| I | Uniform asymptotic / Chester-Friedman-Ursell (coalescing saddles, Airy-type) | MODERATE | If PADE's two singularities at z≈1.016 and z≈1.57 are coalescing in the asymptotic limit, Airy-type expansion would apply; but no evidence of coalescence. |
| J | Darboux's method (coefficient asymptotic from algebraic singularities) | **HIGHEST** | Temme §2.4: directly extracts a_n asymptotic from generating-function f(z) with algebraic singularities. Matches the PADE picture (singularity structure of f(z) = Σ ε_k z^k). Single most direct route. |

## Hypothesis spectrum

- **SELECTED**: Technique applied to T_p sum (or to f(z) of ε_k) gives explicit asymptotic
  formula for ε_k at large k that matches both PADE singularity prediction and Faure spectral
  radius, AND converts via Plancherel/Tao recursion to a polynomial-in-A bound on |μ̂_n(ξ)|.
- **PARTIAL**: Asymptotic matches PADE quantitatively (or Faure quantitatively) but the
  conversion to |μ̂_n(ξ)| bound has a gap, OR matches one prediction but not both.
- **NO_FIT**: Technique categorically inapplicable (rare for Watson/saddle-point — universal tools).
- **BLOCKER**: R78/R79 sum form not extractable (FALSE — extracted above).
- **MODE_H_CIRCULAR**: Technique requires the closure asymptotic as input (very unlikely).

## Decision rule

If saddle-point applied to the R78/R79 bilinear sum produces a leading asymptotic for ε_k that:
- Predicts singularity at |z|=√3≈1.732 → consistent with Faure prediction
- Predicts complex-conjugate-pair structure with θ≈0.68 → consistent with PADE
- Numerically matches ε_k=8..13 within tight tolerance
- Converts to closure on |μ̂_n(ξ)| via existing R75/R76/R77 identities

Then SELECTED. If any one element fails: PARTIAL with precise specification of the failed step.

## Pre-registered priors on outcome

- SELECTED: 25-35%
- PARTIAL: 35%
- NO_FIT: 15%
- BLOCKER: 10% (mostly extracted now; revised down to ~5%)
- MODE_H_CIRCULAR: 5%

## Pre-registered secondary routing if not SELECTED

1. **Multi-singularity Flajolet-Sedgewick VI.4-VI.5** (PADE's recommended next step; complex-pair
   singularity asymptotics).
2. **Direct profinite semiclassical theory construction** (Faure's priority-3; multi-session monograph work).

## Discipline checklist

- [x] Mode E: techniques extracted VERBATIM from Temme 2013, Cohn, Manton, Glasner, Uppsala PDFs
- [x] Pre-reg locked before Phase 2 execution
- [x] R78/R79 sum form fully specified (above)
- [x] PADE + Faure predictions extracted and recorded
- [ ] Phase 2 — execute the asymptotic computation
- [ ] Phase 3 — compare against PADE + Faure + numerical ε_k=8..13
- [ ] Phase 4 — disposition

## Files referenced (all read VERBATIM)

- C:/Collatz/result_78.md (Cochrane attack obstructed at D=0)
- C:/Collatz/result_79.md (van der Corput obstructed)
- C:/Collatz/result_78_extended.md (R78.4-78.6 saddle-point closed form at r=3, q=3)
- C:/Collatz/PATH2_BILINEAR_FROM_CLOSED_FORM.md (family-level bilinear analysis; |T_p| ≤ 2N at r=3)
- C:/Collatz/bilinear_pair_operator.py (T_M definition, M_n bilinear moments)
- C:/Collatz/result_80_bilinear_attack.py (within-saddle-class structure at r=3)
- C:/Collatz/PADE_NUMERICAL_DISPOSITION.md (leading singularity at |z|≈1.57, complex pair, period 9.2)
- C:/Collatz/PADE_NUMERICAL_DATA.md (ε_k k=1..13 table)
- C:/Collatz/FAURE_DISPOSITION.md (Faure 2009 spectral radius 1/√3 prediction)
- C:/Users/Nate/OneDrive/Documents/Collatz Papers/Plancherel saturation/writeups/FHAT_THEOREM_VERIFICATION_RESULTS.md (|F̂_p^full(ξ)| = p^{(r+3)/2} on principal-unit support, 33 cells)
- C:/tmp/watson/Cohn_Watson_Lemma_Notes_Nebraska.txt (Watson lemma verbatim)
- C:/tmp/watson/Manton_Cambridge_DAMTP_Asymptotic_Methods.txt (Watson + Laplace + steepest descent)
- C:/tmp/watson/Temme_2013_Uniform_Asymptotic_Methods_Integrals.txt (Darboux's method §2.4, Mellin-Barnes §2.5, uniform/coalescing §4)
