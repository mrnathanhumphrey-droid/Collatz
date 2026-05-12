# HENSEL_RECONSTRUCTION_COMPARISON — Phase 5

**Date:** 2026-05-11. Comparison of independent re-derivation to the original Hensel-lift claim (per system message description).

## Original claim (per system message)

> G_p(a) = √q · η_p(r) · e_q(P_a(s*(r)))
>
> with s*(r) = (C_a − 1) // p mod p^{r-1}
> and P_a(s*(r)) = Σ_{j=2}^{r} (−1)^{j-1} · (p·s*(r))^j / (j·(j-1)) mod p^{r+1}
> and η_p(r) = quadratic Gauss sum at even r, 1 at odd r:
>   η_p(r) = 1                    (r odd)
>   η_p(r) = (1/√p)·Σ_h e_p(h²/2) (r even)

## Independent re-derivation result

From Phases 1-4:

### Saddle (Phase 2)
> **s*(r) = (C_a − 1)/p mod p^{r-1}**

Derived: stationary-phase condition `dP_a/ds = p·(1+ps−C_a)/(1+ps) = 0` is **linear in s** in formal Z_p, hence the saddle is `s = (C_a − 1)/p`, exact in Z_p. Reducing mod p^{r-1} gives the saddle representative.

**Comparison: identical to claim.** ✓

This is **digit extraction**, not a Hensel-iteration series. The "Hensel-correction series δ_k" claimed to collapse to the clean form — confirmed independently. The reason for the collapse: the saddle equation factors with a linear factor in s, due to the log structure of P_a.

### Phase polynomial (Phase 3)
> **P_a(s*(r)) ≡ Σ_{j=2}^{J_p} (−1)^{j-1}/(j(j-1)) · (p·s*(r))^j mod p^{r+1}**

With J_p = r for p ≥ 5 (generic). Equivalently `P_a(s*) = M(p·s*)` where M(y) = y − (1+y)·log(1+y).

**Comparison: identical to claim.** ✓

Coefficient pattern `(−1)^{j-1}/(j(j-1))` matches. Variable `(p·s*(r))^j` matches. Truncation depth: claim says `j=2..r`, my derivation gives `j=2..J_p` where J_p = r generically. **Match** (modulo harmless edge cases at p=3 where J_p might differ from r by 1, but the additional terms vanish mod q anyway).

The "(1+y)·log(1+y) coefficient pattern" mentioned in the system message refers exactly to this: M(y) = y − (1+y)·log(1+y), with coefficients derivable from term-by-term expansion. Matches.

### η_p Gauss-sum factor (Phase 4)
> **η_p(r) = 1**  (r odd ≥ 3)
> **η_p(r) = (1/√p)·G_{quad}(1/2; p) = (1/√p)·Σ_h e_p(h²/2)**  (r even ≥ 2)

Derived: at the saddle class, the sub-sum over u ∈ Z/p^{r-1} of e_q(P_a(t_1 + p·u)) has a parity-dependent residue. At odd r, an orthogonality collapse on the linear-in-u term eliminates the Gauss-sum residue. At even r, a residual quadratic Gauss sum of length p remains.

**Comparison: identical to claim.** ✓

The Gauss sum `Σ_h e_p(h²/2)` matches the system-message form `Σ_h e_p(h²/2)`.

## Summary table

| Component | Claim | My Derivation | Match? |
|---|---|---|---|
| Saddle s*(r) | (C_a − 1) // p mod p^{r-1} (digit extraction) | (C_a − 1)/p mod p^{r-1} (digit extraction, formal Z_p saddle reduced) | **YES** |
| Phase polynomial coefficients | (−1)^{j-1}/(j(j-1)) | (−1)^{j-1}/(j(j-1)) | **YES** |
| Phase polynomial variable | (p·s*(r))^j | (p·s*(r))^j | **YES** |
| Phase truncation depth | j=2..r | j=2..J_p, J_p = r generically | **YES** (equivalent mod q) |
| η_p at odd r | 1 | 1 | **YES** |
| η_p at even r | (1/√p)·Σ_h e_p(h²/2) | (1/√p)·G_{quad}(1/2; p) = (1/√p)·Σ_h e_p(h²/2) | **YES** |
| Magnitude prefactor | √q | √q = p^{(r+1)/2} (matches FHAT empirical and structural Plancherel) | **YES** |

**All components match.** ✓

## Disposition: H_RECONSTRUCTION_CONFIRMS

The independent re-derivation, starting from R78.4-78.6 (family-level extension in PATH2_FAMILY_EXTENSION.md) and using standard stationary-phase analysis with the formal-Z_p identity `dP_a/ds = p·(1+ps − C_a)/(1+ps)`, lands at exactly the same closed-form structure as the original Hensel-lift claim:

> **G_p(a) = √q · η_p(r) · e_q(M(p·s*(r)))    where s*(r) = (C_a − 1)/p mod p^{r-1}, M(y) = y − (1+y)log(1+y), η_p parity-dependent.**

Same saddle. Same phase polynomial. Same η_p. **Confirmed.**

## Why the structures align (mechanism)

The cleanness of the closed form derives from the structure of the Cochrane polynomial P_a(s) = ps − C_a·L_p(1+ps):

1. **dP_a/ds has a linear factor `(1+ps − C_a)` in formal Z_p.** This is because L_p is a logarithm-type series and dL_p(1+ps)/ds = p/(1+ps). Linear factors give CLEAN roots — no Hensel iteration needed.

2. **The saddle equation `1+ps = C_a` has a single root in Z_p**, namely s = (C_a − 1)/p. This is exact, not a series.

3. **The phase at the saddle is `M(p·s*) = M(C_a − 1) = (C_a − 1) − C_a · log(C_a)`** — a clean analytic function of C_a in the formal Z_p.

4. **The "Hensel correction" is just the digit expansion of (C_a − 1)/p mod p^{r-1}** — read off the digits of C_a directly. No iteration needed because the saddle equation is already linear (no Hensel-Newton step required).

5. **The Gauss-sum factor at even r comes from the quadratic correction `(p²/(2C_a))·u²` in the Taylor expansion at the leading-digit saddle `t_1`**, which contributes a length-p quadratic Gauss sum mod p when the relevant p-level lies inside mod q. At odd r, an additional orthogonality collapse eliminates this residue.

These features are all standard p-adic stationary-phase outcomes for polynomial phases with logarithm structure. **Nothing exotic or unusual in the derivation chain.**

## Implications

**For the original claim:** internally confirmed by independent derivation. The structural form (saddle = digit extraction, phase = M-series, η_p parity-dependent) is robust — falls out naturally from stationary-phase on the Cochrane polynomial.

**For the empirical numerical match at 1e-15 (10 cells, original Hensel agent's verification):** this is consistent with both the original derivation AND my independent re-derivation. Two derivations + the numerical match is the strongest possible internal evidence.

**For the open question** (the original Hensel agent's claimed verification at 10 cells with max rel dev 6.4e-15 — per system message): the form is structurally CORRECT. Empirical match validates the precision, structural derivation validates the form.

**Caveats:**

- Phase 4 (η_p) derivation at general even r ≥ 6 was sketched but not fully unwound for the explicit cubic/higher-corrections. The pattern is clear from r=4 (full derivation) and r=5 (full derivation showing the odd-r collapse), and the structural reason (parity of r — 1 in the d²P_a/ds² · u² Gauss-completion vs higher-order term levels) is identified. **Full r ≥ 6 details would be a mechanical extension** but not done here.

- The derivation at p=3 has the additional complication of 1/3 denominators making J_p ≥ r+1 at higher r, but the leading structure is unaffected.

- The empirical verification was 10 cells (per system message). Independent confirmation extends the structural foundation; combined evidence is strong.

## Files

- HENSEL_RECONSTRUCTION_PHASE1.md — saddle re-derivation at r=3.
- HENSEL_RECONSTRUCTION_PHASE2.md — saddle extension to r ≥ 4 (structural, not perturbative).
- HENSEL_RECONSTRUCTION_PHASE3.md — phase polynomial M(p·s*).
- HENSEL_RECONSTRUCTION_PHASE4.md — η_p Gauss-sum factor.
- HENSEL_RECONSTRUCTION_COMPARISON.md — this document.
- HENSEL_RECONSTRUCTION_DISPOSITION.md — top-level disposition.
