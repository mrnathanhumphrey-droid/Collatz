# RESULT — P-HYDRA FAMILY: unified 2-parameter Plancherel law S_{k+1}/S_k → q(p−1)/(p+1); Collatz is DOUBLY critical (2026-07-28)

**Probe:** `probes/probe_phydra_family.py`. Wilson's Q1/Q2/Q3 on the "class" of q-difference Plancherel constants. The genuine
sibling family varies the **valuation prime p** (Siegel (q,p)-Hydra: `x ↦ (qx+1)/p^v`), not just q. Chain on `(ℤ/q^k)*`:
`r ↦ (qr+1)·p^{−v}`, `v~Geom(1−1/p)`, `S_k = X_k − X_{k−1}`, `X_k = q^k·Σπ²`. Gate: p=2,q=3 reproduces the certified
S-ladder (S_1=2/3, S_2=10/21) to 5.6×10⁻¹⁶.

## The finding — a clean unified two-parameter law
```
   (q,p)   meas S_{k+1}/S_k    q(p-1)/(p+1)     |diff|
   (3,2)        1.00092           1.00000       9e-4    <- Collatz (boundary, slow)
   (5,2)        1.66675           1.66667       8e-5    <- q-family (=q/3)
   (7,2)        2.33321           2.33333       1e-4
   (3,5)        2.00021           2.00000       2e-4    <- p-family (=3(p-1)/(p+1))
   (3,7)        2.24999           2.25000       1e-5
   (5,7)        3.75000           3.75000       4e-9    <- MIXED, both != base
   (7,5)        4.66667           4.66667       3e-8    <- MIXED
```
> **`S_{k+1}^{(q,p)} / S_k^{(q,p)} → q·(p−1)/(p+1)`** — verified across both axes (mixed cases to ~10⁻⁸). Reduces to the
> banked `q/3` at p=2 (q-family) and to `3(p−1)/(p+1)` at q=3 (p-family).

## Collatz is DOUBLY critical
The Plancherel mass **converges (finite S_∞) iff the ratio = 1**, i.e. on the **boundary curve `q(p−1) = p+1`** (equivalently
`q = (p+1)/(p−1)`); off it, `S_k` diverges geometrically and `S_∞ = ∞` (the value must be renormalized, a `c̃(q,p)`).
**Collatz `(q,p)=(3,2)` sits exactly on this curve** (`3·1 = 3 = p+1`), and it is the boundary of *both* one-parameter slices at
once — `q/3 = 1 ⟺ q=3` **and** `3(p−1)/(p+1) = 1 ⟺ p=2`. So `3x+1` with halving is the **unique doubly-critical member** of the
whole (q,p) family: the one point where both the multiplier-family and the valuation-family ratios equal 1. Not a coincidence of
the constant — a structural fact about where the finite-value regime lives.

## What this answers (Q1/Q2/Q3), honestly
- **Q3 (functional equation in the family): YES — but for the associated GRADED, not the value.** `q(p−1)/(p+1)` is a clean
  2-parameter functional law governing the **divergence rate** = the multiplicative / `𝔾_m` (the `D` in `M=D(I+N)`) structure. It
  pins the boundary curve and the leading behavior exactly. It does **not** evaluate the boundary constant: at Collatz the rate is
  1 (marginal), and the *value* still lives in the unipotent `N` part (infinite Mahler depth). The functional equation is the
  gr-level law; the value's functional equation, if any, is the `N`-part *beyond* this — the sharpened pen target.
- **Q1 (finite-depth c̃ vs infinite-depth S_∞ — same web?): different, as expected — it's the extension-vs-graded split.** The
  off-boundary members are divergent; their renormalized `c̃(q,p)` are the "nice"/rational associated-graded (finite-depth); the
  boundary member (Collatz) is the full infinite-depth extension. The relation between them *is* the filtration we already have.
- **Q2 (PSLQ independence of the limit values): infeasible AND moot.** On the boundary we have ~3 digits of `S_∞^{(3,2)}`
  (infinite Mahler order ⟹ no acceleration; ~16 levels/digit), far short of the 30 PSLQ needs; off the boundary the "values" are
  `∞`. So there is no honest digit-hunt to run. Independence of the family is a **Hardouin-criterion PEN result** on the extension
  classes (structural), not a PSLQ on values — exactly as pre-registered.

## Same class (structural), all infinite-depth
Exact `S_k^{(3,p)}` denominators grow doubly-exponentially for **every** p (log_p-den ratio → ~3 > 2, the MAHLER signature):
p=2 (2.77, 3.66, 3.62, 3.47), p=5 (2.53, 5.53, 3.99), p=7 (3.59, 4.21). So all siblings are **infinite Mahler order** — the same
*type* — even though only the boundary member has a finite value.

## Net
- **Unified law `S_{k+1}/S_k → q(p−1)/(p+1)`** (verified both axes + mixed); **boundary curve `q(p−1)=p+1`**; **Collatz is the
  unique doubly-critical point** on it. This is the family's associated-graded functional equation and it hands Wilson's Q3 the
  gr-level law explicitly.
- The **value** at the boundary (S_∞ ≈ 0.475) is untouched — still infinite unipotent/Mahler depth; the gr functional equation
  governs the rate/boundary, not the constant. Same wall, sharper map.
- **Not at stake:** 7/15 (floor 0.473177), MAHLER, GARSIA, DENOM, SOLSTICE, R1–R30. **Newly banked:** the (q,p) Plancherel law +
  Collatz's double-criticality.
