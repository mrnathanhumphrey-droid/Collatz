# RESULT — DENOM OBSTRUCTION: the denominator theorem CANNOT decide the arithmetic type of S∞ (2026-07-27)

**Pen (Wilson) + numerical gate (Claude).** Does the denominator/norm structure of the S-ladder obstruct S∞ from being
rational (or algebraic)? Verdict: **no — provably, and for a precise reason.** A clean, honest closure of the arithmetic route.

## The S-ladder denominators (banked)
`S_r ∈ ℚ` at every level, denominator controlled by `2^M − 1`, `M = 2·3^{r-1}`, up to powers of 3:
- `S₁ = 2/3` (den 3; `2²−1 = 3`).
- `S₂ = 10/21 = 10/(3·7)` (den 21 | `3ᵃ·(2⁶−1) = 3ᵃ·63`). ✓ gated.
- den `S_r | 3ᵃ·(2^{2·3^{r-1}} − 1)`, so **`log₂ D_r ~ M = 2·3^{r-1}` — doubly exponential in r.**
- ⚠️Zsygmondy footnote: `2⁶−1 = 3²·7` is a Zsygmondy *exception* (no primitive prime of order 6; `ord₂(7)=3`), so
  the "new prime has `ord₂ = 2·3^{r-1}`" attribution wrinkles at r=2; it does **not** affect `D_r`'s size. Resumes at
  r=3 (19 is primitive-for-18).

## The obstruction runs the WRONG way (the finding)
Suppose `S∞ = a/b` (lowest terms). Each `a/b − S_r` is a nonzero rational with denominator | `b·D_r`, so by Liouville
`|a/b − S_r| ≥ 1/(b·D_r)`. The measured convergence is geometric, `|S∞ − S_r| ~ C·(0.867)^r`. Combined:
```
1/(b·D_r) ≤ C·(0.867)^r   ⟺   den_bits(r) ≥ conv_bits(r)
den_bits  ~ M = 2·3^{r-1}   (doubly exponential)
conv_bits ~ r·log₂(1/0.867) = 0.206·r   (LINEAR)
```
Gate (r=1..7): gap `den_bits − conv_bits` = **1.8, 5.6, 17.4, 53.2, 161, 485, 1457 bits** — an exploding margin.
**The inequality is satisfied with enormous room ⟹ NO obstruction to rationality.** Same for algebraicity: Roth needs
`(0.867)^r ≳ c/D_r^{2+ε} ~ 2^{−(2+ε)·2·3^{r-1}}`, again satisfied trivially. A Diophantine obstruction bites only when the
*approximation outruns the denominator*; here it is the reverse by a doubly-exponential margin. The `S_r` are so inefficient
as rational approximants (huge denominator, only geometric accuracy) that they carry **no** Diophantine information about S∞.

## Reduction (the honest endpoint)
The denominator/norm structure **cannot** decide whether S∞ is rational, algebraic, or transcendental — the finite-level
denominators are far too large, too fast, to constrain a geometrically-approached limit. So the arithmetic type is **not
accessible from the arithmetic we built**; it reduces to the *analytic* type of ν:
> **S∞ rational is tied to ν being absolutely continuous; ν singular ⟹ S∞ is the Plancherel mass of a singular measure,
> generically irrational (and 7/15 "wrong in kind").** ("tied to / generically" — a strong localization, not a proven iff.)
And ν's singular-vs-a.c. type is the **Bernoulli-convolution question proper** (1/λ = 2 is a Pisot base ⟹ singular is the
default expectation) — a named hard problem, with no special structure on our side that shortcuts it.

## Net (three routes now closed, cleanly)
- **q-family isolated** (q=3 the unique nontrivial contracting odd integer — no neighbors to interpolate). [result: q-fam]
- **base-3 lattice nonexistent** (two incommensurate scales 2,3; oscillation quasi-periodic). [`result_LATTICE.md`]
- **denominator/norm route cannot decide arithmetic type** (this result — doubly-exp denominators vs geometric convergence).

The sign question reduces, with **no further leverage available to us**, to the singular-vs-a.c. type of ν — hard,
field-recognized, not shortcuttable from here. Value stands: **S∞ ≈ 0.475, exact floor `2·T_20 = 0.473177`**, 7/15 excluded
barring a hidden `ρ₃ > 0.999` real mode. The **denominator theorem is a genuine standalone result** (the S-ladder denominator
structure + this "cannot decide arithmetic type" corollary); Wilson writes it up. Not at stake: P6D–P6K identities,
S_{i+1}=2T_i, R1–R30, the i=20 no-crossing observation.
