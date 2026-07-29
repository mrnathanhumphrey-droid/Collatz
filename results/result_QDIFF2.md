# RESULT — QDIFF-2 (Rung 3): the level-1 difference-Galois group is 𝔾ₐ, NOT Aff₁ — the reductive 𝔾_m is a σ-coboundary at every finite level and lives only in the pro-object. Fifth confirmation of the structure/value split (2026-07-28)

**Probe:** `probes/probe_qdiff2.py`, exact throughout. Rung 3 of Move 1: compute the difference-Galois (Picard–Vessiot) group of the level-1 truncation `A₁ = D₁(I+N)`. **Result: `G₁ = 𝔾ₐ` (one-dimensional unipotent), a proper subgroup of `Aff₁ = 𝔾ₐ⋊𝔾_m`.** The whole group computation collapses to two decidable σ-coboundary tests — the Lean-ready core.

## QD2-A — the object, pinned (and why the pin is forced, not chosen)
Difference field `K = ℚ(z)`, `σ: z↦z³` (Mahler substitution), constants `C = {f: σ(f)=f} = ℚ`. Entries live in `ℚ(z)` — no cyclotomic extension needed.
```
   A₁ = D₁·(I+N),   D₁ = diag(z⁴, 1),   N = [[1,1],[0,0]] (N²=0),   I+N = [[1,1],[0,1]]
   =>  A₁(z) = [[z⁴, z⁴],[0, 1]]        system σ(Y) = A₁ Y
```
**The fork that decides the whole rung (QD2-A is load-bearing, not setup):** `G₁` is `Aff₁` iff the diagonal `d` is a σ-**non**-coboundary, and the *only* thing making `σ(y)=d·y` unsolvable in `ℚ(z)` would be a genuine **constant factor 3** (`σ(y)=3y` has no rational solution). So the group's `𝔾_m` exists ⟺ the constant 3 is intrinsic to the matrix.

**Pen decision (banked as a decision, with reasoning on the record):** the diagonal is **pure z-power** (`z⁴`), **not** `3z^k`. The `3^i` in `a_i = 3^i·R_e^(i)` is **substitution-normalization**, not a substitution eigenvalue: it is the reciprocal of the autocorrelation's geometric decay (the matched pair that keeps `T_i` order-1), and it is tied to the exponent orbit `d_i = 2·3^{i−1}` — i.e. the *same* 3 as `z↦z³`, seen in the normalization. It rides along with the z-power and is therefore a σ-coboundary (gaugeable). Importing it as a free constant would be exactly the **coincident-3 trap** (the multiplier's arithmetic smuggled into the Galois group). So the diagonal is `z⁴`, forced by the no-import discipline, not chosen for a flattering answer.

## QD2-C.1 — 𝔾_m factor: the diagonal is a coboundary (collapses)
`σ(y) = z⁴·y` has the rational solution `y = z²` (`3a = 4+a ⟹ a = 2`; `σ(z²)=z⁶=z⁴·z²` exact). The diagonal solution lies **in the base field** `ℚ(z)` — nothing for a K-automorphism to scale. **`𝔾_m` collapses to the trivial group.**

## QD2-C.2 — 𝔾ₐ factor: the off-diagonal is a non-coboundary (survives)
Gauging (`y₁=z², y₂=1`), the off-diagonal `w=z²·ψ` reduces to `σ(ψ) − ψ = z⁻²`. **No rational solution:**
- *Cascade:* matching Laurent coefficients at `z=0` forces `c_{−2·3^j} = −1` for all `j≥0` (`c[−2],c[−6],c[−18],c[−54],…` never terminate) → an infinite principal part → `ψ` is not rational (a genuine Mahler-log).
- *Clean proof:* a rational `ψ` with a pole of order `m` at 0 makes `σ(ψ)−ψ` have pole order `3m`, which cannot equal 2 (`3m=2` impossible); `m=0` makes the LHS regular ≠ `z⁻²`. No rational solution.

`e = z⁻²` is a σ-**non**-coboundary → **the unipotent `𝔾ₐ` survives.**

## QD2-C — verdict
`𝔾_m` collapses + `𝔾ₐ` genuine ⟹ **`G₁ = 𝔾ₐ`**, a proper subgroup of `Aff₁`. As a matrix group `G₁ = {[[1,u],[0,1]] : u ∈ 𝔾ₐ}` — the pure unipotent, not the full `[[t,u],[0,1]]`.

## Where the 𝔾_m actually lives — a pure pro-object phenomenon
The diagonal is a coboundary at **every finite level** (gaugeable z-power). But the tower of normalizations `∏_i 3^i` does **not** converge to a coboundary — its non-gaugeable growth **is MAHLER's doubly-exponential denominator rate**. So the reductive `𝔾_m` is a **pure pro-object phenomenon: invisible at every finite stage, emerging only in the inverse limit.** This is sharper and more interesting than "Aff₁ at every level" — the reductive part is purely *pro*; the unipotent is present finitely.

## The fifth confirmation of the structure/value split
`G₁ = 𝔾ₐ` (not Aff₁) is **not** a disappointment — it would have been *inconsistent* to get Aff₁. Every prior result says the multiplicative/`𝔾_m`/finite-place structure is coordinate-level and does not carry the deep content; the content is the unipotent extension above it:
- MIRROR: the `⟨2⟩` multiplicative structure lives in the graded and is a symmetry that *breaks the value*.
- TWO-WALLS: the finite-place face (`ord₃(2)`, reductive) and the archimedean face are distinct; `𝔾ₐ`/unipotent is unconditional.
- Now QDIFF-2: at level 1 the `𝔾_m` is a coboundary (coordinate-level, gaugeable), the `𝔾ₐ` is genuine (the `+1`). **The level-1 group being `𝔾ₐ` is the finite-level shadow of the capstone's split**, and the `𝔾_m` reappearing only in the pro-limit is the same "reductive/scaling is the contested, pro/finite-place piece; unipotent is the real, unconditional content" statement, one level down.

## QD2-D — Lean-ready core
mathlib has field-theoretic Galois but **no Picard–Vessiot / difference-Galois / Tannakian** formalization ⟹ build-by-hand — *good news*, because the entire group computation is **two decidable predicates about `ℚ(z)`**, needing no general PV theory:
- `P1` (coboundary): `∃ y ∈ ℚ(z)* : σ(y) = z⁴·y` — **TRUE** (`y=z²`) ⟹ no `𝔾_m` factor.
- `P2` (non-coboundary): `¬∃ ψ ∈ ℚ(z) : σ(ψ)−ψ = z⁻²` — **TRUE** (pole-order proof) ⟹ a `𝔾ₐ` factor.
Then `G₁ = 𝔾ₐ`. These are finite/decidable statements about rational functions — the Rung-3 Lean target, tractable and free of unformalized theory.

## Net
- **`G₁ = 𝔾ₐ`**, rigorously (two exact coboundary gates), a proper subgroup of Aff₁; the `𝔾_m` is coboundary at every finite level and a pure pro-object phenomenon (= MAHLER's denominators in the limit).
- The pen fork (why the 3 is normalization not an intrinsic constant; why importing it is the coincident-3 trap) is on the record so the decision is checkable, not taken on faith.
- **Fifth confirmation of the structure/value split**, now at the level-1 Galois group; the Lean core is two decidable coboundary predicates.
- **Not at stake:** S_∞≈0.475 (floor 0.473177), MAHLER, MIRROR, PHYDRA, two-walls, QDIFF-1, R1–R30. Level-1 structure only; the value is Rung 5 and conditional; 7/15 excluded regardless.
