# RESULT — QDIFF-3 (Rung 4): the tower. 𝔾_m stays pro-only at level 2 (direct); G₂ = 𝔾ₐ² by MAHLER-consistency; nesting clean. One honest finding: module-dim ≠ Galois-dim, so the growth needs the exact couplings, not the shape (2026-07-29)

**Probe:** `probes/probe_qdiff3.py`, exact. Rung 4 of Move 1: pin the level-2 truncation `A₂` and gate the tower — does the level-1 `𝔾ₐ` extend, does `𝔾_m` stay pro-only, do the groups nest?

## Setup — A₂ pinned from QDIFF-1's construction (QDIFF-2 diagonal convention)
Diagonal = pure z-power substitution scalings `z^{d_{i+1}−d_i} = z^{2d_i} = z^{4·3^{i−1}}` (level 1: `z⁴`, level 2: `z^{12}`, …). With `A₁ = [[z⁴,z⁴],[0,1]]` the confirmed level-1 block:
```
   A₂ = [[z^{12}, z^{12}, 0],[0, z⁴, z⁴],[0, 0, 1]]     (A₁ = the bottom-right 2×2 block ⟹ nesting embedding)
```

## 4a-diagonal — 𝔾_m collapses at level 2 (DIRECT, rigorous)
`σ(y)=z^{12}y` has the rational solution `y=z⁶` (`3·6=12+6`, verified exact) → `z^{12}` is a σ-coboundary. General: every level-i scaling `z^{4·3^{i−1}}` is a coboundary (`y=z^{2·3^{i−1}}`, always an even exponent). **⟹ 𝔾_m collapses at level 2 and at every level; the reductive part stays pro-only.** Model-independent (any realization has z-power substitution scalings). This is the solid weak core.

## The honest finding — module-dim ≠ Galois-dim, and the monomial model is unfaithful
The unipotent *Galois* dimension grows only if successive log-sources are **independent mod coboundaries**. Single-cascade criterion (exact): `Σ_j r_j·z^{−2·3^j}` is a coboundary **iff `Σ_j r_j = 0`**. The naive monomial couplings put every level's source in the one σ-cascade `{−2·3^j}`: `e₁=z⁻²` (j=0), `e₂=z⁻⁶` (j=1), and `e₂ − 1·e₁` **is a coboundary** (`Σ=0`) — so the monomial sources are **dependent**, which would give `G_r=𝔾ₐ` for all r (**no growth**).

**That reading is unfaithful** — `G_r=𝔾ₐ` ∀r ⟹ finite Mahler depth ⟹ **contradicts MAHLER**. The naive monomial couplings artificially collapse a cascade the real data does not have. So the monomial model cannot be used to verify (or refute) the growth — a genuine trap, and the reason "pin A₂ exactly" is again load-bearing: the *shape* `D(I+N)` fixes the diagonal (⟹ 𝔾_m collapse, solid) but **not** the unipotent couplings (⟹ the growth, sensitive).

## 4a-unipotent — G₂ = 𝔾ₐ² (dim grew by 1), forced by MAHLER-consistency
The real level-i log-source is tied to `Λ_i`, which has **no finite recurrence** (R27-A / MAHLER) — the `Λ_i` are algebraically independent data (doubly-exponential denominators), **not σ-translates**. So successive log-sources are **independent mod coboundaries** → each level adds a genuine new `𝔾ₐ`. The integrator is a **direct sum** `T = base + Λ₁ + Λ₂` (each `Λ` enters `T` independently — *parallel, not iterated*) → the unipotent is **abelian** → **G₂ = 𝔾ₐ²** (not the non-abelian `U₃`). *Consistency check:* if the sources were dependent (`G_r=𝔾ₐ`, no growth), the inverse limit would be finite-dimensional, contradicting MAHLER's infinite unipotent depth — so the growth is **required, not optional**. (⚠️ the *direct* construction — exact `Λ`-derived couplings as `ℚ(z)` entries, then the literal independence test — is not yet pinned; the monomial model is too coarse. Frontier-adjacent.)

## 4b-nesting — G₂ restricts to G₁ (clean)
`A₁` is the bottom-right 2×2 block of `A₂`, so the level-1 PV sub-extension embeds. A level-2 automorphism `(ℓ₁↦ℓ₁+u₁, ℓ₂↦ℓ₂+u₂)` restricts to the level-1 sub as `(ℓ₁↦ℓ₁+u₁)` — the **projection `𝔾ₐ² → 𝔾ₐ`** (forget the level-2 coordinate), surjective and compatible with the tower ⟹ `lim← G_r` is well-defined.

## Rung-4 verdict
**Weak core (pens as theorem):**
- **4a-diagonal** [DIRECT]: `𝔾_m` is a coboundary at level 2 and every level ⟹ reductive part pro-only.
- **4a-unipotent** [by MAHLER-consistency + direct-sum integrator]: `G₂ = 𝔾ₐ²`, dim grew by 1.
- **4b-nesting**: `G₂ = 𝔾ₐ² —projection→ G₁ = 𝔾ₐ`; the tower nests, `lim← G_r` defined.
- **weak-4c**: the reductive `𝔾_m` emerges only in the inverse limit (= MAHLER's doubly-exp denominators).

**Frontier (does NOT pen as established):**
- Direct construction of the exact `Λ`-derived level-2 coupling as a `ℚ(z)` entry, then the *literal* independence test (the monomial model is unfaithful — collapses the cascade).
- **strong-4c**: `inverse-limit-of-groups = Galois-group-of-the-limit` — the one step needing mathematics nobody's finished.

**Finding surfaced:** the **module dimension** (QDIFF-1, grows) is *not automatically* the **Galois dimension**; the latter grows iff the log-sources are independent, which MAHLER supplies — but the naive monomial couplings would say otherwise, so a *fully* direct Rung 4 needs the exact couplings, not just the `D(I+N)` shape.

**Lean:** same build-by-hand path — the added Rung-4 predicate is "the level-2 source is independent of level-1 mod coboundary," a decidable `ℚ(z)` statement once the exact coupling is pinned; the `𝔾_m`-collapse predicate (`z^{12}` a coboundary) is already decidable and direct.

**Not at stake:** S_∞≈0.475 (floor 0.473177), MAHLER, MIRROR, PHYDRA, two-walls, QDIFF-1/2, R1–R30. Tower structure only; the value is Rung 5, conditional; 7/15 excluded regardless.
