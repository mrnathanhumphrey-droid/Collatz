# QSC_MOMENT_PREDICTIONS — predicted vs measured Syracuse moments

**Date:** 2026-05-15
**Mode:** E.
**Companion:** `QSC_VERBATIM.md`, `QSC_SYRACUSE_IDENTIFICATION.md`, `QSC_DISPOSITION.md`.

---

## 0. Setup

Syracuse measured rows under Reading B / marginal centering / V_TRUNC = 16 / sum_entries reduction (`D1_DISPOSITION.md`):

- Row (b): `ϕ(X̃_{j_1} · X̃_{j_2})` ≈ 0 (numerically 1.076×10⁻⁷, noise floor)
- Row (d): `ϕ(X̃_{j_1} · X̃_{j_2} · X̃_{j_1})` = 1.078×10⁻¹ (4–7 orders above noise across all 4 scalar reductions)
- Row (f): `ϕ(X̃_{j_1} · X̃_{j_2} · X̃_{j_1} · X̃_{j_2})` = 6.089×10⁻¹

Fubini inner factor `F(v_1, v_1') = E_{(v_2)}[X̃_{j_2} · X̃_{j_1} · X̃_{j_2}] = 6.347×10⁻²` (constant across 12 grid points).

We tabulate predictions under each plausible QSC identification.

---

## 1. Identification A — X̃_j ↔ dA_j + dA_j† (Hermitian quadrature, "Brownian-noise" analog)

**Filtration:** Syracuse `B_j` ↔ HP `F_t = Γ_{(t]}`. (NB: mismatch in §1 of QSC_SYRACUSE_IDENTIFICATION.md — Syracuse B abelian, HP F non-abelian. Force the identification anyway and check moments.)

**Increment:** X̃_j is interpreted as the Hermitian quadrature `Q_j := dA_j + dA_j†` (where dA_j, dA_j† are discrete-time increments at "time" j with `dA_j · dA_j† = δ_{j,j} · 1`, the only non-zero Itô product within this restricted table).

### Row (b)

`ϕ(Q_{j_1} Q_{j_2}) = ϕ((dA_{j_1} + dA_{j_1}†)(dA_{j_2} + dA_{j_2}†))` for j_1 ≠ j_2.

Expansion: `dA_{j_1} dA_{j_2} + dA_{j_1} dA_{j_2}† + dA_{j_1}† dA_{j_2} + dA_{j_1}† dA_{j_2}†`. Under the Itô table at distinct time slices, all four products give 0 (no coincident slices → no `dt` contribution). Vacuum expectation = 0.

**Prediction (b):** 0. **Match Syracuse 0:** ✓

### Row (d)

`ϕ(Q_{j_1} Q_{j_2} Q_{j_1})` with j_1 ≠ j_2.

Expand into 8 monomials. Each has the form `dA^{(ε_1)}_{j_1} dA^{(ε_2)}_{j_2} dA^{(ε_3)}_{j_1}` with ε_i ∈ {0, *} (creation/annihilation). By time-locality of the Itô table, the only non-trivial reductions happen at COINCIDING time slices. Here `j_1` appears twice (positions 1 and 3) but the j_2 in the middle is not coincident with either, so we can't reduce `dA^{ε_1}_{j_1}` with `dA^{ε_3}_{j_1}` directly past `dA^{ε_2}_{j_2}` — they're separated by a different-slice operator.

In HP / AP, distinct-slice operators COMMUTE (the filtration tensor splits as `Γ_{(t]} ⊗ Γ_{(t}`). So we can re-order: `dA^{ε_1}_{j_1} · dA^{ε_3}_{j_1} · dA^{ε_2}_{j_2}` (modulo sign — for bosonic, the order is preserved on adapted-process commuting). Now the j_1-j_1 product reduces by Itô. Then the j_2-factor's vacuum expectation is 0 (since each of `dA_{j_2}, dA_{j_2}†` has vacuum expectation 0 — `dA` annihilates vacuum, `dA†` creates a 1-particle state which is orthogonal to vacuum).

**Prediction (d):** 0. **Match Syracuse 0.108:** ✗ MISMATCH BY 5+ ORDERS OF MAGNITUDE.

### Row (f)

`ϕ(Q_{j_1} Q_{j_2} Q_{j_1} Q_{j_2})` with j_1 ≠ j_2.

By commuting same-slice factors together (since distinct slices commute):

`= ϕ(Q_{j_1}^2 · Q_{j_2}^2)` (after re-ordering)

`Q_{j_1}^2 = (dA + dA†)^2 = dA^2 + dA dA† + dA† dA + (dA†)^2`. By HP Itô: `dA · dA† = dt` (the only non-zero), all others 0. So `Q_{j_1}^2 = dt_{j_1}` (using discrete-time Itô analog `dA_j · dA_j† = δ_{j_1} dt_{j_1}` = the indicator at slice j_1).

`Q_{j_2}^2 = dt_{j_2}` likewise.

`Q_{j_1}^2 Q_{j_2}^2 = dt_{j_1} · dt_{j_2}`. In discrete-time interpretation, this is the product of two slice indicators at distinct slices, which equals... the indicator at both slices, with value `h · h = h²` for time-step h. Taking vacuum expectation: `dt_{j_1} dt_{j_2} → 0` as h → 0 (a second-order infinitesimal).

In the discrete (AP) version: `dt_{j_1}` is the constant operator `h · I` at slice j_1, and `dt_{j_1} · dt_{j_2} = h² · I`. Vacuum expectation `h²`.

**Prediction (f):** O(h²) → 0 as h → 0; in the discrete-time fixed-h setting, the value is `h²`. **Match Syracuse 0.609:** ✗ MISMATCH (the second-order infinitesimal nature of `Q^2 Q^2` doesn't produce a fixed O(1) value).

But there's a subtlety. Under non-commutative re-ordering at the SAME slice, `Q_{j_1} Q_{j_2} Q_{j_1} Q_{j_2}` is NOT the same as `Q_{j_1}^2 Q_{j_2}^2`. It's:

`Q_{j_1} Q_{j_2} Q_{j_1} Q_{j_2}`. Since j_1 ≠ j_2 and distinct slices commute, we CAN re-order to `Q_{j_1} Q_{j_1} Q_{j_2} Q_{j_2} = Q_{j_1}^2 Q_{j_2}^2`. Same result.

**Prediction (f) under Identification A:** O(h²). **Mismatch** Syracuse's 0.609.

### Identification A summary

| Row | Syracuse | Predicted (A) | Match |
|---|---|---|---|
| (b) | 0 | 0 | ✓ |
| (d) | 0.108 | 0 | ✗ |
| (f) | 0.609 | O(h²) → 0 | ✗ |

**Identification A FAILS at rows (d) and (f).** Standard HP quadrature does not reproduce Syracuse's non-zero cross-step moments.

---

## 2. Identification B — X̃_j ↔ dΛ_j (number / gauge)

**Increment:** X̃_j ↔ `dΛ_j = a_j^{11}` in AP discrete-time (counting / number differential at slice j).

### Row (b)

`ϕ(dΛ_{j_1} dΛ_{j_2})` for j_1 ≠ j_2. Distinct slices commute and have independent vacuum expectations. `ϕ(dΛ_j) = ϕ(a_j^{11}) = ⟨Ω, a_j^{11} Ω⟩ = ⟨e_0, |e_1⟩⟨e_1| e_0⟩_j = 0`. So `ϕ(dΛ_{j_1} dΛ_{j_2}) = ϕ(dΛ_{j_1}) · ϕ(dΛ_{j_2}) = 0`.

**Prediction (b):** 0. **Match:** ✓

### Row (d)

`ϕ(dΛ_{j_1} dΛ_{j_2} dΛ_{j_1})`. By distinct-slice commuting, `= ϕ(dΛ_{j_1}^2 · dΛ_{j_2})`. By matrix-unit composition `a^{11} · a^{11} = a^{11}`, so `dΛ^2_{j_1} = dΛ_{j_1}`. Then `ϕ(dΛ_{j_1} · dΛ_{j_2}) = 0` per row (b).

**Prediction (d):** 0. **Match Syracuse 0.108:** ✗

### Row (f)

`ϕ(dΛ_{j_1}^2 · dΛ_{j_2}^2) = ϕ(dΛ_{j_1} · dΛ_{j_2}) = 0`.

**Prediction (f):** 0. **Match Syracuse 0.609:** ✗

### Identification B summary

| Row | Syracuse | Predicted (B) | Match |
|---|---|---|---|
| (b) | 0 | 0 | ✓ |
| (d) | 0.108 | 0 | ✗ |
| (f) | 0.609 | 0 | ✗ |

**Identification B FAILS at rows (d) and (f).** Same structural reason as A: time-local Itô table + vacuum expectation kills cross-step moments.

---

## 3. Identification C — X̃_j ↔ adapted-process F_j · (quantum increment)

**Setup:** Try `X̃_j = F_j · dQ_j` where `F_j ∈ B_{j-1}` (the prior-measurable adapted integrand) and `dQ_j` is a "pure quantum increment" at slice j (an HP differential).

For Syracuse, `F_j` could encode the χ_j phase factor (with `b_{[1, j-1]}`-dependence) and `dQ_j` the within-pair (v_{2j-1}, v_{2j}) randomness.

### Structural issue (per §3 of QSC_SYRACUSE_IDENTIFICATION.md)

The Syracuse off-diagonal operator:

> `Off_j(ξ) = Σ_{v ≠ v'} 2^{-v} 2^{-v'} · χ_j(b_{[1, j-1]}; v, v', ξ) · σ_{-(v + v')}`

does NOT factor as `F(b_{[1, j-1]}) · G(v, v', ξ)` because χ_j has `b_{[1, j-1]}` and `(v, v')` coupled multiplicatively inside the exponential. The B-measurable content is mixed inside the Σ_{v, v'}.

Even if we forced an identification `X̃_j = F_j · dQ_j`, the F_j and dQ_j would themselves be sums of products, and the moment calculation would not simplify cleanly.

### Best-case prediction

If `F_j` and `dQ_j` were cleanly separable, then `ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1}) = ϕ(F_{j_1} dQ_{j_1} F_{j_2} dQ_{j_2} F_{j_1} dQ_{j_1})`. Moving the B-measurable F's around via B-linearity of E_B, and using independence of dQ_j at distinct slices, this STILL gives 0 (because dQ_{j_2} alone has E[dQ_{j_2}] = 0 from centering).

**Prediction (d) under Identification C (forced):** 0. **Match Syracuse 0.108:** ✗

### Identification C summary

Even with the adapted-integrand structure, **vacuum expectation of the middle dQ_{j_2} kills the row (d) prediction**. The fundamental issue: in any QSC framework where the increments are time-local AND vacuum-centered, the cross-step triple product `X_{j_1} X_{j_2} X_{j_1}` vanishes because of the centered middle factor.

To get a NON-ZERO row (d), Syracuse needs the middle X̃_{j_2} to NOT vanish in expectation conditional on the surrounding context. The way this works in Syracuse: the surrounding X̃_{j_1}'s have "imprinted" the b_{[1, j_2]}-accumulator data into the trace structure, and the X̃_{j_2}-expectation conditional on this is NOT zero (it picks up the phase χ_{j_2}(b_{[1, j_2 - 1]}) which is non-trivially imprinted by the prior X̃_{j_1}'s).

This is a level-graded conditional expectation behavior that QSC's standard adapted-process integrand structure does NOT encode.

---

## 4. The Fubini constant factor — what would QSC predict?

Syracuse measures `F(v_1, v_1') = E_{(v_2)}[X̃_{j_2} · X̃_{j_1} · X̃_{j_2}] = 6.347×10⁻²` CONSTANT in (v_1, v_1').

### QSC prediction for analogous Fubini object

Under any of identifications A, B, C, the inner expectation `E_{(v_2)}[X̃_{j_2} X̃_{j_1} X̃_{j_2}]` involves the squared X̃_{j_2} sandwiching X̃_{j_1}. In HP/AP, squared annihilation `dA_{j_2}^2 = 0`, squared creation `(dA_{j_2}†)^2 = 0`, squared number `dΛ_{j_2}^2 = dΛ_{j_2}`. The only non-vanishing squared product is `dA_{j_2} · dA_{j_2}† = dt_{j_2}`.

If we identify X̃_{j_2} = `dA_{j_2} + dA_{j_2}†`, then `X̃_{j_2}^2 = dA dA† + dA† dA = dt + 0 = dt` (using HP Itô). So `X̃_{j_2} X̃_{j_1} X̃_{j_2}` has structure that involves `dt_{j_2}` (a scalar in the j_2-slice) sandwiching X̃_{j_1}. After taking E_{(v_2)} (analogous to integrating over the j_2-fiber), this reduces to `c · X̃_{j_1}` for some scalar c.

**QSC prediction for F(v_1, v_1'):** `c · X̃_{j_1}` after E_{(v_2)} — this IS a function of (v_1, v_1') via the X̃_{j_1} factor, so it's NOT constant in (v_1, v_1').

**Syracuse measures F constant in (v_1, v_1').** This is the STRONGER claim:

> The inner factor `E_{(v_2)}[X̃_{j_2} X̃_{j_1} X̃_{j_2}]` is INDEPENDENT of the j_1-fiber data — purely a scalar at the (j_1)-position.

QSC's adapted-process structure would predict the inner factor depends on the j_1-fiber via the X̃_{j_1} content. Syracuse's CONSTANCY means the j_1-fiber dependence has been COMPLETELY ABSORBED into the inner-pair integration at j_2. This is a TIGHTER algebraic property than QSC delivers.

### Mode-E note

The constancy of F(v_1, v_1') is "structural cleanness beyond QSC Itô" — Syracuse's inner factor is **more constrained** than what HP or AP would naturally produce. This is consistent with the diagnosis that Syracuse needs a NEW framework (Outcome C of the brief), not a QSC fit.

---

## 5. Summary table

| Framework / identification | Row (b) match | Row (d) match | Row (f) match | F = const match |
|---|---|---|---|---|
| HP quadrature `dA + dA†` | ✓ (0) | ✗ (0 vs 0.108) | ✗ (h² vs 0.609) | ✗ |
| HP number `dΛ` | ✓ (0) | ✗ (0 vs 0.108) | ✗ (0 vs 0.609) | ✗ |
| AP matrix-unit `a_j^{ij}` (any single) | ✓ (0) | ✗ (0 vs 0.108) | ✗ (0 vs 0.609) | ✗ |
| Adapted integrand `F_j · dQ_j` | ✓ (0) | ✗ (0 vs 0.108) | ✗ (0 vs 0.609) | ✗ |

**No QSC identification matches Syracuse beyond row (b) (which is also trivially matched by every block-factorization framework already eliminated in D1/D2).**

The row (d) and row (f) non-zero values are structurally INCOMPATIBLE with QSC's time-local Itô calculus + vacuum-centered increments.

---

## 6. Mode-E gap

The above predictions assume the "natural" / canonical identifications. A bespoke QSC-like framework that adds **non-time-local Itô products driven by B-measurable phase coupling** could in principle reproduce Syracuse's moments — but this is not HP, not AP, not Köstler-Speicher. It would be a NEW framework that BORROWS the adapted-process language of QSC but rewrites the Itô table.

**See `QSC_DISPOSITION.md` for the verdict.**
