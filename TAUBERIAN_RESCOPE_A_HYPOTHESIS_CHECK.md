# TAUBERIAN_RESCOPE_A_HYPOTHESIS_CHECK (FS Ch. VI × inputs (1)-(4))

**Date:** 2026-05-13.

Read against Theorem VI.4 (single singularity) and VI.5 (multiple singularities).

---

## h × I matrix — Theorem VI.4 (single dominant singularity)

Reading f(z) = Σ ε_n z^n.

| Hypothesis | (1) ε_k | (2) C1 | (3) C2 | (4) BT |
|---|---|---|---|---|
| h_1: f analytic at 0 | SATISFIED (8 nonzero Taylor coeffs, finite radius) | N/A | N/A | N/A |
| h_2: f has singularity at ζ on circle of convergence | UNVERIFIABLE — need full ε_n sequence; ratio test on first 8 ratios gives ratios ≈ 0.43-0.53 for k=2..6, suggesting radius ≈ 2, but k=7 jump (ratio 2.36) breaks the pattern. ζ unknown. | N/A | N/A | N/A |
| h_3: analytic continuation to ζ·Δ_0 | UNVERIFIABLE — Δ-analyticity requires proof of continuation past |z| = radius; not available from inputs. | N/A | N/A | N/A |
| h_4: singular expansion f(z) = σ(z/ζ) + O(τ(z/ζ)) in S = {(1-z)^{-α} λ(z)^β} | UNVERIFIABLE — empirical fit gives α-1 ≈ -3.35 for k=3..6 (i.e. α ≈ -2.35), but k=7,8 breaks this. The standard scale (1-z)^{-α} predicts coefficient n^{α-1}; observed n^{-3.35} decay then resurgence at k=7 is incompatible with a single (1-z)^{-α} term. | N/A | N/A | N/A |

**Theorem VI.4 disposition: BLOCKER** — all the load-bearing hypotheses h_2, h_3, h_4 require analytic-continuation knowledge of g that the 8-coefficient input does not supply.

---

## h × I matrix — Theorem VI.5 (multiple dominant singularities)

This is more interesting: the k=7 jump in |ε_k|·2^k is empirically consistent with **two** (or more) dominant singularities ζ_1, ζ_2 on the boundary of convergence, each contributing a term to the asymptotic, with the slower-decay term taking over at higher k.

| Hypothesis | (1) ε_k | (2) C1 | (3) C2 | (4) BT |
|---|---|---|---|---|
| h_1: f analytic in |z| < ρ with finitely many dominant singularities on |z| = ρ | UNVERIFIABLE — would need to identify the multiple singularities. | N/A | N/A | (4) ARCHIMEDEAN finding bears on this: if the c=7/45 closure requires *adelic* substrate (BT disposition), then the dominant-singularities structure of any single p-adic generating function will NOT capture the full closure target. The number of dominant singularities on the unit circle would need to encode adelic data. |
| h_2-h_4: Δ-analyticity through indented disc with cones at each ζ_j | UNVERIFIABLE | N/A | N/A | N/A |

**Theorem VI.5 disposition: BLOCKER** — same load-bearing failures as VI.4, compounded by the unknown number/locations of dominant singularities. Input (4) actively suggests this is the *wrong category* of object (need adelic, not single-place generating series).

---

## Aggregate disposition for A: BLOCKER

The Flajolet-Sedgewick singularity-analysis machinery requires:
- An explicit generating function (or at least its analytic continuation), AND
- Knowledge of the location and nature of dominant singularities on the boundary of convergence.

Inputs (1)-(4) supply neither. We have 8 numerical coefficients with an empirical 2-regime pattern (geometric-like decay for k=2..6, jump at k=7,8), and structural descriptions of the substrate (renewal-walk, BMP support, archimedean place) that do not specify the singular structure of any generating series.

**The empirical 2-regime pattern in (1) is consistent with Theorem VI.5's multi-singularity setup, but identifying the specific σ_j(z), τ_j(z) for each ζ_j requires either analytic theory we don't have or many more coefficients.**

Disposition: **BLOCKER**.
