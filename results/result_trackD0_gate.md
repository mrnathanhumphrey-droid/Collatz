# GATE G-D0 — M_tower as a gauge-fibered weighted shift (Wilson's D-0 frozen). ALL 5 POINTS PASS, L=2,3.

**Date:** 2026-07-18  Probe `probes/probe_trackD0_gate.py`. Direct/exact, exhaustive at L=2,3. Gates the
frozen decomposition M_tower = S + F, the three new lemmas (D0.1–D0.3), and the four specializations.

## 1. KERNEL IDENTITY — EXACT
The (u′, s | tape)-form carry **γ′ = ⌊γ/3⌋ + T̂ + c** (T̂=(T−T₀)/3, c=(d₀+T₀)/3) rebuilds M_tower
**entry-identical** to build_M_gen: max|diff| = 0.0, nnz equal (3240 @L2, 892296 @L3).

## 2. LEMMA D0.1 (carry-bit law) — HOLDS
**c = 1_{v₃(γ)=0}** on every surviving branch: **0 violations**, exhaustive, both L. The depth class *is* the
carry-in bit (unit-carry always injects c=1; divisible never, c=0).

## 3. LEMMA D0.2 (position-space triangularity) — HOLDS
Exhaustive at L=3, j=0,1,2: with (γ mod 3^{j+1}, u′ mod 3^{j+1}, e′ mod 2·3^j) fixed, varying the higher
digits leaves **γ′ mod 3^j invariant** — 0 violations at every j. The addition is depth-triangular; the
dynamics never reads upward. (This is the character-side rule 3^{L−1−j}|κ transported to position space.)

## 4. SPECIALIZATIONS — byte-agreement with the banked tables
- **4a cell row-sums** = {2/9, 5/18, 5/9, 4/9} (exact set match, both L).
- **4b cascade marginal** = 2·3^{−(j+1)} + tail: L=3 → {2/3, 2/9, 1/9} exact; L=2 → {2/3, 1/3} (the v=L−1
  truncation folds the 1/9 tail — expected, finite-section at the tape top).
- **4c collective 2×2** = (1/27)[[5,4],[4,5]]: the source-uniform **4×4 cell transfer** is block-anti-diagonal
  (the v0↔odd / v≥1↔even parity-flip 2-cycle),
  `27·T = [[0,0,4,2],[10,5,0,0],[0,0,5,5/2],[8,4,0,0]]`, with spectrum **{1/3, 1/27, 0, 0}** (EXACT, both L,
  byte-identical). Its non-null core has trace 10/27 and det 1/81 = exactly the invariants of (1/27)[[5,4],[4,5]].
- **4d k=0 kernel = E-FORM**: the 2c0-G2 sector form B̂[k_out,k_in]=R_{k_in}(s)·N_{k_in−k_out}/D closes (k=0
  slice = R₀(s)·N₀/D = the E-FORM chain).

## 5. D0.3 (critical weighting) — HOLDS, and the smoke test is decisive
**ρ(S) = 1/3 to machine precision, both L.** Subcriticality smoke test (pre-registered, ALGEBRAIC):

| | λ=0.4 | λ=0.5 | λ=0.6 |
|---|---|---|---|
| ρ(S), L=2 | 0.3333333333 | 0.3333333333 | 0.3333333333 |
| ρ(S), L=3 | 0.3333333333 | 0.3333333333 | 0.3333333333 |

**ρ(S) = 1/3 EXACT at every λ ⟹ the criticality is WEIGHT-FREE — structural, living in S's
routing/cascade, NOT in the halving weight λ=½.** The mean-field shift is critically weighted for its own
structural reasons at all λ. Consequence for D-5: the (q≥5, λ≠½) **subcriticality cannot come from S** — it
must be carried by **F** (the fluctuation sectors) or the q-arithmetic. This localizes the boundary/gap
mechanism precisely, before D-1's fork or D-2's symbol.

## Status
D-0 frozen object **gated exact**: kernel identity, D0.1, D0.2, all four specializations, and D0.3 (with the
weight-free criticality finding). No definitional drift. **Cleared for D-1** (the H1/H2 fork).

Probe `probes/probe_trackD0_gate.py`.
