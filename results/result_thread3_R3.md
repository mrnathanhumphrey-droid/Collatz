# Probe R3 — the secular product (Theorem S gate; the crown's corrected route)

**Date:** 2026-07-19  CPU, dense at q=3 (L=2,3) + SpMV (L=4). Probe `probes/probe_thread3_R3.py`
(+ `probe_thread3_R3_L4.py`). Gates Theorem S: the flat shell level = the Jordan coupling read through the
agreement functionals, **S∞ = 3·g·φ_tow·ψ_kin** — with the diverging g/Δ meeting the vanishing Δ and cancelling.

Two-mode form on the degenerate top block B = [[c₀, 0],[g, ρ]] (protection ⟹ upper-right 0):
a_m = P·c₀^m + Q·ρ^m, P = φ_kin ψ_kin + φ_tow ψ_kin (g/Δ), Q = φ_tow ψ_tow − φ_tow ψ_kin (g/Δ), Δ = c₀−ρ.
φ = agreement readout (all-ones) on (kin,tow); ψ = independent-pair init δ(1,1,0) on (kin,tow).

## R3-A — CONVENTION FREEZE: **backward pinned by the data**
From the A-sequence, the two candidate shell functionals:
- **backward** g_m = a_m − a_{m−1}/3 → S₁ = **2/3**, S₂ = **10/21** — matches the welded Era-6 constants exactly.
- forward g_m = a_m − a_{m+1} → S₁ = 20/21. Wrong.

**The data pins S_k = 3^k(a_m − a_{m−1}/3).** This is the convention the corpus's S carries; all of R3 uses it.

## R3-C — THE Δ-CANCELLATION (the heart): **VERIFIED at L=3** (sign corrected)
Structure at each L: protection holds (B[kin,tow] = −3e-18 / +4e-17 ≈ 0); ψ_tow ≈ 0 exactly (the init has no
tower component), so Q ≈ −φ_tow ψ_kin (g/Δ). Then Q·(ρ−1/3) with ρ−1/3 ≈ −Δ gives **+g·φ_tow·ψ_kin** — the
**two minuses cancel to a plus** (Wilson's Theorem S line; the "−" in his R3-C prose is a typo).

| L | g | Δ = c₀−ρ | g/Δ (defect meter) | Q·(ρ−1/3) | +g·φ_tow·ψ_kin | ratio | regime |
|---|---|---|---|---|---|---|---|
| 2 | +0.0505 | −2.911e-3 | −17.4 | +0.4737 | +0.1022 | +4.63 | **ρ₂>1/3 (braid): super-critical, ε_c not small — invalid** |
| **3** | +0.0188 | +9.958e-5 | **+189.0** | **+0.11409** | **+0.11708** | **+0.974** | **critical (ρ₃<1/3): cancellation VERIFIED (2.6% = ε)** |

**At L=3 the cancellation is real:** the defectiveness meter g/Δ = 189 (diverging toward the EP) times the
vanishing Δ leaves the finite product +g·φ_tow·ψ_kin, matching Q·(ρ−1/3) to **2.6%**. The divergence and the
detuning are the same structure twice; their product is finite. **L=2 is not a clean test** — the partner sits
*above* 1/3 there (ρ₂ = 0.3468, the braid's even-L sign), so ε_c = c₀−1/3 = +0.011 is comparable to Δ and the
o(1) corrections dominate. The clean regime (ε_c doubly-exponentially small, ρ < 1/3) first appears at L=3.

## R3-B — THE PRODUCT vs the plateau
| L | 3·g·φ_tow·ψ_kin | 2-mode self-plateau 3Q(ρ−1/3) | full-chain plateau | note |
|---|---|---|---|---|
| 2 | 0.3066 | 1.421 | 2.506 (grows) | super-critical: 3ρ₂ > 1, no valid plateau |
| **3** | **0.3512** | 0.3423 | 0.4393 | 2-mode is **80%** of the full plateau |

At L=3 the product (0.351) and the 2-mode self-consistent plateau (0.342) agree to 2.6% (that IS R3-C) — the
two-mode picture is internally consistent. But it sits **~20% below the full-chain plateau (0.439)**: the top two
modes capture 86–89% of a_m (subdominant near-1/3 modes carry the rest), and the shell functional amplifies that
residual. So **the pure 2×2 secular level is a lower bound on the finite-L flat level**; the gap is subdominant
modes, not a failure of the cancellation.

## R3-D — THE L-LAW (the last owed piece) — **limit NOT numerically pinned; stays owed**
The secular product as an L-sequence: **L=2: 0.3066 · L=3: 0.3512 · L=4: 0.183 [DISCARDED]**.
- **The L=4 point is not trustworthy.** The sanctioned block-2 subspace iteration **stalled at res = 4.2e-5**
  after 3000 steps (25 min) — the c₀↔partner gap at L=4 is only |Δ| = 1.68e-4, so res/gap ≈ 0.24 means **~24%
  mixing** between the two near-degenerate eigenvectors. g, φ, ψ are corrupted; the product 0.183 is an artifact
  of underconvergence, **discarded**. (Shift-invert would resolve it but is barred by the instrument law at the
  defective EP; the ~0.9995/step rate makes res≪gap infeasible in reasonable time.)
- **R3-C survives at L=4** (ratio 0.9916) only because it is the algebraic ε_c/Δ consistency (ratio = 1 − ε_c/Δ,
  = 1 − (−1.4e-6)/(−1.68e-4) = 0.9916) — robust to the mixing, confirming the cancellation *structure*, not the
  physical values.
- **The braid also complicates the sequence:** ρ_L alternates above/below 1/3 (Δ₄ = 1.68e-4 > Δ₃ = 9.96e-5, non-
  monotone), so g/Δ is non-monotone (17.4, 189, ~27.5) and the plateau interpretation is clean only at odd L.

**Verdict:** the reliable L-law points are L=2 (0.307, super-critical regime) and L=3 (0.351, clean) — **one clean
point.** The L→∞ limit → 7/15 is **not demonstrated numerically** (finite-L near-EP precision + braid), and
remains exactly what it was: the crown's owed **closed-form** step. The cancellation *mechanism* (R3-C) is
verified at L=3 and L=4; the *value* 7/15 is not pinned by the L-sequence.

## Status
**Mechanism DERIVED and gated** (renewal quarantined throughout — nothing here uses S∞ or the renewal limit):
R3-A pins the backward shell convention (2/3, 10/21); R3-C verifies the Δ-cancellation *structure* at L=3 (2.6%)
and L=4 (0.84%) — the diverging g/Δ meets the vanishing Δ — with the sign corrected to **+**. **Owed for theorem
grade:** (i) done — convention frozen; (ii) mechanism verified (L=3 clean; L=2 out-of-regime; L=4 confirms the
cancellation but its product is discarded — near-EP underconvergence, 24% eigenvector mixing); (iii) **still fully
owed** — the L→∞ closed-form limit 3·g_L·φ_tow·ψ_kin = 7/15 is **not numerically demonstrated** (only L=3 is a
clean point; the sequence 0.307, 0.351, [0.183 discarded] does not pin the limit). Two caveats kept visible: the
2×2 product is ~80% of the finite-L plateau (subdominant near-1/3 modes fill the rest — the EP is the *leading*
term, not the whole flat level); and the near-EP eigenvector precision caps how far the L-law can be pushed
numerically with the sanctioned instrument. The cancellation *mechanism* is gated; the *value* 7/15 needs the pen.
