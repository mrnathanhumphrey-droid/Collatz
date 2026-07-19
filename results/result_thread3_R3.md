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

## R3-D — THE L-LAW (the last owed piece)
The secular product as an L-sequence, heading to 7/15 = 0.46667:
- **L=2: 0.3066** · **L=3: 0.3512** · **L=4: pending** (`probe_thread3_R3_L4.py`, appended on landing).

Climbing toward 7/15. Note the braid: ρ_L alternates above/below 1/3 (even/odd L), so the *plateau
interpretation* is only clean at odd L, but the *product* 3·g_L·φ_tow·ψ_kin is well-defined at every L and is the
sequence whose L→∞ limit the crown must pin to 7/15 in closed form.

## Status
**Mechanism DERIVED and gated** (renewal quarantined throughout — nothing here uses S∞ or the renewal limit):
R3-A pins the backward shell convention (2/3, 10/21); R3-C verifies the Δ-cancellation at L=3 to 2.6%
(g/Δ = 189 × Δ → finite g·φ·ψ), with the sign corrected to **+**. **Owed for theorem grade:** (i) done — convention
frozen; (ii) partially — verified clean at L=3, L=2 outside the regime, L=4 pending; (iii) the L→∞ closed-form
limit of 3·g_L·φ_tow·ψ_kin = 7/15 (the product 0.307 → 0.351 → … climbs toward it; the closed form is the pen's
final step). Caveat kept visible: the 2×2 product is ~80% of the finite-L plateau (subdominant modes fill the
rest), so "just the exceptional point" is the *leading* term, not the whole flat level at finite L.
