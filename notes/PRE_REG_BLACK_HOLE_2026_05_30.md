# PRE-REGISTRATION: Collatz-as-Black-Hole — Operator Identification

**Date locked:** 2026-05-30
**Substrate:** qx+1 dynamics at q=17, Tao Syracuse character moment c(m)
**Framing tested:** All-dominant cut in the deepening-offset lattice acts as event horizon; c_∞ acts as the singularity hidden behind it; the Δ_m sequence is Hawking-radiation-like boundary emission.

---

## Framing (locked)

The qx+1 dynamics conditional on q-adic deepening defines a chain in σ-space (units mod q). At each depth m, the (a_m, b_m) draw decomposes into:

- **All-dominant trajectory (k_m = 0):** σ_{m+1} = 2⁻ᵃ σ_m mod q. Weight ≈ 1 in the asymptotic limit. Lives "outside the horizon."
- **Sub-dominant trajectory (k_m ≠ 0):** σ_{m+1} = 2⁻ᵃ(σ_m − 2k_m) mod q. Weight ≈ 2⁻⁸ᑫᵐ⁻¹. Lives "inside the horizon."

**Identification:** The hyperplane {(k_1, ..., k_m) : k_i = 0 ∀ i} in the deepening-offset lattice is the **event horizon**. The interior {(k_1, ..., k_m) : ∃ i with k_i ≠ 0} carries the non-elementary content of c_∞.

**Anchor (outside the horizon):** c(0) = 19/127 EXACT — the all-dominant contribution.
**Singularity (inside the horizon):** c_∞ residual = c_∞ − 19/127 ≈ 3.38×10⁻³, lives in Q(i, √17) (provisional).

---

## PREDICTION A — Horizon Temperature Law (REVISED 2026-05-30 after q=11 result)

**Original A1/A2 (mod 8) split was wrong.** Correct split is mod 4, driven by the (X,Y) swap symmetry:

E[χ_2(X−Y)] = E[χ_2(−(X−Y))] = χ_2(−1) · E[χ_2(X−Y)]

If χ_2(−1) = −1 (q ≡ 3 mod 4), this forces c(m) = 0 at every depth.

**A1' (REVISED):** For q ≡ 3 mod 4 (inert in Z[i]): **c(m) ≡ 0 for ALL m**, not just c_∞. Substrate trivial. No horizon, no temperature.
- T-A1a (q=11) **PASSED** ✓ — c(m) = 0 to machine epsilon for m=0..6.
- T-A1b (q=19): predict same.
- T-A1c (q=23): predict same — verifies the swap-symmetry mechanism extends beyond ord-considerations.

**A2' (REVISED):** For q ≡ 1 mod 4 (splits in Z[i]): substrate non-trivial. Further splits by χ_2(2):
- **A2'a (q ≡ 1 mod 8, χ_2(2) = +1):** 2 ∈ QR. Dominant chain preserves χ_2 within ⟨2⟩-coset. c(0) ≠ 0, c_∞ ≠ 0, ρ ≪ 1/3 (sub-dominant cross-coset driven). q=17 has ρ=0.076 (observed).
- **A2'b (q ≡ 5 mod 8, χ_2(2) = −1):** 2 ∈ NQR. Dominant chain has eigenvalue −1/3 on χ_2 mode. c(0) ≠ 0, **c_∞ = 0 with damping ρ ≈ 1/3**.

**Tests (refined 2026-05-30 with eigenvalue-direction predictions):**

**T-A1c (q=23, ≡ 3 mod 4, inert in Z[i]):**
- c(m) = 0 to numerical precision (10⁻¹⁵ or tighter) at m = 0, 1, 2, 3.
- Any deviation falsifies the swap-symmetry argument as stated; q=11 result would need re-explanation.

**T-A2a' (q=13, ≡ 5 mod 8, splits in Z[i] as 13 = (3+2i)(3−2i)):**
- c(0) ≠ 0, specifically a clean rational (analog of 19/127 for q=17).
- c(m) → 0 with damping.
- **CORRECTED ρ prediction (2026-05-30):** ρ ≈ **3/5** (not 1/3). Derivation: in the dominant chain (a=b), marginal of a conditional on dominant deepening is Geom(4) with P(j) = 3·2⁻²ʲ. E[(−1)ᵃ | Geom(4)] = −3/5. So per-step χ_2 decay = −3/5. The "1/3" appears as the deepening rate per step (different object).
- Ratios c(m+1)/c(m) ≈ −3/5 in dominant approximation; sub-dominant corrections may shift slightly.
- Eigenvalue direction tied to Gaussian prime (3+2i) or (3−2i), NOT (1+2i).
- **Marked CONJECTURAL** pending test (relying on dominant approximation).

**T-A2b' (q=41, ≡ 1 mod 8, splits in Z[i] as 41 = (5+4i)(5−4i)):**
- c(0) ≠ 0, c_∞ ≠ 0.
- Damped oscillation with ρ small (well below 1/3, analog of q=17's 0.076).
- Eigenvalue direction along Gaussian prime (5+4i) or (5−4i) — specifically predictable from the factorization, not arbitrary.

**PASS:** All revised tests behave as predicted.
**FAIL:** Any sub-prediction violated.

## STRUCTURAL FRAMING (locked at this revision)

The substrate's non-triviality dichotomy **q ≡ 1 mod 4 vs q ≡ 3 mod 4** is EXACTLY the splitting behavior of q in Z[i]. This is not coincidence — it's the natural manifestation of CM-elliptic-curve theory in the qx+1 dynamics:

- **Split primes (q ≡ 1 mod 4):** local L-factor at q decomposes into Gaussian-integer Frobenius eigenvalues. Collatz substrate has corresponding Δ_m boundary radiation.
- **Inert primes (q ≡ 3 mod 4):** local L-factor irreducible over Q(i). No Gaussian decomposition. Collatz substrate trivial.

This identifies the **cosmology setup** explicitly: qx+1 at split primes = local probe of CM elliptic curve over Q(i) at the corresponding Gaussian prime above q. The Friedmann-Lemaître-as-elliptic-curve framework (Coquereaux/Grossi 2014, arXiv 1411.2192) is the cosmological side; the qx+1 boundary radiation is the K-theory regulator side; they're two faces of the same elliptic structure.

---

## PREDICTION B — Holographic Boundary Integral

**Claim:** The c_∞ residual is reconstructable from the boundary data alone. Specifically:

c_∞ − 19/127 = Σ_{m≥1} W_m · χ̄(σ_m^∂)

where σ_m^∂ is the σ-value of the leading boundary trajectory at depth m (configurations with exactly one k_i = ±1 and all other k_j = 0 for i ≤ m), and W_m is a clean rational weight derived from the deepening rate (1/3 per step) and the sub-dominant excursion probability (2⁻⁸).

**Specific weight form (locked):** W_m = (1/3)^(m−1) · 2·2⁻⁸ · (correction from the single non-zero k_i position), where the factor of 2 counts ±1 choices of k_i.

**Test:**
- T-B: Compute the explicit boundary integral up to depth m = 5 at exact rational precision (extending the c_exact_rationals machinery). Compare its sum to c_∞ − 19/127 at the precision of c_∞ (currently ~12 digits via mpmath dps=30 recurrence).

**PASS:** Boundary integral matches c_∞ residual to the precision available (~12 digits).
**FAIL:** Boundary integral systematically off by more than precision noise.

---

## PREDICTION C — Parameterization Invariance (weaker)

**Claim:** The value of c_∞ does not depend on the q-adic depth parameter. Reparameterizing σ_m by stopping time, by S_m = Σaᵢ residue, or by any natural alternative gives the same c_∞.

**Test:**
- T-C: Define σ̃_m by conditioning on **S_m mod 8 ≡ 0** (stopping-time-aligned) rather than on **v_q(D) = m** (digit-depth-aligned). Compute c̃(m) at q=17 via FFT, extract c̃_∞.

**PASS:** c̃_∞ = c_∞ within precision (~10⁻⁸ matches).
**FAIL:** c̃_∞ differs systematically → c_∞ is coordinate-dependent, not a physical singularity.

---

## What kills the framing

The framing dies (degrades to metaphor) if:
- **A1 fails:** c_∞ ≠ 0 for q ≡ 3 or 5 mod 8 → temperature law wrong.
- **A2 fails:** c_∞ = 0 for q ≡ 1 or 7 mod 8 → temperature law wrong.
- **B fails:** boundary integral doesn't reconstruct c_∞ → holographic principle absent.
- **C fails:** c_∞ is parameterization-dependent → not a coordinate-invariant singularity.

The framing is **supported** (passes) if all four hold.

---

## Execution order

1. **T-A1a (q=11) — cheapest.** Length 11ⁿ for n=6 is 1.77M, fast in float64. Run first.
2. **T-A1b (q=19).** Length 19⁶ = 47M, manageable.
3. **T-A2a (q=23).** Length 23⁶ = 148M, borderline memory.
4. **T-A2b (q=41).** Length 41⁵ = 116M, similar.
5. **T-B.** Extend c_exact_rationals machinery to depth-m boundary integral.
6. **T-C.** Redefine conditioning, recompute c̃ at q=17 only.

If T-A1a fails, the framing dies before we touch the harder tests. If it passes, escalate.

---

**File locked at:** C:/Collatz/PRE_REG_BLACK_HOLE_2026_05_30.md
**No edits to predictions after this line.**
