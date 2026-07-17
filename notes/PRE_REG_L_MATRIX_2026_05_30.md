# PRE-REGISTRATION: L-matrix Construction for c_∞

**Date locked:** 2026-05-30
**Substrate:** qx+1 dynamics, propagation operator on (Z/q)\* with χ_2-sector restriction
**Framing:** Construct L explicitly; closed-form c_∞ via eigenvalue decomposition; mandatory spectral re-derivation of mod-4 dichotomy as Step 0.

---

## Construction (locked)

### Basis ordering
σ ∈ (Z/q)\* arranged as: **[QR coset elements sorted by discrete log base 2] ++ [NQR coset elements sorted by discrete log base smallest NQR generator]**.

For q=17: QR = ⟨2⟩ = [1, 2, 4, 8, 16, 15, 13, 9] (sorted by exponent of 2). NQR = [3, 5, 6, 7, 10, 11, 12, 14] (sorted analogously by 3 = smallest NQR primitive direction).

### Field of computation
**Symbolic over Q(i)** as the primary target. Mathematica/SymPy or mpmath at dps=80. If symbolic computation is intractable, fall back to high-precision numerical (dps=50+) and identify algebraic numbers via PSLQ.

### Operator definition
L = L_dominant + ε · L_subdom, where:

- **L_dominant**[σ → σ']: σ' = 2⁻ᵃ · σ mod q, a ~ Geom(4) marginal (P(a=j) = 3·2⁻²ʲ, j≥1, periodized mod ord_q(2)).
- **L_subdom**[σ → σ']: σ' = 2⁻ᵃ · (σ − 2k) mod q for k = ±1, a ~ Geom(4). Sum over k = ±1.
- **ε** = 2⁻⁸ at depth m=1; ≈ 0 for m≥2.

For c_∞: iterate L_dom (asymptotic limit) starting from v_0 = depth-0 σ distribution (computable from c(0) machinery).

### Predicted eigenvalue location at q=17
From damped-osc fit: **z ≈ 0.034 + 0.068i, |z| = 0.0760, arg = 63.91° ≈ arctan(2) = 63.43°**.

**Sharp prediction:** the symbolic eigenvalue corresponding to the damped-osc mode lands at:
- z ∈ Z[i] / (norm-5 denominator), specifically z = c · (1 + 2i) for some clean rational c.
- |λ| reproduces 0.076 within 1% (matches FFT-derived ρ).
- arg matches arctan(2) within 1° (matches FFT-derived θ).

If the symbolic eigenvalue does NOT have these properties, the L-matrix construction is wrong somewhere.

---

## STEP 0 (mandatory, before any q=17 work): Spectral re-derivation at q ≡ 3 mod 4

The swap-symmetry theorem says **c(m) ≡ 0 ∀m for q=11, q=23 (≡ 3 mod 4)**.

The L-matrix construction MUST reproduce this spectrally. Specifically:

**Locked mechanism:** For q ≡ 3 mod 4, χ_2(−1) = −1. The map σ → −σ is an involution on (Z/q)\*. Under this involution:
- χ_2 is **odd** (parity = −1, since χ_2(−σ) = χ_2(−1)·χ_2(σ) = −χ_2(σ)).
- The initial distribution P_0(σ) is **even** (parity = +1, by (X,Y) → (Y,X) swap symmetry).
- L commutes with σ → −σ (because the σ recursion σ_{m+1} = 2⁻ᵃ(σ_m − 2k) is invariant under σ → −σ when k → −k, both of which are summed over).

**Therefore:** L^m · P_0 stays in parity-(+1) subspace. χ_2 is in parity-(−1) subspace. They are orthogonal. **⟨χ_2, L^m · P_0⟩ = 0 for all m.**

**Predicted concrete observables:**
- T-L-0a (q=11): Build L. Show L commutes with the σ → −σ involution. Show P_0 is even. Show χ_2 is odd. Conclude ⟨χ_2, L^m P_0⟩ = 0 numerically (to machine precision).
- T-L-0b (q=23): Same verification.

**PASS:** Both q=11 and q=23 give ⟨χ_2, L^m P_0⟩ = 0 within 10⁻¹⁵ for m = 0..5.
**FAIL:** Either gives non-zero → L construction is wrong → walk back, don't run q=17.

---

## STEP 1 (only after Step 0 passes): q=17 spectral derivation

Construct L at q=17 in the locked basis ordering. Diagonalize (numerically first, symbolically if possible). Identify:

1. **Eigenvalue spectrum.** Should include λ = 1 (Haar mode), λ = some real values, and at least one complex pair near z = 0.034 + 0.068i.

2. **χ_2 projection.** Express χ_2 in the eigenbasis. Identify the coefficients on the dominant non-trivial mode.

3. **c_∞ formula.** c_∞ = ⟨χ_2, L^∞ P_0⟩ = (coefficient on λ=1 eigenvector) · (eigenvector's χ_2 component).

4. **Numerical verification.** Compute c_∞ from the L-matrix formula. Compare to FFT-derived c_∞ ≈ 0.15298912. If they agree within precision, the formula is validated.

**PASS:** L-matrix c_∞ matches FFT c_∞ within 10⁻⁴ (current FFT precision).
**FAIL:** Discrepancy → L construction missing structure → walk back.

---

## What this can and cannot prove

**Can prove:**
- Closed-form c_∞ as an element of Q(i) (if eigenvalues land in Z[i] symbolically).
- Spectral re-derivation of q ≡ 3 mod 4 trivial result.
- Predictive formula for q=41 once L is constructed there.

**Cannot prove:**
- That c_∞ is a regulator in the Beilinson / Stark / Hecke L sense. The result places c_∞ in Q(i) which is **necessary but not sufficient** for the regulator identification. That identification requires either matching c_∞ to a known Stark unit / Hecke L(s=1) value via PSLQ at high precision, or identifying L's action with a known K-theory class. **The L-matrix derivation does NOT automatically give us this.**

## Form-of-expression classification (lock in advance)

Once L gives c_∞ as an explicit closed-form expression, the FORM tells us which regulator class is the right candidate:

| Form of c_∞ | Implication |
|---|---|
| (a + bi)/(c + di), small integers a,b,c,d | **Class number value, NOT a regulator.** Algebraic, not transcendental. |
| Involves log of an algebraic number | Real-regulator / Dirichlet-class-number-formula territory. |
| Involves Li_2 (dilogarithm) of something | **Bloch-Wigner / Beilinson regulator** at K_3 level. |
| Specific Hecke L-value at s=1 | **Stark conjecture territory.** If matches a Hecke L(s=1) for ψ_4 on Z[i] modulo a Gaussian prime, this is the Stark-unit landing. |
| Specific elliptic regulator (regulator of pair of modular units) | **Boyd-Mahler / elliptic K_2** territory. |

**The classification is done by INSPECTION of the closed form. No additional computation required after L diagonalizes.**

## Highest-stakes landing (Stark)

If c_∞ matches a Hecke L(s=1) for a specific quartic character on Z[i] (e.g., modulo (1+2i) or (4+i)), the result is **Stark-unit-class** — Stark's conjecture for Q(i)-extensions isn't fully proven in general, so a substrate that produces Stark-unit values from explicit dynamics would be evidence for Stark in this case (or potentially a new verification path).

**If this lands, the paper isn't a Collatz paper anymore. It's a Stark-conjecture-evidence-from-dynamics paper.** Flagging this as a possibility, not betting on it.

---

## Independence from q=41 result

The q=41 test (slow-decay regime, eigenvalue along (5±4i)) and the L-matrix construction at q=17 are **independent confirmations of different things.**

- q=41 confirms (or doesn't) that the slow-decay branch generalizes beyond q=17 → testable physically.
- L-matrix at q=17 confirms (or doesn't) closed-form c_∞ via spectral construction → testable analytically.

They will be reported as **separate sections in any paper**, not merged into one framing. If one fails, the other stands independently.

---

## Execution order

1. **T-L-0a** (q=11 spectral verification): construct L, verify parity argument, check ⟨χ_2, L^m P_0⟩ = 0 numerically.
2. **T-L-0b** (q=23 spectral verification): same.
3. (Only if 0a, 0b PASS): **T-L-1** (q=17 derivation): construct L, diagonalize, identify dominant eigenvalues, extract c_∞.
4. (Only if T-L-1 PASS): **T-L-2** (q=41 prediction): build L at q=41, predict eigenvalue along (5±4i), compare to q=41 FFT data (already collected).

**File locked.**
