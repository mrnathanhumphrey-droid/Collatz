"""
result_77_T_diagonal.py — derive the 2x2 diagonal-approximation operator T_diag
acting on (P_+, P_-) from Tao's recursion + R66 chain rule + class symmetry.

DERIVATION (rigorous):

Tao's recursion: μ̂_{n+1}(ξ) = Σ_v 2^{-v} A_v(ξ) μ̂_n(ξ·2^{-v} mod 3^n)

R66 chain rule: v even → r' ≡ 1 mod 3 (class +); v odd → r' ≡ 2 mod 3 (class -)

So:
  μ̂_{n+1}^+(ξ) = Σ_{v even} 2^{-v} A_v(ξ) [μ̂_n^+(ξ·2^{-v}) + μ̂_n^-(ξ·2^{-v})]   = E^+ + E^-
  μ̂_{n+1}^-(ξ) = Σ_{v odd}  2^{-v} A_v(ξ) [μ̂_n^+(ξ·2^{-v}) + μ̂_n^-(ξ·2^{-v})]   = O^+ + O^-

For the diagonal v = v' contribution to P_{n+1}^{ab}(c):
- ξ → u = ξ · 2^{-v} mod 3^{n+1}, a unit shuffle with class behavior:
    v even: u ≡ ξ mod 3 (class preserved), so sum stays in class c
    v odd:  u ≡ 2ξ mod 3 (class flipped), so sum moves to class 2c mod 3
- Project u to u mod 3^n: 3-to-1 cover.
- Each level-n moment scales by 3 (cover) × Σ_v 4^{-v} (where v has the given parity).

Coefficients:
  Σ_{v even} 4^{-v} = 1/16 + 1/256 + ... = (1/16)/(1 - 1/16) = 1/15
  Σ_{v odd}  4^{-v} = 1/4 + 1/64 + ...   = (1/4)/(1 - 1/16) = 4/15
  Total: 1/3 (matches E[2^{-2v}] for Geom(2))

Diagonal contribution to P_{n+1}^{ab}(c) = (3 · Σ_{v ∈ parity_a∩parity_b}) · P_n^{ab}(class flow)

Using class-c symmetry (P_n^{ab}(1) = P_n^{ab}(2) = P_n^{ab} for n ≥ 2):
  P_{n+1}^{++}(c)_diag = 3 · (1/15) · [P_n^{++} + P_n^{--}] = (1/5) · (P_n^+ + P_n^-)
  P_{n+1}^{--}(c)_diag = 3 · (4/15) · [P_n^{++} + P_n^{--}] = (4/5) · (P_n^+ + P_n^-)

(Cross terms P_n^{+-} = 0 for n ≥ 2 don't contribute.)

So in the diagonal approximation:
  (P_+, P_-)_{n+1} ≈ T_diag · (P_+, P_-)_n
  T_diag = (1/5) [[1, 1], [4, 4]]

Eigenvalues: λ_1 = 1 (eigenvector (1, 4)), λ_2 = 0 (eigenvector (1, -1))

CRITICAL OBSERVATION: T_diag has λ_1 = 1, meaning S preservation (S_{n+1} = S_n in diag approx).
This is consistent with conservation laws: in the diagonal approximation, the "Plancherel mass"
S is exactly preserved.

The actual rate-1/2 comes from the OFF-DIAGONAL (v ≠ v') corrections, which inject perturbation
into the (1, 4) deviation direction. The off-diagonal analysis is the remaining work.

This script verifies T_diag's spectrum and traces the off-diagonal correction empirically.
"""
import sys
import os
import time
from fractions import Fraction
import cmath
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def main():
    print("# Result 77: T_diag derivation + spectrum + off-diagonal trace")
    print()
    print("## T_diag (rigorous from Tao + R66):")
    print("  T_diag = (1/5) · [[1, 1], [4, 4]]")
    print()

    # T_diag matrix
    T_diag = (Fraction(1, 5)) * np.array([[Fraction(1), Fraction(1)],
                                           [Fraction(4), Fraction(4)]], dtype=object)

    # Eigenvalues over Q via characteristic polynomial:
    # det(T - λI) = (1/5 - λ)(4/5 - λ) - (1/5)(4/5) = λ² - λ
    # Roots: λ = 0 or λ = 1
    print("  Characteristic polynomial: λ² - λ = 0")
    print("  Eigenvalues: λ_1 = 1, λ_2 = 0")
    print()
    print("  Eigenvector at λ=1: (1, 4)  (preserves S)")
    print("  Eigenvector at λ=0: (1, -1) (zero direction)")
    print()

    # Verify on level-2 data
    print("## Verification at finite levels:")
    print()

    # Level 2: P_+ = 1/21, P_- = 4/21
    P2 = (Fraction(1, 21), Fraction(4, 21))
    # Apply T_diag
    P3_diag_pred = (Fraction(1, 5) * (P2[0] + P2[1]),
                    Fraction(4, 5) * (P2[0] + P2[1]))
    print(f"  Level 2: P_+ = {P2[0]} = {float(P2[0]):.10f}, P_- = {P2[1]} = {float(P2[1]):.10f}")
    print(f"  T_diag · P_2 = ({P3_diag_pred[0]} = {float(P3_diag_pred[0]):.10f},")
    print(f"                  {P3_diag_pred[1]} = {float(P3_diag_pred[1]):.10f})")
    # P3_diag_pred should equal (P_+^2 + P_-^2)/5, 4(P_+^2 + P_-^2)/5 = (1/21)/(...
    # Actually (1/21 + 4/21)/5 = (5/21)/5 = 1/21
    # So T_diag·P_2 = (1/21, 4/21) = P_2 (eigenvalue 1, since (1/21, 4/21) is on (1, 4) direction)

    # Compare to actual level 3: P_+ ≈ 0.04616, P_- ≈ 0.18463
    print(f"  Actual level 3:  P_+ = 0.046157, P_- = 0.184630")
    print(f"  Diff: ΔP_+ = {0.046157 - float(P3_diag_pred[0]):+.6e}, ΔP_- = {0.184630 - float(P3_diag_pred[1]):+.6e}")
    print()

    # The diff is the OFF-DIAGONAL contribution to P_3
    delta_P_plus = 0.046157 - float(P3_diag_pred[0])
    delta_P_minus = 0.184630 - float(P3_diag_pred[1])
    print(f"  Off-diagonal correction at n=2→3:")
    print(f"    ΔP_+ = {delta_P_plus:+.6e}, ΔP_- = {delta_P_minus:+.6e}")
    print(f"    Ratio ΔP_-/ΔP_+ = {delta_P_minus/delta_P_plus:.6f}  (should be 4 if on (1,4))")
    print()

    # Now: the OFF-DIAGONAL correction is what gives rate 1/2.
    # Track its decay across levels.
    print("## Off-diagonal correction sequence O_n := P_+^{n+1} - (P_+^n + P_-^n)/5:")
    print()
    # From data:
    P_plus = {2: 1/21, 3: 0.04615746803, 4: 0.04642144084, 5: 0.04655149198, 6: 0.04661687610}
    P_minus = {2: 4/21, 3: 0.18462987212, 4: 0.18568576337, 5: 0.18620596793, 6: 0.18646750441}
    # Verify class symmetry assumption: P_-^n = 4 P_+^n? not quite
    # Actually let's verify: 4/21 = 4 · (1/21)? Yes
    # P_-^3 = 4 · P_+^3? 4 · 0.04616 = 0.18464 ≈ 0.18463. YES (within rounding)
    # P_-^n / P_+^n stays close to 4 due to (1, 4) eigenvector

    print(f"  {'n':>3}  {'P_+^n':>14}  {'P_-^n':>14}  {'P_-/P_+':>10}  {'O_{n+1}':>14}  {'O ratio':>10}")
    O_seq = {}
    prev_O = None
    for n in [2, 3, 4, 5]:
        np1 = n + 1
        if np1 in P_plus:
            O = P_plus[np1] - (P_plus[n] + P_minus[n])/5
            ratio_pm = P_minus[n] / P_plus[n]
            O_ratio = O / prev_O if prev_O is not None else float('nan')
            print(f"  {n:>3}  {P_plus[n]:>14.10f}  {P_minus[n]:>14.10f}  {ratio_pm:>10.6f}  "
                  f"{O:>+14.6e}  {O_ratio:>+10.6f}")
            O_seq[n] = O
            prev_O = O
    print()

    # The O_n sequence is the off-diagonal correction. Its rate of decay = rate of S → 7/15.
    # Empirically O ratio → 1/2, confirming rate-1/2.
    print("## Asymptotic structure:")
    print("  P_+^n = 7/150 + δ_n,  P_-^n = 14/75 + 4δ_n (always on (1,4) line)")
    print("  S_n = 7/15 + 10·δ_n  (so ε_n = 10·δ_n)")
    print()
    # Empirical δ_n
    print(f"  δ_n = P_+^n - 7/150:")
    for n in [2, 3, 4, 5, 6]:
        delta = P_plus[n] - 7/150
        print(f"    δ_{n} = {delta:+.6e}")
    print()

    # The rate of δ_n → 0 is the rate-1/2 we seek.
    # In T_diag approx, δ would be eigenvalue-1 (preserved). Off-diag changes this to 1/2.
    print("## Why T_diag eigenvalue 1 lifts to actual eigenvalue 1/2:")
    print("  T_diag's null mode (1, -1) has eigenvalue 0 (instantly killed).")
    print("  T_diag's slow mode (1, 4) has eigenvalue 1 (would preserve δ).")
    print("  But OFF-DIAG injection adds bilinear terms μ̂_n^+(ξ_1) μ̂_n^-*(ξ_2) for ξ_1 ≠ ξ_2.")
    print("  These contribute to (P_+, P_-) via cross-frequency interference.")
    print("  Effective lift: λ_2 = 1 - (off-diag contraction) = 1/2.")
    print()
    print("  The factor 1/2 emerges from P(v=1)·1 = 1/2 (Geom(2) probability of v=1).")
    print("  Specifically: at level n+1, the v=1 contribution gives a 'fresh' deviation")
    print("  scaled by 2^{-1} = 1/2 from the previous level. This is the dominant")
    print("  off-diagonal coupling, giving the asymptotic decay rate 1/2.")
    print()

    # Final asymptotic
    # ε_n ~ -1/30 · (1/2)^n (from k=6 fit)
    # 1/30 = 7/15 / 14 = (S_∞)/14
    print("## Sharpened closed form (numerical fit through k=6):")
    print("  ε_n = S_n - 7/15 ~ -(1/30) · (1/2)^n + O((1/4)^n)")
    print("  Equivalently: c = 7/45 + (-1/90) · (1/2)^k · 3^k = 7/45 - (1/90)·(3/2)^k... hmm")
    print()
    print("  More precisely: ‖d_{k+1}‖² · 3^k = c + ε_{k+1}/3 = 7/45 - (1/90) · (1/2)^{k+1} + ...")
    print("    = 7/45 - (1/180) · (1/2)^k + lower order")
    print()
    print("  The 1/30 = S_∞/14 = 7/(15·14) = 7/210 — the 14 likely comes from")
    print("  some combinatorial factor in the off-diagonal expansion (open).")
    print()


if __name__ == "__main__":
    main()
