"""
path_B_explicit_phase.py — derive ψ(a) analytically from Cochrane Proposition 4 + truncated
3-adic log. Verify against empirical F̂(3a) values.

By Cochrane Prop 4: χ(1+3s) = e_q(-C_χ · L(1+3s)) where L(1+3s) = Σ_{j=1}^J (-1)^{j-1}/j · (3s)^j.

For our character χ_a (defined by χ_a(4) = e_period(-a)):
  χ_a(4) = e_q(-3a)  (since q/period = 3)
  And by Prop 4 with v = 4, s = 1: χ_a(4) = e_q(-C_a · L(4))
  So C_a · L(4) ≡ 3a mod q
  ⟹ C_a ≡ 3a / L(4) mod q

L(4) is the truncated p-adic log of 4, computable explicitly.

Then F̂(3a) = 3 · Σ_{v ∈ 1+3Z/q} χ_a(v) · e_q(v)
           = 3 · e_q(1) · Σ_{s=0}^{period-1} e_q(3s - C_a · L(1+3s))

Define P_a(s) := 3s - C_a · L(1+3s) — explicit polynomial in s of degree J.

Then F̂(3a) = 3 · e_q(1) · G(a) where G(a) = Σ_s e_q(P_a(s)) is a finite Gauss sum.

|G(a)| = √q (standard Gauss-sum theorem on cyclic principal-unit subgroup).
"""
import sys
import os
import math
import cmath
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def truncated_3adic_log(s, J):
    """L(1+3s) = Σ_{j=1}^J (-1)^{j-1}/j · (3s)^j as a Fraction."""
    L = Fraction(0)
    for j in range(1, J + 1):
        L += Fraction((-1)**(j-1), j) * (Fraction(3 * s)) ** j
    return L


def J_for_p3(m):
    """Compute J_{3,1,m} = max j with j - v_3(j) < m."""
    j = 1
    while True:
        # v_3(j+1)
        x = j + 1
        v = 0
        while x % 3 == 0:
            x //= 3
            v += 1
        if (j + 1) - v >= m:
            return j
        j += 1


def F_hat_via_explicit(r, a, c=1):
    """F̂(3a) via the explicit Cochrane formula:
    F̂(3a) = 3 · e_q(c) · Σ_{s=0}^{period-1} e_q(3cs - C_a · L(1+3s))

    Note: the c·v term where v = 1+3s gives c·(1+3s) = c + 3cs, so we get e_q(c) · e_q(3cs).
    """
    q = 3**(r+1)
    period = 3**r
    m = r + 1
    J = J_for_p3(m)

    # Compute L(4) = L(1+3·1) for s = 1
    L4 = truncated_3adic_log(1, J)
    # L4 should be a Fraction; reduce mod q
    L4_mod = (L4.numerator * pow(L4.denominator, -1, q)) % q

    # C_a satisfies C_a · L(4) ≡ 3a mod q
    # L(4) has v_3 = 1 (since L(1+3) starts with 3·1), so factor out:
    # L(4) = 3 · L̃ with L̃ unit mod p^{m-1}
    # 3·L̃·C_a ≡ 3a mod q ⟹ L̃·C_a ≡ a mod q/3 = p^{m-1}
    # C_a ≡ a · L̃^{-1} mod p^{m-1}
    p_mm1 = 3**(m-1)  # p^{m-1}
    L_tilde = L4_mod // 3  # L̃ = L(4)/3, integer since L(4) divisible by 3
    L_tilde_inv = pow(L_tilde, -1, p_mm1)
    C_a = (a * L_tilde_inv) % p_mm1

    # Compute Σ_s e_q(3cs - C_a · L(1+3s))
    total = complex(0, 0)
    for s in range(period):
        L_val = truncated_3adic_log(s, J)
        L_mod = (L_val.numerator * pow(L_val.denominator, -1, q) if s != 0 else 0) % q if s != 0 else 0
        # Phase = 3cs - C_a · L_mod
        phase_int = (3 * c * s - C_a * L_mod) % q
        total += cmath.exp(2j * cmath.pi * phase_int / q)

    # Multiply by 3 · e_q(c)
    F_hat = 3 * cmath.exp(2j * cmath.pi * c / q) * total
    return F_hat, C_a, L4_mod


def F_hat_direct(r, a, c=1):
    """Direct computation: F̂(3a) = Σ_{u=0}^{q-1} e_q(c·4^u - 3au)."""
    q = 3**(r+1)
    total = complex(0, 0)
    pw = 1
    for u in range(q):
        phase_int = (c * pw - 3 * a * u) % q
        total += cmath.exp(2j * cmath.pi * phase_int / q)
        pw = (pw * 4) % q
    return total


def main():
    print("# Path B: explicit phase formula via Cochrane Prop 4")
    print()

    for r in [2, 3]:
        q = 3**(r+1)
        period = 3**r
        m = r + 1
        J = J_for_p3(m)
        print(f"## r = {r}: q = {q}, period = {period}, J_{{3,1,{m}}} = {J}")
        print()

        L4 = truncated_3adic_log(1, J)
        print(f"  L(4) = L(1+3) (truncated to J = {J} terms) = {L4} = {float(L4):.6f}")
        L4_mod = (L4.numerator * pow(L4.denominator, -1, q)) % q
        print(f"  L(4) mod {q} = {L4_mod}")
        print()

        # For each a in support, compute via formula vs direct
        print(f"  {'a':>4}  {'C_a':>6}  {'F̂_formula':>30}  {'F̂_direct':>30}  {'diff':>10}")
        for a in range(period):
            if a % 3 != 1:
                continue
            F_form, C_a, _ = F_hat_via_explicit(r, a)
            F_dir = F_hat_direct(r, a)
            diff = abs(F_form - F_dir)
            match = "✓" if diff < 1e-9 else "✗"
            print(f"  {a:>4}  {C_a:>6}  {F_form.real:>+12.4f}{F_form.imag:>+12.4f}i   "
                  f"{F_dir.real:>+12.4f}{F_dir.imag:>+12.4f}i  {diff:>10.2e}  {match}")
        print()

    # Now the key insight: the phase ψ(a) depends on C_a and the Gauss-sum structure
    # of Σ_s e_q(3s - C_a · L(1+3s)).
    # For specific structure, compute the saddle-point / explicit formula.
    print("# Closed-form structure of ψ(a)")
    print()
    print("F̂(3a) = 3 · e_q(1) · G(a)  where G(a) = Σ_s e_q(P_a(s))")
    print("with P_a(s) = 3s - C_a · L(1+3s) and C_a = 3a · L(4)^{-1} mod q.")
    print()
    print("L(1+3s) = 3s - 9s²/2 + 9s³ - ... (truncated to J terms)")
    print()
    print("So P_a(s) = 3s - C_a · (3s - 9s²/2 + 9s³ - ...) = 3s(1 - C_a) + (C_a · 9/2) s² - C_a · 9 s³ + ...")
    print()

    # For r=2, J=3: P_a(s) = 3(1-C_a)s + (9C_a/2)s² - 9·C_a·s³ + 0·s^4+...
    # mod 27: depends on C_a mod 27.
    print("# At r=2 (q=27, J=3): P_a(s) = 3(1-C_a)s + (9C_a/2)s² - 9C_a·s³ mod 27")
    print()
    for a in [1, 4, 7]:
        F_form, C_a, _ = F_hat_via_explicit(2, a)
        print(f"  a = {a}: C_a = {C_a}")
        # Compute P_a coefficients mod 27
        # 3(1-C_a) mod 27, 9C_a/2 mod 27, -9C_a mod 27
        coef_1 = (3 * (1 - C_a)) % 27
        # 9·C_a / 2 mod 27: need 2^{-1} mod 27 = 14
        coef_2 = (9 * C_a * 14) % 27
        coef_3 = (-9 * C_a) % 27
        print(f"    P_a coefficients mod 27: u^1 → {coef_1}, u^2 → {coef_2}, u^3 → {coef_3}")
    print()


if __name__ == "__main__":
    main()
