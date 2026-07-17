"""
path_B_saddle_point.py — derive saddle-point closed form for the Gauss sum G(a) and
extract explicit ψ(a) closed-form.

For r=2: P_a(s) = 3(1-C_a)s + 18·s² + 18·s³ mod 27 (where C_a ≡ 1 mod 3)
        = 3(1-C_a)s + 18·s²·(1+s) mod 27

Note: at higher s³ coefficient mod 27, we need to track 3-adic structure carefully.

The Gauss sum G(a) = Σ_{s=0}^{period-1} e_q(P_a(s)).
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
    L = Fraction(0)
    for j in range(1, J + 1):
        L += Fraction((-1)**(j-1), j) * (Fraction(3 * s)) ** j
    return L


def J_for_p3(m):
    j = 1
    while True:
        x = j + 1
        v = 0
        while x % 3 == 0:
            x //= 3
            v += 1
        if (j + 1) - v >= m:
            return j
        j += 1


def G_a(r, C_a):
    """Compute G(a) = Σ_{s=0}^{period-1} e_q(P_a(s)) where P_a(s) = 3s - C_a · L(1+3s).
    period = 3^r, q = 3^{r+1}.
    """
    q = 3**(r+1)
    period = 3**r
    J = J_for_p3(r + 1)
    total = complex(0, 0)
    for s in range(period):
        L_val = truncated_3adic_log(s, J)
        L_mod = ((L_val.numerator * pow(L_val.denominator, -1, q)) % q) if s != 0 else 0
        phase_int = (3 * s - C_a * L_mod) % q
        total += cmath.exp(2j * cmath.pi * phase_int / q)
    return total


def main():
    print("# Path B saddle-point analysis: extract ψ(a) closed form")
    print()

    for r in [2, 3, 4]:
        q = 3**(r+1)
        period = 3**r
        J = J_for_p3(r + 1)
        p_mm1 = 3**r
        L4_frac = truncated_3adic_log(1, J)
        L4_mod = (L4_frac.numerator * pow(L4_frac.denominator, -1, q)) % q
        L_tilde = L4_mod // 3
        L_tilde_inv = pow(L_tilde, -1, p_mm1)

        print(f"## r = {r}: q = {q}, period = {period}, J = {J}")
        print(f"  L̃ = L(4)/3 mod {p_mm1} = {L_tilde}")
        print(f"  L̃^(-1) mod {p_mm1} = {L_tilde_inv}")
        print()

        # For a ∈ support, compute C_a, then G(a), then ψ(a) = G(a)/√q
        sqrt_q = math.sqrt(q)
        print(f"  {'a':>4}  {'C_a':>6}  {'s*':>6}  {'G(a)':>30}  {'|G|/√q':>10}  {'ψ(a) phase / 2π':>20}")
        psi_data = []
        for a in range(period):
            if a % 3 != 1:
                continue
            C_a = (a * L_tilde_inv) % p_mm1
            G = G_a(r, C_a)
            mag_norm = abs(G) / sqrt_q
            phase_norm = (cmath.phase(G) / (2 * cmath.pi)) % 1.0

            # Saddle point: dP/ds = 0 mod q (modular saddle)
            # P_a(s) = 3s - C_a · L(1+3s) → dP/ds = 3 - 3C_a + 9C_a·s + ...
            # For r=2: dP/ds = 3 - 3C_a + 9C_a·s - 27C_a·s² ≡ 3(1 - C_a) + 9C_a·s mod 27
            # Setting = 0: 9C_a·s ≡ 3(C_a - 1) mod 27, s ≡ (C_a - 1)/(3C_a) mod 3
            # Since C_a ≡ 1 mod 3: s* = (C_a - 1)/3 mod 3
            s_star = ((C_a - 1) // 3) % 3

            psi_data.append((a, C_a, s_star, G))
            print(f"  {a:>4}  {C_a:>6}  {s_star:>6}  {G.real:>+12.4f}{G.imag:>+12.4f}i   "
                  f"{mag_norm:>10.4f}  {phase_norm:>20.6f}")
        print()

        # Now: hypothesize ψ(a) phase has saddle-point form e_q(P_a(s*))
        # where s* is the saddle point. Compute predicted phase:
        print(f"  Saddle-point phase prediction P_a(s*) mod q vs empirical phase × q:")
        print(f"  {'a':>4}  {'C_a':>6}  {'s*':>6}  {'P_a(s*) mod q':>16}  {'empirical phase × q':>22}  {'match?':>8}")
        for (a, C_a, s_star, G) in psi_data:
            # Compute P_a(s*) mod q
            L_val = truncated_3adic_log(s_star, J)
            if s_star == 0:
                L_mod = 0
            else:
                L_mod = (L_val.numerator * pow(L_val.denominator, -1, q)) % q
            P_at_saddle = (3 * s_star - C_a * L_mod) % q

            empirical_phase = (cmath.phase(G) / (2 * cmath.pi)) % 1.0
            empirical_phase_q = (empirical_phase * q) % q
            # Saddle phase prediction (offset by some constant maybe)
            match = "✓" if abs((P_at_saddle - empirical_phase_q) % q) < 0.5 or \
                          abs((P_at_saddle - empirical_phase_q) % q - q) < 0.5 else "✗"
            print(f"  {a:>4}  {C_a:>6}  {s_star:>6}  {P_at_saddle:>16}  {empirical_phase_q:>22.4f}  {match:>8}")
        print()
        print()


if __name__ == "__main__":
    main()
