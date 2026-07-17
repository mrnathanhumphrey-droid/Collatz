"""
result_78_polynomial_identification.py — Step 1 of R78: identify polynomials f(u), g(u)
in Q(u) for Kalafatelis's S_{r,ℓ,ε}(m) sum, then apply Cochrane Theorem 2.

Kalafatelis's sum (from brief):
  S_{r,ℓ,ε}(m) = Σ_{u=0}^{3^{r-1}-1} e_{3^{r+1}}(c_{ℓ,ε} · 4^u - 9mu)

where c_{ℓ,ε} = 2^ε · ω_r^ℓ ∈ (Z/3^{r+1})*, ω_r is a primitive cube root of unity in (Z/3^{r+1})*.

Strategy:
1. Convert 4^u to polynomial in u via binomial expansion modulo 3^{r+1}:
   4^u = (1+3)^u = Σ_{k=0}^r C(u, k) · 3^k mod 3^{r+1}
   (terms with k ≥ r+1 vanish mod 3^{r+1})

2. Set f(u) = 1 (constant — character is trivial since we have only phase factors)
   Set g(u) = c_{ℓ,ε} · Σ_{k=0}^r C(u, k) · 3^k - 9m·u

3. Compute h(u) = g'(u) (Cochrane Theorem 2 with C·f'/f = 0 since f=const)

4. Factor h(u) = 3^τ · H(u) with H ≢ 0 mod 3. Compute D = degp H+.

5. Apply Theorem 2 bound: |S/p^m| ≤ D · p^{-(m-τ-1)/(2D)}.

PROBLEM IDENTIFIED:
  - g(u) is polynomial of degree r in u
  - g'(u) is polynomial of degree r-1 in u
  - Coefficient of u^j in g'(u) has v_3 ≥ j+1 - v_3((j+1)!)
  - In particular v_3 of constant term = 1, all other terms v_3 ≥ 2
  - Factoring τ=1: H(u) = g'(u)/3 has constant term ≠ 0 mod 3, all other terms ≡ 0 mod 3
  - **D = degp H+ = 0** (constant polynomial mod 3)

CONSEQUENCE:
  - Cochrane Theorem 2 with D=0: bound becomes "vanishing sum" by H(a) ≢ 0 mod p^{m-ℓ-τ}
  - But S_{r,ℓ,ε}(m) is non-vanishing (Kalafatelis treats it as a non-trivial obstruction)
  - **Cochrane Theorem 2 DOES NOT directly apply because our sum is INCOMPLETE**
    (range u ∈ [0, 3^{r-1}-1], modulus 3^{r+1}; Theorem 2 requires complete summation mod p^m)

This script verifies the identification numerically and documents the obstruction.
"""
import sys
import os
from fractions import Fraction
import math

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def binomial_poly_coeffs(k):
    """Return coefficients [a_0, a_1, ..., a_k] of C(u, k) = u(u-1)...(u-k+1)/k! as polynomial in u.
    Returns Fraction coefficients."""
    if k == 0:
        return [Fraction(1)]
    poly = [Fraction(1)]  # represents constant 1
    # multiply by (u - i)/(i+1) iteratively for i = 0, 1, ..., k-1, then divide by k!
    # Actually C(u, k) = product_{i=0}^{k-1} (u-i) / k!
    # Start with 1, multiply by (u - i) k times, then divide by k!
    poly = [Fraction(1)]  # 1
    for i in range(k):
        # multiply poly by (u - i)
        new_poly = [Fraction(0)] * (len(poly) + 1)
        for j in range(len(poly)):
            new_poly[j+1] += poly[j]    # u term
            new_poly[j] += -i * poly[j]  # constant term
        poly = new_poly
    # Divide by k!
    fact = Fraction(math.factorial(k))
    poly = [c / fact for c in poly]
    return poly


def derivative(poly):
    """Polynomial derivative: poly = [a_0, a_1, ...] -> [a_1, 2a_2, 3a_3, ...]."""
    return [Fraction(j) * poly[j] for j in range(1, len(poly))]


def v_3(q):
    """3-adic valuation of a Fraction. Returns ∞ if q == 0."""
    if q == 0:
        return float('inf')
    num, den = q.numerator, q.denominator
    v = 0
    while num % 3 == 0:
        num //= 3
        v += 1
    while den % 3 == 0:
        den //= 3
        v -= 1
    return v


def main():
    print("# Result 78 Step 1: Polynomial identification for Cochrane Theorem 2")
    print()

    # For r = 2, 3, 4: build g(u), compute g'(u), find τ and D
    for r in [2, 3, 4]:
        N_phase = 3**(r+1)  # phase modulus
        N_range = 3**(r-1)  # summation range
        print(f"## Level r = {r}: phase modulus 3^{r+1} = {N_phase}, range 3^{r-1} = {N_range}")
        print()

        # g(u) = c · Σ_{k=0}^r C(u, k) · 3^k - 9m·u
        # where c = c_{ℓ,ε} is a 3-adic unit (we don't fix it; just track 3-adic structure)
        # Build g_coeffs as polynomial in u, with c left symbolic (just compute c=1 case)

        # Actually we want g(u) - 9m·u, but 9m is just shift in constant of g(u). We can ignore m
        # for the H structure analysis (m only shifts constant of g, so shifts constant of g'... wait
        # 9mu is linear in u, so its derivative is 9m which is constant. So 9m shifts the constant of g'.)

        # For c = 1: g(u) = Σ_{k=0}^r C(u, k) · 3^k
        # = 1 + 3u + 9·u(u-1)/2 + 27·u(u-1)(u-2)/6 + ... + 3^r·C(u, r)

        g_coeffs = [Fraction(0)] * (r + 1)  # polynomial of degree r in u
        for k in range(r + 1):
            c_uk = binomial_poly_coeffs(k)  # C(u, k) as poly
            scale = Fraction(3**k)
            for j, coef in enumerate(c_uk):
                g_coeffs[j] += scale * coef

        # g'(u) = derivative
        gp_coeffs = derivative(g_coeffs)
        print(f"  g(u) = c · Σ_{{k=0}}^{r} C(u, k) · 3^k  (degree {len(g_coeffs)-1})")
        print(f"  g'(u) coefficients (c = 1 case):")
        for j, coef in enumerate(gp_coeffs):
            v = v_3(coef)
            print(f"    u^{j}: {coef} = {float(coef):.6f},  v_3 = {v}")
        print()

        # Find τ = min v_3 of coefficients of g'(u)
        v_3_list = [v_3(c) for c in gp_coeffs if c != 0]
        tau = min(v_3_list) if v_3_list else float('inf')
        print(f"  τ = min v_3 of g'(u) coefficients = {tau}")

        # H(u) = g'(u) / 3^τ
        H_coeffs = [c / Fraction(3**tau) for c in gp_coeffs]
        # H mod 3
        H_mod3 = [c.numerator * pow(c.denominator, -1, 3) % 3 if c.denominator % 3 != 0 else 'div0'
                  for c in H_coeffs]
        print(f"  H(u) = g'(u) / 3^{tau}  coefficients mod 3: {H_mod3}")

        # D = deg of H mod 3 (highest non-zero coefficient mod 3)
        D = -1
        for j in range(len(H_mod3) - 1, -1, -1):
            if H_mod3[j] != 0 and H_mod3[j] != 'div0':
                D = j
                break
        print(f"  D = degp H+ = {D}")
        print()

        # Apply Cochrane Theorem 2 if D > 0
        m = r + 1  # Cochrane's m
        if D > 0:
            bound_log = math.log10(D) + (-(m - tau - 1) / (2 * D)) * math.log10(3)
            print(f"  Cochrane Theorem 2: |sum/p^m| ≤ D · p^{{-(m-τ-1)/(2D)}}")
            print(f"  Numerical: {D} · 3^{{{-(m-tau-1)/(2*D):.4f}}}  ≈ 10^{bound_log:.3f}")
            print(f"  Hence |sum| ≤ {D} · 3^{{{(m - (m-tau-1)/(2*D)):.4f}}} (× p^m factor)")
        else:
            print(f"  ⚠ D = {D}: Cochrane Theorem 2 trivializes (sum vanishes if H(a) ≢ 0 mod p)")
            print(f"  But sum is INCOMPLETE (range 3^{r-1} < modulus 3^{r+1}), so trivial vanishing doesn't apply.")
            print(f"  → Cochrane Theorem 2 doesn't directly close eq 190.")
        print()
    print()

    # Document the obstruction
    print("# Obstruction analysis:")
    print()
    print("Cochrane Theorem 2 requires COMPLETE summation Σ_{n mod p^m}.")
    print("Our S_{r,ℓ,ε}(m) is INCOMPLETE: range u ∈ [0, 3^{r-1}-1], modulus 3^{r+1}.")
    print()
    print("The polynomial identification gives g(u) of degree r with:")
    print("  - g(u) = c·Σ_{k=0}^r C(u,k)·3^k - 9m·u")
    print("  - g'(u) has v_3 of coefficients = j+1 - v_3((j+1)!) for u^j (≥ 1 for j ≥ 0)")
    print("  - Factoring τ = min v_3 = 1: H(u) = g'(u)/3 has degp H+ = 0 (constant mod 3)")
    print("  - Cochrane Theorem 2 with D = 0 gives 'sum vanishes' (trivial bound)")
    print()
    print("The INCOMPLETE summation means Cochrane Theorem 2's conclusion doesn't apply.")
    print("Need an INCOMPLETE-sum analog: Cochrane Theorem 1 (smooth modulus, doesn't apply to 3^{r+1})")
    print("or alternative bridging technique.")
    print()

    # Save analysis
    out = os.path.join(OUTDIR, "r78_polynomial_identification.txt")
    with open(out, 'w') as f:
        f.write("Result 78 Step 1: Polynomial identification analysis\n\n")
        f.write("Identification: f(u) = 1, g(u) = c · Σ_{k=0}^r C(u, k) · 3^k - 9m·u\n")
        f.write("h(u) = g'(u), τ = 1, D = degp H+ = 0\n\n")
        f.write("Cochrane Theorem 2 doesn't apply: requires complete sum mod p^m, ours is incomplete.\n")
    print(f"[save] {out}")


if __name__ == "__main__":
    main()
