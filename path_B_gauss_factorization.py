"""
path_B_gauss_factorization.py — compute F̂(3a) exactly for r = 2, 3 and look for
closed-form phase pattern.

F̂(3a) = Σ_{u=0}^{q-1} e_q(c·4^u - 3au) where q = 3^{r+1}.

By substitution v = 4^u (with u = log_4(v)), v ranges over principal units {1+3Z/q}
(period 3^r), each visited 3 times since u ranges over Z/q with period 3^r in 4^u:

F̂(3a) = 3 · Σ_{v ∈ 1+3Z/q} e_q(c·v - 3a·log_4(v))

We computed |F̂(3a)| = 3√q on its support. The phase φ(a) such that F̂(3a) = 3√q · e^{iφ(a)}
is what we want to extract — if it has clean closed form, eq 190 closure follows.
"""
import sys
import os
import math
import cmath

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def discrete_log_4(v, q, period):
    """Compute u s.t. 4^u ≡ v mod q, u ∈ {0, 1, ..., period-1}."""
    pw = 1
    for u in range(period):
        if pw == v:
            return u
        pw = (pw * 4) % q
    return None


def F_hat_via_gauss(r, c, a):
    """Compute F̂(3a) via Gauss-sum decomposition.

    F̂(3a) = 3 · Σ_{v ∈ principal units} e_q(c·v - 3a·log_4(v))
    """
    q = 3**(r+1)
    period = 3**r
    # Compute log_4 table for principal units
    pu = []  # list of principal unit values, indexed by log
    pw = 1
    for _ in range(period):
        pu.append(pw)
        pw = (pw * 4) % q

    total = complex(0, 0)
    for u in range(period):
        v = pu[u]
        log_v = u  # since pu[u] = 4^u
        phase = ((c * v - 3 * a * log_v) % q) / q
        total += cmath.exp(2j * cmath.pi * phase)
    return 3 * total


def F_hat_direct(r, c, a):
    """Direct: F̂(3a) = Σ_{u=0}^{q-1} e_q(c·4^u - 3a·u)."""
    q = 3**(r+1)
    total = complex(0, 0)
    pw = 1
    for u in range(q):
        phase = ((c * pw - 3 * a * u) % q) / q
        total += cmath.exp(2j * cmath.pi * phase)
        pw = (pw * 4) % q
    return total


def main():
    print("# Path B: explicit Gauss-sum factorization of F̂(3a)")
    print()

    for r in [2, 3]:
        q = 3**(r+1)
        period = 3**r
        sqrt_3q = math.sqrt(9 * q)  # = 3√q
        print(f"## r = {r}: q = {q}, period = {period}, |F̂| = 3√q = {sqrt_3q:.6f}")
        print()

        # Verify equivalence of F_hat_via_gauss and F_hat_direct
        c = 1
        for a in [1, 4, 7]:
            via = F_hat_via_gauss(r, c, a)
            direct = F_hat_direct(r, c, a)
            print(f"  a = {a}: via_gauss = {via:.4f}, direct = {direct:.4f}, "
                  f"diff = {abs(via - direct):.2e}")
        print()

        # For a in support {a ≡ 1 mod 3}, compute F̂(3a) and extract phase φ(a)
        # F̂(3a) = 3√q · e^{iφ(a)}
        # Note: only a ≡ 1 mod 3 in Z/period are in support (others give F̂ = 0)
        print(f"  F̂(3a) for a ≡ 1 mod 3, a ∈ Z/{period} (the support):")
        print(f"    {'a':>4}  {'F̂_real':>12}  {'F̂_imag':>12}  {'|F̂|':>10}  {'phase φ(a)/2π':>16}  {'(period·a) mod q':>18}")
        phase_values = []
        for a in range(period):
            if a % 3 != 1:
                continue
            F = F_hat_via_gauss(r, c=1, a=a)
            mag = abs(F)
            if mag > 0.01:
                phi_norm = (cmath.phase(F) / (2 * cmath.pi)) % 1.0
            else:
                phi_norm = -1
            phase_values.append((a, F, phi_norm))
            # Suggested closed form for phi: φ(a)/2π = some function of a, q
            # First guess: φ(a)/2π ≡ (something · a) / period mod 1
            print(f"    {a:>4}  {F.real:>12.6f}  {F.imag:>12.6f}  {mag:>10.4f}  {phi_norm:>16.6f}")

        # Look for linear pattern in phi: φ(a)/2π = α · a + β mod 1
        # Differences:
        if len(phase_values) >= 2:
            print(f"\n  Phase differences (φ(a_{{i+1}}) - φ(a_i))/2π:")
            for i in range(len(phase_values) - 1):
                a1, _, p1 = phase_values[i]
                a2, _, p2 = phase_values[i+1]
                diff = (p2 - p1) % 1.0
                if diff > 0.5:
                    diff -= 1
                print(f"    a={a1} → a={a2}: Δφ/2π = {diff:+.6f},  (a2 - a1)·1/q = {(a2-a1)/q:+.6f}")
        print()

    # Try to identify the closed form. Hypothesis: φ(a)/2π = (something) · a / q, possibly + constant
    print("# Hypothesis testing: closed-form φ(a)")
    print()

    # For the principal-unit Gauss sum, the standard formula involves:
    # G(χ) := Σ_{v ∈ 1+3Z/q} χ(v) e_q(v) where χ is a multiplicative character on principal units.
    # Our sum has e_q(c·v - 3a·log_4(v)) = e_q(c·v) · χ_a(v) where χ_a(v) = e_q(-3a·log_4(v)) = e_period(-a·log_4(v)).
    # χ_a is a character on principal units (cyclic of order 3^r), with χ_a(4) = e_period(-a).
    # Specifically χ_a is determined by χ_a(4^u) = e_period(-au).
    # Which is exactly the character e_period(-a · log_4(·)).
    #
    # So F̂(3a) = 3 · Σ_v χ_a(v) e_q(c·v) = 3 · G(χ_a, c)  (Gauss sum twisted by c)
    # where G(χ, c) := Σ_v χ(v) e_q(c·v).

    # Standard Gauss-sum theory: |G(χ, c)|² = q for primitive χ. With q = 3^{r+1} and our χ_a primitive
    # iff a coprime to period = 3^r... but we have a ≡ 1 mod 3, so 3 ∤ a, so a coprime to 3 (yes).
    # However: a ∈ {1, 4, 7, ..., period - 2} — does this make χ_a primitive?
    # χ_a is primitive on (1+3Z/q) iff a generates the dual group of (1+3Z/q) ≅ Z/period.
    # χ_a(4^u) = e_period(-au), so χ_a generates the dual iff a ∈ (Z/period)*.
    # For a ≡ 1 mod 3 (so 3 ∤ a in Z/period), a is coprime to period = 3^r iff... 3 ∤ a. ✓
    # So χ_a primitive for all a in our support. Then |G(χ_a, c)| = √q.
    #
    # F̂(3a) = 3 · G(χ_a, 1).
    # |F̂(3a)| = 3 · √q = 3√q. ✓ matches Theorem 78.3.
    #
    # The PHASE of G(χ_a, c) is the "Gauss sum phase" — for principal-unit characters on Z/3^{r+1},
    # this has specific closed form (Davenport-Hasse, etc.).

    # Hypothesis: F̂(3a) = 3√q · e^{iπ/4 · sign} · (Legendre-symbol-like factor) · e_q(c · ψ(a))
    # for some inverse-character function ψ(a).

    # Test: compute F̂(3a) / (3√q) and look at the unit-modulus part.
    print(f"  F̂(3a) / (3√q) for r=2 (the unit-modulus phases):")
    r = 2
    q = 3**(r+1)
    period = 3**r
    sqrt_3q_factor = math.sqrt(9 * q)
    print(f"    {'a':>4}  {'F̂/(3√q) = u(a)':>30}  {'|u(a)|':>8}")
    u_values = []
    for a in range(period):
        if a % 3 != 1:
            continue
        F = F_hat_via_gauss(r, c=1, a=a)
        u = F / sqrt_3q_factor
        u_values.append((a, u))
        print(f"    {a:>4}  {u.real:>+12.6f} + {u.imag:>+12.6f}i  {abs(u):>8.6f}")
    print()

    # Look for pattern: u(1)*u(4)*u(7) = ?
    if len(u_values) == 3:
        prod = u_values[0][1] * u_values[1][1] * u_values[2][1]
        print(f"  u(1) · u(4) · u(7) = {prod:.6f}")
        print(f"  u(1) · u(4) = {(u_values[0][1] * u_values[1][1]):.6f}")
        print(f"  Sum u(1) + u(4) + u(7) = {sum(u for _, u in u_values):.6f}")

    # Check c-dependence
    print(f"\n  F̂(3a, c) varying c at r=2:")
    print(f"  {'c':>4}  {'a=1':>20}  {'a=4':>20}  {'a=7':>20}")
    for c_val in [1, 2, 4, 5, 7, 8]:
        row = [f"{c_val:>4}"]
        for a in [1, 4, 7]:
            F = F_hat_via_gauss(r, c=c_val, a=a)
            row.append(f"{abs(F):.4f}")
        print("  " + "  ".join(row))


if __name__ == "__main__":
    main()
