"""
l2_flattening_probe.py — Empirical L^2-flattening probe for the Syracuse stationary measure μ_n.

For the Markov chain on (Z/3^n Z)* with transition K_n[r→s]
the stationary distribution π_n on (Z/3^n Z)* is computed exactly in Q.
We then check whether the convolution π_n * π_n flattens in L^2:

  ||π_n * π_n||_2^2  versus  ||π_n||_2^4 / (effective support).

L^2-flattening (in the Khalil / BKS sense, formulated abstractly) means
that repeated self-convolution drives the L^2 mass toward uniform —
quantitatively, that ||μ * μ||_2 is significantly smaller than ||μ||_2.

In the discrete setting on a finite abelian group G of size N:
  - Uniform measure u: u(x) = 1/N for all x. ||u||_2^2 = 1/N.
  - For a probability measure μ on G, ||μ||_2^2 ≥ 1/N (Cauchy-Schwarz).
  - Convolution u * μ = u for any probability μ.
  - "Flattening" means ||μ^{*k}||_2^2 → 1/N (or close), i.e. equilibrating to uniform.

A natural one-step L^2-flattening estimate at scale 3^{-n} is:
  ||π_n * π_n||_2  ≤  c · ||π_n||_2 · (something < 1)
or in the Khalil formulation,
  ||π_n * π_n||_q ≤ N^{-α} · ||π_n||_q^2-something
for q=2, expressing dimension increase under self-convolution.

We compute:
  E_n   = ||π_n||_2^2 = Σ_r π_n(r)^2           (collision probability / Renyi-2)
  E2_n  = ||π_n * π_n||_2^2 = Σ_s (π_n * π_n)(s)^2
  U_n   = 1 / |(Z/3^n Z)*| = 1 / (2 · 3^(n-1))  (uniform collision probability)

The classical Markov-chain mixing identity (random-walk-on-group L^2 decay):
  ||π_n * π_n - u||_2 ≤ (something < 1) · ||π_n - u||_2
would be the natural L^2-flattening estimate at scale n.

We report:
  - E_n, E2_n, U_n
  - ratio E2_n / E_n^2 (should be ≤ 1, smaller = stronger flattening)
  - ratio (E_n - U_n) / U_n (excess over uniform)
  - ratio (E2_n - U_n) / (E_n - U_n) (one-step convolution flattening factor)
"""

import sys
import os
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from s_infinity_exact import build_markov_rational, stationary_rational


def conv_on_Z3n(pi, coprime_states, N):
    """Convolution of pi (supported on coprime_states ⊂ Z/N) with itself.
    Returns a list of length N giving (pi * pi)(s) for s in 0..N-1, as Fractions.
    """
    out = [Fraction(0)] * N
    n = len(coprime_states)
    for i in range(n):
        if pi[i] == 0:
            continue
        ri = coprime_states[i]
        for j in range(n):
            if pi[j] == 0:
                continue
            rj = coprime_states[j]
            s = (ri + rj) % N
            out[s] += pi[i] * pi[j]
    return out


def main():
    print("# L^2-flattening probe for Syracuse stationary measure π_n on (Z/3^n Z)*")
    print("# Convention: ||μ||_2^2 = Σ_x μ(x)^2  (counting measure on Z/3^n Z)")
    print()
    print(f"{'n':>3}  {'|G_n*|':>8}  {'E_n=||π_n||_2^2':>22}  {'E2_n=||π*π||_2^2':>22}  {'U_n=1/|G*|':>14}  {'E_n/U_n':>12}  {'E2_n/U_n':>12}  {'flat ratio':>12}")
    print()

    for n in [1, 2, 3, 4, 5]:
        N = 3**n
        size_Gstar = 2 * 3**(n - 1)
        U_n = Fraction(1, size_Gstar)

        K, coprime = build_markov_rational(n)
        pi = stationary_rational(K)

        E_n = sum(p * p for p in pi)

        pi_conv = conv_on_Z3n(pi, coprime, N)
        E2_n = sum(x * x for x in pi_conv)

        # Flattening factor: how much does one convolution reduce the L^2 excess over uniform?
        # Note: π * π is supported on Z/3^n Z (not just coprime), and its uniform would be 1/N over Z/N
        # but π is supported on G* of size 2·3^(n-1), so π*π is supported on Z/3^n Z with effective support 3^n.
        # Compare to "uniform on the convolution support".
        # Effective L^2 reduction:
        ratio_flat = float(E2_n) / float(E_n)  # ||π*π||_2^2 / ||π||_2^2 — pure ratio
        excess_E = float(E_n) - float(U_n)
        # For comparison: uniform-on-Z/N collision is 1/N
        U_N = Fraction(1, N)
        excess_E2 = float(E2_n) - float(U_N)
        if excess_E > 0:
            flat_excess = excess_E2 / excess_E
        else:
            flat_excess = float("nan")

        print(
            f"{n:>3}  {size_Gstar:>8}  "
            f"{float(E_n):>22.10e}  "
            f"{float(E2_n):>22.10e}  "
            f"{float(U_n):>14.6e}  "
            f"{float(E_n) / float(U_n):>12.6f}  "
            f"{float(E2_n) / float(U_N):>12.6f}  "
            f"{flat_excess:>12.6f}"
        )

    print()
    print("# Interpretation:")
    print("# E_n/U_n = 1 would mean π_n IS uniform on G*. Larger ratio = farther from uniform.")
    print("# E2_n/U_N = 1 would mean π_n * π_n is uniform on Z/3^n Z. Larger = farther.")
    print("# 'flat ratio' = excess(E2) / excess(E): < 1 indicates flattening under one convolution.")
    print("# An L^2-flattening estimate at scale 3^{-n} would require this ratio to decay quantitatively in n.")


if __name__ == "__main__":
    main()
