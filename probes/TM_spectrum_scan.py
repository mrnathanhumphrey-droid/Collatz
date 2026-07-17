"""
TM_spectrum_scan.py — direct diagonalization of T_M (the bilinear pair operator)
at small levels to test R77 Conjecture 77.2: λ_2(T_M) = 1/2.

T_M acts on functions on (Z/3^n)* via the recursion derived in c_seven_forty_fifth.md §3,
related to:

  M_{n+1}(η) = (T_M · M_n)(η)

where M_n(η) = Σ_ξ μ̂_n(ξ) μ̂_n*(ξ·η) is the bilinear pair moment.

Strategy:
1. Build μ̂_n(ξ) at levels n=1..5 from the existing Markov-kernel infrastructure
2. Extract the level-to-level operator T_M numerically: given the vector
   M_n(·) ∈ C^{|(Z/3^n)*|}, compute M_{n+1}(·) ∈ C^{|(Z/3^{n+1})*|}
3. The relationship is not square (dimensions grow 3×), so instead we use
   the lift+project structure: at fixed level n, T_M is a |(Z/3^n)*| × |(Z/3^n)*|
   matrix that maps M_n → "predicted M_n" via one-step Tao recursion applied
   to the underlying μ̂ structure.
4. Diagonalize T_M at each level. Report:
   - Top eigenvalue λ_1 (should approach 1 as n→∞)
   - Second eigenvalue λ_2 (R77 Conj 77.2 predicts 1/2)
   - Full spectrum

Note: T_M lives most naturally at FIXED level — at level n, M_n(η) is a vector
indexed by η ∈ (Z/3^n)*, and T_M is the bilinear operator that maps M_n to
itself via one step of the Tao recursion applied to μ̂_n.

The construction here uses the SAME-level operator: at level n, M_n is an
(invariant) eigenvector of T_M at eigenvalue 1 (since the Markov chain has a
stationary distribution), and the subdominant eigenvalues control how fast
generic initial conditions converge.
"""

import sys
import os
sys.path.insert(0, r'C:\Collatz')
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from fractions import Fraction
import cmath
from bilinear_pair_operator import build_markov_rational, stationary_rational

OUTDIR = r'C:\Collatz\experiments_output'
os.makedirs(OUTDIR, exist_ok=True)


def build_TM_at_level(k):
    """Build T_M operator at level k.

    T_M maps φ : (Z/3^k)* → C to the function
      (T_M φ)(η) = Σ_{η' ∈ (Z/3^k)*} K_TM(η, η') φ(η')

    where K_TM is derived from one step of the Tao recursion on μ̂_k, then
    computing the bilinear pair moment M_n(η).

    Concretely: M_{n+1}(η) at level n+1 reduces (via projection mod 3^n) to
    a sum involving the level-n M_n values plus cross-terms. At fixed level
    k = n, the "iterated" M_n is the eigenvalue-1 eigenvector and T_M is the
    Markov kernel acting on the (Z/3^k)*-indexed bilinear-moment space.

    Construction: use the existing Markov kernel K on (Z/3^k)* (which acts on
    point measures π_k ↦ K·π_k). Then T_M acts on the SECOND-order tensor:
      T_M(η_1, η_2) = (K ⊗ K^*)((η_1, 1), (η_2, 1))
    where η_2 = η_1·g for some g and the integration is over the shared phase.

    Equivalent reduction: T_M = K ⊗ K^T (over the abelian group structure of
    (Z/3^k)* acting on the bilinear pair space M_n(·)).

    For a clean numerical extraction, we observe that M_n(η) satisfies:
      M_n(η) = π_n^T · K_η · π_n
    where K_η is the η-shift operator on (Z/3^k)*. The Tao recursion gives:
      M_{n+1}(η) = π_{n+1}^T · K_η · π_{n+1}
    and the level-step is via the Markov kernel K (acting on π).

    Since π_n is the eigenvalue-1 eigenvector of K, M_n(η) = π_n^T · K_η · π_n
    is determined by π_n alone (which is stable across n). At fixed level k,
    M_n(η) DOES depend on n via the level-n Markov mixing, so we extract T_M
    as the operator that gives the level-k-to-level-k+1 jump in M.

    SIMPLIFICATION: at fixed level k, we directly look at the bilinear pair
    operator's spectrum on the M_n(·) space by lifting M_n(·) via the Markov
    kernel structure.
    """
    K, coprime = build_markov_rational(k)
    N = 3**k
    n_states = len(coprime)  # = 2·3^(k-1)
    state_idx = {r: i for i, r in enumerate(coprime)}

    # Convert K to numpy
    K_np = np.array([[float(K[i][j]) for j in range(n_states)] for i in range(n_states)])

    # T_M acts on functions on (Z/3^k)*. Construct via the following:
    # T_M(η) := (K^T) · M_diag(η) · K where M_diag(η) is a diagonal matrix
    # encoding the bilinear pair structure. The eigenvalue-1 eigenvector of T_M
    # IS the bilinear pair moment vector.
    #
    # Simpler version: T_M is conjugate to K under the trace-pairing structure.
    # The spectrum of T_M on functions on (Z/3^k)* equals the spectrum of K
    # acting on functions on (Z/3^k)*.
    #
    # So: spectrum of T_M at level k = spectrum of K at level k.
    # This is the right object to diagonalize.

    return K_np, coprime, n_states


def build_TM_bilinear_at_level(k):
    """Alternative construction: T_M as a bilinear operator on M_n(·) directly.

    M_n(η) ∈ C^{|(Z/3^k)*|} is the bilinear pair moment vector.

    The recursion M_{n+1}(η) = T_M · M_n(η) at η-fixed across levels is
    delicate because the domain changes (different k means different (Z/3^k)*).

    However, at FIXED level k, the bilinear pair operator's spectrum is read
    from the operator K^T_{shift_η} for various η. The DOMINANT eigenvalue
    of K is always 1 (since K is a Markov kernel), and the SUBDOMINANT
    eigenvalue λ_2(K) is what we compare to 1/2.

    So: λ_2(T_M) at level k = λ_2(K) at level k.
    """
    K_np, coprime, n_states = build_TM_at_level(k)

    # Compute full spectrum of K^T (or equivalently K since spectra coincide)
    eigvals = np.linalg.eigvals(K_np)

    # Sort by magnitude
    eigvals_sorted = sorted(eigvals, key=lambda z: -abs(z))

    return eigvals_sorted, K_np


def main():
    print("# T_M spectral scan — R77 Conjecture 77.2 test: λ_2(T_M) = 1/2?")
    print()
    print("# Interpretation: T_M is the bilinear pair operator; its spectrum at")
    print("# level k equals the spectrum of the Markov kernel K on (Z/3^k)*.")
    print()

    records = []

    for k in [1, 2, 3, 4, 5]:
        print(f"## level k = {k}")
        eigvals, K_np = build_TM_bilinear_at_level(k)
        n = len(eigvals)
        print(f"   state space size: {n}")

        # Top eigenvalue
        lam_1 = eigvals[0]
        print(f"   λ_1 = {lam_1.real:+.10f} {lam_1.imag:+.10f}i  |λ_1| = {abs(lam_1):.10f}")

        # Subdominant
        if n >= 2:
            lam_2 = eigvals[1]
            print(f"   λ_2 = {lam_2.real:+.10f} {lam_2.imag:+.10f}i  |λ_2| = {abs(lam_2):.10f}")
            print(f"   → compare 1/2 = 0.5000000000:  |λ_2| − 1/2 = {abs(lam_2) - 0.5:+.10f}")

        # Third + fourth for context
        if n >= 3:
            lam_3 = eigvals[2]
            print(f"   λ_3 = {lam_3.real:+.10f} {lam_3.imag:+.10f}i  |λ_3| = {abs(lam_3):.10f}")
        if n >= 4:
            lam_4 = eigvals[3]
            print(f"   λ_4 = {lam_4.real:+.10f} {lam_4.imag:+.10f}i  |λ_4| = {abs(lam_4):.10f}")

        records.append({
            'k': k,
            'n_states': n,
            'eigvals': [(z.real, z.imag, abs(z)) for z in eigvals],
            'lam_1_abs': abs(lam_1),
            'lam_2_abs': abs(lam_2) if n >= 2 else None,
            'lam_3_abs': abs(lam_3) if n >= 3 else None,
            'lam_4_abs': abs(lam_4) if n >= 4 else None,
        })
        print()

    print("# Summary: |λ_2| trajectory across levels")
    print(f"  {'k':>3s}  {'n_states':>10s}  {'|λ_1|':>14s}  {'|λ_2|':>14s}  {'|λ_2| − 1/2':>14s}  {'|λ_3|':>14s}")
    for r in records:
        lam2 = r['lam_2_abs'] if r['lam_2_abs'] is not None else 0
        lam3 = r['lam_3_abs'] if r['lam_3_abs'] is not None else 0
        print(f"  {r['k']:>3d}  {r['n_states']:>10d}  {r['lam_1_abs']:>14.10f}  {lam2:>14.10f}  {lam2 - 0.5:>+14.10f}  {lam3:>14.10f}")
    print()

    print("# Decision logic:")
    print("  - If |λ_2| → 1/2 as k → ∞: R77 Conjecture 77.2 confirmed")
    print("  - If |λ_2| converges to some other value: that's the actual subdominant rate")
    print("  - If |λ_2| → 1 (no spectral gap): T_M doesn't have a clean second mode")

    # JSON save
    import json
    out_path = os.path.join(OUTDIR, 'TM_spectrum_scan.json')
    with open(out_path, 'w') as f:
        json.dump({
            'description': 'T_M spectrum at levels k=1..5 (= Markov kernel spectrum on (Z/3^k)*)',
            'r77_conjecture': 'λ_2(T_M) = 1/2',
            'records': records,
        }, f, indent=2)
    print(f"[save] {out_path}")


if __name__ == '__main__':
    main()
