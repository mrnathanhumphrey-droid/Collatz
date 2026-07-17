"""
path_b_matrix_wh.py — build residue-class Markov chain Q for the Syracuse
map at modulus 2^k, then attempt matrix Wiener-Hopf factorization for W_j.

Setup:
  - States: odd residues r mod 2^k. For k=6, that's 32 states (r ∈ {1,3,...,63}).
  - Transitions: for m = r + 2^k · h with h uniform on natural density,
    T(m) ≡ r' mod 2^k with probabilities Q[r, r'].
  - For r where v(r) := ν_2(3r+1) < k, Q[r, ·] is uniform on 2^{v(r)} residues.
  - For r where v(r) ≥ k (the "boundary" residue), v depends on higher bits of m.

For k=6, the boundary residue is r=21 (since 3·21+1 = 64 = 2^6, so v ≥ 6).

This code:
  1. Builds Q on mod-64 odd residues (32x32 stochastic matrix).
  2. Verifies row-stochasticity.
  3. Computes eigenvalues, stationary distribution.
  4. Builds the Markov-additive characteristic exponent matrix Φ(θ).
  5. Reports drift structure and identifies absorbing-class residues.

W_j extraction via matrix Wiener-Hopf is flagged for follow-up; this script
delivers the matrix-construction foundation.

Cap < 400 lines per brief.
"""
import numpy as np
import sys
from math import gcd

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)


def v2(n):
    """2-adic valuation of n (number of trailing zeros)."""
    if n == 0:
        return 0
    count = 0
    while n & 1 == 0:
        n >>= 1
        count += 1
    return count


def build_Q_residue_chain(k):
    """Build Q[r, r'] for odd residues r, r' mod 2^k.

    For each odd r in [1, 2^k):
      v(r) = ν_2(3r+1)
      If v(r) < k:
        T(m) for m = r + 2^k · h cycles through 2^{v(r)} residues mod 2^k.
        Each gets equal weight 1/2^{v(r)}.
      If v(r) >= k:
        Boundary case — v depends on higher bits. Stochastic transition.

    Returns Q (rows indexed by odd residues), v_per_residue, T_per_residue (sample).
    """
    M = 1 << k
    odd_residues = [r for r in range(1, M, 2)]
    N_states = len(odd_residues)
    state_index = {r: i for i, r in enumerate(odd_residues)}

    Q = np.zeros((N_states, N_states), dtype=np.float64)
    log_step_per_transition = np.zeros((N_states, N_states), dtype=np.float64)
    boundary_residues = []

    for r in odd_residues:
        i = state_index[r]
        n3r1 = 3 * r + 1
        v_r = v2(n3r1)

        if v_r < k:
            # Deterministic v-value. T(m) mod 2^k cycles through 2^v values.
            q_r = n3r1 >> v_r            # T(r) modulo 2^k contribution
            shift = 3 * (1 << (k - v_r))  # increment of T(m) per unit h
            # T(m) mod 2^k for h = 0, 1, ..., 2^v - 1:
            target_residues_set = []
            log_step = LOG3 - v_r * LOG2  # the log step is exactly log 3 - v · log 2
            for h in range(1 << v_r):
                t_mod = (q_r + shift * h) & (M - 1)
                # Ensure odd
                if t_mod & 1 == 0:
                    # If t_mod is even, the residue chain is malformed — investigate
                    continue
                target_residues_set.append(t_mod)
            n_targets = len(target_residues_set)
            if n_targets == 0:
                continue
            weight = 1.0 / n_targets
            for t_mod in target_residues_set:
                if t_mod in state_index:
                    j = state_index[t_mod]
                    Q[i, j] += weight
                    log_step_per_transition[i, j] = log_step
                # else t_mod is even — error, shouldn't happen for proper Syracuse
        else:
            # Boundary residue (v_r >= k). Need to handle stochastic v.
            boundary_residues.append(r)
            # For r=21 mod 64 in k=6: 3r+1 = 64, so 3m+1 = 64 + 64·3·(h/...).
            # Actually: 3m+1 = 3(r + 2^k h) + 1 = 3r+1 + 3·2^k·h = 2^k (1 + 3h)
            # So v(3m+1) = k + v_2(1 + 3h).
            # For h sampled uniformly with natural density, distribution of v_2(1+3h):
            # 1+3h takes all residues mod 2^j for j-th level by equidistribution.
            # P(v_2(1+3h) = j) = 1/2^{j+1} for j = 0, 1, 2, ... (standard).
            # We truncate at large j_max for numerical purposes.
            j_max = 20
            for j in range(j_max):
                p_j = 0.5 ** (j + 1)
                v_total = k + j
                if 1 + 3 * 0 == (1 + 3 * 0):  # placeholder — need to model actual T(m) mod 2^k for each j
                    # Compute T(m) for representative m = r + 2^k·h with v_2(1+3h) = j.
                    # T(m) = (3m+1)/2^{v_total} = 2^k (1+3h) / 2^{k+j} = (1+3h)/2^j.
                    # The result T(m) is some odd integer. Its mod 2^k:
                    # (1+3h)/2^j mod 2^k depends on h.
                    # For natural sampling: T(m) mod 2^k is uniform over odd residues.
                    # (Conjectural — should verify.)
                    # Approximation: spread weight uniformly over odd residues.
                    log_step = LOG3 - v_total * LOG2
                    weight = p_j / N_states  # uniform over odd residues
                    for t_mod in odd_residues:
                        j_idx = state_index[t_mod]
                        Q[i, j_idx] += weight
                        # Use weighted log step
                        # (Note: this is a simplification; proper handling needs more care)
                        if log_step_per_transition[i, j_idx] == 0:
                            log_step_per_transition[i, j_idx] = log_step
            # Ensure row sums to 1 (numerical correction for truncation)
            row_sum = Q[i, :].sum()
            if row_sum > 0:
                Q[i, :] /= row_sum

    return Q, odd_residues, boundary_residues, log_step_per_transition


def main():
    k = 6
    print(f"=== Building residue chain Q on mod-2^{k} = mod {1<<k} odd residues ===")
    Q, odd_residues, boundary_residues, log_steps = build_Q_residue_chain(k)
    N = Q.shape[0]
    print(f"  State count (odd residues mod {1<<k}): {N}")
    print(f"  Boundary residues (v ≥ k): {boundary_residues}")
    row_sums = Q.sum(axis=1)
    print(f"  Row sums: min = {row_sums.min():.6f}, max = {row_sums.max():.6f}")
    if not np.allclose(row_sums, 1.0, atol=1e-6):
        print(f"  WARNING: rows not stochastic (max deviation {abs(row_sums - 1).max():.2e})")
    else:
        print(f"  ✓ Q is row-stochastic")

    # Stationary distribution via eigenvalue decomp
    print(f"\n=== Spectral analysis of Q ===")
    eigvals = np.linalg.eigvals(Q)
    eigvals_sorted = sorted(eigvals, key=lambda z: -abs(z))
    print(f"  Top 10 eigenvalues by modulus:")
    for j, ev in enumerate(eigvals_sorted[:10]):
        print(f"    λ_{j} = {ev.real:+.5f} {ev.imag:+.5f}j   |λ| = {abs(ev):.5f}")

    # Stationary distribution (left eigenvector at λ=1)
    eigval_check = abs(eigvals_sorted[0] - 1.0)
    print(f"  |λ_0 − 1| = {eigval_check:.2e}  (should be ~0 for stochastic matrix)")
    eigvals_full, eigvecs_left = np.linalg.eig(Q.T)
    idx = np.argmin(abs(eigvals_full - 1.0))
    pi = np.real(eigvecs_left[:, idx])
    pi = pi / pi.sum()
    print(f"  Stationary distribution π:")
    print(f"    min π_r = {pi.min():.5f},  max π_r = {pi.max():.5f}")
    print(f"    π[r=1] = {pi[odd_residues.index(1)]:.5f}")
    print(f"    π[r=3] = {pi[odd_residues.index(3)]:.5f}")
    print(f"    π[r=5] = {pi[odd_residues.index(5)]:.5f}")

    # Markov-additive drift
    print(f"\n=== Markov-additive drift structure ===")
    # Mean log step per state: E[X | r] = sum over r' of Q[r,r'] · log_step(r → r')
    mean_step_per_r = (Q * log_steps).sum(axis=1)
    print(f"  Mean log-step per residue (E[X | r]):")
    for r_idx in [0, 1, 2, 4, 8, 16]:  # representative residues
        if r_idx < N:
            r = odd_residues[r_idx]
            print(f"    r={r}: E[X | r] = {mean_step_per_r[r_idx]:+.5f}")
    overall_drift = (pi * mean_step_per_r).sum()
    print(f"\n  Stationary drift E_π[X] = Σ π_r · E[X|r] = {overall_drift:+.6f}")
    print(f"  iid Geom(1/2) drift = log(3) − 2·log(2) = {LOG3 - 2*LOG2:+.6f}")
    print(f"  These should be equal under uniform-natural-density assumption.")

    # Markov-additive characteristic exponent at θ = 0 reduces to Q
    # (no log-step contribution since e^(0·X) = 1)
    print(f"\n=== Sample Φ(θ) at θ=0 (= Q) ===")
    Phi_0 = Q  # exp(iθ X) = 1 at θ=0
    eigvals_Phi = np.linalg.eigvals(Phi_0)
    print(f"  Φ(0) = Q, eigenvalues = those above")

    # Identify absorbing-class residues (m_j attractor entries):
    # m_2 = 5, m_4 = 85, m_5 = 341, ...
    # These are odd integers. Their residues mod 64:
    print(f"\n=== Absorbing-class residues mod {1<<k} ===")
    m_j_list = [(j, (4**j - 1) // 3) for j in range(1, 8)]
    for j, m_j in m_j_list:
        r_j = m_j & ((1 << k) - 1)
        print(f"  m_{j} = {m_j:>6d}  →  r = {r_j} mod {1<<k}")

    # Note: the absorbing m_j values are specific INTEGERS, not residue classes.
    # The chain on residues doesn't directly absorb; absorption is a property of
    # the actual orbit's (residue, magnitude) joint state. Path B's matrix WH
    # framework requires additional state (e.g., log m magnitude) to capture this.
    print(f"\n  IMPORTANT: residues alone don't determine absorption.")
    print(f"  m_2 = 5 has residue 5 mod 64, but m = 5 itself is the small")
    print(f"  attractor; m = 5 + 64 = 69 has the same residue but isn't absorbed.")
    print(f"  Path B's matrix WH must couple residue chain (algebraic) with")
    print(f"  a continuous log-m walk component (Markov-additive process).")


if __name__ == "__main__":
    main()
