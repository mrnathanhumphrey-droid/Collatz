"""
forward_q_eigvec_test.py

Test the v3.5 framing hypothesis: BOTH the inverse-tree stationary measure AND
the forward-orbit trajectory measure are leading left-eigenvectors of their
respective natural-density transition matrices on residues mod 2^k.

Inverse-tree:  M_closed (built in inverse_tree_residue_build.py)
Forward-orbit: Q (Syracuse residue chain Path B uses)

If both match their empirical equilibria, the v3.5 framing holds and the
structural anti-correlation (r ≈ -0.20 between the two ratio profiles)
reduces to "different transition matrices on the same state space."
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")


def v2(n):
    if n == 0:
        return 0
    c = 0
    while n & 1 == 0:
        n >>= 1
        c += 1
    return c


def build_Q_forward_syracuse(k):
    """Build Q[r1, r2] on ODD residues r1, r2 mod 2^k for forward Syracuse map.

    For odd m = r + 2^k * h with h uniform-natural-density:
      v(m) = nu_2(3m+1)
      T(m) = (3m+1) / 2^v
    Standard fact: v has Geom(1/2) distribution, P(v=j) = 2^(-j) for j=1,2,...
    Conditional on v=j, T(m) mod 2^k for natural-density h:
      - If v < k: T(m) cycles deterministically through 2^v specific residues
        as h varies through 2^v values; then repeats. Each cycle-residue gets
        equal weight.
      - If v >= k: T(m) = (3m+1)/2^v is dominated by high bits of h; mod 2^k
        is uniform over odd residues.
    """
    M = 1 << k
    odds = [r for r in range(1, M, 2)]
    N = len(odds)
    idx = {r: i for i, r in enumerate(odds)}
    Q = np.zeros((N, N))
    for r in odds:
        i = idx[r]
        n3r1 = 3 * r + 1
        v_r = v2(n3r1)
        if v_r < k:
            # Deterministic v at this resolution.
            q_r = n3r1 >> v_r
            shift = 3 * (1 << (k - v_r))
            for h in range(1 << v_r):
                t = (q_r + shift * h) & (M - 1)
                if t & 1:  # only count odd targets
                    Q[i, idx[t]] += 1.0 / (1 << v_r)
        else:
            # Boundary: v_r >= k. v stochastic.
            # P(v = v_r + j) = 2^(-(j+1)) for j = 0, 1, 2, ...
            # T(m) mod 2^k uniform over odd residues for j large enough.
            # Truncate at j_max = 20.
            for j in range(20):
                p_j = 2.0 ** (-(j + 1))
                # Uniform on odd residues at this approximation level
                Q[i, :] += p_j / N
            # Ensure rows sum to 1
            s = Q[i, :].sum()
            if s > 0:
                Q[i, :] /= s
    return Q, odds


def leading_left_eigvec(Q):
    """Stationary distribution: left eigenvector of Q at lambda=1."""
    eigvals, eigvecs = np.linalg.eig(Q.T)
    idx = np.argmin(np.abs(eigvals - 1.0))
    pi = np.real(eigvecs[:, idx])
    pi = np.abs(pi)
    pi = pi / pi.sum()
    return pi, complex(eigvals[idx])


def main():
    print(f"=" * 70)
    print(f"v3.5 framing test: forward Syracuse Q leading eigvec vs trajectory measure")
    print(f"=" * 70)

    # Test at multiple k to see if spikes emerge at finer resolution
    for k_test in [5, 6, 8]:
        Q_test, odds_test = build_Q_forward_syracuse(k_test)
        pi_test, lam_test = leading_left_eigvec(Q_test)
        pi_ratio_test = pi_test / (1.0 / len(odds_test))
        print(f"\n[k={k_test}, mod {1<<k_test}] Q stationary range: "
              f"min={pi_ratio_test.min():.6f}, max={pi_ratio_test.max():.6f}, "
              f"std={pi_ratio_test.std():.6f}")

    k = 5  # mod 32
    M = 1 << k

    Q, odds = build_Q_forward_syracuse(k)
    print(f"\n[forward Q] mod {M}, {len(odds)} odd states, row sums:")
    rs = Q.sum(axis=1)
    print(f"  min={rs.min():.6f}, max={rs.max():.6f} (should be 1.0)")

    pi, lam = leading_left_eigvec(Q)
    print(f"\n[forward Q stationary] leading left-eigvec at lambda={lam.real:.6f}+{lam.imag:.4f}i")
    print(f"  This is the predicted forward-orbit trajectory measure on odd residues mod {M}.")

    # Express as ratio vs uniform-on-odd-residues (1/16 each)
    pi_ratio = pi / (1.0 / len(odds))

    print(f"\n[predicted vs forward-orbit empirical from agent2_findings.md]")
    fwd_emp = {1: 0.930, 3: 0.981, 5: 1.232, 7: 1.051, 9: 0.981, 11: 0.996,
               13: 0.911, 15: 1.045, 17: 1.082, 19: 0.890, 21: 0.887,
               23: 1.052, 25: 0.896, 27: 1.040, 29: 1.081, 31: 0.946}
    print(f"  {'r':>4} {'Q_eigvec':>10} {'fwd_emp':>10} {'diff':>10}")
    pred = []
    emp = []
    for r in sorted(fwd_emp.keys()):
        i = odds.index(r)
        p = pi_ratio[i]
        e = fwd_emp[r]
        pred.append(p); emp.append(e)
        print(f"  {r:>4} {p:>10.4f} {e:>10.4f} {p-e:>+10.4f}")

    pred = np.array(pred); emp = np.array(emp)
    r_corr = float(np.corrcoef(pred, emp)[0, 1])
    print(f"\n  Pearson r(Q_eigvec, fwd_emp) = {r_corr:+.4f}")
    max_abs_diff = float(np.max(np.abs(pred - emp)))
    print(f"  Max |Q_eigvec - fwd_emp| = {max_abs_diff:.4f}")
    rmse = float(np.sqrt(np.mean((pred - emp) ** 2)))
    print(f"  RMSE = {rmse:.4f}")

    print(f"\n[interpretation]")
    if r_corr > 0.9 and max_abs_diff < 0.1:
        print(f"  Q's leading eigvec MATCHES forward-orbit trajectory measure tightly.")
        print(f"  v3.5 framing CONFIRMED: both inverse-tree and forward-orbit measures")
        print(f"  are leading eigvecs of natural-density transition matrices on the")
        print(f"  same residue space. Their structural anti-correlation reduces to")
        print(f"  M_closed (inverse) vs Q (forward) being structurally different matrices.")
    elif r_corr > 0.5:
        print(f"  Partial match (Pearson r={r_corr:.2f}). v3.5 framing partially holds")
        print(f"  but the predicted Q-eigvec doesn't perfectly recover the forward-orbit")
        print(f"  trajectory measure. Document the gap.")
    else:
        print(f"  No match (Pearson r={r_corr:.2f}). The forward-orbit trajectory measure")
        print(f"  is NOT the leading eigvec of Q at this resolution. v3.5 framing")
        print(f"  fails — forward-orbit needs a different structural object.")


if __name__ == "__main__":
    main()
