"""
cross_freq_compute.py — derivation probe for cross-frequency bilinear closure.

Question: does the off-diagonal v ≠ v' contribution of Tao's recursion to
P_{n+1}^{ab}(c) reduce to a linear combination of {P_n^{ab}(c)} (same-frequency
class-resolved moments at level n), or does it require an enlarged span?

Step 1 (analytic, done in the writeup): the cross-frequency contribution at fixed
pair (v, v') reduces via lift-fiber orthogonality (j=0,1,2 sum) and unit shuffle
to a function only of g = v' - v (for v < v', both even). Specifically:

    Q_n^{++}(c; v, v') = 3 · X_n^{++}(c; g),    g = v' - v
    X_n^{ab}(c; g) := Σ_{s ∈ (Z/3^n)× ∩ c-class}
                          e^{-2πi s · ẽ_g / 3^n}
                          · μ̂_n^a(s) · μ̂_n^{b*}(s · 2^{-g} mod 3^n)
    ẽ_g := (1 - 2^{-g}) / 3   mod 3^n
           (well-defined because 3 | (1 - 2^{-g}) for g even)

Step 2 (computational, this script): pick small n (n = 2, n = 3), compute X_n^{ab}(c; g)
for the leading non-trivial pairs (g = 2, g = 4, ...) and check whether each lies in
the linear span of {P_n^{ab}(c)} entries.

If X spans an enlarged space, we have H_CROSS_CLOSES_ON_ENLARGED_SPAN.
If X happens to lie inside span{P_n}, then H_CROSS_CLOSES_ON_SAME_FREQ (unexpected).

ALL ARITHMETIC OVER Q[ω] (ω = e^{2πi/3}); the partial sum reduces to
3^{n-1} · Σ ... e^{-2πi c·h/3^n} = ω^{-c·m}-weighted lattice sums when h ≡ 0 mod 3^{n-1}.

Setup follows R66's class flow: r' ≡ 2^{-v} mod 3 ⇒ v even gives +; v odd gives -.
Hence for P^{++} computation only v, v' both ≥ 2 even contribute; (v, v') = (2, 4), (2, 6), (4, 6), ...
"""
from __future__ import annotations
import sys
from fractions import Fraction
import cmath
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")


# ====================================================================
# Markov chain π_n exactly over Q
# ====================================================================
def build_markov_rational(k: int):
    N = 3 ** k
    M = 2 * 3 ** (k - 1)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime_states = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    n = len(coprime_states)
    K = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    Z_v = Fraction(2 ** M - 1, 2 ** M)
    for r in coprime_states:
        for r_v in range(1, M + 1):
            p = Fraction(1, 2 ** r_v) / Z_v
            target = ((3 * r + 1) * powers_inv2[r_v - 1]) % N
            K[state_idx[r]][state_idx[target]] += p
    return K, coprime_states


def stationary_rational(K):
    n = len(K)
    A = [[K[j][i] - (Fraction(1) if i == j else Fraction(0)) for j in range(n)] for i in range(n)]
    A[n - 1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n - 1] = Fraction(1)
    for col in range(n):
        pivot = -1
        for row in range(col, n):
            if A[row][col] != 0:
                pivot = row
                break
        if pivot == -1:
            raise ValueError(f"Singular at col {col}")
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            b[col], b[pivot] = b[pivot], b[col]
        piv = A[col][col]
        for j in range(col, n):
            A[col][j] /= piv
        b[col] /= piv
        for row in range(n):
            if row != col and A[row][col] != 0:
                factor = A[row][col]
                for j in range(col, n):
                    A[row][j] -= factor * A[col][j]
                b[row] -= factor * b[col]
    return b


# ====================================================================
# class-resolved Fourier characters μ̂_n^a(s) over float-complex
# (the algebraic structure of the question doesn't need Q[ω]; we will
#  use float to numerically diagnose whether the cross-frequency objects
#  lie in span{P_n^{ab}(c)}, which is a linear-algebra question)
# ====================================================================
def char_funcs_class_resolved(pi_q, coprime, k):
    """Compute μ̂_n^+(s) and μ̂_n^-(s) for all s ∈ Z/3^k, as floats."""
    N = 3 ** k
    pi_f = {coprime[i]: float(pi_q[i]) for i in range(len(coprime))}
    mu_plus = {s: complex(0) for s in range(N)}
    mu_minus = {s: complex(0) for s in range(N)}
    for s in range(N):
        for r in coprime:
            phase = cmath.exp(-2j * cmath.pi * r * s / N)
            if r % 3 == 1:
                mu_plus[s] += pi_f[r] * phase
            else:
                mu_minus[s] += pi_f[r] * phase
    return mu_plus, mu_minus


def compute_P_table(mu_plus, mu_minus, coprime):
    """Compute P^{ab}(c) for a,b in {+,-}, c in {1,2}.
    P^{ab}(c) := Σ_{ξ ∈ (Z/3^n)*, ξ ≡ c mod 3} μ̂_n^a(ξ) μ̂_n^{b*}(ξ).
    Returns dict keyed by (a,b,c)."""
    P = {(a, b, c): complex(0) for a in '+-' for b in '+-' for c in (1, 2)}
    for s in coprime:
        c = s % 3
        for a, mu_a in [('+', mu_plus), ('-', mu_minus)]:
            for b, mu_b in [('+', mu_plus), ('-', mu_minus)]:
                P[(a, b, c)] += mu_a[s] * mu_b[s].conjugate()
    return P


def compute_X_cross(mu_plus, mu_minus, coprime, k, g, c_target):
    """Compute the cross-frequency object
       X_n^{ab}(c; g) := Σ_{s ∈ (Z/3^n)*, s ≡ c_target mod 3}
                            e^{-2πi s·ẽ_g/3^n}
                            · μ̂_n^a(s) · μ̂_n^{b*}(s·2^{-g} mod 3^n)
    for (a, b) ∈ {+,-}², g even ≥ 2.

    ẽ_g := (1 - 2^{-g}) / 3 mod 3^n  — well-defined since 3 | (1 - 2^{-g}).
    """
    N = 3 ** k
    # Compute 2^{-g} mod 3^n, then (1 - 2^{-g})/3 mod 3^n
    inv2 = pow(2, -1, N)
    inv2_g = pow(inv2, g, N)
    one_minus = (1 - inv2_g) % N
    # one_minus should be divisible by 3
    assert one_minus % 3 == 0, f"ẽ_g not integer: g={g}, 1-2^-g={one_minus} mod {N}"
    e_tilde = (one_minus // 3) % N

    X = {(a, b): complex(0) for a in '+-' for b in '+-'}
    for s in coprime:
        if s % 3 != c_target:
            continue
        # phase e^{-2πi s · ẽ_g / 3^n}
        phase = cmath.exp(-2j * cmath.pi * s * e_tilde / N)
        s_shift = (s * inv2_g) % N
        for a, mu_a in [('+', mu_plus), ('-', mu_minus)]:
            for b, mu_b in [('+', mu_plus), ('-', mu_minus)]:
                X[(a, b)] += phase * mu_a[s] * mu_b[s_shift].conjugate()
    return X, e_tilde


def basis_P_real(P):
    """Build the real R^8 basis from P:
    [P^++(1), P^++(2), P^--(1), P^--(2),
     Re P^+-(1), Im P^+-(1), Re P^+-(2), Im P^+-(2)]
    (P^-+ is conjugate of P^+-, automatic.)"""
    return np.array([
        P[('+', '+', 1)].real,
        P[('+', '+', 2)].real,
        P[('-', '-', 1)].real,
        P[('-', '-', 2)].real,
        P[('+', '-', 1)].real,
        P[('+', '-', 1)].imag,
        P[('+', '-', 2)].real,
        P[('+', '-', 2)].imag,
    ])


def X_to_real8(X):
    """Likewise pack X^{ab} into an 8-real vector."""
    return np.array([
        X[('+', '+')].real,
        X[('+', '+')].imag,
        X[('-', '-')].real,
        X[('-', '-')].imag,
        X[('+', '-')].real,
        X[('+', '-')].imag,
        X[('-', '+')].real,
        X[('-', '+')].imag,
    ])


def main():
    print("# cross_freq_compute.py — diagnostic for cross-frequency closure")
    print()
    print("# Question: do X_n^{ab}(c; g) (cross-frequency bilinear residues, parameterized")
    print("#           by g = v' - v for v, v' both even) lie in span{P_n^{ab}(c)}?")
    print()

    for k in [2, 3]:
        print(f"## Level n = {k}")
        print("-" * 60)
        K, coprime = build_markov_rational(k)
        pi_q = stationary_rational(K)
        mu_plus, mu_minus = char_funcs_class_resolved(pi_q, coprime, k)
        P = compute_P_table(mu_plus, mu_minus, coprime)

        # Print P table
        print(f"  P^{{ab}}(c) table (n={k}):")
        for c in (1, 2):
            for a in '+-':
                for b in '+-':
                    val = P[(a, b, c)]
                    print(f"    P^{{{a}{b}}}({c}) = {val.real:+.10e}  {val.imag:+.4e}i")
        print()

        # Structural collapse check (R76 §11): P^{+-} = 0 and P^++(1) = P^++(2), P^--(1) = P^--(2)
        for c in (1, 2):
            for a in '+-':
                val = P[('+', '-', c)]
                if abs(val) > 1e-8:
                    print(f"  WARNING: P^+-({c}) not zero: {val}")
        same_plus = abs(P[('+', '+', 1)] - P[('+', '+', 2)])
        same_minus = abs(P[('-', '-', 1)] - P[('-', '-', 2)])
        print(f"  Class-c symmetry: |P^++(1) - P^++(2)| = {same_plus:.2e}, "
              f"|P^--(1) - P^--(2)| = {same_minus:.2e}")
        print()

        # Try the cross-frequency objects for g = 2, 4, ... (even non-zero)
        # We're computing X for v, v' both even, so g = v' - v ∈ {2, 4, 6, ...}
        print(f"  Cross-frequency moments X_n^{{ab}}(c; g) for g ∈ {{2, 4, 6}}:")
        X_data = {}
        for g in [2, 4, 6]:
            try:
                for c_t in (1, 2):
                    X, e_tilde = compute_X_cross(mu_plus, mu_minus, coprime, k, g, c_t)
                    X_data[(g, c_t)] = X
                    print(f"    g={g}, c={c_t}, ẽ_g={e_tilde} mod {3**k}:")
                    for ab in [('+', '+'), ('-', '-'), ('+', '-'), ('-', '+')]:
                        v = X[ab]
                        print(f"      X^{{{ab[0]}{ab[1]}}} = {v.real:+.6e}  {v.imag:+.4e}i")
            except AssertionError as e:
                print(f"    g={g}: SKIP ({e})")
        print()

        # Now the linear-algebra question: do X_n^{ab}(c; g) live in span{P_n^{ab}(c)}?
        # We have 8 P entries (from class-resolved + c). For each (g, c_t) we have 8
        # X entries. Stack them as columns of a matrix; rank > 8 ⇒ enlarged span needed.
        # Note span{P_n^{ab}(c)}_real is itself ≤ 8-dimensional.

        # Build augmented matrix [P_vector | X_g=2,c=1 | X_g=2,c=2 | ...]
        # Each as a real 8-vector (real and imag parts as needed).
        cols = []
        col_names = []
        cols.append(basis_P_real(P))
        col_names.append("P-basis (single vector)")
        for (g, c_t), X in X_data.items():
            cols.append(X_to_real8(X))
            col_names.append(f"X_g={g}_c={c_t}")
        M = np.column_stack(cols)
        # Rank diagnostics
        rank_full = np.linalg.matrix_rank(M, tol=1e-6)
        rank_P_only = np.linalg.matrix_rank(cols[0].reshape(-1, 1), tol=1e-6)
        print(f"  Augmented [P | X_g_c columns]: shape={M.shape}, rank={rank_full}")
        print(f"  P-only rank: {rank_P_only}")
        print()

        # Better diagnostic: span of {P^{ab}(c) entries} as an 8-dim object — the
        # 'span' here is the linear span over R of these 8 real numbers, which is at
        # most 1-dim. So we should compare span of {entries of P, indexed by (a,b,c)}
        # as basis vectors of R^... actually, the P-basis is naturally 8-dim:
        # the 8 'coordinates' of the moment tensor. So the natural ambient space is R^8.
        # The question becomes: does each X (cross-frequency object) project to a
        # vector that lies in span{the 4 basis P^{ab}(c) tensors at level n}?
        #
        # Reframing: the question is whether the LINEAR MAP (mu_n^+, mu_n^-) ⇒ Q-thing
        # only sees the 8 P-coordinates, or whether it sees additional data.
        # Since X is a different functional of (mu_n^+, mu_n^-), they are generally
        # *different functionals*, hence not in the P-span as functionals on the
        # underlying ((Z/3^n)^×)_+- structure. Direct check: vary π_n while preserving
        # the 8 P-values and check if X changes.
        # Practical implementation deferred.

        # However we can at this n compute the magnitudes — if X is non-trivial in
        # magnitude relative to P, that's evidence the cross-frequency object carries
        # independent information.
        if X_data:
            P_mag = np.linalg.norm(basis_P_real(P))
            X_mags = {key: np.linalg.norm(X_to_real8(X)) for key, X in X_data.items()}
            print(f"  |P-vector| = {P_mag:.4e}")
            for key, mag in X_mags.items():
                print(f"  |X_{key}| = {mag:.4e}  (relative to |P|: {mag/P_mag:.3e})")
        print()
        print()


if __name__ == "__main__":
    main()
