"""
dark_subspace_probe.py — Dark-subspace classification of Syracuse's DWM quantum trajectory
per Benoist-Pellegrini-Szczepanek 2024 (arXiv 2409.18655).

A dark subspace D ⊂ H = L²((Z/3^n)*) is a subspace invariant under ALL Kraus operators
in the support of µ:
  M_v^{(j, b_prior)} D ⊂ D
for all v ∈ {1, 2, 3, ...}, j ∈ {1, 2, ...}, b_prior ∈ {0, 1, 2, ...}.

Equivalently (when D and D^perp are both invariant), the projection P_D onto D
commutes with all M_v^{(j, b_prior)}:
  [P_D, M_v^{(j, b_prior)}] = 0  for all (v, j, b_prior).

The COMMUTANT A' := {X : [X, M] = 0 for all M in the family} is the natural object.
- If A' = C·I (one-dimensional), Syracuse is FULLY IRREDUCIBLE — no non-trivial dark
  subspaces; unique invariant measure.
- If dim(A') > 1, non-trivial dark subspaces exist; the invariant-measure classification
  is non-trivial.

This probe:
1. Builds the Kraus operator family at n=2, 3 (Z/9, Z/27 coprime states).
2. Computes the joint commutant by finding null space of:
     T(X) := stacked [X, M_v^{(j, b_prior)}] vectorized for all (v, j, b_prior).
3. Reports dim(commutant) and identifies any non-trivial central projections.

If dim(commutant) > 1: factorize via spectral decomposition of central projections,
identify dark subspaces explicitly.
"""
import sys, os, json, cmath, time
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def build_kraus_family(n, V_MAX=16, J_MAX=3):
    """Build all Kraus operators M_v^{(j, b_prior)} acting on V_n = L²((Z/3^n)*).
    Returns list of (label, M_matrix) tuples.

    M_v^{(j, b_prior)}(f)(xi) = 2^{-v/2} · A_v^{(j)}(xi, b_prior) · f(xi · 2^{-v} mod 3^n)

    where A_v^{(j)}(xi, b_prior) = e^{-2*pi*i * xi * x_j(b_prior) * 2^{-v} / 3^n}
    and x_j(b_prior) = 3^{2j-2} * 2^{-b_prior} mod 3^n.
    """
    N = 3 ** n
    inv2 = pow(2, -1, N)
    coprime = [r for r in range(N) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(coprime)}
    n_dim = len(coprime)
    two_pi_over_N = 2.0 * np.pi / N

    # Enumerate distinct x_j(b_prior) values
    # j=1: x_1(0) = 1 (no past)
    # j>=2: x_j(b_prior) = 3^{2j-2} * 2^{-b_prior} mod 3^n
    # b_prior takes values determined by past v's; just enumerate small range
    family = []
    seen_keys = set()
    for j in range(1, J_MAX + 1):
        b_max = 3 * V_MAX  # enough to cycle through all possible x_j values
        for b_prior in range(0, b_max + 1):
            x_j = pow(3, 2 * j - 2, N) * pow(inv2, b_prior, N) % N
            for v in range(1, V_MAX + 1):
                pow_inv2_v = pow(inv2, v, N)
                # x_j_eff = x_j * 2^{-v} mod N (in the phase)
                x_phase = (x_j * pow_inv2_v) % N
                # Key by (x_phase, v) since these determine M_v fully
                # (the shift sigma_{-v} is determined by v; phase depends on x_j*2^{-v} mod N)
                key = (x_phase, v)
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                M = np.zeros((n_dim, n_dim), dtype=complex)
                weight = 2.0 ** (-v / 2)
                for xi in coprime:
                    target_xi = (xi * pow_inv2_v) % N
                    if target_xi % 3 == 0 or target_xi not in idx:
                        continue
                    phase = cmath.exp(-1j * two_pi_over_N * xi * x_phase)
                    M[idx[xi], idx[target_xi]] = weight * phase
                family.append((f"M_v={v}_xphase={x_phase}_j={j}_b={b_prior}", M))
    return family, coprime


def compute_commutant_dim(family, n_dim, tol=1e-9):
    """Compute the dimension of the commutant A' = {X : [X, M] = 0 for all M}.

    Stack the linearization X -> [X, M] for each M as rows of a big matrix T.
    The commutant is the null space of T.
    """
    K = len(family)
    # X is n_dim x n_dim. Vectorize as n_dim^2 vector.
    # For each M, [X, M] = X M - M X is a linear map n_dim^2 -> n_dim^2.
    # T_M(vec X) = vec(X M - M X).
    # Stack T_M for all M: T is (K * n_dim^2) x n_dim^2.
    T_rows = []
    for label, M in family:
        # Vec(XM) = (M^T ⊗ I) vec(X), Vec(MX) = (I ⊗ M) vec(X).
        # So vec([X, M]) = (M^T ⊗ I − I ⊗ M) vec(X).
        Mt = M.T
        I = np.eye(n_dim, dtype=complex)
        T_M = np.kron(Mt, I) - np.kron(I, M)
        T_rows.append(T_M)
    T = np.vstack(T_rows)

    # Compute null space via SVD
    U, sigma, Vh = np.linalg.svd(T, full_matrices=False)
    # Null space: rows of Vh corresponding to sigma < tol
    null_mask = sigma < tol * max(sigma.max(), 1.0)
    null_dim = int(np.sum(null_mask))
    null_basis = Vh[null_mask, :].conj().T  # columns are null vectors
    return null_dim, null_basis, sigma


def analyze_commutant(null_basis, n_dim):
    """Given basis for the commutant (each column is a n_dim^2 vector = vectorized matrix),
    find the central projections / dark subspaces.
    """
    K = null_basis.shape[1]
    print(f"  commutant dim = {K}")
    # Reshape each null vector to a n_dim x n_dim matrix
    matrices = [null_basis[:, k].reshape(n_dim, n_dim) for k in range(K)]
    # Check: do they include the identity? (Should always be in the commutant)
    I_norm = np.linalg.norm(np.eye(n_dim))
    for k, A in enumerate(matrices):
        # Find component along identity
        c_I = np.trace(A) / n_dim
        A_minus_cI = A - c_I * np.eye(n_dim)
        # Print info
        print(f"  commutant basis [{k}]: trace/n = {c_I:.4f}, ||A - (tr/n)·I|| = {np.linalg.norm(A_minus_cI):.4e}")
    if K == 1:
        return [], "fully irreducible: commutant = C·I"
    # For K > 1, find central elements (= elements commuting with each other)
    # Simplification: find common eigenspaces of (any) non-identity commutant element
    # to identify dark subspaces.
    return matrices, f"commutant dim {K} > 1: non-trivial dark subspace structure"


def main():
    out = {}
    for n in (2, 3):
        print(f"\n=== n={n} (V_n = L²((Z/{3**n})*), dim {2*3**(n-1)}) ===")
        t0 = time.time()
        family, coprime = build_kraus_family(n, V_MAX=12, J_MAX=3)
        n_dim = len(coprime)
        print(f"  built {len(family)} distinct Kraus operators in {time.time()-t0:.2f}s")

        null_dim, null_basis, sigma = compute_commutant_dim(family, n_dim)
        print(f"  smallest 5 singular values: {sigma[-5:]}")
        print(f"  largest 5 singular values: {sigma[:5]}")

        matrices, verdict = analyze_commutant(null_basis, n_dim)
        print(f"  VERDICT: {verdict}")

        out[f"n={n}"] = {
            "n": n,
            "n_dim": int(n_dim),
            "n_kraus": len(family),
            "commutant_dim": int(null_dim),
            "smallest_sigmas": [float(s) for s in sigma[-10:]],
            "verdict": verdict,
        }

    with open(os.path.join(OUTDIR, "dark_subspace_probe.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote", os.path.join(OUTDIR, "dark_subspace_probe.json"))


if __name__ == "__main__":
    main()
