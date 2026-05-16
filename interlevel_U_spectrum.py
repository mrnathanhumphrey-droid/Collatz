"""
interlevel_U_spectrum.py — Singular value structure of the Fourier-side Tao transfer operator
  U_n : V_n^Fourier -> V_{n+1}^Fourier

defined by

  (U_n f)(xi') := (1/Z) sum_{v >= 1} 2^{-v} * exp(-2*pi*i * xi' * 2^{-v} / 3^{n+1}) * f(xi' * 2^{-v} mod 3^n)

for xi' in (Z/3^{n+1})*, with f a function on (Z/3^n)*.

By my K_k lemma, P_fiber_sum @ U_n collapses to K_n on V_n (trivial spectrum {1, 0, ...}).
This probe instead looks at the SINGULAR VALUE STRUCTURE of U_n directly:
  - sigma(U_n): singular values, real
  - eigs(U_n^* U_n) on V_n: same as sigma^2, with eigenvectors carrying rate info
  - eigs(U_n U_n^*) on V_{n+1}: ditto, on the bigger space

Plus a "shift-iterate" endomorphism:
  T_shift_n := T_sum @ U_n : V_n -> V_n  (should recover K_n by marginal consistency)
  T_W_n     := P_W @ U_n   : V_n -> W_n  (the W_n-component of inter-level transfer; non-square)

The interesting non-trivial structure (if any) lives in singular values of U_n that are NOT 1 or 0.
If sigma(U_n) is just {1, 0, ..., 0}, the inter-level operator is also trivial and the period-9 CC pair
must live in an infinite-dimensional / continuous-spectrum object.
"""
import sys, os, json, time, cmath
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def build_U(n, V_MAX=30):
    """U_n : V_n -> V_{n+1} (Fourier side, complex).

    Domain  basis: coprime ξ in Z/3^n   (size 2*3^{n-1})
    Codomain basis: coprime ξ' in Z/3^{n+1} (size 2*3^n)
    Matrix entry  U[ξ', ξ] = (1/Z) sum_{v>=1, ξ'*2^{-v} mod 3^n == ξ} 2^{-v} * e^{-2πi ξ' 2^{-v} / 3^{n+1}}
    """
    N_dom = 3 ** n
    N_cod = 3 ** (n + 1)
    coprime_dom = [r for r in range(N_dom) if r % 3 != 0]
    coprime_cod = [r for r in range(N_cod) if r % 3 != 0]
    idx_dom = {r: i for i, r in enumerate(coprime_dom)}
    idx_cod = {r: i for i, r in enumerate(coprime_cod)}
    n_dom = len(coprime_dom)
    n_cod = len(coprime_cod)

    inv2_cod = pow(2, -1, N_cod)
    pow_inv2_cod = [pow(inv2_cod, v, N_cod) for v in range(1, V_MAX + 1)]
    Z = 1.0 - 2.0 ** (-V_MAX)
    two_pi_over_3pp = 2.0 * np.pi / N_cod

    U = np.zeros((n_cod, n_dom), dtype=complex)
    for xi_prime in coprime_cod:
        i = idx_cod[xi_prime]
        for v in range(1, V_MAX + 1):
            xi_prime_inv2v = (xi_prime * pow_inv2_cod[v - 1]) % N_cod  # in Z/3^{n+1}
            # Reduce to level n: xi = xi_prime_inv2v mod 3^n. Must be coprime to 3.
            xi = xi_prime_inv2v % N_dom
            if xi == 0 or xi % 3 == 0:
                continue
            j = idx_dom[xi]
            weight = (2.0 ** (-v)) / Z
            phase = cmath.exp(-1j * two_pi_over_3pp * xi_prime * pow_inv2_cod[v - 1])
            U[i, j] += weight * phase
    return U, coprime_dom, coprime_cod


def build_T_sum(n):
    """T_sum : V_{n+1} -> V_n via fiber-summing.
    T_sum[xi, xi'] = 1 if xi' lifts xi (i.e., xi' mod 3^n = xi), else 0.
    Returns shape (n_dom, n_cod).
    """
    N_dom = 3 ** n
    N_cod = 3 ** (n + 1)
    coprime_dom = [r for r in range(N_dom) if r % 3 != 0]
    coprime_cod = [r for r in range(N_cod) if r % 3 != 0]
    idx_dom = {r: i for i, r in enumerate(coprime_dom)}
    idx_cod = {r: i for i, r in enumerate(coprime_cod)}
    T = np.zeros((len(coprime_dom), len(coprime_cod)))
    for xi_prime in coprime_cod:
        xi = xi_prime % N_dom
        if xi % 3 == 0:
            continue
        T[idx_dom[xi], idx_cod[xi_prime]] = 1.0
    return T


def build_W_basis_at_n_plus_1(n):
    """Orthonormal basis for W_n inside V_{n+1} (3-fiber-zero-mean functions)."""
    N_cod = 3 ** (n + 1)
    N_dom = 3 ** n
    coprime_dom = [r for r in range(N_dom) if r % 3 != 0]
    coprime_cod = [r for r in range(N_cod) if r % 3 != 0]
    idx_cod = {r: i for i, r in enumerate(coprime_cod)}
    cols = []
    for r in coprime_dom:
        lifts = [r, r + N_dom, r + 2 * N_dom]
        idxs = [idx_cod[l] for l in lifts]
        v1 = np.zeros(len(coprime_cod), dtype=complex)
        v1[idxs[0]] = 1.0 / np.sqrt(2.0)
        v1[idxs[1]] = -1.0 / np.sqrt(2.0)
        v2 = np.zeros(len(coprime_cod), dtype=complex)
        v2[idxs[0]] = 1.0 / np.sqrt(6.0)
        v2[idxs[1]] = 1.0 / np.sqrt(6.0)
        v2[idxs[2]] = -2.0 / np.sqrt(6.0)
        cols.append(v1)
        cols.append(v2)
    return np.column_stack(cols)


def run(n, V_MAX=30):
    print(f"\n=== n={n} (U_n : V_{n} ({2*3**(n-1)}-dim) -> V_{n+1} ({2*3**n}-dim)) ===")
    t0 = time.time()
    U, _, _ = build_U(n, V_MAX=V_MAX)
    T_sum = build_T_sum(n)
    print(f"  built U_n and T_sum in {time.time()-t0:.2f}s")

    # 1. Singular values of U_n
    sigmas = np.linalg.svd(U, compute_uv=False)
    print(f"  sigma(U_n) top 8: {sigmas[:8]}")
    print(f"  sigma(U_n) bottom 5: {sigmas[-5:]}")

    # 2. Eigenvalues of U_n^* U_n (square Hermitian, eigenvalues = sigma^2)
    UsU = U.conj().T @ U
    eigs_UsU = np.linalg.eigvalsh(UsU)
    eigs_UsU_sorted = np.sort(eigs_UsU)[::-1]
    print(f"  eig(U_n^* U_n) top 8: {eigs_UsU_sorted[:8]}")

    # 3. T_shift_n := T_sum @ U_n : V_n -> V_n. Should be K_n (trivial spectrum) by marginal consistency.
    T_shift = T_sum @ U
    eigs_T_shift = np.linalg.eigvals(T_shift)
    abs_eigs_T_shift = np.sort(np.abs(eigs_T_shift))[::-1]
    print(f"  eig(T_sum @ U_n) top 8 |.|: {abs_eigs_T_shift[:8]}")
    # Also report the complex eigenvalues themselves
    top_complex = sorted(eigs_T_shift, key=lambda z: -abs(z))[:6]
    print(f"  eig(T_sum @ U_n) top 6 complex: {[f'{z.real:+.4f}{z.imag:+.4f}j' for z in top_complex]}")

    # 4. P_W @ U_n : V_n -> W_n (non-square; singular values)
    W_basis = build_W_basis_at_n_plus_1(n)  # shape (n_cod, dim_W)
    PW_U = W_basis.conj().T @ U  # shape (dim_W, n_dom)
    sigmas_PW_U = np.linalg.svd(PW_U, compute_uv=False)
    print(f"  sigma(P_W @ U_n) top 8: {sigmas_PW_U[:8]}")
    print(f"  sigma(P_W @ U_n) bottom 5: {sigmas_PW_U[-5:]}")

    # 5. ENDOMORPHISM on V_n built from P_W_residue then re-project back via U_n^* P_W
    # T_residue := U_n^* @ P_W @ U_n : V_n -> V_n. Square Hermitian-positive-semidef.
    T_residue = U.conj().T @ (W_basis @ W_basis.conj().T) @ U
    eigs_T_residue = np.linalg.eigvalsh(T_residue)
    eigs_T_residue_sorted = np.sort(eigs_T_residue)[::-1]
    print(f"  eig(U^* P_W U) top 8 (real, hermitian): {eigs_T_residue_sorted[:8]}")

    return {
        "n": n,
        "dim_Vn": int(2 * 3 ** (n - 1)),
        "dim_Vn1": int(2 * 3 ** n),
        "dim_Wn": int(4 * 3 ** (n - 1)),
        "V_MAX": V_MAX,
        "sigma_U_top": [float(x) for x in sigmas[:10]],
        "eig_UsU_top": [float(x) for x in eigs_UsU_sorted[:10]],
        "abs_eig_T_shift_top": [float(x) for x in abs_eigs_T_shift[:10]],
        "complex_eig_T_shift_top": [[float(z.real), float(z.imag)] for z in top_complex],
        "sigma_PW_U_top": [float(x) for x in sigmas_PW_U[:10]],
        "eig_T_residue_top": [float(x) for x in eigs_T_residue_sorted[:10]],
    }


def main():
    out = {}
    for n in (1, 2, 3, 4):
        out[f"n={n}"] = run(n, V_MAX=30 if n <= 3 else 25)
    with open(os.path.join(OUTDIR, "interlevel_U_spectrum.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote", os.path.join(OUTDIR, "interlevel_U_spectrum.json"))

    # Summary table: top singular values across n
    print("\n--- Summary: top 5 sigma(U_n) ---")
    print(f"{'n':>3} {'dim_Vn':>7} {'dim_Vn1':>8} " + "".join(f"{'σ'+str(i+1):>9}" for i in range(5)))
    for n in (1, 2, 3, 4):
        r = out[f"n={n}"]
        s = r["sigma_U_top"]
        print(f"{n:>3} {r['dim_Vn']:>7} {r['dim_Vn1']:>8} " + "".join(f"{x:>9.5f}" for x in s[:5]))

    print("\n--- Summary: top 5 sigma(P_W @ U_n) (W-component of inter-level transfer) ---")
    print(f"{'n':>3} {'dim_Vn':>7} {'dim_Wn':>7} " + "".join(f"{'σ'+str(i+1):>9}" for i in range(5)))
    for n in (1, 2, 3, 4):
        r = out[f"n={n}"]
        s = r["sigma_PW_U_top"]
        print(f"{n:>3} {r['dim_Vn']:>7} {r['dim_Wn']:>7} " + "".join(f"{x:>9.5f}" for x in s[:5]))


if __name__ == "__main__":
    main()
