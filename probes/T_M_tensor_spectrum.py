"""
T_M_tensor_spectrum.py — Bilinear pair-correlation operator U_n ⊗ conj(U_n) acting on
V_n ⊗ V_n* (dim N_n², where N_n = 2*3^(n-1)).

This is the natural setting for T_M without V_M closure issues. The pair-correlation
M_n(η) = Σ_ξ μ̂_n(ξ) μ̂_n*(ξη) corresponds to a SPECIFIC element of V_n ⊗ V_n*; the
recursion μ̂_n → μ̂_{n+1} via U_n lifts to (μ̂_n ⊗ μ̂_n*) → U_n μ̂_n ⊗ conj(U_n) μ̂_n*
which is U_n ⊗ conj(U_n) at the tensor level.

The "natural endomorphism on V_n ⊗ V_n*" via inter-level transfer is:
   Phi_M := (P_back ⊗ P_back) ∘ (U_n ⊗ conj(U_n)) : V_n ⊗ V_n* → V_n ⊗ V_n*
where P_back: V_{n+1} → V_n is some "level reduction" operator.

But P_fiber_sum @ U_n = 0 (cube-root cancellation), so trivial fiber-sum gives zero.
Need twisted projection. The natural choice for pair-correlation is the pair (T^omega @ U_n)
on both factors.

PROBE: compute spectrum of (T^ω1 ⊗ T^{ω2}) ∘ (U_n ⊗ conj(U_n)) for ω1, ω2 ∈ {ω_3, ω_3²}.
Looking for:
1. The 43/45 eigenvalue (T_lead reproduction) somewhere
2. CC pair eigenvalues matching empirical period-9.2
3. Sub-leading structure

Note: V_n ⊗ V_n* has dim N_n² which grows fast:
  n=2: 36, n=3: 324, n=4: 2916.
At n=3, full eigenvalue computation of 324x324 complex matrix is fast (~seconds).
At n=4, 2916² ≈ 8.5M entries; eigenvalue computation ~minute.
"""
import sys, os, json, time, cmath
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"


def build_U(n, V_MAX=30):
    """U_n : V_n -> V_{n+1} (Fourier side, complex)."""
    N_dom = 3 ** n
    N_cod = 3 ** (n + 1)
    coprime_dom = [r for r in range(N_dom) if r % 3 != 0]
    coprime_cod = [r for r in range(N_cod) if r % 3 != 0]
    idx_dom = {r: i for i, r in enumerate(coprime_dom)}
    idx_cod = {r: i for i, r in enumerate(coprime_cod)}
    inv2_cod = pow(2, -1, N_cod)
    pow_inv2_cod = [pow(inv2_cod, v, N_cod) for v in range(1, V_MAX + 1)]
    Z = 1.0 - 2.0 ** (-V_MAX)
    two_pi_over_3pp = 2.0 * np.pi / N_cod
    U = np.zeros((len(coprime_cod), len(coprime_dom)), dtype=complex)
    for xi_prime in coprime_cod:
        i = idx_cod[xi_prime]
        for v in range(1, V_MAX + 1):
            xi_prime_inv2v = (xi_prime * pow_inv2_cod[v - 1]) % N_cod
            xi = xi_prime_inv2v % N_dom
            if xi == 0 or xi % 3 == 0:
                continue
            j = idx_dom[xi]
            U[i, j] += (2.0 ** (-v)) / Z * cmath.exp(-1j * two_pi_over_3pp * xi_prime * pow_inv2_cod[v - 1])
    return U


def build_T_omega(n, omega):
    """T^omega : V_{n+1} -> V_n, twisted fiber selection."""
    N_dom = 3 ** n
    N_cod = 3 ** (n + 1)
    coprime_dom = [r for r in range(N_dom) if r % 3 != 0]
    coprime_cod = [r for r in range(N_cod) if r % 3 != 0]
    idx_dom = {r: i for i, r in enumerate(coprime_dom)}
    idx_cod = {r: i for i, r in enumerate(coprime_cod)}
    T = np.zeros((len(coprime_dom), len(coprime_cod)), dtype=complex)
    for xi in coprime_dom:
        for a in (0, 1, 2):
            xi_prime = xi + a * N_dom
            if xi_prime % 3 == 0:
                continue
            T[idx_dom[xi], idx_cod[xi_prime]] = (1.0 / 3.0) * (omega ** a)
    return T


def run(n, V_MAX=30):
    print(f"\n=== n={n} (V_n^bilin dim = {(2*3**(n-1))**2}) ===")
    t0 = time.time()
    U = build_U(n, V_MAX=V_MAX)
    omega3 = cmath.exp(2j * np.pi / 3)
    T1 = build_T_omega(n, 1.0)
    T_w = build_T_omega(n, omega3)
    T_w2 = build_T_omega(n, omega3 ** 2)
    print(f"  built U, T^ω, T^ω² in {time.time()-t0:.2f}s")

    # Phi_omega = T^omega @ U_n on V_n (we already have this)
    Phi_w = T_w @ U
    Phi_w2 = T_w2 @ U

    # Bilinear tensor operators on V_n ⊗ V_n via (T^ω1, T^ω2_conjugate) projection:
    # The pair-correlation Tao recursion is U_n ⊗ conj(U_n) at the tensor level.
    # For the back-projection, use T^omega ⊗ conj(T^omega) at the dual:
    # Phi_M^{(ω,ω')} := (T^ω ⊗ T^ω'_conj) ∘ (U_n ⊗ conj(U_n)) = (T^ω @ U_n) ⊗ conj(T^ω' @ U_n)
    #
    # Note conj(T^ω) = T^conj(ω). For ω = ω_3, conj = ω_3². For ω = ω_3², conj = ω_3.
    #
    # By the tensor structure:
    # spec(Phi_M^{ω,ω'}) = { μ · conj(ν) : μ ∈ spec(Phi_ω), ν ∈ spec(Phi_ω') }

    # Compute spec of each Phi_ω first
    eigs_Phi_w = np.linalg.eigvals(Phi_w)
    eigs_Phi_w2 = np.linalg.eigvals(Phi_w2)

    # Top moduli
    abs_w = sorted([abs(z) for z in eigs_Phi_w], reverse=True)
    abs_w2 = sorted([abs(z) for z in eigs_Phi_w2], reverse=True)
    print(f"  spec(Phi_ω) top |z|: {abs_w[:5]}")
    print(f"  spec(Phi_ω²) top |z|: {abs_w2[:5]}")

    # Tensor products: Phi_ω ⊗ conj(Phi_ω) means eigenvalues are mu_i * conj(mu_j).
    # But conjugating eigenvalues: conj(Phi_omega's eigvals). Then tensor with Phi_omega:
    # eigvals = {mu_i * conj(mu_j)}.
    # We want the "natural" pair-correlation operator on V_n ⊗ V_n* under inter-level transfer.
    #
    # The pair-correlation Tao recursion (M_n -> M_{n+1}) lifts the structure
    # μ̂_n(ξ) μ̂_n*(ξ') at the level-n bilinear, applies U_n on first factor and
    # conj(U_n) on second factor, gets level-(n+1) bilinear.
    # Project back via fiber-summing on both factors with twists ω, ω' (with ω·ω' = 1 to
    # preserve the "trace" / autocorrelation structure where η = ξ' ξ^{-1} carries through):
    # Effective projection: pair (T^ω, T^{ω^-1}) for each ω.
    # This corresponds to spec(Phi_omega) ⊗ spec(conj(Phi_{omega})) = spec(Phi_omega ⊗ conj(Phi_omega))
    # whose eigenvalues are { mu_i * conj(mu_j) : mu_i, mu_j ∈ spec(Phi_omega) }.

    print(f"\n  Tensor product spectra (M_n ⊗ conj_M_n bilinears):")
    for label, eigs in [("Phi_ω ⊗ conj(Phi_ω)", eigs_Phi_w), ("Phi_ω² ⊗ conj(Phi_ω²)", eigs_Phi_w2)]:
        # tensor eigs = mu_i * conj(mu_j)
        tensor_eigs = np.array([m * np.conj(n_e) for m in eigs for n_e in eigs])
        # Drop those with |z| < 0.1 (small/numerical)
        big = tensor_eigs[np.abs(tensor_eigs) > 0.1]
        big_sorted = sorted(big, key=lambda z: -abs(z))
        print(f"    {label}: top 10 (by |z|)")
        for i, z in enumerate(big_sorted[:10]):
            mod = abs(z)
            arg = np.angle(z)
            per = 2 * np.pi / abs(arg) if abs(arg) > 1e-9 else float('inf')
            print(f"      [{i}] {z.real:+.6f}{z.imag:+.6f}j  |z|={mod:.6f}  arg={arg:+.4f}  period={per:.3f}")

    # Mixed tensor: Phi_omega ⊗ conj(Phi_omega²)
    # eigs = { mu_i (in Phi_omega) * conj(nu_j) (in Phi_omega^2) }
    mixed_eigs = np.array([m * np.conj(n_e) for m in eigs_Phi_w for n_e in eigs_Phi_w2])
    big_mixed = mixed_eigs[np.abs(mixed_eigs) > 0.1]
    big_mixed_sorted = sorted(big_mixed, key=lambda z: -abs(z))
    print(f"\n  Mixed tensor Phi_ω ⊗ conj(Phi_ω²): top 10 (by |z|)")
    for i, z in enumerate(big_mixed_sorted[:10]):
        mod = abs(z)
        arg = np.angle(z)
        per = 2 * np.pi / abs(arg) if abs(arg) > 1e-9 else float('inf')
        print(f"    [{i}] {z.real:+.6f}{z.imag:+.6f}j  |z|={mod:.6f}  arg={arg:+.4f}  period={per:.3f}")

    # 43/45 check: does 43/45 ≈ 0.956 appear anywhere?
    target = 43.0 / 45.0
    print(f"\n  Target 43/45 = {target:.6f}; closest eigenvalue across all bilinear tensors:")
    all_tensor = np.concatenate([
        np.array([m * np.conj(n_e) for m in eigs_Phi_w for n_e in eigs_Phi_w]),
        np.array([m * np.conj(n_e) for m in eigs_Phi_w2 for n_e in eigs_Phi_w2]),
        np.array([m * np.conj(n_e) for m in eigs_Phi_w for n_e in eigs_Phi_w2]),
    ])
    closest = min(all_tensor, key=lambda z: abs(abs(z) - target))
    print(f"    closest: {closest.real:+.6f}{closest.imag:+.6f}j  |z|={abs(closest):.6f}  arg={np.angle(closest):+.4f}")

    return {
        "n": n,
        "spec_Phi_omega_abs_top5": abs_w[:5],
        "spec_Phi_omega2_abs_top5": abs_w2[:5],
        "tensor_w_w_top10": [[float(z.real), float(z.imag)] for z in sorted([m*np.conj(n_e) for m in eigs_Phi_w for n_e in eigs_Phi_w], key=lambda z: -abs(z))[:10]],
        "tensor_w2_w2_top10": [[float(z.real), float(z.imag)] for z in sorted([m*np.conj(n_e) for m in eigs_Phi_w2 for n_e in eigs_Phi_w2], key=lambda z: -abs(z))[:10]],
        "tensor_w_w2_mixed_top10": [[float(z.real), float(z.imag)] for z in big_mixed_sorted[:10]],
        "closest_to_43_45": [float(closest.real), float(closest.imag), float(abs(closest))],
    }


def main():
    out = {}
    for n in (2, 3, 4, 5):
        out[f"n={n}"] = run(n, V_MAX=30 if n <= 3 else 25)
    with open(os.path.join(OUTDIR, "T_M_tensor_spectrum.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote", os.path.join(OUTDIR, "T_M_tensor_spectrum.json"))


if __name__ == "__main__":
    main()
