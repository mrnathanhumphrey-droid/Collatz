"""
T_M_truncated_spectrum.py — Numerical spectrum of T_M acting on V_n^M, truncated by
dropping phase-twist contributions that exit V_M (per T_V_DISPOSITION's
H_M_RECURSION_UNDERSPECIFIED finding).

V_n^M is the space of η-indexed pair-correlation moments,
  M_n(η) := Σ_{ξ ∈ (Z/3^n)*} μ̂_n(ξ) · μ̂_n*(ξ · η),  for η in (Z/3^n)*.
dim V_n^M = N_n = 2 * 3^(n-1).

Tao recursion gives:
  M_{n+1}(η') = Σ_{v, v'} (3 · 2^{-v-v'} / Z^2) * 1[η' ≡ 2^{v'-v} (mod 3)] *
                          M_n_with_phase(η · 2^{v-v'}, δ'(v,v',η'))
where δ'(v,v',η') := (1 - η' · 2^{v-v'}) / 3 is generically non-zero,
and M_n_with_phase(η, δ') := Σ_ξ e^{-2πi ξ δ' / 3^n} μ̂_n(ξ) μ̂_n*(ξη).

The δ' ≠ 0 terms generate phase-twisted moments NOT in V_n^M (per T_V_RECURSION
phase obstruction). Truncated T_M sets δ' = 0 in the projection back to V_n^M.

This script:
1. Builds T_M_trunc as a linear operator on V_n^M (square, dim N_n).
   Specifically T_M_trunc: V_{n+1}^M → V_n^M via the recursion + truncation.
   For an endomorphism, take T_M_trunc on V_n^M via composition with the "lift" or via
   self-iteration at level n.
2. The simpler thing to compute first: T_M_trunc as the operator from M_n (level n moments)
   into M_{n+1} (level n+1 moments), then project M_{n+1} back to "level n moments" via
   restriction (η' mod 3^n -> η) keeping only δ' = 0 terms.
3. Compute spectrum at n=2, 3, 4. Look for: 43/45 embedded (consistency check),
   plus other eigenvalues — real, complex CC pair, etc.

Cross-check: T_M_trunc's (1, 4)-eigenvalue on the (P_+, P_-) class-resolved 2x2 subspace
should equal T_lead's 43/45 by construction.
"""
import sys, os, json, time, cmath
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"


def build_T_M_trunc(n, V_MAX=20):
    """Build T_M_trunc: V_n^M -> V_n^M as an N_n x N_n complex matrix.

    Recursion structure:
      M_{n+1}(η') ← Σ_(v,v') with η' ≡ 2^{v'-v} (mod 3) and δ' = (1 - η'·2^{v-v'})/3
                    of 3·2^{-v-v'}/Z² · M_n_with_phase(η·2^{v-v'}, δ')
    Truncation: keep only δ' = 0 terms (i.e., 1 ≡ η'·2^{v-v'} mod 3^{n+1}, NOT just mod 3).
    Then project to level n via η = η' mod 3^n.

    The truncation gives a CLEAN linear operator on V_n^M (no phase-twist exit).
    """
    N_dom = 3 ** n
    N_cod = 3 ** (n + 1)
    coprime_n = [r for r in range(N_dom) if r % 3 != 0]
    coprime_np = [r for r in range(N_cod) if r % 3 != 0]
    idx_n = {r: i for i, r in enumerate(coprime_n)}
    N_n = len(coprime_n)  # dim V_n^M
    inv2_np = pow(2, -1, N_cod)
    Z = 1.0 - 2.0 ** (-V_MAX)

    # For each η' ∈ (Z/3^{n+1})*, compute the M_{n+1}(η') contribution as a linear
    # combination of M_n(η) for various η ∈ (Z/3^n)*.
    # Then reduce η' -> η_red := η' mod 3^n to get back to V_n^M.
    # KEEP only δ' = 0 contributions in this reduction.

    T = np.zeros((N_n, N_n), dtype=complex)  # T[i, j] = coeff of M_n(eta_j) in M_{n+1, trunc}(eta_red_i)
    count_table = np.zeros(N_n, dtype=int)  # how many η' lift each η

    for eta_prime in coprime_np:
        eta_red = eta_prime % N_dom
        if eta_red % 3 == 0:
            continue
        i = idx_n[eta_red]
        count_table[i] += 1
        # Loop (v, v')
        for v in range(1, V_MAX + 1):
            for vp in range(1, V_MAX + 1):
                # Parity check at mod 3:
                # eta' * 2^{v - vp} mod 3 should be 1 (i.e., eta' ≡ 2^{vp - v} mod 3)
                two_vmvp_mod3 = pow(2, v - vp, 3)
                if (eta_prime * two_vmvp_mod3) % 3 != 1:
                    continue
                # δ' check: we want δ = (1 - η' · 2^{v-vp}) ≡ 0 (mod 3^{n+1})
                # NOT just mod 3. Truncated to δ' = 0.
                two_vmvp_full = pow(2, v - vp, N_cod) if v >= vp else pow(inv2_np, vp - v, N_cod)
                if (eta_prime * two_vmvp_full) % N_cod != 1:
                    continue
                # δ' = 0 case: phase trivial.
                # New η for M_n: η_new = (eta_red) · 2^{v - vp} mod 3^n
                two_vmvp_n = pow(2, v - vp, N_dom) if v >= vp else pow(pow(2, -1, N_dom), vp - v, N_dom)
                eta_new = (eta_red * two_vmvp_n) % N_dom
                if eta_new == 0 or eta_new % 3 == 0:
                    continue
                j = idx_n[eta_new]
                weight = 3.0 * (2.0 ** (-v - vp)) / (Z * Z)
                T[i, j] += weight

    # Average over the 3 lifts (count_table should be 3 for each i)
    for i in range(N_n):
        if count_table[i] > 0:
            T[i, :] /= count_table[i]
    return T, coprime_n


def project_class_resolved(T, coprime_n):
    """Project T to the (P_+, P_-) 2x2 class-resolved subspace.
    P_+ = average over η ≡ 1 (mod 3); P_- = average over η ≡ 2 (mod 3)."""
    N_n = T.shape[0]
    plus = [i for i, r in enumerate(coprime_n) if r % 3 == 1]
    minus = [i for i, r in enumerate(coprime_n) if r % 3 == 2]
    # Project: average over each class
    P_plus = np.zeros(N_n, dtype=complex)
    P_minus = np.zeros(N_n, dtype=complex)
    P_plus[plus] = 1.0 / len(plus)
    P_minus[minus] = 1.0 / len(minus)
    # 2x2 matrix in (P_+, P_-) basis
    # M_{++} = (Σ_{i in plus} Σ_{j in plus} T[i,j]) / |plus|
    M22 = np.zeros((2, 2), dtype=complex)
    M22[0, 0] = sum(T[i, j] for i in plus for j in plus) / (len(plus))
    M22[0, 1] = sum(T[i, j] for i in plus for j in minus) / (len(plus))
    M22[1, 0] = sum(T[i, j] for i in minus for j in plus) / (len(minus))
    M22[1, 1] = sum(T[i, j] for i in minus for j in minus) / (len(minus))
    return M22


def run(n, V_MAX=20):
    print(f"\n=== n={n} (V_n^M dim = {2*3**(n-1)}) ===")
    t0 = time.time()
    T, coprime = build_T_M_trunc(n, V_MAX=V_MAX)
    print(f"  built T_M_trunc in {time.time()-t0:.2f}s")

    # Sanity: T should be approximately real (since δ' = 0 truncation kills imaginary phases)
    imag_norm = np.linalg.norm(T.imag)
    real_norm = np.linalg.norm(T.real)
    print(f"  ||Im T|| / ||Re T|| = {imag_norm/max(real_norm,1e-30):.2e}")

    # Compute full spectrum
    eigs = np.linalg.eigvals(T)
    abs_eigs = sorted(np.abs(eigs), reverse=True)
    print(f"  top 10 |eigs|: {[f'{x:.6f}' for x in abs_eigs[:10]]}")

    # Group by complex value
    eigs_by_mod = sorted(eigs, key=lambda z: -abs(z))[:15]
    print(f"  top 15 complex eigs:")
    for i, z in enumerate(eigs_by_mod):
        arg = np.angle(z)
        per = 2 * np.pi / abs(arg) if abs(arg) > 1e-9 else float('inf')
        print(f"    [{i}] {z.real:+.6f}{z.imag:+.6f}j  |z|={abs(z):.6f}  arg={arg:+.4f}  period={per:.3f}")

    # Project to (P_+, P_-) class-resolved 2x2 — should give T_lead's 43/45 if construction is right
    T22 = project_class_resolved(T, coprime)
    eigs22 = np.linalg.eigvals(T22)
    print(f"  (P_+, P_-) class-resolved 2x2:")
    print(f"    matrix = {T22.real}")
    print(f"    eigs = {eigs22}")
    print(f"    expected T_lead's (43/45 = {43/45:.6f}, 0)")

    return {
        "n": n,
        "N_n": len(coprime),
        "imag_norm_over_real": float(imag_norm / max(real_norm, 1e-30)),
        "abs_eigs_top10": [float(x) for x in abs_eigs[:10]],
        "complex_eigs_top15": [[float(z.real), float(z.imag)] for z in eigs_by_mod],
        "class_resolved_eigs": [[float(z.real), float(z.imag)] for z in eigs22],
    }


def main():
    out = {}
    for n in (2, 3, 4):
        out[f"n={n}"] = run(n, V_MAX=20 if n <= 3 else 15)
    with open(os.path.join(OUTDIR, "T_M_truncated_spectrum.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote", os.path.join(OUTDIR, "T_M_truncated_spectrum.json"))


if __name__ == "__main__":
    main()
