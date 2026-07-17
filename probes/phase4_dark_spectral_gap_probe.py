"""
phase4_dark_spectral_gap_probe.py — Spectral gap of Syracuse's per-step DWM channel
restricted to D_W under the j ≥ 2 sub-family.

Per PHASE2_APPROX_DARK_RESULT.md: D_W (3-fiber-zero-mean, dim n_dim − 2) is EXACTLY
preserved by j ≥ 2 Kraus operators. The asymptotic c=7/45 rate (if a spectral-gap
phenomenon on D_W) should be the second-largest |eigenvalue| of L|_{D_W} where
  L(ρ) = Σ_{v=1}^∞ M_v^{(j, b_prior)} · ρ · (M_v^{(j, b_prior)})†   (fixed j ≥ 2, b_prior)

Build L|_{D_W} as a superoperator matrix (d_W² × d_W²), eigendecompose, extract
{|λ_1|, |λ_2|, |λ_3|, ...}. λ_1 = 1 (trace preservation). λ_2 is the per-step mixing rate.

Test for c=7/45 closure:
  - λ_2 ≈ 0.984 (empirical PADE Hadamard, n=10..13): direct match.
  - λ_2 ≈ 43/45 = 0.95556: T_lead-anchored rate.
  - λ_2 << 0.95: Phase 5 inverse-limit machinery needed.

Scan over (j, b_prior) ∈ {2, 3} × {0, 1, 2} and V_MAX = 16 to check stability.
"""
import sys, os, json, cmath, time
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def build_kraus_family_for_fixed_j_bprior(n, j, b_prior, V_MAX):
    N = 3 ** n
    inv2 = pow(2, -1, N)
    coprime = [r for r in range(N) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(coprime)}
    n_dim = len(coprime)
    two_pi_over_N = 2.0 * np.pi / N
    x_j = pow(3, 2 * j - 2, N) * pow(inv2, b_prior, N) % N

    family = []
    for v in range(1, V_MAX + 1):
        pow_inv2_v = pow(inv2, v, N)
        x_phase = (x_j * pow_inv2_v) % N
        M = np.zeros((n_dim, n_dim), dtype=complex)
        weight = 2.0 ** (-v / 2)
        for xi in coprime:
            target_xi = (xi * pow_inv2_v) % N
            if target_xi % 3 == 0 or target_xi not in idx:
                continue
            phase = cmath.exp(-1j * two_pi_over_N * xi * x_phase)
            M[idx[xi], idx[target_xi]] = weight * phase
        family.append({"v": v, "M": M, "x_phase": x_phase})
    return family, coprime, n_dim


def build_DW_basis(coprime):
    n_dim = len(coprime)
    class_plus = [i for i, xi in enumerate(coprime) if xi % 3 == 1]
    class_minus = [i for i, xi in enumerate(coprime) if xi % 3 == 2]
    e_plus = np.zeros(n_dim, dtype=complex)
    e_plus[class_plus] = 1.0 / np.sqrt(len(class_plus))
    e_minus = np.zeros(n_dim, dtype=complex)
    e_minus[class_minus] = 1.0 / np.sqrt(len(class_minus))
    P_class = np.outer(e_plus, e_plus.conj()) + np.outer(e_minus, e_minus.conj())
    P_W = np.eye(n_dim, dtype=complex) - P_class
    eigvals, eigvecs = np.linalg.eigh(P_W)
    mask = eigvals > 0.5
    B_W = eigvecs[:, mask]
    return B_W, P_W


def channel_superop_restricted(family, B_W):
    d_W = B_W.shape[1]
    super_op = np.zeros((d_W * d_W, d_W * d_W), dtype=complex)
    for op in family:
        M = op["M"]
        M_W = B_W.conj().T @ M @ B_W
        super_op += np.kron(M_W.conj(), M_W)
    return super_op


def channel_superop_full(family, n_dim):
    super_op = np.zeros((n_dim * n_dim, n_dim * n_dim), dtype=complex)
    for op in family:
        M = op["M"]
        super_op += np.kron(M.conj(), M)
    return super_op


def main():
    out = {}
    for n in (2, 3):
        print(f"\n=== n={n} (H_n = L²((Z/{3**n})*), dim {2*3**(n-1)}) ===")
        N = 3 ** n
        coprime = [r for r in range(N) if r % 3 != 0]
        n_dim = len(coprime)
        B_W, P_W = build_DW_basis(coprime)
        d_W = B_W.shape[1]
        print(f"  d_W (= dim D_W) = {d_W}")

        results_at_n = {}
        for j in (2, 3):
            for b_prior in (0, 1, 2):
                family, _, _ = build_kraus_family_for_fixed_j_bprior(n, j, b_prior, V_MAX=16)
                # Verify family preserves D_W exactly
                max_leak = max(
                    np.linalg.norm((np.eye(n_dim) - P_W) @ op["M"] @ P_W, ord=2)
                    for op in family
                )
                # POVM resolution check
                S = sum(op["M"].conj().T @ op["M"] for op in family)
                povm_dev = np.linalg.norm(S - np.eye(n_dim), ord=2)

                L_DW = channel_superop_restricted(family, B_W)
                eigvals = np.linalg.eigvals(L_DW)
                eigvals_sorted = sorted(eigvals, key=lambda x: -abs(x))
                d_W = B_W.shape[1]
                n_print = min(20, d_W * d_W)
                top = eigvals_sorted[:n_print]
                # Count eigenvalues near 1.0 (within 1e-3 of unity in magnitude)
                n_near_1 = sum(1 for e in eigvals_sorted if abs(abs(e) - 1.0) < 1e-3)
                # Largest |eigenvalue| strictly below the 1.0-cluster (~0.99 cutoff)
                below_1 = [abs(e) for e in eigvals_sorted if abs(e) < 0.99]
                lambda_first_below = below_1[0] if below_1 else 0.0

                print(f"  j={j}, b_prior={b_prior}: W-leak={max_leak:.2e}, POVM-dev={povm_dev:.2e}")
                print(f"    n_near_1 = {n_near_1} (commutant-dim on D_W) ; first |λ| below cluster = {lambda_first_below:.6f}")
                print(f"    top {n_print} |λ|: " + ", ".join(f"{abs(e):.4f}" for e in top))
                print(f"    compare: 43/45={43/45:.6f}, empirical 0.984, 1/√3={1/np.sqrt(3):.6f}")

                key = f"j={j},b_prior={b_prior}"
                results_at_n[key] = {
                    "j": j, "b_prior": b_prior,
                    "top20_abs": [float(abs(e)) for e in top],
                    "top20_complex": [[float(e.real), float(e.imag)] for e in top],
                    "n_near_1": int(n_near_1),
                    "lambda_first_below_cluster": float(lambda_first_below),
                    "max_leak": float(max_leak),
                    "povm_dev": float(povm_dev),
                }
        out[f"n={n}"] = {"d_W": int(d_W), "results": results_at_n}

    print("\n--- SUMMARY (commutant on D_W + first |λ| below cluster) ---")
    for nk, d in out.items():
        print(f"{nk}: d_W={d['d_W']}")
        for key, r in d["results"].items():
            print(f"  {key}: n_near_1 = {r['n_near_1']:2d} ; first |λ| below = {r['lambda_first_below_cluster']:.6f}")

    with open(os.path.join(OUTDIR, "phase4_dark_spectral_gap_probe.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote", os.path.join(OUTDIR, "phase4_dark_spectral_gap_probe.json"))


if __name__ == "__main__":
    main()
