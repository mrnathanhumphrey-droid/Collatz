"""
phase4_n4_verification.py — Direct verification of Phase 4 closed form at n=4.

R3's Theorem 4.1 predicts: at level n, the first below-commutant eigenvalue of
L|_{D_W} for the j ≥ 2 sub-family is λ_below(n) = 0.5/|1 − 0.5·e^{iπ/3^{n-1}}|.

Predictions:
  n=2 (dim H=6, d_W=4): 1/√3 = 0.577350  ✓ (verified in phase4_dark_spectral_gap_probe.py)
  n=3 (dim H=18, d_W=16): 0.897582 ✓ (verified)
  n=4 (dim H=54, d_W=52): 0.986745 ← NEW VERIFICATION
  n=5 (dim H=162, d_W=160): 0.998499 (predicted, computational scale-up needed)

This probe verifies n=4. Super-operator dim = 52² = 2704; SVD is feasible.
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


def main():
    out = {}
    for n in (4,):
        print(f"\n=== n={n} (H_n = L²((Z/{3**n})*), dim {2*3**(n-1)}) ===")
        N = 3 ** n
        coprime = [r for r in range(N) if r % 3 != 0]
        n_dim = len(coprime)
        B_W, P_W = build_DW_basis(coprime)
        d_W = B_W.shape[1]
        print(f"  d_W (dim D_W) = {d_W}, super-op dim = {d_W*d_W}")

        # Predicted top below-commutant eigenvalue
        angle = np.pi / 3 ** (n - 1)
        lam_below_predicted = 0.5 / abs(1 - 0.5 * np.exp(1j * angle))
        period_predicted = 2 * np.pi / angle
        print(f"  Predicted λ_below(n={n}) = {lam_below_predicted:.6f}, period = {period_predicted:.2f}")

        # Scan over j, b_prior
        results = {}
        for j in (2, 3):
            for b_prior in (0, 1, 2):
                t0 = time.time()
                family, _, _ = build_kraus_family_for_fixed_j_bprior(n, j, b_prior, V_MAX=16)
                max_leak = max(
                    np.linalg.norm((np.eye(n_dim) - P_W) @ op["M"] @ P_W, ord=2)
                    for op in family
                )
                L_DW = channel_superop_restricted(family, B_W)
                t_super = time.time() - t0
                t0 = time.time()
                eigvals = np.linalg.eigvals(L_DW)
                t_eig = time.time() - t0
                eigvals_sorted = sorted(eigvals, key=lambda x: -abs(x))
                # count near-1
                n_near_1 = sum(1 for e in eigvals_sorted if abs(abs(e) - 1.0) < 1e-3)
                below_1 = [e for e in eigvals_sorted if abs(e) < 0.99]
                lambda_below_obs = below_1[0] if below_1 else None
                lam_below_mag = abs(lambda_below_obs) if lambda_below_obs is not None else None
                lam_below_arg = float(np.angle(lambda_below_obs)) if lambda_below_obs is not None else None
                lam_below_period = (2 * np.pi / abs(lam_below_arg)) if (lam_below_arg is not None and abs(lam_below_arg) > 1e-10) else None
                print(f"  j={j}, b_prior={b_prior}: W-leak={max_leak:.2e}, super-op build {t_super:.1f}s, eig {t_eig:.1f}s")
                print(f"    n_near_1 = {n_near_1}, λ_below |·| = {lam_below_mag:.6f}, period = {lam_below_period:.4f}")
                print(f"    match to prediction: {abs(lam_below_mag - lam_below_predicted):.6f} (rel err {abs(lam_below_mag - lam_below_predicted)/lam_below_predicted:.2e})")
                results[f"j={j},b_prior={b_prior}"] = {
                    "j": j, "b_prior": b_prior,
                    "n_near_1": n_near_1,
                    "lambda_below_mag": float(lam_below_mag),
                    "lambda_below_period": float(lam_below_period) if lam_below_period else None,
                    "predicted": float(lam_below_predicted),
                    "rel_err": float(abs(lam_below_mag - lam_below_predicted) / lam_below_predicted),
                    "max_leak": float(max_leak),
                }
        out[f"n={n}"] = {"d_W": int(d_W), "predicted_lambda_below": float(lam_below_predicted), "predicted_period": float(period_predicted), "results": results}

    with open(os.path.join(OUTDIR, "phase4_n4_verification.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote phase4_n4_verification.json")


if __name__ == "__main__":
    main()
