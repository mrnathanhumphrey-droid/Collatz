"""
phase4_n4_true_DW_verification.py — Re-verify Phase 4 closed form at n=4
using the TRUE D_W = 3-fiber-zero-mean (dim 36), NOT class^⊥ (dim 52).

R3 §3 claims D_W = 3-fiber-zero-mean is EXACTLY dark under j ≥ 2. The Phase 2
probe used class^⊥ as a proxy for D_W. At n=2, 3 these coincide. At n=4 they
DIVERGE: class^⊥ has dim 52, but true D_W has dim 36.

This probe:
1. Builds the TRUE D_W basis at n=4 via null-space of fiber-sum constraints.
2. Verifies j ≥ 2 family preserves true D_W exactly (machine epsilon leakage).
3. Computes L|_{D_W} spectrum on the d_W = 36 subspace.
4. Compares first below-commutant eigenvalue to predicted 0.5/|1−0.5·e^{iπ/27}| = 0.987.
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
    return family, coprime, idx, n_dim


def build_TRUE_DW_basis(coprime, idx, n, N):
    """True D_W = 3-fiber-zero-mean subspace.
    Conditions: for each ξ_0 in (Z/3^{n-1})* coprime to 3, Σ_{a=0,1,2} f(ξ_0 + a·3^{n-1}) = 0.
    Number of conditions = number of distinct 3-fibers = 2·3^{n-2}.
    """
    n_dim = len(coprime)
    fiber_step = 3 ** (n - 1)
    fibers = []
    seen = set()
    for r in coprime:
        fiber = tuple(sorted([(r + a * fiber_step) % N for a in range(3)]))
        if fiber not in seen and all(f % 3 != 0 for f in fiber):
            seen.add(fiber)
            fibers.append(fiber)
    n_fibers = len(fibers)
    print(f"  3-fibers: {n_fibers} (expected 2·3^{n-2} = {2*3**(n-2)})")
    # Build constraint matrix C: each row = indicator of fiber
    C = np.zeros((n_fibers, n_dim), dtype=complex)
    for i, fiber in enumerate(fibers):
        for xi in fiber:
            C[i, idx[xi]] = 1.0
    # D_W = orthogonal complement of row span of C
    # Use SVD: D_W is span of right singular vectors with singular value 0
    U, S, Vh = np.linalg.svd(C, full_matrices=True)
    rank_C = int(np.sum(S > 1e-9 * max(S.max(), 1.0)))
    print(f"  rank(C) = {rank_C} (expected {n_fibers}, fibers linearly independent)")
    B_W = Vh.conj().T[:, rank_C:]
    print(f"  d_W (TRUE D_W dim) = {B_W.shape[1]}")
    # Build projection P_W = B_W · B_W^†
    P_W = B_W @ B_W.conj().T
    return B_W, P_W, fibers


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
    for n in (3, 4):
        print(f"\n=== n={n} (H_n = L²((Z/{3**n})*), dim {2*3**(n-1)}) ===")
        N = 3 ** n
        coprime = [r for r in range(N) if r % 3 != 0]
        idx_map = {r: i for i, r in enumerate(coprime)}
        n_dim = len(coprime)
        print(f"  n_dim = {n_dim}, TRUE D_W = 3-fiber-zero-mean")
        B_W, P_W, fibers = build_TRUE_DW_basis(coprime, idx_map, n, N)
        d_W = B_W.shape[1]

        # Predicted top below-commutant eigenvalue (from R3 formula)
        angle = np.pi / 3 ** (n - 1)
        lam_below_predicted = 0.5 / abs(1 - 0.5 * np.exp(1j * angle))
        period_predicted = 2 * np.pi / angle
        print(f"  Predicted λ_below(n={n}) = {lam_below_predicted:.6f}, period = {period_predicted:.2f}")

        # Verify j ≥ 2 family preserves TRUE D_W
        results = {}
        for j in (2, 3):
            for b_prior in (0,):
                t0 = time.time()
                family, _, _, _ = build_kraus_family_for_fixed_j_bprior(n, j, b_prior, V_MAX=16)
                max_leak = max(
                    np.linalg.norm((np.eye(n_dim) - P_W) @ op["M"] @ P_W, ord=2)
                    for op in family
                )
                L_DW = channel_superop_restricted(family, B_W)
                eigvals = np.linalg.eigvals(L_DW)
                eigvals_sorted = sorted(eigvals, key=lambda x: -abs(x))
                n_near_1 = sum(1 for e in eigvals_sorted if abs(abs(e) - 1.0) < 1e-3)
                below_1 = [e for e in eigvals_sorted if abs(e) < 0.99]
                lambda_below_obs = below_1[0] if below_1 else None
                if lambda_below_obs is not None:
                    lam_below_mag = abs(lambda_below_obs)
                    lam_below_arg = float(np.angle(lambda_below_obs))
                    lam_below_period = (2 * np.pi / abs(lam_below_arg)) if abs(lam_below_arg) > 1e-10 else None
                else:
                    lam_below_mag = None
                    lam_below_period = None
                t1 = time.time() - t0
                print(f"  j={j}, b_prior={b_prior}: TRUE D_W W-leak = {max_leak:.3e}, n_near_1 = {n_near_1}")
                if lam_below_mag is not None:
                    print(f"    λ_below |·| = {lam_below_mag:.6f} (predicted {lam_below_predicted:.6f})")
                    print(f"    rel err = {abs(lam_below_mag - lam_below_predicted)/lam_below_predicted:.3e}")
                    print(f"    period = {lam_below_period}")
                results[f"j={j},b_prior={b_prior}"] = {
                    "max_leak": float(max_leak),
                    "n_near_1": n_near_1,
                    "lambda_below_mag": float(lam_below_mag) if lam_below_mag else None,
                    "lambda_below_period": float(lam_below_period) if lam_below_period else None,
                    "predicted": float(lam_below_predicted),
                    "rel_err": float(abs(lam_below_mag - lam_below_predicted) / lam_below_predicted) if lam_below_mag else None,
                }
        out[f"n={n}"] = {
            "d_W_true": int(d_W),
            "predicted_lambda_below": float(lam_below_predicted),
            "results": results,
        }

    with open(os.path.join(OUTDIR, "phase4_n4_true_DW_verification.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote phase4_n4_true_DW_verification.json")


if __name__ == "__main__":
    main()
