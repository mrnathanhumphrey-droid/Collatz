"""
phase4_formula_grid_verification.py — Efficient verification of the corrected R3 formula
across a (n, j) grid using matrix-free ARPACK on the super-operator.

Tests: λ_below(n, j) = 0.5 / |1 − 0.5·e^{iπ/3^{min(n−1, 2j−1)}}|

Grid: (n, j) ∈ {2..6} × {2..4}. Crucial discriminator: at fixed n=5:
  j=2 (saturated, 2j-1=3 < n-1=4): 0.987
  j=3 (unsaturated, 2j-1=5 > n-1=4): 0.9985
If the j=3 case gives 0.9985 (not 0.987), the j-saturation formula is confirmed.

Matrix-free trick: super-op `L(ρ) = Σ_v M_v · ρ · M_v†` applied to vectorized ρ via
direct Kraus action — no dense super-op storage. ARPACK eigs gives top-K efficiently.
"""
import sys, os, json, cmath, time
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigs
sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def build_kraus_family(n, j, b_prior, V_MAX=16):
    N = 3 ** n
    inv2 = pow(2, -1, N)
    coprime = [r for r in range(N) if r % 3 != 0]
    idx_map = {r: i for i, r in enumerate(coprime)}
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
            target = (xi * pow_inv2_v) % N
            if target % 3 == 0:
                continue
            phase = cmath.exp(-1j * two_pi_over_N * xi * x_phase)
            M[idx_map[xi], idx_map[target]] = weight * phase
        family.append(M)
    return family, coprime, idx_map, n_dim


def build_TRUE_DW_basis(coprime, idx_map, n, N):
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
    C = np.zeros((n_fibers, n_dim), dtype=complex)
    for i, fb in enumerate(fibers):
        for xi in fb:
            C[i, idx_map[xi]] = 1.0
    _, S, Vh = np.linalg.svd(C, full_matrices=True)
    rank_C = int(np.sum(S > 1e-9 * max(S.max(), 1.0)))
    B_W = Vh.conj().T[:, rank_C:]
    return B_W


def superop_action(family, B_W):
    """Return a LinearOperator that applies super-op L|_{D_W} to vectorized rho."""
    d_W = B_W.shape[1]
    # Pre-compute M_v in D_W basis
    M_W_list = [B_W.conj().T @ M @ B_W for M in family]
    def matvec(vec_rho):
        # vec_rho is column-major vec of a (d_W × d_W) matrix
        rho = vec_rho.reshape(d_W, d_W, order='F')
        result = np.zeros_like(rho)
        for M_W in M_W_list:
            result += M_W @ rho @ M_W.conj().T
        return result.reshape(-1, order='F')
    return LinearOperator((d_W * d_W, d_W * d_W), matvec=matvec, dtype=complex)


def verify_point(n, j, b_prior=0, V_MAX=16, k_top=10):
    N = 3 ** n
    coprime = [r for r in range(N) if r % 3 != 0]
    idx_map = {r: i for i, r in enumerate(coprime)}
    n_dim = len(coprime)
    family, _, _, _ = build_kraus_family(n, j, b_prior, V_MAX=V_MAX)
    B_W = build_TRUE_DW_basis(coprime, idx_map, n, N)
    d_W = B_W.shape[1]

    # Verify D_W preservation
    P_W = B_W @ B_W.conj().T
    max_leak = max(np.linalg.norm((np.eye(n_dim) - P_W) @ M @ P_W, ord=2) for M in family)

    # Predicted λ_below (CORRECTED: uses min(n-1, 2j-1))
    m_exp = min(n - 1, 2 * j - 1)
    angle = np.pi / 3 ** m_exp
    lam_pred = 0.5 / abs(1 - 0.5 * np.exp(1j * angle))

    L_op = superop_action(family, B_W)
    # ARPACK eigs — use sigma shift near predicted to converge faster
    try:
        eigvals, _ = eigs(L_op, k=k_top, which='LM', maxiter=5000, tol=1e-10)
    except Exception:
        # retry with sigma shift near predicted
        try:
            eigvals, _ = eigs(L_op, k=k_top, sigma=complex(lam_pred, 0), maxiter=5000, tol=1e-8)
        except Exception:
            eigvals = np.array([np.nan])
    eigvals_sorted = sorted(eigvals, key=lambda x: -abs(x))
    n_near_1 = sum(1 for e in eigvals_sorted if abs(abs(e) - 1.0) < 1e-3)
    below_1 = [e for e in eigvals_sorted if abs(e) < 0.999]
    if below_1:
        lam_obs = below_1[0]
        lam_mag = abs(lam_obs)
        lam_arg = float(np.angle(lam_obs))
        period = 2 * np.pi / abs(lam_arg) if abs(lam_arg) > 1e-10 else None
        rel_err = abs(lam_mag - lam_pred) / lam_pred
    else:
        # Need more eigs to get below 1.0
        eigvals, _ = eigs(L_op, k=k_top + 20, which='LM', maxiter=3000, tol=1e-10)
        eigvals_sorted = sorted(eigvals, key=lambda x: -abs(x))
        below_1 = [e for e in eigvals_sorted if abs(e) < 0.999]
        lam_obs = below_1[0] if below_1 else None
        lam_mag = abs(lam_obs) if lam_obs is not None else None
        lam_arg = float(np.angle(lam_obs)) if lam_obs is not None else None
        period = 2 * np.pi / abs(lam_arg) if (lam_arg is not None and abs(lam_arg) > 1e-10) else None
        rel_err = abs(lam_mag - lam_pred) / lam_pred if lam_mag else None

    return {
        "n": n, "j": j, "d_W": int(d_W),
        "max_leak": float(max_leak),
        "lam_predicted": float(lam_pred),
        "lam_observed": float(lam_mag) if lam_mag else None,
        "period": float(period) if period else None,
        "rel_err": float(rel_err) if rel_err else None,
        "n_near_1": int(n_near_1),
    }


def main():
    out = []
    grid = [
        (3, 2), (3, 3),
        (4, 2), (4, 3), (4, 4),
        (5, 2), (5, 3), (5, 4),
    ]
    print(f"{'n':>3} {'j':>3} {'d_W':>5} {'W-leak':>10} {'predicted':>10} {'observed':>10} {'rel_err':>10} {'period':>8} {'n_near_1':>10}")
    for (n, j) in grid:
        t0 = time.time()
        r = verify_point(n, j)
        elapsed = time.time() - t0
        marker = "  OK" if (r["rel_err"] is not None and r["rel_err"] < 0.01) else "  ?"
        print(f"{r['n']:>3} {r['j']:>3} {r['d_W']:>5} {r['max_leak']:>10.2e} {r['lam_predicted']:>10.6f} {r['lam_observed'] or 0:>10.6f} {r['rel_err'] or 0:>10.2e} {r['period'] or 0:>8.2f} {r['n_near_1']:>10}  ({elapsed:.1f}s){marker}")
        out.append({**r, "elapsed_sec": float(elapsed)})

    with open(os.path.join(OUTDIR, "phase4_formula_grid_verification.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote phase4_formula_grid_verification.json")


if __name__ == "__main__":
    main()
