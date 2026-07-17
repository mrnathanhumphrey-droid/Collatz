"""
result_renormalization_spectrum.py
==================================
Two complementary extractions of the slow-mode spectrum governing ε_k:

  Part 1 — Linear-recurrence fit on the ε_k scalar sequence.
    Fit ε_{k+1} = α₁ ε_k + α₂ ε_{k-1} + ... at orders r = 1..4.
    Roots of the characteristic polynomial are the effective eigenvalues
    of the renormalization in the ε projection.
    Predicts ε_12, ε_13.

  Part 2 — Inter-level renormalization operator R_k at k = 3..6.
    R_k acts on the wavelet residual subspace W_k = V_{k+1} ⊖ T(V_k),
    where T is the level-k → k+1 measure-lift.
    R_k = P_W K_{k+1} P_W^T   (block of K_{k+1} in W_k orthonormal basis).
    Eigenvalues of R_k at small k are the spectrum of the structural
    renormalization. If a complex pair of magnitude ~constant emerges
    across k=3..6, the slow oscillating mode is structurally identified.
"""
from __future__ import annotations

import csv
import os
import sys
import time
from fractions import Fraction

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz"
ENV_CSV = os.path.join(OUT_DIR, "result_q_sweep_test_1_envelope.csv")
OUT_MD = os.path.join(OUT_DIR, "result_renormalization_spectrum.md")
OUT_RECUR_CSV = os.path.join(OUT_DIR,
                              "result_renormalization_recurrence_fits.csv")
OUT_RSPEC_CSV = os.path.join(OUT_DIR,
                              "result_renormalization_R_spectrum.csv")

# Cached eps values from prior verified runs
EPS_CACHED = {
    6: -4.9790566522e-04,
    7: -1.1752368304e-03,
    8: -7.4554636729e-04,
    9: -7.5202571564e-06,
    10: +7.2075091711e-04,
    11: +1.5019670121e-03,
}


def order_of_two(N):
    assert N % 2 == 1
    m = 1
    v = 2 % N
    while v != 1:
        v = (v * 2) % N
        m += 1
    return m


def build_K_dense(q, k):
    N = q ** k
    M = order_of_two(N)
    inv2 = pow(2, -1, N)
    powers_inv2 = np.empty(M, dtype=np.int64)
    pi = inv2
    for v in range(M):
        powers_inv2[v] = pi
        pi = (pi * inv2) % N
    coprime = np.array([r for r in range(N) if r % q != 0], dtype=np.int64)
    n = len(coprime)
    state_idx = -np.ones(N, dtype=np.int64)
    for i, r in enumerate(coprime):
        state_idx[r] = i
    Z_v = 1.0 - 2.0 ** (-M)
    weights = (2.0 ** (-np.arange(1, M + 1))) / Z_v
    K = np.zeros((n, n), dtype=np.float64)
    for i_r in range(n):
        r = int(coprime[i_r])
        base = (q * r + 1) % N
        targets = (base * powers_inv2) % N
        js = state_idx[targets]
        K[i_r] = np.bincount(js, weights=weights, minlength=n)
    return K, coprime, state_idx


def build_T_lift(q, k_lower):
    """T: V_k → V_{k+1} as a sparse-ish matrix S with shape (n_k, n_{k+1}).
       u @ S gives the natural mass-preserving lift of u from V_k to V_{k+1}.
       S[r_idx, r'_idx] = 1/q if r' ≡ r mod q^k_lower (and both coprime to q).
       Each row of S has exactly q nonzero entries summing to 1.
    """
    k = k_lower
    N_k = q ** k
    N_kp1 = q ** (k + 1)
    coprime_k = [r for r in range(N_k) if r % q != 0]
    coprime_kp1 = [r for r in range(N_kp1) if r % q != 0]
    n_k = len(coprime_k)
    n_kp1 = len(coprime_kp1)
    state_idx_k = {r: i for i, r in enumerate(coprime_k)}
    state_idx_kp1 = {r: i for i, r in enumerate(coprime_kp1)}
    S = np.zeros((n_k, n_kp1), dtype=np.float64)
    inv_q = 1.0 / q
    for r_p in coprime_kp1:
        i_p = state_idx_kp1[r_p]
        r = r_p % N_k
        i = state_idx_k[r]
        S[i, i_p] = inv_q
    return S, coprime_k, coprime_kp1


def build_W_basis(S):
    """Given S of shape (n_k, n_{k+1}), return orthonormal basis P_W
    (shape (n_{k+1} - n_k, n_{k+1})) for W_k = orthogonal complement of
    row space of S in R^{n_{k+1}}."""
    Q, _ = np.linalg.qr(S.T, mode="complete")
    n_k = S.shape[0]
    P_W = Q[:, n_k:].T  # rows are orthonormal basis for W_k
    return P_W


def project_R(K_kp1, P_W):
    """R_k = P_W @ K_{k+1} @ P_W^T."""
    return P_W @ K_kp1 @ P_W.T


def stationary(K, tol=1e-14, max_iter=10000):
    n = K.shape[0]
    pi = np.full(n, 1.0 / n)
    for it in range(max_iter):
        pi_new = pi @ K
        pi_new /= pi_new.sum()
        r = float(np.linalg.norm(pi_new - pi, ord=1))
        pi = pi_new
        if r < tol:
            return pi, it + 1, r
    return pi, max_iter, r


def fit_recurrence(eps_values, ks, order):
    """Fit eps_{k+order} = sum_{j=0}^{order-1} alpha_j eps_{k+order-1-j}
    by ordinary least squares on consecutive (k, eps_k) pairs.

    Returns: alpha_vec (length order), residuals, predicted vs actual table.
    """
    y = np.array([eps_values[k] for k in ks], dtype=np.float64)
    n_y = len(y)
    if n_y < order + 1:
        return None
    n_eq = n_y - order
    A = np.zeros((n_eq, order))
    b = np.zeros(n_eq)
    for i in range(n_eq):
        b[i] = y[order + i]
        for j in range(order):
            A[i, j] = y[order + i - 1 - j]
    coeffs, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    pred = A @ coeffs
    resid = b - pred
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((b - b.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    # Characteristic polynomial: z^order - α_1 z^(order-1) - α_2 z^(order-2) - ...
    poly_coefs = np.array([1.0] + list(-coeffs))
    roots = np.roots(poly_coefs)
    return {
        "order": order,
        "alphas": coeffs,
        "roots": roots,
        "ss_res": ss_res,
        "r2": r2,
        "predicted": pred,
        "actual": b,
        "n_eq": n_eq,
    }


def predict_forward(eps_values, ks, alphas, order, k_targets):
    """Roll the recurrence forward. eps_values is dict, alphas is the fitted
    coefficients in the order eps_{k+order} = α_1 eps_{k+order-1} + ..."""
    eps_extended = dict(eps_values)
    for k_target in k_targets:
        # eps_{k_target} = sum_{j=0}^{order-1} alphas[j] * eps_{k_target-1-j}
        val = 0.0
        for j in range(order):
            val += alphas[j] * eps_extended[k_target - 1 - j]
        eps_extended[k_target] = val
    return eps_extended


def main():
    print("=" * 78)
    print("RENORMALIZATION SPECTRUM — recurrence + R operator")
    print("=" * 78)

    # Load eps values
    S_cache = {}
    with open(ENV_CSV) as f:
        for row in csv.DictReader(f):
            if int(row["q"]) != 3:
                continue
            S_cache[int(row["n"])] = Fraction(int(row["S_n_num"]),
                                                int(row["S_n_den"]))
    eps_all = {k: float(S_cache[k] - Fraction(7, 15))
               for k in S_cache}
    eps_all.update(EPS_CACHED)
    ks_sorted = sorted(eps_all.keys())
    print(f"\nLoaded eps_k for k = {ks_sorted}")
    for k in ks_sorted:
        print(f"  eps_{k} = {eps_all[k]:+.10e}")

    # ==============================================================
    # PART 1 — Linear recurrence fit
    # ==============================================================
    print(f"\n{'=' * 78}")
    print("PART 1 — Linear-recurrence fit on epsilon_k sequence")
    print(f"{'=' * 78}")

    fits = {}
    for order in [1, 2, 3, 4]:
        fit = fit_recurrence(eps_all, ks_sorted, order)
        if fit is None:
            continue
        fits[order] = fit
        print(f"\nOrder {order}: {fit['n_eq']} equations on {len(ks_sorted)} points")
        print(f"  alphas = {fit['alphas']}")
        print(f"  R^2 = {fit['r2']:.6f}, ss_res = {fit['ss_res']:.4e}")
        print(f"  characteristic roots:")
        for j, root in enumerate(fit["roots"]):
            mag = abs(root)
            arg = np.angle(root)
            print(f"    z_{j+1} = {root.real:+.4f} {'+' if root.imag >= 0 else '-'} "
                  f"{abs(root.imag):.4f}i  "
                  f"|z| = {mag:.4f}, arg = {arg:.4f} rad "
                  f"({np.degrees(arg):+.1f} deg)")

    # Predict forward using the highest-order well-conditioned fit
    # Choose the fit with best R^2 that's not over-fit (use BIC-like check)
    print(f"\n--- Forward predictions for eps_12, eps_13, eps_14, eps_15 ---")
    pred_table = {}
    for order in sorted(fits.keys()):
        fit = fits[order]
        eps_ext = predict_forward(eps_all, ks_sorted, fit["alphas"],
                                    order, [12, 13, 14, 15])
        pred_table[order] = eps_ext
        print(f"\nOrder {order} predictions:")
        for k_p in [12, 13, 14, 15]:
            print(f"  eps_{k_p} = {eps_ext[k_p]:+.4e}")

    # ==============================================================
    # PART 2 — R operator spectrum at k=3..6
    # ==============================================================
    print(f"\n{'=' * 78}")
    print("PART 2 — Renormalization operator R_k spectrum at k = 3..6")
    print(f"{'=' * 78}")
    print("\nR_k acts on W_k = V_{k+1} ⊖ T(V_k), where T is the measure-lift.")
    print("R_k = P_W K_{k+1} P_W^T   (orthonormal projection)\n")

    K_cache = {}  # k → K_k
    R_results = {}  # k → eigenvalues of R_k
    for k_lower in [3, 4, 5, 6]:
        k_upper = k_lower + 1
        print(f"--- R_{k_lower} (acts on W_{k_lower}, K_{k_upper} restricted) ---")
        t0 = time.time()
        if k_upper not in K_cache:
            K_kp1, _, _ = build_K_dense(3, k_upper)
            K_cache[k_upper] = K_kp1
        else:
            K_kp1 = K_cache[k_upper]
        n_kp1 = K_kp1.shape[0]
        S, _, _ = build_T_lift(3, k_lower)
        n_k = S.shape[0]
        P_W = build_W_basis(S)
        dim_W = P_W.shape[0]
        t_setup = time.time() - t0
        print(f"  K_{k_upper}: {n_kp1} states, T_lift S: ({n_k}, {n_kp1}), "
              f"dim(W_{k_lower}) = {dim_W}, setup {t_setup:.2f}s")

        t0 = time.time()
        R_k = project_R(K_kp1, P_W)
        eigs = np.linalg.eigvals(R_k)
        t_eig = time.time() - t0

        # Sort by magnitude descending
        order_idx = np.argsort(-np.abs(eigs))
        eigs_sorted = eigs[order_idx]
        R_results[k_lower] = eigs_sorted

        print(f"  R_{k_lower} computed + eigvals in {t_eig:.2f}s")
        print(f"  top 10 eigenvalues by |·| (truncated; full list in CSV):")
        for j, e in enumerate(eigs_sorted[:10]):
            mag = abs(e)
            arg = np.angle(e)
            print(f"    {j+1:>2}: {e.real:+.6f} "
                  f"{'+' if e.imag >= 0 else '-'} {abs(e.imag):.6f}i  "
                  f"|·| = {mag:.6f}, arg = {arg:+.4f} ({np.degrees(arg):+6.1f}°)")
        # Identify oscillating pairs (complex with non-trivial argument)
        complex_pairs = [(j, e) for j, e in enumerate(eigs_sorted)
                          if abs(e.imag) > 1e-8]
        if complex_pairs:
            print(f"  largest |·| complex eigenvalue: ", end="")
            j, e = complex_pairs[0]
            print(f"|z| = {abs(e):.6f}, arg = {np.angle(e):+.4f} rad "
                  f"({np.degrees(np.angle(e)):+.1f}°), "
                  f"period = {2*np.pi/abs(np.angle(e)):.2f} steps")
        print()

    # Track top-magnitude eigenvalue across k
    print(f"--- Cross-k top eigenvalue tracking ---")
    print(f"{'k':>3} {'top |z|':>10} {'2nd |z|':>10} {'top complex |z|':>18} "
          f"{'top complex arg':>18}")
    for k in sorted(R_results.keys()):
        eigs = R_results[k]
        top_mag = abs(eigs[0])
        second_mag = abs(eigs[1])
        complex_eigs = [e for e in eigs if abs(e.imag) > 1e-8]
        if complex_eigs:
            top_c = max(complex_eigs, key=lambda x: abs(x))
            top_c_mag = abs(top_c)
            top_c_arg = np.angle(top_c)
            top_c_period = 2 * np.pi / abs(top_c_arg) if abs(top_c_arg) > 1e-12 else float('inf')
        else:
            top_c_mag = float('nan')
            top_c_arg = float('nan')
            top_c_period = float('nan')
        print(f"{k:>3} {top_mag:>10.6f} {second_mag:>10.6f} "
              f"{top_c_mag:>18.6f} {top_c_arg:>+18.4f}")

    # =================== Outputs ===================
    print(f"\nWriting outputs...")
    with open(OUT_RECUR_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["order", "n_eq", "r2", "ss_res",
                    "alphas (comma-sep)", "roots (comma-sep)"])
        for order, fit in sorted(fits.items()):
            alphas_str = ",".join(f"{a:.10f}" for a in fit["alphas"])
            roots_str = ",".join(f"({r.real:+.6f}{r.imag:+.6f}j)"
                                  for r in fit["roots"])
            w.writerow([order, fit["n_eq"], f"{fit['r2']:.6f}",
                        f"{fit['ss_res']:.4e}", alphas_str, roots_str])
    print(f"  saved {OUT_RECUR_CSV}")

    with open(OUT_RSPEC_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["k_lower", "rank", "real", "imag", "magnitude", "arg_rad",
                    "arg_deg", "period_if_complex"])
        for k, eigs in sorted(R_results.items()):
            for j, e in enumerate(eigs):
                mag = abs(e)
                arg = np.angle(e)
                period = (2 * np.pi / abs(arg)) if abs(arg) > 1e-12 \
                    else float("inf")
                w.writerow([k, j + 1, f"{e.real:.10f}", f"{e.imag:.10f}",
                            f"{mag:.10f}", f"{arg:.10f}",
                            f"{np.degrees(arg):.6f}",
                            f"{period:.4f}"])
    print(f"  saved {OUT_RSPEC_CSV}")

    # ==================== Markdown ====================
    L = []
    L.append("# Renormalization Spectrum — recurrence fit + R operator")
    L.append("")
    L.append("**Date:** 2026-05-05.  Two complementary extractions of the "
             "slow-mode spectrum governing ε_k.")
    L.append("")
    L.append("## Part 1 — Linear-recurrence fit on ε_k")
    L.append("")
    L.append("Fitting ε_{k+order} = α_1 ε_{k+order-1} + ... + α_order ε_k by "
             "OLS on consecutive (k, ε_k) pairs from k=2..11 (10 data points).")
    L.append("")
    L.append("| order | n_eq | R² | top \\|root\\| | next \\|root\\| | "
             "complex pair? |")
    L.append("|---:|---:|---:|---:|---:|---|")
    for order in sorted(fits.keys()):
        fit = fits[order]
        roots = fit["roots"]
        mags = sorted([abs(r) for r in roots], reverse=True)
        complex_pairs = [r for r in roots if abs(r.imag) > 1e-8]
        cstr = "—"
        if complex_pairs:
            top_c = max(complex_pairs, key=lambda x: abs(x))
            cstr = (f"\\|z\\| = {abs(top_c):.4f}, "
                    f"arg = {np.angle(top_c):+.3f} rad "
                    f"(period {2*np.pi/abs(np.angle(top_c)):.2f})")
        L.append(f"| {order} | {fit['n_eq']} | {fit['r2']:.5f} | "
                 f"{mags[0]:.4f} | "
                 f"{mags[1] if len(mags) > 1 else float('nan'):.4f} | "
                 f"{cstr} |")
    L.append("")
    L.append("### Forward predictions ε_12..ε_15")
    L.append("")
    L.append("| k | order=1 | order=2 | order=3 | order=4 |")
    L.append("|---:|---:|---:|---:|---:|")
    for k_p in [12, 13, 14, 15]:
        row = [str(k_p)]
        for order in [1, 2, 3, 4]:
            if order in pred_table:
                row.append(f"{pred_table[order][k_p]:+.4e}")
            else:
                row.append("—")
        L.append("| " + " | ".join(row) + " |")
    L.append("")

    # Per-order details
    for order in sorted(fits.keys()):
        fit = fits[order]
        L.append(f"### Order {order} details")
        L.append("")
        L.append(f"- α coefficients: " +
                 ", ".join(f"α_{j+1} = {a:+.6f}"
                            for j, a in enumerate(fit['alphas'])))
        L.append(f"- R² = {fit['r2']:.6f}, ss_res = {fit['ss_res']:.4e}")
        L.append(f"- Characteristic roots:")
        for j, root in enumerate(fit["roots"]):
            mag = abs(root)
            arg = np.angle(root)
            L.append(f"  - z_{j+1} = {root.real:+.6f} "
                     f"{'+' if root.imag >= 0 else '-'} "
                     f"{abs(root.imag):.6f}i,  |z| = {mag:.6f}, "
                     f"arg = {arg:+.4f} rad ({np.degrees(arg):+.1f}°)")
        L.append("")

    # Part 2
    L.append("## Part 2 — R operator at k = 3..6")
    L.append("")
    L.append("R_k acts on W_k = V_{k+1} ⊖ T(V_k), the wavelet residual "
             "subspace at level k. R_k = P_W K_{k+1} P_W^T  (orthonormal "
             "basis projection).")
    L.append("")
    L.append("### Cross-k tracking of leading R_k eigenvalues")
    L.append("")
    L.append("| k | dim(W_k) | top \\|z\\| | 2nd \\|z\\| | "
             "top complex \\|z\\| | top complex arg | period |")
    L.append("|---:|---:|---:|---:|---:|---:|---:|")
    for k in sorted(R_results.keys()):
        eigs = R_results[k]
        dim_W = len(eigs)
        top_mag = abs(eigs[0])
        second_mag = abs(eigs[1])
        complex_eigs = [e for e in eigs if abs(e.imag) > 1e-8]
        if complex_eigs:
            top_c = max(complex_eigs, key=lambda x: abs(x))
            top_c_mag = abs(top_c)
            top_c_arg = np.angle(top_c)
            top_c_period = (2 * np.pi / abs(top_c_arg)
                            if abs(top_c_arg) > 1e-12 else float("inf"))
            L.append(f"| {k} | {dim_W} | {top_mag:.6f} | "
                     f"{second_mag:.6f} | {top_c_mag:.6f} | "
                     f"{top_c_arg:+.4f} ({np.degrees(top_c_arg):+.1f}°) | "
                     f"{top_c_period:.2f} |")
        else:
            L.append(f"| {k} | {dim_W} | {top_mag:.6f} | "
                     f"{second_mag:.6f} | — | — | — |")
    L.append("")

    L.append("### Top 10 eigenvalues of R_k by |·|")
    for k in sorted(R_results.keys()):
        L.append("")
        L.append(f"#### R_{k}  (dim = {len(R_results[k])})")
        L.append("")
        L.append("| rank | Re | Im | |z| | arg (rad) | arg (deg) | "
                 "period |")
        L.append("|---:|---:|---:|---:|---:|---:|---:|")
        for j, e in enumerate(R_results[k][:10]):
            mag = abs(e)
            arg = np.angle(e)
            period = (2 * np.pi / abs(arg)) if abs(arg) > 1e-12 \
                else float("inf")
            L.append(f"| {j+1} | {e.real:+.6f} | {e.imag:+.6f} | "
                     f"{mag:.6f} | {arg:+.4f} | "
                     f"{np.degrees(arg):+.1f} | {period:.2f} |")
    L.append("")

    L.append("## Files")
    L.append("")
    L.append("- `result_renormalization_recurrence_fits.csv` — per-order "
             "fits with roots")
    L.append("- `result_renormalization_R_spectrum.csv` — full R_k spectra "
             "at k=3..6")
    L.append("- `result_renormalization_spectrum.md` — this writeup")

    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    print(f"  saved {OUT_MD}")
    print("\nDone.")


if __name__ == "__main__":
    main()
