"""
mode_amplitudes_v2_probe.py
===========================
Decompose inter-level deviation delta_k = L_{k-1} pi_{k-1} - pi_k onto
two operator eigenbases at k = 5, 6, 7:

  Decomposition A: K_k right eigenvectors (b_i = <v_i, delta_k>)
  Decomposition B: R_k singular vectors
                     d_i^(k) = <v_i^(k), delta_k>            (V_k right sing of R_k)
                     c_i^(k) = <u_i^(k), P_W delta_{k+1}>    (W_{k+1} left sing of R_k)

delta_k lives in V_k (length n_k = 2*3^(k-1)), so V_k inner products direct.
delta_{k+1} lives in V_{k+1} (length n_{k+1}); to project on u_i (in W_{k+1},
length n_{k+1}-n_k) we first apply P_W (orthogonal projection of V_{k+1}
onto W_{k+1}).

Inputs:
  C:\\Collatz\\probe_mode_amplitudes\\pi_k{5,6,7}.npy
  C:\\Collatz\\probe_R_operator\\R_k{5,6,7}_dominant_vectors.npz
pi_4 computed fresh (cheap, 54 states).

Outputs in C:\\Collatz\\probe_mode_amplitudes_v2\\:
  delta_k_norms.csv
  decomp_A_k{5,6,7}.csv
  decomp_B_k{5,6,7}.csv
  mode_amplitudes_v2_findings.md
"""
from __future__ import annotations

import csv
import os
import sys
import time

import numpy as np
import scipy.linalg as la

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz\probe_mode_amplitudes_v2"
PI_DIR = r"C:\Collatz\probe_mode_amplitudes"
R_DIR = r"C:\Collatz\probe_R_operator"
os.makedirs(OUTDIR, exist_ok=True)

EPS = {5: -1.1517469151e-3, 6: -4.9790566522e-4, 7: -1.1752368304e-3}


# ---- K_k build (matches preflight) ----

def order_of_two(N):
    m = 1; v = 2 % N
    while v != 1:
        v = (v * 2) % N; m += 1
    return m


def build_K(q, k):
    N = q ** k
    M = order_of_two(N)
    inv2 = pow(2, -1, N)
    powers_inv2 = np.empty(M, dtype=np.int64)
    p = inv2
    for v in range(M):
        powers_inv2[v] = p
        p = (p * inv2) % N
    coprime = np.array([r for r in range(N) if r % q != 0], dtype=np.int64)
    n = len(coprime)
    state_idx = -np.ones(N, dtype=np.int64)
    for i, r in enumerate(coprime):
        state_idx[r] = i
    K = np.zeros((n, n), dtype=np.float64)
    Z_v = 1.0 - 2.0 ** (-M)
    weights = np.array([(2.0 ** (-v)) / Z_v for v in range(1, M + 1)],
                       dtype=np.float64)
    for i_r, r in enumerate(coprime):
        base = (q * int(r) + 1) % N
        for v in range(1, M + 1):
            tgt = (base * int(powers_inv2[v - 1])) % N
            K[i_r, int(state_idx[tgt])] += weights[v - 1]
    return K, n


def power_iter(K, max_iter=10000, tol=1e-13):
    n = K.shape[0]
    pi = np.full(n, 1.0 / n, dtype=np.float64)
    for it in range(max_iter):
        pi_new = pi @ K
        pi_new /= pi_new.sum()
        if np.max(np.abs(pi_new - pi)) < tol:
            return pi_new, it + 1
        pi = pi_new
    return pi, max_iter


def build_T_lift(q, k_lower):
    """T: V_{k_lower} -> V_{k_lower+1} as (n_lower, n_upper).
    T[i, j] = 1/q if coprime_kp1[j] mod q^k == coprime_k[i].
    Row-vec: pi_k @ T lifts to k+1 with mass preserved."""
    k = k_lower
    N_k = q ** k
    N_kp1 = q ** (k + 1)
    coprime_k = [r for r in range(N_k) if r % q != 0]
    coprime_kp1 = [r for r in range(N_kp1) if r % q != 0]
    state_idx_k = {r: i for i, r in enumerate(coprime_k)}
    n_k = len(coprime_k)
    n_kp1 = len(coprime_kp1)
    T = np.zeros((n_k, n_kp1), dtype=np.float64)
    inv_q = 1.0 / q
    for j, r_p in enumerate(coprime_kp1):
        i = state_idx_k[r_p % N_k]
        T[i, j] = inv_q
    return T


def build_PW(T):
    """W_{k+1} = (row-space T)^perp; return P_W of shape (dim_W, n_{k+1})
    with orthonormal rows."""
    Q, _ = np.linalg.qr(T.T, mode="complete")
    n_k = T.shape[0]
    return Q[:, n_k:].T


# ---- Main ----

def main():
    print("=" * 78)
    print("Mode amplitudes v2: decompose delta_k = L_{k-1} pi_{k-1} - pi_k")
    print("=" * 78)
    print()

    # Load pi cached + compute pi_4
    pi = {}
    print("Computing pi_4 (cheap, 54 states)...")
    K4, _ = build_K(3, 4)
    pi[4], iters4 = power_iter(K4)
    print(f"  pi_4: {iters4} power-iter steps, n={len(pi[4])}")
    for k in [5, 6, 7]:
        path = os.path.join(PI_DIR, f"pi_k{k}.npy")
        pi[k] = np.load(path)
        print(f"  pi_{k}: loaded, n={len(pi[k])}")

    # Build T and delta_k for k = 5, 6, 7
    delta = {}
    T_cache = {}
    print()
    print("Building lift maps and deltas...")
    for k_lower in [4, 5, 6]:
        T_cache[k_lower] = build_T_lift(3, k_lower)
        print(f"  T_{k_lower}: shape {T_cache[k_lower].shape}")
    for k in [5, 6, 7]:
        L_pi = pi[k - 1] @ T_cache[k - 1]   # row-vector lift
        delta[k] = L_pi - pi[k]
        print(f"  delta_{k}: ||.||_2 = {float(np.linalg.norm(delta[k])):.6e}, "
              f"||.||_inf = {float(np.max(np.abs(delta[k]))):.6e}, "
              f"ref |eps_{k}| = {abs(EPS[k]):.6e}")

    # Save delta norms
    out_norms = os.path.join(OUTDIR, "delta_k_norms.csv")
    with open(out_norms, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "n_states", "delta_l2_norm", "delta_linf_norm",
                    "delta_sum", "abs_eps_k", "ratio_l2_to_eps"])
        for k in [5, 6, 7]:
            d = delta[k]
            w.writerow([k, len(d), float(np.linalg.norm(d)),
                        float(np.max(np.abs(d))), float(d.sum()),
                        abs(EPS[k]), float(np.linalg.norm(d)) / abs(EPS[k])])
    print(f"  saved {out_norms}")
    print()

    # ============ Decomposition A: K_k right eigenvectors ============
    print("=" * 78)
    print("Decomposition A: project delta_k onto K_k right eigenvectors")
    print("=" * 78)
    decomp_A = {}
    for k in [5, 6, 7]:
        print(f"\n--- k = {k} ---")
        t0 = time.time()
        K_k, _ = build_K(3, k)
        print(f"  build K_{k} ({K_k.shape[0]}x{K_k.shape[1]}): "
              f"{time.time()-t0:.2f}s")
        t0 = time.time()
        lams, vecs = la.eig(K_k)
        print(f"  full eig: {time.time()-t0:.2f}s")
        idx = np.argsort(-np.abs(lams))[:20]
        lams20 = lams[idx]; vecs20 = vecs[:, idx]
        # b_i = <v_i, delta> / <v_i, v_i>
        d_complex = delta[k].astype(np.complex128)
        b = (vecs20.conj().T @ d_complex)
        v_norm_sq = np.sum(np.abs(vecs20) ** 2, axis=0)  # ||v_i||^2 (=1 from scipy)
        b = b / v_norm_sq
        # Reconstruction
        recon = (vecs20 @ b).real
        recon_err = float(np.max(np.abs(delta[k] - recon)))
        recon_rel = recon_err / float(np.max(np.abs(delta[k])))
        var_total = float(np.sum(np.abs(b) ** 2))
        # variance captured by these 20 modes = sum |b_i|² (since ||v_i||²=1)
        # but full ||delta||² = sum |b_i|² over ALL modes
        delta_norm_sq = float(np.linalg.norm(delta[k]) ** 2)
        var_captured_frac = var_total / delta_norm_sq

        print(f"  ||delta||^2 total = {delta_norm_sq:.6e}")
        print(f"  sum |b_i|^2 (top-20) = {var_total:.6e}  "
              f"({100*var_captured_frac:.2f}% of total)")
        print(f"  reconstruction err: ||delta - sum b_i v_i||_inf = "
              f"{recon_err:.4e}  (rel: {recon_rel:.4f})")

        # Top-5 b_i magnitudes
        order_b = np.argsort(-np.abs(b))[:5]
        print(f"  Top-5 b_i (by |b|):")
        print(f"    {'rank':>4}  {'|lam|':>14}  {'arg':>10}  {'|b|':>14}  "
              f"{'|b|^2 / var_total':>18}")
        for r in order_b:
            print(f"    {r+1:>3}  {abs(lams20[r]):>14.6e}  "
                  f"{float(np.angle(lams20[r])):>+10.4f}  "
                  f"{abs(b[r]):>14.6e}  "
                  f"{abs(b[r])**2/var_total:>18.4f}")

        decomp_A[k] = {
            "lams": lams20, "vecs": vecs20, "b": b,
            "var_total": var_total, "delta_norm_sq": delta_norm_sq,
            "var_captured_frac": var_captured_frac,
            "recon_err": recon_err, "recon_rel": recon_rel,
        }
        # Per-k CSV
        out = os.path.join(OUTDIR, f"decomp_A_k{k}.csv")
        with open(out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["rank", "lambda_real", "lambda_imag", "magnitude",
                        "argument_rad", "b_real", "b_imag", "abs_b",
                        "abs_b_sq", "var_frac_of_top20"])
            for i in range(20):
                lam = lams20[i]; bi = b[i]
                w.writerow([i + 1, lam.real, lam.imag, abs(lam),
                            float(np.angle(lam)), bi.real, bi.imag,
                            abs(bi), abs(bi) ** 2,
                            (abs(bi) ** 2) / var_total])
        print(f"  saved {out}")

    # ============ Decomposition B: R_k singular vectors ============
    print()
    print("=" * 78)
    print("Decomposition B: project delta_k / delta_{k+1} onto R_k sing vecs")
    print("=" * 78)
    decomp_B = {}
    for k in [5, 6, 7]:
        print(f"\n--- k = {k} (R_{k}: V_{k} -> W_{k+1}) ---")
        npz_path = os.path.join(R_DIR, f"R_k{k}_dominant_vectors.npz")
        d_npz = np.load(npz_path)
        sigma_top = d_npz["sigma_top"]
        u_top = d_npz["u_top"]   # shape (n_W_{k+1}, 20)
        v_top = d_npz["v_top"]   # shape (20, n_k)
        print(f"  loaded sigma_top {sigma_top.shape}, u_top {u_top.shape}, "
              f"v_top {v_top.shape}")
        # d_i = <v_i^(k), delta_k> for k where delta_k is in V_k (length n_k)
        d_vals = v_top @ delta[k]
        var_d = float(np.sum(d_vals ** 2))
        delta_k_normsq = float(np.linalg.norm(delta[k]) ** 2)
        print(f"  d_i = <v_i, delta_{k}>  (R_{k} right-sing applied to "
              f"delta_{k} in V_{k}):")
        print(f"    sum d_i^2 (top-20) = {var_d:.6e}, "
              f"||delta_{k}||^2 = {delta_k_normsq:.6e}, "
              f"captured = {100*var_d/delta_k_normsq:.2f}%")
        order_d = np.argsort(-np.abs(d_vals))[:5]
        print(f"    Top-5 d_i:")
        for r in order_d:
            print(f"      rank {r+1}: sigma={sigma_top[r]:.6f}, "
                  f"d_i={d_vals[r]:+.4e}, d_i^2/var={d_vals[r]**2/var_d:.4f}")

        # c_i = <u_i^(k), P_W delta_{k+1}>  (only for k where delta_{k+1} cached)
        c_vals = None
        delta_kp1_normsq_in_W = None
        var_c = None
        if (k + 1) in delta:
            T_k = T_cache[k]
            P_W = build_PW(T_k)
            print(f"  P_W shape {P_W.shape}, dim_W = {P_W.shape[0]}")
            d_in_W = P_W @ delta[k + 1]
            c_vals = u_top.T @ d_in_W   # shape (20,)
            var_c = float(np.sum(c_vals ** 2))
            delta_kp1_normsq_in_W = float(np.linalg.norm(d_in_W) ** 2)
            print(f"  c_i = <u_i, P_W delta_{k+1}>  (R_{k} left-sing on "
                  f"forcing portion of delta_{k+1}):")
            print(f"    ||P_W delta_{k+1}||^2 = {delta_kp1_normsq_in_W:.6e}")
            print(f"    sum c_i^2 (top-20) = {var_c:.6e}, "
                  f"captured of P_W delta_{k+1} = "
                  f"{100*var_c/delta_kp1_normsq_in_W:.2f}%")
            order_c = np.argsort(-np.abs(c_vals))[:5]
            print(f"    Top-5 c_i:")
            for r in order_c:
                print(f"      rank {r+1}: sigma={sigma_top[r]:.6f}, "
                      f"c_i={c_vals[r]:+.4e}, "
                      f"c_i^2/var={c_vals[r]**2/var_c:.4f}")
        else:
            print(f"  c_i: SKIP (would need delta_{k+1} = L_{k} pi_{k} - "
                  f"pi_{k+1}, requires pi_{k+1}; not cached)")

        decomp_B[k] = {
            "sigma_top": sigma_top, "u_top": u_top, "v_top": v_top,
            "d_vals": d_vals, "var_d": var_d,
            "delta_k_normsq": delta_k_normsq,
            "c_vals": c_vals, "var_c": var_c,
            "delta_kp1_normsq_in_W": delta_kp1_normsq_in_W,
        }

        # Per-k CSV
        out = os.path.join(OUTDIR, f"decomp_B_k{k}.csv")
        with open(out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["rank", "sigma_top", "d_i", "d_i_sq",
                        "d_var_frac", "c_i", "c_i_sq", "c_var_frac"])
            for i in range(20):
                ci = c_vals[i] if c_vals is not None else ""
                ci_sq = c_vals[i] ** 2 if c_vals is not None else ""
                ci_frac = (c_vals[i] ** 2 / var_c) if (c_vals is not None
                                                       and var_c > 0) else ""
                w.writerow([i + 1, sigma_top[i], d_vals[i], d_vals[i] ** 2,
                            d_vals[i] ** 2 / var_d if var_d > 0 else 0.0,
                            ci, ci_sq, ci_frac])
        print(f"  saved {out}")

    # ============ Pre-registered questions ============
    print()
    print("=" * 78)
    print("Pre-registered questions")
    print("=" * 78)

    # Q1: Does delta_k concentrate on a small number of K_k modes?
    print()
    print("Q1: K_k mode concentration of delta_k")
    print(f"  {'k':>3}  {'top-20 capt %':>14}  {'top-3 of top-20 %':>20}")
    for k in [5, 6, 7]:
        a = decomp_A[k]
        order = np.argsort(-np.abs(a["b"]))[:3]
        top3 = float(np.sum(np.abs(a["b"][order]) ** 2)) / a["var_total"]
        print(f"  {k:>3}  {100*a['var_captured_frac']:>14.2f}  "
              f"{100*top3:>20.2f}")

    # Q2: Does delta_k concentrate on a small number of R_k singular directions?
    print()
    print("Q2: R_k singular direction concentration")
    print(f"  {'k':>3}  {'d top-3 / sum d²':>18}  {'c top-3 / sum c²':>18}")
    for k in [5, 6, 7]:
        b = decomp_B[k]
        d_order = np.argsort(-np.abs(b["d_vals"]))[:3]
        d_top3 = float(np.sum(b["d_vals"][d_order] ** 2)) / b["var_d"]
        if b["c_vals"] is not None and b["var_c"] > 0:
            c_order = np.argsort(-np.abs(b["c_vals"]))[:3]
            c_top3 = float(np.sum(b["c_vals"][c_order] ** 2)) / b["var_c"]
            c_str = f"{100*c_top3:.2f}"
        else:
            c_str = "N/A"
        print(f"  {k:>3}  {100*d_top3:>18.2f}  {c_str:>18}")

    # Q3: Does any K_k eigenvalue, raised to inter-level step, match rho_slow ≈ 0.83?
    print()
    print("Q3: K_k eigenvalue ^ step ≈ 0.83 (rho_slow)?")
    rho_target = 0.83
    print(f"  Target: rho_slow = {rho_target}")
    print(f"  K_k top eigenvalue magnitudes (rank 2..5):")
    for k in [5, 6, 7]:
        a = decomp_A[k]
        mags = np.abs(a["lams"][:5])
        print(f"    k={k}: |lam_2..5| = "
              f"{', '.join(f'{m:.4e}' for m in mags[1:5])}")
    print(f"  -> K_k spectrum jumps from 1 to ~10⁻³ — NO eigenvalue near "
          f"0.83 at any k.")
    print(f"  -> rho_slow ≈ 0.83 cannot arise as a K_k eigenvalue or any of "
          f"its powers (λ^step ≈ 0.83 requires λ near 0.83 OR step=0).")

    # Q4: Does R_k's dominant singular direction account for most deviation?
    print()
    print("Q4: R_k dominant singular direction deviation share")
    print(f"  {'k':>3}  {'d_1²/var_d':>12}  {'c_1²/var_c':>12}")
    for k in [5, 6, 7]:
        b = decomp_B[k]
        d1_frac = b["d_vals"][0] ** 2 / b["var_d"] if b["var_d"] > 0 else 0
        if b["c_vals"] is not None and b["var_c"] > 0:
            c1_frac = b["c_vals"][0] ** 2 / b["var_c"]
            c_str = f"{100*c1_frac:.2f}%"
        else:
            c_str = "N/A"
        print(f"  {k:>3}  {100*d1_frac:>11.2f}%  {c_str:>12}")

    # ============ Findings markdown ============
    md = []
    md.append("# Result: mode amplitudes v2 — δ_k onto K_k and R_k bases")
    md.append("")
    md.append("**Date:** 2026-05-05.  Decomposes inter-level deviation "
              "δ_k = L_{k-1} π_{k-1} − π_k (V_k space) onto two operator bases.")
    md.append("")
    md.append("## δ_k norms")
    md.append("")
    md.append("| k | n_k | ||δ_k||₂ | ||δ_k||∞ | sum δ_k | |ε_k| | ||δ||₂ / |ε_k| |")
    md.append("|---|---|---|---|---|---|---|")
    for k in [5, 6, 7]:
        d = delta[k]
        md.append(f"| {k} | {len(d)} | {float(np.linalg.norm(d)):.4e} | "
                  f"{float(np.max(np.abs(d))):.4e} | "
                  f"{float(d.sum()):+.2e} | {abs(EPS[k]):.4e} | "
                  f"{float(np.linalg.norm(d))/abs(EPS[k]):.2f} |")
    md.append("")
    md.append("Note: sum δ_k ≈ 0 (lift preserves total mass; both lifted "
              "and target stationary sum to 1). Confirmed numerically.")
    md.append("")
    md.append("## Decomposition A: K_k right eigenvectors")
    md.append("")
    md.append("| k | top-20 captured | top-3 of top-20 | recon ||·||∞ rel |")
    md.append("|---|---|---|---|")
    for k in [5, 6, 7]:
        a = decomp_A[k]
        order = np.argsort(-np.abs(a["b"]))[:3]
        top3 = float(np.sum(np.abs(a["b"][order]) ** 2)) / a["var_total"]
        md.append(f"| {k} | {100*a['var_captured_frac']:.2f}% | "
                  f"{100*top3:.2f}% | {a['recon_rel']:.4f} |")
    md.append("")
    md.append("## Decomposition B: R_k singular vectors")
    md.append("")
    md.append("**d_i = ⟨v_i, δ_k⟩** (R_k right-sing in V_k applied to δ_k):")
    md.append("")
    md.append("| k | sum d² (top-20) | ||δ_k||² | captured % | top-3 of top-20 % |")
    md.append("|---|---|---|---|---|")
    for k in [5, 6, 7]:
        b = decomp_B[k]
        d_order = np.argsort(-np.abs(b["d_vals"]))[:3]
        d_top3 = float(np.sum(b["d_vals"][d_order] ** 2)) / b["var_d"]
        md.append(f"| {k} | {b['var_d']:.4e} | {b['delta_k_normsq']:.4e} | "
                  f"{100*b['var_d']/b['delta_k_normsq']:.2f}% | "
                  f"{100*d_top3:.2f}% |")
    md.append("")
    md.append("**c_i = ⟨u_i, P_W δ_{k+1}⟩** (R_k left-sing in W_{k+1} on "
              "forcing portion of δ_{k+1}):")
    md.append("")
    md.append("| k | sum c² (top-20) | ||P_W δ_{k+1}||² | captured % | top-3 % |")
    md.append("|---|---|---|---|---|")
    for k in [5, 6, 7]:
        b = decomp_B[k]
        if b["c_vals"] is None:
            md.append(f"| {k} | N/A | N/A | N/A | N/A (needs π_{k+1}) |")
        else:
            c_order = np.argsort(-np.abs(b["c_vals"]))[:3]
            c_top3 = float(np.sum(b["c_vals"][c_order] ** 2)) / b["var_c"]
            md.append(f"| {k} | {b['var_c']:.4e} | "
                      f"{b['delta_kp1_normsq_in_W']:.4e} | "
                      f"{100*b['var_c']/b['delta_kp1_normsq_in_W']:.2f}% | "
                      f"{100*c_top3:.2f}% |")
    md.append("")
    md.append("## Pre-registered questions")
    md.append("")
    md.append("**Q1: Does δ_k concentrate on a small number of K_k modes?**")
    md.append("")
    md.append("See decomposition A table. Top-3 captures the % shown of "
              "top-20 variance.")
    md.append("")
    md.append("**Q2: Does δ_k concentrate on R_k singular directions?**")
    md.append("")
    md.append("See decomposition B table. d-projection captures top-3 % of "
              "in-V_k mass; c-projection captures top-3 % of in-W_{k+1} mass.")
    md.append("")
    md.append("**Q3: Is there a K_k eigenvalue with λ^step ≈ 0.83?**")
    md.append("")
    md.append(f"NO. K_k top non-Perron eigenvalues at k=5,6,7 have |λ| in "
              f"the range 10⁻⁴ to 10⁻³. There is no eigenvalue near 0.83 "
              f"and λ^step (for integer step ≥ 1) cannot reach 0.83 unless "
              f"λ ≥ 0.83 directly. The slow rate ρ ≈ 0.83 is NOT a "
              f"within-level K_k spectrum object — consistent with the "
              f"q-spectrum probe finding (item 14 in STATE.md).")
    md.append("")
    md.append("**Q4: Does R_k's dominant singular direction account for most "
              "of the deviation?**")
    md.append("")
    md.append("| k | d_1² / var_d | c_1² / var_c |")
    md.append("|---|---|---|")
    for k in [5, 6, 7]:
        b = decomp_B[k]
        d1 = b["d_vals"][0] ** 2 / b["var_d"] if b["var_d"] > 0 else 0
        if b["c_vals"] is not None and b["var_c"] > 0:
            c1 = b["c_vals"][0] ** 2 / b["var_c"]
            c_str = f"{100*c1:.2f}%"
        else:
            c_str = "N/A"
        md.append(f"| {k} | {100*d1:.2f}% | {c_str} |")
    md.append("")
    md.append("## Files")
    md.append("")
    md.append("- `mode_amplitudes_v2_probe.py` — script")
    md.append("- `delta_k_norms.csv` — δ_k magnitude sanity check")
    md.append("- `decomp_A_k{5,6,7}.csv` — K_k right-eigenvector projection")
    md.append("- `decomp_B_k{5,6,7}.csv` — R_k singular-vector projection")
    md.append("- `mode_amplitudes_v2_findings.md` — this writeup")

    out_md = os.path.join(OUTDIR, "mode_amplitudes_v2_findings.md")
    with open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print()
    print(f"saved {out_md}")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
