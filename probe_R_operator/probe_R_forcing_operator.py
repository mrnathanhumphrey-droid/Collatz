"""
probe_R_forcing_operator.py
============================
Inter-level forcing operator R_k: V_k -> W_{k+1}.

Setup:
  V_k         = level-k coprime function space (dim n_k = 2*3^{k-1})
  L_k         = uniform lift V_k -> V_{k+1} via L_k[r, r'] = 1/3 for the
                3 coprime preimages of r in (Z/3^{k+1})*; row-vector convention
                u @ L_k preserves total mass.
  K_{k+1}     = Tao-Syracuse Markov kernel on V_{k+1}
  W_{k+1}     = V_{k+1} \ominus L_k(V_k), orthogonal complement of the lift
                image inside V_{k+1} (dim = n_{k+1} - n_k = 2*n_k)
  P_{W_{k+1}} = orthonormal projection V_{k+1} -> W_{k+1}
  R_k         = P_{W_{k+1}} ∘ K_{k+1} ∘ L_k

Brief said "V_k = span(pi_k)" (1-dim) but then requested top-20 singular values;
the only consistent reading uses V_k = full level-k function space (n_k-dim),
matching the prior R operator probe's notation. Documented in the markdown.

Levels k = 5, 6, 7. K_8 (n=4374) is the largest dense matrix; SVD at
k=7 expected ~30-60s on this hardware.

Outputs (to C:\\Collatz\\probe_R_operator\\):
- R_k{k}_singular_values.csv     (rank, sigma per k, top 20 plus tail summary)
- R_k{k}_dominant_vectors.npz    (u_1, v_1, sigma_1)
- R_operator_findings.md         (writeup with verdict)
"""
from __future__ import annotations

import csv
import math
import os
import sys
import time

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz\probe_R_operator"
os.makedirs(OUT_DIR, exist_ok=True)

LEVELS = [5, 6, 7]
TOP_N_REPORT = 20
Q = 3  # 3x+1 case


# ---------- kernel + lift + projection (from prior probe, verified) ----------

def order_of_two(N: int) -> int:
    assert N % 2 == 1
    m = 1
    v = 2 % N
    while v != 1:
        v = (v * 2) % N
        m += 1
    return m


def build_K_dense(q: int, k: int):
    """Tao Markov kernel K_k as right-stochastic dense matrix on coprime
    classes of Z/q^k. Row r sums to 1; (u @ K) is the next-step row vector."""
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
    Z_v = 1.0 - 2.0 ** (-min(M, 1074))
    weights = np.zeros(M, dtype=np.float64)
    for vv in range(min(M, 1074)):
        weights[vv] = (2.0 ** -(vv + 1)) / Z_v
    K = np.zeros((n, n), dtype=np.float64)
    for i_r in range(n):
        r = int(coprime[i_r])
        base = (q * r + 1) % N
        targets = (base * powers_inv2) % N
        js = state_idx[targets]
        K[i_r] = np.bincount(js, weights=weights, minlength=n)
    return K, coprime, state_idx


def build_T_lift(q: int, k_lower: int):
    """T: V_k -> V_{k+1} as (n_k, n_{k+1}) matrix.
    T[i, j] = 1/q if coprime_kp1[j] mod q^k == coprime_k[i], else 0.
    Each row has q nonzero entries summing to 1; row-vector u @ T is the
    uniform lift preserving total mass."""
    k = k_lower
    N_k = q ** k
    N_kp1 = q ** (k + 1)
    coprime_k = [r for r in range(N_k) if r % q != 0]
    coprime_kp1 = [r for r in range(N_kp1) if r % q != 0]
    n_k = len(coprime_k)
    n_kp1 = len(coprime_kp1)
    state_idx_k = {r: i for i, r in enumerate(coprime_k)}
    state_idx_kp1 = {r: i for i, r in enumerate(coprime_kp1)}
    T = np.zeros((n_k, n_kp1), dtype=np.float64)
    inv_q = 1.0 / q
    for r_p in coprime_kp1:
        i_p = state_idx_kp1[r_p]
        r = r_p % N_k
        i = state_idx_k[r]
        T[i, i_p] = inv_q
    return T, coprime_k, coprime_kp1


def build_W_basis(T: np.ndarray):
    """Given T of shape (n_k, n_{k+1}), return P_W of shape
    (dim_W = n_{k+1} - n_k, n_{k+1}) with orthonormal rows spanning
    W_{k+1} = (row-space of T)^perp inside R^{n_{k+1}}."""
    Q, _ = np.linalg.qr(T.T, mode="complete")
    n_k = T.shape[0]
    P_W = Q[:, n_k:].T
    return P_W


# ---------- per-k probe ----------

def probe_one(k_lower: int):
    """Compute R_k = P_W K_{k+1} L_k as forcing operator V_k -> W_{k+1}.
    Returns dict with singular values, dominant u/v vectors, timings."""
    k_upper = k_lower + 1

    t0 = time.time()
    K_kp1, _, _ = build_K_dense(Q, k_upper)
    t_K = time.time() - t0
    n_kp1 = K_kp1.shape[0]
    print(f"  K_{k_upper} built: shape {K_kp1.shape}, t={t_K:.1f}s")

    t0 = time.time()
    T, _, _ = build_T_lift(Q, k_lower)
    n_k = T.shape[0]
    P_W = build_W_basis(T)
    dim_W = P_W.shape[0]
    t_TP = time.time() - t0
    print(f"  T={T.shape}, dim(V_{k_lower})={n_k}, dim(W_{k_upper})={dim_W}, "
          f"t={t_TP:.1f}s")

    # Sanity check: P_W @ T.T should be ~ 0 (W is orthogonal to row-space of T)
    res_orth = float(np.max(np.abs(P_W @ T.T)))
    print(f"  orthogonality check max|P_W @ T^T|: {res_orth:.2e}")

    # Build R_k as column-vector linear map V_k -> W_{k+1}
    # Row-vector form u @ T @ K @ P_W^T;  column-vector form is its transpose:
    # R_k_col = P_W @ K^T @ T^T, shape (dim_W, n_k)
    t0 = time.time()
    R_k = P_W @ K_kp1.T @ T.T
    t_build_R = time.time() - t0
    print(f"  R_k (col-vec form) shape {R_k.shape}, t_build={t_build_R:.1f}s")

    # SVD: R_k = U @ diag(sigma) @ Vh.  Thin SVD since dim_W >= n_k.
    t0 = time.time()
    U, sigma, Vh = np.linalg.svd(R_k, full_matrices=False)
    t_svd = time.time() - t0
    print(f"  SVD: {len(sigma)} singular values, "
          f"sigma_1={sigma[0]:.6f}, sigma_min={sigma[-1]:.6e}, t={t_svd:.1f}s")

    # Sanity: sigma_max should be <= 1 (K is stochastic, T sums to 1, P_W is
    # orthogonal projection).
    return {
        "k_lower": k_lower,
        "n_k": n_k,
        "dim_W": dim_W,
        "sigma": sigma,
        "U_top": U[:, :TOP_N_REPORT].copy(),  # dominant left singular vecs in W
        "Vh_top": Vh[:TOP_N_REPORT, :].copy(),  # dominant right singular vecs in V
        "u_1": U[:, 0].copy(),
        "v_1": Vh[0, :].copy(),
        "sigma_1": float(sigma[0]),
        "res_orth": res_orth,
        "t_K": t_K, "t_TP": t_TP, "t_build_R": t_build_R, "t_svd": t_svd,
    }


# ---------- write outputs ----------

def write_singular_values_csv(res: dict, k: int):
    path = os.path.join(OUT_DIR, f"R_k{k}_singular_values.csv")
    sigma = res["sigma"]
    n = len(sigma)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "sigma", "log10_sigma"])
        for i, s in enumerate(sigma):
            ls = math.log10(max(s, 1e-300))
            w.writerow([i + 1, f"{s:.15e}", f"{ls:.6f}"])
    print(f"  csv: {path} ({n} singular values)")


def write_dominant_vectors_npz(res: dict, k: int):
    path = os.path.join(OUT_DIR, f"R_k{k}_dominant_vectors.npz")
    np.savez(
        path,
        sigma_top=res["sigma"][:TOP_N_REPORT],
        u_top=res["U_top"],
        v_top=res["Vh_top"],
        sigma_full=res["sigma"],
    )
    print(f"  npz: {path}")


def write_findings_md(all_res: dict):
    path = os.path.join(OUT_DIR, "R_operator_findings.md")
    rho_slow = 0.826934  # from order-3 recurrence (result_renormalization_recurrence_fits.csv)
    order3_complex_mag = 0.192080
    order2_top = 0.312245
    rate_one_half = 0.5

    def closest_to(x, candidates):
        best = min(candidates.items(), key=lambda kv: abs(x - kv[1]))
        return best[0], best[1], abs(x - best[1])

    references = {
        "rho_slow (order-3 real root)": rho_slow,
        "rate-1/2": rate_one_half,
        "order-3 complex magnitude": order3_complex_mag,
        "order-2 top root": order2_top,
        "1.0 (broken projection)": 1.0,
        "0.0 (trivial forcing)": 0.0,
    }

    lines = []
    lines.append("# R_k forcing operator — singular value spectrum")
    lines.append("")
    lines.append("**Date:** 2026-05-05 (post-compact). Probe target: forcing "
                 "operator R_k = P_{W_{k+1}} ∘ K_{k+1} ∘ L_k (V_k -> W_{k+1}), "
                 "where W_{k+1} = V_{k+1} \\ominus L_k(V_k). Tests whether the "
                 "rate-determining structure missing from the prior R-operator "
                 "probe (which used the W -> W self-map and returned zero) "
                 "lives in the V -> W forcing block.")
    lines.append("")
    lines.append("## Convention note")
    lines.append("")
    lines.append("Brief specified `V_k = span(pi_k)` (1-dim) but also requested "
                 "top-20 singular values — those two are inconsistent (a 1-dim "
                 "input gives at most one nonzero singular value). Adopted the "
                 "interpretation V_k = full level-k coprime function space "
                 "(n_k = 2·3^{k-1}-dim), consistent with the prior R-operator "
                 "probe's V_k notation. R_k is therefore an n_k-dim → 2·n_k-dim "
                 "linear map and yields n_k singular values.")
    lines.append("")
    lines.append("## Construction")
    lines.append("")
    lines.append("- Row-vector convention: `u @ L_k @ K_{k+1} @ P_{W_{k+1}}^T` "
                 "transports a level-k function through the lift, the next-level "
                 "Markov dynamics, and the orthogonal projection onto W_{k+1}.")
    lines.append("- L_k is the uniform fiber lift L_k[i, j] = 1/3 if "
                 "coprime_{k+1}[j] mod 3^k = coprime_k[i] else 0; row-stochastic.")
    lines.append("- W_{k+1} basis P_W obtained from QR of L_k^T (orthonormal "
                 "complement of L_k's row-space inside V_{k+1}). Sanity check "
                 "max|P_W @ L_k^T| reported per k below — should be ~ machine eps.")
    lines.append("- R_k as a column-vector linear map V_k -> W_{k+1} is "
                 "`R_k = P_W @ K_{k+1}^T @ L_k^T`, shape (dim(W_{k+1}), n_k). "
                 "SVD via `numpy.linalg.svd(R_k, full_matrices=False)`.")
    lines.append("")
    lines.append("## Pre-registered reference rates")
    lines.append("")
    lines.append("| label | value | source |")
    lines.append("|---|---|---|")
    lines.append(f"| rho_slow (order-3 real root) | {rho_slow:.6f} | "
                 f"`result_renormalization_recurrence_fits.csv` |")
    lines.append(f"| order-3 complex pair magnitude | {order3_complex_mag:.6f} | "
                 f"same |")
    lines.append(f"| order-2 top root | {order2_top:.6f} | same |")
    lines.append(f"| rate-1/2 (legacy walked-back claim) | {rate_one_half:.4f} | "
                 f"R75/R76 |")
    lines.append(f"| 1.0 | broken projection diagnostic | walk-back gate |")
    lines.append(f"| 0.0 | trivial forcing diagnostic | walk-back gate |")
    lines.append("")

    # Per-k summary
    lines.append("## Per-k summary")
    lines.append("")
    lines.append("| k | dim(V_k) | dim(W_{k+1}) | sigma_1 | sigma_2 | sigma_3 "
                 "| sigma_20 | sigma_min | orth_resid |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for k in sorted(all_res.keys()):
        r = all_res[k]
        s = r["sigma"]
        n_top = lambda i: f"{s[i]:.6e}" if i < len(s) else "--"
        lines.append(f"| {k} | {r['n_k']} | {r['dim_W']} | "
                     f"{n_top(0)} | {n_top(1)} | {n_top(2)} | {n_top(19)} | "
                     f"{s[-1]:.6e} | {r['res_orth']:.2e} |")
    lines.append("")

    # Top 20 singular values per k
    lines.append("## Top 20 singular values per k")
    lines.append("")
    for k in sorted(all_res.keys()):
        r = all_res[k]
        lines.append(f"### k = {k}  (V_k dim = {r['n_k']}, "
                     f"W_{k+1} dim = {r['dim_W']})")
        lines.append("")
        lines.append("| rank | sigma | log10(sigma) |")
        lines.append("|---|---|---|")
        for i in range(min(TOP_N_REPORT, len(r["sigma"]))):
            s = r["sigma"][i]
            lines.append(f"| {i+1} | {s:.10e} | {math.log10(max(s, 1e-300)):.4f} |")
        lines.append("")

    # Cross-k tracking + closeness to reference rates
    lines.append("## sigma_1 vs reference rates")
    lines.append("")
    lines.append("| k | sigma_1 | closest reference | distance |")
    lines.append("|---|---|---|---|")
    for k in sorted(all_res.keys()):
        s1 = all_res[k]["sigma_1"]
        label, ref_val, dist = closest_to(s1, references)
        lines.append(f"| {k} | {s1:.6f} | {label} ({ref_val:.4f}) | {dist:.4f} |")
    lines.append("")

    # Cross-k decay/growth
    lines.append("## Cross-k decay/growth of sigma_1")
    lines.append("")
    lines.append("| k -> k+1 | sigma_1(k) | sigma_1(k+1) | ratio | log_3 ratio |")
    lines.append("|---|---|---|---|---|")
    ks_sorted = sorted(all_res.keys())
    for i in range(len(ks_sorted) - 1):
        ka, kb = ks_sorted[i], ks_sorted[i+1]
        sa = all_res[ka]["sigma_1"]
        sb = all_res[kb]["sigma_1"]
        ratio = sb / sa if sa > 0 else float("inf")
        lr = math.log(ratio, 3) if ratio > 0 else float("-inf")
        lines.append(f"| {ka} -> {kb} | {sa:.6f} | {sb:.6f} | "
                     f"{ratio:.6f} | {lr:.4f} |")
    lines.append("")

    # Verdict
    lines.append("## Verdict")
    lines.append("")
    sigma1s = [all_res[k]["sigma_1"] for k in ks_sorted]
    near_one = all(abs(s - 1.0) < 0.05 for s in sigma1s)
    near_zero = all(s < 0.01 for s in sigma1s)
    near_rho = any(abs(s - rho_slow) < 0.05 for s in sigma1s)
    near_half = any(abs(s - 0.5) < 0.05 for s in sigma1s)

    if near_one:
        lines.append("**Walk-back gate fires (sigma_1 ≈ 1):** the projection "
                     "is not orthogonalizing as intended, or R_k has a "
                     "preserved direction. Construction needs review.")
    elif near_zero:
        lines.append("**Walk-back gate fires (sigma_1 ≈ 0):** the V -> W "
                     "forcing is trivial; rate-determining structure is not "
                     "located in this block. Look at W -> W self-map (already "
                     "null), or higher-order forcing across multiple levels.")
    elif near_rho:
        lines.append(f"**Pre-registered match: sigma_1 hits rho_slow ≈ "
                     f"{rho_slow:.3f}.** Strong evidence that R_k carries "
                     "the rate-determining structure of ε_k convergence. The "
                     "slow mode in the order-3 recurrence is an SVD direction "
                     "of the V -> W forcing operator.")
    elif near_half:
        lines.append("**sigma_1 near 1/2:** consistent with the legacy "
                     "rate-1/2 envelope claim that has been empirically "
                     "walked back at the ε_k scalar level. Whether R_k's "
                     "spectrum is structurally 1/2 or this is a finite-k "
                     "artifact is the next question.")
    else:
        lines.append("**No clean match to any pre-registered reference.** "
                     "Singular values cluster but do not match rho_slow ≈ 0.83, "
                     "rate 1/2, or any order-2/3 recurrence root within 0.05. "
                     "Reporting and pausing for analysis as per the brief's "
                     "third walk-back gate.")
    lines.append("")
    lines.append(f"sigma_1 across k = {ks_sorted}: " +
                 ", ".join(f"{s:.4f}" for s in sigma1s))
    lines.append("")

    # Timings
    lines.append("## Timings")
    lines.append("")
    lines.append("| k | t_K | t_lift+P_W | t_build_R | t_SVD |")
    lines.append("|---|---|---|---|---|")
    for k in ks_sorted:
        r = all_res[k]
        lines.append(f"| {k} | {r['t_K']:.1f}s | {r['t_TP']:.1f}s | "
                     f"{r['t_build_R']:.1f}s | {r['t_svd']:.1f}s |")
    lines.append("")

    lines.append("## Files")
    lines.append("")
    for k in ks_sorted:
        lines.append(f"- `R_k{k}_singular_values.csv` — full singular spectrum")
        lines.append(f"- `R_k{k}_dominant_vectors.npz` — top-{TOP_N_REPORT} u/v vectors")
    lines.append("- `R_operator_findings.md` — this writeup")
    lines.append("")

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print(f"\n[md: {path}]")


def main():
    print("=" * 78)
    print("R_k forcing operator probe (V_k -> W_{k+1})")
    print(f"Levels: {LEVELS}")
    print("=" * 78)

    all_res = {}
    for k in LEVELS:
        print(f"\n--- k = {k} ---")
        t_start = time.time()
        res = probe_one(k)
        all_res[k] = res
        write_singular_values_csv(res, k)
        write_dominant_vectors_npz(res, k)
        print(f"  total t={time.time()-t_start:.1f}s")

    write_findings_md(all_res)
    print("\nDone.")


if __name__ == "__main__":
    main()
