"""
Experiment 48 — followup to exp 47: re-fit (K_h, b) from body, refit ξ_X.

Result 12 found ξ_X ≈ +0.095 with user-specified K_h=10.4282, b=2.275 — but
those came from σ-std-vs-log(n), a different observable than K_eff_band(q).
At N=2^27, the body K_eff_band(q) values for q ∈ {0.375, 0.625} suggest a
different intercept and slope than the user-specified ones. Re-fit the body
empirically and check whether the tail ξ aligns with σ-records ξ=-0.075
under self-consistent body calibration.

Adds finer tail bands (midpoints at q ∈ {0.825, 0.925, 0.97, 0.995}) for
multi-point GPD fit. Runs at N ∈ {2^24, 2^25, 2^27} to test N-invariance.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from scipy.stats import norm
from scipy.optimize import brentq, minimize

sys.stdout.reconfigure(encoding="utf-8")

K_H_USER = 3.0 / (np.log(4.0) - np.log(3.0))   # 10.4282 user reference
B_USER = 2.275
GPD_XI_REF = -0.075  # σ-records reference


@njit(parallel=True, cache=True)
def first_passage_4thresh(starts, thresholds_per_n):
    n = len(starts)
    s_arr = np.full((n, 4), -1, dtype=np.int32)
    max_steps = 20000
    for i in prange(n):
        m = starts[i]
        steps = 0
        next_idx = 0
        while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
            s_arr[i, next_idx] = 0
            next_idx += 1
        while next_idx < 4 and steps < max_steps:
            if m & 1:
                m = 3 * m + 1
            else:
                m = m >> 1
            steps += 1
            while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
                s_arr[i, next_idx] = steps
                next_idx += 1
    return s_arr


def build_thresholds(n_arr):
    log_n = np.log(n_arr.astype(np.float64))
    sqrt_n = np.sqrt(n_arr.astype(np.float64))
    n_two_third = n_arr.astype(np.float64) ** (2.0 / 3.0)
    sqrt_log = sqrt_n * log_n
    sqrt_div_log = sqrt_n / np.maximum(log_n, 1.0)
    raw = np.column_stack([n_two_third, sqrt_log, sqrt_n, sqrt_div_log])
    sort_idx = np.argsort(-raw, axis=1)
    thresh_sorted = np.take_along_axis(raw, sort_idx, axis=1).astype(np.int64)
    return thresh_sorted, sort_idx, raw


def K_eff_band(sigma, s_per_thresh, thresh_phys):
    log_f = np.log(np.maximum(thresh_phys, 1.0))
    R = sigma[:, None] - s_per_thresh.astype(np.float64)
    log_f_mean = log_f.mean(axis=0)
    R_mean = R.mean(axis=0)
    x = log_f_mean
    y = R_mean
    x_c = x - x.mean()
    y_c = y - y.mean()
    slope = (x_c * y_c).sum() / (x_c * x_c).sum()
    return float(slope)


def gpd_factor(p, xi):
    if abs(xi) < 1e-9:
        return -np.log(1.0 - p)
    return ((1.0 - p) ** (-xi) - 1.0) / xi


def fit_xi_least_squares(q_thresh, q_tail, excess_X):
    """Fit (ξ, σ_GPD) by NLS on tail excess data.
    excess_X[i] = X_q[i] - X_thresh; q_tail[i] in (q_thresh, 1).
    Model: excess = σ_GPD · g(p_i; ξ) where p_i = (q_i - q_thresh)/(1-q_thresh).
    """
    p = (q_tail - q_thresh) / (1.0 - q_thresh)
    def loss(params):
        xi, sigma = params
        if sigma <= 0: return 1e20
        if abs(xi) > 0.5: return 1e20  # bounds
        pred = sigma * np.array([gpd_factor(pi, xi) for pi in p])
        return float(np.sum((pred - excess_X) ** 2))
    res = minimize(loss, x0=[0.0, np.median(excess_X)], method='Nelder-Mead',
                   options={'xatol': 1e-6, 'fatol': 1e-9})
    return float(res.x[0]), float(res.x[1]), float(res.fun)


def analyze_N(N_val, data_dir, n_sample=500_000, rng=None):
    print(f"\n========== N = {N_val:,} = 2^{int(np.log2(N_val))} ==========")
    df = pl.read_parquet(data_dir / f"main_N{N_val}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    sigma_arr = df["sigma"].to_numpy().astype(np.int32)
    print(f"  loaded {len(n_arr):,} odd rows")

    if len(n_arr) > n_sample:
        idx = rng.choice(len(n_arr), size=n_sample, replace=False)
        n_arr = n_arr[idx]; sigma_arr = sigma_arr[idx]

    t0 = time.perf_counter()
    thresh_sorted, sort_idx, raw_thresh = build_thresholds(n_arr)
    s_arr = first_passage_4thresh(n_arr, thresh_sorted)
    inv_sort = np.argsort(sort_idx, axis=1)
    s_phys = np.zeros((len(n_arr), 4), dtype=np.int32)
    for col in range(4):
        rev = inv_sort[:, col]
        s_phys[:, col] = np.take_along_axis(s_arr, rev[:, None], axis=1).flatten()
    print(f"  walked in {time.perf_counter()-t0:.1f}s")

    # Bands at midpoints {0.125, 0.375, 0.625, 0.825, 0.925, 0.97, 0.995}
    # Edges at {25, 50, 75, 90, 95, 99, 100} percentiles
    edges_pct = [25, 50, 75, 90, 95, 99]
    edges = np.percentile(sigma_arr, edges_pct)
    band_specs = [
        (0.125, sigma_arr <= edges[0]),
        (0.375, (sigma_arr > edges[0]) & (sigma_arr <= edges[1])),
        (0.625, (sigma_arr > edges[1]) & (sigma_arr <= edges[2])),
        (0.825, (sigma_arr > edges[2]) & (sigma_arr <= edges[3])),  # 75-90%
        (0.925, (sigma_arr > edges[3]) & (sigma_arr <= edges[4])),  # 90-95%
        (0.97,  (sigma_arr > edges[4]) & (sigma_arr <= edges[5])),  # 95-99%
        (0.995, sigma_arr > edges[5]),                               # 99-100%
    ]

    K_emp_all = []
    z_q_all = []
    n_band_all = []
    for q, mask in band_specs:
        K = K_eff_band(sigma_arr[mask].astype(np.float64),
                        s_phys[mask], raw_thresh[mask])
        K_emp_all.append(K)
        z_q_all.append(float(norm.ppf(q)))
        n_band_all.append(int(mask.sum()))
    K_emp = np.array(K_emp_all)
    z_q = np.array(z_q_all)
    q_arr = np.array([q for q, _ in band_specs])

    print(f"  band K_eff:")
    print(f"    {'q':>6} {'z_q':>7} {'K_emp':>9} {'n_band':>9}")
    for i in range(len(q_arr)):
        print(f"    {q_arr[i]:>6.3f} {z_q[i]:>+7.3f} {K_emp[i]:>9.4f} {n_band_all[i]:>9,}")

    # Body fit (q ∈ {0.375, 0.625}): two-point line through body interior
    # (q=0.125 is in left tail, possibly skewed; use the "body-body" pair only)
    body_mask = (q_arr >= 0.30) & (q_arr <= 0.70)
    z_b = z_q[body_mask]
    K_b = K_emp[body_mask]
    if len(z_b) >= 2:
        z_c = z_b - z_b.mean()
        K_c = K_b - K_b.mean()
        b_fit = float((z_c * K_c).sum() / (z_c * z_c).sum())
        K_h_fit = float(K_b.mean() - b_fit * z_b.mean())
    else:
        b_fit, K_h_fit = B_USER, K_H_USER

    print(f"\n  Body fit (q ∈ [0.30, 0.70]): K_h_fit = {K_h_fit:.4f}, b_fit = {b_fit:.4f}")
    print(f"  User reference:                K_h     = {K_H_USER:.4f}, b     = {B_USER:.4f}")

    # X-frame translation (USER): X_q^user = (K_emp - K_H_USER)/B_USER
    X_user = (K_emp - K_H_USER) / B_USER
    # X-frame translation (BODY): X_q^body = (K_emp - K_h_fit)/b_fit
    X_body = (K_emp - K_h_fit) / b_fit

    # Right-tail fit for both X-frames
    q_thresh = 0.75
    X_thresh = norm.ppf(q_thresh)  # 0.674

    # Tail q's: {0.825, 0.925, 0.97, 0.995}
    tail_mask = q_arr >= q_thresh
    q_tail = q_arr[tail_mask]
    excess_user = X_user[tail_mask] - X_thresh
    excess_body = X_body[tail_mask] - X_thresh

    print(f"\n  Tail q values: {q_tail}")
    print(f"  Excess in X^user (with K_h={K_H_USER:.3f}, b={B_USER:.3f}):")
    print(f"    {[f'{v:+.4f}' for v in excess_user]}")
    print(f"  Excess in X^body (with K_h={K_h_fit:.3f}, b={b_fit:.3f}):")
    print(f"    {[f'{v:+.4f}' for v in excess_body]}")

    # Multi-point GPD fit (NLS)
    xi_user, sigma_user, sse_user = fit_xi_least_squares(q_thresh, q_tail, excess_user)
    xi_body, sigma_body, sse_body = fit_xi_least_squares(q_thresh, q_tail, excess_body)
    print(f"\n  Multi-point GPD fit on {len(q_tail)} tail points:")
    print(f"    USER frame: ξ = {xi_user:+.4f}, σ_GPD = {sigma_user:.4f}, SSE = {sse_user:.4e}")
    print(f"    BODY frame: ξ = {xi_body:+.4f}, σ_GPD = {sigma_body:.4f}, SSE = {sse_body:.4e}")
    print(f"    σ-records reference: ξ_σrec ≈ {GPD_XI_REF}")
    print(f"    USER-frame ξ vs σrec ξ: gap = {xi_user - GPD_XI_REF:+.4f}")
    print(f"    BODY-frame ξ vs σrec ξ: gap = {xi_body - GPD_XI_REF:+.4f}")

    # Predict from each fit and compare
    K_pred_user = np.zeros_like(K_emp)
    K_pred_body = np.zeros_like(K_emp)
    for i, q in enumerate(q_arr):
        if q < q_thresh:
            K_pred_user[i] = K_H_USER + B_USER * z_q[i]
            K_pred_body[i] = K_h_fit + b_fit * z_q[i]
        else:
            p = (q - q_thresh) / (1.0 - q_thresh)
            X_user_pred = X_thresh + sigma_user * gpd_factor(p, xi_user)
            X_body_pred = X_thresh + sigma_body * gpd_factor(p, xi_body)
            K_pred_user[i] = K_H_USER + B_USER * X_user_pred
            K_pred_body[i] = K_h_fit + b_fit * X_body_pred

    res_user = K_emp - K_pred_user
    res_body = K_emp - K_pred_body
    print(f"\n  Predictions:")
    print(f"  {'q':>6} {'K_emp':>9} {'K_user':>9} {'res_user':>9} {'K_body':>9} {'res_body':>9}")
    for i in range(len(q_arr)):
        print(f"  {q_arr[i]:>6.3f} {K_emp[i]:>9.4f} {K_pred_user[i]:>9.4f} {res_user[i]:>+9.4f} "
              f"{K_pred_body[i]:>9.4f} {res_body[i]:>+9.4f}")
    print(f"\n  Max |res_user| = {np.max(np.abs(res_user)):.3f}, "
          f"max |res_body| = {np.max(np.abs(res_body)):.3f}")

    return {
        'N': N_val,
        'log2N': int(np.log2(N_val)),
        'K_h_fit': K_h_fit, 'b_fit': b_fit,
        'xi_user': xi_user, 'sigma_user': sigma_user,
        'xi_body': xi_body, 'sigma_body': sigma_body,
        'q_arr': q_arr, 'K_emp': K_emp,
        'K_pred_user': K_pred_user, 'K_pred_body': K_pred_body,
        'res_user': res_user, 'res_body': res_body,
        'tail_q': q_tail, 'excess_user': excess_user, 'excess_body': excess_body,
    }


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(exist_ok=True)

    rng = np.random.default_rng(42)
    Ns = [1 << 24, 1 << 25, 1 << 27]
    results = []
    for N in Ns:
        r = analyze_N(N, data_dir, n_sample=500_000, rng=rng)
        results.append(r)

    # N-invariance check
    print(f"\n\n========== N-invariance of fitted ξ ==========")
    print(f"  {'N':>10}  {'K_h_fit':>9}  {'b_fit':>8}  "
          f"{'ξ_user':>9}  {'ξ_body':>9}")
    for r in results:
        print(f"  2^{r['log2N']:<8}  {r['K_h_fit']:>9.4f}  {r['b_fit']:>8.4f}  "
              f"{r['xi_user']:>+9.4f}  {r['xi_body']:>+9.4f}")
    xi_body_arr = np.array([r['xi_body'] for r in results])
    xi_user_arr = np.array([r['xi_user'] for r in results])
    print(f"\n  ξ_body spread: {np.ptp(xi_body_arr):.4f}, mean = {xi_body_arr.mean():+.4f}")
    print(f"  ξ_user spread: {np.ptp(xi_user_arr):.4f}, mean = {xi_user_arr.mean():+.4f}")
    print(f"  σ-records ξ:   {GPD_XI_REF}")

    # Save CSV
    rows = []
    for r in results:
        for i, q in enumerate(r['q_arr']):
            rows.append({
                'N': r['N'], 'log2N': r['log2N'],
                'q': float(q),
                'K_emp': float(r['K_emp'][i]),
                'K_pred_user': float(r['K_pred_user'][i]),
                'K_pred_body': float(r['K_pred_body'][i]),
                'res_user': float(r['res_user'][i]),
                'res_body': float(r['res_body'][i]),
                'K_h_fit': r['K_h_fit'], 'b_fit': r['b_fit'],
                'xi_user': r['xi_user'], 'xi_body': r['xi_body'],
            })
    out_csv = out_dir / "48_piecewise_followup.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
