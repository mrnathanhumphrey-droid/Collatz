"""
Experiment 47 — piecewise body-Gaussian + tail-GPD test for K_eff_band(q).

Hypothesis: standardized residual X = (σ - K_h·log n)/scale has Gaussian body
(q ∈ [0.25, 0.75]) and GPD right tail (q ≥ 0.75) with shape ξ = -0.075
(k-invariant from σ-records work).

Model:
  K_eff_band(q) = K_h + b·X_q
where K_h = 10.4282 (asymptotic), b = 2.275 (σ-std-vs-log(n) slope).

Body (q ∈ [0.25, 0.75]):  X_q = z_q (standard Gaussian quantile)
Right tail (q ≥ 0.75):    X_q = X_thresh + F⁻¹_GPD(q'; ξ, σ_GPD)
                          q' = (q - 0.75)/0.25, X_thresh = z_0.75 = 0.674
                          F⁻¹_GPD(p; ξ, σ_GPD) = σ_GPD · ((1-p)^(-ξ) - 1)/ξ

Test data: K_eff_band at q ∈ {0.125, 0.375, 0.625, 0.875, 0.975}.
First 4 from existing exp 45 (quartile midpoints); 5th computed here as
top-5% band on N=2^27 first-passage R(f) regression.

Calibration: σ_GPD matched at q=0.875 (one tail data point); validate at q=0.975.

Verdict criterion (from prompt):
  - residual ≤ 0.5 at all tail q  → piecewise closed form CONFIRMED
  - residual ≤ 1.0 partial valid  → boundary not sharp; partial mechanism
  - residual > 1.5                → tail isn't GPD with ξ=-0.075 (falsified)
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from scipy.stats import norm

sys.stdout.reconfigure(encoding="utf-8")

K_H = 3.0 / (np.log(4.0) - np.log(3.0))
B_SLOPE = 2.275
GPD_XI = -0.075
Q_THRESH = 0.75


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


def build_thresholds_and_logf(n_arr):
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
    """Slope of (sigma_bar - s_bar(f)) on log(f) using 4 thresholds, OLS."""
    log_f = np.log(np.maximum(thresh_phys, 1.0))  # n x 4
    R = sigma[:, None] - s_per_thresh.astype(np.float64)  # n x 4
    log_f_mean = log_f.mean(axis=0)  # length 4
    R_mean = R.mean(axis=0)
    x = log_f_mean
    y = R_mean
    x_c = x - x.mean()
    y_c = y - y.mean()
    slope = (x_c * y_c).sum() / (x_c * x_c).sum()
    return float(slope)


def gpd_inverse(p, xi, sigma):
    """F⁻¹_GPD(p; ξ, σ) for ξ ≠ 0."""
    return sigma * ((1.0 - p) ** (-xi) - 1.0) / xi


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "experiments_output"

    N_val = 1 << 27
    print(f"# N = {N_val:,} = 2^27")
    t0 = time.perf_counter()
    df = pl.read_parquet(data_dir / f"main_N{N_val}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    sigma_arr = df["sigma"].to_numpy().astype(np.int32)
    print(f"  loaded {len(n_arr):,} odd rows in {time.perf_counter()-t0:.1f}s")

    # Subsample for speed
    rng = np.random.default_rng(42)
    n_sample = 500_000
    if len(n_arr) > n_sample:
        idx = rng.choice(len(n_arr), size=n_sample, replace=False)
        n_arr = n_arr[idx]; sigma_arr = sigma_arr[idx]
        print(f"  subsampled to {n_sample:,}")

    print(f"  building thresholds & walking ...")
    t0 = time.perf_counter()
    thresh_sorted, sort_idx, raw_thresh = build_thresholds_and_logf(n_arr)
    s_arr = first_passage_4thresh(n_arr, thresh_sorted)
    inv_sort = np.argsort(sort_idx, axis=1)
    s_phys = np.zeros((len(n_arr), 4), dtype=np.int32)
    for col in range(4):
        rev = inv_sort[:, col]
        s_phys[:, col] = np.take_along_axis(s_arr, rev[:, None], axis=1).flatten()
    print(f"    walked in {time.perf_counter()-t0:.1f}s")

    # σ-quantile bands at user's 5 quantile midpoints
    q_edges = np.percentile(sigma_arr, [25, 50, 75, 95])
    print(f"\n  σ quantile edges: q25={q_edges[0]:.0f}, q50={q_edges[1]:.0f}, "
          f"q75={q_edges[2]:.0f}, q95={q_edges[3]:.0f}")

    bands = [
        (0.125, sigma_arr <= q_edges[0]),
        (0.375, (sigma_arr > q_edges[0]) & (sigma_arr <= q_edges[1])),
        (0.625, (sigma_arr > q_edges[1]) & (sigma_arr <= q_edges[2])),
        (0.875, (sigma_arr > q_edges[2]) & (sigma_arr <= q_edges[3])),
        (0.975, sigma_arr > q_edges[3]),
    ]

    K_emp = []
    z_q_list = []
    for q, mask in bands:
        n_band = int(mask.sum())
        K = K_eff_band(sigma_arr[mask].astype(np.float64),
                       s_phys[mask], raw_thresh[mask])
        K_emp.append(K)
        z_q_list.append(float(norm.ppf(q)))
        sigma_mean = sigma_arr[mask].mean()
        print(f"    q={q:.3f}  z_q={z_q_list[-1]:+.3f}  K_emp={K:>8.4f}  "
              f"σ̄={sigma_mean:>7.2f}  n={n_band:,}")

    K_emp = np.array(K_emp)
    z_q_arr = np.array(z_q_list)
    q_arr = np.array([q for q, _ in bands])

    # Gaussian prediction
    K_gauss = K_H + B_SLOPE * z_q_arr

    # Piecewise: Gaussian body + GPD right tail (with calibration at q=0.875)
    X_thresh = norm.ppf(Q_THRESH)  # ≈ 0.674
    # Calibrate σ_GPD at q=0.875: empirical excess X = (K_emp[3] - K_h)/b - X_thresh
    X_at_0875 = (K_emp[3] - K_H) / B_SLOPE
    excess_emp = X_at_0875 - X_thresh
    q_prime_0875 = (0.875 - Q_THRESH) / (1.0 - Q_THRESH)  # = 0.5
    gpd_factor_0875 = ((1.0 - q_prime_0875) ** (-GPD_XI) - 1.0) / GPD_XI
    sigma_GPD = excess_emp / gpd_factor_0875

    print(f"\n  Calibration:")
    print(f"    X_thresh = z_0.75 = {X_thresh:.4f}")
    print(f"    Empirical X at q=0.875: ({K_emp[3]:.3f} - {K_H:.3f})/{B_SLOPE} = {X_at_0875:.4f}")
    print(f"    Excess at q=0.875: {excess_emp:.4f}")
    print(f"    GPD inverse factor at q'=0.5: {gpd_factor_0875:.4f}")
    print(f"    Calibrated σ_GPD = {sigma_GPD:.4f}")

    K_piecewise = np.zeros(5)
    for i, q in enumerate(q_arr):
        if q < Q_THRESH:
            # Body: Gaussian
            K_piecewise[i] = K_H + B_SLOPE * z_q_arr[i]
        else:
            # Tail: Gaussian threshold + GPD excess
            q_prime = (q - Q_THRESH) / (1.0 - Q_THRESH)
            excess = sigma_GPD * gpd_inverse(q_prime, GPD_XI, 1.0)  # (already in sigma_GPD units)
            X_q = X_thresh + excess
            K_piecewise[i] = K_H + B_SLOPE * X_q

    # Residuals
    resid_gauss = K_emp - K_gauss
    resid_piecewise = K_emp - K_piecewise

    print(f"\n  ===== Predictions vs empirical (N=2^27) =====")
    print(f"  {'q':>6} {'z_q':>7} {'K_emp':>9} {'K_gauss':>9} {'K_piecew':>9} "
          f"{'res_G':>8} {'res_PW':>8}")
    for i, q in enumerate(q_arr):
        print(f"  {q:>6.3f} {z_q_arr[i]:>+7.3f} {K_emp[i]:>9.4f} {K_gauss[i]:>9.4f} "
              f"{K_piecewise[i]:>9.4f} {resid_gauss[i]:>+8.3f} {resid_piecewise[i]:>+8.3f}")

    print(f"\n  Max |residual| Gaussian:   {np.max(np.abs(resid_gauss)):.3f}")
    print(f"  Max |residual| piecewise:  {np.max(np.abs(resid_piecewise)):.3f}")
    print(f"  Tail residuals (q ∈ {{0.875, 0.975}}):")
    print(f"    Gaussian:   {resid_gauss[3]:+.3f}, {resid_gauss[4]:+.3f}")
    print(f"    Piecewise:  {resid_piecewise[3]:+.3f}, {resid_piecewise[4]:+.3f}")

    # Verdict
    max_tail_pw = max(abs(resid_piecewise[3]), abs(resid_piecewise[4]))
    print(f"\n  ===== Verdict =====")
    if max_tail_pw <= 0.5:
        verdict = "CONFIRMED — piecewise body+GPD-tail closed form within ±0.5 at tails."
    elif max_tail_pw <= 1.0:
        verdict = "PARTIAL — boundary not sharp at q=0.75; piecewise improves over Gaussian but residual structure remains."
    else:
        verdict = f"FALSIFIED — piecewise residual {max_tail_pw:.3f} > 1.5 at tails. Tail isn't GPD with ξ=-0.075 OR b/K_h calibration is wrong for this regime."
    print(f"  {verdict}")

    rows = []
    for i, q in enumerate(q_arr):
        rows.append({
            'N': N_val,
            'q': float(q),
            'z_q': float(z_q_arr[i]),
            'K_eff_emp': float(K_emp[i]),
            'K_eff_gauss': float(K_gauss[i]),
            'K_eff_piecewise': float(K_piecewise[i]),
            'resid_gauss': float(resid_gauss[i]),
            'resid_piecewise': float(resid_piecewise[i]),
            'sigma_GPD_calibrated': float(sigma_GPD),
            'xi_GPD': float(GPD_XI),
        })
    out_csv = out_dir / "47_piecewise_body_tail.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")

    # Solve for exact ξ that fits both tail points (q=0.875 calibration + q=0.975 validation)
    excess_0875 = (K_emp[3] - K_H) / B_SLOPE - X_thresh
    excess_0975 = (K_emp[4] - K_H) / B_SLOPE - X_thresh
    target_ratio = excess_0975 / excess_0875

    def g(p, xi):
        if abs(xi) < 1e-9:
            return -np.log(1.0 - p)
        return ((1.0 - p) ** (-xi) - 1.0) / xi

    def ratio_at_xi(xi):
        return g(0.9, xi) / g(0.5, xi)

    from scipy.optimize import brentq
    xi_fit = brentq(lambda x: ratio_at_xi(x) - target_ratio, -0.3, 0.3)
    sigma_GPD_fit = excess_0875 / g(0.5, xi_fit)
    print(f"\n  ===== Two-point GPD fit (q=0.875 + q=0.975) =====")
    print(f"    excess at q=0.875: {excess_0875:.4f}")
    print(f"    excess at q=0.975: {excess_0975:.4f}")
    print(f"    ratio = {target_ratio:.4f}")
    print(f"    fitted ξ = {xi_fit:+.4f}")
    print(f"    fitted σ_GPD = {sigma_GPD_fit:.4f}")
    print(f"  σ-records reference: ξ_σrec ≈ -0.075 (Weibull domain, bounded)")
    print(f"  X residual fitted ξ = {xi_fit:+.4f}: "
          f"{'positive (Pareto, heavy tail)' if xi_fit > 0.01 else 'near zero (exponential)' if abs(xi_fit) < 0.01 else 'negative (Weibull, bounded)'}")

    # Plot
    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                                  gridspec_kw={'height_ratios': [2, 1]})
        axes[0].plot(q_arr, K_emp, 'ko-', label='empirical', markersize=8)
        axes[0].plot(q_arr, K_gauss, 'b--', label='Gaussian (K_h + b·z_q)', alpha=0.7)
        axes[0].plot(q_arr, K_piecewise, 'r-.',
                     label=f'piecewise (Gauss body + GPD tail ξ={GPD_XI})', alpha=0.7)
        axes[0].axhline(K_H, color='gray', linestyle=':', alpha=0.5, label=f'K_h={K_H:.3f}')
        axes[0].axvline(Q_THRESH, color='gray', linestyle=':', alpha=0.3)
        axes[0].set_ylabel('K_eff_band(q)')
        axes[0].set_title(f'Piecewise body+tail closed form for K_eff_band at N=2^27')
        axes[0].legend(fontsize=9)
        axes[0].grid(alpha=0.3)
        axes[1].plot(q_arr, resid_gauss, 'b--o', label='Gaussian residual', alpha=0.7)
        axes[1].plot(q_arr, resid_piecewise, 'r-.s', label='piecewise residual', alpha=0.7)
        axes[1].axhline(0, color='black', linewidth=0.7)
        axes[1].axhline(0.5, color='green', linestyle=':', alpha=0.5, label='±0.5 success')
        axes[1].axhline(-0.5, color='green', linestyle=':', alpha=0.5)
        axes[1].set_xlabel('quantile q')
        axes[1].set_ylabel('residual')
        axes[1].legend(fontsize=9)
        axes[1].grid(alpha=0.3)
        plt.tight_layout()
        plot_path = out_dir / "47_piecewise_body_tail.png"
        plt.savefig(plot_path, dpi=120)
        print(f"[save] {plot_path}")
    except Exception as e:
        print(f"  (plot skipped: {e})")


if __name__ == "__main__":
    main()
