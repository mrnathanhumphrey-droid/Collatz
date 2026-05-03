"""
Experiment 41 — MSE / scale-dependent structure of σ residuals along orbits.

Tests whether the slope K_h − K_eff = 0.486 on log(threshold) lives in the
variance / MSE structure of σ residuals, binned by orbit-value scale.

Hypothesis (user, 2026-05-02): two prior mechanism candidates ruled out
(trajectory-weighted σ-vs-log-v slope and E[v]=1.995 correction). Try
MSE-by-scale: if var(resid_partial(v)) is scale-dependent in the right way,
its integral against the trajectory-visit measure could reproduce 0.486.

Method (per prompt at experiments/41_mse_scale_residuals_PROMPT.md):
  1. Sample ≥200K odd starts from N=2^25 parquet.
  2. Compute α_det(r) for r = n mod 2^k (k=8). Sanity-check the
     start-point residual resid(n) = σ(n) − α_det(r) − K_h·log(n).
  3. Walk each orbit; at each visit (v, steps_remaining), record
     resid_partial(v) = steps_remaining(v) − K_h·log(v) into 20 log-v bins.
  4. Per-bin mean and var. Slopes above/below √N. Integral of d(var)/d(log v)
     against trajectory weights, compared to 0.486.

Decisive criteria:
  - mean has clean kink at √N with slopes K_h above and K_eff below
    → two-regime descent identified.
  - integral of d(var)/d(log v) · weight = 0.486 ± 5%
    → MSE-by-scale identified.
  - else → ruled out, look elsewhere.

Usage:
    python 41_mse_scale_residuals.py
"""
import gc
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

K_H = 3.0 / (np.log(4.0) - np.log(3.0))
LOG3 = np.log(3.0)
LOG2 = np.log(2.0)


def compute_alpha_det(k):
    """For each odd residue r mod 2^k, compute α_det(r) =
    prefix_steps(r) + K_h · (j(r)·log 3 − k·log 2),
    where j(r) is the number of odd Collatz steps in the first k T-map
    applications starting from a representative n = 2^k + r.
    Returns alpha_det[0..2^k − 1] with even indices NaN."""
    M = 1 << k
    alpha = np.full(M, np.nan, dtype=np.float64)
    j_arr = np.zeros(M, dtype=np.int32)
    pre_arr = np.zeros(M, dtype=np.int32)
    for r in range(1, M, 2):
        n = M + r
        odd_count = 0
        steps = 0
        for _ in range(k):
            if n & 1:
                n = (3 * n + 1) // 2
                odd_count += 1
                steps += 2
            else:
                n = n // 2
                steps += 1
        j_arr[r] = odd_count
        pre_arr[r] = steps
        alpha[r] = steps + K_H * (odd_count * LOG3 - k * LOG2)
    return alpha, j_arr, pre_arr


@njit(parallel=True, cache=True)
def walk_and_bin(starts, sigmas, log_v_min, log_v_max, n_bins):
    """For each odd start, walk Collatz orbit. At every value v on the orbit,
    accumulate resid_partial(v) = (σ_total − steps_done) − K_h·log(v)
    into log-v bins.

    Returns (counts, sum_r, sum_r2) per bin."""
    n = len(starts)
    chunk = max(1, n // 256)
    n_chunks = (n + chunk - 1) // chunk
    counts_local = np.zeros((n_chunks, n_bins), dtype=np.int64)
    sum_r_local = np.zeros((n_chunks, n_bins), dtype=np.float64)
    sum_r2_local = np.zeros((n_chunks, n_bins), dtype=np.float64)
    K_h = 3.0 / (np.log(4.0) - np.log(3.0))
    bin_width = (log_v_max - log_v_min) / n_bins
    max_steps = 20000

    for c in prange(n_chunks):
        start_i = c * chunk
        end_i = min(start_i + chunk, n)
        for i in range(start_i, end_i):
            m = starts[i]
            sigma_n = sigmas[i]
            steps_done = 0
            while m > 1 and steps_done < max_steps:
                log_v = np.log(np.float64(m))
                if log_v >= log_v_min and log_v < log_v_max:
                    steps_remaining = sigma_n - steps_done
                    resid = float(steps_remaining) - K_h * log_v
                    bin_idx = int((log_v - log_v_min) / bin_width)
                    if 0 <= bin_idx < n_bins:
                        counts_local[c, bin_idx] += 1
                        sum_r_local[c, bin_idx] += resid
                        sum_r2_local[c, bin_idx] += resid * resid
                if m & 1:
                    m = 3 * m + 1
                else:
                    m = m >> 1
                steps_done += 1
            # final visit at m=1: log(v) = 0
            if m == 1:
                log_v = 0.0
                if log_v >= log_v_min and log_v < log_v_max:
                    steps_remaining = sigma_n - steps_done
                    resid = float(steps_remaining) - K_h * log_v
                    bin_idx = int((log_v - log_v_min) / bin_width)
                    if 0 <= bin_idx < n_bins:
                        counts_local[c, bin_idx] += 1
                        sum_r_local[c, bin_idx] += resid
                        sum_r2_local[c, bin_idx] += resid * resid

    counts = counts_local.sum(axis=0)
    sum_r = sum_r_local.sum(axis=0)
    sum_r2 = sum_r2_local.sum(axis=0)
    return counts, sum_r, sum_r2


def fit_slope(x, y):
    m = ~np.isnan(y) & np.isfinite(y)
    if m.sum() < 3:
        return float("nan"), float("nan")
    x_, y_ = x[m], y[m]
    xm, ym = x_.mean(), y_.mean()
    num = ((x_ - xm) * (y_ - ym)).sum()
    den = ((x_ - xm) ** 2).sum()
    if den < 1e-20:
        return float("nan"), float("nan")
    slope = num / den
    intercept = ym - slope * xm
    return float(slope), float(intercept)


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "outputs" / "exp41"
    out_dir.mkdir(parents=True, exist_ok=True)

    N = 1 << 25
    print(f"[load] N = 2^25 = {N:,}", flush=True)
    df = pl.read_parquet(data_dir / f"main_N{N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    sigma_arr = df["sigma"].to_numpy().astype(np.int64)
    print(f"  {len(n_arr):,} odd rows loaded", flush=True)

    rng = np.random.default_rng(42)
    n_samples = 200_000
    idx = rng.choice(len(n_arr), size=n_samples, replace=False)
    starts = n_arr[idx]
    sigmas = sigma_arr[idx]
    print(f"  sampled {n_samples:,} odd starts (seed=42)", flush=True)

    # ---- Step 1: α_det at start, sanity-check resid(n) ----
    k_mod = 8
    print(f"\n[step 1] α_det(r) for k = {k_mod}", flush=True)
    alpha_det, j_arr, pre_arr = compute_alpha_det(k_mod)
    M = 1 << k_mod
    r_starts = (starts & (M - 1)).astype(np.int64)
    alpha_at_starts = alpha_det[r_starts]
    log_starts = np.log(starts.astype(np.float64))
    resid_n = sigmas.astype(np.float64) - alpha_at_starts - K_H * log_starts
    print(f"  α_det range: [{np.nanmin(alpha_det[1::2]):.3f}, {np.nanmax(alpha_det[1::2]):.3f}]")
    print(f"  resid(n) at start: mean = {resid_n.mean():.4f}, "
          f"std = {resid_n.std():.4f}, "
          f"median = {np.median(resid_n):.4f}")
    # Same residual without the α_det adjustment (pure resid_partial(v=n))
    resid_n_noalpha = sigmas.astype(np.float64) - K_H * log_starts
    print(f"  resid(n) no-α (= resid_partial at start): "
          f"mean = {resid_n_noalpha.mean():.4f}, "
          f"std = {resid_n_noalpha.std():.4f}")

    # ---- Step 2 & 3: walk orbits, bin by log(v) ----
    log_v_min = 0.0
    log_v_max = np.log(float(N))
    n_bins = 20
    bin_width = (log_v_max - log_v_min) / n_bins
    bin_lo = log_v_min + np.arange(n_bins) * bin_width
    bin_hi = bin_lo + bin_width
    bin_centers = bin_lo + bin_width / 2

    print(f"\n[step 2+3] walking {n_samples:,} orbits; binning resid_partial(v) "
          f"into {n_bins} log-v bins from {log_v_min:.2f} to {log_v_max:.2f}",
          flush=True)
    t0 = time.perf_counter()
    counts, sum_r, sum_r2 = walk_and_bin(
        starts, sigmas, log_v_min, log_v_max, n_bins)
    elapsed = time.perf_counter() - t0
    total_visits = int(counts.sum())
    print(f"  done in {elapsed:.1f}s   total visits binned: {total_visits:,}",
          flush=True)

    valid = counts > 100
    mean_r = np.full(n_bins, np.nan)
    var_r = np.full(n_bins, np.nan)
    mean_r[valid] = sum_r[valid] / counts[valid]
    var_r[valid] = sum_r2[valid] / counts[valid] - mean_r[valid] ** 2

    sqrt_N = np.sqrt(float(N))
    log_sqrt_N = np.log(sqrt_N)

    print(f"\n  log(√N) = {log_sqrt_N:.4f}   (between bins "
          f"{int((log_sqrt_N - log_v_min) / bin_width)} and "
          f"{int((log_sqrt_N - log_v_min) / bin_width) + 1})")
    print(f"  Bin   log(v)_lo   log(v)_hi          count  "
          f"mean(resid)   var(resid)")
    for b in range(n_bins):
        marker = "  <-- sqrt(N)" if bin_lo[b] <= log_sqrt_N < bin_hi[b] else ""
        m_str = f"{mean_r[b]:>10.3f}" if valid[b] else "       N/A"
        v_str = f"{var_r[b]:>10.3f}" if valid[b] else "       N/A"
        print(f"  {b:>3}   {bin_lo[b]:>9.3f}   {bin_hi[b]:>9.3f}   "
              f"{counts[b]:>12,d}   {m_str}   {v_str}{marker}")

    # ---- Slope fits: above / below √N ----
    above = bin_centers > log_sqrt_N
    below = ~above
    s_above, i_above = fit_slope(bin_centers[above], mean_r[above])
    s_below, i_below = fit_slope(bin_centers[below], mean_r[below])
    s_full, i_full = fit_slope(bin_centers, mean_r)

    print(f"\n=== Slope of mean(resid_partial) vs log(v) ===")
    print(f"  full-range slope = {s_full:>+.4f}    intercept = {i_full:>+.4f}")
    print(f"  above sqrt(N):     slope = {s_above:>+.4f}    intercept = {i_above:>+.4f}")
    print(f"  below sqrt(N):     slope = {s_below:>+.4f}    intercept = {i_below:>+.4f}")
    print(f"  K_h               = {K_H:.4f}    (expected slope above if 2-regime)")
    print(f"  K_eff (calibrated) = {K_H - 0.486:.4f}  (expected slope below if 2-regime)")
    print(f"  slope_above - slope_below = {s_above - s_below:>+.4f}")
    print(f"  target K_h - K_eff = +0.486")

    # ---- MSE integral test ----
    print(f"\n=== MSE-by-scale integral test ===")
    dvar_dlog = np.full(n_bins, np.nan)
    for b in range(n_bins - 1):
        if valid[b] and valid[b + 1]:
            dvar_dlog[b] = (var_r[b + 1] - var_r[b]) / bin_width

    print(f"  Per-bin d(var)/d(log v):")
    for b in range(n_bins - 1):
        if not np.isnan(dvar_dlog[b]):
            print(f"    bin {b:>2} -> {b+1:>2}: "
                  f"log(v) ~ {bin_centers[b]:>5.2f} : "
                  f"d(var)/d(log v) = {dvar_dlog[b]:>+9.3f} "
                  f"(weight = {counts[b]:>10,d})")

    weights = counts.astype(np.float64)
    fwd_idx = np.arange(n_bins - 1)
    valid_int = ~np.isnan(dvar_dlog[:n_bins - 1])
    if valid_int.sum() > 0:
        w = weights[fwd_idx][valid_int]
        d = dvar_dlog[fwd_idx][valid_int]
        if w.sum() > 0:
            integral = (w * d).sum() / w.sum()
            print(f"\n  Trajectory-weighted integral:")
            print(f"    integral d(var)/d(log v) * weight(v) dv / integral weight(v) dv = {integral:>+8.4f}")
            print(f"    target K_h - K_eff = +0.486")
            diff = abs(integral - 0.486)
            rel = diff / 0.486
            if rel < 0.05:
                verdict = "MSE-by-scale identified (within 5%)"
            elif rel < 0.20:
                verdict = "Partial mechanism (within 20%)"
            else:
                verdict = "Ruled out (>20% off)"
            print(f"    relative error = {rel*100:.1f}%   verdict: {verdict}")

    # Variance ratio above/below √N (scale-dependence indicator)
    var_above_mean = np.nanmean(var_r[above])
    var_below_mean = np.nanmean(var_r[below])
    print(f"\n  mean(var) above sqrt(N) = {var_above_mean:.3f}")
    print(f"  mean(var) below sqrt(N) = {var_below_mean:.3f}")
    print(f"  ratio (below / above)  = {var_below_mean / var_above_mean:.3f}")

    # ---- Plots ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        print(f"\n[plots] writing to {out_dir}", flush=True)

        # Plot 1: mean(resid) vs log(v)
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(bin_centers, mean_r, "o-", label="mean(resid_partial)")
        ax.axvline(log_sqrt_N, color="red", linestyle="--", alpha=0.5,
                   label=f"log(sqrt(N)) = {log_sqrt_N:.2f}")
        # K_h slope reference (anchored at right edge)
        anchor_x = log_v_max
        anchor_y = mean_r[valid][-1] if valid.any() else 0.0
        ref_x = np.array([log_v_min, log_v_max])
        ref_y_kh = anchor_y + K_H * (ref_x - anchor_x)
        ax.plot(ref_x, ref_y_kh, "g:", alpha=0.6,
                label=f"slope K_h = {K_H:.2f} (reference)")
        ax.set_xlabel("log(v)")
        ax.set_ylabel("mean(steps_remaining - K_h * log v)")
        ax.set_title("Plot 1: mean(resid_partial) vs log(v)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / "plot1_mean_resid.png", dpi=120)
        plt.close()

        # Plot 2: var(resid) vs log(v)
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(bin_centers, var_r, "o-")
        ax.axvline(log_sqrt_N, color="red", linestyle="--", alpha=0.5,
                   label=f"log(sqrt(N)) = {log_sqrt_N:.2f}")
        ax.set_xlabel("log(v)")
        ax.set_ylabel("var(steps_remaining - K_h * log v)")
        ax.set_title("Plot 2: var(resid_partial) vs log(v)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / "plot2_var_resid.png", dpi=120)
        plt.close()

        # Plot 3: cumulative var (running sum from low log v upward)
        cum_var = np.full(n_bins, np.nan)
        running = 0.0
        for b in range(n_bins):
            if valid[b]:
                running += var_r[b]
                cum_var[b] = running
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(bin_centers, cum_var, "o-")
        ax.axvline(log_sqrt_N, color="red", linestyle="--", alpha=0.5,
                   label=f"log(sqrt(N)) = {log_sqrt_N:.2f}")
        ax.set_xlabel("log(v)")
        ax.set_ylabel("cumulative var(resid_partial)")
        ax.set_title("Plot 3: running sum of var(resid_partial) vs log(v)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / "plot3_cum_var.png", dpi=120)
        plt.close()
        print(f"  3 plots saved to {out_dir}", flush=True)
    except Exception as e:
        print(f"  [warn] plotting failed: {e}", flush=True)

    print(f"\n[done]")


if __name__ == "__main__":
    main()
