"""
Experiment 45 — first-passage K_eff with σ-quantile-band stratification.

Re-runs the finite-N boundary decomposition from exp 44 using Result 6's
K_eff definition (slope of R(f) = σ̄ − s̄(f) on log(f) from 4-threshold
first-passage), which has the genuine 1.35 spread across N. Adds σ-quantile
stratification to test whether boundary effects shift band-weighting.

Tests:

  Effect A (sample-truncation bias):
    K_eff_full     = first-passage K_eff on all odd n in [1, N].
    K_eff_interior = same but on n < N/8 (drops top 3 octaves).
    If K_int N-stable while K_full oscillates, Effect A confirmed.

  Effect B (orbit-mass concentration):
    K_eff_above = first-passage K_eff restricted to n > √N.
    K_eff_below = same restricted to n ≤ √N.
    If these differ, Effect B real.

  σ-band stratification:
    Within each N, partition starts by σ-quantile (q1=0-25%, q2=25-50%,
    q3=50-75%, q4=75-100%). Compute K_eff per band.
    If K_eff per band is N-stable but K_eff_full oscillates → boundary
    shifts band weights = mechanism via band reweighting.
    If K_eff per band oscillates with N → boundary not the mechanism.

Output: CSV (N, K_eff variants, per-band K_eff) + plots.

Usage:
    python 45_first_passage_K_eff_band_decomp.py
"""
import gc
import sys
import time
import csv
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

K_H = 3.0 / (np.log(4.0) - np.log(3.0))


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
    thresholds_sorted = np.take_along_axis(raw, sort_idx, axis=1).astype(np.int64)
    return raw, sort_idx, thresholds_sorted


def get_first_passage_data(n_arr, sigma_arr):
    """Walk all starts to first-passage of 4 thresholds. Returns (s_phys,
    thresh_phys, valid mask) where columns are [n^(2/3), √n·log n, √n,
    √n/log n] in physical (de-sorted) order."""
    raw, sort_idx, thresh_sorted = build_thresholds(n_arr)
    s_sorted = first_passage_4thresh(n_arr, thresh_sorted)
    inv_sort = np.argsort(sort_idx, axis=1)
    s_phys = np.zeros_like(s_sorted)
    thresh_phys = np.zeros_like(thresh_sorted)
    for col in range(4):
        rev = inv_sort[:, col]
        s_phys[:, col] = np.take_along_axis(s_sorted, rev[:, None], axis=1).flatten()
        thresh_phys[:, col] = np.take_along_axis(thresh_sorted, rev[:, None], axis=1).flatten()
    valid = np.all(s_phys >= 0, axis=1)
    return s_phys, thresh_phys, valid


def compute_K_eff(sigma_arr, s_phys, thresh_phys):
    """K_eff = slope of R(f) = σ̄ − s̄(f) on log(f) over the 4 thresholds."""
    if len(sigma_arr) < 10:
        return float("nan")
    sigma_mean = sigma_arr.mean()
    s_mean = s_phys.mean(axis=0).astype(np.float64)
    log_f_mean = np.log(thresh_phys.astype(np.float64)).mean(axis=0)
    R = sigma_mean - s_mean
    # OLS on 4 points
    xm, ym = log_f_mean.mean(), R.mean()
    den = ((log_f_mean - xm) ** 2).sum()
    if den < 1e-12:
        return float("nan")
    slope = ((log_f_mean - xm) * (R - ym)).sum() / den
    return float(slope)


def bootstrap_K_eff(sigma_arr, s_phys, thresh_phys, n_boot=1000, rng=None):
    rng = rng or np.random.default_rng(42)
    n = len(sigma_arr)
    if n < 10:
        return float("nan"), float("nan"), float("nan")
    boot_K = np.zeros(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        boot_K[b] = compute_K_eff(sigma_arr[idx], s_phys[idx], thresh_phys[idx])
    valid = ~np.isnan(boot_K)
    if valid.sum() < 10:
        return float("nan"), float("nan"), float("nan")
    point = compute_K_eff(sigma_arr, s_phys, thresh_phys)
    ci_lo, ci_hi = np.percentile(boot_K[valid], [2.5, 97.5])
    return float(point), float(ci_lo), float(ci_hi)


def analyze_one_N(data_dir, N, n_sample=500_000, rng=None):
    rng = rng or np.random.default_rng(42)
    df = pl.read_parquet(data_dir / f"main_N{N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    sigma_arr = df["sigma"].to_numpy().astype(np.int64)
    log2N = int(round(np.log2(N)))
    sqrt_N = int(np.sqrt(float(N)))
    cutoff_int = N // 8

    n_total = len(n_arr)
    if n_total > n_sample:
        idx = rng.choice(n_total, n_sample, replace=False)
        n_arr = n_arr[idx]
        sigma_arr = sigma_arr[idx]
    n_used = len(n_arr)

    print(f"\n=== N = 2^{log2N} = {N:,}  (using {n_used:,} of {n_total:,} odd rows) ===",
          flush=True)
    print(f"  cutoff_interior (n < N/8): {cutoff_int:,}")
    print(f"  cutoff_above/below (√N):   {sqrt_N:,}")

    t0 = time.perf_counter()
    s_phys, thresh_phys, valid = get_first_passage_data(n_arr, sigma_arr)
    n_arr = n_arr[valid]
    sigma_arr = sigma_arr[valid]
    s_phys = s_phys[valid]
    thresh_phys = thresh_phys[valid]
    print(f"  walked {len(n_arr):,} valid orbits in {time.perf_counter()-t0:.1f}s",
          flush=True)

    out = {"N": N, "log2N": log2N, "log_N": np.log(float(N)), "n_used": len(n_arr)}

    # K_eff_full
    pt, lo, hi = bootstrap_K_eff(sigma_arr, s_phys, thresh_phys, rng=rng)
    out["K_eff_full"] = pt
    out["K_eff_full_lo"] = lo
    out["K_eff_full_hi"] = hi
    out["K_eff_full_n"] = len(n_arr)
    print(f"  K_eff_full     = {pt:.4f}  CI [{lo:.4f}, {hi:.4f}]  (n={len(n_arr):,})")

    # K_eff_interior (n < N/8)
    mask = n_arr < cutoff_int
    pt, lo, hi = bootstrap_K_eff(sigma_arr[mask], s_phys[mask], thresh_phys[mask], rng=rng)
    out["K_eff_interior"] = pt
    out["K_eff_interior_lo"] = lo
    out["K_eff_interior_hi"] = hi
    out["K_eff_interior_n"] = int(mask.sum())
    print(f"  K_eff_interior = {pt:.4f}  CI [{lo:.4f}, {hi:.4f}]  (n={int(mask.sum()):,})")

    # K_eff_above (n > √N)
    mask = n_arr > sqrt_N
    pt, lo, hi = bootstrap_K_eff(sigma_arr[mask], s_phys[mask], thresh_phys[mask], rng=rng)
    out["K_eff_above"] = pt
    out["K_eff_above_lo"] = lo
    out["K_eff_above_hi"] = hi
    out["K_eff_above_n"] = int(mask.sum())
    print(f"  K_eff_above    = {pt:.4f}  CI [{lo:.4f}, {hi:.4f}]  (n={int(mask.sum()):,})")

    # K_eff_below (n ≤ √N)
    mask = n_arr <= sqrt_N
    pt, lo, hi = bootstrap_K_eff(sigma_arr[mask], s_phys[mask], thresh_phys[mask], rng=rng)
    out["K_eff_below"] = pt
    out["K_eff_below_lo"] = lo
    out["K_eff_below_hi"] = hi
    out["K_eff_below_n"] = int(mask.sum())
    print(f"  K_eff_below    = {pt:.4f}  CI [{lo:.4f}, {hi:.4f}]  (n={int(mask.sum()):,})")

    # σ-quantile bands
    q_edges = np.percentile(sigma_arr, [25, 50, 75])
    print(f"  σ quantile edges: q25={q_edges[0]:.0f}, q50={q_edges[1]:.0f}, q75={q_edges[2]:.0f}")
    print(f"  σ band K_eff:")
    bands = [
        ("q1_0_25", sigma_arr <= q_edges[0]),
        ("q2_25_50", (sigma_arr > q_edges[0]) & (sigma_arr <= q_edges[1])),
        ("q3_50_75", (sigma_arr > q_edges[1]) & (sigma_arr <= q_edges[2])),
        ("q4_75_100", sigma_arr > q_edges[2]),
    ]
    for name, mask in bands:
        n_band = int(mask.sum())
        sigma_band_mean = sigma_arr[mask].mean() if n_band > 0 else float("nan")
        pt, lo, hi = bootstrap_K_eff(sigma_arr[mask], s_phys[mask], thresh_phys[mask],
                                     rng=rng)
        out[f"K_eff_{name}"] = pt
        out[f"K_eff_{name}_lo"] = lo
        out[f"K_eff_{name}_hi"] = hi
        out[f"K_eff_{name}_n"] = n_band
        out[f"sigma_mean_{name}"] = sigma_band_mean
        out[f"weight_{name}"] = n_band / len(n_arr)
        print(f"    {name}: K_eff = {pt:.4f}  CI [{lo:.4f}, {hi:.4f}]  "
              f"(n={n_band:,}, σ̄_band={sigma_band_mean:.2f}, w={n_band/len(n_arr):.4f})")

    return out


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "outputs" / "exp45"
    out_dir.mkdir(parents=True, exist_ok=True)

    N_values = [
        1 << 20, 1 << 22, 1 << 23, 1 << 24, 1 << 25, 1 << 27,
    ]

    print(f"First-passage K_eff finite-N decomposition + σ-band stratification")
    print(f"K_h = {K_H:.6f}")
    print(f"K_eff defined as slope of R(f) = σ̄ − s̄(f) on log(f), 4 thresholds.")
    print(f"Sample 500K starts per N (or all if smaller).")
    print(f"Bootstrap: 1000 resamples.")

    rng = np.random.default_rng(42)
    results = []
    for N in N_values:
        r = analyze_one_N(data_dir, N, n_sample=500_000, rng=rng)
        results.append(r)
        gc.collect()

    # CSV output
    csv_path = out_dir / "K_eff_band_decomp.csv"
    keys = list(results[0].keys())
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for r in results:
            w.writerow([r[k] for k in keys])
    print(f"\n[csv] wrote {csv_path}")

    # Summary tables
    print(f"\n=== Cross-N: A/B start-range decomp ===")
    print(f"  {'N':>8}  {'K_full':>10}  {'K_int':>10}  {'K_above':>10}  {'K_below':>10}  "
          f"{'K_h-K_full':>10}")
    for r in results:
        print(f"  2^{r['log2N']:<6}  {r['K_eff_full']:>10.4f}  "
              f"{r['K_eff_interior']:>10.4f}  {r['K_eff_above']:>10.4f}  "
              f"{r['K_eff_below']:>10.4f}  {K_H - r['K_eff_full']:>+10.4f}")

    print(f"\n=== Cross-N: σ-band K_eff ===")
    print(f"  {'N':>8}  {'K_q1':>10}  {'K_q2':>10}  {'K_q3':>10}  {'K_q4':>10}  "
          f"({'q1 σ̄':>6}, {'q4 σ̄':>6})")
    for r in results:
        print(f"  2^{r['log2N']:<6}  {r['K_eff_q1_0_25']:>10.4f}  "
              f"{r['K_eff_q2_25_50']:>10.4f}  {r['K_eff_q3_50_75']:>10.4f}  "
              f"{r['K_eff_q4_75_100']:>10.4f}  ({r['sigma_mean_q1_0_25']:>6.1f}, "
              f"{r['sigma_mean_q4_75_100']:>6.1f})")

    # Spread analysis
    K_full = np.array([r["K_eff_full"] for r in results])
    K_int = np.array([r["K_eff_interior"] for r in results])
    K_above = np.array([r["K_eff_above"] for r in results])
    K_below = np.array([r["K_eff_below"] for r in results])
    K_q1 = np.array([r["K_eff_q1_0_25"] for r in results])
    K_q2 = np.array([r["K_eff_q2_25_50"] for r in results])
    K_q3 = np.array([r["K_eff_q3_50_75"] for r in results])
    K_q4 = np.array([r["K_eff_q4_75_100"] for r in results])

    print(f"\n=== Spread across N ∈ [2²⁰, 2²⁷] ===")
    print(f"  K_eff_full:     spread = {np.ptp(K_full):.4f}")
    print(f"  K_eff_interior: spread = {np.ptp(K_int):.4f}")
    print(f"  K_eff_above:    spread = {np.ptp(K_above):.4f}")
    print(f"  K_eff_below:    spread = {np.ptp(K_below):.4f}")
    print(f"  K_eff_q1:       spread = {np.ptp(K_q1):.4f}    (lowest σ band)")
    print(f"  K_eff_q2:       spread = {np.ptp(K_q2):.4f}")
    print(f"  K_eff_q3:       spread = {np.ptp(K_q3):.4f}")
    print(f"  K_eff_q4:       spread = {np.ptp(K_q4):.4f}    (highest σ band)")

    # Effect A verdict
    print(f"\n=== Effect A (sample-truncation bias) verdict ===")
    print(f"  K_full spread = {np.ptp(K_full):.4f}, K_int spread = {np.ptp(K_int):.4f}")
    if np.ptp(K_int) < 0.2 and np.ptp(K_full) > 0.4:
        a_v = "K_int N-stable while K_full oscillates → Effect A CONFIRMED"
    elif np.ptp(K_int) < np.ptp(K_full) * 0.5:
        a_v = "K_int much more stable than K_full → Effect A LIKELY contributes"
    elif np.ptp(K_int) > np.ptp(K_full) * 0.8:
        a_v = "K_int and K_full have similar spreads → Effect A NOT primary"
    else:
        a_v = "Mixed signal"
    print(f"  Verdict: {a_v}")

    # Effect B verdict
    print(f"\n=== Effect B (orbit-mass concentration) verdict ===")
    print(f"  K_above − K_below per N:")
    n_sig_diffs = 0
    for r in results:
        d = r["K_eff_above"] - r["K_eff_below"]
        below_w = r["K_eff_below_hi"] - r["K_eff_below_lo"]
        sig = abs(d) > below_w
        n_sig_diffs += int(sig)
        marker = " *" if sig else ""
        print(f"    2^{r['log2N']}: above={r['K_eff_above']:.3f}, below={r['K_eff_below']:.3f}, "
              f"diff={d:+.3f} (below CI width {below_w:.3f}){marker}")
    if n_sig_diffs >= 4:
        b_v = "K_above ≠ K_below at most N → Effect B CONFIRMED"
    elif n_sig_diffs >= 2:
        b_v = f"K_above ≠ K_below at {n_sig_diffs}/{len(results)} N → Effect B PARTIAL"
    else:
        b_v = "K_above ≈ K_below within CI → Effect B RULED OUT"
    print(f"  Verdict: {b_v}")

    # σ-band verdict
    print(f"\n=== σ-band K_eff stratification verdict ===")
    print(f"  Hypothesis: if K_eff per band is N-stable but K_eff_full oscillates,")
    print(f"              boundary effects shift band weights = mechanism.")
    print(f"  Per-band spreads: q1={np.ptp(K_q1):.3f}, q2={np.ptp(K_q2):.3f}, "
          f"q3={np.ptp(K_q3):.3f}, q4={np.ptp(K_q4):.3f}")
    avg_band_spread = np.mean([np.ptp(K_q1), np.ptp(K_q2), np.ptp(K_q3), np.ptp(K_q4)])
    print(f"  Avg per-band spread = {avg_band_spread:.4f}")
    print(f"  K_full spread       = {np.ptp(K_full):.4f}")
    if avg_band_spread < np.ptp(K_full) * 0.5:
        band_v = "Per-band K_eff is N-stable; K_full oscillation is from band-reweighting"
    elif avg_band_spread < np.ptp(K_full) * 0.8:
        band_v = "Per-band K_eff partially stable; band-reweighting contributes"
    else:
        band_v = "Per-band K_eff oscillates as much as K_full → reweighting NOT mechanism"
    print(f"  Verdict: {band_v}")

    # Asymptotic fit
    print(f"\n=== Asymptotic fit K_full(N) = K_h − C / (log N)^α ===")
    log_Ns = np.array([r["log_N"] for r in results])
    gap = K_H - K_full
    print(f"  log N values:    {log_Ns}")
    print(f"  K_h − K_full:    {gap}")
    if (gap > 0).all():
        x = np.log(log_Ns)
        y = np.log(gap)
        xm, ym = x.mean(), y.mean()
        slope_log = ((x - xm) * (y - ym)).sum() / ((x - xm) ** 2).sum()
        intercept_log = ym - slope_log * xm
        alpha_est = -slope_log
        C_est = np.exp(intercept_log)
        print(f"  α ≈ {alpha_est:.4f}    C ≈ {C_est:.4f}")
        print(f"  Tao bound exponent: 0.6")
    else:
        print(f"  NOT ALL GAPS POSITIVE — sign-changing. Cannot fit Tao-style form.")
        print(f"  → K_full oscillates around K_h, not monotone-converges.")

    # Plots
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        log2Ns = [r["log2N"] for r in results]

        # Plot 1: K_full vs K_interior
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.errorbar(log2Ns, K_full,
                    yerr=[K_full - [r["K_eff_full_lo"] for r in results],
                          [r["K_eff_full_hi"] for r in results] - K_full],
                    fmt="o-", label="K_eff_full")
        ax.errorbar(log2Ns, K_int,
                    yerr=[K_int - [r["K_eff_interior_lo"] for r in results],
                          [r["K_eff_interior_hi"] for r in results] - K_int],
                    fmt="s-", label="K_eff_interior (n<N/8)")
        ax.axhline(K_H, color="black", linestyle="--", alpha=0.5, label=f"K_h = {K_H:.3f}")
        ax.set_xlabel("log_2(N)")
        ax.set_ylabel("K_eff (first-passage R-slope)")
        ax.set_title("Plot 1: Effect A test (K_full vs K_interior)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / "plot1_full_vs_interior.png", dpi=120)
        plt.close()

        # Plot 2: K_above vs K_below
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.errorbar(log2Ns, K_above,
                    yerr=[K_above - [r["K_eff_above_lo"] for r in results],
                          [r["K_eff_above_hi"] for r in results] - K_above],
                    fmt="o-", label="K_eff_above (n > √N)")
        ax.errorbar(log2Ns, K_below,
                    yerr=[K_below - [r["K_eff_below_lo"] for r in results],
                          [r["K_eff_below_hi"] for r in results] - K_below],
                    fmt="s-", label="K_eff_below (n ≤ √N)")
        ax.axhline(K_H, color="black", linestyle="--", alpha=0.5, label=f"K_h = {K_H:.3f}")
        ax.set_xlabel("log_2(N)")
        ax.set_ylabel("K_eff (first-passage R-slope)")
        ax.set_title("Plot 2: Effect B test (K_above vs K_below)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / "plot2_above_vs_below.png", dpi=120)
        plt.close()

        # Plot 3: σ-band K_eff
        fig, ax = plt.subplots(figsize=(8, 5))
        for K_q, lab in [(K_q1, "q1 (0-25%)"), (K_q2, "q2 (25-50%)"),
                         (K_q3, "q3 (50-75%)"), (K_q4, "q4 (75-100%)")]:
            ax.plot(log2Ns, K_q, "o-", label=lab)
        ax.plot(log2Ns, K_full, "k--", linewidth=2, label="K_full (overall)")
        ax.axhline(K_H, color="black", linestyle=":", alpha=0.5, label=f"K_h = {K_H:.3f}")
        ax.set_xlabel("log_2(N)")
        ax.set_ylabel("K_eff per σ-quantile band")
        ax.set_title("Plot 3: σ-band K_eff stratification")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / "plot3_sigma_band.png", dpi=120)
        plt.close()
        print(f"\n[plots] saved to {out_dir}")
    except Exception as e:
        print(f"  [warn] plotting failed: {e}")

    print(f"\n[done]")


if __name__ == "__main__":
    main()
