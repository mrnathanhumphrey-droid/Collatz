"""
Experiment 44 — discriminate Effect A (sample-truncation bias) vs Effect B
(orbit-mass concentration) as the K_eff finite-N mechanism.

Per user prompt:
  Part 1 — sample-truncation bias.
    K_eff_full(N)     = slope σ vs log(n) on all odd n in [1, N].
    K_eff_interior(N) = same but restricted to n < N/8 (drops top 3 octaves).
    If interior is N-stable while full oscillates, Effect A confirmed.

  Part 2 — orbit-mass concentration.
    K_eff_above(N) = slope σ vs log(n) restricted to n > √N.
    K_eff_below(N) = slope σ vs log(n) restricted to n ≤ √N.
    If these differ substantially at any N, Effect B is real.

  Part 3 — asymptotic extrapolation:
    Fit K_eff(N) = K_h − C/(log N)^α if a mechanism is confirmed.
    α ≈ 0.6 would match Tao's O((log N)^0.6) bound.

Constraints (per prompt):
  - Use existing parquet only.
  - Cap ≤15 min wall.
  - Bootstrap CIs (1000 resamples).
  - Output CSV (N, four K_effs, CIs).

NOTE: This experiment's "K_eff" = slope of σ vs log(n) on a sample subset.
This is the *μ_β analog* (global stopping-time slope), NOT the same quantity
as Result 6's K_eff which is the slope of R(f) on log(f) from first-passage.
The two are related but not identical. The prompt frames both as testing the
same finite-N mechanism, so we report this μ_β-style slope variant.

Usage:
    python 44_finite_N_boundary_decomp.py
"""
import gc
import sys
import time
import csv
from pathlib import Path

import numpy as np
import polars as pl

sys.stdout.reconfigure(encoding="utf-8")

K_H = 3.0 / (np.log(4.0) - np.log(3.0))


def fit_slope_with_bootstrap(n_arr, sigma_arr, n_boot=1000, n_subsample=200_000, rng=None):
    """OLS slope of sigma on log(n), with bootstrap CI.

    For large arrays, subsample to n_subsample first (without replacement);
    bootstrap is done with replacement on the subsample (or full if smaller).

    Returns (slope, ci_lo, ci_hi, n_used).
    """
    rng = rng or np.random.default_rng(42)
    n = len(n_arr)
    if n < 5:
        return float("nan"), float("nan"), float("nan"), n
    if n > n_subsample:
        idx = rng.choice(n, n_subsample, replace=False)
        x = np.log(n_arr[idx].astype(np.float64))
        y = sigma_arr[idx].astype(np.float64)
    else:
        x = np.log(n_arr.astype(np.float64))
        y = sigma_arr.astype(np.float64)
    nu = len(x)
    # Direct OLS
    xm, ym = x.mean(), y.mean()
    den = ((x - xm) ** 2).sum()
    if den < 1e-12:
        return float("nan"), float("nan"), float("nan"), nu
    slope = ((x - xm) * (y - ym)).sum() / den
    # Bootstrap
    boot_slopes = np.zeros(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, nu, nu)
        xb = x[idx]
        yb = y[idx]
        xbm, ybm = xb.mean(), yb.mean()
        denb = ((xb - xbm) ** 2).sum()
        if denb < 1e-12:
            boot_slopes[b] = float("nan")
        else:
            boot_slopes[b] = ((xb - xbm) * (yb - ybm)).sum() / denb
    valid = ~np.isnan(boot_slopes)
    if valid.sum() > 10:
        ci_lo, ci_hi = np.percentile(boot_slopes[valid], [2.5, 97.5])
    else:
        ci_lo = ci_hi = float("nan")
    return float(slope), float(ci_lo), float(ci_hi), nu


def analyze_one_N(data_dir, N, rng=None):
    rng = rng or np.random.default_rng(42)
    df = pl.read_parquet(data_dir / f"main_N{N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    sigma_arr = df["sigma"].to_numpy().astype(np.int64)
    log2N = int(round(np.log2(N)))
    sqrt_N = int(np.sqrt(float(N)))
    cutoff_full_interior = N // 8  # drops top 3 octaves

    print(f"\n=== N = 2^{log2N} = {N:,}  ({len(n_arr):,} odd rows) ===", flush=True)
    print(f"  cutoff for interior (n < N/8 = 2^{log2N-3}): {cutoff_full_interior:,}")
    print(f"  cutoff for above/below (sqrt(N) ≈ 2^{log2N/2:.1f}): {sqrt_N:,}")

    out = {"N": N, "log2N": log2N, "log_N": np.log(float(N))}

    # K_eff_full
    t0 = time.perf_counter()
    s, lo, hi, nu = fit_slope_with_bootstrap(n_arr, sigma_arr, rng=rng)
    out["K_eff_full"] = s
    out["K_eff_full_lo"] = lo
    out["K_eff_full_hi"] = hi
    out["K_eff_full_n"] = nu
    print(f"  K_eff_full     = {s:.4f}  [95% CI {lo:.4f}, {hi:.4f}]  "
          f"(n={nu:,}, {time.perf_counter()-t0:.1f}s)")

    # K_eff_interior (n < N/8)
    mask = n_arr < cutoff_full_interior
    t0 = time.perf_counter()
    s, lo, hi, nu = fit_slope_with_bootstrap(n_arr[mask], sigma_arr[mask], rng=rng)
    out["K_eff_interior"] = s
    out["K_eff_interior_lo"] = lo
    out["K_eff_interior_hi"] = hi
    out["K_eff_interior_n"] = nu
    print(f"  K_eff_interior = {s:.4f}  [95% CI {lo:.4f}, {hi:.4f}]  "
          f"(n={nu:,}, {time.perf_counter()-t0:.1f}s)")

    # K_eff_above (n > √N)
    mask = n_arr > sqrt_N
    t0 = time.perf_counter()
    s, lo, hi, nu = fit_slope_with_bootstrap(n_arr[mask], sigma_arr[mask], rng=rng)
    out["K_eff_above"] = s
    out["K_eff_above_lo"] = lo
    out["K_eff_above_hi"] = hi
    out["K_eff_above_n"] = nu
    print(f"  K_eff_above    = {s:.4f}  [95% CI {lo:.4f}, {hi:.4f}]  "
          f"(n={nu:,}, {time.perf_counter()-t0:.1f}s)")

    # K_eff_below (n ≤ √N)
    mask = n_arr <= sqrt_N
    t0 = time.perf_counter()
    s, lo, hi, nu = fit_slope_with_bootstrap(n_arr[mask], sigma_arr[mask], rng=rng)
    out["K_eff_below"] = s
    out["K_eff_below_lo"] = lo
    out["K_eff_below_hi"] = hi
    out["K_eff_below_n"] = nu
    print(f"  K_eff_below    = {s:.4f}  [95% CI {lo:.4f}, {hi:.4f}]  "
          f"(n={nu:,}, {time.perf_counter()-t0:.1f}s)")

    return out


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "outputs" / "exp44"
    out_dir.mkdir(parents=True, exist_ok=True)

    N_values = [
        1 << 20, 1 << 22, 1 << 23, 1 << 24, 1 << 25, 1 << 27,
    ]

    print(f"K_eff finite-N boundary decomposition")
    print(f"K_h = {K_H:.6f}")
    print(f"K_eff defined here as slope of σ vs log(n) on the indicated subset.")
    print(f"Bootstrap: 1000 resamples, percentile CIs.")
    print(f"For full/interior: subsample to 200K when |subset| > 200K.")

    rng = np.random.default_rng(42)
    results = []
    for N in N_values:
        r = analyze_one_N(data_dir, N, rng=rng)
        results.append(r)
        gc.collect()

    # CSV output
    csv_path = out_dir / "K_eff_finite_N.csv"
    keys = [
        "N", "log2N", "log_N",
        "K_eff_full", "K_eff_full_lo", "K_eff_full_hi", "K_eff_full_n",
        "K_eff_interior", "K_eff_interior_lo", "K_eff_interior_hi", "K_eff_interior_n",
        "K_eff_above", "K_eff_above_lo", "K_eff_above_hi", "K_eff_above_n",
        "K_eff_below", "K_eff_below_lo", "K_eff_below_hi", "K_eff_below_n",
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for r in results:
            w.writerow([r[k] for k in keys])
    print(f"\n[csv] wrote {csv_path}")

    # Summary table
    print(f"\n=== Cross-N summary ===")
    print(f"  {'N':>8} {'K_full':>8} {'K_int':>8} {'K_above':>9} {'K_below':>9}  "
          f"{'full_CI':>16} {'int_CI':>16} {'above_CI':>16}")
    for r in results:
        print(f"  2^{r['log2N']:<6} "
              f"{r['K_eff_full']:>8.4f} {r['K_eff_interior']:>8.4f} "
              f"{r['K_eff_above']:>9.4f} {r['K_eff_below']:>9.4f}  "
              f"[{r['K_eff_full_lo']:.3f},{r['K_eff_full_hi']:.3f}]  "
              f"[{r['K_eff_interior_lo']:.3f},{r['K_eff_interior_hi']:.3f}]  "
              f"[{r['K_eff_above_lo']:.3f},{r['K_eff_above_hi']:.3f}]")

    # Spread analysis
    K_full = np.array([r["K_eff_full"] for r in results])
    K_int = np.array([r["K_eff_interior"] for r in results])
    K_above = np.array([r["K_eff_above"] for r in results])
    K_below = np.array([r["K_eff_below"] for r in results])

    print(f"\n=== Spread across N ∈ [2²⁰, 2²⁷] ===")
    print(f"  K_eff_full:     spread = {K_full.max() - K_full.min():.4f}")
    print(f"  K_eff_interior: spread = {K_int.max() - K_int.min():.4f}")
    print(f"  K_eff_above:    spread = {K_above.max() - K_above.min():.4f}")
    print(f"  K_eff_below:    spread = {K_below.max() - K_below.min():.4f}")

    # Effect A: K_eff_interior should be N-stable IF effect A is mechanism
    print(f"\n=== Effect A (sample-truncation bias) verdict ===")
    if K_int.max() - K_int.min() < K_full.max() - K_full.min() / 2:
        a_verdict = "K_int much more stable than K_full → Effect A LIKELY contributing"
    elif K_int.max() - K_int.min() < 0.05:
        a_verdict = "K_int N-stable (<0.05 spread) → Effect A CONFIRMED"
    else:
        a_verdict = "K_int still oscillates with similar amplitude → Effect A NOT primary"
    print(f"  K_full spread = {K_full.max() - K_full.min():.4f}")
    print(f"  K_int  spread = {K_int.max() - K_int.min():.4f}")
    print(f"  Verdict: {a_verdict}")

    # Effect B: K_eff_above vs K_eff_below should differ at any single N
    print(f"\n=== Effect B (orbit-mass concentration) verdict ===")
    diffs = K_above - K_below
    print(f"  K_above - K_below per N:")
    for r, d in zip(results, diffs):
        print(f"    2^{r['log2N']}: above={r['K_eff_above']:.3f}, below={r['K_eff_below']:.3f}, diff={d:+.3f}  "
              f"(below CI width = {r['K_eff_below_hi']-r['K_eff_below_lo']:.3f})")
    sig_diffs = []
    for r in results:
        below_width = r["K_eff_below_hi"] - r["K_eff_below_lo"]
        diff = r["K_eff_above"] - r["K_eff_below"]
        if abs(diff) > below_width:
            sig_diffs.append(True)
        else:
            sig_diffs.append(False)
    if all(sig_diffs):
        b_verdict = "K_above ≠ K_below at all N (above below CI width) → Effect B LIKELY"
    elif any(sig_diffs):
        b_verdict = f"K_above ≠ K_below at {sum(sig_diffs)}/{len(sig_diffs)} N → Effect B PARTIAL"
    else:
        b_verdict = "K_above ≈ K_below within CI → Effect B RULED OUT"
    print(f"  Verdict: {b_verdict}")

    # Part 3: asymptotic fit
    print(f"\n=== Asymptotic fit K_full(N) = K_h − C / (log N)^α ===")
    log_Ns = np.array([r["log_N"] for r in results])
    gap = K_H - K_full
    print(f"  log N values:    {log_Ns}")
    print(f"  K_h − K_full:    {gap}")
    if (gap > 0).all():
        # Fit log(gap) = log(C) − α · log(log N)
        x = np.log(log_Ns)
        y = np.log(gap)
        xm, ym = x.mean(), y.mean()
        slope_log = ((x - xm) * (y - ym)).sum() / ((x - xm) ** 2).sum()
        intercept_log = ym - slope_log * xm
        alpha_est = -slope_log
        C_est = np.exp(intercept_log)
        print(f"  All gaps positive → log-log fit feasible.")
        print(f"  α ≈ {alpha_est:.4f}    C ≈ {C_est:.4f}")
        print(f"  Tao bound exponent: 0.6   (match if α ≈ 0.6)")
    else:
        print(f"  NOT ALL GAPS POSITIVE — sign-changing K_h − K_full.")
        print(f"  Cannot fit K_h − C/(log N)^α (model can't accommodate sign changes).")
        print(f"  → Suggests K_full oscillates around K_h, not monotone-converges.")

    # Plots
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 5))
        log2Ns = [r["log2N"] for r in results]
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
        ax.set_ylabel("K_eff = slope of σ vs log(n)")
        ax.set_title("Plot 1: K_eff_full vs K_eff_interior (Effect A test)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / "plot1_full_vs_interior.png", dpi=120)
        plt.close()

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
        ax.set_ylabel("K_eff = slope of σ vs log(n)")
        ax.set_title("Plot 2: K_eff_above vs K_eff_below (Effect B test)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / "plot2_above_vs_below.png", dpi=120)
        plt.close()
        print(f"\n[plots] saved to {out_dir}")
    except Exception as e:
        print(f"  [warn] plotting failed: {e}")

    print(f"\n[done]")


if __name__ == "__main__":
    main()
