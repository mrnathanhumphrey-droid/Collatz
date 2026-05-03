"""
Experiment 09 — Multi-statistic prefix decomposition

The prefix decomposition was developed for sigma (total stopping time). This
script tests whether the same a_final covariate parameterizes other Collatz
trajectory functionals — i.e., is a_final the right structural axis for
trajectory statistics in general, or specific to sigma?

For each outcome y in {sigma, syracuse, odd_steps, even_steps, log(max_excursion)}:
  1. Per-class OLS slope and intercept on log(n), at k=6 (32 odd residue classes mod 64)
  2. Moment-corrected per-class slope variance: tau_beta_sq = max(0, var(slope) - mean(SE^2))
  3. Predicted alpha_det from deterministic prefix algebra:
        alpha_det(r) = prefix_steps(r) + slope_heuristic * log(a_final(r) / 2^k)
     R^2 of alpha_actual vs alpha_det quantifies a_final's predictive power
  4. AD 2-sample test on per-class residual distributions, pairs grouped by
     same-a_final vs different-a_final. Ratio of median AD statistics quantifies
     a_final's distributional clustering power.

If a_final captures structure for all 5 outcomes:
  - tau_beta_sq -> 0 (slope universal across classes given prefix correction)
  - R^2 of alpha prediction high (~0.99+)
  - AD ratio (diff/same) >> 1 (same-a_final pairs come from same distribution)

Usage:
    python 09_multi_stat_decomposition.py --N 10000000 --k 6
"""
import argparse
import sys
from pathlib import Path
from itertools import combinations

import numpy as np
import polars as pl
from scipy.stats import anderson_ksamp

sys.stdout.reconfigure(encoding="utf-8")


def deterministic_prefix(r, a0, max_steps=400):
    """Return (prefix_steps, a_final, c_final) for residue r mod a0."""
    a, c = a0, r
    steps = 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3 * c + 1
        steps += 1
    return steps, a, c


def heuristic_slope(outcome_name):
    """Asymptotic slope of E[outcome | n] vs log(n) for odd n."""
    L4_3 = np.log(4.0) - np.log(3.0)
    if outcome_name == "sigma":
        return 3.0 / L4_3                 # ~10.4282
    if outcome_name == "syracuse":
        return 1.0 / L4_3                 # ~3.4761 (Syracuse compresses each odd step into one)
    if outcome_name == "odd_steps":
        return 1.0 / L4_3                 # same: each Syracuse step = one odd step
    if outcome_name == "even_steps":
        return 2.0 / L4_3                 # ~6.9521 (even/odd ratio is 2:1 in heuristic)
    if outcome_name == "log_max_excursion":
        # log(max_excursion) is bounded by sigma * log(3/2) above and by log(n) below.
        # Empirically scales with log(n); slope is unknown a priori, fit it.
        return None
    return None


def analyze_outcome(name, y, log_n, class_idx, K, a_final_per, c_final_per, prefix_per):
    """Return dict of analyses for one outcome variable."""
    # Per-class OLS
    alpha_actual = np.zeros(K)
    beta_actual = np.zeros(K)
    se_beta = np.zeros(K)
    eps_per = {}

    for k in range(K):
        m = class_idx == k
        if m.sum() < 50:
            continue
        x = log_n[m]
        yi = y[m]
        bk, ak = np.polyfit(x, yi, 1)
        beta_actual[k] = bk
        alpha_actual[k] = ak
        # SE of slope from residuals
        resid = yi - (ak + bk * x)
        eps_per[k] = resid
        s2 = np.sum(resid**2) / (len(yi) - 2)
        sxx = np.sum((x - x.mean())**2)
        se_beta[k] = np.sqrt(s2 / sxx) if sxx > 0 else np.nan

    # Slope universality: between-class variance vs within-class SE
    var_beta = np.var(beta_actual, ddof=1)
    mean_se_sq = np.mean(se_beta**2)
    tau_beta_sq_corrected = max(0.0, var_beta - mean_se_sq)
    mu_beta = float(np.mean(beta_actual))

    # Heuristic slope (where defined)
    slope_h = heuristic_slope(name)

    # Predicted alpha_det
    M = 2 ** int(np.log2(2 * (K + 1) - 2)) if K > 0 else 64  # fallback; we know K=32 → M=64
    M = 64  # fixed for k=6
    if slope_h is not None:
        # alpha_det(r) = prefix_steps(r) + slope_h * log(a_final(r) / M)
        alpha_det = np.zeros(K)
        for k in range(K):
            alpha_det[k] = prefix_per[k] + slope_h * np.log(a_final_per[k] / float(M))
    else:
        # log_max_excursion: fit slope_h from data via R^2 maximization?
        # Simpler: use mu_beta as the empirical slope.
        alpha_det = np.zeros(K)
        for k in range(K):
            alpha_det[k] = prefix_per[k] + mu_beta * np.log(a_final_per[k] / float(M))

    # R^2 of alpha_actual vs alpha_det
    if np.var(alpha_det) > 0:
        b_fit, a_fit = np.polyfit(alpha_det, alpha_actual, 1)
        fitted = a_fit + b_fit * alpha_det
        ss_total = np.sum((alpha_actual - alpha_actual.mean())**2)
        ss_resid = np.sum((alpha_actual - fitted)**2)
        r_sq = 1.0 - ss_resid / ss_total if ss_total > 0 else float("nan")
    else:
        r_sq = float("nan")

    # AD test: pair up classes, group by same-a_final vs different-a_final
    rng = np.random.default_rng(42)
    eps_sub = {}
    n_per_class = 20_000
    for k, eps in eps_per.items():
        if len(eps) > n_per_class:
            idx = rng.choice(len(eps), n_per_class, replace=False)
            eps_sub[k] = eps[idx]
        else:
            eps_sub[k] = eps

    same_a_stats, diff_a_stats = [], []
    for ki, kj in combinations(range(K), 2):
        if ki not in eps_sub or kj not in eps_sub:
            continue
        try:
            res = anderson_ksamp([eps_sub[ki], eps_sub[kj]])
            stat = res.statistic
        except Exception:
            stat = np.nan
        if a_final_per[ki] == a_final_per[kj]:
            same_a_stats.append(stat)
        else:
            diff_a_stats.append(stat)
    same_med = float(np.nanmedian(same_a_stats)) if same_a_stats else float("nan")
    diff_med = float(np.nanmedian(diff_a_stats)) if diff_a_stats else float("nan")
    ad_ratio = diff_med / same_med if abs(same_med) > 1e-6 else float("inf")

    return dict(
        outcome=name,
        K=K,
        mu_beta=mu_beta,
        slope_heuristic=slope_h,
        var_beta=var_beta,
        mean_se_sq=mean_se_sq,
        tau_beta_sq_corrected=tau_beta_sq_corrected,
        alpha_actual_sd=float(np.std(alpha_actual, ddof=1)),
        alpha_det_R2=r_sq,
        ad_same_med=same_med,
        ad_diff_med=diff_med,
        ad_ratio=ad_ratio,
        n_same_pairs=len(same_a_stats),
        n_diff_pairs=len(diff_a_stats),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10_000_000)
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--data", type=str, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data) if args.data else here.parent / "data"
    out_dir = Path(args.out) if args.out else here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] N={args.N:,}, k={args.k}", flush=True)
    df = pl.read_parquet(data_dir / f"main_N{args.N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n = df["n"].to_numpy().astype(np.int64)
    log_n = np.log(n.astype(np.float64))
    M = 1 << args.k
    res = (n % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)  # odd residues only, indexed 0..K-1
    K = M // 2  # number of odd classes
    print(f"        odd-n rows={len(n):,}, K={K} odd classes mod {M}")

    # Per-class deterministic prefix lookup
    a_final_per = np.zeros(K, dtype=np.int64)
    c_final_per = np.zeros(K, dtype=np.int64)
    prefix_per = np.zeros(K, dtype=np.int64)
    for k in range(K):
        r = 2 * k + 1
        s, a, c = deterministic_prefix(r, M)
        prefix_per[k] = s
        a_final_per[k] = a
        c_final_per[k] = c

    # Outcomes to analyze
    outcomes = {
        "sigma": df["sigma"].to_numpy().astype(np.float64),
        "syracuse": df["syracuse"].to_numpy().astype(np.float64),
        "odd_steps": df["odd_steps"].to_numpy().astype(np.float64),
        "even_steps": df["even_steps"].to_numpy().astype(np.float64),
        "log_max_excursion": np.log(df["max_excursion"].to_numpy().astype(np.float64)),
    }

    results = []
    for name, y in outcomes.items():
        print(f"\n=== {name} ===", flush=True)
        r = analyze_outcome(name, y, log_n, class_idx, K, a_final_per, c_final_per, prefix_per)
        results.append(r)
        print(f"  mu_beta={r['mu_beta']:.4f}  heuristic={r['slope_heuristic']}", flush=True)
        print(f"  var(beta)={r['var_beta']:.6f}  mean(SE^2)={r['mean_se_sq']:.6f}  "
              f"tau^2 corrected={r['tau_beta_sq_corrected']:.6f}", flush=True)
        print(f"  alpha_det R^2={r['alpha_det_R2']:.4f}  "
              f"SD(alpha_actual)={r['alpha_actual_sd']:.3f}", flush=True)
        print(f"  AD same-a_final median={r['ad_same_med']:.3f}  "
              f"AD diff-a_final median={r['ad_diff_med']:.3f}  "
              f"ratio={r['ad_ratio']:.1f}x", flush=True)
        print(f"  pairs: same={r['n_same_pairs']}, diff={r['n_diff_pairs']}", flush=True)

    # Summary table
    print("\n=== Summary (k=6, N={:,}) ===".format(args.N))
    cols = ["outcome", "mu_beta", "tau^2_corr", "alpha_R^2", "AD_same", "AD_diff", "AD_ratio"]
    print(f"{cols[0]:<20} {cols[1]:>10} {cols[2]:>12} {cols[3]:>10} {cols[4]:>10} {cols[5]:>10} {cols[6]:>10}")
    for r in results:
        print(f"{r['outcome']:<20} {r['mu_beta']:>10.4f} {r['tau_beta_sq_corrected']:>12.6f} "
              f"{r['alpha_det_R2']:>10.4f} {r['ad_same_med']:>10.3f} "
              f"{r['ad_diff_med']:>10.3f} {r['ad_ratio']:>10.1f}")

    # Save CSV
    summary = pl.DataFrame([
        {
            "outcome": r["outcome"],
            "mu_beta": r["mu_beta"],
            "slope_heuristic": r["slope_heuristic"],
            "var_beta": r["var_beta"],
            "mean_se_sq": r["mean_se_sq"],
            "tau_beta_sq_corrected": r["tau_beta_sq_corrected"],
            "alpha_actual_sd": r["alpha_actual_sd"],
            "alpha_det_R2": r["alpha_det_R2"],
            "ad_same_med": r["ad_same_med"],
            "ad_diff_med": r["ad_diff_med"],
            "ad_ratio": r["ad_ratio"],
            "n_same_pairs": r["n_same_pairs"],
            "n_diff_pairs": r["n_diff_pairs"],
        }
        for r in results
    ])
    out_csv = out_dir / f"09_multi_stat_decomposition_k{args.k}_N{args.N}.csv"
    summary.write_csv(out_csv)
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
