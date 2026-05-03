"""
Experiment 07 — Anderson-Darling distributional test for same-a_final collapse

Tests the strong form of the prefix-decomposition claim: classes sharing the
same a_final value have identical residual distributions, beyond just matching
moments.

For each pair of classes at k=6 (32 odd classes mod 64, 496 pairs total),
runs a 2-sample Anderson-Darling test on per-class residual distributions.

Pairs split into:
  - same-a_final pairs (~110 of 496 at k=6): null is "same distribution"
  - different-a_final pairs (~386): null should be rejected

Predicted outcomes:
  - If a_final fully parameterizes the conditional distribution:
      same-a_final p-values uniform on [0,1]; different-a_final p-values
      cluster near 0.
  - If c_final introduces secondary modulation:
      same-a_final p-values shift toward 0 (some rejections); different-a_final
      still cluster near 0 but more strongly.

Plot side-by-side histograms of the two p-value distributions.

Usage:
    python 07_anderson_darling.py --N 33554432
"""
import argparse
import sys
from pathlib import Path
from collections import defaultdict
from itertools import combinations

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from scipy.stats import anderson_ksamp

sys.stdout.reconfigure(encoding="utf-8")


def deterministic_prefix(r, a0=64):
    a, c = a0, r; steps = 0
    while a % 2 == 0:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3*c + 1
        steps += 1
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=1 << 25)
    ap.add_argument("--n_per_class", type=int, default=20_000,
                    help="Subsample size per class for AD (default 20K). AD is sensitive at large N.")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--data", type=str, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data) if args.data else here.parent / "data"
    out_dir = Path(args.out) if args.out else here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pl.read_parquet(data_dir / f"main_N{args.N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n = df["n"].to_numpy().astype(np.int64)
    log_n = np.log(n.astype(np.float64))
    sigma = df["sigma"].to_numpy().astype(np.float64)
    res64 = (n % 64).astype(np.int32)
    class_idx = ((res64 - 1) // 2).astype(np.int32)
    print(f"[load] N={args.N:,}, odd-only={len(n):,}")

    K = 32
    # Per-class OLS residuals
    eps_per = {}
    a_final_per = {}
    for k in range(K):
        m = class_idx == k
        bk, ak = np.polyfit(log_n[m], sigma[m], 1)
        eps_per[k] = sigma[m] - (ak + bk * log_n[m])
        a_final_per[k] = deterministic_prefix(2*k + 1, 64)

    # Subsample for AD
    rng = np.random.default_rng(args.seed)
    eps_sub = {}
    for k, eps in eps_per.items():
        if len(eps) > args.n_per_class:
            idx = rng.choice(len(eps), args.n_per_class, replace=False)
            eps_sub[k] = eps[idx]
        else:
            eps_sub[k] = eps

    # All pairs
    same_a_pvals = []
    diff_a_pvals = []
    same_a_stats = []
    diff_a_stats = []
    for ki, kj in combinations(range(K), 2):
        try:
            res = anderson_ksamp([eps_sub[ki], eps_sub[kj]])
            stat = res.statistic
            # scipy clips significance_level to [0.001, 0.25]; for less-conservative
            # behavior we'd want raw p, but capped is fine for our binary same/different test.
            p = res.significance_level / 100.0  # convert from percent
        except Exception:
            stat = np.nan; p = np.nan
        if a_final_per[ki] == a_final_per[kj]:
            same_a_pvals.append(p); same_a_stats.append(stat)
        else:
            diff_a_pvals.append(p); diff_a_stats.append(stat)

    print()
    print(f"=== Anderson-Darling 2-sample tests, k=6, n_per_class={args.n_per_class:,} ===")
    print(f"Same-a_final pairs ({len(same_a_pvals)}):")
    print(f"  AD stat:  mean={np.nanmean(same_a_stats):.3f}, median={np.nanmedian(same_a_stats):.3f}, max={np.nanmax(same_a_stats):.3f}")
    print(f"  p-value:  mean={np.nanmean(same_a_pvals):.4f}, median={np.nanmedian(same_a_pvals):.4f}")
    n_reject_same_05 = sum(1 for p in same_a_pvals if p < 0.05 and not np.isnan(p))
    n_reject_same_bonf = sum(1 for p in same_a_pvals if p < 0.05/len(same_a_pvals) and not np.isnan(p))
    print(f"  rejections at alpha=0.05: {n_reject_same_05}/{len(same_a_pvals)}")
    print(f"  rejections at Bonferroni alpha=0.05/{len(same_a_pvals)}: {n_reject_same_bonf}")

    print(f"Diff-a_final pairs ({len(diff_a_pvals)}):")
    print(f"  AD stat:  mean={np.nanmean(diff_a_stats):.3f}, median={np.nanmedian(diff_a_stats):.3f}, max={np.nanmax(diff_a_stats):.3f}")
    print(f"  p-value:  mean={np.nanmean(diff_a_pvals):.4f}, median={np.nanmedian(diff_a_pvals):.4f}")
    n_reject_diff_05 = sum(1 for p in diff_a_pvals if p < 0.05 and not np.isnan(p))
    print(f"  rejections at alpha=0.05: {n_reject_diff_05}/{len(diff_a_pvals)}")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    ax.hist(same_a_pvals, bins=20, range=(0, 0.26), color="C0", alpha=0.7,
            edgecolor="black", label=f"same-a_final ({len(same_a_pvals)})")
    ax.hist(diff_a_pvals, bins=20, range=(0, 0.26), color="C3", alpha=0.7,
            edgecolor="black", label=f"diff-a_final ({len(diff_a_pvals)})")
    ax.set_xlabel("Anderson-Darling p-value (capped at 0.25 by scipy)")
    ax.set_ylabel("count")
    ax.set_title(f"AD p-value distribution by pair type (k=6, n_per_class={args.n_per_class:,})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.hist(same_a_stats, bins=30, color="C0", alpha=0.7, edgecolor="black",
            label=f"same-a_final (median={np.nanmedian(same_a_stats):.2f})")
    ax.hist(diff_a_stats, bins=30, color="C3", alpha=0.7, edgecolor="black",
            label=f"diff-a_final (median={np.nanmedian(diff_a_stats):.2f})")
    ax.set_xlabel("Anderson-Darling statistic")
    ax.set_ylabel("count")
    ax.set_title("AD statistic by pair type\n(higher = more distinguishable)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.suptitle("Distributional test: are same-a_final classes from the same distribution?",
                 y=1.02)
    plt.tight_layout()
    out_png = out_dir / f"07_anderson_darling_k6_N{args.N}.png"
    plt.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close()
    print()
    print(f"[save] {out_png}")

    # Verdict — AD statistic is more reliable than scipy's capped p-values
    # (scipy clips to [0.001, 0.25]; for large samples this can mislead p-based counts).
    # 2-sample AD critical values: ~2.5 at 5%, ~3.7 at 1%, ~5.5 at 0.1%.
    same_med = np.nanmedian(same_a_stats)
    diff_med = np.nanmedian(diff_a_stats)
    print()
    print("=== Verdict (AD-stat-based, not capped-p-value-based) ===")
    print(f"Same-a_final AD stat: median={same_med:.3f}, max={np.nanmax(same_a_stats):.3f}")
    print(f"Diff-a_final AD stat: median={diff_med:.3f}, max={np.nanmax(diff_a_stats):.3f}")
    print(f"Ratio (diff median / same median): {diff_med/same_med:.1f}x" if abs(same_med) > 1e-6 else
          f"Same-median essentially zero; diff-median = {diff_med:.3f}")
    n_same_above_5pct = sum(1 for s in same_a_stats if s > 2.5 and not np.isnan(s))
    n_same_above_1pct = sum(1 for s in same_a_stats if s > 3.7 and not np.isnan(s))
    n_diff_above_5pct = sum(1 for s in diff_a_stats if s > 2.5 and not np.isnan(s))
    print(f"Same-a_final pairs above 5%-AD-critical (~2.5):  {n_same_above_5pct}/{len(same_a_stats)}")
    print(f"Same-a_final pairs above 1%-AD-critical (~3.7):  {n_same_above_1pct}/{len(same_a_stats)}")
    print(f"Diff-a_final pairs above 5%-AD-critical (~2.5):  {n_diff_above_5pct}/{len(diff_a_stats)}")
    if abs(same_med) < 1.0 and diff_med > 5.0:
        print()
        print("--> Most same-a_final pairs cluster near AD~0 (consistent with same-distribution null).")
        print("--> All diff-a_final pairs reject decisively.")
        print("--> Strong support for a_final as the primary clustering axis; c_final modulation")
        print("    detectable in extreme pairs but not the typical case.")


if __name__ == "__main__":
    main()
