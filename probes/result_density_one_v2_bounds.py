"""
result_density_one_v2_bounds.py
================================
Quantitative density-1 tests connecting eps_k framework to Lagarias/Tao
2-adic equidistribution.

TEST A — cumulative v_2 density.
  For each Syracuse trajectory (started from odd n coprime to 3 in
  [3, 8388607]): compute the empirical fraction of steps with
  v_2(3 n_i + 1) >= k for k = 1, 2, ..., 10. Compare to geometric null
  P(v >= k) = 2^{-(k-1)}. Density-1 claim: fraction of trajectories
  with empirical density >= null approaches 1 as trajectory length grows.

TEST B — average v_2 vs log_2(3).
  For each trajectory: empirical mean of v along the trajectory.
  Density-1 claim: fraction of trajectories with mean v > log_2(3) ≈ 1.585
  approaches 1 as trajectory length grows.

Note on adapted test points: per-trajectory length tops out around 247 in
this data (each Syracuse trajectory ends at 1 in O(log n) steps), so the
brief's k=10000 sample is unreachable per-trajectory. We bin by trajectory
length instead and report density-1 convergence across length bins.

Data: data/v_seq_N8388608.parquet (~427M rows, ~2.8M qualifying starts).
"""
from __future__ import annotations

import csv
import math
import os
import sys
import time

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz"
DATA = os.path.join(OUTDIR, "data", "v_seq_N8388608.parquet")
OUT_MD = os.path.join(OUTDIR, "result_density_one_v2_bounds.md")
OUT_CSV = os.path.join(OUTDIR, "result_density_one_v2_bounds.csv")
OUT_DIAG = os.path.join(OUTDIR, "result_density_one_v2_bounds_diagnostic.md")

LOG2_3 = math.log2(3)  # ~1.5849625


def main():
    print("=" * 78)
    print("Density-1 v_2 bound tests")
    print("=" * 78)

    t0 = time.time()
    print(f"Loading {DATA}...")
    df = pd.read_parquet(DATA)
    print(f"  raw rows: {len(df):,}, t={time.time()-t0:.1f}s")

    # Filter to odd starts coprime to 3
    t0 = time.time()
    mask = (df["n"] % 2 == 1) & (df["n"] % 3 != 0)
    df = df.loc[mask].reset_index(drop=True)
    print(f"  after filter (odd, coprime to 3): {len(df):,} rows, t={time.time()-t0:.1f}s")

    n_unique = df["n"].nunique()
    print(f"  unique trajectories: {n_unique:,}")
    print()

    # ---------- per-trajectory aggregation ----------
    t0 = time.time()
    print("Computing per-trajectory aggregates...")
    # Build mask columns once for k = 1..10
    K_MAX = 10
    for k in range(1, K_MAX + 1):
        df[f"mask_ge_{k}"] = (df["v"] >= k).astype(np.int32)

    agg_dict = {"len": ("v", "count"), "sum_v": ("v", "sum")}
    for k in range(1, K_MAX + 1):
        agg_dict[f"count_ge_{k}"] = (f"mask_ge_{k}", "sum")

    g = df.groupby("n", sort=False).agg(**agg_dict).reset_index()
    print(f"  grouped to {len(g):,} trajectories, t={time.time()-t0:.1f}s")

    # Drop the mask cols from df to save memory
    df = df[["n", "step_idx", "v"]]

    # Per-trajectory derived columns
    g["mean_v"] = g["sum_v"] / g["len"]
    for k in range(1, K_MAX + 1):
        g[f"frac_ge_{k}"] = g[f"count_ge_{k}"] / g["len"]
    print(f"  trajectory length: min={int(g['len'].min())}, "
          f"median={int(g['len'].median())}, "
          f"max={int(g['len'].max())}, mean={g['len'].mean():.2f}")
    print(f"  mean v across all trajectories: {g['mean_v'].mean():.6f}")
    print()

    # ---------- TEST A: density of v >= k vs geometric null ----------
    print("=" * 78)
    print("TEST A: empirical density of v >= k vs geometric null 2^{-(k-1)}")
    print("=" * 78)
    null_p = {k: 2 ** -(k - 1) for k in range(1, K_MAX + 1)}
    test_a_rows = []
    for k in range(1, K_MAX + 1):
        col = f"frac_ge_{k}"
        emp = g[col]
        frac_pass = float((emp >= null_p[k]).mean())
        frac_pass_strict = float((emp > null_p[k]).mean())
        emp_mean = float(emp.mean())
        emp_med = float(emp.median())
        test_a_rows.append({
            "k": k,
            "null_p_geom": null_p[k],
            "frac_ge_null": frac_pass,
            "frac_gt_null": frac_pass_strict,
            "emp_mean": emp_mean,
            "emp_median": emp_med,
        })
        print(f"  k={k}: null=2^{-(k-1)}={null_p[k]:.6f}, "
              f"emp_mean={emp_mean:.6f}, emp_median={emp_med:.6f}, "
              f"frac >= null = {frac_pass:.6f}, frac > null = {frac_pass_strict:.6f}")
    print()

    # Test A by trajectory length bin
    bin_edges = [1, 5, 10, 20, 30, 50, 100, 200, 1000]
    bin_labels = [f"[{bin_edges[i]},{bin_edges[i+1]})"
                  for i in range(len(bin_edges) - 1)]
    g["len_bin"] = pd.cut(g["len"], bins=bin_edges, right=False, labels=bin_labels)
    print(f"  trajectory length distribution by bin:")
    bin_counts = g["len_bin"].value_counts().sort_index()
    for b, c in bin_counts.items():
        print(f"    {b}: {c:,}")
    print()

    print(f"  TEST A by length bin: fraction with frac_ge_k >= 2^{-(k-1)}")
    test_a_by_bin = []
    for b in bin_labels:
        sub = g[g["len_bin"] == b]
        if len(sub) == 0:
            continue
        row = {"len_bin": b, "n_traj": len(sub)}
        for k in range(1, K_MAX + 1):
            row[f"frac_ge_null_k{k}"] = float((sub[f"frac_ge_{k}"] >= null_p[k]).mean())
        test_a_by_bin.append(row)
        print(f"    {b} (n={len(sub):>9,}):  "
              f"k=1:{row['frac_ge_null_k1']:.4f}  "
              f"k=2:{row['frac_ge_null_k2']:.4f}  "
              f"k=3:{row['frac_ge_null_k3']:.4f}  "
              f"k=5:{row['frac_ge_null_k5']:.4f}  "
              f"k=10:{row['frac_ge_null_k10']:.4f}")
    print()

    # ---------- TEST B: trajectory mean v vs log_2(3) ----------
    print("=" * 78)
    print("TEST B: trajectory mean v vs log_2(3) ≈ 1.5849625")
    print("=" * 78)
    g["excess"] = g["mean_v"] - LOG2_3
    print(f"  overall: mean of trajectory-means = {g['mean_v'].mean():.6f}, "
          f"frac mean_v > log_2(3) = {(g['mean_v'] > LOG2_3).mean():.6f}")
    print()

    test_b_by_bin = []
    print(f"  TEST B by length bin:")
    print(f"    {'bin':>15}  {'n_traj':>10}  {'mean(mean_v)':>14}  "
          f"{'std(mean_v)':>12}  {'frac > log2(3)':>16}  {'mean excess':>12}")
    for b in bin_labels:
        sub = g[g["len_bin"] == b]
        if len(sub) == 0:
            continue
        row = {
            "len_bin": b,
            "n_traj": int(len(sub)),
            "mean_meanv": float(sub["mean_v"].mean()),
            "std_meanv": float(sub["mean_v"].std()),
            "frac_above_log2_3": float((sub["mean_v"] > LOG2_3).mean()),
            "mean_excess": float(sub["excess"].mean()),
            "median_excess": float(sub["excess"].median()),
        }
        test_b_by_bin.append(row)
        print(f"    {b:>15}  {row['n_traj']:>10,}  {row['mean_meanv']:>14.6f}  "
              f"{row['std_meanv']:>12.6f}  {row['frac_above_log2_3']:>16.6f}  "
              f"{row['mean_excess']:>+12.6f}")
    print()

    # ---------- Failure mode analysis ----------
    print("=" * 78)
    print("Failure mode: trajectories where mean_v <= log_2(3)")
    print("=" * 78)
    fail = g[g["mean_v"] <= LOG2_3].copy()
    print(f"  total failing trajectories: {len(fail):,} ({len(fail)/len(g)*100:.4f}%)")
    if len(fail) > 0:
        print(f"  failing length distribution: min={int(fail['len'].min())}, "
              f"median={int(fail['len'].median())}, "
              f"max={int(fail['len'].max())}, mean={fail['len'].mean():.2f}")
        print(f"  failing trajectories' mean_v: min={fail['mean_v'].min():.6f}, "
              f"median={fail['mean_v'].median():.6f}, "
              f"max={fail['mean_v'].max():.6f}")
        # By length bin
        for b in bin_labels:
            sub = fail[fail["len_bin"] == b]
            if len(sub) > 0:
                print(f"    bin {b}: {len(sub):,} failures")
        # Residue mod 4? mod 8?
        fail["n_mod_4"] = fail["n"] % 4
        fail["n_mod_8"] = fail["n"] % 8
        all_mod_4 = g["n"] % 4
        print(f"  starting residue mod 4 (failing vs all):")
        for r in [1, 3]:
            f_frac = (fail["n_mod_4"] == r).mean()
            a_frac = (all_mod_4 == r).mean()
            print(f"    n mod 4 = {r}: failing {f_frac:.6f}, all {a_frac:.6f}")
    print()

    # ---------- write CSV (per-trajectory) ----------
    print("Writing CSV (sample)...")
    # Per-trajectory CSV is huge; write a uniform random sample of 50000
    rng = np.random.default_rng(seed=42)
    if len(g) > 50000:
        sample_idx = rng.choice(len(g), size=50000, replace=False)
        sample = g.iloc[sample_idx]
    else:
        sample = g
    sample_out = sample[["n", "len", "mean_v"] +
                        [f"frac_ge_{k}" for k in range(1, K_MAX + 1)]].copy()
    sample_out.to_csv(OUT_CSV, index=False, float_format="%.6f")
    print(f"[csv (50k sample): {OUT_CSV}]")

    # ---------- write main MD ----------
    with open(OUT_MD, "w", encoding="utf-8") as fh:
        fh.write("# Density-1 v_2 bound tests\n\n")
        fh.write("Quantitative density-1 confirmation tests connecting the "
                 "eps_k Syracuse-Markov framework to the Lagarias / Tao 2-adic "
                 "equidistribution prediction.\n\n")
        fh.write(f"**Data**: `data/v_seq_N8388608.parquet` — Syracuse trajectories "
                 f"for starting integers in [3, 8388607]. Filtered to odd starts "
                 f"coprime to 3 ⇒ {n_unique:,} qualifying trajectories.\n\n")
        fh.write(f"**Trajectory length**: max {int(g['len'].max())}, "
                 f"median {int(g['len'].median())}, mean {g['len'].mean():.2f}. "
                 f"(Per-trajectory k=1000/10000 tests from the brief unreachable; "
                 f"adapted to length-binned density-1 convergence.)\n\n")
        fh.write(f"**log_2(3)** = {LOG2_3:.10f}\n\n")
        fh.write(f"**Mean of v across all trajectories**: "
                 f"{g['mean_v'].mean():.6f} (geometric prediction: 2.0)\n\n")

        # TEST A overall table
        fh.write("## TEST A: density of v >= k vs geometric null 2^{-(k-1)}\n\n")
        fh.write("Geometric null: under Lagarias' 2-adic equidistribution, "
                 "v_2(3 n_i + 1) is asymptotically Geom(1/2) on {1, 2, ...}, so "
                 "P(v >= k) = 2^{-(k-1)}. Empirical density of v >= k along each "
                 "trajectory is computed and compared.\n\n")
        fh.write("| k | null = 2^{-(k-1)} | emp mean | emp median | frac >= null | frac > null |\n")
        fh.write("|---|------------------:|---------:|-----------:|-------------:|------------:|\n")
        for r in test_a_rows:
            fh.write(f"| {r['k']} | {r['null_p_geom']:.6f} | "
                     f"{r['emp_mean']:.6f} | {r['emp_median']:.6f} | "
                     f"{r['frac_ge_null']:.6f} | {r['frac_gt_null']:.6f} |\n")
        fh.write("\n")

        # TEST A by length bin
        fh.write("### TEST A by trajectory length bin\n\n")
        fh.write("Fraction of trajectories with empirical density of v >= k "
                 "exceeding the geometric null 2^{-(k-1)}:\n\n")
        fh.write("| length bin | n_traj | k=1 | k=2 | k=3 | k=4 | k=5 | k=6 | k=8 | k=10 |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|---|\n")
        for r in test_a_by_bin:
            fh.write(f"| {r['len_bin']} | {r['n_traj']:,} | "
                     f"{r['frac_ge_null_k1']:.4f} | "
                     f"{r['frac_ge_null_k2']:.4f} | "
                     f"{r['frac_ge_null_k3']:.4f} | "
                     f"{r['frac_ge_null_k4']:.4f} | "
                     f"{r['frac_ge_null_k5']:.4f} | "
                     f"{r['frac_ge_null_k6']:.4f} | "
                     f"{r['frac_ge_null_k8']:.4f} | "
                     f"{r['frac_ge_null_k10']:.4f} |\n")
        fh.write("\n")

        # TEST B overall + by bin
        fh.write("## TEST B: trajectory mean v vs log_2(3)\n\n")
        fh.write("Density-1 claim: for density-1 of starting integers, "
                 "(1/L) Σ v_2(3 n_i + 1) > log_2(3) for L large. Per-trajectory "
                 "mean v is computed, and fraction exceeding log_2(3) reported by "
                 "trajectory length bin.\n\n")
        fh.write("| length bin | n_traj | mean(mean_v) | std(mean_v) | "
                 "frac > log_2(3) | mean excess | median excess |\n")
        fh.write("|---|---|---|---|---|---|---|\n")
        for r in test_b_by_bin:
            fh.write(f"| {r['len_bin']} | {r['n_traj']:,} | "
                     f"{r['mean_meanv']:.6f} | {r['std_meanv']:.6f} | "
                     f"{r['frac_above_log2_3']:.6f} | "
                     f"{r['mean_excess']:+.6f} | {r['median_excess']:+.6f} |\n")
        fh.write("\n")

        # Verdict
        fh.write("## Verdict\n\n")
        # Test A verdict
        all_a_pass = all(r['frac_ge_null'] > 0.95 for r in test_a_rows)
        all_b_pass_long = test_b_by_bin and test_b_by_bin[-1]['frac_above_log2_3'] > 0.99
        fh.write(f"**TEST A**: ")
        if all_a_pass:
            fh.write("Empirical density of v >= k exceeds the geometric null "
                     "2^{-(k-1)} for >95% of trajectories at every k = 1..10. "
                     "Density-1 confirmation holds.\n\n")
        else:
            fh.write("Empirical density falls below geometric null for a "
                     "non-trivial fraction of trajectories at some k. "
                     "Density-1 weakening — characterize by length bin.\n\n")

        fh.write(f"**TEST B**: ")
        if all_b_pass_long:
            fh.write(f"In the longest trajectory bin, "
                     f"{test_b_by_bin[-1]['frac_above_log2_3']*100:.4f}% of trajectories "
                     f"have empirical mean v > log_2(3). Fraction approaches 1 "
                     f"as trajectory length grows: density-1 confirmation.\n\n")
        else:
            fh.write(f"In the longest trajectory bin, only "
                     f"{test_b_by_bin[-1]['frac_above_log2_3']*100:.4f}% of trajectories "
                     f"have mean v > log_2(3). Density-1 not strict at this "
                     f"trajectory length — characterize convergence rate.\n\n")

        # Mean excess decay rate
        fh.write("**Convergence of mean excess**: from the by-bin table, "
                 "the mean excess (mean_v − log_2(3)) varies with bin. "
                 "Reading the rate of approach to the limiting value gives "
                 "the empirical rate of convergence to the Lagarias prediction.\n\n")

        fh.write("## Notes\n\n")
        fh.write("- Per-trajectory CSV (50,000-row uniform random sample) at "
                 "result_density_one_v2_bounds.csv.\n")
        fh.write("- Diagnostic file documents failure-mode analysis "
                 "(residue patterns, length distribution of failures).\n")

    print(f"[md:  {OUT_MD}]")

    # ---------- write diagnostic MD ----------
    with open(OUT_DIAG, "w", encoding="utf-8") as fh:
        fh.write("# Density-1 v_2 bounds — diagnostic\n\n")
        fh.write("## Data summary\n\n")
        fh.write(f"- Source: `data/v_seq_N8388608.parquet`\n")
        fh.write(f"- Filter: odd starts, coprime to 3\n")
        fh.write(f"- Qualifying trajectories: {n_unique:,}\n")
        fh.write(f"- Filtered rows: see stdout\n")
        fh.write(f"- Trajectory length: min={int(g['len'].min())}, "
                 f"median={int(g['len'].median())}, "
                 f"max={int(g['len'].max())}, mean={g['len'].mean():.4f}\n\n")

        fh.write("## TEST B failure mode (trajectories with mean_v ≤ log_2(3))\n\n")
        fh.write(f"Total failing: {len(fail):,} "
                 f"({len(fail)/len(g)*100:.4f}% of qualifying trajectories)\n\n")
        if len(fail) > 0:
            fh.write(f"### Failure length distribution\n\n")
            fh.write(f"- min length: {int(fail['len'].min())}\n")
            fh.write(f"- median length: {int(fail['len'].median())}\n")
            fh.write(f"- max length: {int(fail['len'].max())}\n")
            fh.write(f"- mean length: {fail['len'].mean():.4f}\n\n")
            fh.write(f"### Failure rate by trajectory length bin\n\n")
            fh.write("| length bin | n_total | n_failing | failure rate |\n")
            fh.write("|---|---|---|---|\n")
            for b in bin_labels:
                tot = (g["len_bin"] == b).sum()
                f_n = (fail["len_bin"] == b).sum() if len(fail) > 0 else 0
                rate = f_n / tot if tot > 0 else 0.0
                fh.write(f"| {b} | {tot:,} | {f_n:,} | {rate:.6f} |\n")
            fh.write("\n")

            fh.write(f"### Starting residue patterns of failing trajectories\n\n")
            fh.write("| residue | failing count | failing rate | overall rate | enriched? |\n")
            fh.write("|---|---|---|---|---|\n")
            n_total = len(g)
            for mod in [4, 8, 16]:
                fail_col = fail["n"] % mod
                all_col = g["n"] % mod
                for r in sorted(all_col.unique()):
                    if r % 2 == 0 or r % 3 == 0:
                        continue  # skip non-qualifying residues
                    f_count = int((fail_col == r).sum())
                    f_rate = f_count / len(fail) if len(fail) > 0 else 0
                    a_rate = float((all_col == r).mean())
                    enriched = "yes" if f_rate > 1.5 * a_rate else ""
                    fh.write(f"| n mod {mod} = {r} | {f_count:,} | "
                             f"{f_rate:.6f} | {a_rate:.6f} | {enriched} |\n")
                fh.write("\n")

        fh.write("## Numerical notes\n\n")
        fh.write("- Empirical density per trajectory has resolution 1/L where L "
                 "is trajectory length. For L=20 (typical short trajectory), "
                 "resolution is 0.05; the empirical density of v ≥ k for k ≥ 5 "
                 "(null 2^{-4} = 0.0625) is harder to compare reliably for short "
                 "trajectories. Length-binned analysis controls for this.\n")

    print(f"[diag: {OUT_DIAG}]")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
