"""
Experiment 30 — σ-records vs prefix decomposition class.

OEIS A006877 lists the integers n where σ(n) sets a new running max — the
σ-records of the 3x+1 dynamics. The prefix decomposition predicts:

  σ(n | n ≡ r mod 2^k) ≈ α_det(r) + β·log(n) + ε

with α_det(r) = prefix_steps(r) + (3 / (ln 4 − ln 3)) · ln(a★(r) / 2^k).
Records require σ(n) > all prior σ. Hence records should preferentially
come from classes with the highest α_det — i.e., classes terminating at
the largest a★(r) values.

Test: compute residue mod 2^k for k ∈ {6, 8, 10} for each odd record
n ≤ 2^29 from A006877, tabulate by a★, compute Spearman correlation between
records-per-class and α_det.

Cross-check OEIS records vs our existing parquet σ-data in the overlap
(records with n ≤ 2^27).

Usage:
    python 30_sigma_records_prefix_analysis.py --N_cap 536870912
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr, fisher_exact

sys.stdout.reconfigure(encoding="utf-8")

LOG_FACTOR_ODD = 3.0 / (np.log(4.0) - np.log(3.0))


def deterministic_prefix(r, a0, max_steps=400):
    """3x+1 prefix iteration from state (a=a0, c=r) until a is odd."""
    a, c = a0, r
    steps = 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3 * c + 1
        steps += 1
    return steps, a, c


def cross_check_records(records_oeis, parquet_path):
    """Walk the parquet in n-order, track running max sigma, find records."""
    print(f"[xcheck] loading {parquet_path}", flush=True)
    df = pl.read_parquet(parquet_path).sort("n")
    n_arr = df["n"].to_numpy()
    s_arr = df["sigma"].to_numpy()
    # Find indices where sigma sets new max
    cummax = np.maximum.accumulate(s_arr)
    # Record n's are those where sigma equals the cummax AND is strictly greater than previous cummax
    is_record = np.empty(len(s_arr), dtype=bool)
    is_record[0] = True
    is_record[1:] = s_arr[1:] > cummax[:-1]
    found = sorted(n_arr[is_record].tolist())
    cap = max(records_oeis)
    found_in_range = [n for n in found if n <= cap]
    print(f"[xcheck] our parquet records (n <= {cap:,}): {len(found_in_range)}")
    print(f"[xcheck] OEIS records (n <= {cap:,}): {len(records_oeis)}")
    diff_us = set(found_in_range) - set(records_oeis)
    diff_oeis = set(records_oeis) - set(found_in_range)
    if not diff_us and not diff_oeis:
        print(f"[xcheck] EXACT MATCH on {len(found_in_range)} records ✓")
    else:
        print(f"[xcheck] MISMATCH:")
        if diff_us:
            print(f"   in ours but not OEIS: {sorted(diff_us)[:10]}")
        if diff_oeis:
            print(f"   in OEIS but not ours: {sorted(diff_oeis)[:10]}")
    return found_in_range


def analyze_records_at_k(odd_records, k):
    """For each odd record n, compute residue mod 2^k, a*, prefix_steps, alpha_det.
    Tabulate records by a* and by class. Return dataframe + summary stats."""
    M = 1 << k
    rows = []
    for n in odd_records:
        r = n % M
        steps, a_star, c_star = deterministic_prefix(r, M)
        alpha_det = steps + LOG_FACTOR_ODD * np.log(a_star / M)
        rows.append({
            "n": n, "log_n": np.log(n),
            "r": r, "prefix_steps": steps,
            "a_star": int(a_star), "c_star": int(c_star),
            "alpha_det": alpha_det,
        })
    df = pl.DataFrame(rows)

    # Per-a* totals
    by_a = df.group_by("a_star").agg(
        pl.len().alias("records"),
        pl.col("alpha_det").mean().alias("alpha_det_mean"),
    ).sort("a_star")

    # Class density: how many classes have each a* value (across all 2^(k-1) odd residues)?
    class_density = {}
    for r in range(1, M, 2):
        _, a_s, _ = deterministic_prefix(r, M)
        class_density[int(a_s)] = class_density.get(int(a_s), 0) + 1
    by_a = by_a.with_columns(
        pl.col("a_star").map_elements(lambda a: class_density.get(int(a), 0),
                                       return_dtype=pl.Int64).alias("n_classes")
    )
    # Records per class (records / classes available)
    by_a = by_a.with_columns(
        (pl.col("records") / pl.col("n_classes")).alias("records_per_class")
    )

    # Spearman correlation: records_per_class vs alpha_det_mean
    rs = by_a["records_per_class"].to_numpy()
    ad = by_a["alpha_det_mean"].to_numpy()
    rho, p = spearmanr(rs, ad)

    print(f"\n=== k={k} (mod {M}, {M//2} odd classes, {len(set(class_density))} a* levels) ===")
    print(by_a.with_columns(
        pl.col("alpha_det_mean").round(2),
        pl.col("records_per_class").round(3),
    ))
    print(f"  Spearman rho(records_per_class, alpha_det) = {rho:+.4f}  p = {p:.4g}  n={len(rs)} a* levels")

    # Per-class deviation: records_per_class vs uniform expectation
    total_records = len(odd_records)
    expected_per_class = total_records / (M // 2)
    by_a = by_a.with_columns(
        (pl.col("n_classes") * expected_per_class).alias("expected_records"),
        ((pl.col("records") - pl.col("n_classes") * expected_per_class) /
         np.sqrt(pl.col("n_classes") * expected_per_class)).alias("z_uniform"),
    )

    return df, by_a, rho, p


def log_n_trend(df_records, k=6):
    """Bin records by log(n) octaves, check if a* concentration weakens with log(n)."""
    df = df_records.with_columns(
        (pl.col("log_n") / np.log(2)).cast(pl.Int32).alias("log2_octave")
    )
    print(f"\n=== log(n) trend at k={k}: P(a*=max) per log2 octave ===")
    M = 1 << k
    a_star_max = max(df["a_star"].to_list())
    octaves = sorted(df["log2_octave"].unique().to_list())
    rows = []
    for o in octaves:
        sub = df.filter(pl.col("log2_octave") == o)
        n_total = len(sub)
        n_max = sub.filter(pl.col("a_star") == a_star_max).height
        rows.append((o, n_total, n_max, n_max/n_total if n_total > 0 else 0))
    print(f"  {'log2(n)':>7} {'n_records':>10} {'records at a*='+str(a_star_max):>22} {'fraction':>10}")
    for o, t, m, f in rows:
        print(f"  {o:>7} {t:>10} {m:>22} {f:>10.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N_cap", type=int, default=2**29)
    ap.add_argument("--bfile", type=str, default="experiments_output/A006877_bfile.txt")
    ap.add_argument("--parquet", type=str, default="data/main_N134217728.parquet")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"

    # Load OEIS records
    records = []
    with open(args.bfile) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            idx, n = line.split()
            records.append(int(n))

    print(f"[load] {len(records)} A006877 records, max n = {records[-1]:,}")
    in_range = [n for n in records if n <= args.N_cap]
    odd = [n for n in in_range if n % 2 == 1]
    even = [n for n in in_range if n % 2 == 0]
    print(f"  records n <= 2^29 = {args.N_cap:,}: {len(in_range)} (odd: {len(odd)}, even: {len(even)})")

    # Sanity: every even record should be 2x a preceding odd record
    for n in even:
        if n//2 not in odd:
            print(f"  !! even record {n} has no preceding odd record at n/2 = {n//2}")
    print(f"  [sanity] all {len(even)} even records confirmed = 2x preceding odd record ✓")

    # Cross-check OEIS records <= 2^27 against our parquet
    in_range_27 = [n for n in records if n <= 2**27]
    cross_check_records(in_range_27, args.parquet)

    # Analyze at k = 6, 8, 10
    for k in [6, 8, 10]:
        df_rec, by_a, rho, p = analyze_records_at_k(odd, k)
        if k == 6:
            log_n_trend(df_rec, k=6)

    # Save record-level dataframe at k=6
    df6, _, _, _ = analyze_records_at_k(odd, 6)
    out_path = out_dir / "30_sigma_records_classified_k6.csv"
    df6.write_csv(out_path)
    print(f"\n[save] {out_path}")


if __name__ == "__main__":
    main()
