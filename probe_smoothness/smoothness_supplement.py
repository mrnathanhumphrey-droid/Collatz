"""Supplementary analysis: correlation of pi(r) with max_prime(r), extremes."""
import csv
import math
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
from scipy.stats import pearsonr, spearmanr

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = Path("C:/Collatz/probe_smoothness")

# Load partition CSV
data_by_k = defaultdict(list)
with open(OUTDIR / "result_smooth_rough_partition.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        data_by_k[int(row["k"])].append({
            "r": int(row["r"]),
            "max_prime": int(row["max_prime"]),
            "orbit_length": int(row["orbit_length"]),
            "pi": float(row["pi"]),
        })

print(f"{'k':>3} {'n':>5} {'corr(pi, max_p)':>17} {'corr(pi, log_max_p)':>20} "
      f"{'spearman(pi, max_p)':>20} {'corr(pi^2, max_p)':>19}")
print("-" * 95)
for k, rows in sorted(data_by_k.items()):
    pis = np.array([r["pi"] for r in rows])
    pis_sq = pis ** 2
    mps = np.array([r["max_prime"] for r in rows], dtype=float)
    log_mps = np.log10(np.maximum(mps, 1.0))
    pearson_pi_mp = pearsonr(pis, mps).statistic
    pearson_pi_log = pearsonr(pis, log_mps).statistic
    spearman_pi_mp = spearmanr(pis, mps).statistic
    pearson_pi2_mp = pearsonr(pis_sq, mps).statistic
    print(f"{k:>3} {len(rows):>5} {pearson_pi_mp:>+17.4f} {pearson_pi_log:>+20.4f} "
          f"{spearman_pi_mp:>+20.4f} {pearson_pi2_mp:>+19.4f}")

# Top-10 highest pi residues per k, with their max_prime
print()
print("Top-10 highest pi residues per k:")
for k, rows in sorted(data_by_k.items()):
    sorted_rows = sorted(rows, key=lambda r: -r["pi"])
    print(f"  k={k}:")
    for r in sorted_rows[:10]:
        print(f"    r={r['r']:>5}  pi={r['pi']:.4e}  max_prime={r['max_prime']:>7}  orbit_len={r['orbit_length']}")

print()
print("Top-10 lowest pi residues per k:")
for k, rows in sorted(data_by_k.items()):
    sorted_rows = sorted(rows, key=lambda r: r["pi"])
    print(f"  k={k}:")
    for r in sorted_rows[:10]:
        print(f"    r={r['r']:>5}  pi={r['pi']:.4e}  max_prime={r['max_prime']:>7}  orbit_len={r['orbit_length']}")

# Mass ratio rough explicit
print()
print(f"{'k':>3} {'B':>6} {'count_share_rough':>18} {'mass_share_rough':>17} {'mass_ratio_rough':>17}")
print("-" * 70)
with open(OUTDIR / "result_S_k_conditional.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        k = int(row["k"])
        B = int(row["B"])
        cs_smooth = float(row["count_share_smooth"])
        ms_smooth = float(row["mass_share_smooth"])
        cs_rough = 1.0 - cs_smooth
        ms_rough = 1.0 - ms_smooth
        mr_rough = ms_rough / cs_rough if cs_rough > 0 else float('nan')
        print(f"{k:>3} {B:>6} {cs_rough:>18.4f} {ms_rough:>17.4f} {mr_rough:>17.4f}")
