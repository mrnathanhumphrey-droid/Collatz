"""β estimates using K_c1m0 (canonical (c=1, m=0)) — eliminates (c,m)-resonance confound."""
import sys, csv, math
import numpy as np
from collections import defaultdict
sys.stdout.reconfigure(encoding="utf-8")

rows_by_p = defaultdict(list)
with open(r"C:\Collatz\q_sweep_data.csv") as f:
    for row in csv.DictReader(f):
        rows_by_p[int(row["p"])].append(row)

print(f"{'p':>3}  {'beta_max (c,m sweep)':>22}  {'beta_canonical (c=1,m=0)':>26}  {'diff':>8}")
print("-" * 72)

results = {}
for p in sorted(rows_by_p):
    rows = sorted(rows_by_p[p], key=lambda r: int(r["r"]))
    log_N = np.array([math.log(int(r["N"])) for r in rows])
    log_K_max = np.array([math.log(float(r["K_max_abs"])) for r in rows])
    log_K_c1m0 = np.array([math.log(float(r["K_c1m0_abs"])) for r in rows])

    b_max, _ = np.polyfit(log_N, log_K_max, 1)
    b_c1m0, _ = np.polyfit(log_N, log_K_c1m0, 1)

    print(f"{p:>3}  {b_max:>22.4f}  {b_c1m0:>26.4f}  {b_max-b_c1m0:>+8.4f}")
    results[p] = (b_max, b_c1m0)

print()
b_canonicals = [results[p][1] for p in sorted(results)]
print(f"β_canonical range: [{min(b_canonicals):.4f}, {max(b_canonicals):.4f}]  spread = {max(b_canonicals)-min(b_canonicals):.4f}")

# Per-prime canonical prefactor C_canonical = K_c1m0 / sqrt(N)
print()
print(f"{'p':>3}  {'C_c1m0 mean':>13}  {'C_c1m0 @ rmax':>15}")
for p in sorted(rows_by_p):
    rows = sorted(rows_by_p[p], key=lambda r: int(r["r"]))
    cs = [float(r["K_c1m0_abs"]) / math.sqrt(int(r["N"])) for r in rows]
    print(f"{p:>3}  {np.mean(cs):>13.4f}  {cs[-1]:>15.4f}")
