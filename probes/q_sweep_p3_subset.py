"""Sanity check: how does p=3 β depend on r-window? Compare full r=8..20 vs r=12..20."""
import csv, math
import numpy as np

with open(r"C:\Collatz\q_sweep_data.csv") as f:
    p3 = [r for r in csv.DictReader(f) if int(r["p"]) == 3]

p3 = sorted(p3, key=lambda x: int(x["r"]))
log_N = np.array([math.log(int(r["N"])) for r in p3])
log_K = np.array([math.log(float(r["K_max_abs"])) for r in p3])
rs = [int(r["r"]) for r in p3]

def fit(idx, label):
    if len(idx) < 2: return
    x = log_N[idx]; y = log_K[idx]
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = ((y - y_pred) ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot
    print(f"  {label}: r in {[rs[i] for i in idx]}  n={len(idx)}  beta={slope:.4f}  R2={r2:.4f}")

fit(list(range(len(rs))), "full r=8..20")
fit([i for i, r in enumerate(rs) if r >= 10], "r=10..20")
fit([i for i, r in enumerate(rs) if r >= 12], "r=12..20")
fit([i for i, r in enumerate(rs) if r >= 14], "r=14..20")
fit([i for i, r in enumerate(rs) if r >= 16], "r=16..20")
