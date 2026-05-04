"""q_sweep_analyze.py — fit beta_p, C_p; classify Pattern alpha/beta/gamma."""
import sys, csv, math
import numpy as np
from collections import defaultdict
sys.stdout.reconfigure(encoding="utf-8")

rows_by_p = defaultdict(list)
with open(r"C:\Collatz\q_sweep_data.csv") as f:
    for row in csv.DictReader(f):
        rows_by_p[int(row["p"])].append(row)

print(f"{'p':>3}  {'n_r':>4}  {'r-range':>8}  {'β (OLS)':>10}  {'σ_β':>7}  {'R²':>7}  "
      f"{'C @ rmax':>9}  {'C mean':>8}")
print("-" * 78)

results = {}
for p in sorted(rows_by_p):
    rows = sorted(rows_by_p[p], key=lambda r: int(r["r"]))
    rs = [int(r["r"]) for r in rows]
    log_N = np.array([math.log(int(r["N"])) for r in rows])
    log_K = np.array([math.log(float(r["K_max_abs"])) for r in rows])
    ratio = np.array([float(r["ratio"]) for r in rows])

    # OLS via polyfit + R²
    n = len(rs)
    if n >= 2:
        coeffs, residuals, rank, sv, rcond = np.polyfit(log_N, log_K, 1, full=True)
        slope, intercept = coeffs
        # R²
        y_pred = slope * log_N + intercept
        ss_res = ((log_K - y_pred) ** 2).sum()
        ss_tot = ((log_K - log_K.mean()) ** 2).sum()
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
        # σ_β (standard error of slope)
        if n > 2:
            sigma2 = ss_res / (n - 2)
            x_var = ((log_N - log_N.mean()) ** 2).sum()
            se_slope = math.sqrt(sigma2 / x_var) if x_var > 0 else float("nan")
        else:
            se_slope = float("nan")
    else:
        slope, r2, se_slope = float("nan"), float("nan"), float("nan")

    C_rmax = ratio[-1]
    C_mean = ratio.mean()
    results[p] = {"beta": slope, "se_beta": se_slope, "R2": r2,
                  "C_rmax": C_rmax, "C_mean": C_mean,
                  "rs": rs, "n": n}
    print(f"{p:>3}  {n:>4}  {rs[0]:>3}..{rs[-1]:<3}  {slope:>10.4f}  {se_slope:>7.4f}  "
          f"{r2:>7.4f}  {C_rmax:>9.3f}  {C_mean:>8.3f}")

print()
betas = [results[p]["beta"] for p in sorted(results)]
beta_min, beta_max = min(betas), max(betas)
print(f"β range across primes: [{beta_min:.4f}, {beta_max:.4f}]  spread = {beta_max-beta_min:.4f}")

C_means = [results[p]["C_mean"] for p in sorted(results)]
C_min, C_max = min(C_means), max(C_means)
print(f"C_mean range: [{C_min:.3f}, {C_max:.3f}]  ratio max/min = {C_max/C_min:.3f}")

C_rmaxes = [results[p]["C_rmax"] for p in sorted(results)]
C_max_min = min(C_rmaxes); C_max_max = max(C_rmaxes)
print(f"C @ rmax range: [{C_max_min:.3f}, {C_max_max:.3f}]  ratio = {C_max_max/C_max_min:.3f}")

# Pattern classification (per brief)
all_beta_in_band = all(0.48 <= b <= 0.52 for b in betas)
beta_close = max(betas) - min(betas) < 0.04
C_spread_ok = (C_max / C_min - 1) <= 0.30   # 30% variation threshold
C_rmax_spread_ok = (C_max_max / C_max_min - 1) <= 0.30

if not all_beta_in_band:
    out_band = [(p, results[p]["beta"]) for p in sorted(results) if not (0.48 <= results[p]["beta"] <= 0.52)]
    print(f"\nPattern γ candidate: primes outside [0.48, 0.52]: {out_band}")
elif not C_spread_ok:
    print(f"\n→ Pattern β (universal β ≈ 0.5; prefactor C varies > 30%)")
else:
    print(f"\n→ Pattern α (universal β AND universal prefactor within 30%)")

print()
print("Per-prime detail:")
for p in sorted(results):
    r = results[p]
    print(f"  p={p:>2}: β={r['beta']:.4f} ± {r['se_beta']:.4f}  R²={r['R2']:.4f}  "
          f"C_rmax={r['C_rmax']:.3f}  C_mean={r['C_mean']:.3f}")
