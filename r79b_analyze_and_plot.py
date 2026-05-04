"""R79b analysis: linear fit of log|K| vs log N, log q; comparison to predicted curves."""
import sys
import os
import math
import numpy as np
import csv

sys.stdout.reconfigure(encoding="utf-8")

CSV_PATH = r"C:\Collatz\r79b_S_partial_data.csv"


def read_data(csv_path: str):
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "r": int(row["r"]),
                "q": int(row["q"]),
                "N": int(row["N"]),
                "K_c1m0": float(row["K_c1m0_abs"]),
                "K_max": float(row["K_max_abs"]),
                "rho": float(row["rho"]),
                "delta_emp": float(row["delta_emp"]),
                "baseline": float(row["random_baseline_abs"]),
                "elapsed": float(row["elapsed_s"]),
            })
    rows.sort(key=lambda r: r["r"])
    return rows


def linear_fit(xs, ys):
    """OLS y = a + b·x. Returns (a, b, R², se_b, n)."""
    n = len(xs)
    if n < 2:
        return None, None, None, None, n
    xs = np.array(xs)
    ys = np.array(ys)
    x_mean = xs.mean()
    y_mean = ys.mean()
    ssxy = ((xs - x_mean) * (ys - y_mean)).sum()
    ssxx = ((xs - x_mean) ** 2).sum()
    if ssxx == 0:
        return None, None, None, None, n
    b = ssxy / ssxx
    a = y_mean - b * x_mean
    y_pred = a + b * xs
    ss_res = ((ys - y_pred) ** 2).sum()
    ss_tot = ((ys - y_mean) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    # Standard error of slope: se_b = sqrt(MSE / SSXX), MSE = ss_res / (n-2)
    if n > 2:
        mse = ss_res / (n - 2)
        se_b = math.sqrt(mse / ssxx)
    else:
        se_b = 0.0
    return a, b, r2, se_b, n


def main():
    rows = read_data(CSV_PATH)
    if len(rows) < 2:
        print("Not enough data yet")
        return
    print(f"# Loaded {len(rows)} rows from {CSV_PATH}")
    print()

    # Headline rate fit: log|K|_max vs log N (so slope = empirical exponent w.r.t. N)
    log_N = [math.log(r["N"]) for r in rows]
    log_q = [math.log(r["q"]) for r in rows]
    log_K_max = [math.log(r["K_max"]) for r in rows]
    log_K_c1 = [math.log(r["K_c1m0"]) for r in rows]
    log_base = [math.log(r["baseline"]) if r["baseline"] > 0 else math.log(1e-9) for r in rows]

    # Fit log|K| ~ a + b · log(N)
    a, b, r2, se_b, n = linear_fit(log_N, log_K_max)
    print(f"# Fit log|K_max| = a + b · log(N), full range r = {rows[0]['r']}..{rows[-1]['r']}")
    print(f"  slope b = {b:.4f} ± {se_b:.4f}  (rate-w.r.t.-N saving)")
    print(f"  intercept a = {a:.4f}  (≈ log of constant: {math.exp(a):.3f})")
    print(f"  R² = {r2:.4f}  (n = {n})")
    print()
    print(f"  Interpretation: |K_max| ≈ {math.exp(a):.3f} · N^{b:.4f}")
    print(f"    rate 0.50 = pure square-root cancellation (target)")
    print(f"    rate 1.00 = trivial (no cancellation)")
    print(f"    rate < 0.50 = sub-square-root (Bourgain-Konyagin level)")
    print()

    # High-r fit: r ≥ max-5
    high_rows = rows[-min(6, len(rows)):]
    if len(high_rows) >= 3:
        a_h, b_h, r2_h, se_b_h, n_h = linear_fit(
            [math.log(r["N"]) for r in high_rows],
            [math.log(r["K_max"]) for r in high_rows]
        )
        print(f"# High-r fit (last {n_h} points, r = {high_rows[0]['r']}..{high_rows[-1]['r']})")
        print(f"  slope b = {b_h:.4f} ± {se_b_h:.4f}")
        print(f"  R² = {r2_h:.4f}")
        print()

    # Brief's δ-mapping (against q^{1/2})
    print("# Brief's δ-saving against q^{1/2}:")
    print(f"  ρ(r) := log|K_max|/log(q^{{1/2}}) − 1, δ_brief = −ρ/2")
    print()
    print(f"  {'r':>3} {'|K|_max':>12} {'log|K|/log√N':>14} {'δ_brief':>10} {'δ_pred_sqrt':>13}")
    for row in rows:
        rate_N = math.log(row["K_max"]) / math.log(row["N"])
        delta_brief = row["delta_emp"]
        delta_pred_sqrt = (math.log(row["q"]) - 2 * math.log(row["K_max"])) / (2 * math.log(row["q"]))
        print(f"  {row['r']:>3} {row['K_max']:>12.2f} {rate_N:>14.4f} {delta_brief:>10.4f} {delta_pred_sqrt:>13.4f}")

    # Fit log|K| ~ a + b · log(q^{1/2})
    print()
    log_sqrt_q = [0.5 * lq for lq in log_q]
    a2, b2, r22, se_b2, _ = linear_fit(log_sqrt_q, log_K_max)
    print(f"# Fit log|K_max| = a + b · log(q^{{1/2}}):")
    print(f"  slope = {b2:.4f} ± {se_b2:.4f}  (b = 1 means |K| ∝ q^{{1/2}})")
    print(f"  intercept = {a2:.4f}  (constant ≈ {math.exp(a2):.3f})")
    print(f"  R² = {r22:.4f}")
    print(f"  Brief's δ = (1 − slope)/2 = {(1 - b2)/2:.4f} ± {se_b2/2:.4f}")
    print()

    # Random-baseline check
    print("# C3: Random-phase baseline comparison")
    print(f"  Expected: random sum has |.|² ≈ N, so |.| ≈ √N. Compare to actual |K_max|.")
    print()
    print(f"  {'r':>3} {'√N':>10} {'|K|_max':>10} {'baseline':>10} {'|K|/√N':>10} {'base/√N':>10}")
    for row in rows:
        sqrt_N = math.sqrt(row["N"])
        ratio_K = row["K_max"] / sqrt_N
        ratio_base = row["baseline"] / sqrt_N
        print(f"  {row['r']:>3} {sqrt_N:>10.2f} {row['K_max']:>10.2f} {row['baseline']:>10.2f} {ratio_K:>10.3f} {ratio_base:>10.3f}")

    print()
    print("# C2: r mod 3 oscillation check")
    print()
    print(f"  {'r':>3} {'r mod 3':>8} {'|K|/√N':>10}")
    for row in rows:
        sqrt_N = math.sqrt(row["N"])
        ratio = row["K_max"] / sqrt_N
        print(f"  {row['r']:>3} {row['r'] % 3:>8} {ratio:>10.4f}")

    by_mod = {0: [], 1: [], 2: []}
    for row in rows:
        by_mod[row['r'] % 3].append(row["K_max"] / math.sqrt(row["N"]))
    print()
    print("  per-class average |K|/√N:")
    for k in [0, 1, 2]:
        if by_mod[k]:
            print(f"    r ≡ {k} mod 3: mean = {np.mean(by_mod[k]):.3f}, std = {np.std(by_mod[k]):.3f}, n = {len(by_mod[k])}")


if __name__ == "__main__":
    main()
