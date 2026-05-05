"""
esscher_tilt_r58_closure_v2_log_ratio.py
========================================
R58 closure attempt #2: Esscher tilt at log(R58 weight / R60 weight) per R69's
mechanism identification.

Prior attempt (esscher_tilt_r58_closure.md) tilted by σ_orbit (depth) and
failed: best λ ≈ -0.01, improvement only +0.014. Uniform depth-tilt cannot
fix opposite-sign residuals at r=5/r=23 (over-predicted enhancement) vs
r=13 (over-predicted depletion).

This attack uses R69's mechanism: R58 weights inverse-tree nodes by subtree
size; R60 weights by survivor-conditioned forward orbit visits (size-stratified
Markov stationary on (residue, log-size)). The closure observable is per-node
log-ratio of the two weightings.

Geometric-interpolation framing (numerically stable):
  w_λ(m) = w_R58(m)^(1-λ) · w_R60(m)^λ
  λ = 0 -> R58, λ = 1 -> R60, intermediate λ -> hybrid

Per-node R60 weight (proxy):
  w_R60(m) = π_R60(r, b) / N_cell(r, b)
  where (r, b) = (m mod 32, floor(log_2 m)) is m's size-stratified Markov
  cell, π_R60(r, b) is the Markov stationary mass on that cell, and
  N_cell(r, b) is the count of inverse-tree nodes in that cell. Normalization
  ensures D_1(r) = R60's residue-marginal prediction by construction.

Output:
  esscher_tilt_v2_log_ratio.md
  esscher_tilt_v2_log_ratio_residuals.csv  -- per-r at λ=0 vs λ=best
  esscher_tilt_v2_log_ratio_lambda_sweep.csv
  esscher_tilt_v2_log_ratio_perclass_logratio.csv  -- mean/std/skew of r(m) per residue
"""
from __future__ import annotations

import csv
import io
import math
import sys
import time
from collections import defaultdict, deque
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

OUT = Path(r"C:\Collatz")


# --------------------------------------------------------------------------- #
# Inverse tree + subtree sizes (R58)                                          #
# --------------------------------------------------------------------------- #

def build_inverse_tree_from_one(max_value):
    """BFS-build inverse Syracuse tree from m=1 up to max_value. Odd integers only."""
    tree = {1: {'parent': None, 'depth': 0}}
    q = deque([1])
    while q:
        m = q.popleft()
        d = tree[m]['depth']
        if m % 3 == 0:
            continue
        v_start = 2 if (m % 3 == 1) else 1
        for v in range(v_start, 64, 2):
            num = m * (1 << v) - 1
            if num <= 0:
                continue
            if num % 3 != 0:
                continue
            pred = num // 3
            if pred > max_value:
                break
            if pred & 1 == 0:
                continue
            if pred == m:
                continue
            if pred not in tree:
                tree[pred] = {'parent': m, 'depth': d + 1}
                q.append(pred)
    return tree


def compute_subtree_sizes(tree):
    children = defaultdict(list)
    for m, info in tree.items():
        if info['parent'] is not None:
            children[info['parent']].append(m)
    by_depth = sorted(tree.keys(), key=lambda m: -tree[m]['depth'])
    size = {m: 1 for m in tree}
    for m in by_depth:
        for c in children[m]:
            size[m] += size[c]
    return size


# --------------------------------------------------------------------------- #
# R60 size-stratified per-(r, b) stationary loader                            #
# --------------------------------------------------------------------------- #

def load_r60_stationary(path=r"C:\Collatz\size_stratified_eigvec.csv"):
    """Return dict (r, b) -> v_PF (Markov stationary mass)."""
    pi = {}
    with open(path, encoding="utf-8") as fh:
        r = csv.DictReader(fh)
        for row in r:
            rr = int(row["r"])
            bb = int(row["b"])
            v = float(row["v_PF"])
            pi[(rr, bb)] = v
    return pi


# --------------------------------------------------------------------------- #
# D_emp loader                                                                 #
# --------------------------------------------------------------------------- #

def load_empirical_D(t=90):
    path = OUT / "experiments_output" / "chang_qsd_test.csv"
    with open(path, encoding="utf-8") as fh:
        r = csv.DictReader(fh)
        for row in r:
            if int(row["t"]) == t:
                D = {}
                for k, v in row.items():
                    if k.startswith("D_r"):
                        D[int(k[3:])] = float(v)
                return D
    return None


# --------------------------------------------------------------------------- #
# Aggregation                                                                  #
# --------------------------------------------------------------------------- #

def predicted_D_mod32(weights):
    """For each odd r in {1, 3, ..., 31}, sum weights[m] over m ≡ r (mod 32),
    normalize so mean over 16 odd residues = 1."""
    by_r = defaultdict(float)
    total = 0.0
    for m, w in weights.items():
        if m & 1:
            by_r[m % 32] += w
            total += w
    if total == 0:
        return {r: 0.0 for r in range(1, 32, 2)}
    n_residues = 16
    return {r: by_r[r] / total * n_residues for r in range(1, 32, 2)}


def pearson(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mx) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (den_x * den_y) if den_x * den_y > 0 else 0.0


def mae(xs, ys):
    return sum(abs(x - y) for x, y in zip(xs, ys)) / len(xs)


# --------------------------------------------------------------------------- #
# Main                                                                         #
# --------------------------------------------------------------------------- #

def make_w_R60_per_node(tree, pi_R60_by_rb, B=64):
    """For each odd m in tree, w_R60(m) = π_R60(r, b) / N_cell(r, b).
    (r, b) = (m mod 32, min(floor(log2 m), B-1)). Cells with π = 0 get
    a tiny floor (1e-30) to avoid log(0); those nodes contribute negligibly
    to the geometric interpolation at any λ.
    """
    cell_count = defaultdict(int)
    for m in tree:
        if not (m & 1):
            continue
        r = m % 32
        b = min(int(math.log2(m)) if m > 0 else 0, B - 1)
        cell_count[(r, b)] += 1

    w_R60 = {}
    floor = 1e-30
    for m in tree:
        if not (m & 1):
            continue
        r = m % 32
        b = min(int(math.log2(m)) if m > 0 else 0, B - 1)
        pi_rb = pi_R60_by_rb.get((r, b), 0.0)
        n_rb = cell_count[(r, b)]
        if pi_rb > 0 and n_rb > 0:
            w_R60[m] = pi_rb / n_rb
        else:
            w_R60[m] = floor
    return w_R60, cell_count


def run_at_N(N, pi_R60, label):
    print(f"\n=== {label}: build inverse tree at N={N:,} ===", flush=True)
    t0 = time.time()
    tree = build_inverse_tree_from_one(N)
    print(f"  tree size: {len(tree):,} odd nodes  (built in {time.time() - t0:.1f}s)", flush=True)

    sizes = compute_subtree_sizes(tree)
    w_R58 = {m: float(sizes[m]) for m in tree if m & 1}
    w_R60, cell_count = make_w_R60_per_node(tree, pi_R60)

    # Log-ratio per node
    log_ratio = {}
    for m in w_R58:
        a = w_R58[m]
        b = w_R60[m]
        if a > 0 and b > 0:
            log_ratio[m] = math.log(a) - math.log(b)
        else:
            log_ratio[m] = 0.0

    # D_emp at t=90
    D_emp = load_empirical_D(t=90)
    odd_residues = list(range(1, 32, 2))
    emp_xs = [D_emp[r] for r in odd_residues]

    # Baseline R58 (λ=0): D = predicted_D_mod32(w_R58)
    D_R58 = predicted_D_mod32(w_R58)
    r58_xs = [D_R58[r] for r in odd_residues]
    pearson_R58 = pearson(r58_xs, emp_xs)
    mae_R58 = mae(r58_xs, emp_xs)
    print(f"  Baseline R58 Pearson: {pearson_R58:+.4f}  MAE: {mae_R58:.4f}", flush=True)

    # Test endpoint R60 (λ=1): D = predicted_D_mod32(w_R60)
    D_R60 = predicted_D_mod32(w_R60)
    r60_xs = [D_R60[r] for r in odd_residues]
    pearson_R60 = pearson(r60_xs, emp_xs)
    mae_R60 = mae(r60_xs, emp_xs)
    print(f"  Pure R60 (λ=1) Pearson: {pearson_R60:+.4f}  MAE: {mae_R60:.4f}", flush=True)

    return tree, sizes, w_R58, w_R60, log_ratio, D_emp, D_R58, D_R60


def lambda_sweep(w_R58, w_R60, D_emp, lambdas):
    """Geometric interpolation w_λ(m) = w_R58(m)^(1-λ) · w_R60(m)^λ. Pearson per λ."""
    odd_residues = list(range(1, 32, 2))
    emp_xs = [D_emp[r] for r in odd_residues]
    rows = []
    for lam in lambdas:
        # Compute w_λ in log space for numerical stability
        # w_λ(m) = exp((1-λ) log w_R58(m) + λ log w_R60(m))
        w_lam = {}
        log_max = -1e300
        for m in w_R58:
            a = w_R58[m]; b = w_R60[m]
            if a <= 0 or b <= 0:
                continue
            lw = (1.0 - lam) * math.log(a) + lam * math.log(b)
            w_lam[m] = lw
            if lw > log_max:
                log_max = lw
        # Subtract log_max for numerical stability and exponentiate
        w_lam_actual = {m: math.exp(lw - log_max) for m, lw in w_lam.items()}
        D_lam = predicted_D_mod32(w_lam_actual)
        lam_xs = [D_lam[r] for r in odd_residues]
        p = pearson(lam_xs, emp_xs)
        m_ = mae(lam_xs, emp_xs)
        rows.append({"lambda": lam, "pearson": p, "mae": m_, "D_lambda": D_lam})
    return rows


def per_residue_log_ratio_stats(log_ratio, w_R58, tree):
    """For each odd r in {1,...,31}, distribution stats of log_ratio over m ≡ r (mod 32)."""
    by_r = defaultdict(list)
    for m, lr in log_ratio.items():
        by_r[m % 32].append(lr)
    out = {}
    for r in range(1, 32, 2):
        vals = np.array(by_r[r], dtype=np.float64)
        if len(vals) == 0:
            out[r] = {"n": 0, "mean": np.nan, "std": np.nan, "skew": np.nan}
            continue
        m = vals.mean()
        s = vals.std()
        sk = ((vals - m) ** 3).mean() / (s ** 3) if s > 0 else 0.0
        out[r] = {"n": len(vals), "mean": float(m), "std": float(s), "skew": float(sk)}
    return out


def main():
    print("=" * 78)
    print("R58 closure attempt #2: Esscher tilt at log(R58/R60) per R69 mechanism")
    print("=" * 78)

    # Load R60 stationary
    pi_R60 = load_r60_stationary()
    print(f"\n  Loaded R60 stationary: {len(pi_R60)} (r, b) cells")
    nonzero = sum(1 for v in pi_R60.values() if v > 0)
    print(f"  Nonzero π cells: {nonzero}")

    # Step 1+2: Build at N=2^22 and reproduce R58 baseline
    N = 1 << 22
    tree, sizes, w_R58, w_R60, log_ratio, D_emp, D_R58, D_R60 = run_at_N(N, pi_R60, "Step 1+2 (test)")

    # Step 3: Per-residue log-ratio distribution
    print(f"\n=== Step 3: Per-residue log-ratio r(m) = log(w_R58/w_R60) ===")
    stats = per_residue_log_ratio_stats(log_ratio, w_R58, tree)
    print(f"  {'r':>3}  {'n':>8}  {'mean':>9}  {'std':>9}  {'skew':>9}  D_emp - D_R58")
    odd_residues = list(range(1, 32, 2))
    for r in odd_residues:
        s = stats[r]
        residual = D_emp[r] - D_R58[r]
        print(f"  {r:>3}  {s['n']:>8}  {s['mean']:>+9.3f}  {s['std']:>9.3f}  "
              f"{s['skew']:>+9.3f}  {residual:>+12.4f}")

    # Step 4: Lambda sweep
    print(f"\n=== Step 4: λ sweep on geometric interpolation ===")
    lambdas = np.arange(0.0, 1.01, 0.05).tolist()
    sweep = lambda_sweep(w_R58, w_R60, D_emp, lambdas)
    print(f"  {'λ':>5}  {'Pearson':>10}  {'MAE':>8}")
    for row in sweep:
        print(f"  {row['lambda']:>5.2f}  {row['pearson']:>+10.4f}  {row['mae']:>8.4f}")
    best = max(sweep, key=lambda x: x["pearson"])
    print(f"\n  Best λ: {best['lambda']:.2f}  Pearson: {best['pearson']:+.4f}  MAE: {best['mae']:.4f}")
    D_best = best["D_lambda"]

    # Step 5: Per-residue check at best λ
    print(f"\n=== Step 5: Per-residue residuals at λ_best = {best['lambda']:.2f} ===")
    print(f"  {'r':>3}  {'D_emp':>8}  {'D_R58':>8}  {'D_R60':>8}  "
          f"{'D_best':>8}  {'res_R58':>8}  {'res_best':>9}  change")
    for r in odd_residues:
        e = D_emp[r]; r58 = D_R58[r]; r60 = D_R60[r]; b = D_best[r]
        res_R58 = e - r58
        res_best = e - b
        change = abs(res_R58) - abs(res_best)
        marker = " ← QSD-extreme" if r in {5, 13, 23} else ""
        print(f"  {r:>3}  {e:>8.4f}  {r58:>8.4f}  {r60:>8.4f}  "
              f"{b:>8.4f}  {res_R58:>+8.4f}  {res_best:>+9.4f}  {change:>+8.4f}{marker}")

    # Step 6: Train-test
    print(f"\n=== Step 6: Train-test (fit at 2^21, test at 2^22) ===")
    N_train = 1 << 21
    tree_tr, sizes_tr, w_R58_tr, w_R60_tr, lr_tr, D_emp_tr, D_R58_tr, D_R60_tr = run_at_N(N_train, pi_R60, "train")
    sweep_tr = lambda_sweep(w_R58_tr, w_R60_tr, D_emp_tr, lambdas)
    best_tr = max(sweep_tr, key=lambda x: x["pearson"])
    print(f"  Best λ_train: {best_tr['lambda']:.2f}  Pearson_train: {best_tr['pearson']:+.4f}")

    # Apply λ_train to N=2^22 (already built above)
    sweep_tt = lambda_sweep(w_R58, w_R60, D_emp, [best_tr["lambda"]])
    test_at_lam_train = sweep_tt[0]
    print(f"  Apply λ_train at N=2^22: Pearson_test = {test_at_lam_train['pearson']:+.4f}  "
          f"MAE = {test_at_lam_train['mae']:.4f}")
    print(f"  vs. λ optimized at N=2^22: Pearson = {best['pearson']:+.4f}")
    print(f"  Train-test gap: {abs(test_at_lam_train['pearson'] - best['pearson']):.4f}")

    # CSV outputs
    out_dir = OUT
    with open(out_dir / "esscher_tilt_v2_log_ratio_lambda_sweep.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["lambda", "pearson", "mae"])
        for row in sweep:
            w.writerow([row["lambda"], row["pearson"], row["mae"]])

    with open(out_dir / "esscher_tilt_v2_log_ratio_residuals.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["r", "D_emp", "D_R58", "D_R60", "D_best", "res_R58", "res_best", "is_qsd_extreme"])
        for r in odd_residues:
            e = D_emp[r]; r58 = D_R58[r]; r60 = D_R60[r]; b = D_best[r]
            w.writerow([r, e, r58, r60, b, e - r58, e - b, r in {5, 13, 23}])

    with open(out_dir / "esscher_tilt_v2_log_ratio_perclass_logratio.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["r", "n", "mean_log_ratio", "std_log_ratio", "skew_log_ratio"])
        for r in odd_residues:
            s = stats[r]
            w.writerow([r, s["n"], s["mean"], s["std"], s["skew"]])

    print(f"\n[csv outputs written to {out_dir}]")


if __name__ == "__main__":
    main()
