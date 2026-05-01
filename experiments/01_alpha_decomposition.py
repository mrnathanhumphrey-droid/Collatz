"""
Experiment 01 — Alpha decomposition

For each odd residue r mod 2^k, compute the deterministic Collatz prefix:
track state = a*m + c symbolically, applying Collatz rules whenever a is even
(parity of state forced by c alone). Stop when a becomes odd (parity now
depends on m). Record (prefix_steps, a_final, c_final).

Predict per-class intercept:
    alpha_det(r) = prefix_steps(r) + (3/(ln4-ln3)) * ln(a_final(r) / 2^k)

Compare to per-class OLS intercept of sigma vs ln(n) on data restricted to
n in residue class r mod 2^k.

Output: figure (predicted vs actual alpha) and R^2 + per-class residuals.

Usage:
    python 01_alpha_decomposition.py --N 33554432 --k 6
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")

LOG_FACTOR_ODD = 3.0 / (np.log(4.0) - np.log(3.0))   # ~10.4282


def deterministic_prefix(r, a0, max_steps=400):
    """Return (steps, a_final, c_final) for the symbolic prefix from state (a=a0, c=r)."""
    a, c = a0, r
    steps = 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3 * c + 1
        steps += 1
    return steps, a, c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=1 << 25)
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--data", type=str, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data) if args.data else here.parent / "data"
    out_dir = Path(args.out) if args.out else here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    M = 1 << args.k
    K = M // 2

    # Load and filter to odd n
    df = pl.read_parquet(data_dir / f"main_N{args.N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n = df["n"].to_numpy().astype(np.int64)
    log_n = np.log(n.astype(np.float64))
    sigma = df["sigma"].to_numpy().astype(np.float64)
    res = (n % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)
    print(f"[load] N={args.N:,}, odd-only={len(n):,}, k={args.k}, K={K}")

    # Per-class OLS
    alpha_actual = np.zeros(K)
    beta_actual = np.zeros(K)
    se_alpha = np.zeros(K)
    n_per = np.zeros(K, dtype=np.int64)
    for kk in range(K):
        m = class_idx == kk
        if m.sum() < 50:
            continue
        bk, ak = np.polyfit(log_n[m], sigma[m], 1)
        rk = sigma[m] - (ak + bk * log_n[m])
        alpha_actual[kk] = ak
        beta_actual[kk] = bk
        n_per[kk] = m.sum()
        var_x = log_n[m].var()
        mean_x = log_n[m].mean()
        se_alpha[kk] = rk.std() * np.sqrt(1.0/n_per[kk] + mean_x**2 / (n_per[kk] * var_x))

    # Deterministic prefix predictions
    alpha_pred = np.zeros(K)
    prefix_arr = np.zeros(K, dtype=int)
    a_final_arr = np.zeros(K, dtype=int)
    c_final_arr = np.zeros(K, dtype=int)
    for kk in range(K):
        r = 2 * kk + 1
        steps, a_f, c_f = deterministic_prefix(r, M)
        alpha_pred[kk] = steps + LOG_FACTOR_ODD * np.log(a_f / float(M))
        prefix_arr[kk] = steps
        a_final_arr[kk] = a_f
        c_final_arr[kk] = c_f

    # Linear fit
    b_fit, a_fit = np.polyfit(alpha_pred, alpha_actual, 1)
    fitted = a_fit + b_fit * alpha_pred
    resid = alpha_actual - fitted
    r_squared = 1 - np.sum(resid**2) / np.sum((alpha_actual - alpha_actual.mean())**2)

    print(f"[fit] alpha_actual = {a_fit:.4f} + {b_fit:.4f} * alpha_pred")
    print(f"[fit] R^2 = {r_squared:.6f}")
    print(f"[fit] SD(alpha_actual) = {alpha_actual.std():.3f}")
    print(f"[fit] SD(alpha_stoch residual) = {resid.std():.3f}")
    print(f"[fit] mean(per-class SE alpha) = {se_alpha.mean():.3f}")
    print(f"[fit] ratio SD(resid)/mean(SE) = {resid.std()/se_alpha.mean():.3f}  "
          f"(>>1 = real residual structure, ~1 = at noise floor)")
    print(f"[fit] unique a_final values: {sorted(set(a_final_arr.tolist()))}")

    # Figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax = axes[0]
    ax.errorbar(alpha_pred, alpha_actual, yerr=1.96*se_alpha, fmt="o", markersize=5,
                lw=0.8, color="C0", alpha=0.8)
    xs = np.linspace(alpha_pred.min(), alpha_pred.max(), 100)
    ax.plot(xs, a_fit + b_fit * xs, "r--", lw=1.2,
            label=f"fit: y = {a_fit:.2f} + {b_fit:.3f}x  (R^2={r_squared:.6f})")
    ax.set_xlabel("predicted alpha = prefix_steps + 10.4282 * ln(a_final / 2^k)")
    ax.set_ylabel("alpha actual (per-class OLS intercept)")
    ax.set_title(f"Alpha decomposition  (N={args.N:,}, k={args.k}, K={K})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.errorbar(np.arange(K), resid, yerr=1.96*se_alpha, fmt="o", markersize=4,
                lw=0.6, color="C0", alpha=0.8)
    ax.axhline(0, color="red", linestyle="--", alpha=0.6)
    ax.set_xlabel("class index k (residue 2k+1 mod 2^k)")
    ax.set_ylabel("alpha_stoch (actual - prediction)")
    ax.set_title(f"Per-class residual after prefix prediction\n"
                 f"SD(resid)={resid.std():.3f}  ratio to noise={resid.std()/se_alpha.mean():.3f}")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = out_dir / f"01_alpha_decomposition_k{args.k}_N{args.N}.png"
    plt.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"[save] {out_png}")

    # CSV output
    df_out = pl.DataFrame({
        "class_idx": np.arange(K),
        "residue": np.arange(1, M, 2),
        "prefix_steps": prefix_arr,
        "a_final": a_final_arr,
        "c_final": c_final_arr,
        "alpha_pred": alpha_pred,
        "alpha_actual": alpha_actual,
        "alpha_resid": resid,
        "se_alpha": se_alpha,
        "n_per_class": n_per,
    })
    out_csv = out_dir / f"01_alpha_decomposition_k{args.k}_N{args.N}.csv"
    df_out.write_csv(out_csv)
    print(f"[save] {out_csv}")


if __name__ == "__main__":
    main()
