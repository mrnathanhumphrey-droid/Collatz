"""
Experiment 24 — alpha-decomposition k-sweep at N = 2^27.

Existing writeup table covers k in {6, 7, 8, 9}. The reviewer noted that k = 4
and k = 5 were skipped (small but non-trivial gap), and that k in {10, 11, 12}
extends the universality claim further if per-class N stays adequate.

Sweep: for k in {4, 5, 6, 7, 8, 9, 10, 11, 12}, compute on the same N = 2^27
parquet:
  - per-class OLS alpha_actual via sigma vs ln(n)
  - per-class deterministic-prefix prediction alpha_det
  - linear fit alpha_actual ~ a + b * alpha_det, R^2
  - noise-floor ratio: SD(residual) / mean(per-class SE on alpha)

The expected signature for the universality claim: noise-floor ratio in
[0.91, 0.99] band across all k, with no breakaway at higher k.

Usage:
    python 24_k_sweep_alpha_decomposition.py --N 134217728
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

sys.stdout.reconfigure(encoding="utf-8")

LOG_FACTOR_ODD = 3.0 / (np.log(4.0) - np.log(3.0))


def deterministic_prefix(r, a0, max_steps=400):
    a, c = a0, r
    steps = 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3 * c + 1
        steps += 1
    return steps, a, c


def compute_table_for_k(n, log_n, sigma, k):
    """Return (R^2, SD_resid, mean_SE, ratio, n_classes_with_data, min_per_class)."""
    M = 1 << k
    K = M // 2
    res = (n % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)

    alpha_actual = np.zeros(K)
    se_alpha = np.zeros(K)
    n_per = np.zeros(K, dtype=np.int64)

    # Build per-class OLS via groupby-like aggregation. To keep memory bounded
    # at high k, we loop classes; numpy fancy-indexing per class is fine
    # since class_idx is small int32.
    sort_idx = np.argsort(class_idx)
    class_idx_s = class_idx[sort_idx]
    log_n_s = log_n[sort_idx]
    sigma_s = sigma[sort_idx]

    # Find run starts
    boundaries = np.concatenate([[0], np.where(np.diff(class_idx_s) != 0)[0] + 1, [len(class_idx_s)]])

    have_data = np.zeros(K, dtype=bool)
    for bi in range(len(boundaries) - 1):
        start, end = boundaries[bi], boundaries[bi + 1]
        kk = int(class_idx_s[start])
        cnt = end - start
        if cnt < 50:
            continue
        log_x = log_n_s[start:end]
        y = sigma_s[start:end]
        bk, ak = np.polyfit(log_x, y, 1)
        rk = y - (ak + bk * log_x)
        alpha_actual[kk] = ak
        n_per[kk] = cnt
        var_x = log_x.var()
        mean_x = log_x.mean()
        se_alpha[kk] = rk.std() * np.sqrt(1.0 / cnt + mean_x ** 2 / (cnt * var_x))
        have_data[kk] = True

    # Deterministic prefix predictions
    alpha_pred = np.zeros(K)
    a_final_arr = np.zeros(K, dtype=np.int64)
    for kk in range(K):
        r = 2 * kk + 1
        steps, a_f, _ = deterministic_prefix(r, M)
        alpha_pred[kk] = steps + LOG_FACTOR_ODD * np.log(a_f / float(M))
        a_final_arr[kk] = a_f

    # Filter to classes that had data
    mask = have_data
    if mask.sum() < 4:
        return None
    ap = alpha_pred[mask]
    aa = alpha_actual[mask]
    se = se_alpha[mask]
    np_per = n_per[mask]

    b_fit, a_fit = np.polyfit(ap, aa, 1)
    fitted = a_fit + b_fit * ap
    resid = aa - fitted
    sst = np.sum((aa - aa.mean()) ** 2)
    r2 = 1.0 - np.sum(resid ** 2) / sst
    sd_resid = float(resid.std())
    mean_se = float(se.mean())
    ratio = sd_resid / mean_se if mean_se > 0 else float("nan")
    n_unique_a_final = len(set(a_final_arr[mask].tolist()))

    return dict(
        k=k,
        mod=M,
        n_classes=int(mask.sum()),
        K_total=K,
        n_a_final=n_unique_a_final,
        per_class_min=int(np_per.min()),
        per_class_median=int(np.median(np_per)),
        per_class_max=int(np_per.max()),
        r2=r2,
        sd_resid=sd_resid,
        mean_se=mean_se,
        ratio=ratio,
        b_fit=float(b_fit),
        a_fit=float(a_fit),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=1 << 27)
    ap.add_argument("--ks", type=str, default="4,5,6,7,8,9,10,11,12")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] N={args.N:,}", flush=True)
    t0 = time.perf_counter()
    df = pl.read_parquet(data_dir / f"main_N{args.N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n = df["n"].to_numpy().astype(np.int64)
    log_n = np.log(n.astype(np.float64))
    sigma = df["sigma"].to_numpy().astype(np.float64)
    print(f"        odd rows: {len(n):,}  load+convert: {time.perf_counter()-t0:.1f}s",
          flush=True)
    print()

    ks = [int(k) for k in args.ks.split(",")]

    print(f"{'k':>3} {'mod':>6} {'classes':>8} {'a_final_lvls':>13} "
          f"{'per_class_min':>14} {'per_class_med':>14} "
          f"{'R^2':>9} {'SD_resid':>10} {'mean_SE':>10} {'ratio':>7}")
    rows = []
    for k in ks:
        t0 = time.perf_counter()
        r = compute_table_for_k(n, log_n, sigma, k)
        if r is None:
            continue
        rows.append(r)
        print(f"{k:>3} {r['mod']:>6} {r['n_classes']:>8} {r['n_a_final']:>13} "
              f"{r['per_class_min']:>14,} {r['per_class_median']:>14,} "
              f"{r['r2']:>9.6f} {r['sd_resid']:>10.4f} {r['mean_se']:>10.4f} "
              f"{r['ratio']:>7.3f}    [{time.perf_counter()-t0:.1f}s]")

    print()
    print("Notes:")
    print("  ratio in [0.91, 0.99] = residuals scale with sampling SE noise floor")
    print("  ratio >> 1            = real residual structure beyond prefix prediction")
    print("  a_final levels: theoretical max is k for k >= 1 (one per 3^j, j=1..k);")
    print("                  actual count = number of distinct a_final encountered")

    out = pl.DataFrame(rows)
    out_csv = out_dir / f"24_k_sweep_alpha_decomposition_N{args.N}.csv"
    out.write_csv(out_csv)
    print()
    print(f"[save] {out_csv}")


if __name__ == "__main__":
    main()
