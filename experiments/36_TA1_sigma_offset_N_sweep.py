"""
Experiment 36 (TA.1) — N-stability of sigma structural offset.

Tests whether the sigma gap of approximately -2.4 nats is constant in N or
drifts. Build sigma cache at N in {2^25, 2^27, 2^28, 2^30, 2^32}, compute
per-class mean sigma at k=8, compare to alpha_det + K_h * log_n_class_mean.

Reports:
  - slope of sigma_mean ~ alpha_det + K_h * log_n at each N
  - offset (= the gap)
  - per-class SE
  - whether variation across N exceeds variation across k at fixed N

Memory-conscious: at N=2^32 the sigma cache is ~17 GB.

Usage:
    python 36_TA1_sigma_offset_N_sweep.py
"""
import gc
import sys
import time
from pathlib import Path

import numpy as np
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

K_H = 3.0 / (np.log(4.0) - np.log(3.0))
PATH_MAX = 200_000


@njit(cache=True, boundscheck=False)
def compute_sigma_cache(N):
    sigma = np.full(N + 1, -1, dtype=np.int32)
    sigma[1] = 0
    path = np.empty(PATH_MAX, dtype=np.int64)
    for n in range(2, N + 1):
        if sigma[n] != -1:
            continue
        m = np.int64(n)
        path_len = 0
        overflow = False
        while True:
            if m <= N and sigma[m] != -1:
                break
            if path_len >= PATH_MAX:
                overflow = True
                break
            path[path_len] = m
            path_len += 1
            if m & 1:
                m = 3 * m + 1
            else:
                m = m >> 1
            if m == 1:
                break
        if overflow:
            m = np.int64(n)
            steps = 0
            while m > 1:
                if m & 1:
                    m = 3 * m + 1
                else:
                    m = m >> 1
                steps += 1
            sigma[n] = np.int32(steps)
            continue
        cur = np.int32(0) if m == 1 else sigma[m]
        for i in range(path_len - 1, -1, -1):
            v = path[i]
            cur += 1
            if v <= N:
                sigma[v] = cur
    return sigma


@njit(cache=True, parallel=True, boundscheck=False)
def per_class_sums(sigma, N, k):
    """Per-class sum of sigma, sum of log(n), and count, on odd n in [3, N]."""
    M = 1 << k
    K_cls = M // 2
    n_threads = 64
    chunk = max(2, N // n_threads)
    if chunk % 2 == 0:
        chunk += 1
    n_chunks = (N + chunk - 1) // chunk
    sigma_sum = np.zeros((n_chunks, K_cls), dtype=np.float64)
    sigma_sum_sq = np.zeros((n_chunks, K_cls), dtype=np.float64)
    log_sum = np.zeros((n_chunks, K_cls), dtype=np.float64)
    cnt = np.zeros((n_chunks, K_cls), dtype=np.int64)
    mask_M = M - 1
    for c in prange(n_chunks):
        lo = c * chunk
        hi = min(lo + chunk, N + 1)
        if lo < 3:
            lo = 3
        if lo % 2 == 0:
            lo += 1
        for n in range(lo, hi, 2):
            r = n & mask_M
            kk = (r - 1) >> 1
            s = float(sigma[n])
            sigma_sum[c, kk] += s
            sigma_sum_sq[c, kk] += s * s
            log_sum[c, kk] += np.log(np.float64(n))
            cnt[c, kk] += 1
    return sigma_sum.sum(axis=0), sigma_sum_sq.sum(axis=0), log_sum.sum(axis=0), cnt.sum(axis=0)


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


def compute_alpha_det(k):
    M = 1 << k
    K_cls = M // 2
    alpha = np.zeros(K_cls)
    for kk in range(K_cls):
        r = 2 * kk + 1
        s, a_f, _ = deterministic_prefix(r, M)
        alpha[kk] = s + K_H * np.log(a_f / float(M))
    return alpha


def fit_slope_offset(x, y, weights=None):
    if weights is None:
        weights = np.ones_like(x)
    w_sum = weights.sum()
    x_bar = (weights * x).sum() / w_sum
    y_bar = (weights * y).sum() / w_sum
    x_c = x - x_bar
    y_c = y - y_bar
    slope = (weights * x_c * y_c).sum() / (weights * x_c * x_c).sum()
    offset = y_bar - slope * x_bar
    pred = slope * x + offset
    ss_res = (weights * (y - pred) ** 2).sum()
    ss_tot = (weights * y_c ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return float(slope), float(offset), float(r2)


def process_N(N, ks=(8, 10, 12)):
    M_max = 1 << max(ks)
    print()
    print("#" * 100)
    print(f"# N = {N:,} = 2^{int(round(np.log2(N)))}")
    print("#" * 100)
    print(f"  sigma cache memory: {(N + 1) * 4 / 1e9:.2f} GB")
    t0 = time.perf_counter()
    sigma = compute_sigma_cache(N)
    print(f"  cache built in {time.perf_counter()-t0:.1f}s   max(sigma)={int(sigma[1:].max())}", flush=True)

    # Sanity check
    known = {1: 0, 27: 111, 703: 170, 871: 178, 6171: 261, 77031: 350, 837799: 524}
    for nk, sk in known.items():
        if nk <= N:
            assert int(sigma[nk]) == sk, f"sigma({nk}) != {sk}"
    print(f"  sigma checkpoints OK", flush=True)

    rows = []
    for k in ks:
        t0 = time.perf_counter()
        sig_sum, sig_sq, log_sum, cnt = per_class_sums(sigma, N, k)
        cnt_safe = np.maximum(cnt, 1)
        sigma_mean = sig_sum / cnt_safe
        log_n_mean = log_sum / cnt_safe
        sigma_var = sig_sq / cnt_safe - sigma_mean * sigma_mean
        # Per-class SE on the mean
        sigma_se = np.sqrt(np.maximum(sigma_var, 0) / cnt_safe)

        alpha_det = compute_alpha_det(k)
        # y = sigma_mean - K_h * log_n_mean (the "Tao-residual")
        y = sigma_mean - K_H * log_n_mean
        x = alpha_det
        slope, offset, r2 = fit_slope_offset(x, y)

        # Also report y bar and x bar at <log N>
        log_n_grand = log_n_mean.mean()
        gap_at_mean = sigma_mean.mean() - alpha_det.mean() - K_H * log_n_grand

        # Median per-class SE
        median_SE = float(np.median(sigma_se))

        print(f"  k={k:>2}: slope={slope:.4f} offset={offset:>+.4f} R^2={r2:.6f} "
              f"<log N>={log_n_grand:.4f} gap_at_mean={gap_at_mean:>+.4f} "
              f"median_class_SE={median_SE:.4f}    [{time.perf_counter()-t0:.1f}s]", flush=True)
        rows.append((N, k, slope, offset, r2, gap_at_mean, log_n_grand, median_SE))

    # Free
    del sigma
    gc.collect()
    return rows


def main():
    Ns = [1 << 25, 1 << 27, 1 << 28, 1 << 30, 1 << 32]
    all_rows = []
    for N in Ns:
        all_rows.extend(process_N(N, ks=(8, 10, 12)))

    # Cross-N comparison
    print()
    print("=" * 100)
    print("Summary: σ offset (gap from K_h * <log N>) vs N at each k")
    print("=" * 100)
    print(f"{'N':>14} {'log2(N)':>9} {'<log N>':>10} {'k':>3} "
          f"{'slope':>9} {'offset':>10} {'gap_at_mean':>13} {'med_class_SE':>13}")
    for r in all_rows:
        N, k, slope, offset, r2, gap, log_n, se = r
        print(f"{N:>14,} {int(round(np.log2(N))):>9} {log_n:>10.4f} {k:>3} "
              f"{slope:>9.4f} {offset:>+10.4f} {gap:>+13.4f} {se:>13.4f}")

    # Variation across N at fixed k
    print()
    print("Variation in gap_at_mean across N, at each k:")
    for k in (8, 10, 12):
        gaps = [r[5] for r in all_rows if r[1] == k]
        if len(gaps) > 1:
            print(f"  k={k}: gaps = {[f'{g:>+.4f}' for g in gaps]}   "
                  f"range = {max(gaps) - min(gaps):.4f}   "
                  f"std = {np.std(gaps):.4f}")

    # Variation across k at fixed N
    print()
    print("Variation in gap_at_mean across k, at each N:")
    for N in Ns:
        gaps = [r[5] for r in all_rows if r[0] == N]
        if len(gaps) > 1:
            print(f"  N=2^{int(round(np.log2(N)))}: gaps = {[f'{g:>+.4f}' for g in gaps]}   "
                  f"range = {max(gaps) - min(gaps):.4f}   "
                  f"std = {np.std(gaps):.4f}")

    # Save
    import polars as pl
    df = pl.DataFrame(
        all_rows,
        schema=["N", "k", "slope", "offset", "R2", "gap_at_mean", "log_n_mean", "median_class_SE"],
        orient="row",
    )
    out = Path(__file__).resolve().parent.parent / "experiments_output"
    out.mkdir(parents=True, exist_ok=True)
    df.write_csv(out / "36_TA1_sigma_offset_N_sweep.csv")
    print()
    print(f"[save] {out / '36_TA1_sigma_offset_N_sweep.csv'}")


if __name__ == "__main__":
    main()
