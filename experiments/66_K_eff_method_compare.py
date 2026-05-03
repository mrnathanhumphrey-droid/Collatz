"""
Compare K_eff_OLS (within-band σ on log_n regression) to K_eff_first_passage
(per-band slope of mean R = σ - σ_threshold on mean log threshold) on the
SAME orbits, with matched bootstrap.

If they differ systematically: ΔK_band U-shape is in methodology difference,
not a single structural quantity. Constant 4 "boundary correction" depends on
which definition you use.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

import io
sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

LOG_2 = np.log(2.0); LOG_3 = np.log(3.0)
K_H = 3.0 / np.log(4.0/3.0)
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_with_thresholds(starts, thresh_arr, max_value, max_steps):
    """Walk on UNFOLDED Collatz (3n+1, m/2 atomic steps), capture σ at 4 thresholds."""
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    s_at_thresh = np.full((n, 4), -1, dtype=np.int32)
    T_arr = np.zeros(n, dtype=np.int32)  # # odd steps for σ-identity decomposition
    sumv_arr = np.zeros(n, dtype=np.int64)
    ok_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = np.int64(starts[i])
        steps = 0; T = 0; sumv = 0; next_idx = 0
        # Initial threshold passes
        while next_idx < 4 and m <= thresh_arr[i, next_idx]:
            s_at_thresh[i, next_idx] = 0; next_idx += 1
        failed = False
        while m != 1 and steps < max_steps:
            if m & 1:
                if m > max_value // 3:
                    failed = True; break
                # 3m+1 followed by halvings: count v
                x = 3*m + 1; v = 0
                while (x & 1) == 0:
                    x >>= 1; v += 1
                m = x
                steps += 1 + v  # σ unfolded: 1 odd step + v halvings
                T += 1
                sumv += v
            else:
                m = m >> 1
                steps += 1
            while next_idx < 4 and m <= thresh_arr[i, next_idx]:
                s_at_thresh[i, next_idx] = steps; next_idx += 1
        if not failed and m == 1:
            sigma_arr[i] = steps
            T_arr[i] = T
            sumv_arr[i] = sumv
            ok_arr[i] = True
    return sigma_arr, s_at_thresh, T_arr, sumv_arr, ok_arr


def K_eff_first_passage(sigma_band, s_band, log_thresh_band):
    """4-point regression: mean R vs mean log threshold, restricted to band orbits."""
    R = sigma_band[:, None] - s_band.astype(np.float64)
    x = log_thresh_band.mean(axis=0)
    y = R.mean(axis=0)
    xc = x - x.mean(); yc = y - y.mean()
    return float((xc*yc).sum() / (xc*xc).sum())


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"

    log2N = 36
    N = 1 << log2N
    n_per_seed = 100_000
    seeds = [42, 137, 271, 314, 1729]

    print(f"# N = 2^{log2N}, {len(seeds)} seeds × {n_per_seed:,} = {len(seeds)*n_per_seed:,} orbits", flush=True)

    band_defs = [
        (0.125, 0.0, 0.25),
        (0.375, 0.25, 0.50),
        (0.625, 0.50, 0.75),
        (0.875, 0.75, 1.00),
        (0.975, 0.95, 1.00),
    ]

    K_ols_per_seed = {q: [] for q, _, _ in band_defs}
    K_fp_per_seed = {q: [] for q, _, _ in band_defs}

    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1)//2, size=n_per_seed, dtype=np.int64) + 1
        starts = starts.astype(np.int64)
        # 4 thresholds per orbit: N^(2/3), √n·log n, √n, √n/log n
        n_f = starts.astype(np.float64)
        thresh = np.column_stack([
            np.power(n_f, 2.0/3.0),
            np.sqrt(n_f) * np.log(n_f),
            np.sqrt(n_f),
            np.sqrt(n_f) / np.log(n_f),
        ]).astype(np.int64)
        # Sort each row descending so threshold[0] is highest (passed first as orbit descends)
        thresh = -np.sort(-thresh, axis=1)
        log_thresh = np.log(np.maximum(thresh.astype(np.float64), 1.0))

        sigma, s_at_thresh, T, sumv, ok = walk_with_thresholds(starts, thresh, MAX_VAL, 1_000_000)
        ok_mask = ok & (s_at_thresh > 0).all(axis=1) if False else ok  # simpler: just ok
        sigma = sigma[ok_mask]; s_at_thresh = s_at_thresh[ok_mask]
        T = T[ok_mask]; sumv = sumv[ok_mask]
        log_thresh = log_thresh[ok_mask]
        starts_ok = starts[ok_mask]
        log_n = np.log(starts_ok.astype(np.float64))

        sigma = sigma.astype(np.float64); T = T.astype(np.float64); sumv = sumv.astype(np.float64)
        sigma_resid = sigma - K_H * log_n

        for q, lo, hi in band_defs:
            if lo == 0.0:
                lo_v = -np.inf
            else:
                lo_v = float(np.percentile(sigma_resid, lo*100))
            if hi == 1.0:
                hi_v = np.inf
            else:
                hi_v = float(np.percentile(sigma_resid, hi*100))
            mask = (sigma_resid > lo_v) & (sigma_resid <= hi_v)
            if mask.sum() < 100: continue

            sig_b = sigma[mask]; logn_b = log_n[mask]
            # K_ols: within-band σ-on-log_n
            xc = logn_b - logn_b.mean(); yc = sig_b - sig_b.mean()
            K_ols = float((xc*yc).sum() / (xc*xc).sum())
            K_ols_per_seed[q].append(K_ols)

            # K_first_passage on same band orbits
            valid = (s_at_thresh[mask] > 0).all(axis=1)
            if valid.sum() < 100: continue
            K_fp = K_eff_first_passage(sig_b[valid], s_at_thresh[mask][valid], log_thresh[mask][valid])
            K_fp_per_seed[q].append(K_fp)

    # Bootstrap statistics
    print(f"\n=== K_eff comparison: within-band OLS vs first-passage (5-seed bootstrap, N=2^{log2N}) ===", flush=True)
    print(f"  {'q':>6}  {'K_ols mean':>10}  {'K_ols sd':>10}  {'K_fp mean':>10}  {'K_fp sd':>10}  "
          f"{'K_ols-K_fp':>11}  {'sd(diff)':>10}", flush=True)
    rows = []
    for q, _, _ in band_defs:
        ols = np.array(K_ols_per_seed[q]); fp = np.array(K_fp_per_seed[q])
        if len(ols) < 2 or len(fp) < 2: continue
        diff = ols - fp
        print(f"  {q:>6.3f}  {ols.mean():>10.4f}  {ols.std():>10.4f}  "
              f"{fp.mean():>10.4f}  {fp.std():>10.4f}  "
              f"{diff.mean():>+10.4f}  {diff.std():>10.4f}", flush=True)
        rows.append({
            'q': q,
            'K_ols_mean': float(ols.mean()), 'K_ols_sd': float(ols.std()),
            'K_fp_mean': float(fp.mean()), 'K_fp_sd': float(fp.std()),
            'gap_mean': float(diff.mean()), 'gap_sd': float(diff.std()),
        })

    pl.DataFrame(rows).write_csv(out_dir / "66_K_eff_method_compare.csv")
    print(f"\n[save] CSV", flush=True)


if __name__ == "__main__":
    main()
