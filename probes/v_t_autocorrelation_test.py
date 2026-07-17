"""
v_t_autocorrelation_test.py

Conditional autocorrelation of v_t given σ-band at q=3, N=2^36.

For each σ-band q ∈ {0.125, 0.375, 0.625, 0.875, 0.975}:
  - Walk orbits, track v_t for t = 0..T_track-1
  - Compute lag-k autocorr ρ_k(band) = Corr(v_t, v_{t+k} | band) for k ∈ {1, 2, 3, 5, 10}
  - 5-seed bootstrap for CIs

Reuses walker logic from experiments/60_per_band_esscher_verify.py.
"""
import sys
import time
from pathlib import Path
import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")
log_lines = []
def log(s): print(s, flush=True); log_lines.append(s)


MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_syracuse_with_v_seq(starts, max_value, max_syr_steps, T_track):
    n = len(starts)
    sigma_arr = np.zeros(n, dtype=np.int64)
    n_odd_arr = np.zeros(n, dtype=np.int64)
    v_seq = np.full((n, T_track), -1, dtype=np.int8)
    ok_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = starts[i]
        sigma_total = 0
        syr_steps = 0
        failed = False
        while m != 1 and syr_steps < max_syr_steps:
            if (m & 1) == 0:
                m = m >> 1
                sigma_total += 1
                continue
            if m > max_value // 3:
                failed = True; break
            threex_p1 = 3 * m + 1
            v = 0
            tmp = threex_p1
            while (tmp & 1) == 0:
                tmp >>= 1
                v += 1
            if syr_steps < T_track:
                v_seq[i, syr_steps] = v
            m = tmp
            sigma_total += 1 + v
            syr_steps += 1
        if not failed and m == 1:
            sigma_arr[i] = sigma_total
            n_odd_arr[i] = syr_steps
            ok_arr[i] = True
    return sigma_arr, n_odd_arr, v_seq, ok_arr


def compute_lag_autocorr(v_seq, n_odd, valid_mask, k):
    """For orbits with valid_mask, compute Corr(v_t, v_{t+k}) pooled across all valid (t, t+k) pairs."""
    pairs_x = []
    pairs_y = []
    valid_idx = np.where(valid_mask)[0]
    T_max = v_seq.shape[1]
    for i in valid_idx:
        end_t = min(int(n_odd[i]), T_max) - k
        if end_t <= 0:
            continue
        x = v_seq[i, :end_t].astype(np.int32)
        y = v_seq[i, k:end_t + k].astype(np.int32)
        if len(x) > 0:
            pairs_x.append(x); pairs_y.append(y)
    if not pairs_x:
        return float("nan"), 0
    X = np.concatenate(pairs_x).astype(np.float64)
    Y = np.concatenate(pairs_y).astype(np.float64)
    if len(X) < 50 or X.std() < 1e-9 or Y.std() < 1e-9:
        return float("nan"), len(X)
    rho = float(np.corrcoef(X, Y)[0, 1])
    return rho, len(X)


def main():
    log("=" * 78)
    log("v_t conditional autocorrelation test at q=3, N=2^36")
    log("=" * 78)

    log2N = 36
    N = 1 << log2N
    n_per_seed = 100_000
    seeds = [42, 137, 271, 314, 1729]
    T_track = 28
    LAGS = [1, 2, 3, 5, 10]
    BANDS = [(0.125, 0.0, 0.25), (0.375, 0.25, 0.50),
             (0.625, 0.50, 0.75), (0.875, 0.75, 1.00),
             (0.975, 0.95, 1.00)]

    # Per-seed, per-band, per-lag autocorrelation values
    results = {(b[0], k): [] for b in BANDS for k in LAGS}
    n_pairs = {(b[0], k): [] for b in BANDS for k in LAGS}

    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1) // 2, size=n_per_seed, dtype=np.int64) + 1
        t0 = time.perf_counter()
        sigma, n_odd, v_seq, ok = walk_syracuse_with_v_seq(starts, MAX_VAL, 1_000_000, T_track)
        log(f"  seed={seed}: walked in {time.perf_counter()-t0:.1f}s, ok={int(ok.sum()):,}")

        starts_ok = starts[ok]
        sigma_ok = sigma[ok].astype(np.float64)
        n_odd_ok = n_odd[ok]
        v_seq_ok = v_seq[ok]

        log_n = np.log(starts_ok.astype(np.float64))
        # OLS sigma ~ a + b * log_n to get residuals
        log_n_c = log_n - log_n.mean()
        sigma_c = sigma_ok - sigma_ok.mean()
        beta = float((log_n_c * sigma_c).sum() / (log_n_c * log_n_c).sum())
        alpha = float(sigma_ok.mean() - beta * log_n.mean())
        sigma_resid = sigma_ok - (alpha + beta * log_n)
        q_thresh = {p: float(np.percentile(sigma_resid, 100 * p))
                    for p in [0.25, 0.50, 0.75, 0.95]}

        for q, lo, hi in BANDS:
            if hi == 1.00 and lo == 0.95:
                mask = sigma_resid > q_thresh[0.95]
            elif lo == 0.75:
                mask = sigma_resid > q_thresh[0.75]
            elif lo == 0.50:
                mask = (sigma_resid > q_thresh[0.50]) & (sigma_resid <= q_thresh[0.75])
            elif lo == 0.25:
                mask = (sigma_resid > q_thresh[0.25]) & (sigma_resid <= q_thresh[0.50])
            else:
                mask = sigma_resid <= q_thresh[0.25]
            for k in LAGS:
                rho, n_p = compute_lag_autocorr(v_seq_ok, n_odd_ok, mask, k)
                results[(q, k)].append(rho)
                n_pairs[(q, k)].append(n_p)

    # Aggregate across seeds
    log(f"\n=== Conditional autocorrelation ρ_k(band) — 5-seed bootstrap ===")
    log(f"  {'band':>7} {'lag':>4} {'mean ρ':>9} {'SD':>8} {'n_pairs':>10} {'verdict':>15}")
    rows_out = []
    for q, _, _ in BANDS:
        for k in LAGS:
            rhos = np.array(results[(q, k)])
            n_p_avg = np.mean(n_pairs[(q, k)])
            mean_rho = float(np.nanmean(rhos))
            sd_rho = float(np.nanstd(rhos))
            verdict = "INDEP (|ρ|<0.02)" if abs(mean_rho) < 0.02 else (
                      "WEAK (0.02-0.05)" if abs(mean_rho) < 0.05 else "STRUCTURE")
            log(f"  {q:>7.3f} {k:>4d} {mean_rho:>+9.4f} {sd_rho:>8.4f} {n_p_avg:>10.0f} {verdict:>15}")
            rows_out.append({"band": q, "lag": k, "mean_rho": mean_rho,
                             "sd_rho": sd_rho, "n_pairs_mean": float(n_p_avg)})

    # Cross-band comparison
    log(f"\n=== Cross-band comparison ===")
    for k in LAGS:
        rhos_across = [float(np.nanmean(results[(q, k)])) for q, _, _ in BANDS]
        log(f"  lag={k}: ρ across bands = "
            f"[{', '.join(f'{r:+.4f}' for r in rhos_across)}], "
            f"max|ρ|={max(abs(r) for r in rhos_across):.4f}")

    # Verdict per brief outcomes
    all_rhos = [results[(q, k)] for q, _, _ in BANDS for k in LAGS]
    max_abs_mean_rho = max(abs(float(np.nanmean(r))) for r in all_rhos)
    log(f"\n=== Decisive outcome ===")
    log(f"  max |mean ρ_k(band)| across all (band, lag) = {max_abs_mean_rho:.4f}")
    if max_abs_mean_rho < 0.02:
        log(f"  --> OUTCOME (a): conditional independence verified within ±0.02.")
        log(f"      Three-slice characterization is sufficient.")
    elif max_abs_mean_rho < 0.05:
        log(f"  --> WEAK structure (max |ρ| < 0.05). Likely outcome (b) — small Markov-like decay.")
    else:
        log(f"  --> SUBSTANTIAL structure (max |ρ| > 0.05). Outcome (b)/(c)/(d). Characterize decay.")
        # Check if decay across lags is consistent across bands
        log(f"\n  Decay structure across lags per band:")
        for q, _, _ in BANDS:
            decay = [float(np.nanmean(results[(q, k)])) for k in LAGS]
            log(f"    band {q}: ρ = {[f'{r:+.4f}' for r in decay]}")

    pl.DataFrame(rows_out).write_csv(Path(r"C:\Collatz") / "v_t_autocorrelation_test.csv")
    Path(r"C:\Collatz/v_t_autocorrelation_test_log.txt").write_text("\n".join(log_lines), encoding="utf-8")
    log(f"\n[wrote] v_t_autocorrelation_test.csv, v_t_autocorrelation_test_log.txt")


if __name__ == "__main__":
    main()
