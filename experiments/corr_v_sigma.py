"""
Approach 1 — derive E[v]_q125 = 2.216 via joint-Gaussian conditional moment.

For each orbit i: V_i = n_even_i / n_odd_i (per-orbit mean v_t).
                  Σ_i = total stopping time σ_i.

If (V_i, Σ_i) is bivariate Gaussian under the trajectory measure:
  E[V | Σ < q_25] = E[V] + ρ · SD[V] · E[Z_Σ | Z_Σ < z_25]
  E[Z|Z < -0.6745] = -φ(0.6745)/Φ(-0.6745) = -0.3178/0.25 = -1.2729

Computes per-orbit V and Σ from existing parquet (N=2^27 = 134M, full
enumeration), measures Corr(V, Σ), plugs into truncated-Gaussian formula.
Compares prediction to empirical E[V|q_125] at this N AND extrapolates via
the predicted-vs-empirical relationship.

Also walks 500K orbits at higher N (2^36, 2^40, 2^42) to see if Corr(V, Σ)
and the joint-Gaussian prediction stabilize toward 2.216.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from scipy.stats import norm

sys.stdout.reconfigure(encoding="utf-8")

K_H = 3.0 / (np.log(4.0) - np.log(3.0))
LOG_2 = np.log(2.0)
LOG_3 = np.log(3.0)
MAX_VAL = np.int64(2**62)

# Truncated-Gaussian: E[Z | Z < z_25] for z_25 = norm.ppf(0.25)
Z_25 = float(norm.ppf(0.25))                       # -0.6745
PHI_AT_Z25 = float(norm.pdf(Z_25))                 # 0.3178
TRUNC_E_Z = -PHI_AT_Z25 / 0.25                     # -1.2729


def K_of_v(v):
    return (1 + v) / (v * LOG_2 - LOG_3)


def analyze_orbits(n_arr, sigma_arr, n_odd, n_even, label):
    """Per-orbit V = n_even/n_odd; compute trunc-Gaussian prediction for E[V|q_125]."""
    log_n = np.log(n_arr.astype(np.float64))
    V_per_orbit = n_even.astype(np.float64) / np.maximum(n_odd.astype(np.float64), 1)
    sigma_f = sigma_arr.astype(np.float64)

    # Detrend σ by log_n (we care about RESIDUAL σ, not absolute, since absolute
    # σ varies hugely with log n in a heterogeneous-N parquet).
    # Within a single-N walker (orbits in [3, N]), log_n still varies a lot, so detrend.
    # Use linear regression σ = α + β·log_n; work with residuals σ_resid.
    log_n_c = log_n - log_n.mean()
    sigma_c = sigma_f - sigma_f.mean()
    beta = float((log_n_c * sigma_c).sum() / (log_n_c * log_n_c).sum())
    alpha = float(sigma_f.mean() - beta * log_n.mean())
    sigma_resid = sigma_f - (alpha + beta * log_n)

    # Per-orbit moments
    E_V = float(V_per_orbit.mean())
    SD_V = float(V_per_orbit.std(ddof=1))
    E_sigma = float(sigma_f.mean())
    SD_sigma = float(sigma_f.std(ddof=1))
    SD_sigma_resid = float(sigma_resid.std(ddof=1))

    # Two correlations: raw Corr(V, σ) and Corr(V, σ_resid) [σ_resid removes log_n trend]
    rho_raw = float(np.corrcoef(V_per_orbit, sigma_f)[0, 1])
    rho_resid = float(np.corrcoef(V_per_orbit, sigma_resid)[0, 1])

    # Truncated-Gaussian prediction (per-orbit conditioning on σ_resid quartile)
    # The bottom σ-quartile band selects orbits with smallest σ AT SAME log_n,
    # i.e., bottom quartile of σ_resid. Use rho_resid.
    pred_E_V_q125 = E_V + rho_resid * SD_V * TRUNC_E_Z

    # Empirical: bottom-σ-quartile conditional E[V]
    # Within this single-N parquet, bottom σ-quartile = bottom quartile of σ_f among orbits
    # in this parquet. But because σ scales with log_n, this is a non-uniform mix.
    # Strict version: bottom quartile of σ_resid (σ at same log_n).
    q25_resid = np.percentile(sigma_resid, 25)
    mask_low_resid = sigma_resid <= q25_resid
    E_V_q125_resid = float(V_per_orbit[mask_low_resid].mean())

    # Loose version: bottom quartile of σ regardless of log_n (matches the K_eff_band
    # measurement in prior exps which used σ directly, not σ_resid).
    q25_raw = np.percentile(sigma_f, 25)
    mask_low_raw = sigma_f <= q25_raw
    E_V_q125_raw = float(V_per_orbit[mask_low_raw].mean())

    print(f"\n=== {label} ===")
    print(f"  N orbits = {len(n_arr):,}")
    print(f"  Per-orbit moments:")
    print(f"    E[V]         = {E_V:.4f}")
    print(f"    SD[V]        = {SD_V:.4f}")
    print(f"    E[σ]         = {E_sigma:.2f}")
    print(f"    SD[σ]        = {SD_sigma:.2f}")
    print(f"    SD[σ_resid]  = {SD_sigma_resid:.2f}")
    print(f"  Correlations:")
    print(f"    Corr(V, σ)         = {rho_raw:+.4f}")
    print(f"    Corr(V, σ_resid)   = {rho_resid:+.4f}")
    print(f"  Truncated-Gaussian prediction E[V | σ_resid < q_25]:")
    print(f"    = {E_V:.4f} + {rho_resid:+.4f} · {SD_V:.4f} · {TRUNC_E_Z:.4f}")
    print(f"    = {pred_E_V_q125:.4f}")
    print(f"  Empirical:")
    print(f"    E[V | σ_resid < q_25]   = {E_V_q125_resid:.4f}")
    print(f"    E[V | σ_raw < q_25]     = {E_V_q125_raw:.4f}")
    print(f"  Asymptote target: 2.2161")
    print(f"  Gap to target (resid version): {E_V_q125_resid - 2.2161:+.4f}")

    return {
        'label': label, 'n_orbits': len(n_arr),
        'E_V': E_V, 'SD_V': SD_V,
        'E_sigma': E_sigma, 'SD_sigma': SD_sigma, 'SD_sigma_resid': SD_sigma_resid,
        'rho_raw': rho_raw, 'rho_resid': rho_resid,
        'pred_E_V_q125': pred_E_V_q125,
        'E_V_q125_resid': E_V_q125_resid,
        'E_V_q125_raw': E_V_q125_raw,
    }


@njit(parallel=True, cache=True)
def walk_orbit_v(starts, max_value, max_steps):
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    n_odd = np.zeros(n, dtype=np.int32)
    n_even = np.zeros(n, dtype=np.int32)
    ok_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = np.int64(starts[i])
        steps = 0
        odd_count = 0; even_count = 0
        failed = False
        while m != 1 and steps < max_steps:
            if m & 1:
                if m > max_value // 3:
                    failed = True; break
                m = 3 * m + 1
                odd_count += 1
            else:
                m = m >> 1
                even_count += 1
            steps += 1
        if not failed and m == 1:
            sigma_arr[i] = steps
            ok_arr[i] = True
            n_odd[i] = odd_count
            n_even[i] = even_count
    return sigma_arr, ok_arr, n_odd, n_even


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "experiments_output"

    results = []

    # 1. From existing parquet at N=2^27 (134M rows, full enumeration)
    N_val = 1 << 27
    print(f"# Loading parquet at N=2^27 ({N_val:,})")
    df = pl.read_parquet(data_dir / f"main_N{N_val}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    sigma_arr = df["sigma"].to_numpy().astype(np.int32)
    odd_steps = df["odd_steps"].to_numpy().astype(np.int32)
    even_steps = df["even_steps"].to_numpy().astype(np.int32)
    print(f"  Loaded {len(n_arr):,} odd orbits")

    r = analyze_orbits(n_arr, sigma_arr, odd_steps, even_steps,
                        f"Full enumeration N=2^27")
    results.append(r)

    # Subsample the parquet to match 500K-orbit walker conventions
    rng = np.random.default_rng(42)
    idx = rng.choice(len(n_arr), size=500_000, replace=False)
    r2 = analyze_orbits(n_arr[idx], sigma_arr[idx], odd_steps[idx], even_steps[idx],
                         f"Subsample 500K from N=2^27")
    results.append(r2)

    # 2. Walk fresh orbits at higher N (sampling)
    for log2N in [30, 32, 34, 36, 38, 40, 42]:
        N = 1 << log2N
        print(f"\n# Walking {500_000:,} orbits at N=2^{log2N}")
        rng = np.random.default_rng(42)
        starts = 2 * rng.integers(1, (N - 1) // 2, size=500_000, dtype=np.int64) + 1
        t0 = time.perf_counter()
        sigma, ok, n_odd, n_even = walk_orbit_v(starts, MAX_VAL, 500_000)
        print(f"  Walk: {time.perf_counter()-t0:.1f}s, ok={ok.sum():,}, failed={(~ok).sum()}")
        r = analyze_orbits(starts[ok], sigma[ok], n_odd[ok], n_even[ok],
                            f"Sample N=2^{log2N}")
        results.append(r)

    # Trajectory across N
    print(f"\n=== Joint-Gaussian prediction vs empirical, across N ===")
    print(f"  {'label':<30}  {'rho_resid':>10}  {'SD[V]':>9}  {'pred':>9}  {'emp(resid)':>11}  {'emp(raw)':>10}")
    for r in results:
        print(f"  {r['label']:<30}  {r['rho_resid']:>+10.4f}  {r['SD_V']:>9.4f}  "
              f"{r['pred_E_V_q125']:>9.4f}  {r['E_V_q125_resid']:>11.4f}  "
              f"{r['E_V_q125_raw']:>10.4f}")

    # Save
    rows = [{'label': r['label'], 'n_orbits': r['n_orbits'],
             'E_V': r['E_V'], 'SD_V': r['SD_V'],
             'rho_raw': r['rho_raw'], 'rho_resid': r['rho_resid'],
             'pred_E_V_q125': r['pred_E_V_q125'],
             'E_V_q125_resid': r['E_V_q125_resid'],
             'E_V_q125_raw': r['E_V_q125_raw']} for r in results]
    out_csv = out_dir / "corr_v_sigma.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}")


if __name__ == "__main__":
    main()
