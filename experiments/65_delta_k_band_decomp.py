"""
ΔK_band algebraic decomposition test.

Within σ-resid band q, K_eff_band(q) = slope(σ on log_n | band).
By identity σ = T·(1+V_orbit), this slope decomposes as:

  slope(σ, log_n | band) = (1+E[V|band])·slope(T, log_n | band)
                         + E[T|band]·slope(V_orbit, log_n | band)
                         + κ_TVZ / Var(log_n | band)

where κ_TVZ = E[(T-E[T])(V-E[V])(log_n-E[log_n]) | band] is the joint third
central moment.

K_bulk(E_band) = (1+E_band)/(E_band·log2 - log3) is the asymptotic random-walk slope.
ΔK_band = K_eff_band − K_bulk = slope_decomp − K_bulk.

Test: does the algebraic decomposition + asymptotic K_bulk reproduce empirical
ΔK_band, OR does it require new structure beyond the σ-identity?
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
K_H = 3.0 / np.log(4.0/3.0)  # ~ 10.4185
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_v_orbit_stats(starts, max_value, max_syr):
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    T_arr = np.zeros(n, dtype=np.int32)
    sumv_arr = np.zeros(n, dtype=np.int64)
    ok_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = np.int64(starts[i])
        sigma_total = 0; T = 0; sumv = 0
        failed = False
        while m != 1 and T < max_syr:
            if (m & 1) == 0:
                m = m >> 1; sigma_total += 1; continue
            if m > max_value // 3:
                failed = True; break
            x = 3*m + 1; v = 0
            while (x & 1) == 0:
                x >>= 1; v += 1
            m = x
            sigma_total += 1 + v
            sumv += v
            T += 1
        if not failed and m == 1:
            sigma_arr[i] = sigma_total
            T_arr[i] = T
            sumv_arr[i] = sumv
            ok_arr[i] = True
    return sigma_arr, T_arr, sumv_arr, ok_arr


def K_bulk(E_band):
    """Asymptotic random-walk slope from per-step Esscher mean."""
    return (1.0 + E_band) / (E_band * LOG_2 - LOG_3)


def slope_OLS(X, Y):
    """OLS slope of Y on X."""
    Xc = X - X.mean(); Yc = Y - Y.mean()
    return float((Xc*Yc).sum() / (Xc*Xc).sum())


def analyze_N(N, n_per_seed=100_000, seeds=(42,137,271,314,1729)):
    log2N = int(np.log2(N))
    print(f"\n# N = 2^{log2N} = {N:,}", flush=True)
    t0 = time.time()

    sigma_all = []; T_all = []; sumv_all = []; logn_all = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1)//2, size=n_per_seed, dtype=np.int64) + 1
        sigma, T, sumv, ok = walk_v_orbit_stats(starts, MAX_VAL, 1_000_000)
        sigma_all.append(sigma[ok].astype(np.float64))
        T_all.append(T[ok].astype(np.float64))
        sumv_all.append(sumv[ok].astype(np.float64))
        logn_all.append(np.log(starts[ok].astype(np.float64)))

    sigma = np.concatenate(sigma_all); T = np.concatenate(T_all)
    sumv = np.concatenate(sumv_all); log_n = np.concatenate(logn_all)
    V_orbit = sumv / np.maximum(T, 1)

    print(f"  total ok: {len(sigma):,}  walk-time: {time.time()-t0:.1f}s", flush=True)

    # Standardize σ via theoretical K_h drift → bands defined in σ_resid quantiles
    sigma_resid = sigma - K_H * log_n

    band_defs = [
        (0.125, 0.0,  0.25),
        (0.375, 0.25, 0.50),
        (0.625, 0.50, 0.75),
        (0.875, 0.75, 1.00),
        (0.975, 0.95, 1.00),
    ]

    rows = []
    for q, lo, hi in band_defs:
        if lo == 0.0:
            lo_val = -np.inf
        else:
            lo_val = float(np.percentile(sigma_resid, lo*100))
        if hi == 1.0:
            hi_val = np.inf
        else:
            hi_val = float(np.percentile(sigma_resid, hi*100))

        mask = (sigma_resid > lo_val) & (sigma_resid <= hi_val)
        n_band = int(mask.sum())
        if n_band < 100: continue

        sigma_b = sigma[mask]; T_b = T[mask]; V_b = V_orbit[mask]; logn_b = log_n[mask]

        E_V_band = float(V_b.mean())
        E_T_band = float(T_b.mean())
        var_logn_band = float(logn_b.var())
        E_logn_band = float(logn_b.mean())

        # K_eff_band via OLS slope(σ on log_n | band)
        K_eff_band_ols = slope_OLS(logn_b, sigma_b)

        # Decomposition: slope_T = slope(T on log_n | band), slope_V = slope(V on log_n | band)
        slope_T = slope_OLS(logn_b, T_b)
        slope_V = slope_OLS(logn_b, V_b)

        # Candidate decomposition (linearization, ignoring cubic term)
        candidate_no_cubic = (1.0 + E_V_band) * slope_T + E_T_band * slope_V

        # Cubic third joint cumulant
        Tc = T_b - E_T_band; Vc = V_b - E_V_band; logn_c = logn_b - E_logn_band
        kappa_TVZ = float((Tc * Vc * logn_c).mean())
        candidate_with_cubic = candidate_no_cubic + kappa_TVZ / var_logn_band

        # K_bulk from Esscher mean
        K_bulk_E = K_bulk(E_V_band)

        # ΔK_band
        Delta_K_emp = K_eff_band_ols - K_bulk_E
        Delta_K_decomp = candidate_with_cubic - K_bulk_E

        # Asymptotic dT/dlogn = 1/(E·log2 - log3)
        dT_dlogn_asymp = 1.0 / (E_V_band * LOG_2 - LOG_3)

        rows.append({
            'log2N': log2N, 'q': q, 'n_band': n_band,
            'E_V_band': E_V_band, 'E_T_band': E_T_band,
            'var_logn_band': var_logn_band,
            'K_eff_ols': K_eff_band_ols,
            'K_bulk_E': K_bulk_E,
            'slope_T': slope_T, 'slope_V': slope_V,
            'dT_dlogn_asymp': dT_dlogn_asymp,
            'kappa_TVZ': kappa_TVZ,
            'cubic_contribution': kappa_TVZ / var_logn_band,
            'candidate_no_cubic': candidate_no_cubic,
            'candidate_with_cubic': candidate_with_cubic,
            'Delta_K_emp': Delta_K_emp,
            'Delta_K_decomp': Delta_K_decomp,
            'gap_decomp': K_eff_band_ols - candidate_with_cubic,  # should be ~0
        })

    return rows


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(exist_ok=True)

    print(f"K_h = {K_H:.6f}", flush=True)
    print(f"K_bulk(E) = (1+E)/(E·log2 − log3) — asymptotic random-walk slope", flush=True)

    all_rows = []
    for log2N in [32, 34, 36, 38]:
        N = 1 << log2N
        rows = analyze_N(N)
        all_rows.extend(rows)

    df = pl.DataFrame(all_rows)

    # Sanity: decomposition should match empirical K_eff_ols (all algebra)
    print(f"\n\n=== Sanity check: decomposition algebraic identity ===", flush=True)
    print(f"  K_eff_ols vs (1+E_V)·slope_T + E_T·slope_V + κ_TVZ/Var(log_n)", flush=True)
    print(f"  {'log2N':>6}  {'q':>6}  {'K_eff_ols':>10}  {'cand_w/cubic':>13}  {'gap':>10}", flush=True)
    for r in all_rows:
        print(f"  {r['log2N']:>6}  {r['q']:>6.3f}  {r['K_eff_ols']:>10.4f}  "
              f"{r['candidate_with_cubic']:>13.4f}  {r['gap_decomp']:>+10.5f}", flush=True)

    # Per-N display: empirical ΔK_band, K_bulk, components
    print(f"\n\n=== K_eff_ols vs K_bulk; ΔK_band per band ===", flush=True)
    print(f"  {'log2N':>6}  {'q':>6}  {'E_V':>7}  {'K_bulk':>8}  {'K_eff_ols':>10}  "
          f"{'ΔK_band':>9}  {'slope_T':>9}  {'slope_V':>9}  {'cubic':>10}", flush=True)
    for r in all_rows:
        print(f"  {r['log2N']:>6}  {r['q']:>6.3f}  {r['E_V_band']:>7.4f}  "
              f"{r['K_bulk_E']:>8.3f}  {r['K_eff_ols']:>10.3f}  "
              f"{r['Delta_K_emp']:>+9.3f}  {r['slope_T']:>9.4f}  {r['slope_V']:>+9.5f}  "
              f"{r['cubic_contribution']:>+9.5f}", flush=True)

    # Decomposition of ΔK_band into pieces:
    # ΔK_band = (1+E_V)·[slope_T - dT/dlogn_asymp] + E_T·slope_V + cubic
    print(f"\n\n=== ΔK_band decomposition pieces ===", flush=True)
    print(f"  ΔK_band = (1+E_V)·(slope_T − 1/(E·log2−log3))  +  E_T·slope_V  +  κ_TVZ/Var(log_n)", flush=True)
    print(f"  {'log2N':>6}  {'q':>6}  {'piece_T':>10}  {'piece_V':>10}  {'piece_cubic':>11}  "
          f"{'sum':>9}  {'ΔK_emp':>9}  {'gap':>9}", flush=True)
    for r in all_rows:
        piece_T = (1.0 + r['E_V_band']) * (r['slope_T'] - r['dT_dlogn_asymp'])
        piece_V = r['E_T_band'] * r['slope_V']
        piece_cubic = r['cubic_contribution']
        sum_pieces = piece_T + piece_V + piece_cubic
        print(f"  {r['log2N']:>6}  {r['q']:>6.3f}  {piece_T:>+9.4f}  {piece_V:>+9.4f}  "
              f"{piece_cubic:>+10.5f}  {sum_pieces:>+8.4f}  {r['Delta_K_emp']:>+8.4f}  "
              f"{sum_pieces - r['Delta_K_emp']:>+8.5f}", flush=True)

    # Verdict on whether decomposition closes
    print(f"\n\n=== VERDICT ===", flush=True)
    max_gap = float(np.abs([r['gap_decomp'] for r in all_rows]).max())
    print(f"  max |K_eff_ols - candidate_with_cubic|: {max_gap:.5f}", flush=True)
    if max_gap < 0.01:
        print(f"  → algebraic identity exact; decomposition is tautological from σ = T·(1+V)", flush=True)

    # The substantive question: are slope_T, slope_V, κ_TVZ derivable from Esscher framework?
    # Test: is slope_T close to dT/dlogn_asymp = 1/(E_V·log2 - log3)?
    print(f"\n=== Substantive test: slope_T vs Esscher-asymptotic 1/(E·log2 − log3) ===", flush=True)
    print(f"  {'log2N':>6}  {'q':>6}  {'slope_T':>9}  {'asymptotic':>11}  {'ratio':>8}  {'gap':>9}", flush=True)
    for r in all_rows:
        ratio = r['slope_T'] / r['dT_dlogn_asymp']
        gap = r['slope_T'] - r['dT_dlogn_asymp']
        print(f"  {r['log2N']:>6}  {r['q']:>6.3f}  {r['slope_T']:>9.5f}  "
              f"{r['dT_dlogn_asymp']:>11.5f}  {ratio:>8.4f}  {gap:>+9.5f}", flush=True)

    print(f"\n=== slope_V vs zero (asymptotic, V=E_band constant) ===", flush=True)
    for r in all_rows:
        print(f"  log2N={r['log2N']}  q={r['q']:.3f}  slope_V = {r['slope_V']:+.6f}  "
              f"E_T·slope_V = {r['E_T_band']*r['slope_V']:+.4f}", flush=True)

    out_csv = out_dir / "65_delta_k_band_decomp.csv"
    df.write_csv(out_csv)
    print(f"\n[save] {out_csv}", flush=True)


if __name__ == "__main__":
    main()
