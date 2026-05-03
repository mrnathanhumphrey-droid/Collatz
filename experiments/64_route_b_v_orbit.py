"""
Route B: Direct E[V_orbit | σ-band] derivation from Esscher-tilted joint structure.

Per Result 25 (verified), per-step v_t conditional on σ-band is Esscher-tilted
Geom(1/2) with band-specific w_q. Stationary across t. So E[v_t | band] = E_band
(known per-band).

Question: Does E[V_orbit | band] = E_band, where V_orbit = (1/T)·Σ_t v_t?

If yes (within ±0.005): V_orbit | band closes via Esscher directly. Outcome (a).
If no: characterize correction. Test functional forms.

Also: identity σ_orbit = T·(1 + V_orbit) constrains the (T, V_orbit) joint
within a band. Compute Cov[V_orbit, T | band], E[V_orbit | T, band] vs T.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from scipy.stats import norm

import io
sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

LOG_2 = np.log(2.0); LOG_3 = np.log(3.0)
K_H = 3.0 / np.log(4.0/3.0)  # ~ 10.4185
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_v_orbit_stats(starts, max_value, max_syr):
    """Walk orbits; return σ, T (=n_odd), and per-orbit v stats."""
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    T_arr = np.zeros(n, dtype=np.int32)
    sumv_arr = np.zeros(n, dtype=np.int64)  # Σ_t v_t
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


@njit(parallel=True, cache=True)
def walk_with_v_seq(starts, max_value, max_syr, T_track):
    """Walk orbits and track v_seq for first T_track Syracuse steps (for E[v_t|T,band] check)."""
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    T_arr = np.zeros(n, dtype=np.int32)
    v_seq = np.full((n, T_track), -1, dtype=np.int8)
    ok_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = np.int64(starts[i])
        sigma_total = 0; T = 0
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
            if T < T_track:
                v_seq[i, T] = v
            T += 1
        if not failed and m == 1:
            sigma_arr[i] = sigma_total
            T_arr[i] = T
            ok_arr[i] = True
    return sigma_arr, T_arr, v_seq, ok_arr


def trunc_z_moments(a, b):
    Phi_a = norm.cdf(a); Phi_b = norm.cdf(b)
    phi_a = norm.pdf(a); phi_b = norm.pdf(b)
    P = Phi_b - Phi_a
    if P < 1e-15: return 0.0, 0.0
    EZ = (phi_a - phi_b) / P
    a_phi_a = a * phi_a if np.isfinite(a) else 0.0
    b_phi_b = b * phi_b if np.isfinite(b) else 0.0
    EZ2 = 1.0 + (a_phi_a - b_phi_b) / P
    return EZ, EZ2 - 1.0


def analyze_N(N, n_per_seed=100_000, seeds=(42,137,271,314,1729), do_v_seq=False):
    log2N = int(np.log2(N))
    print(f"\n# N = 2^{log2N} = {N:,}", flush=True)
    t0 = time.time()

    if do_v_seq:
        T_track = 32
        v_seq_all = []

    sigma_all = []; T_all = []; sumv_all = []; logn_all = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1)//2, size=n_per_seed, dtype=np.int64) + 1
        if do_v_seq:
            sigma, T, v_seq, ok = walk_with_v_seq(starts, MAX_VAL, 1_000_000, T_track)
            v_seq_all.append(v_seq[ok])
        else:
            sigma, T, sumv, ok = walk_v_orbit_stats(starts, MAX_VAL, 1_000_000)
            sumv_all.append(sumv[ok].astype(np.float64))
        sigma_all.append(sigma[ok].astype(np.float64))
        T_all.append(T[ok].astype(np.float64))
        logn_all.append(np.log(starts[ok].astype(np.float64)))

    sigma = np.concatenate(sigma_all); T = np.concatenate(T_all); log_n = np.concatenate(logn_all)
    if do_v_seq:
        v_seq = np.concatenate(v_seq_all)
        sumv = (v_seq * (v_seq >= 0)).sum(axis=1).astype(np.float64)
        # NOTE: for v_seq mode we use only first T_track steps' sum, not full sum
        # but V_orbit = sum_v / T_orbit only valid if T <= T_track
        # We'll handle both: full V_orbit needs full sumv, so use non-v_seq mode for V_orbit metric
    else:
        sumv = np.concatenate(sumv_all)
    print(f"  total ok: {len(sigma):,}  walk-time: {time.time()-t0:.1f}s", flush=True)

    V_orbit = sumv / np.maximum(T, 1)
    # Sanity: σ = T·(1 + V_orbit)?
    sigma_check = T * (1 + V_orbit)
    diff = float(np.abs(sigma - sigma_check).max())
    print(f"  σ - T·(1+V_orbit) max abs diff: {diff:.2e}  (expect 0 — exact identity)", flush=True)

    # Standardize σ via theoretical K_h drift
    sigma_resid = sigma - K_H * log_n
    mu_S = float(sigma_resid.mean()); sd_S = float(sigma_resid.std())

    # Build σ-quartile bands per Result 25 definitions
    band_defs = [
        (0.125, 0.0,  0.25),
        (0.375, 0.25, 0.50),
        (0.625, 0.50, 0.75),
        (0.875, 0.75, 1.00),
        (0.975, 0.95, 1.00),
    ]

    # μ_V global
    mu_V_global = float(V_orbit.mean()); sd_V_global = float(V_orbit.std())
    print(f"  global μ_V_orbit = {mu_V_global:.4f}  σ_V_orbit = {sd_V_global:.4f}", flush=True)

    band_results = []
    for q, lo, hi in band_defs:
        if lo == 0.0:
            lo_val = -np.inf; a_z = -np.inf
        else:
            lo_val = float(np.percentile(sigma_resid, lo*100))
            a_z = (lo_val - mu_S) / sd_S
        if hi == 1.0:
            hi_val = np.inf; b_z = np.inf
        else:
            hi_val = float(np.percentile(sigma_resid, hi*100))
            b_z = (hi_val - mu_S) / sd_S

        mask = (sigma_resid > lo_val) & (sigma_resid <= hi_val)
        n_band = int(mask.sum())
        if n_band < 100: continue

        # E[V_orbit | band]
        V_in = V_orbit[mask]; T_in = T[mask]; sumv_in = sumv[mask]
        E_V_orbit_band = float(V_in.mean())
        SE_V_orbit_band = float(V_in.std() / np.sqrt(n_band))

        # E_band (per-step): Σ all v over Σ all T = total sumv / total T in band
        # This is the per-step average: average v across all (orbit, t) pairs in band
        E_band_per_step = float(sumv_in.sum() / T_in.sum())
        SE_band_per_step = float(np.sqrt(((sumv_in/np.maximum(T_in,1) - E_band_per_step)**2 * T_in / T_in.mean()).mean() / n_band))

        # Conditional moments
        cov_VT_band = float(((V_in - V_in.mean()) * (T_in - T_in.mean())).mean())
        var_T_band = float(T_in.var())
        var_V_band = float(V_in.var())
        rho_VT_band = cov_VT_band / np.sqrt(var_V_band * var_T_band) if var_V_band > 1e-15 and var_T_band > 1e-15 else 0.0

        # E_band - E[V_orbit] correction
        correction = E_V_orbit_band - E_band_per_step

        # Truncated Gaussian moments for Edgeworth-shape comparison
        EZ_band, EZ2m1_band = trunc_z_moments(a_z, b_z)

        band_results.append({
            'q': q, 'lo': lo, 'hi': hi, 'a_z': a_z, 'b_z': b_z,
            'n_band': n_band,
            'E_V_orbit_band': E_V_orbit_band,
            'SE_V_orbit_band': SE_V_orbit_band,
            'E_band_per_step': E_band_per_step,
            'correction': correction,
            'cov_VT': cov_VT_band, 'rho_VT': rho_VT_band,
            'mean_T': float(T_in.mean()), 'sd_T': float(T_in.std()),
            'mean_V': float(V_in.mean()), 'sd_V': float(V_in.std()),
            'EZ_band': EZ_band, 'EZ2m1_band': EZ2m1_band,
        })

    return {
        'log2N': log2N, 'mu_V_global': mu_V_global, 'sd_V_global': sd_V_global,
        'mu_S': mu_S, 'sd_S': sd_S,
        'bands': band_results,
        'V_orbit': V_orbit, 'sigma_resid': sigma_resid, 'T': T, 'sumv': sumv,
    }


def stratify_v_orbit_by_T(r):
    """For the q=0.125 band, stratify orbits by T-deciles and compute E[V_orbit | T-decile]."""
    log2N = r['log2N']
    sigma_resid = r['sigma_resid']; V_orbit = r['V_orbit']; T = r['T']
    mu_S = r['mu_S']; sd_S = r['sd_S']

    print(f"\n  --- T-stratification within q=0.125 band at N=2^{log2N} ---", flush=True)
    # bottom 25% σ-resid band
    cut = float(np.percentile(sigma_resid, 25))
    mask = sigma_resid <= cut
    V_in = V_orbit[mask]; T_in = T[mask]
    print(f"  band size: {mask.sum():,}  T range: [{T_in.min():.0f}, {T_in.max():.0f}]  mean T: {T_in.mean():.1f}", flush=True)

    # 10 T-deciles
    T_qs = np.percentile(T_in, np.linspace(10, 90, 9))
    T_edges = np.concatenate([[T_in.min()-0.5], T_qs, [T_in.max()+0.5]])
    print(f"  decile  T_lo   T_hi   n   E[V_orbit|T-decile]   E_band-per-step|T-decile", flush=True)
    for i in range(len(T_edges)-1):
        m = (T_in > T_edges[i]) & (T_in <= T_edges[i+1])
        n = int(m.sum())
        if n < 50: continue
        E_V = float(V_in[m].mean())
        # per-step E[v_t | T-decile, band]: should equal E_V if v_t stationary within (T,band) cell
        # For each orbit, V_orbit = sumv / T. So mean over orbits of V_orbit weighted by T = E[v_t | band, T-decile]
        # Without v_seq we can compute E[v_t | (T-decile, band)] = sum(sumv[m]) / sum(T[m])
        sumv_in = r['sumv'][mask]
        Esep = float(sumv_in[m].sum() / T_in[m].sum())
        print(f"     {i+1:>2}  {T_edges[i]:>5.0f}  {T_edges[i+1]:>5.0f}  {n:>6,}  "
              f"{E_V:>20.4f}  {Esep:>22.4f}", flush=True)


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(exist_ok=True)

    print(f"K_h (theoretical) = 3/log(4/3) = {K_H:.6f}", flush=True)
    print(f"Identity check: σ_orbit = T·(1 + V_orbit) where V_orbit = Σv_t / T", flush=True)

    all_results = {}
    for log2N in [32, 34, 36, 38]:
        N = 1 << log2N
        r = analyze_N(N)
        all_results[log2N] = r

    # Per-N per-band comparison
    print(f"\n\n=== E[V_orbit | band] vs E_band-per-step (Esscher prediction) ===", flush=True)
    print(f"  {'log2N':>6}  {'q':>6}  {'n_band':>8}  "
          f"{'E[V_orbit|band]':>16}  {'E_band-per-step':>16}  {'correction':>11}  {'SE':>9}", flush=True)
    for log2N, r in all_results.items():
        for b in r['bands']:
            print(f"  {log2N:>6}  {b['q']:>6.3f}  {b['n_band']:>8,}  "
                  f"{b['E_V_orbit_band']:>16.4f}  {b['E_band_per_step']:>16.4f}  "
                  f"{b['correction']:>+10.4f}  {b['SE_V_orbit_band']:>9.5f}", flush=True)

    # Conditional Cov[V, T | band]
    print(f"\n\n=== Conditional Cov[V_orbit, T | band] ===", flush=True)
    print(f"  {'log2N':>6}  {'q':>6}  {'mean_T':>9}  {'sd_T':>9}  "
          f"{'mean_V':>9}  {'sd_V':>9}  {'ρ(V,T)':>10}", flush=True)
    for log2N, r in all_results.items():
        for b in r['bands']:
            print(f"  {log2N:>6}  {b['q']:>6.3f}  {b['mean_T']:>9.2f}  {b['sd_T']:>9.2f}  "
                  f"{b['mean_V']:>9.4f}  {b['sd_V']:>9.4f}  {b['rho_VT']:>+9.4f}", flush=True)

    # T-stratification within q=0.125 (most extreme band) at N=2^36
    stratify_v_orbit_by_T(all_results[36])

    # Correction analysis: is correction structurally functional in band?
    print(f"\n\n=== Correction (E[V_orbit|band] − E_band-per-step) — fit candidates ===", flush=True)
    for log2N, r in all_results.items():
        bands = r['bands']
        EZ2m1 = np.array([b['EZ2m1_band'] for b in bands])
        EZ = np.array([b['EZ_band'] for b in bands])
        corr = np.array([b['correction'] for b in bands])
        sd_V_g = r['sd_V_global']

        # Through-origin fit corr = c·EZ2m1
        denom = float((EZ2m1**2).sum())
        c_EZ2 = float((EZ2m1 * corr).sum() / denom) if denom > 1e-15 else 0.0
        pred_EZ2 = c_EZ2 * EZ2m1
        ss_res = float(((corr - pred_EZ2)**2).sum())
        ss_tot = float(((corr - corr.mean())**2).sum())
        R2_EZ2 = 1 - ss_res/ss_tot if ss_tot > 1e-15 else 1.0

        # corr = c·EZ
        denom = float((EZ**2).sum())
        c_EZ = float((EZ * corr).sum() / denom) if denom > 1e-15 else 0.0
        pred_EZ = c_EZ * EZ
        ss_res2 = float(((corr - pred_EZ)**2).sum())
        R2_EZ = 1 - ss_res2/ss_tot if ss_tot > 1e-15 else 1.0

        max_abs = float(np.abs(corr).max())
        print(f"  N=2^{log2N}  σ_V_g={sd_V_g:.4f}  max|corr|={max_abs:.4f}", flush=True)
        print(f"    fit corr = c·E[Z²-1|band]: c = {c_EZ2:+.5f}  R² = {R2_EZ2:.4f}  "
              f"c/σ_V = {c_EZ2/sd_V_g:+.4f}", flush=True)
        print(f"    fit corr = c·E[Z|band]:   c = {c_EZ:+.5f}  R² = {R2_EZ:.4f}", flush=True)

    # Compare correction to Result 29 C ≈ 0.21·σ_V "phenomenological"
    print(f"\n=== Compare correction-shape coeff to Result 29 'C' empirical (σ_V-rescaled) ===", flush=True)
    for log2N, r in all_results.items():
        bands = r['bands']
        EZ2m1 = np.array([b['EZ2m1_band'] for b in bands])
        corr = np.array([b['correction'] for b in bands])
        denom = float((EZ2m1**2).sum())
        c_EZ2 = float((EZ2m1 * corr).sum() / denom) if denom > 1e-15 else 0.0
        sd_V_g = r['sd_V_global']
        print(f"  log2N={log2N}  c (Route B) = {c_EZ2:+.5f}  σ_V = {sd_V_g:.4f}  c/σ_V = {c_EZ2/sd_V_g:+.4f}", flush=True)

    # ============ VERDICT ============
    print(f"\n=== VERDICT ===", flush=True)
    max_corr_per_N = {}
    for log2N, r in all_results.items():
        max_abs = float(np.abs([b['correction'] for b in r['bands']]).max())
        max_corr_per_N[log2N] = max_abs
        print(f"  N=2^{log2N}  max |correction| = {max_abs:.4f}", flush=True)
    overall_max = max(max_corr_per_N.values())
    if overall_max < 0.005:
        print("  → outcome (a): E[V_orbit|band] = E_band within ±0.005. Esscher closes V_orbit|band exactly.", flush=True)
    elif overall_max < 0.05:
        print("  → outcome (b): correction exists but is small (<0.05). Functional form testable.", flush=True)
    else:
        print("  → outcome (c) or (d): correction ≥ 0.05. Needs full conditional structure analysis.", flush=True)

    # Save
    rows = []
    for log2N, r in all_results.items():
        for b in r['bands']:
            rows.append({
                'log2N': log2N, 'q': b['q'],
                'n_band': b['n_band'],
                'E_V_orbit_band': b['E_V_orbit_band'],
                'E_band_per_step': b['E_band_per_step'],
                'correction': b['correction'],
                'SE_V_orbit_band': b['SE_V_orbit_band'],
                'mean_T': b['mean_T'], 'sd_T': b['sd_T'],
                'mean_V': b['mean_V'], 'sd_V': b['sd_V'],
                'rho_VT_band': b['rho_VT'],
                'EZ_band': b['EZ_band'], 'EZ2m1_band': b['EZ2m1_band'],
                'sd_V_global': r['sd_V_global'],
            })
    df = pl.DataFrame(rows)
    out_csv = out_dir / "64_route_b_v_orbit.csv"
    df.write_csv(out_csv)
    print(f"\n[save] {out_csv}", flush=True)


if __name__ == "__main__":
    main()
