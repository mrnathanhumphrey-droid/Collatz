"""
Two tests:
  T1: w_q / z_q exactly 0.13 across extended bands q ∈ {0.05, 0.125, 0.20,
       0.375, 0.50, 0.625, 0.80, 0.875, 0.95, 0.975}?
  T2: ΔK_band = K_eff_band − K_bulk per band; does it correlate with the
       entry-class structure (P(j=2)·W_2 weighted by tilt)?

Walks 500K orbits at N=2^36, tracks per-orbit (σ, n_odd, v_seq, last_v before
terminate, first-passage steps to 4 thresholds). Per band:
  - E[v]_band → w_q, K_bulk = K(E[v])
  - K_eff_band via first-passage regression
  - ΔK_band = K_eff − K_bulk
  - P(j=2 | band) = P(orbit's last v == 4) — last Syracuse step lands at 1 from m=5
    (since 3·5+1 = 16 = 2^4, so v_last = 4 ⇔ m_last = 5 = m_j for j=2)
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
K_H = 3.0 / (LOG_2*2 - LOG_3)
MAX_VAL = np.int64(2**62)


def E_v_from_w(w):
    return 1.0 / (1.0 - 2.0**(-(1.0 + w)))

def w_from_E_v(E_v):
    if E_v <= 1: return float('inf')
    r = 1.0 - 1.0/E_v
    return -np.log2(r) - 1.0

def K_of_v(v):
    return (1 + v) / (v * LOG_2 - LOG_3)


@njit(parallel=True, cache=True)
def walk_orbit_full(starts, thresholds_per_n, max_value, max_syr_steps, T_track):
    """Walks Syracuse, records:
      - sigma, n_odd, v_seq[t<T_track], last_v_before_terminate
      - first-passage step counts to 4 thresholds (DESCENDING from start)
      - ok flag
    """
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    n_odd_arr = np.zeros(n, dtype=np.int32)
    v_seq = np.full((n, T_track), -1, dtype=np.int8)
    last_v_arr = np.full(n, -1, dtype=np.int8)
    s_arr = np.full((n, 4), -1, dtype=np.int32)
    ok_arr = np.zeros(n, dtype=np.bool_)

    for i in prange(n):
        m = np.int64(starts[i])
        sigma_total = 0
        syr_steps = 0
        last_v = np.int8(-1)
        next_idx = 0
        # Initial threshold check (orbit starts above all thresholds usually)
        while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
            s_arr[i, next_idx] = 0
            next_idx += 1
        failed = False
        while m != 1 and syr_steps < max_syr_steps:
            if (m & 1) == 0:
                m = m >> 1
                sigma_total += 1
                # Threshold checks
                while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
                    s_arr[i, next_idx] = sigma_total
                    next_idx += 1
                continue
            if m > max_value // 3:
                failed = True; break
            threex_p1 = 3 * m + 1
            v = np.int8(0)
            tmp = threex_p1
            while (tmp & 1) == 0:
                tmp >>= 1
                v += 1
            last_v = v
            if syr_steps < T_track:
                v_seq[i, syr_steps] = v
            m = tmp
            sigma_total += 1 + v
            syr_steps += 1
            while next_idx < 4 and m <= thresholds_per_n[i, next_idx]:
                s_arr[i, next_idx] = sigma_total
                next_idx += 1
        if not failed and m == 1:
            sigma_arr[i] = sigma_total
            n_odd_arr[i] = syr_steps
            last_v_arr[i] = last_v
            ok_arr[i] = True
    return sigma_arr, n_odd_arr, v_seq, last_v_arr, s_arr, ok_arr


def K_eff_slope(sigma, s_per_thresh, thresh_phys):
    log_f = np.log(np.maximum(thresh_phys, 1.0))
    R = sigma[:, None] - s_per_thresh.astype(np.float64)
    log_f_mean = log_f.mean(axis=0)
    R_mean = R.mean(axis=0)
    x = log_f_mean - log_f_mean.mean()
    y = R_mean - R_mean.mean()
    return float((x * y).sum() / (x * x).sum())


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"

    log2N = 36
    N = 1 << log2N
    n_per_seed = 100_000
    seeds = [42, 137, 271, 314, 1729]
    T_track = 28

    # Extended quantile bands: midpoints
    q_list = [0.05, 0.125, 0.20, 0.375, 0.50, 0.625, 0.80, 0.875, 0.95, 0.975]
    # Each band defined by quantile interval; we'll select from σ_resid quantile cuts
    # For midpoint q with width Δq centered at q (use Δq=0.10 unless at boundary):
    # Or use disjoint quartile-style edges:
    # 0.05: [0, 0.10), 0.125: [0.10, 0.15], ... need careful design.
    # Actually simplest: each midpoint q corresponds to interval [q-Δ/2, q+Δ/2) for some Δ.
    # Use Δq = 0.10 except at edges.
    # Equivalently: take percentile cuts {0.0, 0.10, 0.15, 0.25, 0.30, 0.45, 0.55, 0.70, 0.85, 0.95, 1.0}
    #
    # Cleaner: use disjoint bands at percentiles that match the midpoints
    # Bands and edges (low_q, high_q):
    band_defs = [
        (0.05,  0.00, 0.10),  # 0-10% midpoint 0.05
        (0.125, 0.10, 0.15),  # 10-15% midpoint 0.125 (width 5%)
        (0.20,  0.15, 0.25),  # 15-25% midpoint 0.20
        (0.375, 0.25, 0.50),  # 25-50% midpoint 0.375
        (0.50,  0.45, 0.55),  # 45-55% midpoint 0.50 (width 10%)
        (0.625, 0.50, 0.75),  # 50-75% midpoint 0.625
        (0.80,  0.75, 0.85),  # 75-85% midpoint 0.80
        (0.875, 0.85, 0.90),  # 85-90% midpoint 0.875 (note: this overlaps with q4 quartile)
        (0.95,  0.90, 1.00),  # 90-100% midpoint 0.95
        (0.975, 0.95, 1.00),  # 95-100% midpoint 0.975 (sub-band of 0.95)
    ]

    print(f"# Walking 500K orbits at N=2^{log2N}", flush=True)

    # Aggregation per band: (E[v], K_eff_band, P(last_v=4), P(last_v=8), counts)
    K_max = 12
    band_data = {q: {'sum_v': 0.0, 'sum_v2': 0.0, 'count_v': 0,
                     'sigma': [], 's_phys': [], 'raw_thresh': [],
                     'last_v_count': np.zeros(20, dtype=np.int64),
                     'orbit_count': 0,
                     'P_v_t': np.zeros((T_track, K_max+1), dtype=np.int64),
                     'count_t': np.zeros(T_track, dtype=np.int64)}
                 for q, _, _ in band_defs}

    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1) // 2, size=n_per_seed, dtype=np.int64) + 1
        log_n = np.log(starts.astype(np.float64))
        sqrt_n = np.sqrt(starts.astype(np.float64))
        nt = starts.astype(np.float64) ** (2.0 / 3.0)
        sl = sqrt_n * log_n
        sdl = sqrt_n / np.maximum(log_n, 1.0)
        raw = np.column_stack([nt, sl, sqrt_n, sdl])
        sort_idx = np.argsort(-raw, axis=1)
        th = np.take_along_axis(raw, sort_idx, axis=1).astype(np.int64)

        t0 = time.perf_counter()
        sigma, n_odd, v_seq, last_v, s_arr, ok = walk_orbit_full(
            starts, th, MAX_VAL, 1_000_000, T_track)
        elapsed = time.perf_counter() - t0
        print(f"  seed={seed}: walked in {elapsed:.1f}s, ok={int(ok.sum()):,}", flush=True)

        starts_ok = starts[ok]; sigma_ok = sigma[ok].astype(np.float64)
        n_odd_ok = n_odd[ok]; v_seq_ok = v_seq[ok]; last_v_ok = last_v[ok]
        s_arr_ok = s_arr[ok]; sort_idx_ok = sort_idx[ok]; raw_ok = raw[ok]

        inv_sort = np.argsort(sort_idx_ok, axis=1)
        s_phys = np.zeros((len(starts_ok), 4), dtype=np.int32)
        for col in range(4):
            s_phys[:, col] = np.take_along_axis(s_arr_ok, inv_sort[:, col:col+1], axis=1).flatten()

        log_n_ok = np.log(starts_ok.astype(np.float64))
        log_n_c = log_n_ok - log_n_ok.mean()
        sigma_c = sigma_ok - sigma_ok.mean()
        beta = float((log_n_c * sigma_c).sum() / (log_n_c * log_n_c).sum())
        alpha = float(sigma_ok.mean() - beta * log_n_ok.mean())
        sigma_resid = sigma_ok - (alpha + beta * log_n_ok)

        # For each band, get the σ_resid percentile cuts and select orbits
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
            if n_band == 0: continue

            # Aggregate v_per_orbit (sum-of-v / n_odd)
            n_odd_band = n_odd_ok[mask]
            # Per-orbit V = (sum of v across orbit) / n_odd. We can sum v_seq or use n_odd.
            # Better: sum v_seq for each orbit (need to mask -1) — but for orbits with n_odd > T_track,
            # we don't have all v's. Use n_odd-derived: total v in orbit = sigma - n_odd
            # (since sigma = n_odd + sum_v over orbit)
            sum_v_orbit = sigma_ok[mask] - n_odd_band.astype(np.float64)
            V_orbit = sum_v_orbit / np.maximum(n_odd_band.astype(np.float64), 1)
            band_data[q]['sum_v'] += float(V_orbit.sum())
            band_data[q]['sum_v2'] += float((V_orbit**2).sum())
            band_data[q]['count_v'] += n_band
            band_data[q]['orbit_count'] += n_band
            # First-passage: collect for K_eff_band regression
            band_data[q]['sigma'].append(sigma_ok[mask])
            band_data[q]['s_phys'].append(s_phys[mask])
            band_data[q]['raw_thresh'].append(raw_ok[mask])
            # Last v
            for lv in range(20):
                band_data[q]['last_v_count'][lv] += int((last_v_ok[mask] == lv).sum())
            # v_t per t
            for t in range(T_track):
                valid_t = (n_odd_band > t) & np.ones(n_band, dtype=bool)
                v_at_t = v_seq_ok[mask][valid_t, t].astype(np.int32)
                band_data[q]['count_t'][t] += len(v_at_t)
                for k in range(1, K_max+1):
                    band_data[q]['P_v_t'][t, k] += int((v_at_t == k).sum())

    # Compute per-band E[v], K_eff_band, K_bulk, ΔK
    print(f"\n=== Per-band test results ===\n", flush=True)
    print(f"  {'q':>6}  {'z_q':>9}  {'E[v]':>8}  {'w_q':>10}  {'w_q/z_q':>9}  "
          f"{'K_bulk':>8}  {'K_eff_band':>11}  {'ΔK':>10}  {'P(j=2)':>9}", flush=True)

    rows = []
    for q, lo, hi in band_defs:
        d = band_data[q]
        if d['count_v'] == 0: continue
        E_v = d['sum_v'] / d['count_v']
        # K_eff_band via combined first-passage regression
        sigma_b = np.concatenate(d['sigma'])
        s_phys_b = np.concatenate(d['s_phys'])
        raw_b = np.concatenate(d['raw_thresh'])
        K_eff = K_eff_slope(sigma_b, s_phys_b, raw_b)
        K_bulk = K_of_v(E_v)
        delta_K = K_eff - K_bulk
        w_q = w_from_E_v(E_v)
        z_q = float(norm.ppf(q))
        ratio = w_q / z_q if abs(z_q) > 0.01 else float('nan')
        # P(j=2): orbit landed via last_v=4 (since 3·5+1=16=2^4)
        P_j2 = d['last_v_count'][4] / d['orbit_count']

        print(f"  {q:>6.3f}  {z_q:>+8.3f}  {E_v:>8.4f}  {w_q:>+9.4f}  {ratio:>9.4f}  "
              f"{K_bulk:>8.3f}  {K_eff:>11.3f}  {delta_K:>+9.3f}  {P_j2:>9.5f}", flush=True)

        rows.append({
            'q': q, 'z_q': z_q,
            'E_v_band': float(E_v),
            'w_q': float(w_q),
            'w_q_over_z_q': float(ratio),
            'K_bulk': float(K_bulk),
            'K_eff_band': float(K_eff),
            'delta_K': float(delta_K),
            'P_last_v_eq_4': float(P_j2),
        })

    # T1 verdict: w_q / z_q across bands
    print(f"\n\n=== T1: w_q / z_q ratio across bands ===\n", flush=True)
    ratios = [r['w_q_over_z_q'] for r in rows if abs(r['z_q']) > 0.01]
    if ratios:
        ratios = np.array(ratios)
        print(f"  Mean ratio: {ratios.mean():.4f}", flush=True)
        print(f"  SD ratio:   {ratios.std(ddof=1):.4f}", flush=True)
        print(f"  Range:      [{ratios.min():.4f}, {ratios.max():.4f}]", flush=True)
        if ratios.std(ddof=1) < 0.01:
            print(f"  → CONSTANT — w_q ≈ {ratios.mean():.3f} · z_q is closed-form", flush=True)
        elif ratios.std(ddof=1) < 0.03:
            print(f"  → Near-constant; small structural correction may be present", flush=True)
        else:
            print(f"  → NOT a constant; w_q(q) is not simply linear in z_q", flush=True)

    # T2 verdict: ΔK_band correlation with P(j=2)·W_2 weighted
    # W_2 ≈ 7.156 from compute_threads_findings.md
    W_2 = 7.156
    print(f"\n=== T2: ΔK_band correlation with P(j=2 | band) · W_2 = {W_2} ===\n", flush=True)
    print(f"  {'q':>6}  {'ΔK_band':>10}  {'P(j=2|band)':>12}  {'P·W_2':>10}  {'ratio ΔK/(P·W_2)':>20}", flush=True)
    for r in rows:
        pred = r['P_last_v_eq_4'] * W_2
        ratio_dk = r['delta_K'] / pred if abs(pred) > 1e-3 else float('nan')
        r['P_j2_W2'] = pred
        r['delta_K_over_pred'] = ratio_dk
        print(f"  {r['q']:>6.3f}  {r['delta_K']:>+9.3f}  {r['P_last_v_eq_4']:>12.5f}  "
              f"{pred:>+9.3f}  {ratio_dk:>20.4f}", flush=True)

    # Correlation
    dks = np.array([r['delta_K'] for r in rows])
    preds = np.array([r['P_j2_W2'] for r in rows])
    corr = float(np.corrcoef(dks, preds)[0, 1])
    print(f"\n  Pearson correlation (ΔK_band vs P(j=2)·W_2): {corr:+.4f}", flush=True)
    # OLS slope
    pc = preds - preds.mean(); dc = dks - dks.mean()
    slope = float((pc*dc).sum() / (pc*pc).sum())
    intercept = float(dks.mean() - slope*preds.mean())
    pred_dk = intercept + slope*preds
    R2 = 1 - float(np.sum((dks - pred_dk)**2) / np.sum((dks - dks.mean())**2))
    print(f"  OLS fit: ΔK_band = {intercept:+.3f} + {slope:+.3f} · (P·W_2)  R² = {R2:.4f}", flush=True)

    pl.DataFrame(rows).write_csv(out_dir / "61_w_q_linearity_and_boundary.csv")
    print(f"\n[save] CSV", flush=True)


if __name__ == "__main__":
    main()
