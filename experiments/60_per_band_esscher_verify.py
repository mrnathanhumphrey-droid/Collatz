"""
Verify partial-Esscher-tilt framework across all σ-quartile bands at N=2^36.

For each band q ∈ {0.125, 0.375, 0.625, 0.875, 0.975}:
  - Track v_t at each step t=0..27 (cap to avoid survivor bias)
  - Compute E[v_t | band] across t — verify stationarity
  - Solve E_w[v] = E[v]_band for w_q (Esscher tilt parameter, Path C convention)
  - Compare predicted P_w(v=k) to empirical P(v_t=k | band)
  - Predict K_band from E[v]_band via K(v) = (1+v)/(v·log2 - log3)
  - Compare to empirical K_eff_band asymptote

Esscher tilt of Geom(1/2) on {1, 2, ...} at parameter w (Path C):
  P_w(v=k) = (1-r) · r^(k-1)  where r = 2^(-(1-w))
  Equivalently: r = 2^(w-1)
  E_w[v] = 1/(1-r) = 1/(1 - 2^(w-1))
  Inverse: r = 1 - 1/E[v], w = 1 + log_2(r)

  Sanity: w=0 → r=1/2 → E[v]=2.0 (Geom(1/2)) ✓
          w=1 → r=1 → E[v]=∞ (degenerate); actually w<1 required
          For E[v]=4/3 (Path C full tilt): r = 1 - 3/4 = 1/4, w = 1 + log_2(1/4) = -1
          But user said Path C w*=1. Convention difference.

  Use convention E[v]_target → solve r = 1 - 1/E[v], w_path_c = 1 - 2·r ...
  Actually I'll just report (E_w[v], r, w_implied) and check consistency.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

import io
sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

LOG_2 = np.log(2.0)
LOG_3 = np.log(3.0)
K_H = 3.0 / (LOG_2 * 2 - LOG_3)
MAX_VAL = np.int64(2**62)

# Empirical K_eff_band asymptotes from prior work
EMPIR_K_BAND = {
    0.125: 7.50,
    0.375: 9.05,
    0.625: 10.66,
    0.875: 14.50,  # peak then plateau, "asymptote" estimate
    0.975: 18.88,
}


def E_v_from_w_pathc(w):
    """E[v] under Esscher tilt at Path C w (where w=0 unconditional, w=1 full tilt giving 4/3).
    P_w(v=k) ∝ 2^(-k) · 2^(-w·k) · ... actually let me use the simpler convention.
    Tilt on log m drop: P_w(v=k) ∝ Geom(1/2)(v) · e^{-w · (v·log2 - log3)}
                              = 2^(-k) · 2^(-w·k) · 2^(w·log_2(3))
                              ∝ 2^(-k(1+w))
    Z_w = sum_{k≥1} 2^(-k(1+w)) = 2^(-(1+w)) / (1 - 2^(-(1+w)))
    P_w(v=k) = 2^(-k(1+w)) · (1 - 2^(-(1+w))) / 2^(-(1+w))
            = (1 - 2^(-(1+w))) · 2^(-(k-1)(1+w))
    E_w[v] = 1 / (1 - 2^(-(1+w)))

    Sanity:
      w=0: E[v] = 1/(1 - 1/2) = 2 ✓
      w=1: E[v] = 1/(1 - 1/4) = 4/3 ✓
    """
    return 1.0 / (1.0 - 2.0**(-(1.0 + w)))


def w_from_E_v_pathc(E_v):
    """Inverse: solve E[v] = 1/(1 - 2^(-(1+w))) for w."""
    if E_v <= 1: return float('inf')
    r = 1.0 - 1.0/E_v   # = 2^(-(1+w))
    if r <= 0: return float('inf')
    one_plus_w = -np.log2(r)
    return one_plus_w - 1.0


def P_v_geom_tilt(k, w):
    """P_w(v=k) for Esscher-tilted Geom(1/2) at Path C parameter w."""
    one_plus_w = 1.0 + w
    p = 1.0 - 2.0**(-one_plus_w)
    return p * (1.0 - p)**(k - 1) if k >= 1 else 0.0


def K_of_v(v):
    return (1 + v) / (v * LOG_2 - LOG_3)


@njit(parallel=True, cache=True)
def walk_syracuse_with_v_seq(starts, max_value, max_syr_steps, T_track):
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    n_odd_arr = np.zeros(n, dtype=np.int32)
    v_seq = np.full((n, T_track), -1, dtype=np.int8)
    ok_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = np.int64(starts[i])
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


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"

    log2N = 36
    N = 1 << log2N
    n_per_seed = 100_000
    seeds = [42, 137, 271, 314, 1729]
    T_track = 28  # cap below survivor-bias range

    # Sanity check on the formula
    print(f"# Sanity check on E_w[v] formula:", flush=True)
    for w_test in [0.0, 1.0, -0.135]:
        ev = E_v_from_w_pathc(w_test)
        w_back = w_from_E_v_pathc(ev)
        print(f"  w={w_test:+.4f} → E[v]={ev:.4f} → w_back={w_back:+.4f}", flush=True)
    print(f"# Sanity: q=0.125 expected w_q ≈ -0.136, E[v]=2.216:", flush=True)
    print(f"  E_w_pathc(-0.136) = {E_v_from_w_pathc(-0.136):.4f}  (expect 2.216)", flush=True)
    print(f"  w_from_E_v(2.216) = {w_from_E_v_pathc(2.216):+.4f}  (expect -0.136)", flush=True)
    print(f"  K(2.216) = {K_of_v(2.216):.4f}  (expect 7.35)", flush=True)
    print(flush=True)

    print(f"# Walking 500K orbits at N=2^{log2N}, tracking T={T_track} steps", flush=True)

    # Storage: per-band counts and v-distribution
    # Bands: q1=0-25, q2=25-50, q3=50-75, q4=75-100, q5=95-100
    # Midpoints: 0.125, 0.375, 0.625, 0.875, 0.975
    band_labels = [(0.125, 'q1: 0-25%'),
                   (0.375, 'q2: 25-50%'),
                   (0.625, 'q3: 50-75%'),
                   (0.875, 'q4: 75-100%'),
                   (0.975, 'q5: 95-100%')]

    K_max = 12  # track P(v=k) for k=1..K_max
    P_vt_per_band = {q: np.zeros((T_track, K_max+1), dtype=np.int64) for q, _ in band_labels}
    counts_per_band = {q: np.zeros(T_track, dtype=np.int64) for q, _ in band_labels}
    sum_v_per_band = {q: np.zeros(T_track, dtype=np.float64) for q, _ in band_labels}
    sumsq_v_per_band = {q: np.zeros(T_track, dtype=np.float64) for q, _ in band_labels}

    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1) // 2, size=n_per_seed, dtype=np.int64) + 1
        t0 = time.perf_counter()
        sigma, n_odd, v_seq, ok = walk_syracuse_with_v_seq(starts, MAX_VAL, 1_000_000, T_track)
        elapsed = time.perf_counter() - t0
        print(f"  seed={seed}: walked in {elapsed:.1f}s, ok={int(ok.sum()):,}", flush=True)

        starts_ok = starts[ok]
        sigma_ok = sigma[ok].astype(np.float64)
        n_odd_ok = n_odd[ok]
        v_seq_ok = v_seq[ok]

        log_n = np.log(starts_ok.astype(np.float64))
        log_n_c = log_n - log_n.mean()
        sigma_c = sigma_ok - sigma_ok.mean()
        beta = float((log_n_c * sigma_c).sum() / (log_n_c * log_n_c).sum())
        alpha = float(sigma_ok.mean() - beta * log_n.mean())
        sigma_resid = sigma_ok - (alpha + beta * log_n)
        q25 = float(np.percentile(sigma_resid, 25))
        q50 = float(np.percentile(sigma_resid, 50))
        q75 = float(np.percentile(sigma_resid, 75))
        q95 = float(np.percentile(sigma_resid, 95))

        # Per-band masks
        masks = {
            0.125: sigma_resid <= q25,
            0.375: (sigma_resid > q25) & (sigma_resid <= q50),
            0.625: (sigma_resid > q50) & (sigma_resid <= q75),
            0.875: sigma_resid > q75,
            0.975: sigma_resid > q95,
        }

        for q, mask in masks.items():
            for t in range(T_track):
                valid = (n_odd_ok > t) & mask
                v_at_t = v_seq_ok[valid, t].astype(np.int32)
                counts_per_band[q][t] += len(v_at_t)
                sum_v_per_band[q][t] += float(v_at_t.sum())
                sumsq_v_per_band[q][t] += float((v_at_t.astype(np.int64) ** 2).sum())
                for k in range(1, K_max + 1):
                    P_vt_per_band[q][t, k] += int((v_at_t == k).sum())

    # Compute E[v_t | band] and P(v=k | band) per t
    print(f"\n=== Per-band E[v_t] and stationarity check (t = 0..{T_track-1}) ===", flush=True)

    band_results = {}
    for q, label in band_labels:
        E_vt = sum_v_per_band[q] / np.maximum(counts_per_band[q], 1)
        Var_vt = sumsq_v_per_band[q] / np.maximum(counts_per_band[q], 1) - E_vt**2
        # P(v=k | band, t)
        P_vt = P_vt_per_band[q] / np.maximum(counts_per_band[q][:, None], 1)
        # Average across t (stationarity check)
        E_v_band = float(E_vt.mean())
        SD_E_vt = float(E_vt.std(ddof=1))
        # Average P(v=k) across t
        P_band = P_vt.mean(axis=0)

        band_results[q] = {
            'label': label,
            'E_vt_per_t': E_vt,
            'E_v_band': E_v_band,
            'SD_E_vt': SD_E_vt,
            'P_band': P_band,
            'P_vt': P_vt,
            'counts': counts_per_band[q],
        }

        print(f"\n  {label}  (mid q = {q})", flush=True)
        print(f"    E[v_t] across t (first 28 + survivor-bias check):", flush=True)
        print(f"      {'t':>3}", end='', flush=True)
        for t in [0, 1, 5, 10, 15, 20, 25, 27]:
            print(f"  {'t='+str(t):>7}", end='', flush=True)
        print('', flush=True)
        print(f"      {'E':>3}", end='', flush=True)
        for t in [0, 1, 5, 10, 15, 20, 25, 27]:
            print(f"  {E_vt[t]:>7.4f}", end='', flush=True)
        print('', flush=True)
        print(f"    Mean E[v_t] (stationary range): {E_v_band:.4f}, SD across t: {SD_E_vt:.4f}", flush=True)

    # Solve for w_q per band
    print(f"\n\n=== Esscher-tilt parameters w_q per band (Path C convention) ===", flush=True)
    print(f"  {'band':>20}  {'E[v]_band':>10}  {'w_q':>9}  {'E_w[v] check':>13}  {'gap':>9}", flush=True)
    for q, label in band_labels:
        E_v_band = band_results[q]['E_v_band']
        w_q = w_from_E_v_pathc(E_v_band)
        E_v_back = E_v_from_w_pathc(w_q)
        gap = E_v_band - E_v_back
        band_results[q]['w_q'] = w_q
        band_results[q]['E_v_back'] = E_v_back
        print(f"  {label:>20}  {E_v_band:>10.4f}  {w_q:>+9.4f}  "
              f"{E_v_back:>13.4f}  {gap:>+9.5f}", flush=True)

    # Verify P_w(v=k) prediction matches empirical
    print(f"\n\n=== P(v=k | band) — empirical vs Esscher-tilt prediction ===", flush=True)
    for q, label in band_labels:
        E_v_band = band_results[q]['E_v_band']
        w_q = band_results[q]['w_q']
        P_band_emp = band_results[q]['P_band']
        print(f"\n  {label}:  E[v]={E_v_band:.4f}, w_q={w_q:+.4f}", flush=True)
        print(f"    {'k':>3}  {'P_w(v=k) pred':>14}  {'P emp':>9}  {'gap':>9}", flush=True)
        max_gap = 0.0
        for k in range(1, 7):
            pred = P_v_geom_tilt(k, w_q)
            emp = P_band_emp[k]
            gap = emp - pred
            max_gap = max(max_gap, abs(gap))
            print(f"    {k:>3}  {pred:>14.5f}  {emp:>9.5f}  {gap:>+9.5f}", flush=True)
        band_results[q]['max_pred_gap'] = max_gap
        verdict = ("MATCH (≤0.005)" if max_gap < 0.005 else
                   "approx (≤0.02)" if max_gap < 0.02 else
                   "FAIL (>0.02)")
        print(f"    Max |gap|: {max_gap:.5f} → {verdict}", flush=True)

    # Compute K_band predicted vs empirical
    print(f"\n\n=== K_band: predicted from E[v]_band via K(v) formula vs empirical asymptote ===", flush=True)
    print(f"  {'band':>20}  {'E[v]_band':>10}  {'K(E[v]) pred':>14}  {'K emp':>10}  {'gap':>10}", flush=True)
    for q, label in band_labels:
        E_v_band = band_results[q]['E_v_band']
        K_pred = K_of_v(E_v_band)
        K_emp = EMPIR_K_BAND.get(q, float('nan'))
        gap = K_pred - K_emp if not np.isnan(K_emp) else float('nan')
        band_results[q]['K_pred'] = K_pred
        band_results[q]['K_emp'] = K_emp
        print(f"  {label:>20}  {E_v_band:>10.4f}  {K_pred:>14.4f}  "
              f"{K_emp:>10.4f}  {gap:>+10.4f}", flush=True)

    # Aggregate K_full prediction (4 quartiles, equal weight)
    K_band_q1234 = [band_results[q]['K_pred'] for q in [0.125, 0.375, 0.625, 0.875]]
    K_full_pred = float(np.mean(K_band_q1234))
    print(f"\n  K_full predicted = mean of quartile K_band's = {K_full_pred:.4f}", flush=True)
    print(f"  K_h reference                                  = {K_H:.4f}", flush=True)
    print(f"  Gap to K_h: {K_full_pred - K_H:+.4f}", flush=True)

    # w_q vs q smooth function?
    print(f"\n\n=== w_q vs q trajectory ===", flush=True)
    print(f"  {'q':>6}  {'w_q':>9}  {'E[v]_band':>10}", flush=True)
    for q, _ in band_labels:
        print(f"  {q:>6.3f}  {band_results[q]['w_q']:>+9.4f}  {band_results[q]['E_v_band']:>10.4f}", flush=True)

    # Save CSV
    rows = []
    for q, label in band_labels:
        r = band_results[q]
        for k in range(1, 7):
            rows.append({
                'q_band': q, 'label': label,
                'E_v_band': r['E_v_band'],
                'w_q': r['w_q'],
                'k': k,
                'P_pred': float(P_v_geom_tilt(k, r['w_q'])),
                'P_emp': float(r['P_band'][k]),
                'gap': float(r['P_band'][k] - P_v_geom_tilt(k, r['w_q'])),
                'K_pred': r['K_pred'],
                'K_emp': r['K_emp'],
                'K_gap': r['K_pred'] - r['K_emp'] if not np.isnan(r['K_emp']) else float('nan'),
                'SD_E_vt': r['SD_E_vt'],
            })
    pl.DataFrame(rows).write_csv(out_dir / "60_per_band_esscher_verify.csv")
    print(f"\n[save] {out_dir / '60_per_band_esscher_verify.csv'}", flush=True)


if __name__ == "__main__":
    main()
