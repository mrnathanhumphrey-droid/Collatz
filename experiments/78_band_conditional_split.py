"""
Verify band-conditional E_V deviations decompose as predicted m ≡ 21/53 mod 64 split.

Per Result 50/52:
  m ≡ 53 mod 64 (V=5 deterministic) and m ≡ 21 mod 64 (V≥6 shifted Geom(1/2))
  E[V | r=21 mod 32, band] = 7·p_21(B) + 5·p_53(B) = 5 + 2·p_21(B)

Test: does empirical E_V match this prediction across bands?
If yes (outcome a), open piece reduces to 1D function p_21(B).
"""
import sys
import io
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

LOG_2 = np.log(2.0); LOG_3 = np.log(3.0)
K_H = 3.0 / np.log(4.0/3.0)
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_capture(starts, max_T, max_visits):
    """Walk orbits, record (T_step, V, m_mod_64, gap) per r=21 mod 32 visit."""
    n = len(starts)
    visits = np.full((n, max_visits, 4), -1, dtype=np.int32)
    visit_count = np.zeros(n, dtype=np.int32)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    T_arr = np.zeros(n, dtype=np.int32)

    for i in prange(n):
        m = np.int64(starts[i])
        sigma = 0; T = 0
        last_visit_T = -1
        vc = 0
        failed = False
        while m != 1 and T < max_T:
            if (m & 1) == 0:
                m = m >> 1; sigma += 1; continue
            if m > MAX_VAL // 3:
                failed = True; break
            r32 = m & 31
            x = 3*m + 1; v = 0
            while (x & 1) == 0:
                x >>= 1; v += 1
            if r32 == 21 and vc < max_visits:
                gap = T - last_visit_T if last_visit_T >= 0 else -1
                visits[i, vc, 0] = T
                visits[i, vc, 1] = v
                visits[i, vc, 2] = m & 63  # m mod 64 (must be 21 or 53)
                visits[i, vc, 3] = gap
                last_visit_T = T
                vc += 1
            sigma += 1 + v
            T += 1
            m = x
        if not failed and m == 1:
            sigma_arr[i] = sigma
            T_arr[i] = T
            visit_count[i] = vc
    return visits, visit_count, sigma_arr, T_arr


def main():
    out_dir = Path("C:/Collatz/experiments_output")

    log2N = 36
    N = 1 << log2N
    n_per_seed = 100_000
    seeds = [42, 137, 271, 314, 1729]

    print(f"# Band-conditional m ≡ 21/53 mod 64 split at N=2^{log2N}", flush=True)

    all_visits = []; all_vc = []; all_sigma = []; all_logn = []
    t0 = time.time()
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1)//2, size=n_per_seed, dtype=np.int64) + 1
        visits, vc, sigma, T = walk_capture(starts, 600, 60)
        ok = sigma > 0
        all_visits.append(visits[ok])
        all_vc.append(vc[ok])
        all_sigma.append(sigma[ok].astype(np.float64))
        all_logn.append(np.log(starts[ok].astype(np.float64)))

    visits = np.concatenate(all_visits, axis=0)
    vc = np.concatenate(all_vc)
    sigma = np.concatenate(all_sigma)
    log_n = np.concatenate(all_logn)
    n_orbits = len(sigma)
    print(f"  walked in {time.time()-t0:.1f}s, {n_orbits:,} orbits", flush=True)

    # σ-bands
    sigma_resid = sigma - K_H * log_n
    edges = np.percentile(sigma_resid, [25, 50, 75, 95])
    band = np.digitize(sigma_resid, edges)

    # Build flat visit table
    flat_band = []; flat_V = []; flat_mmod64 = []; flat_gap = []
    for i in range(n_orbits):
        for k in range(vc[i]):
            flat_band.append(band[i])
            flat_V.append(visits[i, k, 1])
            flat_mmod64.append(visits[i, k, 2])
            flat_gap.append(visits[i, k, 3])
    flat_band = np.array(flat_band, dtype=np.int8)
    flat_V = np.array(flat_V, dtype=np.int8)
    flat_mmod64 = np.array(flat_mmod64, dtype=np.int8)
    flat_gap = np.array(flat_gap, dtype=np.int32)
    n_visits = len(flat_band)
    print(f"  visits: {n_visits:,}", flush=True)

    # Sanity: m mod 64 ∈ {21, 53}
    assert set(np.unique(flat_mmod64)).issubset({21, 53}), f"unexpected residues: {np.unique(flat_mmod64)}"

    # ============= Step 1+3: per-band p_21 and predicted vs empirical E_V =============
    print(f"\n=== Step 1-3: Band-conditional p_21 vs empirical E_V ===", flush=True)
    print(f"  Predicted E_V = 5 + 2·p_21(B) (under closed-form Result 50)", flush=True)
    print(f"\n  {'band':>7}  {'n':>10}  {'p_21':>7}  {'p_53':>7}  "
          f"{'pred_E_V':>9}  {'emp_E_V':>9}  {'gap':>9}", flush=True)
    band_names = ['0-25', '25-50', '50-75', '75-95', '95-100']
    rows_main = []
    for b in range(5):
        mask = flat_band == b
        n_b = int(mask.sum())
        if n_b < 100: continue
        m64 = flat_mmod64[mask]
        v_b = flat_V[mask]
        p_21 = float((m64 == 21).mean())
        p_53 = float((m64 == 53).mean())
        pred_E_V = 5 + 2*p_21
        emp_E_V = float(v_b.mean())
        gap = pred_E_V - emp_E_V
        # bootstrap SE on emp_E_V
        se_emp = float(v_b.std() / np.sqrt(n_b))
        print(f"  {band_names[b]:>7}  {n_b:>10,}  {p_21:>7.4f}  {p_53:>7.4f}  "
              f"{pred_E_V:>9.4f}  {emp_E_V:>9.4f}  {gap:>+9.5f}", flush=True)
        rows_main.append({
            'band': band_names[b], 'n_visits': n_b,
            'p_21': p_21, 'p_53': p_53,
            'pred_E_V': pred_E_V, 'emp_E_V': emp_E_V,
            'gap': gap, 'se_emp': se_emp,
        })

    # Verdict on accuracy
    max_gap = max(abs(r['gap']) for r in rows_main)
    print(f"\n  Max |pred − emp|: {max_gap:.5f}", flush=True)
    if max_gap < 0.005:
        print(f"  → outcome (a) within 0.005 — mechanism fully verified", flush=True)
    elif max_gap < 0.05:
        print(f"  → outcome (b) within 0.05 — mostly verified, small residual", flush=True)
    else:
        print(f"  → outcome (c) — mechanism doesn't fully account for E_V deviations", flush=True)

    # ============= Step 4: Verify within m ≡ 21 mod 64, V−6 ~ Geom(1/2) =============
    print(f"\n=== Step 4: Within m ≡ 21 mod 64, verify V−6 ~ Geom(1/2) ===", flush=True)
    print(f"  Predicted P(V=6+k | m ≡ 21 mod 64) = 2^(-(k+1))", flush=True)
    print(f"\n  {'band':>7}  {'P(V=6)':>8}  {'P(V=7)':>8}  {'P(V=8)':>8}  {'P(V=9)':>8}  "
          f"{'P(V≥10)':>9}  {'⟨V⟩-6':>7}", flush=True)
    print(f"  {'pred':>7}  {'0.5000':>8}  {'0.2500':>8}  {'0.1250':>8}  {'0.0625':>8}  "
          f"{'0.0625':>9}  {'1.0000':>7}", flush=True)
    for b in range(5):
        mask = (flat_band == b) & (flat_mmod64 == 21)
        n_b = int(mask.sum())
        if n_b < 100: continue
        v_b = flat_V[mask]
        p6 = (v_b == 6).mean(); p7 = (v_b == 7).mean(); p8 = (v_b == 8).mean()
        p9 = (v_b == 9).mean(); p10p = (v_b >= 10).mean()
        mean_minus_6 = float(v_b.mean()) - 6
        print(f"  {band_names[b]:>7}  {p6:>8.4f}  {p7:>8.4f}  {p8:>8.4f}  {p9:>8.4f}  "
              f"{p10p:>9.4f}  {mean_minus_6:>7.4f}", flush=True)

    # ============= Step 5: p_21 across bands — open piece signature =============
    print(f"\n=== Step 5: p_21(B) curve — the 1D open-piece signature ===", flush=True)
    p_21_vals = [r['p_21'] for r in rows_main]
    print(f"  Uniform prediction: p_21 = 0.5 in every band", flush=True)
    print(f"  Empirical range: {min(p_21_vals):.4f} to {max(p_21_vals):.4f}", flush=True)
    for r in rows_main:
        deviation = r['p_21'] - 0.5
        print(f"  {r['band']:>7}  p_21 = {r['p_21']:.4f}  (deviation {deviation:+.4f})", flush=True)

    # ============= Step 6: Test p_21(B, G) — does gap add information beyond band? =============
    print(f"\n=== Step 6: p_21 vs gap (G-coupling test) ===", flush=True)
    valid_g = flat_gap > 0
    # Stratify by gap deciles within each band
    print(f"  For each band, stratify visits by gap; p_21 should DECREASE with G", flush=True)
    print(f"  (long gap → more bit-mixing → land on r=53 mod 64 → V=5)", flush=True)
    print(f"\n  band 50-75 (representative middle band):", flush=True)
    mask_b = (flat_band == 2) & valid_g
    g_in_band = flat_gap[mask_b]
    m64_in_band = flat_mmod64[mask_b]
    g_quantiles = np.percentile(g_in_band, [0, 20, 40, 60, 80, 100])
    print(f"  {'gap range':>15}  {'n':>8}  {'p_21':>7}", flush=True)
    for q_idx in range(5):
        lo, hi = g_quantiles[q_idx], g_quantiles[q_idx+1]
        m_g = (g_in_band >= lo) & (g_in_band <= hi)
        n_q = int(m_g.sum())
        if n_q == 0: continue
        p21_q = float((m64_in_band[m_g] == 21).mean())
        print(f"  [{lo:>4.0f}, {hi:>4.0f}]    {n_q:>8,}  {p21_q:>7.4f}", flush=True)

    # ============= Step 6b: pooled p_21 vs gap across all bands =============
    print(f"\n=== Step 6b: Pooled p_21 vs gap (across bands) ===", flush=True)
    valid_g = flat_gap > 0
    g_all = flat_gap[valid_g]
    m64_all = flat_mmod64[valid_g]
    g_qs = np.percentile(g_all, [0, 20, 40, 60, 80, 100])
    print(f"  {'gap range':>15}  {'n':>10}  {'p_21':>7}", flush=True)
    for q_idx in range(5):
        lo, hi = g_qs[q_idx], g_qs[q_idx+1]
        m_g = (g_all >= lo) & (g_all <= hi)
        n_q = int(m_g.sum())
        p21_q = float((m64_all[m_g] == 21).mean())
        print(f"  [{lo:>4.0f}, {hi:>4.0f}]    {n_q:>10,}  {p21_q:>7.4f}", flush=True)

    # Save
    pl.DataFrame(rows_main).write_csv(out_dir / "78_band_conditional_split.csv")
    print(f"\n[save] CSV", flush=True)


if __name__ == "__main__":
    main()
