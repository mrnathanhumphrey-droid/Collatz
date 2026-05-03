"""
Return-time distribution to r ≡ 21 mod 32.

For each orbit at N=2^34, identify Syracuse-step indices where m mod 32 == 21.
Compute inter-visit gaps G_n = T_{n+1} − T_n, v_2 at each visit V_n.

Tests:
1. Marginal P(G) — geometric or other?
2. Per-σ-band P(G|band)
3. Autocorrelation ρ(G_n, G_{n+1})
4. Coupling ρ(G_n, V_n) and ρ(V_n, G_{n+1})
5. First-passage T_1 by starting residue
"""
import sys
import io
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from scipy import stats

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

LOG2 = np.log(2.0); LOG3 = np.log(3.0)
K_H = 3.0 / np.log(4.0/3.0)
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_capture_visits(starts, max_T, max_visits):
    """Walk orbits, record visits to r ≡ 21 mod 32.
    Returns flat (orbit_id, visit_idx, T_step, V_at_visit, m_mod_64)."""
    n = len(starts)
    visits = np.full((n, max_visits, 4), -1, dtype=np.int32)
    # Columns: T_step, V_at_visit, m_mod_64, prev_gap
    visit_count = np.zeros(n, dtype=np.int32)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    T_arr = np.zeros(n, dtype=np.int32)
    sumv_arr = np.zeros(n, dtype=np.int64)

    for i in prange(n):
        m = np.int64(starts[i])
        sigma = 0; T = 0; sumv = 0
        last_visit_T = -1
        vc = 0
        failed = False
        while m != 1 and T < max_T:
            if (m & 1) == 0:
                m = m >> 1; sigma += 1; continue
            if m > MAX_VAL // 3:
                failed = True; break
            r32 = m & 31  # m mod 32
            x = 3*m + 1; v = 0
            while (x & 1) == 0:
                x >>= 1; v += 1

            if r32 == 21 and vc < max_visits:
                gap = T - last_visit_T if last_visit_T >= 0 else -1
                visits[i, vc, 0] = T
                visits[i, vc, 1] = v
                visits[i, vc, 2] = m & 63  # m mod 64
                visits[i, vc, 3] = gap
                last_visit_T = T
                vc += 1

            sigma += 1 + v
            sumv += v
            T += 1
            m = x

        if not failed and m == 1:
            sigma_arr[i] = sigma
            T_arr[i] = T
            sumv_arr[i] = sumv
            visit_count[i] = vc
    return visits, visit_count, sigma_arr, T_arr, sumv_arr


def main():
    out_dir = Path("C:/Collatz/experiments_output")

    log2N = 34
    N = 1 << log2N
    n_per_seed = 50_000  # 250k orbits
    seeds = [42, 137, 271, 314, 1729]
    max_T = 500
    max_visits = 60

    print(f"# r=21 mod 32 return-time analysis at N=2^{log2N}", flush=True)
    print(f"# {len(seeds)} seeds × {n_per_seed:,} = {len(seeds)*n_per_seed:,} orbits", flush=True)

    all_visits = []; all_vc = []; all_sigma = []; all_T = []; all_logn = []; all_starts = []
    t0 = time.time()
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1)//2, size=n_per_seed, dtype=np.int64) + 1
        visits, vc, sigma, T, sumv = walk_capture_visits(starts, max_T, max_visits)
        ok = sigma > 0
        all_visits.append(visits[ok])
        all_vc.append(vc[ok])
        all_sigma.append(sigma[ok].astype(np.float64))
        all_T.append(T[ok].astype(np.float64))
        all_logn.append(np.log(starts[ok].astype(np.float64)))
        all_starts.append(starts[ok])

    visits = np.concatenate(all_visits, axis=0)
    vc = np.concatenate(all_vc)
    sigma = np.concatenate(all_sigma)
    T_total = np.concatenate(all_T)
    log_n = np.concatenate(all_logn)
    starts = np.concatenate(all_starts)
    n_orbits = len(starts)
    print(f"  walked in {time.time()-t0:.1f}s, {n_orbits:,} orbits", flush=True)

    # σ-bands
    sigma_resid = sigma - K_H * log_n
    edges = np.percentile(sigma_resid, [25, 50, 75, 95])
    band = np.digitize(sigma_resid, edges)  # 0..4

    # Build flat visit table
    print(f"  building flat visit table...", flush=True)
    flat_orbit = []; flat_visit_idx = []; flat_T = []; flat_V = []; flat_gap = []; flat_band = []; flat_mmod64 = []
    flat_logn = []; flat_starts_mod64 = []
    for i in range(n_orbits):
        nv = vc[i]
        for k in range(nv):
            flat_orbit.append(i)
            flat_visit_idx.append(k)
            flat_T.append(visits[i, k, 0])
            flat_V.append(visits[i, k, 1])
            flat_mmod64.append(visits[i, k, 2])
            flat_gap.append(visits[i, k, 3])
            flat_band.append(band[i])
            flat_logn.append(log_n[i])
            flat_starts_mod64.append(starts[i] & 63)

    flat_orbit = np.array(flat_orbit, dtype=np.int32)
    flat_visit_idx = np.array(flat_visit_idx, dtype=np.int32)
    flat_T = np.array(flat_T, dtype=np.int32)
    flat_V = np.array(flat_V, dtype=np.int32)
    flat_gap = np.array(flat_gap, dtype=np.int32)
    flat_band = np.array(flat_band, dtype=np.int32)
    flat_mmod64 = np.array(flat_mmod64, dtype=np.int32)
    flat_logn = np.array(flat_logn, dtype=np.float32)
    flat_starts_mod64 = np.array(flat_starts_mod64, dtype=np.int32)

    n_visits_total = len(flat_orbit)
    print(f"  total r=21 visits: {n_visits_total:,}", flush=True)

    # ============= Step 2: Marginal P(G) =============
    print(f"\n=== Step 2: Marginal P(G) — return-time distribution ===", flush=True)
    valid_gap = flat_gap > 0
    G = flat_gap[valid_gap]
    print(f"  Total inter-visit gaps: {len(G):,}", flush=True)
    mean_G = G.mean(); std_G = G.std()
    median_G = np.median(G); skew_G = stats.skew(G); kurt_G = stats.kurtosis(G)
    print(f"  ⟨G⟩ = {mean_G:.4f}  SD = {std_G:.4f}  median = {median_G}  skew = {skew_G:.3f}  excess kurt = {kurt_G:.3f}", flush=True)

    # P(G=g) for g=1..30
    print(f"\n  {'g':>4}  {'P(G=g)':>10}  {'count':>10}  {'log P':>10}", flush=True)
    for g in [1, 2, 3, 5, 8, 10, 15, 20, 30, 50]:
        c = int((G == g).sum())
        p = c / len(G)
        print(f"  {g:>4}  {p:>10.6f}  {c:>10,}  {np.log(p) if p>0 else float('nan'):>10.3f}", flush=True)

    # Tail decay rate: log P(G=g) vs g, fit linear
    g_vals = np.arange(1, 30)
    p_vals = np.array([(G==g).sum() / len(G) for g in g_vals])
    valid = p_vals > 1/len(G)/10  # only fit g where we have ≥0.1 expected count
    if valid.sum() > 5:
        slope, intercept, r_val, _, _ = stats.linregress(g_vals[valid], np.log(p_vals[valid]))
        print(f"\n  Linear fit log P(G=g) ≈ {intercept:+.3f} + {slope:+.4f}·g", flush=True)
        print(f"  R² = {r_val**2:.4f}", flush=True)
        print(f"  Implied geometric rate λ ≈ {-slope:.4f}", flush=True)
        print(f"  Predicted Geom mean = 1/λ = {-1/slope:.2f}  vs empirical {mean_G:.2f}", flush=True)
        # Test geometric: variance should equal mean²·(1-λ)/λ ≈ mean²(1−1/mean)
        # More precisely Geom(p) on {1, 2, ...} has mean=1/p, var=(1-p)/p²
        # So var/mean² = (1-p)/1 = 1 − p = 1 − 1/mean
        var_over_meansq = (std_G/mean_G)**2
        geom_pred = 1 - 1/mean_G
        print(f"  Var/Mean² empirical = {var_over_meansq:.4f}  vs Geom prediction (1−1/⟨G⟩) = {geom_pred:.4f}", flush=True)

    # Best-fit Geometric KS test
    p_geom = 1 / mean_G
    cdf_geom = lambda g: 1 - (1-p_geom)**g
    g_max = G.max()
    emp_cdf = np.array([(G <= g).sum() / len(G) for g in range(1, g_max+1)])
    th_cdf = np.array([cdf_geom(g) for g in range(1, g_max+1)])
    ks_geom = float(np.abs(emp_cdf - th_cdf).max())
    print(f"\n  KS distance to Geom({p_geom:.4f}): {ks_geom:.4f}", flush=True)

    # ============= Step 3: Per-σ-band P(G|band) =============
    print(f"\n=== Step 3: Per-σ-band P(G|band) ===", flush=True)
    band_names = ['0-25', '25-50', '50-75', '75-95', '95-100']
    print(f"  {'band':>6}  {'n':>10}  {'⟨G⟩':>8}  {'SD':>8}  {'P(G=1)':>9}  {'fit λ':>8}", flush=True)
    for b in range(5):
        mask = (flat_band == b) & (flat_gap > 0)
        n_b = int(mask.sum())
        if n_b < 100: continue
        G_b = flat_gap[mask]
        mean_b = G_b.mean(); sd_b = G_b.std()
        p1_b = (G_b == 1).mean()
        # Geometric rate
        lam_b = 1/mean_b
        print(f"  {band_names[b]:>6}  {n_b:>10,}  {mean_b:>8.4f}  {sd_b:>8.4f}  {p1_b:>9.4f}  {lam_b:>8.4f}", flush=True)

    # ============= Step 4: Autocorrelation ρ(G_n, G_{n+1}) =============
    print(f"\n=== Step 4: Autocorrelation of consecutive gaps ===", flush=True)
    # For each orbit, take pairs (G_n, G_{n+1}) where both are valid
    # Sort by orbit_id, visit_idx
    order = np.lexsort((flat_visit_idx, flat_orbit))
    o_sorted = flat_orbit[order]; v_sorted = flat_visit_idx[order]; g_sorted = flat_gap[order]
    # Pairs: same orbit, consecutive visit_idx, both gaps > 0
    same_orbit = o_sorted[1:] == o_sorted[:-1]
    consec = v_sorted[1:] == v_sorted[:-1] + 1
    both_gap = (g_sorted[1:] > 0) & (g_sorted[:-1] > 0)
    pair_mask = same_orbit & consec & both_gap
    G_n = g_sorted[:-1][pair_mask].astype(np.float64)
    G_np1 = g_sorted[1:][pair_mask].astype(np.float64)
    rho_lag1 = np.corrcoef(G_n, G_np1)[0, 1] if len(G_n) > 100 else float('nan')
    print(f"  Lag-1 pairs: {len(G_n):,}", flush=True)
    print(f"  ρ(G_n, G_{{n+1}}) = {rho_lag1:+.5f}", flush=True)

    # Lag 2
    same_o2 = o_sorted[2:] == o_sorted[:-2]
    consec2 = v_sorted[2:] == v_sorted[:-2] + 2
    both_g2 = (g_sorted[2:] > 0) & (g_sorted[:-2] > 0)
    pair_mask2 = same_o2 & consec2 & both_g2
    G_n_2 = g_sorted[:-2][pair_mask2].astype(np.float64)
    G_np2 = g_sorted[2:][pair_mask2].astype(np.float64)
    rho_lag2 = np.corrcoef(G_n_2, G_np2)[0, 1] if len(G_n_2) > 100 else float('nan')
    print(f"  Lag-2 pairs: {len(G_n_2):,}", flush=True)
    print(f"  ρ(G_n, G_{{n+2}}) = {rho_lag2:+.5f}", flush=True)

    # ============= Step 5: Coupling ρ(G_n, V_n) =============
    print(f"\n=== Step 5: Coupling between gap and v_2 draw ===", flush=True)
    # ρ(G_n, V_n) — does the gap leading TO visit n correlate with v_2 AT visit n?
    valid = flat_gap > 0
    G_a = flat_gap[valid].astype(np.float64); V_a = flat_V[valid].astype(np.float64)
    rho_GV = np.corrcoef(G_a, V_a)[0, 1] if len(G_a) > 100 else float('nan')
    print(f"  ρ(G_n, V_n) = {rho_GV:+.5f}  (n={len(G_a):,})", flush=True)

    # ρ(V_n, G_{n+1}) — does v_2 at visit n correlate with gap to visit n+1?
    # Use sorted arrays
    v_sorted_arr = flat_V[order]
    same_o = o_sorted[1:] == o_sorted[:-1]
    consec = v_sorted[1:] == v_sorted[:-1] + 1
    g_next_valid = g_sorted[1:] > 0
    pair_mask3 = same_o & consec & g_next_valid
    V_n_seq = v_sorted_arr[:-1][pair_mask3].astype(np.float64)
    G_np1_seq = g_sorted[1:][pair_mask3].astype(np.float64)
    rho_VG = np.corrcoef(V_n_seq, G_np1_seq)[0, 1] if len(V_n_seq) > 100 else float('nan')
    print(f"  ρ(V_n, G_{{n+1}}) = {rho_VG:+.5f}  (n={len(V_n_seq):,})", flush=True)

    # ============= Step 6: First-passage T_1 by starting residue =============
    print(f"\n=== Step 6: First-passage T_1 by starting residue mod 64 ===", flush=True)
    # For each orbit, T_1 = step of first visit (visit_idx=0)
    first_visit_mask = flat_visit_idx == 0
    T1_arr = flat_T[first_visit_mask]
    starts_mod64 = flat_starts_mod64[first_visit_mask]
    print(f"  {'r_0 mod 64':>10}  {'n':>8}  {'⟨T_1⟩':>9}  {'median T_1':>10}  {'P(T_1=0)':>10}", flush=True)
    chang_I2 = {7, 27, 31, 59, 63}
    for r0 in [3, 7, 21, 27, 31, 53, 59, 63]:
        mask = starts_mod64 == r0
        n_r = int(mask.sum())
        if n_r < 50: continue
        t1 = T1_arr[mask]
        m_t1 = float(t1.mean()); med = int(np.median(t1))
        p0 = (t1 == 0).mean()  # already at r=21 immediately
        marker = " (Chang I_2)" if r0 in chang_I2 else (" (r=21!)" if r0 == 21 else "")
        print(f"  {r0:>10}  {n_r:>8,}  {m_t1:>9.3f}  {med:>10}  {p0:>10.4f}{marker}", flush=True)

    # ============= VERDICT =============
    print(f"\n=== VERDICT ===", flush=True)
    print(f"  Marginal: ⟨G⟩ = {mean_G:.3f}, ", flush=True)
    print(f"  Geometric KS distance: {ks_geom:.4f} ({'pass' if ks_geom<0.02 else 'fail'} <0.02 threshold)", flush=True)
    print(f"  Autocorrelation ρ(G_n,G_{{n+1}}) = {rho_lag1:+.5f}", flush=True)
    print(f"  Coupling ρ(G_n,V_n) = {rho_GV:+.5f}, ρ(V_n,G_{{n+1}}) = {rho_VG:+.5f}", flush=True)
    if abs(rho_lag1) < 0.05 and ks_geom < 0.02 and abs(rho_GV) < 0.05 and abs(rho_VG) < 0.05:
        print(f"  → outcome (a): clean renewal structure", flush=True)
    elif abs(rho_lag1) < 0.05 and ks_geom > 0.02:
        print(f"  → outcome (b): independent gaps but non-Geometric — band-mixing", flush=True)
    elif abs(rho_lag1) > 0.05:
        print(f"  → outcome (c): non-renewal, gaps have memory", flush=True)
    if abs(rho_GV) > 0.05 or abs(rho_VG) > 0.05:
        print(f"  → outcome (d) flag: V/G coupling present", flush=True)

    # Save
    pl.DataFrame({
        'g': list(range(1, 51)),
        'P_G': [(G==g).sum()/len(G) for g in range(1, 51)],
    }).write_csv(out_dir / "74_r21_return_time_marginal.csv")
    print(f"\n[save] CSVs written", flush=True)


if __name__ == "__main__":
    main()
