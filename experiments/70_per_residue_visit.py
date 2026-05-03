"""
Per-residue-visit dynamics test.

For each Syracuse step, record (residue_in mod 2^k, v, visit_count_so_far,
position_in_orbit). Aggregate across 200K orbits at N=2^36.

Tests:
1. P(v | r) per residue mod 32 — does it match Geom(1/2)?
2. ⟨v | r, visit=k⟩ vs k — is dynamic visit-number-independent (Markov)?
3. P(v | r, position) — is dynamic stationary along orbit?
4. Autocorrelation Cov(v_visit_i, v_visit_{i+1} | same r same orbit)
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
MAX_VAL = np.int64(2**62)

K_MOD = 32  # residue modulus = 2^5


@njit(cache=True)
def walk_one_orbit(start, max_steps, residue_buf, v_buf):
    """Walk one orbit, fill residue_buf[t] = m mod K_MOD, v_buf[t] = v at step t.
    Returns T = number of Syracuse steps if converged, -1 otherwise."""
    m = np.int64(start)
    T = 0
    while m != 1 and T < max_steps:
        if (m & 1) == 0:
            m = m >> 1
            continue
        if m > MAX_VAL // 3:
            return -1
        residue_buf[T] = m % K_MOD
        x = 3*m + 1; v = 0
        while (x & 1) == 0:
            x >>= 1; v += 1
        v_buf[T] = v if v < 127 else 127
        m = x
        T += 1
    if m == 1:
        return T
    return -1


@njit(parallel=True, cache=True)
def walk_all_orbits(starts, max_T_per_orbit):
    """Walk all orbits; return flat arrays of (orbit_id, t, T, residue, v)."""
    n = len(starts)
    # Per-orbit residue/v buffers (allocated max_T_per_orbit each)
    all_residues = np.full((n, max_T_per_orbit), -1, dtype=np.int8)
    all_v = np.full((n, max_T_per_orbit), -1, dtype=np.int8)
    T_arr = np.zeros(n, dtype=np.int32)
    for i in prange(n):
        T = walk_one_orbit(starts[i], max_T_per_orbit, all_residues[i], all_v[i])
        if T > 0:
            T_arr[i] = T
        else:
            T_arr[i] = 0  # failed
    return all_residues, all_v, T_arr


def main():
    out_dir = Path("C:/Collatz/experiments_output")

    log2N = 36
    N = 1 << log2N
    n_per_seed = 50_000  # 200k total
    seeds = [42, 137, 271, 314]
    max_T = 600  # safe upper bound

    print(f"# Per-residue-visit dynamics at N=2^{log2N}, {len(seeds)} seeds × {n_per_seed:,} = {len(seeds)*n_per_seed:,} orbits", flush=True)
    print(f"# Residue modulus = 2^5 = {K_MOD} (16 odd residues)", flush=True)

    t0 = time.time()
    all_starts = []
    all_residues = []
    all_v = []
    all_T = []

    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1)//2, size=n_per_seed, dtype=np.int64) + 1
        all_starts.append(starts)
        residues, v_arr, T_arr = walk_all_orbits(starts, max_T)
        all_residues.append(residues)
        all_v.append(v_arr)
        all_T.append(T_arr)

    starts = np.concatenate(all_starts)
    residues = np.concatenate(all_residues, axis=0)
    v_arr = np.concatenate(all_v, axis=0)
    T_arr = np.concatenate(all_T)
    print(f"  walked in {time.time()-t0:.1f}s", flush=True)

    # Filter converged orbits
    ok = T_arr > 0
    starts = starts[ok]; residues = residues[ok]; v_arr = v_arr[ok]; T_arr = T_arr[ok]
    n_orbits = len(starts)
    print(f"  ok orbits: {n_orbits:,}", flush=True)

    # ============== Build flat per-visit table ==============
    # For each (orbit, t): orbit_id, t, T, residue, v, visit_count, position
    print(f"\n# Building per-visit aggregates...", flush=True)
    t0 = time.time()
    total_visits = int(T_arr.sum())
    print(f"  total visits: {total_visits:,}", flush=True)

    # Build arrays via flat concatenation per orbit
    # For visit-count we need per-orbit running counter per residue
    flat_orbit_id = np.empty(total_visits, dtype=np.int32)
    flat_t = np.empty(total_visits, dtype=np.int16)
    flat_T = np.empty(total_visits, dtype=np.int16)
    flat_r = np.empty(total_visits, dtype=np.int8)
    flat_v = np.empty(total_visits, dtype=np.int8)
    flat_visit = np.empty(total_visits, dtype=np.int16)

    @njit(cache=True)
    def fill_flat(residues, v_arr, T_arr, flat_orbit_id, flat_t, flat_T, flat_r, flat_v, flat_visit, K_MOD):
        idx = 0
        n = len(T_arr)
        visit_count_buf = np.zeros(K_MOD, dtype=np.int16)
        for i in range(n):
            T = T_arr[i]
            for k in range(K_MOD):
                visit_count_buf[k] = 0
            for t in range(T):
                r = residues[i, t]
                v = v_arr[i, t]
                if r < 0: continue
                vc = visit_count_buf[r]
                visit_count_buf[r] = vc + 1
                flat_orbit_id[idx] = i
                flat_t[idx] = t
                flat_T[idx] = T
                flat_r[idx] = r
                flat_v[idx] = v
                flat_visit[idx] = vc + 1  # 1-indexed visit number
                idx += 1
        return idx

    actual = fill_flat(residues, v_arr, T_arr, flat_orbit_id, flat_t, flat_T, flat_r, flat_v, flat_visit, K_MOD)
    flat_orbit_id = flat_orbit_id[:actual]
    flat_t = flat_t[:actual]
    flat_T = flat_T[:actual]
    flat_r = flat_r[:actual]
    flat_v = flat_v[:actual]
    flat_visit = flat_visit[:actual]
    print(f"  built flat table in {time.time()-t0:.1f}s, rows={actual:,}", flush=True)

    flat_position = flat_t.astype(np.float32) / flat_T.astype(np.float32)

    # ============== Test 1: P(v | r) per residue ==============
    print(f"\n=== Test 1: P(v | r) per odd residue mod 32 ===", flush=True)
    odd_residues = sorted([r for r in range(1, K_MOD, 2)])
    print(f"  Geom(1/2) baseline: P(v=k) = 2^(-k), E[v]=2.0", flush=True)
    print(f"\n  {'r':>3}  {'n_visits':>9}  {'⟨v⟩':>7}  {'P(v=1)':>8}  {'P(v=2)':>8}  {'P(v=3)':>8}  "
          f"{'P(v=4)':>8}  {'gap_⟨v⟩':>9}", flush=True)

    rows_t1 = []
    for r in odd_residues:
        mask = flat_r == r
        n_v = int(mask.sum())
        if n_v < 100: continue
        v_at_r = flat_v[mask].astype(np.float64)
        mean_v = float(v_at_r.mean())
        p1 = float((v_at_r == 1).mean())
        p2 = float((v_at_r == 2).mean())
        p3 = float((v_at_r == 3).mean())
        p4 = float((v_at_r == 4).mean())
        gap = mean_v - 2.0
        print(f"  {r:>3}  {n_v:>9,}  {mean_v:>7.4f}  {p1:>8.4f}  {p2:>8.4f}  {p3:>8.4f}  "
              f"{p4:>8.4f}  {gap:>+8.4f}", flush=True)
        rows_t1.append({'r': r, 'n_visits': n_v, 'mean_v': mean_v,
                        'p1': p1, 'p2': p2, 'p3': p3, 'p4': p4, 'gap_to_2': gap})

    # Highlight {m_j} class: residues 5, 21, 85 mod ... mod 32: 5, 21
    # m_j = (4^j-1)/3: j=1 → m=1 (excluded), j=2 → m=5, j=3 → m=21, j=4 → m=85 → 85 mod 32 = 21
    # So in mod 32, m_j collapses to {5, 21}
    print(f"\n  Notable residues:", flush=True)
    for r in [5, 21]:
        rows_r = [row for row in rows_t1 if row['r'] == r]
        if rows_r:
            print(f"    r={r} (m_j sub-stratum): ⟨v⟩={rows_r[0]['mean_v']:.4f}", flush=True)

    # ============== Test 2: visit-number dependence ==============
    print(f"\n=== Test 2: ⟨v | r, visit=k⟩ across visit numbers ===", flush=True)
    print(f"  Markov on residues: ⟨v | r, visit=k⟩ should be invariant in k.", flush=True)
    print(f"\n  {'r':>3}  " + "  ".join(f"⟨v⟩@v={k}" for k in [1,2,3,4,5]) +
          "    drift v1→v5", flush=True)

    rows_t2 = []
    for r in odd_residues:
        line = f"  {r:>3}  "
        means = []
        for k in [1, 2, 3, 4, 5]:
            mask = (flat_r == r) & (flat_visit == k)
            n_k = int(mask.sum())
            if n_k < 50:
                line += f"  ----     "
                means.append(np.nan)
            else:
                m_v = float(flat_v[mask].mean())
                line += f"{m_v:>7.4f}  "
                means.append(m_v)
        if not np.isnan(means[0]) and not np.isnan(means[-1]):
            drift = means[-1] - means[0]
            line += f"   {drift:>+.4f}"
            rows_t2.append({'r': r, 'mean_v_v1': means[0], 'mean_v_v5': means[-1], 'drift': drift})
        print(line, flush=True)

    # Aggregate drift statistic
    drifts = [row['drift'] for row in rows_t2 if not np.isnan(row['drift'])]
    print(f"\n  Across all residues: mean drift v1→v5 = {np.mean(drifts):+.5f}, SD = {np.std(drifts):.5f}", flush=True)
    print(f"  Max |drift|: {max(abs(d) for d in drifts):.4f}", flush=True)

    # ============== Test 3: position-within-orbit ==============
    print(f"\n=== Test 3: ⟨v | r, position⟩ across orbit fraction ===", flush=True)
    print(f"  Stationarity: ⟨v | r, position⟩ should be invariant in position.", flush=True)
    pos_bins = [(0.0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.0)]
    print(f"\n  {'r':>3}  " + "  ".join(f"⟨v⟩@[{lo:.2f}-{hi:.2f}]" for lo,hi in pos_bins) +
          "    drift early→late", flush=True)

    rows_t3 = []
    for r in odd_residues:
        line = f"  {r:>3}  "
        means = []
        for lo, hi in pos_bins:
            mask = (flat_r == r) & (flat_position >= lo) & (flat_position < hi)
            n_p = int(mask.sum())
            if n_p < 50:
                line += f"  ----            "
                means.append(np.nan)
            else:
                m_v = float(flat_v[mask].mean())
                line += f"      {m_v:>6.4f}      "
                means.append(m_v)
        if not np.isnan(means[0]) and not np.isnan(means[-1]):
            drift = means[-1] - means[0]
            line += f"   {drift:>+.4f}"
            rows_t3.append({'r': r, 'mean_v_early': means[0], 'mean_v_late': means[-1], 'drift_pos': drift})
        print(line, flush=True)

    drifts_p = [row['drift_pos'] for row in rows_t3 if not np.isnan(row['drift_pos'])]
    print(f"\n  Across all residues: mean position drift = {np.mean(drifts_p):+.5f}, SD = {np.std(drifts_p):.5f}", flush=True)
    print(f"  Max |drift|: {max(abs(d) for d in drifts_p):.4f}", flush=True)

    # ============== Test 4: autocorrelation across visits to same residue ==============
    print(f"\n=== Test 4: autocorrelation Cov(v_i, v_{{i+1}} | same r same orbit) ===", flush=True)
    # For each orbit, partition by residue, take consecutive v pairs
    print(f"  Building visit-pairs per (orbit, residue) ...", flush=True)
    t0 = time.time()

    # Group by (orbit_id, residue), within each take v sorted by visit_count, form pairs
    # Use polars for efficient groupby
    df = pl.DataFrame({
        'orbit_id': flat_orbit_id,
        'r': flat_r,
        'v': flat_v,
        'visit': flat_visit,
    })
    # Sort by orbit_id, r, visit then take diff/lag
    df_sorted = df.sort(['orbit_id', 'r', 'visit'])
    # For pairs: shift v by 1 within (orbit_id, r) groups
    df_lag = df_sorted.with_columns([
        pl.col('v').shift(1).over(['orbit_id', 'r']).alias('v_prev'),
        pl.col('visit').shift(1).over(['orbit_id', 'r']).alias('visit_prev'),
    ])
    pairs = df_lag.filter(pl.col('v_prev').is_not_null())

    print(f"  built {len(pairs):,} consecutive-visit pairs in {time.time()-t0:.1f}s", flush=True)

    # Per-residue autocorrelation
    print(f"\n  {'r':>3}  {'n_pairs':>9}  {'⟨v⟩':>8}  {'⟨v_prev⟩':>10}  {'cov':>9}  {'corr':>9}", flush=True)
    rows_t4 = []
    for r in odd_residues:
        sub = pairs.filter(pl.col('r') == r)
        n = len(sub)
        if n < 100: continue
        v = sub['v'].to_numpy().astype(np.float64)
        vp = sub['v_prev'].to_numpy().astype(np.float64)
        cov = float(((v - v.mean()) * (vp - vp.mean())).mean())
        corr = cov / (v.std() * vp.std()) if (v.std() * vp.std()) > 1e-9 else 0.0
        print(f"  {r:>3}  {n:>9,}  {v.mean():>8.4f}  {vp.mean():>10.4f}  "
              f"{cov:>+9.5f}  {corr:>+9.5f}", flush=True)
        rows_t4.append({'r': r, 'n_pairs': n, 'mean_v': float(v.mean()),
                        'mean_v_prev': float(vp.mean()), 'cov': cov, 'corr': corr})

    # Aggregate
    corrs = [row['corr'] for row in rows_t4]
    print(f"\n  Across residues: mean autocorr = {np.mean(corrs):+.5f}, SD = {np.std(corrs):.5f}", flush=True)
    print(f"  Max |autocorr|: {max(abs(c) for c in corrs):.4f}", flush=True)

    # ============== Verdict ==============
    print(f"\n\n=== VERDICT ===", flush=True)
    max_v_drift = max(abs(d) for d in drifts)
    max_p_drift = max(abs(d) for d in drifts_p)
    max_autocorr = max(abs(c) for c in corrs)
    max_v_gap = max(abs(row['gap_to_2']) for row in rows_t1)
    print(f"  Max |⟨v|r⟩ - 2.0|       (per-residue v-deviation from Geom(1/2)): {max_v_gap:.4f}", flush=True)
    print(f"  Max |⟨v|r,v5⟩ − ⟨v|r,v1⟩| (visit-number drift):                  {max_v_drift:.4f}", flush=True)
    print(f"  Max |⟨v|r,late⟩ − ⟨v|r,early⟩| (position drift):                  {max_p_drift:.4f}", flush=True)
    print(f"  Max |Cov(v_i, v_{{i+1}} | r)| (autocorrelation):                  {max_autocorr:.4f}", flush=True)

    print(f"\n  Outcomes:", flush=True)
    if max_v_drift < 0.02 and max_p_drift < 0.02 and max_autocorr < 0.05:
        print(f"    (a) Per-residue dynamics ARE Markov (visit & position invariant, no autocorr)", flush=True)
    if max_v_drift > 0.05 or max_p_drift > 0.05:
        print(f"    (b) Visit-number/position dependence exists → hidden state", flush=True)
    if max_v_gap > 0.05:
        print(f"    (c) Per-residue ⟨v|r⟩ deviates from Geom(1/2) baseline (2.0) by {max_v_gap:.3f}", flush=True)
    if max_autocorr > 0.05:
        print(f"    (d) Autocorrelation across visits → hidden state at per-residue level", flush=True)

    # Save
    pl.DataFrame(rows_t1).write_csv(out_dir / "70_test1_p_v_given_r.csv")
    pl.DataFrame(rows_t2).write_csv(out_dir / "70_test2_visit_number.csv")
    pl.DataFrame(rows_t3).write_csv(out_dir / "70_test3_position.csv")
    pl.DataFrame(rows_t4).write_csv(out_dir / "70_test4_autocorr.csv")
    print(f"\n[save] CSVs written", flush=True)


if __name__ == "__main__":
    main()
