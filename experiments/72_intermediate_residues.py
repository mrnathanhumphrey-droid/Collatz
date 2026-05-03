"""
Boundary-non-determinism property across the 14 intermediate residues mod 64.

THEORETICAL PREDICTION (proven before computation):

For m ≡ r (mod 64) with v₀ = v_2(3r+1):
  3m+1 = 2^v₀ · u + 2⁶ · 3k   (u odd, k = (m-r)/64)
  3m+1 = 2^v₀ · (u + 2^(6−v₀) · 3k)
For v₀ < 6: 2^(6−v₀) · 3k is even, u odd → sum is odd → v_2(3m+1) = v₀ EXACTLY.
For v₀ = 6 (r=21 only): 3m+1 = 64·(u + 3k). u + 3k can be odd or even depending
on k → v_2 NON-deterministic, depends on bits of k beyond mod 64.

So: outcome (a) BY PURE ARITHMETIC. The 14 intermediate residues (v₀ ∈ {2,3,4})
plus r=53 (v₀=5) are all FULLY DETERMINISTIC at mod 64. Only r=21 is singular.

This script:
  1. Verifies the arithmetic claim numerically across mod 64, 256, 1024
  2. Computes empirical orbit observables (Δ⟨σ⟩, Δ⟨V⟩) per residue mod 64
"""
import sys
import io
import time
from pathlib import Path

import numpy as np
from numba import njit, prange
import polars as pl

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

LOG_2 = np.log(2.0); LOG_3 = np.log(3.0)
K_H = 3.0 / np.log(4.0/3.0)
MAX_VAL = np.int64(2**62)


def v2_of(n):
    """v_2 of integer n."""
    if n == 0: return -1
    v = 0
    while (n & 1) == 0:
        n >>= 1; v += 1
    return v


def step1_classify_mod64():
    """List all 32 odd residues mod 64 with v_2(3r+1)."""
    print(f"\n=== Step 1: Classification of 32 odd residues mod 64 by v_2(3r+1) ===", flush=True)
    classes = {}
    for r in range(1, 64, 2):
        v = v2_of(3*r + 1)
        classes.setdefault(v, []).append(r)
    print(f"  {'v_2':>4}  {'count':>6}  residues", flush=True)
    for v in sorted(classes.keys()):
        rs = classes[v]
        print(f"  {v:>4}  {len(rs):>6}  {rs}", flush=True)
    return classes


def step2_determinism_check(classes):
    """For each residue r mod 64, check if v_2(3m+1) is constant across mod 256 (4 lifts)
    and mod 1024 (16 lifts). PURE ARITHMETIC."""
    print(f"\n=== Step 2: Determinism check at mod 256 and mod 1024 ===", flush=True)
    print(f"  For each r mod 64, lift to m ≡ r + 64j mod 256 (j=0..3) and check v_2(3m+1).", flush=True)
    print(f"\n  {'r':>4}  {'v_2(3r+1)':>10}  {'mod 256 v_2 set':>18}  {'mod 1024 v_2 set':>22}  determinism", flush=True)
    rows = []
    for v_class in sorted(classes.keys()):
        for r in classes[v_class]:
            # mod 256 lifts: m = r + 64j, j=0..3
            v_set_256 = set()
            for j in range(4):
                m = r + 64*j
                v_set_256.add(v2_of(3*m + 1))
            # mod 1024 lifts: m = r + 64j, j=0..15
            v_set_1024 = set()
            for j in range(16):
                m = r + 64*j
                v_set_1024.add(v2_of(3*m + 1))
            det256 = "YES" if len(v_set_256) == 1 else "NO"
            det1024 = "YES" if len(v_set_1024) == 1 else "NO"
            verdict = f"{det256}/{det1024}"
            rows.append({
                'r': r, 'v0': v_class,
                'v_set_256': sorted(v_set_256), 'v_set_1024': sorted(v_set_1024),
                'det_256': det256, 'det_1024': det1024,
            })
            v256_str = str(sorted(v_set_256))
            v1024_str = str(sorted(v_set_1024))
            print(f"  {r:>4}  {v_class:>10}  {v256_str:>18}  {v1024_str:>22}  {verdict}", flush=True)
    return rows


def step3_classification_summary(rows, classes):
    print(f"\n=== Step 3: Classification summary ===", flush=True)
    intermediate_set = []
    for v_class in [2, 3, 4]:  # the 14 intermediate residues
        intermediate_set.extend(classes.get(v_class, []))
    # plus v=5 cluster (just r=53)
    v5_set = classes.get(5, [])
    v6_set = classes.get(6, [])
    chang_I2 = [7, 27, 31, 59, 63]

    print(f"\n  Chang I_2 = {chang_I2} (subset of v=1, 5 residues)", flush=True)
    print(f"  Intermediate (v_2 ∈ {{2,3,4}}) = {sorted(intermediate_set)} ({len(intermediate_set)} residues)", flush=True)
    print(f"  v_2 = 5 cluster = {sorted(v5_set)} ({len(v5_set)} residue)", flush=True)
    print(f"  v_2 = 6 cluster = {sorted(v6_set)} ({len(v6_set)} residue, BOUNDARY)", flush=True)

    # Determinism check
    intermediate_rows = [r for r in rows if r['r'] in intermediate_set]
    n_det_256 = sum(1 for r in intermediate_rows if r['det_256'] == 'YES')
    n_det_1024 = sum(1 for r in intermediate_rows if r['det_1024'] == 'YES')
    print(f"\n  Intermediate residues deterministic at mod 256: {n_det_256}/{len(intermediate_rows)}", flush=True)
    print(f"  Intermediate residues deterministic at mod 1024: {n_det_1024}/{len(intermediate_rows)}", flush=True)

    v5_rows = [r for r in rows if r['r'] in v5_set]
    n_det_v5 = sum(1 for r in v5_rows if r['det_256'] == 'YES')
    print(f"  v=5 residue (r=53) deterministic at mod 256: {n_det_v5}/{len(v5_rows)}", flush=True)

    v6_rows = [r for r in rows if r['r'] in v6_set]
    n_det_v6 = sum(1 for r in v6_rows if r['det_256'] == 'YES')
    print(f"  v=6 residue (r=21) deterministic at mod 256: {n_det_v6}/{len(v6_rows)}", flush=True)

    return {
        'chang_I2': chang_I2, 'intermediate_set': sorted(intermediate_set),
        'v5_set': sorted(v5_set), 'v6_set': sorted(v6_set),
        'n_det_intermediate_256': n_det_256, 'n_total_intermediate': len(intermediate_rows),
    }


@njit(parallel=True, cache=True)
def walk_capture_residue(starts, max_steps):
    n = len(starts)
    sigma_arr = np.zeros(n, dtype=np.int32)
    T_arr = np.zeros(n, dtype=np.int32)
    sumv_arr = np.zeros(n, dtype=np.int64)
    ok_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = np.int64(starts[i])
        sigma_total = 0; T = 0; sumv = 0
        failed = False
        while m != 1 and T < max_steps:
            if (m & 1) == 0:
                m = m >> 1; sigma_total += 1; continue
            if m > MAX_VAL // 3:
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


def step5_orbit_observables(classes):
    """At N=2^32, walk orbits stratified by starting residue mod 64.
    Compute mean σ, mean V_orbit per residue."""
    print(f"\n=== Step 5: Empirical orbit observables per residue mod 64 ===", flush=True)
    log2N = 32
    N = 1 << log2N
    n_per_seed = 200_000  # 1M total
    seeds = [42, 137, 271, 314, 1729]

    print(f"  Walking {len(seeds)*n_per_seed:,} orbits at N=2^{log2N}", flush=True)
    t0 = time.time()
    all_starts = []; all_sigma = []; all_V = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1)//2, size=n_per_seed, dtype=np.int64) + 1
        sigma, T, sumv, ok = walk_capture_residue(starts, 600)
        starts_ok = starts[ok]
        sigma_ok = sigma[ok].astype(np.float64)
        T_ok = T[ok].astype(np.float64)
        sumv_ok = sumv[ok].astype(np.float64)
        V = sumv_ok / np.maximum(T_ok, 1)
        all_starts.append(starts_ok); all_sigma.append(sigma_ok); all_V.append(V)
    starts = np.concatenate(all_starts)
    sigma = np.concatenate(all_sigma)
    V = np.concatenate(all_V)
    print(f"  walked in {time.time()-t0:.1f}s", flush=True)

    starts_mod64 = (starts % 64).astype(np.int32)
    n_orbits = len(starts)

    # Per-residue stats
    print(f"\n  {'r':>4}  {'v0':>3}  {'n_orbits':>9}  {'⟨σ⟩':>9}  {'⟨V⟩':>7}  {'group':>15}", flush=True)
    rows = []
    chang_I2 = {7, 27, 31, 59, 63}
    intermediate_set = set()
    for v in [2, 3, 4]:
        intermediate_set.update(classes.get(v, []))
    v5_set = set(classes.get(5, []))
    v6_set = set(classes.get(6, []))

    for v_class in sorted(classes.keys()):
        for r in classes[v_class]:
            mask = starts_mod64 == r
            n_r = int(mask.sum())
            if n_r < 100: continue
            mean_sigma = float(sigma[mask].mean())
            mean_V = float(V[mask].mean())
            if r in chang_I2:
                grp = "Chang I_2"
            elif r in intermediate_set:
                grp = "intermediate"
            elif r in v5_set:
                grp = "v=5"
            elif r in v6_set:
                grp = "v=6 BOUNDARY"
            else:
                grp = "v=1 other"
            rows.append({'r': r, 'v0': v_class, 'n_orbits': n_r,
                         'mean_sigma': mean_sigma, 'mean_V': mean_V, 'group': grp})
            print(f"  {r:>4}  {v_class:>3}  {n_r:>9,}  {mean_sigma:>9.2f}  {mean_V:>7.4f}  {grp:>15}", flush=True)

    # Group means
    print(f"\n=== Group-level summary ===", flush=True)
    groups = {}
    for row in rows:
        groups.setdefault(row['group'], []).append(row)
    print(f"  {'group':>15}  n_residues  {'⟨σ⟩':>9}  {'SD(⟨σ⟩)':>9}  {'⟨V⟩':>7}  {'SD(⟨V⟩)':>9}", flush=True)
    grp_summary = []
    for g, rs in sorted(groups.items()):
        sigmas = np.array([r['mean_sigma'] for r in rs])
        Vs = np.array([r['mean_V'] for r in rs])
        print(f"  {g:>15}  {len(rs):>10}  {sigmas.mean():>9.2f}  {sigmas.std():>9.3f}  "
              f"{Vs.mean():>7.4f}  {Vs.std():>9.4f}", flush=True)
        grp_summary.append({'group': g, 'n_residues': len(rs),
                            'mean_sigma': float(sigmas.mean()), 'sd_sigma': float(sigmas.std()),
                            'mean_V': float(Vs.mean()), 'sd_V': float(Vs.std())})

    # Δ vs generic
    generic_rs = [row for row in rows if row['group'] in ['Chang I_2', 'v=1 other', 'intermediate']]
    overall_sigma = np.mean([r['mean_sigma'] for row in [groups['Chang I_2'] + groups['v=1 other'] + groups['intermediate']] for r in row]) if all(g in groups for g in ['Chang I_2','v=1 other','intermediate']) else None

    return rows, grp_summary


def main():
    out_dir = Path("C:/Collatz/experiments_output")

    classes = step1_classify_mod64()
    rows = step2_determinism_check(classes)
    summary = step3_classification_summary(rows, classes)

    # Save the determinism table
    det_rows = [{'r': r['r'], 'v0': r['v0'],
                 'v_set_256': str(r['v_set_256']), 'v_set_1024': str(r['v_set_1024']),
                 'det_256': r['det_256'], 'det_1024': r['det_1024']} for r in rows]
    pl.DataFrame(det_rows).write_csv(out_dir / "72_determinism_mod64.csv")

    obs_rows, grp_summary = step5_orbit_observables(classes)
    pl.DataFrame(obs_rows).write_csv(out_dir / "72_orbit_observables.csv")
    pl.DataFrame(grp_summary).write_csv(out_dir / "72_group_summary.csv")

    # Verdict
    print(f"\n=== VERDICT ===", flush=True)
    if summary['n_det_intermediate_256'] == summary['n_total_intermediate']:
        print(f"  Outcome (a): ALL 14 intermediate residues deterministic at mod 256 ✓", flush=True)
        print(f"  Boundary-non-determinism is unique to r=21 mod 64 (v_2=6).", flush=True)
        print(f"  The v_2 spectrum is structurally homogeneous (deterministic) except at v_2 ≥ 6 boundary.", flush=True)
    else:
        n_nondet = summary['n_total_intermediate'] - summary['n_det_intermediate_256']
        print(f"  Outcome (b): {n_nondet}/{summary['n_total_intermediate']} intermediate residues non-deterministic", flush=True)

    print(f"\n[save] CSVs written", flush=True)


if __name__ == "__main__":
    main()
