"""
Check renewal decomposition consistency with Result 32's Cov[T, V_orbit | band].

Two independent decompositions of E[V_orbit | band]:

(σ-identity, Result 32 — algebraic, exact):
  E[V_orbit | band] = E_band-per-step − Cov[T, V_orbit | band] / E[T | band]

(renewal model, Result 47/49):
  V_orbit · T = Σ_{cylinder visits} V_n + Σ_{deterministic steps} v(residue)
  Under renewal: V_n i.i.d. given band; deterministic v determined by residue mod 64

  E[V_orbit | band] ≈ μ_d(B) + λ(B)·(E_V(B) − μ_d(B))
  where:
    μ_d(B) = orbit-weighted average v at non-cylinder steps, given band
    E_V(B) = average V at cylinder visits, given band
    λ(B)  = 1/⟨G | band⟩ = visit rate

If both forms give the same E[V_orbit | band] within bootstrap precision,
the renewal model is internally consistent with the σ-identity result.
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
def walk(starts, max_T):
    """Walk orbits; track (σ, T, sumv, K_visits, sumv_visits, sumv_det) per orbit.
    K_visits = count of r=21 mod 32 visits.
    sumv_visits = sum of v at those visits.
    sumv_det = sum of v at non-cylinder steps."""
    n = len(starts)
    sigma_arr = np.zeros(n, dtype=np.int64)
    T_arr = np.zeros(n, dtype=np.int32)
    sumv_arr = np.zeros(n, dtype=np.int64)
    K_arr = np.zeros(n, dtype=np.int32)
    sumv_visits_arr = np.zeros(n, dtype=np.int64)
    sumv_det_arr = np.zeros(n, dtype=np.int64)
    ok_arr = np.zeros(n, dtype=np.bool_)

    for i in prange(n):
        m = np.int64(starts[i])
        sigma = 0; T = 0; sumv = 0
        K = 0; sumv_v = 0; sumv_d = 0
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
            if r32 == 21:
                K += 1
                sumv_v += v
            else:
                sumv_d += v
            sigma += 1 + v
            sumv += v
            T += 1
            m = x
        if not failed and m == 1:
            sigma_arr[i] = sigma
            T_arr[i] = T
            sumv_arr[i] = sumv
            K_arr[i] = K
            sumv_visits_arr[i] = sumv_v
            sumv_det_arr[i] = sumv_d
            ok_arr[i] = True

    return sigma_arr, T_arr, sumv_arr, K_arr, sumv_visits_arr, sumv_det_arr, ok_arr


def main():
    out_dir = Path("C:/Collatz/experiments_output")

    log2N = 36
    N = 1 << log2N
    n_per_seed = 100_000
    seeds = [42, 137, 271, 314, 1729]

    print(f"# Renewal vs σ-identity consistency check at N=2^{log2N}", flush=True)
    print(f"# {len(seeds)} seeds × {n_per_seed:,} = {len(seeds)*n_per_seed:,} orbits", flush=True)

    all_starts = []; all_sigma = []; all_T = []; all_sumv = []
    all_K = []; all_sumv_v = []; all_sumv_d = []
    t0 = time.time()
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1)//2, size=n_per_seed, dtype=np.int64) + 1
        sigma, T, sumv, K_visits, sumv_v, sumv_d, ok = walk(starts, 600)
        all_starts.append(starts[ok])
        all_sigma.append(sigma[ok].astype(np.float64))
        all_T.append(T[ok].astype(np.float64))
        all_sumv.append(sumv[ok].astype(np.float64))
        all_K.append(K_visits[ok].astype(np.float64))
        all_sumv_v.append(sumv_v[ok].astype(np.float64))
        all_sumv_d.append(sumv_d[ok].astype(np.float64))

    starts = np.concatenate(all_starts)
    sigma = np.concatenate(all_sigma)
    T = np.concatenate(all_T)
    sumv = np.concatenate(all_sumv)
    K = np.concatenate(all_K)
    sumv_v = np.concatenate(all_sumv_v)
    sumv_d = np.concatenate(all_sumv_d)
    log_n = np.log(starts.astype(np.float64))
    n_orbits = len(starts)
    print(f"  walked in {time.time()-t0:.1f}s, {n_orbits:,} orbits", flush=True)

    # Per-orbit V_orbit
    V_orbit = sumv / np.maximum(T, 1)

    # σ-identity sanity
    sigma_check = T * (1 + V_orbit)
    diff_max = np.abs(sigma - sigma_check).max()
    print(f"  σ-identity check: |σ - T(1+V)|_max = {diff_max:.1e}", flush=True)

    # σ-resid bands
    sigma_resid = sigma - K_H * log_n
    edges = np.percentile(sigma_resid, [25, 50, 75, 95])
    band = np.digitize(sigma_resid, edges)
    band_names = ['0-25', '25-50', '50-75', '75-95', '95-100']

    print(f"\n=== Per-band: two decompositions of E[V_orbit | band] ===", flush=True)
    print(f"  (1) σ-identity (Result 32): E_band-per-step − Cov[T,V|B]/E[T|B]", flush=True)
    print(f"  (2) renewal:                μ_d(B) + λ(B)·(E_V(B) − μ_d(B))", flush=True)
    print(f"\n  {'band':>7}  {'⟨T⟩':>8}  {'⟨V_o⟩':>7}  {'E_per_step':>10}  "
          f"{'Cov[T,V]':>10}  {'(1) form':>9}  {'⟨K⟩':>7}  {'E_V':>7}  {'μ_d':>7}  "
          f"{'λ=K/T':>7}  {'(2) form':>9}  {'gap':>9}", flush=True)

    rows = []
    for b in range(5):
        mask = band == b
        n_b = int(mask.sum())
        if n_b < 100: continue
        T_b = T[mask]; V_b = V_orbit[mask]; sumv_b = sumv[mask]
        K_b = K[mask]; sumv_v_b = sumv_v[mask]; sumv_d_b = sumv_d[mask]

        E_T = float(T_b.mean())
        E_V_orbit = float(V_b.mean())
        E_band_per_step = float(sumv_b.sum() / T_b.sum())  # weighted average over all steps in band

        # Cov[T, V_orbit | band]
        cov_TV = float(((T_b - E_T) * (V_b - E_V_orbit)).mean())

        # (1) σ-identity form
        E_V_form1 = E_band_per_step - cov_TV / E_T

        # Renewal-model parameters
        E_K = float(K_b.mean())
        # E_V at cylinder visits (per-visit average)
        # = total sumv at visits / total visits
        sum_K = K_b.sum()
        if sum_K > 0:
            E_V_visit = float(sumv_v_b.sum() / sum_K)
        else:
            E_V_visit = 0.0
        # μ_d: average v at non-cylinder steps
        # = total sumv_det / total non-cylinder steps
        sum_det_steps = (T_b - K_b).sum()
        mu_d = float(sumv_d_b.sum() / sum_det_steps) if sum_det_steps > 0 else 0.0
        # λ = average visit rate per step
        lam = float(K_b.sum() / T_b.sum())  # use orbit-pooled λ

        # (2) renewal form
        E_V_form2 = mu_d + lam * (E_V_visit - mu_d)

        gap_12 = E_V_form1 - E_V_form2
        # Direct E[V_orbit] (per-orbit average) - Result 32 form is exact tautology with E[V_orbit]
        # Both forms should agree with E_V_orbit by construction if no error.
        gap_orbit_v_form2 = E_V_orbit - E_V_form2
        gap_orbit_v_form1 = E_V_orbit - E_V_form1

        print(f"  {band_names[b]:>7}  {E_T:>8.2f}  {E_V_orbit:>7.4f}  {E_band_per_step:>10.4f}  "
              f"{cov_TV:>+10.4f}  {E_V_form1:>9.4f}  {E_K:>7.2f}  {E_V_visit:>7.4f}  "
              f"{mu_d:>7.4f}  {lam:>7.4f}  {E_V_form2:>9.4f}  {gap_12:>+9.5f}", flush=True)

        rows.append({
            'band': band_names[b], 'n_orbits': n_b,
            'E_T': E_T, 'E_V_orbit': E_V_orbit, 'E_band_per_step': E_band_per_step,
            'cov_TV': cov_TV,
            'E_V_form1_sigma_id': E_V_form1,
            'E_K': E_K, 'E_V_visit': E_V_visit, 'mu_d': mu_d, 'lambda': lam,
            'E_V_form2_renewal': E_V_form2,
            'gap_form1_form2': gap_12,
            'gap_form1_orbitV': gap_orbit_v_form1,
            'gap_form2_orbitV': gap_orbit_v_form2,
        })

    # ============= Verdict =============
    print(f"\n=== VERDICT ===", flush=True)
    print(f"  Form (1) σ-identity vs Form (2) renewal:", flush=True)
    print(f"    Form (1) = E_band-per-step - Cov/E[T] (Result 32, exact tautology)", flush=True)
    print(f"    Form (2) = μ_d + λ·(E_V_visit - μ_d) (renewal model approximation)", flush=True)
    print(f"\n  All three E[V_orbit | band] estimates per band:", flush=True)
    print(f"  {'band':>7}  {'direct ⟨V⟩':>12}  {'(1) σ-id':>10}  {'(2) renewal':>13}  "
          f"{'(1)-direct':>11}  {'(2)-direct':>11}", flush=True)
    for r in rows:
        d1 = r['E_V_form1_sigma_id'] - r['E_V_orbit']
        d2 = r['E_V_form2_renewal'] - r['E_V_orbit']
        print(f"  {r['band']:>7}  {r['E_V_orbit']:>12.4f}  {r['E_V_form1_sigma_id']:>10.4f}  "
              f"{r['E_V_form2_renewal']:>13.4f}  {d1:>+11.5f}  {d2:>+11.5f}", flush=True)

    # Form (1) should match direct exactly (tautology); form (2) should match if renewal consistent
    max_d1 = max(abs(r['E_V_form1_sigma_id'] - r['E_V_orbit']) for r in rows)
    max_d2 = max(abs(r['E_V_form2_renewal'] - r['E_V_orbit']) for r in rows)
    print(f"\n  Max |form (1) - direct ⟨V⟩|: {max_d1:.6f}  (should be ~0 by tautology)", flush=True)
    print(f"  Max |form (2) - direct ⟨V⟩|: {max_d2:.6f}", flush=True)
    if max_d2 < 0.005:
        print(f"  → outcome (a): renewal form matches direct measurement to within 0.005", flush=True)
    elif max_d2 < 0.05:
        print(f"  → outcome (b): renewal form matches within 0.05; minor discrepancy", flush=True)
    else:
        print(f"  → outcome (c): renewal form deviates by >0.05; structural feature missing", flush=True)

    # Save
    pl.DataFrame(rows).write_csv(out_dir / "77_renewal_cov_check.csv")
    print(f"\n[save] CSV", flush=True)


if __name__ == "__main__":
    main()
