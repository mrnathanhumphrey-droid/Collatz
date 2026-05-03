"""
Approach 2 — derive residual ~0.04 of E[v]_q125 = 2.216 via σ-quartile-
conditional residue density amplification.

Vectorized: uses np.bincount with weights for per-residue v_first averaging.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout = sys.stdout.detach() if hasattr(sys.stdout, 'detach') else sys.stdout
import io
sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)


@njit(parallel=True, cache=True)
def v_2_array(arr):
    n = len(arr)
    out = np.zeros(n, dtype=np.int32)
    for i in prange(n):
        x = arr[i]
        if x == 0:
            out[i] = 0
            continue
        v = 0
        while (x & 1) == 0:
            x >>= 1
            v += 1
        out[i] = v
    return out


def per_residue_means(res, mod, weights):
    """E[weights | res=r] for each r in [0, mod), via bincount."""
    counts = np.bincount(res, minlength=mod)
    sums = np.bincount(res, weights=weights, minlength=mod)
    means = np.where(counts > 0, sums / np.maximum(counts, 1), 0.0)
    return counts, means


def analyze_N(N, data_dir, n_sample=None):
    log2N = int(np.log2(N))
    print(f"\n========== N = 2^{log2N} = {N:,} ==========", flush=True)
    t0 = time.perf_counter()
    df = pl.read_parquet(data_dir / f"main_N{N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    if n_sample is not None and len(df) > n_sample:
        df = df.sample(n=n_sample, seed=42)
    n_arr = df["n"].to_numpy().astype(np.int64)
    sigma_arr = df["sigma"].to_numpy().astype(np.int32)
    odd_steps = df["odd_steps"].to_numpy().astype(np.int32)
    even_steps = df["even_steps"].to_numpy().astype(np.int32)
    print(f"  Loaded {len(n_arr):,} odd orbits in {time.perf_counter()-t0:.1f}s", flush=True)

    log_n = np.log(n_arr.astype(np.float64))
    sigma_f = sigma_arr.astype(np.float64)
    V_per_orbit = even_steps.astype(np.float64) / np.maximum(odd_steps.astype(np.float64), 1)

    log_n_c = log_n - log_n.mean()
    sigma_c = sigma_f - sigma_f.mean()
    beta = float((log_n_c * sigma_c).sum() / (log_n_c * log_n_c).sum())
    alpha = float(sigma_f.mean() - beta * log_n.mean())
    sigma_resid = sigma_f - (alpha + beta * log_n)
    q25 = float(np.percentile(sigma_resid, 25))
    mask_low = sigma_resid <= q25

    E_V_uncond = float(V_per_orbit.mean())
    E_V_q125 = float(V_per_orbit[mask_low].mean())
    print(f"  E[V] unconditional = {E_V_uncond:.4f}", flush=True)
    print(f"  E[V] | q125         = {E_V_q125:.4f}", flush=True)
    print(f"  Empirical V shift   = {E_V_q125 - E_V_uncond:+.4f}", flush=True)

    # First-step v_2(3m+1) per orbit, vectorized via numba
    print(f"  Computing v_2(3m+1) per orbit...", flush=True)
    t0 = time.perf_counter()
    v_first = v_2_array((3 * n_arr + 1).astype(np.int64)).astype(np.float64)
    print(f"    done in {time.perf_counter()-t0:.1f}s; mean v_first = {v_first.mean():.4f}", flush=True)

    results_per_mod = []
    for k in [5, 6, 7, 11, 17]:
        mod = 1 << k
        if mod > 2**18:
            continue
        res = (n_arr % mod).astype(np.int64)
        res_low = res[mask_low]
        v_low = v_first[mask_low]

        counts_uncond = np.bincount(res, minlength=mod)
        counts_cond = np.bincount(res_low, minlength=mod)
        sums_uncond_v = np.bincount(res, weights=v_first, minlength=mod)
        sums_cond_v = np.bincount(res_low, weights=v_low, minlength=mod)

        odd_residues = np.arange(1, mod, 2)
        cu = counts_uncond[odd_residues].astype(np.float64)
        cc = counts_cond[odd_residues].astype(np.float64)
        su = sums_uncond_v[odd_residues]
        sc = sums_cond_v[odd_residues]

        N_u = cu.sum(); N_c = cc.sum()
        P_uncond = cu / N_u
        P_cond = cc / N_c
        v_per_res_u = np.where(cu > 0, su / np.maximum(cu, 1), 0.0)
        v_per_res_c = np.where(cc > 0, sc / np.maximum(cc, 1), 0.0)

        E_v_first_uncond = float((P_uncond * v_per_res_u).sum())
        E_v_first_cond = float((P_cond * v_per_res_c).sum())
        # Decompose into "density-shift" vs "v-amplification":
        # Δ(P · v) = (P_c - P_u) · v_u + P_u · (v_c - v_u) + (P_c - P_u)·(v_c - v_u)
        # First term = pure conditional density shift assuming v unchanged per residue
        E_v_density_only = float((P_cond * v_per_res_u).sum())

        ratios = P_cond / np.maximum(P_uncond, 1e-30)
        valid = (cu >= 100) & (cc >= 50)
        if valid.sum() > 0:
            order_amp = np.argsort(-ratios * valid.astype(float))[:5]
            max_amp = float(ratios[valid].max())
        else:
            order_amp = []
            max_amp = float('nan')

        print(f"\n  --- mod 2^{k} = {mod} ---", flush=True)
        print(f"    E[v_first | uncond]      = {E_v_first_uncond:.4f}", flush=True)
        print(f"    E[v_first | bottom q]    = {E_v_first_cond:.4f}", flush=True)
        print(f"    Density-only conditional = {E_v_density_only:.4f}", flush=True)
        print(f"    Total shift = {E_v_first_cond - E_v_first_uncond:+.5f}", flush=True)
        print(f"    Density-shift contribution = {E_v_density_only - E_v_first_uncond:+.5f}", flush=True)
        print(f"    Top-5 amplified residues:", flush=True)
        print(f"      {'r':>6}  {'v_u':>6}  {'P_uncond':>10}  {'P_cond':>10}  {'ratio':>7}", flush=True)
        for idx in order_amp:
            r = int(odd_residues[idx])
            print(f"      {r:>6}  {v_per_res_u[idx]:>6.3f}  "
                  f"{P_uncond[idx]:>10.6f}  {P_cond[idx]:>10.6f}  {ratios[idx]:>7.4f}", flush=True)

        # Bucket by v
        max_v = int(v_first.max())
        if k == 17:
            print(f"\n    v-bucket decomposition (mod 2^17):", flush=True)
            print(f"      {'v':>3}  {'P_u(v)':>9}  {'P_c(v)':>9}  {'Δ P(v)':>10}  {'v · ΔP':>10}", flush=True)
            tot = 0.0
            for v_val in range(1, max_v + 1):
                mask_v = v_first == v_val
                if mask_v.sum() == 0: continue
                Pu_v = float(mask_v.sum() / len(v_first))
                Pc_v = float((mask_v & mask_low).sum() / mask_low.sum())
                contrib = v_val * (Pc_v - Pu_v)
                tot += contrib
                print(f"      {v_val:>3}  {Pu_v:>9.5f}  {Pc_v:>9.5f}  "
                      f"{Pc_v - Pu_v:>+10.5f}  {contrib:>+10.5f}", flush=True)
            print(f"      Total v_first shift = {tot:+.5f}", flush=True)

        results_per_mod.append({
            'k': k, 'mod': mod,
            'E_v_first_uncond': E_v_first_uncond,
            'E_v_first_cond': E_v_first_cond,
            'shift_first': E_v_first_cond - E_v_first_uncond,
            'shift_density_only': E_v_density_only - E_v_first_uncond,
            'max_amplification': max_amp,
        })

    return {
        'log2N': log2N, 'N': N, 'n_orbits': len(n_arr),
        'E_V_uncond': E_V_uncond, 'E_V_q125': E_V_q125,
        'V_shift_emp': E_V_q125 - E_V_uncond,
        'per_mod': results_per_mod,
    }


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "experiments_output"

    N_list = [1<<25, 1<<27, 1<<28]
    all_results = []
    for N in N_list:
        path = data_dir / f"main_N{N}.parquet"
        if not path.exists():
            print(f"# Skip 2^{int(np.log2(N))} (parquet missing)", flush=True)
            continue
        n_sample = 5_000_000 if N >= (1<<26) else None
        r = analyze_N(N, data_dir, n_sample=n_sample)
        all_results.append(r)

    print(f"\n\n========== Cross-N summary ==========", flush=True)
    print(f"  {'log2N':>6}  {'V shift (emp)':>14}  {'v_first shift mod 2^17':>24}  {'fraction':>10}", flush=True)
    for r in all_results:
        m17 = next((m for m in r['per_mod'] if m['k'] == 17), None)
        if m17 and r['V_shift_emp'] != 0:
            frac = m17['shift_first'] / r['V_shift_emp']
            print(f"  {r['log2N']:>6}  {r['V_shift_emp']:>+14.4f}  {m17['shift_first']:>+24.5f}  "
                  f"{frac*100:>9.1f}%", flush=True)
    print(f"\n  Asymptote V shift target: +0.221 (E[V]_∞ 1.995 → 2.216)", flush=True)
    print(f"  Approach 1 (joint-Gaussian) captures ~0.20 of the shift (~91%)", flush=True)
    print(f"  Approach 2 captures: see fraction above", flush=True)

    rows = []
    for r in all_results:
        for m in r['per_mod']:
            rows.append({
                'log2N': r['log2N'], 'N': r['N'], 'n_orbits': r['n_orbits'],
                'V_shift_emp': r['V_shift_emp'],
                'mod_k': m['k'], 'mod': m['mod'],
                'E_v_first_uncond': m['E_v_first_uncond'],
                'E_v_first_cond': m['E_v_first_cond'],
                'shift_first': m['shift_first'],
                'shift_density_only': m['shift_density_only'],
                'max_amplification': m['max_amplification'],
            })
    out_csv = out_dir / "57_approach2_residue_conditional.csv"
    pl.DataFrame(rows).write_csv(out_csv)
    print(f"\n[save] {out_csv}", flush=True)


if __name__ == "__main__":
    main()
