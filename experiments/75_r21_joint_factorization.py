"""
Joint observable factorization at r=21 visits.

Tests the central question: given σ-band, are V values at successive r=21 visits
independent (cylinder = i.i.d.) or do they have memory (Markov / higher-order)?

Builds on Result 47 (return-time) which established:
  - Renewal at gap level: ρ(G_n, G_{n+1}) ≈ 0
  - Asymmetric G→V coupling: ρ(G_n, V_n) = -0.139

This test computes:
  - V autocorrelation: ρ(V_n, V_{n+1}) — Markov order test
  - Empirical Markov kernel K(V_{n+1} | V_n, B)
  - Mutual information I(V_n, V_{n+1} | B) and I(V_{n+1} | V_n, V_{n-1}, B)
  - Joint factorization I(V_n, m_higher_n | B) vs I(V_n, V_{n-1} | B)
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
    """Walk orbits, record at each r=21 visit: T_step, V at visit, exit residue mod 64,
    m mod 4096 at visit."""
    n = len(starts)
    visits = np.full((n, max_visits, 5), -1, dtype=np.int32)
    # Columns: T_step, V, exit_mod64, m_mod_4096, gap
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
                visits[i, vc, 2] = x & 63   # exit residue mod 64 (m' after this step)
                visits[i, vc, 3] = m & 4095  # m mod 4096 (high bits of m at visit)
                visits[i, vc, 4] = gap
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


def entropy(p):
    p = p[p > 0]
    return float(-np.sum(p * np.log2(p)))


def mutual_info(joint):
    """Mutual information from 2D joint distribution joint[i,j]."""
    joint = joint / joint.sum() if joint.sum() > 0 else joint
    p_x = joint.sum(axis=1)
    p_y = joint.sum(axis=0)
    H_x = entropy(p_x); H_y = entropy(p_y)
    H_xy = entropy(joint.flatten())
    return H_x + H_y - H_xy


def main():
    out_dir = Path("C:/Collatz/experiments_output")

    log2N = 34
    N = 1 << log2N
    n_per_seed = 50_000
    seeds = [42, 137, 271, 314, 1729]

    print(f"# Joint factorization at r=21 visits, N=2^{log2N}", flush=True)
    print(f"# {len(seeds)} seeds × {n_per_seed:,} = {len(seeds)*n_per_seed:,} orbits", flush=True)

    all_visits = []; all_vc = []; all_sigma = []; all_logn = []
    t0 = time.time()
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1)//2, size=n_per_seed, dtype=np.int64) + 1
        visits, vc, sigma, T = walk_capture(starts, 500, 60)
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
    print(f"  walked in {time.time()-t0:.1f}s", flush=True)

    # σ-bands
    sigma_resid = sigma - K_H * log_n
    edges = np.percentile(sigma_resid, [25, 50, 75, 95])
    band = np.digitize(sigma_resid, edges)

    # Build flat visit table with previous-visit info
    print(f"  building flat visit table with predecessor info...", flush=True)
    flat_orbit = []; flat_visit_idx = []; flat_T = []
    flat_V = []; flat_exit = []; flat_mhigh = []; flat_gap = []; flat_band = []
    flat_V_prev = []; flat_V_prev2 = []; flat_gap_prev = []
    for i in range(n_orbits):
        nv = vc[i]
        for k in range(nv):
            flat_orbit.append(i)
            flat_visit_idx.append(k)
            flat_T.append(visits[i, k, 0])
            flat_V.append(visits[i, k, 1])
            flat_exit.append(visits[i, k, 2])
            flat_mhigh.append(visits[i, k, 3])
            flat_gap.append(visits[i, k, 4])
            flat_band.append(band[i])
            flat_V_prev.append(visits[i, k-1, 1] if k >= 1 else -1)
            flat_V_prev2.append(visits[i, k-2, 1] if k >= 2 else -1)
            flat_gap_prev.append(visits[i, k-1, 4] if k >= 1 else -1)

    flat_V = np.array(flat_V, dtype=np.int8)
    flat_V_prev = np.array(flat_V_prev, dtype=np.int8)
    flat_V_prev2 = np.array(flat_V_prev2, dtype=np.int8)
    flat_exit = np.array(flat_exit, dtype=np.int32)
    flat_mhigh = np.array(flat_mhigh, dtype=np.int32)
    flat_gap = np.array(flat_gap, dtype=np.int32)
    flat_gap_prev = np.array(flat_gap_prev, dtype=np.int32)
    flat_band = np.array(flat_band, dtype=np.int8)
    flat_visit_idx = np.array(flat_visit_idx, dtype=np.int32)

    n_visits = len(flat_V)
    print(f"  visits: {n_visits:,}", flush=True)

    # ============= Step 3: V_n autocorrelation =============
    print(f"\n=== Step 3a: V_n autocorrelation across consecutive r=21 visits ===", flush=True)
    valid = (flat_V_prev >= 5) & (flat_V >= 5)  # both must be valid
    V_a = flat_V[valid].astype(np.float64); V_p = flat_V_prev[valid].astype(np.float64)
    rho_VV = np.corrcoef(V_a, V_p)[0, 1] if len(V_a) > 100 else float('nan')
    print(f"  Pairs (V_n, V_{{n-1}}): {len(V_a):,}", flush=True)
    print(f"  ρ(V_n, V_{{n-1}}) = {rho_VV:+.5f}", flush=True)

    # Lag 2
    valid2 = (flat_V_prev2 >= 5) & (flat_V >= 5)
    V_a2 = flat_V[valid2].astype(np.float64); V_p2 = flat_V_prev2[valid2].astype(np.float64)
    rho_VV2 = np.corrcoef(V_a2, V_p2)[0, 1] if len(V_a2) > 100 else float('nan')
    print(f"  ρ(V_n, V_{{n-2}}) = {rho_VV2:+.5f}", flush=True)

    # ============= Step 3b: Per-band V autocorrelation =============
    print(f"\n=== Step 3b: V autocorrelation per σ-band ===", flush=True)
    band_names = ['0-25', '25-50', '50-75', '75-95', '95-100']
    print(f"  {'band':>6}  {'n':>10}  {'ρ(V_n,V_{{n-1}})':>18}  {'I(V_n;V_{{n-1}})':>17}", flush=True)
    for b in range(5):
        mask = (flat_band == b) & valid
        n_b = int(mask.sum())
        if n_b < 100: continue
        v_b = flat_V[mask].astype(np.float64); vp_b = flat_V_prev[mask].astype(np.float64)
        rho_b = np.corrcoef(v_b, vp_b)[0, 1] if v_b.std() > 0 and vp_b.std() > 0 else 0.0
        # MI via histogram (V values 5..16)
        v_clip = np.clip(flat_V[mask], 5, 16) - 5
        vp_clip = np.clip(flat_V_prev[mask], 5, 16) - 5
        joint = np.zeros((12, 12))
        for v_, vp_ in zip(v_clip, vp_clip):
            joint[v_, vp_] += 1
        mi = mutual_info(joint)
        print(f"  {band_names[b]:>6}  {n_b:>10,}  {rho_b:>+17.5f}  {mi:>17.5f}", flush=True)

    # ============= Step 5: I(V_n, m_higher_n | B) =============
    print(f"\n=== Step 5: I(V_n; m_higher_n | B) — V is determined by arithmetic of m's higher bits ===", flush=True)
    # m_higher = m mod 4096 (= 21 + 32k for some k, k in 0..127)
    # V is determined by k via v_2(2 + 3k) + 5
    # So I(V; m_higher) should be high (close to H(V) = ~1.9 bits)
    print(f"  {'band':>6}  {'I(V;m_high)':>14}  {'H(V|band)':>11}  {'fraction explained':>20}", flush=True)
    for b in range(5):
        mask = flat_band == b
        if mask.sum() < 100: continue
        # bin m_higher (mod 4096, but only m ≡ 21 mod 32 → 128 distinct residues)
        # Just use m_higher values as bins
        v_clip = np.clip(flat_V[mask], 5, 16) - 5
        m_clip = (flat_mhigh[mask] % 4096) // 32  # 128 bins (since m ≡ 21 mod 32)
        # Build joint
        n_v = 12; n_m = 128
        joint = np.zeros((n_v, n_m))
        for v_, m_ in zip(v_clip, m_clip):
            if 0 <= m_ < n_m:
                joint[v_, m_] += 1
        I_Vm = mutual_info(joint)
        p_v = joint.sum(axis=1); p_v = p_v / p_v.sum()
        H_V = entropy(p_v)
        frac = I_Vm / H_V if H_V > 0 else 0.0
        print(f"  {band_names[b]:>6}  {I_Vm:>14.5f}  {H_V:>11.5f}  {frac:>20.5f}", flush=True)
    print(f"\n  Expected: I(V; m_high) ≈ H(V|band) (V deterministic from m_high arithmetic)", flush=True)

    # ============= Step 5b: I(V_n; V_{n-1} | B, m_higher_{n-1}) — does V_{n-1} add beyond m_higher_{n-1}? =============
    print(f"\n=== Step 5b: I(V_n; V_{{n-1}} | B) compared to I(V_n; G_n | B) ===", flush=True)
    # V_n autocorrelation we have. Now compare to gap-coupling.
    print(f"  {'band':>6}  {'I(V_n;V_{{n-1}})':>17}  {'I(V_n;G_n)':>14}", flush=True)
    valid_g = flat_gap > 0
    for b in range(5):
        mask = (flat_band == b) & valid
        if mask.sum() < 100: continue
        v_clip = np.clip(flat_V[mask], 5, 16) - 5
        vp_clip = np.clip(flat_V_prev[mask], 5, 16) - 5
        joint_VV = np.zeros((12, 12))
        for v_, vp_ in zip(v_clip, vp_clip):
            joint_VV[v_, vp_] += 1
        I_VV = mutual_info(joint_VV)

        mask_g = (flat_band == b) & valid_g
        v_clip2 = np.clip(flat_V[mask_g], 5, 16) - 5
        g_clip = np.clip(flat_gap[mask_g], 1, 50) - 1  # 50 bins
        joint_VG = np.zeros((12, 50))
        for v_, g_ in zip(v_clip2, g_clip):
            joint_VG[v_, g_] += 1
        I_VG = mutual_info(joint_VG)
        print(f"  {band_names[b]:>6}  {I_VV:>17.5f}  {I_VG:>14.5f}", flush=True)

    # ============= Step 6: Empirical Markov kernel K(V_{n+1} | V_n, B) =============
    print(f"\n=== Step 6: Empirical Markov kernel K(V_{{n+1}} | V_n, B=middle) ===", flush=True)
    # Use middle band (50-75) for clarity
    b = 2  # 50-75
    mask = (flat_band == b) & valid
    v_clip = np.clip(flat_V[mask], 5, 16) - 5
    vp_clip = np.clip(flat_V_prev[mask], 5, 16) - 5
    K = np.zeros((12, 12))
    for v_, vp_ in zip(v_clip, vp_clip):
        K[vp_, v_] += 1  # rows = V_{n-1}, cols = V_n
    # Normalize rows to get K(V_n | V_{n-1})
    K_norm = K / np.maximum(K.sum(axis=1, keepdims=True), 1)
    print(f"  Row = V_{{n-1}}, Col = V_n, K[i,j] = P(V_n=j | V_{{n-1}}=i)", flush=True)
    print(f"  V values (offset +5):", flush=True)
    print(f"  {'V_n-1 \\ V_n':>14}  " + "  ".join(f"V={k+5:>2}" for k in range(7)), flush=True)
    for i in range(7):
        line = f"  V_{{n-1}}={i+5:>2}        "
        for j in range(7):
            line += f"  {K_norm[i,j]:>5.3f}"
        print(line, flush=True)

    # Eigenvalues
    if K.sum() > 0:
        eigvals = np.linalg.eigvals(K_norm)
        eigvals_real = np.sort(np.abs(eigvals))[::-1]
        print(f"\n  Top 5 |eigenvalues| of K_norm: {eigvals_real[:5]}", flush=True)
        if len(eigvals_real) >= 2:
            print(f"  Spectral gap (1 − λ_2): {1 - eigvals_real[1]:.4f}", flush=True)

    # Marginal entropy comparison
    p_V_marg = K.sum(axis=0); p_V_marg /= p_V_marg.sum()
    H_marg = entropy(p_V_marg)
    H_cond = 0.0
    p_Vp = K.sum(axis=1); p_Vp /= p_Vp.sum()
    for i in range(12):
        if K[i].sum() == 0: continue
        p_v_given_vp = K[i] / K[i].sum()
        H_cond += p_Vp[i] * entropy(p_v_given_vp)
    print(f"\n  H(V_n | band=middle): {H_marg:.4f} bits", flush=True)
    print(f"  H(V_n | V_{{n-1}}, band=middle): {H_cond:.4f} bits", flush=True)
    print(f"  Reduction: {H_marg - H_cond:.4f} bits = I(V_n; V_{{n-1}} | band)", flush=True)

    # ============= Verdict =============
    print(f"\n=== VERDICT ===", flush=True)
    print(f"  Pooled ρ(V_n, V_{{n-1}}) = {rho_VV:+.5f}", flush=True)
    print(f"  Pooled ρ(V_n, V_{{n-2}}) = {rho_VV2:+.5f}", flush=True)
    if abs(rho_VV) < 0.05:
        print(f"  → V values at consecutive r=21 visits are essentially uncorrelated (i.i.d. given band)", flush=True)
    else:
        print(f"  → V values have memory across visits", flush=True)


if __name__ == "__main__":
    main()
