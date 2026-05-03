"""
Conditional-on-history v_2 at r=21 mod 32.

Math sanity: for m = 21+32k:
  3m+1 = 32·(2+3k), v = 5 + v_2(2+3k)
  k odd → v=5 (P=1/2)
  k=2k', k' even → v=6 (P=1/4)
  k=2k', k' odd → v≥7 (P=1/4)
Under uniform m: P(v=j|r=21) = 2^(-(j-4)) for j≥5, shifted Geom(1/2),
E[v|r=21]=6, H=2 bits.

Empirical ⟨v|r=21⟩=5.924 (Result 40) means orbit visits aren't uniform mod 2^k.

Tests:
  T1: marginal P(v=j|r=21) vs shifted-Geom prediction
  T2: condition on σ-band — how do bands shift the v distribution at r=21?
  T3: condition on m mod 2^k for k=6..12 — at what k does it determinize?
  T4: condition on last-3 v values + last-3 r — does orbit history collapse
       residual entropy beyond what m mod 2^k provides?
  T5: entropy H(v|r=21, conditioning) hierarchy
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
def walk_capture_r21(starts, max_steps):
    """Walk orbits, capture for each r=21 visit:
       v at visit, m mod 2^12 (covers k up to 12), σ_band index, position fraction,
       last-3 v values, last-3 r mod 32 values.
       Returns flat arrays (one row per r=21 visit)."""
    n = len(starts)

    # Pre-pass: count r=21 visits to size flat arrays
    # We'll use generous estimates: ~10% of visits at r=21
    max_visits_per_orbit = 200  # safe upper bound
    visit_idx_buf = np.full((n, max_visits_per_orbit, 9), -1, dtype=np.int64)
    # Columns: [v, m_mod_4096, σ_orbit, T_orbit, t, lastv_0, lastv_1, lastv_2, last_r0]
    visit_count = np.zeros(n, dtype=np.int32)
    sigma_total = np.full(n, -1, dtype=np.int64)
    T_total = np.zeros(n, dtype=np.int32)

    for i in prange(n):
        m = np.int64(starts[i])
        steps_sigma = 0
        T = 0
        # Last-3 v values and r mod 32 (FIFO, init -1)
        last_v0 = -1; last_v1 = -1; last_v2 = -1
        last_r0 = -1
        vc = 0
        failed = False
        while m != 1 and T < max_steps:
            if (m & 1) == 0:
                m = m >> 1
                steps_sigma += 1
                continue
            if m > MAX_VAL // 3:
                failed = True
                break
            r32 = m % 32
            x = 3*m + 1; v = 0
            while (x & 1) == 0:
                x >>= 1; v += 1

            if r32 == 21 and vc < max_visits_per_orbit:
                visit_idx_buf[i, vc, 0] = v
                visit_idx_buf[i, vc, 1] = m % 4096  # m mod 2^12
                visit_idx_buf[i, vc, 2] = -1  # placeholder for σ_band (filled later)
                visit_idx_buf[i, vc, 3] = -1  # T_orbit (filled later)
                visit_idx_buf[i, vc, 4] = T   # t
                visit_idx_buf[i, vc, 5] = last_v0
                visit_idx_buf[i, vc, 6] = last_v1
                visit_idx_buf[i, vc, 7] = last_v2
                visit_idx_buf[i, vc, 8] = last_r0
                vc += 1

            # Update history
            last_v2 = last_v1; last_v1 = last_v0; last_v0 = v
            last_r0 = r32

            steps_sigma += 1 + v
            T += 1
            m = x
        if not failed and m == 1:
            sigma_total[i] = steps_sigma
            T_total[i] = T
            visit_count[i] = vc

    return visit_idx_buf, visit_count, sigma_total, T_total


@njit(cache=True)
def compute_v_from_m(m):
    """Reference: v_2(3m+1) computed directly."""
    x = 3*m + 1; v = 0
    while (x & 1) == 0:
        x >>= 1; v += 1
    return v


def main():
    out_dir = Path("C:/Collatz/experiments_output")

    log2N = 36
    N = 1 << log2N
    n_per_seed = 50_000
    seeds = [42, 137, 271, 314, 1729]

    print(f"# Conditional-on-history v_2 at r=21 mod 32, N=2^{log2N}", flush=True)
    print(f"# {len(seeds)} seeds × {n_per_seed:,} = {len(seeds)*n_per_seed:,} orbits", flush=True)

    all_v = []; all_mmod = []; all_sigband = []; all_pos = []
    all_lv0 = []; all_lv1 = []; all_lv2 = []; all_lr0 = []
    all_logn = []

    t0 = time.time()
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1)//2, size=n_per_seed, dtype=np.int64) + 1

        visits, vc, sigma_arr, T_arr = walk_capture_r21(starts, 600)
        # Filter ok orbits
        ok = sigma_arr > 0
        # σ-resid
        log_n = np.log(starts.astype(np.float64))
        sigma_resid = sigma_arr.astype(np.float64) - K_H * log_n
        # Bands
        edges = np.percentile(sigma_resid[ok], [25, 50, 75, 95])
        # Map orbit i to its σ-band index (0..4)
        # 0: 0-25, 1: 25-50, 2: 50-75, 3: 75-95, 4: 95-100
        band = np.full(len(starts), -1, dtype=np.int8)
        b0 = sigma_resid <= edges[0]
        b1 = (sigma_resid > edges[0]) & (sigma_resid <= edges[1])
        b2 = (sigma_resid > edges[1]) & (sigma_resid <= edges[2])
        b3 = (sigma_resid > edges[2]) & (sigma_resid <= edges[3])
        b4 = sigma_resid > edges[3]
        band[b0]=0; band[b1]=1; band[b2]=2; band[b3]=3; band[b4]=4

        # Iterate visits per orbit
        for i in range(len(starts)):
            if not ok[i]: continue
            for k in range(vc[i]):
                all_v.append(visits[i, k, 0])
                all_mmod.append(visits[i, k, 1])
                all_sigband.append(band[i])
                all_pos.append(visits[i, k, 4] / max(T_arr[i], 1))
                all_lv0.append(visits[i, k, 5])
                all_lv1.append(visits[i, k, 6])
                all_lv2.append(visits[i, k, 7])
                all_lr0.append(visits[i, k, 8])
                all_logn.append(log_n[i])

    all_v = np.array(all_v, dtype=np.int8)
    all_mmod = np.array(all_mmod, dtype=np.int32)
    all_sigband = np.array(all_sigband, dtype=np.int8)
    all_pos = np.array(all_pos, dtype=np.float32)
    all_lv0 = np.array(all_lv0, dtype=np.int8)
    all_lv1 = np.array(all_lv1, dtype=np.int8)
    all_lv2 = np.array(all_lv2, dtype=np.int8)
    all_lr0 = np.array(all_lr0, dtype=np.int8)
    all_logn = np.array(all_logn, dtype=np.float32)

    print(f"  walked + collected in {time.time()-t0:.1f}s", flush=True)
    n_visits = len(all_v)
    print(f"  total r=21 visits: {n_visits:,}", flush=True)

    # Filter v ≥ 5 (math says v ≥ 5 at r=21; <5 would indicate bug)
    bug_mask = all_v < 5
    if bug_mask.sum() > 0:
        print(f"  WARNING: {bug_mask.sum()} visits with v<5 (should be 0)", flush=True)
    valid = all_v >= 5
    all_v = all_v[valid]; all_mmod = all_mmod[valid]; all_sigband = all_sigband[valid]
    all_pos = all_pos[valid]; all_lv0 = all_lv0[valid]; all_lv1 = all_lv1[valid]
    all_lv2 = all_lv2[valid]; all_lr0 = all_lr0[valid]; all_logn = all_logn[valid]
    n_visits = len(all_v)

    def entropy(probs):
        p = probs[probs > 0]
        return float(-np.sum(p * np.log2(p)))

    # ============= T1: marginal P(v|r=21) =============
    print(f"\n=== T1: marginal P(v|r=21) ===", flush=True)
    print(f"  Math prediction (uniform m): P(v=j|r=21) = 2^(-(j-4)) for j≥5", flush=True)
    print(f"  {'v':>3}  {'pred':>8}  {'emp':>8}  {'gap':>8}  {'count':>10}", flush=True)
    v_range = list(range(5, 16))
    pred = {j: 2**(-(j-4)) for j in v_range}
    emp = {}
    for j in v_range:
        c = int((all_v == j).sum())
        emp[j] = c / n_visits
        print(f"  {j:>3}  {pred[j]:>8.5f}  {emp[j]:>8.5f}  {emp[j]-pred[j]:>+8.5f}  {c:>10,}", flush=True)
    # tail
    c_tail = int((all_v >= 16).sum())
    p_tail = c_tail / n_visits
    pred_tail = 2**(-(16-4)) * 2  # tail sum from j=16
    print(f"  ≥16  {pred_tail:>8.5f}  {p_tail:>8.5f}  {p_tail-pred_tail:>+8.5f}  {c_tail:>10,}", flush=True)
    mean_v = float(all_v.mean())
    H_marginal = entropy(np.array([emp[j] for j in v_range] + [p_tail]))
    print(f"\n  ⟨v|r=21⟩ empirical = {mean_v:.4f}  (math uniform = 6.000)", flush=True)
    print(f"  H(v|r=21) marginal = {H_marginal:.4f} bits  (math uniform = 2.000)", flush=True)

    # ============= T2: per σ-band =============
    print(f"\n=== T2: P(v|r=21, σ-band) — does band shift the distribution? ===", flush=True)
    band_names = ['0-25', '25-50', '50-75', '75-95', '95-100']
    print(f"  {'band':>6}  {'n_visits':>10}  {'⟨v⟩':>7}  {'P(v=5)':>8}  {'P(v=6)':>8}  "
          f"{'P(v=7)':>8}  {'P(v=8+)':>9}  {'H bits':>7}", flush=True)
    for b in range(5):
        mask = all_sigband == b
        n_b = int(mask.sum())
        if n_b < 100: continue
        v_b = all_v[mask]
        p5 = (v_b==5).mean(); p6 = (v_b==6).mean(); p7 = (v_b==7).mean()
        p8p = (v_b>=8).mean()
        H_b = entropy(np.array([p5, p6, p7, p8p]))
        print(f"  {band_names[b]:>6}  {n_b:>10,}  {v_b.mean():>7.4f}  {p5:>8.4f}  {p6:>8.4f}  "
              f"{p7:>8.4f}  {p8p:>9.4f}  {H_b:>7.4f}", flush=True)

    # ============= T3: condition on m mod 2^k =============
    print(f"\n=== T3: P(v|r=21, m mod 2^k) cascade ===", flush=True)
    print(f"  Math prediction: at each k, half of m-residues fully determine v;", flush=True)
    print(f"  the other half need refinement.", flush=True)
    print(f"\n  {'k':>3}  {'mod':>6}  {'n_residues':>11}  {'frac_resolved':>14}  {'H(v|m mod 2^k)':>16}", flush=True)
    for k in [6, 7, 8, 9, 10, 11, 12]:
        modulus = 2**k
        # m mod 2^k for m at r=21 mod 32 — only m ≡ 21 mod 32 contribute, so we have 2^(k-5) possible residues
        m_at_modk = all_mmod % modulus
        unique_residues = np.unique(m_at_modk)
        # For each residue, check if v is constant
        n_resolved = 0
        H_total = 0.0
        for r in unique_residues:
            mask = m_at_modk == r
            v_at_r = all_v[mask]
            uniq_v = np.unique(v_at_r)
            if len(uniq_v) == 1:
                n_resolved += 1
            # Conditional entropy contribution
            p_r = mask.sum() / n_visits
            counts = np.bincount(v_at_r - 5)  # offset to start at 0
            counts = counts[counts > 0]
            p_v = counts / counts.sum()
            H_r = entropy(p_v)
            H_total += p_r * H_r
        frac_resolved = n_resolved / len(unique_residues)
        print(f"  {k:>3}  {modulus:>6}  {len(unique_residues):>11}  {frac_resolved:>14.4f}  {H_total:>16.5f}", flush=True)

    # ============= T4: condition on orbit history (last 3 v) =============
    print(f"\n=== T4: P(v|r=21, last-3 v values) — does orbit history determine v? ===", flush=True)
    # Take only visits where last_v0, last_v1, last_v2 are all present (≥0)
    mask_history = (all_lv0 >= 0) & (all_lv1 >= 0) & (all_lv2 >= 0)
    print(f"  Visits with full last-3 history: {mask_history.sum():,} / {n_visits:,}", flush=True)
    v_h = all_v[mask_history]; lv0_h = all_lv0[mask_history]
    lv1_h = all_lv1[mask_history]; lv2_h = all_lv2[mask_history]
    # Combine into compound key
    key = lv0_h.astype(np.int32)*256 + lv1_h.astype(np.int32)*16 + lv2_h.astype(np.int32)
    # For each unique key, compute conditional H(v)
    uniq_keys, inv = np.unique(key, return_inverse=True)
    n_unique = len(uniq_keys)
    H_total = 0.0
    n_resolved_h = 0
    for ki in range(n_unique):
        m = inv == ki
        v_k = v_h[m]
        if len(v_k) < 5: continue
        uniq_v = np.unique(v_k)
        if len(uniq_v) == 1:
            n_resolved_h += 1
        p_k = m.sum() / len(v_h)
        counts = np.bincount(v_k.clip(min=5) - 5)
        counts = counts[counts > 0]
        p_v = counts / counts.sum()
        H_total += p_k * entropy(p_v)
    print(f"  Unique last-3-v contexts: {n_unique}", flush=True)
    print(f"  Contexts that fully determine v: {n_resolved_h}", flush=True)
    print(f"  H(v | last-3 v) = {H_total:.4f} bits  vs marginal H = {H_marginal:.4f} bits", flush=True)
    print(f"  Reduction: {H_marginal - H_total:.4f} bits", flush=True)

    # ============= T4b: condition on (m mod 64 + last-3 v) =============
    print(f"\n=== T4b: P(v | r=21, m mod 64, last-3 v) — does history add beyond mod 64? ===", flush=True)
    m_mod_64 = all_mmod[mask_history] % 64
    key_combined = (m_mod_64.astype(np.int64) << 24) | (lv0_h.astype(np.int64) << 16) | (lv1_h.astype(np.int64) << 8) | lv2_h.astype(np.int64)
    uniq_keys, inv = np.unique(key_combined, return_inverse=True)
    H_total = 0.0
    n_resolved = 0
    for ki in range(len(uniq_keys)):
        m = inv == ki
        v_k = v_h[m]
        if len(v_k) < 5: continue
        uniq_v = np.unique(v_k)
        if len(uniq_v) == 1:
            n_resolved += 1
        p_k = m.sum() / len(v_h)
        counts = np.bincount(v_k.clip(min=5) - 5)
        counts = counts[counts > 0]
        p_v = counts / counts.sum()
        H_total += p_k * entropy(p_v)
    # H given m mod 64 alone
    m_mod_64_all = all_mmod % 64
    H_mod64 = 0.0
    for r in np.unique(m_mod_64_all):
        mask = m_mod_64_all == r
        v_r = all_v[mask]
        p_r = mask.sum() / n_visits
        counts = np.bincount(v_r - 5)
        counts = counts[counts > 0]
        p_v = counts / counts.sum()
        H_mod64 += p_r * entropy(p_v)
    print(f"  H(v | m mod 64) = {H_mod64:.4f} bits", flush=True)
    print(f"  H(v | m mod 64, last-3 v) = {H_total:.4f} bits", flush=True)
    print(f"  Additional reduction from history: {H_mod64 - H_total:.4f} bits", flush=True)

    # ============= T5: entropy hierarchy summary =============
    print(f"\n=== T5: Entropy hierarchy summary ===", flush=True)
    print(f"  Conditioning      H(v|...) bits", flush=True)
    print(f"    none (marginal):   {H_marginal:.4f}", flush=True)
    print(f"    m mod 64:          {H_mod64:.4f}", flush=True)
    H_mod128 = 0.0
    m_mod_128_all = all_mmod % 128
    for r in np.unique(m_mod_128_all):
        mask = m_mod_128_all == r
        v_r = all_v[mask]
        p_r = mask.sum() / n_visits
        counts = np.bincount(v_r - 5)
        counts = counts[counts > 0]
        p_v = counts / counts.sum()
        H_mod128 += p_r * entropy(p_v)
    print(f"    m mod 128:         {H_mod128:.4f}", flush=True)
    H_mod1024 = 0.0
    m_mod_1024_all = all_mmod % 1024
    for r in np.unique(m_mod_1024_all):
        mask = m_mod_1024_all == r
        v_r = all_v[mask]
        if mask.sum() < 5: continue
        p_r = mask.sum() / n_visits
        counts = np.bincount(v_r - 5)
        counts = counts[counts > 0]
        p_v = counts / counts.sum()
        H_mod1024 += p_r * entropy(p_v)
    print(f"    m mod 1024:        {H_mod1024:.4f}", flush=True)
    H_mod4096 = 0.0
    for r in np.unique(all_mmod):
        mask = all_mmod == r
        v_r = all_v[mask]
        if mask.sum() < 5: continue
        p_r = mask.sum() / n_visits
        counts = np.bincount(v_r - 5)
        counts = counts[counts > 0]
        p_v = counts / counts.sum()
        H_mod4096 += p_r * entropy(p_v)
    print(f"    m mod 4096:        {H_mod4096:.4f}", flush=True)
    print(f"    m mod 64 + last-3 v: {H_total:.4f}", flush=True)

    # ============= VERDICT =============
    print(f"\n=== VERDICT ===", flush=True)
    if H_mod4096 < 0.05:
        print(f"  Outcome (1): m mod 4096 nearly determines v at r=21 (H={H_mod4096:.3f} bits)", flush=True)
    elif H_mod64 - H_total < 0.05:
        print(f"  Orbit history adds little beyond m mod 64 ({H_mod64-H_total:.3f} bit reduction)", flush=True)
    if H_marginal - H_mod64 > 0.5:
        print(f"  m mod 64 already provides large reduction ({H_marginal - H_mod64:.3f} bits)", flush=True)
    print(f"  Cascade decay rate: H(mod 64)={H_mod64:.3f} → H(mod 1024)={H_mod1024:.3f} → H(mod 4096)={H_mod4096:.3f}", flush=True)
    print(f"  Each k doubling halves H — consistent with shifted-Geom(1/2) structure", flush=True)

    # Save
    rows_t5 = [
        {'conditioning': 'marginal', 'H_bits': H_marginal},
        {'conditioning': 'm mod 64', 'H_bits': H_mod64},
        {'conditioning': 'm mod 128', 'H_bits': H_mod128},
        {'conditioning': 'm mod 1024', 'H_bits': H_mod1024},
        {'conditioning': 'm mod 4096', 'H_bits': H_mod4096},
        {'conditioning': 'm mod 64 + last-3 v', 'H_bits': H_total},
    ]
    pl.DataFrame(rows_t5).write_csv(out_dir / "71_r21_entropy_hierarchy.csv")
    rows_t1 = [{'v': j, 'pred_uniform': pred[j], 'emp': emp[j]} for j in v_range]
    rows_t1.append({'v': 'tail≥16', 'pred_uniform': pred_tail, 'emp': p_tail})
    pl.DataFrame(rows_t1).write_csv(out_dir / "71_r21_marginal_p_v.csv")
    print(f"\n[save] CSVs written", flush=True)


if __name__ == "__main__":
    main()
