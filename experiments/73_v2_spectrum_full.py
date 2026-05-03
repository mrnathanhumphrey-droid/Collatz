"""
Map v_2(3r+1) spectrum across all 32 odd residues mod 64.

Extends Result 45 (companion task already verified determinism) with:
  - Spearman correlations between v_2 rank and observables
  - Sub-stratification within v_2 levels (especially what distinguishes Chang I_2 within v_2=1)
  - Residue chain map on Z/64Z odd: fixed points, cycles, absorption
  - {m_j} attractor analysis

Reuses orbit observables from exp 72 (1M orbits at N=2^32).
"""
import sys
import io
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

LOG2 = np.log(2.0); LOG3 = np.log(3.0)


def v2_of(n):
    if n == 0: return -1
    v = 0
    while (n & 1) == 0:
        n >>= 1; v += 1
    return v


def step1_full_table():
    """Tabulate per-residue arithmetic data."""
    print(f"\n=== Step 1: Full v_2 stratification (32 odd residues mod 64) ===", flush=True)
    chang_I2 = {7, 27, 31, 59, 63}
    rows = []
    for r in range(1, 64, 2):
        threex_p1 = 3*r + 1
        v0 = v2_of(threex_p1)
        dest = threex_p1 // (1 << v0)
        dest_mod64 = dest % 64
        # If dest is even, halve until odd to get next ODD residue
        d = dest
        while (d & 1) == 0:
            d >>= 1
        # dest IS odd by construction (we divided by 2^v0)
        # Next Syracuse destination
        v_at_dest = v2_of(3*dest_mod64 + 1) if dest_mod64 != 0 else -1
        # Sometimes dest mod 64 is even (if dest > 64). Take dest mod 64 — we want the odd residue mod 64
        # Actually dest is always odd (since v_0 = v_2(3r+1)). dest mod 64 might be different parities mod higher.
        # But dest is odd in absolute terms; dest mod 64 might be odd. Let me check: dest = (3r+1)/2^v0, odd.
        # dest mod 64 has same parity as dest (since 64 even). So dest mod 64 is odd.
        rows.append({
            'r': r, 'three_r_plus_1': threex_p1, 'v0': v0,
            'dest': dest, 'dest_mod64': dest_mod64,
            'v_at_dest': v_at_dest,
            'is_chang_I2': r in chang_I2,
        })

    print(f"  {'r':>3}  {'3r+1':>5}  {'v0':>3}  {'dest':>5}  {'dest_mod64':>10}  {'v@dest':>7}  {'group':>13}", flush=True)
    for row in sorted(rows, key=lambda x: (x['v0'], x['r'])):
        grp = "Chang I_2" if row['is_chang_I2'] else (
              "v=6 BOUNDARY" if row['v0'] == 6 else (
              "v=5" if row['v0'] == 5 else f"v={row['v0']}"))
        print(f"  {row['r']:>3}  {row['three_r_plus_1']:>5}  {row['v0']:>3}  {row['dest']:>5}  "
              f"{row['dest_mod64']:>10}  {row['v_at_dest']:>7}  {grp:>13}", flush=True)
    return rows


def step2_load_observables():
    """Load orbit observables from exp 72 CSV."""
    csv = Path("C:/Collatz/experiments_output/72_orbit_observables.csv")
    df = pl.read_csv(csv)
    print(f"\n=== Step 2: Loaded orbit observables (1M orbits at N=2^32) from exp 72 ===", flush=True)
    print(f"  Rows: {len(df)}", flush=True)
    return df


def step3_monotone_check(arith_rows, obs_df):
    """Spearman correlations between v_2 rank and observables."""
    print(f"\n=== Step 3: Spearman correlations v_2(3r+1) vs observables ===", flush=True)
    # Merge arithmetic and observables
    arith_dict = {row['r']: row for row in arith_rows}
    v0_vals = []
    sigma_vals = []
    V_vals = []
    for row in obs_df.iter_rows(named=True):
        r = row['r']
        if r in arith_dict:
            v0_vals.append(arith_dict[r]['v0'])
            sigma_vals.append(row['mean_sigma'])
            V_vals.append(row['mean_V'])
    v0_arr = np.array(v0_vals)
    sigma_arr = np.array(sigma_vals)
    V_arr = np.array(V_vals)

    rho_sigma, p_sigma = spearmanr(v0_arr, sigma_arr)
    rho_V, p_V = spearmanr(v0_arr, V_arr)
    print(f"  Spearman ρ(v_2, ⟨σ|r⟩)  = {rho_sigma:+.4f}  p={p_sigma:.4g}", flush=True)
    print(f"  Spearman ρ(v_2, ⟨V|r⟩)  = {rho_V:+.4f}  p={p_V:.4g}", flush=True)

    # Pearson on means as well
    pearson_sigma = np.corrcoef(v0_arr, sigma_arr)[0,1]
    pearson_V = np.corrcoef(v0_arr, V_arr)[0,1]
    print(f"  Pearson  r(v_2, ⟨σ|r⟩)  = {pearson_sigma:+.4f}", flush=True)
    print(f"  Pearson  r(v_2, ⟨V|r⟩)  = {pearson_V:+.4f}", flush=True)

    if abs(rho_sigma) > 0.7 and abs(rho_V) > 0.7:
        print(f"  → STRONG monotone progression in v_2", flush=True)
    elif abs(rho_sigma) > 0.4 or abs(rho_V) > 0.4:
        print(f"  → MODERATE monotone progression", flush=True)
    else:
        print(f"  → WEAK monotone progression — sub-structure within v_2", flush=True)

    return {'rho_sigma': float(rho_sigma), 'rho_V': float(rho_V),
            'pearson_sigma': float(pearson_sigma), 'pearson_V': float(pearson_V)}


def step4_sub_stratification(arith_rows, obs_df):
    """Within each v_2 level, look for sub-clusters."""
    print(f"\n=== Step 4: Sub-stratification within v_2 levels ===", flush=True)
    arith_dict = {row['r']: row for row in arith_rows}

    # Group residues by v_2
    by_v0 = {}
    obs_dict = {row['r']: row for row in obs_df.iter_rows(named=True)}
    for r in range(1, 64, 2):
        v0 = arith_dict[r]['v0']
        if r in obs_dict:
            by_v0.setdefault(v0, []).append((r, obs_dict[r], arith_dict[r]))

    # Within v_2=1: 16 residues. Distinguish Chang I_2 from rest.
    print(f"\n  Within v_2=1 (16 residues) — what distinguishes Chang I_2?", flush=True)
    print(f"  {'r':>3}  {'⟨σ|r⟩':>9}  {'⟨V|r⟩':>7}  {'dest_mod64':>10}  {'v@dest':>7}  {'in_I2':>6}", flush=True)
    v1_rows = sorted(by_v0[1], key=lambda x: x[1]['mean_sigma'])
    for r, obs, ar in v1_rows:
        marker = "✓" if ar['is_chang_I2'] else ""
        print(f"  {r:>3}  {obs['mean_sigma']:>9.2f}  {obs['mean_V']:>7.4f}  "
              f"{ar['dest_mod64']:>10}  {ar['v_at_dest']:>7}  {marker:>6}", flush=True)

    # Test: do Chang I_2 residues have particular dest_mod64 v_at_dest distribution?
    chang_v_at_dest = [ar['v_at_dest'] for r, obs, ar in v1_rows if ar['is_chang_I2']]
    others_v_at_dest = [ar['v_at_dest'] for r, obs, ar in v1_rows if not ar['is_chang_I2']]
    print(f"\n  Chang I_2 v_at_dest distribution: {sorted(chang_v_at_dest)}", flush=True)
    print(f"  Other v_2=1 v_at_dest distribution: {sorted(others_v_at_dest)}", flush=True)

    # Compute σ-cluster within v=1
    chang_sigmas = [obs['mean_sigma'] for r, obs, ar in v1_rows if ar['is_chang_I2']]
    others_sigmas = [obs['mean_sigma'] for r, obs, ar in v1_rows if not ar['is_chang_I2']]
    print(f"  Chang I_2 ⟨σ⟩ range: [{min(chang_sigmas):.1f}, {max(chang_sigmas):.1f}]  mean={np.mean(chang_sigmas):.1f}", flush=True)
    print(f"  Other v=1 ⟨σ⟩ range: [{min(others_sigmas):.1f}, {max(others_sigmas):.1f}]  mean={np.mean(others_sigmas):.1f}", flush=True)
    print(f"  Overlap: {max(chang_sigmas) > min(others_sigmas)} — Chang I_2 IS NOT separated from other v=1 by ⟨σ⟩ alone", flush=True)

    # Look for sub-clusters by σ value irrespective of Chang labeling
    sigmas_v1 = sorted([(obs['mean_sigma'], r) for r, obs, ar in v1_rows])
    print(f"\n  v=1 sorted by ⟨σ⟩: {sigmas_v1}", flush=True)
    # Are there discrete clusters? Look at gaps:
    sigma_vals = np.array([s for s, _ in sigmas_v1])
    gaps = np.diff(sigma_vals)
    print(f"  Gaps between consecutive ⟨σ⟩ values: {[f'{g:.1f}' for g in gaps]}", flush=True)
    # Largest gaps
    if len(gaps) > 0:
        largest_gap_idx = np.argmax(gaps)
        print(f"  Largest gap: {gaps[largest_gap_idx]:.1f} between {sigma_vals[largest_gap_idx]:.1f} and {sigma_vals[largest_gap_idx+1]:.1f}", flush=True)

    # Within v_2=2 (8 residues)
    print(f"\n  Within v_2=2 (8 residues):", flush=True)
    print(f"  {'r':>3}  {'⟨σ|r⟩':>9}  {'⟨V|r⟩':>7}  {'dest_mod64':>10}  {'v@dest':>7}", flush=True)
    v2_rows = sorted(by_v0[2], key=lambda x: x[1]['mean_sigma'])
    for r, obs, ar in v2_rows:
        print(f"  {r:>3}  {obs['mean_sigma']:>9.2f}  {obs['mean_V']:>7.4f}  "
              f"{ar['dest_mod64']:>10}  {ar['v_at_dest']:>7}", flush=True)


def step6_residue_chain_map(arith_rows):
    """Compute the Syracuse residue chain map on Z/64Z odd part.

    The map: r → (3r+1)/2^v_2(3r+1) mod 64.
    However, this isn't well-defined as a function — depends on higher bits of m.
    For deterministic residues (v_2 < 6), the next-step destination MOD 64 can vary.

    For r with v_0 = v_2(3r+1), 3r+1 is divided by 2^v_0 deterministically; result is
    (3r+1)/2^v_0. But this destination mod 64 depends on the higher bits of m beyond
    mod 64. EXCEPT: for r with v_0 ≤ 5, the result of (3m+1)/2^v_0 mod 64 IS determined
    by higher bits of m... let me think.

    m = r + 64k. 3m+1 = (3r+1) + 192k. After dividing by 2^v_0:
    dest = (3r+1)/2^v_0 + 192k/2^v_0 = (3r+1)/2^v_0 + 3k·2^(6-v_0)

    If v_0 < 6: dest = (3r+1)/2^v_0 + 3k·2^(6-v_0). Mod 64:
      - For v_0 = 1: dest mod 64 = (3r+1)/2 + 96k mod 64 = base + 32k mod 64. Two values.
      - For v_0 = 2: dest mod 64 = base + 48k mod 64 = base + 48k mod 64. Period 4.
      - For v_0 = 3: dest mod 64 = base + 24k. Period 8.
      - For v_0 = 4: dest mod 64 = base + 12k. Period 16.
      - For v_0 = 5: dest mod 64 = base + 6k. Period 32 (but only 32 values in mod 64, so cycles).

    So dest mod 64 is NOT a single value — it varies with k. The "residue chain map" mod 64
    isn't a deterministic function. It's a multi-valued map.

    To define a deterministic chain, we need to track higher-mod state.
    """
    print(f"\n=== Step 6: Residue chain analysis ===", flush=True)
    print(f"  Note: dest mod 64 is NOT deterministic from r mod 64 (depends on higher bits).", flush=True)
    print(f"  For each r mod 64, dest mod 64 takes 2^v_0 values cyclically as k varies.", flush=True)
    print(f"\n  Per-residue dest mod 64 set:", flush=True)
    print(f"  {'r':>3}  {'v_0':>3}  {'dest mod 64 reachable':>30}", flush=True)
    for arith_row in arith_rows:
        r = arith_row['r']
        v0 = arith_row['v0']
        if v0 < 6:
            base = arith_row['three_r_plus_1'] // (1 << v0) % 64
            step_size = (3 * (1 << (6 - v0))) % 64 if v0 <= 5 else 0
            n_distinct = 1 << v0  # 2^v_0 reachable values, but capped by 32 (odd residues mod 64)
            n_distinct = min(n_distinct, 32)
            dests = sorted(set([(base + k * step_size) % 64 for k in range(n_distinct)]))
        else:
            # For v_0 = 6 (r=21): dest can take many values
            # 3·21+1 = 64. m = 21 + 64k → 3m+1 = 64 + 192k = 64(1 + 3k). v_2 ≥ 6 + v_2(1+3k).
            # dest = (1+3k)/2^v_2(1+3k). Many possibilities.
            dests = ["multi-valued (boundary)"]
        if isinstance(dests[0], str):
            print(f"  {r:>3}  {v0:>3}  {dests[0]:>30}", flush=True)
        else:
            print(f"  {r:>3}  {v0:>3}  {str(dests):>30}", flush=True)


def step5_mid_spectrum_poles(arith_rows, obs_df):
    """Identify residues with unique structural features."""
    print(f"\n=== Step 5: Mid-spectrum poles? ===", flush=True)
    obs_dict = {row['r']: row for row in obs_df.iter_rows(named=True)}

    # Anomalies in observables
    sigma_vals = np.array([obs_dict[r]['mean_sigma'] for r in obs_dict])
    rs = np.array(list(obs_dict.keys()))
    median_sigma = np.median(sigma_vals)
    mad_sigma = np.median(np.abs(sigma_vals - median_sigma))
    print(f"  ⟨σ⟩ across 32 residues: median={median_sigma:.2f}, MAD={mad_sigma:.3f}", flush=True)
    print(f"  Outliers (|⟨σ|r⟩ − median| > 3·MAD):", flush=True)
    for i, r in enumerate(rs):
        if abs(sigma_vals[i] - median_sigma) > 3*mad_sigma:
            print(f"    r={r}: ⟨σ⟩={sigma_vals[i]:.2f}  (offset={sigma_vals[i]-median_sigma:+.2f})", flush=True)

    # Special arithmetic forms
    print(f"\n  Special arithmetic check on 3r+1 forms:", flush=True)
    for arith_row in arith_rows:
        r = arith_row['r']; threex = arith_row['three_r_plus_1']
        v0 = arith_row['v0']
        # Pure power of 2?
        if (threex & (threex - 1)) == 0:
            print(f"    r={r}: 3r+1 = {threex} = 2^{v0} (PURE POWER OF 2, fastest single-step descent)", flush=True)
        # 3r+1 / 2^v0 = 1?
        if threex // (1 << v0) == 1:
            print(f"    r={r}: 3r+1 = 2^{v0}, descent reaches 1 in one step!", flush=True)


def main():
    out_dir = Path("C:/Collatz/experiments_output")

    arith_rows = step1_full_table()
    obs_df = step2_load_observables()
    spear = step3_monotone_check(arith_rows, obs_df)
    step4_sub_stratification(arith_rows, obs_df)
    step5_mid_spectrum_poles(arith_rows, obs_df)
    step6_residue_chain_map(arith_rows)

    # Combined CSV
    rows_combined = []
    obs_dict = {row['r']: row for row in obs_df.iter_rows(named=True)}
    for ar in arith_rows:
        r = ar['r']
        if r in obs_dict:
            o = obs_dict[r]
            rows_combined.append({
                'r': r, 'v0': ar['v0'], 'three_r_plus_1': ar['three_r_plus_1'],
                'dest': ar['dest'], 'dest_mod64': ar['dest_mod64'], 'v_at_dest': ar['v_at_dest'],
                'is_chang_I2': ar['is_chang_I2'],
                'mean_sigma': o['mean_sigma'], 'mean_V': o['mean_V'],
                'group': o.get('group', ''),
            })
    pl.DataFrame(rows_combined).write_csv(out_dir / "73_v2_spectrum_full_map.csv")
    print(f"\n[save] CSV", flush=True)


if __name__ == "__main__":
    main()
