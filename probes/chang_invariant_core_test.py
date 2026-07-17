"""
chang_invariant_core_test.py — test whether Chang's I_2 = {7, 27, 31, 59, 63}
mod 64 and our r=21 mod 32 are tracking the same structure.

Steps:
  (1) Mod-64 arithmetic table for all 32 odd residues
  (2) Per-residue v_2 deterministic-ness at mod 256
  (3) Per-residue orbit observables (σ, V, P(j)) at N=2^32 stratified by r mod 64
  (4) Joint Chang-singular-set × σ-band test

Skip (5): Hausdorff dim — requires Chang's specific framework.
"""
import sys
import time
import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2 * LOG2 - LOG3
MAX_VAL = np.int64(2**62)

CHANG_I_2 = {7, 27, 31, 59, 63}  # Chang/Quadrium I_2 invariant core mod 64
OUR_SINGULAR = {21, 53}  # r=21 mod 32 lifts to {21, 53} mod 64


def v2(n):
    if n == 0:
        return 0
    c = 0
    while n & 1 == 0:
        n >>= 1
        c += 1
    return c


def step1_mod64_arithmetic():
    """Compute v_2(3r+1), T(r) mod 64 for all 32 odd r mod 64."""
    print("=" * 80)
    print("Step 1: mod-64 arithmetic table for 32 odd residues")
    print("=" * 80)
    print(f"\nChang's I_2 = {sorted(CHANG_I_2)} mod 64")
    print(f"Our singular {{21, 53}} mod 64 (r=21 mod 32 lifts)")
    print(f"\n{'r':>3}  {'3r+1':>5}  {'v_2':>4}  {'T(r) mod 64':>12}  {'Chang I_2':>9}  {'Our sing':>8}  notes")

    rows = []
    for r in range(1, 64, 2):
        n3r1 = 3 * r + 1
        v = v2(n3r1)
        T_r = (n3r1 >> v) & 63
        in_chang = "✓" if r in CHANG_I_2 else " "
        in_ours = "✓" if r in OUR_SINGULAR else " "
        notes = ""
        if v >= 6:
            notes = "BOUNDARY (v ≥ 6, non-deterministic at higher mod)"
        elif n3r1 == (1 << v):
            notes = f"3r+1 = 2^{v} pure power"
        rows.append({
            'r': r, 'n3r1': n3r1, 'v_2': v, 'T_r_mod_64': T_r,
            'in_chang_I2': r in CHANG_I_2,
            'in_our_singular': r in OUR_SINGULAR,
            'notes': notes,
        })
        print(f"{r:>3}  {n3r1:>5}  {v:>4}  {T_r:>12}  {in_chang:>9}  {in_ours:>8}  {notes}")

    # Distribution of v_2 across residues
    print(f"\n# v_2(3r+1) distribution across 32 odd residues mod 64:")
    from collections import Counter
    v_counts = Counter(row['v_2'] for row in rows)
    for v in sorted(v_counts):
        rs = [row['r'] for row in rows if row['v_2'] == v]
        print(f"  v_2 = {v}: {len(rs)} residues: {rs}")

    # Chang's I_2 v_2 values
    print(f"\n# Chang's I_2 v_2 values:")
    for r in sorted(CHANG_I_2):
        row = next(row for row in rows if row['r'] == r)
        print(f"  r={r}: v_2(3·{r}+1) = v_2({row['n3r1']}) = {row['v_2']}")

    # Our singular set v_2 values
    print(f"\n# Our singular {{21, 53}} v_2 values:")
    for r in sorted(OUR_SINGULAR):
        row = next(row for row in rows if row['r'] == r)
        print(f"  r={r}: v_2(3·{r}+1) = v_2({row['n3r1']}) = {row['v_2']}")

    return rows


def step2_higher_mod_determinism(rows):
    """Test whether each residue mod 64's v_2(3m+1) is deterministic when
    we lift to mod 256. r=21 should be the only non-deterministic one."""
    print(f"\n{'=' * 80}")
    print("Step 2: v_2(3m+1) determinism at mod 256")
    print(f"{'=' * 80}")
    print(f"For each r mod 64, check whether v_2(3m+1) is constant across m ≡ r mod 64 with m mod 256 ∈ {{r, r+64, r+128, r+192}}")

    print(f"\n{'r mod 64':>9}  {'v_2(3·r+1)':>11}  {'v_2(3·(r+64)+1)':>16}  {'v_2(3·(r+128)+1)':>17}  {'v_2(3·(r+192)+1)':>17}  {'deterministic?':>15}")
    nondeterministic_residues = []
    for row in rows:
        r = row['r']
        v_values = [v2(3 * (r + k * 64) + 1) for k in range(4)]
        determ = "yes" if len(set(v_values)) == 1 else "NO"
        if determ == "NO":
            nondeterministic_residues.append(r)
        print(f"  {r:>3}    {v_values[0]:>11}    {v_values[1]:>16}    {v_values[2]:>17}    {v_values[3]:>17}    {determ:>15}")

    print(f"\n# Non-deterministic at mod 256: {nondeterministic_residues}")
    print(f"# Chang's I_2 ⊂ non-deterministic: {set(CHANG_I_2) <= set(nondeterministic_residues)}")
    print(f"# Our singular {{21, 53}} ⊂ non-deterministic: {set(OUR_SINGULAR) <= set(nondeterministic_residues)}")


@njit(parallel=True, cache=True)
def walk(starts, max_steps=400000):
    n = len(starts)
    sigma = np.zeros(n, dtype=np.int32)
    j_arr = np.full(n, -1, dtype=np.int8)
    v_sum = np.zeros(n, dtype=np.int64)
    for i in prange(n):
        m = starts[i]
        steps = 0
        v_total = 0
        j_attr = -1
        while m != 1 and steps < max_steps:
            if m > MAX_VAL // 3:
                break
            three_m = 3 * m + 1
            v = 0
            while (three_m & 1) == 0:
                three_m >>= 1
                v += 1
            steps += 1
            v_total += v
            if three_m == 1:
                if v % 2 == 0:
                    j_attr = v // 2
            m = three_m
        sigma[i] = steps
        j_arr[i] = j_attr
        v_sum[i] = v_total
    return sigma, j_arr, v_sum


def step3_per_residue_observables():
    """Walk orbits at N=2^32, stratify by starting residue mod 64,
    compute σ, V, P(j) per residue."""
    print(f"\n{'=' * 80}")
    print("Step 3: Per-residue orbit observables at N=2^32 (5M orbits)")
    print(f"{'=' * 80}")

    N = 1 << 32
    n_starts = 5_000_000
    rng = np.random.default_rng(42)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_starts, dtype=np.int64) + 1

    t0 = time.perf_counter()
    sigma_arr, j_arr, v_sum_arr = walk(starts)
    print(f"# Walk: {time.perf_counter()-t0:.1f}s, {n_starts:,} orbits")

    valid = sigma_arr > 0
    starts = starts[valid]
    sigma_arr = sigma_arr[valid]
    j_arr = j_arr[valid]
    v_sum_arr = v_sum_arr[valid]
    log_m = np.log(starts.astype(np.float64))
    v_avg = v_sum_arr.astype(np.float64) / np.maximum(sigma_arr, 1).astype(np.float64)

    # Stratify by starting residue mod 64
    res = (starts & 63).astype(np.int64)

    print(f"\n{'r':>3}  {'set':>10}  {'n_orbits':>10}  {'⟨σ|r⟩':>8}  {'⟨V|r⟩':>8}  {'P(j=2)':>8}  {'P(j=4)':>8}  {'P(j=5)':>8}")

    rows = []
    for r in range(1, 64, 2):
        mask = res == r
        n_r = int(mask.sum())
        if n_r < 100:
            continue
        sg = sigma_arr[mask]
        v = v_avg[mask]
        j = j_arr[mask]
        sigma_mean = float(sg.mean())
        v_mean = float(v.mean())
        P_j2 = float((j == 2).sum()) / n_r
        P_j4 = float((j == 4).sum()) / n_r
        P_j5 = float((j == 5).sum()) / n_r
        membership = []
        if r in CHANG_I_2:
            membership.append("CHANG")
        if r in OUR_SINGULAR:
            membership.append("OURS")
        ms_str = "|".join(membership) if membership else "generic"
        rows.append({
            'r': r, 'n_orbits': n_r, 'mean_sigma': sigma_mean,
            'mean_V': v_mean, 'P_j2': P_j2, 'P_j4': P_j4, 'P_j5': P_j5,
            'in_chang': r in CHANG_I_2, 'in_ours': r in OUR_SINGULAR,
        })
        print(f"  {r:>3}  {ms_str:>10}  {n_r:>10,}  {sigma_mean:>8.3f}  {v_mean:>8.4f}  {P_j2:>8.4f}  {P_j4:>8.4f}  {P_j5:>8.4f}")

    # Group statistics
    print(f"\n# Group means (3 groups: Chang's I_2, our singular {{21, 53}}, generic):")
    groups = {
        'Chang I_2 ({7,27,31,59,63})': CHANG_I_2,
        'Our singular ({21, 53})': OUR_SINGULAR,
        'Generic (other 25)': set(range(1, 64, 2)) - CHANG_I_2 - OUR_SINGULAR,
    }
    for label, residue_set in groups.items():
        rs = [row for row in rows if row['r'] in residue_set]
        if not rs:
            continue
        sigmas = np.array([r['mean_sigma'] for r in rs])
        Vs = np.array([r['mean_V'] for r in rs])
        Pj2 = np.array([r['P_j2'] for r in rs])
        Pj4 = np.array([r['P_j4'] for r in rs])
        Pj5 = np.array([r['P_j5'] for r in rs])
        print(f"  {label}:")
        print(f"     n_residues = {len(rs)}, total orbits = {sum(r['n_orbits'] for r in rs):,}")
        print(f"     ⟨σ⟩ mean over residues: {sigmas.mean():.3f} ± {sigmas.std(ddof=1) if len(sigmas)>1 else 0:.3f}")
        print(f"     ⟨V⟩ mean over residues: {Vs.mean():.4f} ± {Vs.std(ddof=1) if len(Vs)>1 else 0:.4f}")
        print(f"     P(j=2) mean: {Pj2.mean():.4f}, range [{Pj2.min():.4f}, {Pj2.max():.4f}]")
        print(f"     P(j=4) mean: {Pj4.mean():.4f}, range [{Pj4.min():.4f}, {Pj4.max():.4f}]")
        print(f"     P(j=5) mean: {Pj5.mean():.4f}, range [{Pj5.min():.4f}, {Pj5.max():.4f}]")

    # Comparison test
    print(f"\n# Statistical comparison: are Chang's I_2 and Our singular distinguishable from generic?")
    chang_sigmas = np.array([row['mean_sigma'] for row in rows if row['in_chang']])
    ours_sigmas = np.array([row['mean_sigma'] for row in rows if row['in_ours']])
    generic_sigmas = np.array([row['mean_sigma'] for row in rows
                                if not row['in_chang'] and not row['in_ours']])
    chang_Vs = np.array([row['mean_V'] for row in rows if row['in_chang']])
    ours_Vs = np.array([row['mean_V'] for row in rows if row['in_ours']])
    generic_Vs = np.array([row['mean_V'] for row in rows
                            if not row['in_chang'] and not row['in_ours']])

    print(f"  Chang I_2 ⟨σ⟩ vs Generic ⟨σ⟩:  {chang_sigmas.mean():.3f} vs {generic_sigmas.mean():.3f}, "
          f"diff = {chang_sigmas.mean() - generic_sigmas.mean():+.3f}")
    print(f"  Our singular ⟨σ⟩ vs Generic ⟨σ⟩:  {ours_sigmas.mean():.3f} vs {generic_sigmas.mean():.3f}, "
          f"diff = {ours_sigmas.mean() - generic_sigmas.mean():+.3f}")
    print(f"  Chang I_2 ⟨V⟩ vs Generic ⟨V⟩:  {chang_Vs.mean():.4f} vs {generic_Vs.mean():.4f}, "
          f"diff = {chang_Vs.mean() - generic_Vs.mean():+.4f}")
    print(f"  Our singular ⟨V⟩ vs Generic ⟨V⟩:  {ours_Vs.mean():.4f} vs {generic_Vs.mean():.4f}, "
          f"diff = {ours_Vs.mean() - generic_Vs.mean():+.4f}")

    # Save
    pl.DataFrame(rows).write_csv("experiments_output/chang_invariant_core_test.csv")
    print(f"\n[save] experiments_output/chang_invariant_core_test.csv")

    return rows


def main():
    arith_rows = step1_mod64_arithmetic()
    step2_higher_mod_determinism(arith_rows)
    obs_rows = step3_per_residue_observables()


if __name__ == "__main__":
    main()
