"""
Arithmetic-deterministic re-derivation of |mu_hat(a/3^k)|^2 from empirical
trajectory measure mass profile (no Geom(1/2) heuristic).

Steps:
  1. v-lookup tables: v(r mod 2^k) for r odd, k in {4..12}
  2. Walk many orbits at N=2^32, record m_t visits
  3. Compute mod-2^k mass profile + verify v-distribution matches empirical
  4. Compute mod-3 fractions (a, b, c) from mass profile
  5. Direct Fourier: mu_hat(a/3^k) = (1/Z) sum_visits exp(2*pi*i*a*m/3^k) for k=1..7
  6. Compare to R66 analytical (Geom-based Markov chain stationary)
  7. Verify decay law: <|mu_hat|^2>_a * 3^(k-1) -> S_infinity
  8. Verdict on whether arithmetic-deterministic gives exact predictions

Outputs:
  experiments_output/86_mod_2k_lookup_table.csv
  experiments_output/86_mod_2k_mass_profile.csv
  experiments_output/86_arithmetic_vs_geom_comparison.csv
  experiments_output/86_S_infinity_arithmetic.csv
  experiments_output/86_arithmetic_deterministic_log.txt
"""
import sys
import io
import time
import math
import cmath
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

OUT = Path("C:/Collatz")
EXP_OUT = Path("C:/Collatz/experiments_output")

MAX_VAL = np.int64(2**62)

results_log = []
def log(s):
    print(s, flush=True)
    results_log.append(s)


# ============================================================
# Step 1: v-lookup tables — v(r mod 2^k) for odd r
# ============================================================

def v_for_residue(r, k_max):
    """v_2(3r+1) capped at k_max (return -1 if undetermined at this k)."""
    if r % 2 == 0:
        return 0  # not odd
    val = 3 * r + 1
    v = 0
    while val % 2 == 0 and v < k_max:
        val //= 2
        v += 1
    if v == k_max:
        return -1  # undetermined at this k
    return v


def build_v_lookup(k):
    """Return dict: r -> v(r) for odd r in Z/2^kZ. v=-1 means undetermined."""
    table = {}
    for r in range(1, 2**k, 2):
        table[r] = v_for_residue(r, k)
    return table


# ============================================================
# Step 2: Walker — record m_t visits at N=2^32
# Accumulate Fourier sums per (k, a) on the fly to avoid storing all m_t
# ============================================================

# We compute mu_hat(a/3^k) for k=1..K_MAX, all primitive a coprime to 3.
# State: complex accumulator per (k, a). Update: visit m -> add exp(2pi i a m / 3^k).

@njit(parallel=True, cache=True)
def walk_fourier_accum(starts, max_T, n_chunks, max_k_3, k_max_2):
    """
    Walk orbits, accumulate Fourier coefficients and mass profile.

    Returns:
      fourier_re, fourier_im : (n_chunks, sum_{k=1..max_k_3} 3^k) complex sums per (k, a)
      mod_2k_counts : (n_chunks, sum_{k=4..k_max_2 step 2} 2^(k-1)) for odd residues
      total_visits : (n_chunks,) total step count
      mod_3_counts : (n_chunks, 3) m mod 3 visit counts
    """
    n = len(starts)
    chunk_size = (n + n_chunks - 1) // n_chunks

    # For each k=1..max_k_3, store all 3^k complex accumulators.
    # offsets_3[k-1] = starting index of k's slots; k uses 3^k slots.
    offsets_3 = np.zeros(max_k_3 + 1, dtype=np.int64)
    for k in range(1, max_k_3 + 1):
        offsets_3[k] = offsets_3[k-1] + 3**k  # k uses 3^k slots
    total_3 = offsets_3[max_k_3]
    fourier_re = np.zeros((n_chunks, total_3), dtype=np.float64)
    fourier_im = np.zeros((n_chunks, total_3), dtype=np.float64)

    # mod-2^k counts for odd residues at k in {4, 6, 8, 10, 12}
    # We'll just store full arrays per k, indexed compactly
    # Use single flat array: offsets_2[k_idx] for each k value
    k_2_values = np.array([4, 6, 8, 10, 12], dtype=np.int64)
    offsets_2 = np.zeros(len(k_2_values) + 1, dtype=np.int64)
    for i in range(len(k_2_values)):
        offsets_2[i+1] = offsets_2[i] + (1 << k_2_values[i])  # full 2^k size, only odd entries used
    total_2 = offsets_2[-1]
    mod_2k_counts = np.zeros((n_chunks, total_2), dtype=np.int64)

    mod_3_counts = np.zeros((n_chunks, 3), dtype=np.int64)
    total_visits = np.zeros(n_chunks, dtype=np.int64)

    for chunk in prange(n_chunks):
        i_lo = chunk * chunk_size
        i_hi = min((chunk + 1) * chunk_size, n)
        for i in range(i_lo, i_hi):
            m = np.int64(starts[i])
            T = 0
            while (m & 1) == 0 and m > 1:
                m >>= 1
            if m == 1:
                continue

            failed = False
            while m != 1 and T < max_T:
                if m > MAX_VAL // 3:
                    failed = True; break

                # Record this m_t visit
                total_visits[chunk] += 1
                m_mod3 = m % 3
                mod_3_counts[chunk, m_mod3] += 1

                # Fourier accumulators
                m_f = np.float64(m)
                for k in range(1, max_k_3 + 1):
                    qk = 3**k
                    r = m % qk
                    base_phase = 2.0 * np.pi * np.float64(r) / np.float64(qk)
                    off_k = offsets_3[k - 1]
                    for a in range(qk):
                        # Skip a=0 trivially (constant) — but include for completeness
                        ang = base_phase * np.float64(a)
                        fourier_re[chunk, off_k + a] += np.cos(ang)
                        fourier_im[chunk, off_k + a] += np.sin(ang)

                # mod-2^k count
                for kk in range(len(k_2_values)):
                    k_2 = k_2_values[kk]
                    rr = m & ((1 << k_2) - 1)
                    mod_2k_counts[chunk, offsets_2[kk] + rr] += 1

                # Step
                x = 3 * m + 1
                while (x & 1) == 0:
                    x >>= 1
                T += 1
                m = x
    return fourier_re, fourier_im, mod_2k_counts, mod_3_counts, total_visits, offsets_3, k_2_values, offsets_2


# ============================================================
# Step 6: Load R66 analytical (Geom-based Markov chain) for comparison
# ============================================================

def load_R66_analytical():
    df = pl.read_csv(OUT / "experiments_output" / "decay_law_derivation.csv")
    # cols: k, a, analytical_mu_hat_sq, empirical_mu_hat_sq
    return df


# ============================================================
# Main
# ============================================================

def main():
    log("=" * 80)
    log("ARITHMETIC-DETERMINISTIC RE-DERIVATION OF |mu_hat(a/3^k)|^2")
    log("=" * 80)
    log(f"\n  Approach: walk orbits at N=2^32, accumulate complex Fourier sums")
    log(f"            per (k, a) directly from m_t visits — no Geom(1/2) heuristic.\n")

    # ============ Step 1: v-lookup tables ============
    log("=" * 80)
    log("STEP 1: v(r mod 2^k) lookup tables for k = 4, 6, 8, 10, 12")
    log("=" * 80)

    rows_lookup = []
    log(f"\n  {'k':>3}  {'#odd':>6}  {'#determined':>11}  {'#undetermined':>13}  {'frac_det':>9}")
    for k in [4, 6, 8, 10, 12]:
        table = build_v_lookup(k)
        n_odd = len(table)
        n_det = sum(1 for v in table.values() if v >= 0)
        n_und = n_odd - n_det
        frac = n_det / n_odd
        log(f"  {k:>3}  {n_odd:>6}  {n_det:>11}  {n_und:>13}  {frac:>9.4f}")
        for r, v in table.items():
            rows_lookup.append({'k': k, 'r': r, 'v': v})

    pl.DataFrame(rows_lookup).write_csv(EXP_OUT / "86_mod_2k_lookup_table.csv")
    log(f"  [save] 86_mod_2k_lookup_table.csv ({len(rows_lookup)} rows)")

    # ============ Step 2: Walk orbits, accumulate Fourier + mass profile ============
    log("\n" + "=" * 80)
    log("STEP 2: Walk orbits, accumulate Fourier sums + mod-2^k mass profile")
    log("=" * 80)

    log2N = 32; N = 1 << log2N
    n_orbits_per_seed = 100_000  # smaller for memory budget
    seeds = [42, 137, 271]
    max_k_3 = 5  # 3^5 = 243 — total Fourier size across k=1..5: 3+9+27+81+243 = 363

    log(f"\n  Walking {len(seeds) * n_orbits_per_seed:,} orbits at N=2^{log2N}")
    log(f"  Fourier accumulators: k = 1..{max_k_3} (3^k coefficients each)")

    # We'll do seeds sequentially and accumulate to keep memory manageable
    n_chunks = 8

    fourier_re_total = None
    fourier_im_total = None
    mod_2k_total = None
    mod_3_total = np.zeros(3, dtype=np.int64)
    total_visits = 0
    offsets_3_g = None
    k_2_values_g = None
    offsets_2_g = None

    t0 = time.time()
    for seed_i, seed in enumerate(seeds):
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1) // 2, size=n_orbits_per_seed, dtype=np.int64) + 1
        f_re, f_im, m2k, m3, tv, off_3, k_2, off_2 = walk_fourier_accum(
            starts, 600, n_chunks, max_k_3, 12)
        offsets_3_g = off_3
        k_2_values_g = k_2
        offsets_2_g = off_2
        if fourier_re_total is None:
            fourier_re_total = f_re.sum(axis=0)
            fourier_im_total = f_im.sum(axis=0)
            mod_2k_total = m2k.sum(axis=0)
        else:
            fourier_re_total += f_re.sum(axis=0)
            fourier_im_total += f_im.sum(axis=0)
            mod_2k_total += m2k.sum(axis=0)
        mod_3_total += m3.sum(axis=0)
        total_visits += int(tv.sum())
        log(f"  seed {seed_i+1}/{len(seeds)}: visits = {int(tv.sum()):,}, elapsed {time.time()-t0:.1f}s")

    log(f"\n  Total visits: {total_visits:,}")
    log(f"  Walk time: {time.time()-t0:.1f}s")

    # ============ Step 3: v-distribution from mass profile ============
    log("\n" + "=" * 80)
    log("STEP 3: v-distribution from mod-2^k mass profile + verify against empirical")
    log("=" * 80)

    # Use mod-2^12 profile to get high-fidelity v-distribution (covers v up to 11)
    k_for_v = 12
    kk_idx = list(k_2_values_g).index(k_for_v)
    profile_12 = mod_2k_total[offsets_2_g[kk_idx]:offsets_2_g[kk_idx+1]]
    total_p12 = profile_12.sum()
    log(f"\n  Using mod-2^{k_for_v} profile: {total_p12:,} odd-residue visits")

    # Build v-distribution from lookup
    v_table_12 = build_v_lookup(k_for_v)
    v_counts = {}
    for r, v in v_table_12.items():
        v_counts.setdefault(v, 0)
        v_counts[v] += int(profile_12[r])

    log(f"\n  P(v = j) from mod-2^12 mass profile:")
    log(f"  {'j':>3}  {'count':>12}  {'P(v=j)':>10}  {'Geom(1/2)':>10}  {'ratio_emp/geom':>14}")
    j_max = 11
    sum_known = 0
    for j in range(1, j_max + 1):
        cnt = v_counts.get(j, 0)
        p_emp = cnt / total_p12
        p_geom = 0.5 ** j
        ratio = p_emp / p_geom if p_geom > 0 else float('nan')
        sum_known += cnt
        log(f"  {j:>3}  {cnt:>12,}  {p_emp:>10.6f}  {p_geom:>10.6f}  {ratio:>14.4f}")

    n_undet = v_counts.get(-1, 0)
    log(f"  v >= {j_max+1} (undetermined): {n_undet:,} ({n_undet/total_p12:.6f})")

    # ============ Step 4: mod-3 mass-fractions ============
    log("\n" + "=" * 80)
    log("STEP 4: mod-3 mass-fractions (a, b, c) from orbit visits")
    log("=" * 80)

    a_emp = mod_3_total[0] / total_visits
    b_emp = mod_3_total[1] / total_visits
    c_emp = mod_3_total[2] / total_visits
    log(f"\n  m mod 3 = 0 (a): {mod_3_total[0]:>12,}  fraction = {a_emp:.6f}")
    log(f"  m mod 3 = 1 (b): {mod_3_total[1]:>12,}  fraction = {b_emp:.6f}")
    log(f"  m mod 3 = 2 (c): {mod_3_total[2]:>12,}  fraction = {c_emp:.6f}")

    log(f"\n  Note: Syracuse map sends odd m -> (3m+1)/2^v, which is coprime to 3.")
    log(f"  So ALL trajectory visits are at m coprime to 3 (a should be ~0).")
    log(f"  R65 reported (0.007, 0.347, 0.646) at value-truncation N=2^22.")
    log(f"  Asymptotic (0, 1/3, 2/3) as N -> infinity.")
    log(f"  Our values: a={a_emp:.6f}, b={b_emp:.6f}, c={c_emp:.6f}")

    # ============ Step 5: Direct Fourier from accumulated sums ============
    log("\n" + "=" * 80)
    log("STEP 5: |mu_hat(a/3^k)|^2 from arithmetic-deterministic profile (direct)")
    log("=" * 80)

    rows_compare = []
    log(f"\n  Direct computation: mu_hat(a/3^k) = (1/Z) sum_visits exp(2*pi*i*a*m/3^k)")
    log(f"  Z = {total_visits:,} visits\n")

    fourier_arith = {}  # (k, a) -> |mu_hat|^2
    for k in range(1, max_k_3 + 1):
        qk = 3 ** k
        off_k = offsets_3_g[k - 1]
        log(f"  --- k = {k}, 3^k = {qk}, primitive a (coprime to 3): {sum(1 for a in range(1, qk) if math.gcd(a, qk) == 1)} ---")
        primitives = [a for a in range(1, qk) if math.gcd(a, qk) == 1]
        # Display sample
        if k <= 3:
            log(f"  {'a':>4}  {'|mu_hat|^2':>12}")
        for a in primitives:
            re = fourier_re_total[off_k + a] / total_visits
            im = fourier_im_total[off_k + a] / total_visits
            mag_sq = re*re + im*im
            fourier_arith[(k, a)] = mag_sq
            if k <= 3:
                log(f"  {a:>4}  {mag_sq:>12.6f}")

        avg = sum(fourier_arith[(k, a)] for a in primitives) / len(primitives)
        max_val = max(fourier_arith[(k, a)] for a in primitives)
        min_val = min(fourier_arith[(k, a)] for a in primitives)
        log(f"  k={k}: mean = {avg:.6f}, max = {max_val:.6f}, min = {min_val:.6f}, n_primitive = {len(primitives)}")

    # ============ Step 6: Compare to R66 analytical (Geom-based) ============
    log("\n" + "=" * 80)
    log("STEP 6: Arithmetic-deterministic vs R66 Geom(1/2)-based analytical")
    log("=" * 80)

    df_R66 = load_R66_analytical()
    log(f"\n  Loaded R66 decay_law_derivation.csv: {df_R66.height} rows")

    log(f"\n  Per-coefficient comparison (k=1..4):")
    log(f"  {'k':>3}  {'a':>4}  {'R66_anal':>10}  {'R66_emp_subtree':>16}  {'Arith_det':>10}  {'gap_anal':>9}  {'gap_emp':>9}")
    for k in range(1, min(5, max_k_3 + 1)):
        qk = 3 ** k
        primitives = [a for a in range(1, qk) if math.gcd(a, qk) == 1]
        for a in primitives[:8]:  # display first 8
            r66_row = df_R66.filter((pl.col('k') == k) & (pl.col('a') == a))
            if r66_row.height == 0: continue
            r66_anal = float(r66_row['analytical_mu_hat_sq'][0])
            r66_emp = float(r66_row['empirical_mu_hat_sq'][0])
            arith = fourier_arith.get((k, a), float('nan'))
            gap_a = arith - r66_anal
            gap_e = arith - r66_emp
            log(f"  {k:>3}  {a:>4}  {r66_anal:>10.6f}  {r66_emp:>16.6f}  {arith:>10.6f}  {gap_a:>+9.4f}  {gap_e:>+9.4f}")
            rows_compare.append({'k': k, 'a': a,
                                'R66_analytical_geom': r66_anal,
                                'R66_empirical_subtree_R58': r66_emp,
                                'arithmetic_deterministic_R86': arith,
                                'gap_arith_minus_anal': gap_a,
                                'gap_arith_minus_emp': gap_e})

    # Aggregate stats per k
    log(f"\n  Aggregate per k (mean over primitive a):")
    log(f"  {'k':>3}  {'<R66_anal>':>11}  {'<R66_emp>':>11}  {'<Arith_det>':>11}  {'arith/anal':>10}  {'arith/emp':>10}")
    for k in range(1, max_k_3 + 1):
        qk = 3 ** k
        primitives = [a for a in range(1, qk) if math.gcd(a, qk) == 1]
        r66_anal_vals = [float(df_R66.filter((pl.col('k') == k) & (pl.col('a') == a))['analytical_mu_hat_sq'][0])
                         for a in primitives if df_R66.filter((pl.col('k') == k) & (pl.col('a') == a)).height > 0]
        r66_emp_vals = [float(df_R66.filter((pl.col('k') == k) & (pl.col('a') == a))['empirical_mu_hat_sq'][0])
                         for a in primitives if df_R66.filter((pl.col('k') == k) & (pl.col('a') == a)).height > 0]
        arith_vals = [fourier_arith.get((k, a), float('nan')) for a in primitives]
        if not r66_anal_vals: continue
        ma = np.mean(r66_anal_vals); me = np.mean(r66_emp_vals); mar = np.mean(arith_vals)
        log(f"  {k:>3}  {ma:>11.6f}  {me:>11.6f}  {mar:>11.6f}  {mar/ma:>10.4f}  {mar/me:>10.4f}")

    # Save comparison
    pl.DataFrame(rows_compare).write_csv(EXP_OUT / "86_arithmetic_vs_geom_comparison.csv")
    log(f"  [save] 86_arithmetic_vs_geom_comparison.csv ({len(rows_compare)} rows)")

    # ============ Step 7: Decay law and S_infinity ============
    log("\n" + "=" * 80)
    log("STEP 7: Decay law verification — <|mu_hat(a/3^k)|^2>_a * 3^(k-1) -> S_infinity")
    log("=" * 80)

    log(f"\n  Under arithmetic-deterministic mass profile:")
    log(f"  Formula: <|mu_hat(a/3^k)|^2>_a = S_inf / (2 * 3^(k-1))  =>  S_inf = <...> * 2 * 3^(k-1)")
    log(f"  {'k':>3}  {'<arith>':>11}  {'2*3^(k-1)':>10}  {'S_inf':>10}  {'<R66_anal>':>11}  {'S_inf_R66':>10}")
    s_inf_rows = []
    for k in range(1, max_k_3 + 1):
        qk = 3 ** k
        primitives = [a for a in range(1, qk) if math.gcd(a, qk) == 1]
        arith_mean = np.mean([fourier_arith.get((k, a), 0) for a in primitives])
        factor = 2 * 3 ** (k - 1)
        product = arith_mean * factor
        # also compute from R66_analytical for comparison
        r66_anal_vals = [float(df_R66.filter((pl.col('k') == k) & (pl.col('a') == a))['analytical_mu_hat_sq'][0])
                         for a in primitives if df_R66.filter((pl.col('k') == k) & (pl.col('a') == a)).height > 0]
        r66_mean = np.mean(r66_anal_vals) if r66_anal_vals else float('nan')
        s_inf_R66 = r66_mean * factor
        log(f"  {k:>3}  {arith_mean:>11.6f}  {factor:>10d}  {product:>10.6f}  {r66_mean:>11.6f}  {s_inf_R66:>10.6f}")
        s_inf_rows.append({'k': k, 'mean_arith_mu_hat_sq': arith_mean,
                          'factor_2_times_3_pow_k_minus_1': factor,
                          'S_inf_estimate_arith': product,
                          'mean_R66_analytical': r66_mean,
                          'S_inf_estimate_R66_anal': s_inf_R66})

    pl.DataFrame(s_inf_rows).write_csv(EXP_OUT / "86_S_infinity_arithmetic.csv")
    log(f"\n  [save] 86_S_infinity_arithmetic.csv")

    # Compare to 7/15 conjecture
    target = 7.0 / 15.0
    log(f"\n  Target conjecture: S_infinity = 7/15 = {target:.6f}")
    final_S_inf = s_inf_rows[-1]['S_inf_estimate_arith'] if s_inf_rows else float('nan')
    log(f"  Arithmetic-deterministic last estimate (k={max_k_3}): {final_S_inf:.6f}")
    log(f"  Gap vs 7/15: {final_S_inf - target:+.6f}")
    # Best estimate from middle k (less finite-N noise)
    if len(s_inf_rows) >= 4:
        mid_S = np.mean([s_inf_rows[2]['S_inf_estimate_arith'], s_inf_rows[3]['S_inf_estimate_arith']])
        log(f"  Mean of k=3,4 (less noisy): {mid_S:.6f}, gap vs 7/15: {mid_S - target:+.6f}")

    # ============ Step 8: Verdict ============
    log("\n" + "=" * 80)
    log("STEP 8: VERDICT")
    log("=" * 80)

    # Test outcome alpha: arithmetic matches empirical exactly within finite-N noise
    # Use k=2 as primary test (8 primitives, well-determined)
    qk_test = 9
    primitives_test = [a for a in range(1, qk_test) if math.gcd(a, qk_test) == 1]
    diffs_emp = [abs(fourier_arith[(2, a)] - float(df_R66.filter((pl.col('k') == 2) & (pl.col('a') == a))['empirical_mu_hat_sq'][0]))
                 for a in primitives_test if df_R66.filter((pl.col('k') == 2) & (pl.col('a') == a)).height > 0]
    diffs_anal = [abs(fourier_arith[(2, a)] - float(df_R66.filter((pl.col('k') == 2) & (pl.col('a') == a))['analytical_mu_hat_sq'][0]))
                  for a in primitives_test if df_R66.filter((pl.col('k') == 2) & (pl.col('a') == a)).height > 0]
    mean_diff_emp = np.mean(diffs_emp); max_diff_emp = max(diffs_emp)
    mean_diff_anal = np.mean(diffs_anal); max_diff_anal = max(diffs_anal)
    log(f"\n  At k=2 (8 primitive a values):")
    log(f"    arith vs R66_empirical (R58 inverse-tree weights):")
    log(f"      mean |diff| = {mean_diff_emp:.6f}, max |diff| = {max_diff_emp:.6f}")
    log(f"    arith vs R66_analytical (Geom(1/2) Markov chain):")
    log(f"      mean |diff| = {mean_diff_anal:.6f}, max |diff| = {max_diff_anal:.6f}")

    if mean_diff_emp < 0.005:
        log(f"\n  Arithmetic-deterministic MATCHES R58 empirical to within finite-N noise.")
        log(f"  Both compute the same Fourier coefficient via different proxies for")
        log(f"  the trajectory measure (R58: inverse-tree subtree weights; R86: forward")
        log(f"  orbit visits at N=2^32). Convergence to the same answer confirms both")
        log(f"  measure-the-same-thing.")
    else:
        log(f"\n  Arithmetic-deterministic and R58 empirical DIFFER measurably.")
        log(f"  Indicates the two trajectory measure proxies are not the same object,")
        log(f"  or one (or both) has finite-N truncation effects.")

    if mean_diff_anal < mean_diff_emp:
        log(f"\n  Surprising: Arith closer to Geom(1/2) analytical than to R58 empirical")
    elif mean_diff_emp < 0.005:
        log(f"\n  Verdict (alpha): arithmetic-deterministic gives exact predictions")
        log(f"  within finite-N noise. Geom(1/2) heuristic is a 5%-precision shortcut")
        log(f"  to the same target.")
    else:
        log(f"\n  Verdict (beta): arithmetic-deterministic differs from both Geom and R58")
        log(f"  Need further investigation of which is the 'true' trajectory measure.")

    (EXP_OUT / "86_arithmetic_deterministic_log.txt").write_text(
        "\n".join(results_log), encoding="utf-8")
    log(f"\n  [save] 86_arithmetic_deterministic_log.txt")


if __name__ == "__main__":
    main()
