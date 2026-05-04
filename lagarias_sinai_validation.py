"""
lagarias_sinai_validation.py — Result 68 candidate.

High-precision empirical test of Lagarias-Sinai heuristic v_t ~ Geom(½).

Tests:
  1. Direct P̂(v=j) vs 2^(-j) at N ∈ {2^28, 2^30, 2^32, 2^34}
  3. v-parity P̂(even)/P̂(odd) vs (1/3, 2/3)         [R64 foundation]
  4. v mod 6, v mod 18 vs Geom(1/2) prediction       [R65 foundation]
  5. Higher moments E[v], Var[v], skewness, kurtosis
  6. σ-band-conditional P̂(v=j)
  7. Convergence rate |P̂(v=j) - 2^(-j)| vs N

Sample size scales with N: more orbits at small N (cheap), fewer at large N
(expensive walk). Aggregate v samples per N for direct comparison.
"""
import math
import sys
import time
from pathlib import Path
import numpy as np
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(r"C:\Collatz")
log_lines = []
def log(s):
    print(s, flush=True)
    log_lines.append(s)

LOG_2 = np.log(2.0)
K_H = 3.0 / np.log(4.0/3.0)
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_capture_v(starts, max_T):
    """
    For each orbit, record:
      v_seq[i, t] = v_2 at Syracuse step t  (int8, capped at 127)
      v_count[i] = number of Syracuse steps
      sigma_arr[i] = total Collatz steps (or -1 if failed/incomplete)
      log_n_arr[i] = log(start)
      r32_seq[i, t] = m_t mod 32 at the start of Syracuse step t
    """
    n = len(starts)
    v_seq = np.full((n, max_T), -1, dtype=np.int8)
    r32_seq = np.full((n, max_T), -1, dtype=np.int8)
    v_count = np.zeros(n, dtype=np.int32)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    log_n_arr = np.zeros(n, dtype=np.float64)
    for i in prange(n):
        m = np.int64(starts[i])
        log_n_arr[i] = np.log(np.float64(m))
        sigma = 0
        T = 0
        failed = False
        while (m & 1) == 0 and m > 1:
            m >>= 1
            sigma += 1
        if m == 1:
            sigma_arr[i] = sigma
            v_count[i] = 0
            continue
        while m != 1 and T < max_T:
            if m > MAX_VAL // 3:
                failed = True
                break
            r32_seq[i, T] = m & 31
            x = 3 * m + 1
            v = 0
            while (x & 1) == 0:
                x >>= 1
                v += 1
            v_seq[i, T] = v if v < 127 else 127
            sigma += 1 + v
            T += 1
            m = x
        if not failed and m == 1:
            sigma_arr[i] = sigma
            v_count[i] = T
    return v_seq, r32_seq, v_count, sigma_arr, log_n_arr


def aggregate_v_at_N(log2N, n_orbits, seed):
    """Walk n_orbits orbits at scale N=2^log2N, return aggregated arrays."""
    N = 1 << log2N
    rng = np.random.default_rng(seed)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_orbits, dtype=np.int64) + 1
    max_T = max(200, int(20 * np.log2(max(N, 2))))
    t0 = time.time()
    v_seq, r32_seq, v_count, sigma_arr, log_n_arr = walk_capture_v(
        starts, max_T)
    walk_t = time.time() - t0
    ok = (sigma_arr > 0) & (v_count > 0)
    n_ok = int(ok.sum())
    # Flatten v values across all ok orbits
    v_seq = v_seq[ok]
    r32_seq = r32_seq[ok]
    v_count = v_count[ok]
    sigma_arr = sigma_arr[ok]
    log_n_arr = log_n_arr[ok]
    flat_v = []
    flat_r32 = []
    flat_orbit = []
    for i in range(n_ok):
        c = int(v_count[i])
        flat_v.append(v_seq[i, :c])
        flat_r32.append(r32_seq[i, :c])
        flat_orbit.append(np.full(c, i, dtype=np.int32))
    flat_v = np.concatenate(flat_v).astype(np.int32)
    flat_r32 = np.concatenate(flat_r32).astype(np.int32)
    flat_orbit = np.concatenate(flat_orbit)
    return {
        'log2N': log2N,
        'n_orbits': n_ok,
        'walk_t': walk_t,
        'flat_v': flat_v,
        'flat_r32': flat_r32,
        'flat_orbit': flat_orbit,
        'v_count': v_count,
        'sigma_arr': sigma_arr,
        'log_n_arr': log_n_arr,
    }


def empirical_p_v(v_arr, max_j=30):
    """Returns array P̂(v=j) for j in 1..max_j."""
    counts = np.zeros(max_j + 1, dtype=np.int64)
    for j in range(1, max_j + 1):
        counts[j] = (v_arr == j).sum()
    n = counts[1:].sum()
    return counts[1:] / max(n, 1), int(n)


def main():
    log("=" * 78)
    log("Lagarias-Sinai heuristic validation: v_t ~ Geom(1/2)")
    log("=" * 78)

    # Adjust orbit count by N to stay tractable
    runs = [
        (28, 2_000_000, 11),
        (30, 1_000_000, 22),
        (32,   500_000, 33),
        (34,   200_000, 44),
    ]

    log(f"\nStep 1: walking orbits at multiple N\n")
    aggs = []
    for log2N, n_orbits, seed in runs:
        agg = aggregate_v_at_N(log2N, n_orbits, seed)
        n_v = len(agg['flat_v'])
        log(f"  N=2^{log2N}: {n_orbits:,} orbits → {agg['n_orbits']:,} ok, "
            f"{n_v:,} v samples, walk {agg['walk_t']:.1f}s")
        aggs.append(agg)

    # ------- Step 1: direct distribution -------
    log("\n" + "=" * 78)
    log("Step 1: direct P̂(v=j) vs Geom(1/2) prediction 2^(-j)")
    log("=" * 78)

    log(f"\n  {'j':>3}  {'2^-j':>9}  " +
        "  ".join(f"N=2^{a['log2N']}" for a in aggs))
    distrs = {}
    for a in aggs:
        p, n = empirical_p_v(a['flat_v'], max_j=15)
        distrs[a['log2N']] = (p, n)

    for j in range(1, 16):
        line = f"  {j:>3}  {2**(-j):>9.6f}"
        for a in aggs:
            p, n = distrs[a['log2N']]
            line += f"   {p[j-1]:.6f}"
        log(line)

    # Per-j ratio P̂(v=j) / 2^(-j)
    log(f"\n  Ratio P̂(v=j) / 2^(-j) (should be ≈ 1):")
    log(f"  {'j':>3}  " + "  ".join(f"N=2^{a['log2N']}" for a in aggs))
    for j in range(1, 11):
        line = f"  {j:>3}"
        for a in aggs:
            p, n = distrs[a['log2N']]
            ratio = p[j-1] / (2**(-j))
            line += f"   {ratio:.6f}"
        log(line)

    # Chi-squared GoF at largest N
    a_max = aggs[-1]
    p_max, n_max = distrs[a_max['log2N']]
    chi2 = 0.0
    dof = 0
    for j in range(1, 16):
        expected = (2 ** (-j)) * n_max
        observed = p_max[j-1] * n_max
        if expected > 5:
            chi2 += (observed - expected) ** 2 / expected
            dof += 1
    log(f"\n  Chi-squared GoF at N=2^{a_max['log2N']}: chi^2 = {chi2:.2f}, "
        f"dof = {dof-1}, chi^2/dof = {chi2/(dof-1):.3f}")

    # ------- Step 3: v-parity (R64 foundation) -------
    log("\n" + "=" * 78)
    log("Step 3: v-parity test (R64 |μ̂(1/3)|² depends on this)")
    log("=" * 78)
    log(f"\n  Geom(1/2) prediction: P(v even) = 1/3 = 0.3333..., "
        f"P(v odd) = 2/3 = 0.6666...")
    log(f"\n  {'N':>10}  {'n_v':>13}  {'P̂(even)':>10}  {'P̂(odd)':>10}  "
        f"{'dev_even':>9}  {'dev_odd':>9}")
    for a in aggs:
        v = a['flat_v']
        n = len(v)
        p_even = (v % 2 == 0).mean()
        p_odd = 1 - p_even
        log(f"  N=2^{a['log2N']:>3}  {n:>13,}  {p_even:>10.6f}  "
            f"{p_odd:>10.6f}  {p_even - 1/3:>+9.6f}  "
            f"{p_odd - 2/3:>+9.6f}")

    # ------- Step 4: higher-q (R65 foundation) -------
    log("\n" + "=" * 78)
    log("Step 4: v mod 6 and v mod 18 (R65 framework)")
    log("=" * 78)
    log(f"\n  Geom(1/2) prediction for v mod ord:")
    log(f"     P(v ≡ r mod ord) = 2^(-r) / (1 - 2^(-ord)) for r in 1..ord")

    for ord_ in [2, 6, 18]:
        log(f"\n  v mod {ord_} (Geom(1/2) prediction tabulated):")
        log(f"  {'r':>3}  {'pred':>10}  " +
            "  ".join(f"N=2^{a['log2N']}" for a in aggs))
        # Predictions
        denom = 1 - 2**(-ord_)
        preds = []
        for r in range(1, ord_ + 1):
            preds.append(2**(-r) / denom)
        for r in range(1, ord_ + 1):
            line = f"  {r:>3}  {preds[r-1]:>10.6f}"
            for a in aggs:
                v = a['flat_v']
                # v mod ord, mapping 0 → ord
                mod_v = ((v - 1) % ord_) + 1
                p_r = (mod_v == r).mean()
                line += f"   {p_r:.6f}"
            log(line)

    # ------- Step 5: higher moments -------
    log("\n" + "=" * 78)
    log("Step 5: moments E[v], Var[v], E[(v-2)^3], excess kurtosis")
    log("=" * 78)
    log(f"\n  Geom(1/2) moments: E[v]=2.000, Var[v]=2.000, "
        f"E[(v-2)^3]=6.000, kurt=6.000")
    log(f"\n  {'N':>10}  {'mean':>9}  {'var':>9}  {'skew*sigma3':>11}  "
        f"{'excess_kurt':>11}")
    for a in aggs:
        v = a['flat_v'].astype(np.float64)
        mean = v.mean()
        var = v.var()
        cm3 = ((v - mean) ** 3).mean()
        cm4 = ((v - mean) ** 4).mean()
        excess_kurt = cm4 / var ** 2 - 3
        log(f"  N=2^{a['log2N']:>3}  {mean:>9.5f}  {var:>9.5f}  "
            f"{cm3:>11.5f}  {excess_kurt:>11.5f}")

    # ------- Step 6: σ-band-conditional -------
    log("\n" + "=" * 78)
    log("Step 6: σ-band-conditional P̂(v=j)")
    log("=" * 78)
    a_mid = aggs[2]  # use N=2^32
    sigma_resid = a_mid['sigma_arr'].astype(np.float64) - K_H * a_mid['log_n_arr']
    edges = np.percentile(sigma_resid, [25, 50, 75, 95])
    band_per_orbit = np.digitize(sigma_resid, edges)
    band_names = ['0-25', '25-50', '50-75', '75-95', '95-100']

    # Map flat_v back to orbit → band
    band_per_v = band_per_orbit[a_mid['flat_orbit']]
    log(f"\n  Conditional at N=2^{a_mid['log2N']}:")
    log(f"  {'band':>8}  {'n_v':>10}  " +
        "  ".join(f"j={j}" for j in range(1, 8)))
    for b in range(5):
        mask = band_per_v == b
        v_b = a_mid['flat_v'][mask]
        n_b = len(v_b)
        if n_b == 0:
            continue
        p_b, _ = empirical_p_v(v_b, max_j=10)
        line = f"  {band_names[b]:>8}  {n_b:>10,}"
        for j in range(1, 8):
            line += f"   {p_b[j-1]:.5f}"
        log(line)

    # Chi-squared per band vs Geom(1/2)
    log(f"\n  Per-band chi^2 vs Geom(1/2) (j=1..15, expected ≥ 5):")
    for b in range(5):
        mask = band_per_v == b
        v_b = a_mid['flat_v'][mask]
        n_b = len(v_b)
        if n_b == 0:
            continue
        p_b, _ = empirical_p_v(v_b, max_j=15)
        chi2_b = 0.0
        dof_b = 0
        for j in range(1, 16):
            expected = (2 ** (-j)) * n_b
            observed = p_b[j-1] * n_b
            if expected > 5:
                chi2_b += (observed - expected) ** 2 / expected
                dof_b += 1
        log(f"    {band_names[b]:>8}: chi^2/dof = {chi2_b/(dof_b-1):.3f}")

    # ------- Step 7: convergence rate -------
    log("\n" + "=" * 78)
    log("Step 7: convergence rate |P̂(v=j) - 2^(-j)| vs sample size")
    log("=" * 78)
    log(f"\n  At each N, deviation |P̂(v=j) - 2^(-j)| at small j:")
    log(f"  {'N':>10}  {'n_v':>13}  " +
        "  ".join(f"|dev j={j}|" for j in range(1, 6)))
    for a in aggs:
        p, n = distrs[a['log2N']]
        line = f"  N=2^{a['log2N']:>3}  {n:>13,}"
        for j in range(1, 6):
            dev = abs(p[j-1] - 2**(-j))
            line += f"   {dev:>9.6f}"
        log(line)

    # Test 1/sqrt(n_v) scaling (CLT)
    log(f"\n  Compare to 1/sqrt(n_v) (CLT prediction for binomial spread):")
    log(f"  {'N':>10}  {'1/sqrt(n_v)':>12}  {'avg|dev j=1..5|':>16}")
    for a in aggs:
        p, n = distrs[a['log2N']]
        avg_dev = np.mean([abs(p[j-1] - 2**(-j)) for j in range(1, 6)])
        log(f"  N=2^{a['log2N']:>3}  {1/np.sqrt(n):>12.6f}  "
            f"{avg_dev:>16.6f}")

    # ------- Step 2: r-conditional (smaller pass) -------
    log("\n" + "=" * 78)
    log("Step 2: r mod 32 conditional P̂(v=j)")
    log("=" * 78)
    log(f"\n  At N=2^{a_mid['log2N']}, P̂(v=1) | r mod 32 (Geom(1/2) → 0.5):")
    log(f"  {'r':>3}  {'n_at_r':>10}  {'P̂(v=1|r)':>10}  {'dev':>8}")
    for r in range(1, 32, 2):  # odd residues only (m is odd)
        mask = a_mid['flat_r32'] == r
        v_r = a_mid['flat_v'][mask]
        n_r = len(v_r)
        if n_r < 100:
            continue
        p_v1 = (v_r == 1).mean()
        log(f"  {r:>3}  {n_r:>10,}  {p_v1:>10.6f}  {p_v1 - 0.5:>+8.6f}")

    log(f"\n  Note: r=21 mod 32 expected to deviate (R49 finds V≥5 at r=21 visits).")

    # ------- Step 8: precision for downstream claims -------
    log("\n" + "=" * 78)
    log("Step 8: precision for downstream claims")
    log("=" * 78)

    a_max = aggs[-1]
    v = a_max['flat_v'].astype(np.float64)
    p_even = (v % 2 == 0).mean()
    p_odd = 1 - p_even
    p_even_dev = abs(p_even - 1/3)

    log(f"\n  At N=2^{a_max['log2N']} ({len(v):,} v samples):")
    log(f"")
    log(f"  (a) Tao K_h = 3/log(4/3) = {K_H:.4f}")
    log(f"      Depends on E[v] = 2.000 exactly.")
    log(f"      Empirical E[v] = {v.mean():.6f}, deviation = "
        f"{abs(v.mean() - 2):.6f}")
    log(f"      K_h sensitivity: dK_h/dE[v] ≈ -K_h/E[v] = {-K_H/2:.4f}")
    log(f"      Implied K_h precision: K_h = {K_H:.4f} ± "
        f"{abs(v.mean() - 2) * K_H/2:.4f}")
    log(f"")
    log(f"  (b) R64 |μ̂(1/3)|² depends on P(v even) = 1/3.")
    log(f"      Empirical P(v even) = {p_even:.6f}, deviation = "
        f"{p_even_dev:.6f}")
    log(f"      |μ̂(1/3)|² depends on this via path-counting; deviation "
        f"propagates as O(dev) in the closed form.")
    log(f"")
    log(f"  (c) R65 4^(-k) decay depends on geometric tail P(v=j) ~ 2^(-j).")
    log(f"      Empirical at j=1..5: deviations ≤ "
        f"{max(abs(distrs[a_max['log2N']][0][j-1] - 2**(-j)) for j in range(1,6)):.6f}")

    # ------- Save outputs -------
    log("\n" + "=" * 78)
    log("Saving outputs")
    log("=" * 78)

    with open(OUT / "v_distribution_by_N.csv", "w") as f:
        f.write("N,j,p_emp,p_geom,ratio,deviation,n_samples\n")
        for a in aggs:
            p, n = distrs[a['log2N']]
            for j in range(1, 16):
                f.write(f"2^{a['log2N']},{j},{p[j-1]:.8f},"
                        f"{2**(-j):.8f},"
                        f"{p[j-1]/(2**(-j)):.8f},"
                        f"{p[j-1] - 2**(-j):+.8e},"
                        f"{n}\n")
    log("  [wrote] v_distribution_by_N.csv")

    with open(OUT / "v_higher_moments.csv", "w") as f:
        f.write("N,n_samples,mean,var,central_3rd,excess_kurt\n")
        for a in aggs:
            v = a['flat_v'].astype(np.float64)
            f.write(f"2^{a['log2N']},{len(v)},{v.mean():.6f},"
                    f"{v.var():.6f},"
                    f"{((v-v.mean())**3).mean():.6f},"
                    f"{((v-v.mean())**4).mean()/v.var()**2 - 3:.6f}\n")
    log("  [wrote] v_higher_moments.csv")

    with open(OUT / "v_conditional_distributions.csv", "w") as f:
        f.write("type,key,j,p_emp,n_samples\n")
        for r in range(1, 32, 2):
            mask = a_mid['flat_r32'] == r
            v_r = a_mid['flat_v'][mask]
            n_r = len(v_r)
            if n_r < 100:
                continue
            for j in range(1, 11):
                p = (v_r == j).mean()
                f.write(f"r_mod_32,{r},{j},{p:.8f},{n_r}\n")
        for b in range(5):
            mask = band_per_v == b
            v_b = a_mid['flat_v'][mask]
            n_b = len(v_b)
            if n_b == 0:
                continue
            for j in range(1, 11):
                p = (v_b == j).mean()
                f.write(f"sigma_band,{band_names[b]},{j},{p:.8f},{n_b}\n")
    log("  [wrote] v_conditional_distributions.csv")

    with open(OUT / "convergence_rate.csv", "w") as f:
        f.write("N,n_samples,inv_sqrt_n,avg_dev_j1_to_5,max_dev_j1_to_5\n")
        for a in aggs:
            p, n = distrs[a['log2N']]
            devs = [abs(p[j-1] - 2**(-j)) for j in range(1, 6)]
            f.write(f"2^{a['log2N']},{n},{1/np.sqrt(n):.6f},"
                    f"{np.mean(devs):.6f},{max(devs):.6f}\n")
    log("  [wrote] convergence_rate.csv")

    # ------- Verdict -------
    log("\n" + "=" * 78)
    log("VERDICT")
    log("=" * 78)
    a_max = aggs[-1]
    p_max, n_max = distrs[a_max['log2N']]
    max_dev = max(abs(p_max[j-1] - 2**(-j)) for j in range(1, 11))
    v = a_max['flat_v'].astype(np.float64)
    p_even_dev = abs((v % 2 == 0).mean() - 1/3)
    mean_dev = abs(v.mean() - 2)
    log(f"\n  At N=2^{a_max['log2N']}, {n_max:,} v samples:")
    log(f"    max |P̂(v=j) - 2^(-j)| (j=1..10) = {max_dev:.6f}")
    log(f"    P̂(v even) - 1/3 = {(v%2==0).mean() - 1/3:+.6f}")
    log(f"    E[v] - 2.0 = {mean_dev:+.6f}")
    log(f"")
    if max_dev < 0.001 and mean_dev < 0.001:
        log(f"  Outcome (α): heuristic holds to high precision (≤ 0.001).")
    elif max_dev < 0.01 and p_even_dev < 0.01:
        log(f"  Outcome (β): heuristic holds approximately; ε ~ 0.001-0.01.")
    else:
        log(f"  Outcome (γ): significant deviations.")

    (OUT / "lagarias_sinai_validation_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )
    log("\n  [wrote] lagarias_sinai_validation_log.txt")


if __name__ == "__main__":
    main()
