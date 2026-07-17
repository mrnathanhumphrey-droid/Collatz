"""
prime_vs_all.py — compare Lagarias-class observables (W_j, ⟨v|j⟩, w_q) when
restricted to PRIME starting values vs all starting values.

Hypothesis (null): primes don't carry distinctive Lagarias-class signal — the
trajectory measure deviations come from 2-adic / mod-2^k residue structure,
not from prime/composite distinction.

Test: at N=2^32, sample ~1M PRIME starting values, walk orbits, compute
(W_j, ⟨v|j⟩, P(j)) per j ∈ {2, 4, 5}. Compare to all-starts (Result 30).

Method:
  - Use Miller-Rabin primality test (deterministic for n < 2^32 via {2,7,61})
  - Sample odd integers from [3, N], filter to primes
  - Cache batches to amortize numba startup
  - Walk 1M prime orbits, aggregate per-j
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2 * LOG2 - LOG3
MAX_VAL = np.int64(2**62)


@njit(cache=True)
def modexp(base, exp, mod):
    """Compute (base^exp) mod mod."""
    result = 1
    base = base % mod
    while exp > 0:
        if exp & 1:
            result = (result * base) % mod
        exp >>= 1
        base = (base * base) % mod
    return result


@njit(cache=True)
def miller_rabin(n):
    """Miller-Rabin with witnesses {2, 7, 61} — deterministic for n < 4759123141.
    Our N = 2^32 = 4294967296 < 4759123141, so deterministic."""
    if n < 2:
        return False
    if n == 2 or n == 3 or n == 7 or n == 61:
        return True
    if n % 2 == 0:
        return False
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 7, 61):
        if a >= n:
            continue
        x = modexp(a, d, n)
        if x == 1 or x == n - 1:
            continue
        composite = True
        for _ in range(r - 1):
            x = (x * x) % n
            if x == n - 1:
                composite = False
                break
        if composite:
            return False
    return True


@njit(parallel=True, cache=True)
def filter_primes(candidates):
    """Mark which candidates are prime."""
    n = len(candidates)
    is_prime = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        if miller_rabin(candidates[i]):
            is_prime[i] = True
    return is_prime


@njit(parallel=True, cache=True)
def walk(starts, max_steps=400000):
    n = len(starts)
    sigma_arr = np.zeros(n, dtype=np.int32)
    j_arr = np.full(n, -1, dtype=np.int8)
    v_sum_arr = np.zeros(n, dtype=np.int64)
    for i in prange(n):
        m = starts[i]
        steps = 0
        j_attr = -1
        v_total = 0
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
        sigma_arr[i] = steps
        j_arr[i] = j_attr
        v_sum_arr[i] = v_total
    return sigma_arr, j_arr, v_sum_arr


def sample_primes(N, n_target, seed):
    """Sample n_target primes from [3, N]."""
    rng = np.random.default_rng(seed)
    primes_collected = []
    batch_size = max(int(n_target / 0.04), 1_000_000)  # ~4% prime density at N=2^32
    while len(primes_collected) < n_target:
        candidates = 2 * rng.integers(1, (N - 1) // 2, size=batch_size, dtype=np.int64) + 1
        is_p = filter_primes(candidates)
        primes_batch = candidates[is_p]
        primes_collected.extend(primes_batch.tolist())
        if len(primes_collected) >= n_target:
            break
    return np.array(primes_collected[:n_target], dtype=np.int64)


def analyze_starts(starts):
    sigma, j_arr, v_sum = walk(starts)
    valid = sigma > 0
    sigma = sigma[valid].astype(np.float64)
    j_arr = j_arr[valid]
    v_sum = v_sum[valid].astype(np.float64)
    starts_v = starts[valid]
    log_m = np.log(starts_v.astype(np.float64))
    v_avg = v_sum / np.maximum(sigma, 1.0)

    results = {}
    overall = {
        'mean_log_m': float(log_m.mean()),
        'mean_sigma': float(sigma.mean()),
        'mean_v_avg_per_orbit': float(v_avg.mean()),
        'n_total': int(valid.sum()),
    }
    for jt in [2, 4, 5]:
        mask = j_arr == jt
        n_j = int(mask.sum())
        if n_j < 10:
            continue
        m_j = (4**jt - 1) // 3
        log_m_j = np.log(float(m_j))
        ml = float(log_m[mask].mean())
        ms = float(sigma[mask].mean())
        mv = float(v_avg[mask].mean())
        W_j = ms - ml / LOG43 - 1 + log_m_j / LOG43
        results[jt] = {
            'n': n_j,
            'P_j': n_j / valid.sum(),
            'mean_log_m': ml,
            'mean_sigma': ms,
            'mean_v': mv,
            'W_j': W_j,
        }
    return overall, results


def main():
    log2N = 32
    N = 1 << log2N
    n_starts = 1_000_000
    seeds = [42, 137, 271]
    j_targets = [2, 4, 5]

    print(f"# prime-vs-all comparison at N=2^{log2N}, {n_starts:,} orbits per condition × {len(seeds)} seeds")
    print(f"# Empirical W_j (50M-orbit, all starts at N=2^36): W_2=+7.156, W_4=-4.755, W_5=+4.590\n")

    # Run all-starts and prime-only across seeds
    all_overall = []
    all_per_j = {jt: [] for jt in j_targets}
    prime_overall = []
    prime_per_j = {jt: [] for jt in j_targets}

    for seed in seeds:
        # All-starts
        t0 = time.perf_counter()
        rng = np.random.default_rng(seed)
        all_starts = 2 * rng.integers(1, (N - 1) // 2, size=n_starts, dtype=np.int64) + 1
        all_ov, all_pj = analyze_starts(all_starts)
        t_all = time.perf_counter() - t0
        all_overall.append(all_ov)
        for jt in j_targets:
            if jt in all_pj:
                all_per_j[jt].append(all_pj[jt])

        # Prime-only
        t0 = time.perf_counter()
        prime_starts = sample_primes(N, n_starts, seed + 1000)
        t_prime_sample = time.perf_counter() - t0
        t0 = time.perf_counter()
        prime_ov, prime_pj = analyze_starts(prime_starts)
        t_prime_walk = time.perf_counter() - t0
        prime_overall.append(prime_ov)
        for jt in j_targets:
            if jt in prime_pj:
                prime_per_j[jt].append(prime_pj[jt])

        print(f"  seed={seed:>5}  all={t_all:.1f}s  prime_sample={t_prime_sample:.1f}s  prime_walk={t_prime_walk:.1f}s")

    # Aggregate
    def agg_field(arr_list, field):
        return np.array([r[field] for r in arr_list])

    print(f"\n# Overall observables (3-seed bootstrap mean ± SE):")
    print(f"#   {'observable':>20}  {'all_starts':>20}  {'prime_starts':>20}  {'gap':>10}  {'z':>6}")
    for field, name in [('mean_log_m', '⟨log m⟩'), ('mean_sigma', '⟨σ_S⟩'),
                          ('mean_v_avg_per_orbit', '⟨v_avg⟩')]:
        a = agg_field(all_overall, field)
        p = agg_field(prime_overall, field)
        a_mean = a.mean(); a_se = a.std(ddof=1) / np.sqrt(len(a)) if len(a) > 1 else 0
        p_mean = p.mean(); p_se = p.std(ddof=1) / np.sqrt(len(p)) if len(p) > 1 else 0
        gap = p_mean - a_mean
        se_gap = np.sqrt(a_se**2 + p_se**2)
        z = gap / se_gap if se_gap > 0 else 0
        print(f"   {name:>20}  {a_mean:>+11.4f}±{a_se:.4f}    {p_mean:>+11.4f}±{p_se:.4f}    "
              f"{gap:>+8.4f}  {z:>+6.2f}")

    print(f"\n# Per-j observables (3-seed bootstrap):")
    for jt in j_targets:
        print(f"\n  --- j={jt} ---")
        a_list = all_per_j[jt]
        p_list = prime_per_j[jt]
        print(f"   {'observable':>20}  {'all_starts':>20}  {'prime_starts':>20}  {'gap':>10}  {'z':>6}")
        for field, name in [('P_j', 'P(j)'), ('mean_log_m', '⟨log m|j⟩'),
                            ('mean_sigma', '⟨σ_S|j⟩'), ('mean_v', '⟨v|j⟩'),
                            ('W_j', f'W_{jt}')]:
            a = agg_field(a_list, field)
            p = agg_field(p_list, field)
            a_mean = a.mean(); a_se = a.std(ddof=1)/np.sqrt(len(a)) if len(a) > 1 else 0
            p_mean = p.mean(); p_se = p.std(ddof=1)/np.sqrt(len(p)) if len(p) > 1 else 0
            gap = p_mean - a_mean
            se_gap = np.sqrt(a_se**2 + p_se**2)
            z = gap / se_gap if se_gap > 0 else 0
            print(f"   {name:>20}  {a_mean:>+11.4f}±{a_se:.4f}    {p_mean:>+11.4f}±{p_se:.4f}    "
                  f"{gap:>+8.4f}  {z:>+6.2f}")

    print(f"\n# Verdict:")
    print(f"#   Null (primes don't differ): all gaps |z| < 2 across all observables")
    print(f"#   Signal: at least one observable shows |z| > 3 with structural pattern")


if __name__ == "__main__":
    main()
