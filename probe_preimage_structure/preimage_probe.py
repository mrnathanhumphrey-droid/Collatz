"""
Probe: Preimage structure of the Syracuse step on (Z/3^k)*.

Empirically construct K_k via M=2^20 integer lifts per coprime state x:
  for each lift n = x + i*3^k (i=0..M-1):
    v = v_2(3n+1)
    target = ((3n+1) >> v) mod 3^k
    K[x, target] += 1/M

Then study:
  Phase 1: per-y preimage set, max K_k[x,y], dominant v values.
  Phase 2: scaling of mean/median/max |Preimage(y)| across k=5,6,7.
  Phase 3: empirical v-distribution vs Geom(1/2) prediction P(v=j) = 2^(-j).
  Phase 4: scaling fit + outcome.

Sanity check: at k=5, compare K_emp to algebraic K_k (truncated-Geom) — should
match within sampling noise.
"""

import csv
import sys
import time
from pathlib import Path
from collections import defaultdict

import numpy as np
import numba

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = Path("C:/Collatz/probe_preimage_structure")
OUTDIR.mkdir(exist_ok=True)


@numba.njit(cache=True, parallel=False)
def build_K_emp_kernel(coprime, idx_arr, N, M_lifts):
    """Numba-jit empirical K via lifts. Returns (K, v_counter[j]).
    v_counter[j] = total count of v=j across all (x, lift) pairs.
    """
    n_coprime = len(coprime)
    K = np.zeros((n_coprime, n_coprime), dtype=np.float64)
    v_counter = np.zeros(64, dtype=np.int64)  # support up to v=63

    for i in range(n_coprime):
        x = coprime[i]
        for lift in range(M_lifts):
            n = x + lift * N
            num = 3 * n + 1
            v = 0
            while num & 1 == 0:
                num >>= 1
                v += 1
            target = num % N
            j = idx_arr[target]
            if j >= 0:
                K[i, j] += 1.0
            if v < 64:
                v_counter[v] += 1
    K /= M_lifts
    return K, v_counter


def build_K_emp(k, M_lifts):
    N = 3 ** k
    coprime = np.array([r for r in range(N) if r % 3 != 0], dtype=np.int64)
    n_coprime = len(coprime)
    idx_arr = -np.ones(N, dtype=np.int64)
    for i in range(n_coprime):
        idx_arr[coprime[i]] = i
    K, v_counter = build_K_emp_kernel(coprime, idx_arr, N, M_lifts)
    return K, coprime, idx_arr, v_counter


def order_of_2_mod(N):
    o, p = 1, 2 % N
    while p != 1:
        p = (p * 2) % N
        o += 1
    return o


def build_K_algebraic(k):
    """Truncated-Geom K_k for sanity-checking against empirical."""
    N = 3 ** k
    M = order_of_2_mod(N)
    coprime = [r for r in range(N) if r % 3 != 0]
    n = len(coprime)
    idx = {r: i for i, r in enumerate(coprime)}
    inv2 = pow(2, -1, N)
    inv2_pow = [1, inv2]
    for v in range(2, M + 1):
        inv2_pow.append((inv2_pow[-1] * inv2) % N)
    Z = 1.0 - 2.0 ** (-M)
    weights = [None] + [(2.0 ** (-v)) / Z for v in range(1, M + 1)]
    K = np.zeros((n, n))
    for i, r in enumerate(coprime):
        base_mod = (3 * r + 1) % N
        for v in range(1, M + 1):
            target = (base_mod * inv2_pow[v]) % N
            j = idx[target]
            K[i, j] += weights[v]
    return K, M


def preimage_stats(K, threshold=0.0):
    """For each y (column), compute (n_preimages, max_K, entropy, total_mass)."""
    n = K.shape[0]
    stats = []
    for y in range(n):
        col = K[:, y]
        mask = col > threshold
        n_pre = int(mask.sum())
        if n_pre == 0:
            stats.append((0, 0.0, 0.0, 0.0))
            continue
        max_K = float(col[mask].max())
        total = float(col[mask].sum())
        # Normalize and compute entropy
        probs = col[mask] / total
        entropy = float(-np.sum(probs * np.log(probs + 1e-300)))
        stats.append((n_pre, max_K, entropy, total))
    return stats


def main():
    M_lifts = 1 << 20
    k_list = [5, 6, 7]

    print(f"Empirical K_k via M = 2^20 = {M_lifts} integer lifts per coprime state")
    print()

    summary = []
    full_stats = {}
    full_v_counters = {}

    for k in k_list:
        N = 3 ** k
        n_coprime = 2 * 3 ** (k - 1)
        print(f"=== k = {k}  (N=3^k={N}, n_coprime={n_coprime}) ===")

        # Build K empirically
        t0 = time.time()
        K_emp, coprime, idx_arr, v_counter = build_K_emp(k, M_lifts)
        elapsed = time.time() - t0
        print(f"  K_emp built in {elapsed:.1f}s; total_lifts = {n_coprime * M_lifts:.2e}")

        # Sanity check at k=5: compare to algebraic K_k
        if k == 5:
            t0 = time.time()
            K_alg, M_alg = build_K_algebraic(k)
            elapsed = time.time() - t0
            diff = np.linalg.norm(K_emp - K_alg, ord='fro')
            rel = diff / np.linalg.norm(K_alg, ord='fro')
            print(f"  K_alg built in {elapsed:.1f}s; ord_{{3^k}}(2) = {M_alg}")
            print(f"  ||K_emp - K_alg||_F / ||K_alg||_F = {rel:.4e}")

        # Phase 1 + 2: preimage stats per y
        # Two thresholds: strict (>0) and effective (>= 1/M = 2^-20)
        stats_strict = preimage_stats(K_emp, threshold=0.0)
        stats_eff = preimage_stats(K_emp, threshold=1.0 / M_lifts * 0.5)  # any non-zero count

        n_pre_strict = np.array([s[0] for s in stats_strict])
        n_pre_eff = np.array([s[0] for s in stats_eff])
        max_K = np.array([s[1] for s in stats_strict])
        entropies = np.array([s[2] for s in stats_strict])
        total_masses = np.array([s[3] for s in stats_strict])

        print(f"  |Preimage(y)| (strict, K>0):    "
              f"min={n_pre_strict.min()}, median={int(np.median(n_pre_strict))}, "
              f"mean={n_pre_strict.mean():.2f}, max={n_pre_strict.max()}")
        print(f"  |Preimage(y)| (effective, K>=1/M): "
              f"min={n_pre_eff.min()}, median={int(np.median(n_pre_eff))}, "
              f"mean={n_pre_eff.mean():.2f}, max={n_pre_eff.max()}")
        print(f"  max K_k[x,y] over (x,y):  median={float(np.median(max_K)):.4f}, "
              f"mean={max_K.mean():.4f}, max={max_K.max():.4f}")
        print(f"  entropy of preimage dist:  median={float(np.median(entropies)):.4f}, "
              f"mean={entropies.mean():.4f}")
        print(f"  log(|Preimage|) baseline:  median={float(np.log(np.median(n_pre_strict))):.4f}")
        # entropy / log|preimage| ratio: 1.0 = uniform spread, 0 = concentrated
        with np.errstate(divide='ignore', invalid='ignore'):
            ent_ratio = np.where(n_pre_strict > 1, entropies / np.log(n_pre_strict), 0)
        print(f"  entropy / log|Preimage|:   median={float(np.median(ent_ratio)):.4f}, "
              f"mean={ent_ratio.mean():.4f}")

        # Save per-y CSV
        with open(OUTDIR / f"preimage_map_k{k}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["y_idx", "y", "n_preimage_strict", "n_preimage_eff",
                        "max_K", "entropy", "total_mass", "entropy_over_log_npre"])
            for y in range(K_emp.shape[0]):
                er = (entropies[y] / np.log(n_pre_strict[y])) if n_pre_strict[y] > 1 else 0.0
                w.writerow([y, int(coprime[y]), int(n_pre_strict[y]), int(n_pre_eff[y]),
                            f"{max_K[y]:.6e}", f"{entropies[y]:.6f}",
                            f"{total_masses[y]:.6e}", f"{er:.6f}"])

        # Phase 3: empirical v-distribution
        total_v = v_counter.sum()
        print(f"\n  Empirical P(v=j) (top 12) vs Geom(1/2) = 2^(-j):")
        print(f"  {'j':>3} {'count':>12} {'P(v=j)_emp':>14} {'P(v=j)_geom':>14} {'ratio':>10}")
        v_rows = []
        for j in range(1, 21):
            p_emp = v_counter[j] / total_v if total_v > 0 else 0.0
            p_geom = 2.0 ** (-j)
            ratio = p_emp / p_geom if p_geom > 0 else float('nan')
            if j <= 12:
                print(f"  {j:>3} {v_counter[j]:>12} {p_emp:>14.6e} {p_geom:>14.6e} {ratio:>10.4f}")
            v_rows.append([j, int(v_counter[j]), p_emp, p_geom, ratio])

        with open(OUTDIR / f"v_distribution_empirical_k{k}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["v", "count", "P_emp", "P_geom", "ratio"])
            for row in v_rows:
                w.writerow([row[0], row[1], f"{row[2]:.6e}", f"{row[3]:.6e}", f"{row[4]:.6f}"])

        summary.append({
            "k": k,
            "n_coprime": n_coprime,
            "preimage_strict_mean": float(n_pre_strict.mean()),
            "preimage_strict_median": float(np.median(n_pre_strict)),
            "preimage_strict_max": int(n_pre_strict.max()),
            "preimage_strict_min": int(n_pre_strict.min()),
            "preimage_eff_mean": float(n_pre_eff.mean()),
            "preimage_eff_median": float(np.median(n_pre_eff)),
            "preimage_eff_max": int(n_pre_eff.max()),
            "max_K_median": float(np.median(max_K)),
            "entropy_mean": float(entropies.mean()),
            "entropy_over_log_npre_median": float(np.median(ent_ratio)),
        })
        full_stats[k] = (n_pre_strict, n_pre_eff, max_K, entropies)
        full_v_counters[k] = v_counter
        print()

    # Phase 4: scaling fits
    print("=" * 70)
    print("PHASE 4: scaling fits")
    print("=" * 70)

    ks = np.array([s["k"] for s in summary], dtype=float)
    means_strict = np.array([s["preimage_strict_mean"] for s in summary])
    means_eff = np.array([s["preimage_eff_mean"] for s in summary])
    maxes_strict = np.array([s["preimage_strict_max"] for s in summary])

    # Linear fit: |Preimage|(k) ≈ a·k + b
    a_strict, b_strict = np.polyfit(ks, means_strict, 1)
    a_eff, b_eff = np.polyfit(ks, means_eff, 1)
    a_max, b_max = np.polyfit(ks, maxes_strict, 1)

    # Log-log fit: |Preimage|(k) ≈ c · k^p
    log_k = np.log(ks)
    log_m_strict = np.log(means_strict)
    log_m_eff = np.log(means_eff)
    p_strict, log_c_strict = np.polyfit(log_k, log_m_strict, 1)
    p_eff, log_c_eff = np.polyfit(log_k, log_m_eff, 1)

    # Constant fit: |Preimage|(k) ≈ const
    const_strict = means_strict.mean()
    const_eff = means_eff.mean()

    # Compute residuals
    def rss(predicted, actual):
        return float(np.sum((np.array(predicted) - np.array(actual)) ** 2))
    rss_lin_strict = rss(a_strict * ks + b_strict, means_strict)
    rss_const_strict = rss([const_strict] * len(ks), means_strict)
    rss_lin_eff = rss(a_eff * ks + b_eff, means_eff)
    rss_const_eff = rss([const_eff] * len(ks), means_eff)

    print()
    print("Fit comparisons (RSS smaller = better):")
    print(f"  Strict (>0):    linear  a·k+b: a={a_strict:.4f}, b={b_strict:.4f}, RSS={rss_lin_strict:.4e}")
    print(f"                  log-log: c·k^p: p={p_strict:.4f}, c={np.exp(log_c_strict):.4f}")
    print(f"                  constant:      ≈{const_strict:.2f}, RSS={rss_const_strict:.4e}")
    print(f"  Effective (≥1/M): linear: a={a_eff:.4f}, b={b_eff:.4f}, RSS={rss_lin_eff:.4e}")
    print(f"                    log-log: p={p_eff:.4f}, c={np.exp(log_c_eff):.4f}")
    print(f"                    constant:    ≈{const_eff:.2f}, RSS={rss_const_eff:.4e}")
    print(f"  max_strict:     linear  a·k+b: a={a_max:.4f}, b={b_max:.4f}")

    # Save scaling CSV
    with open(OUTDIR / "scaling_analysis.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "n_coprime",
                    "preimage_strict_mean", "preimage_strict_median",
                    "preimage_strict_max", "preimage_strict_min",
                    "preimage_eff_mean", "preimage_eff_median", "preimage_eff_max",
                    "max_K_median", "entropy_mean", "entropy_over_log_npre_median"])
        for s in summary:
            w.writerow([s["k"], s["n_coprime"],
                        f"{s['preimage_strict_mean']:.4f}", f"{s['preimage_strict_median']:.0f}",
                        s["preimage_strict_max"], s["preimage_strict_min"],
                        f"{s['preimage_eff_mean']:.4f}", f"{s['preimage_eff_median']:.0f}",
                        s["preimage_eff_max"],
                        f"{s['max_K_median']:.6f}", f"{s['entropy_mean']:.6f}",
                        f"{s['entropy_over_log_npre_median']:.6f}"])
        w.writerow([])
        w.writerow(["fit_type", "param_a_or_p", "param_b_or_c", "rss"])
        w.writerow(["strict_linear", f"{a_strict:.6f}", f"{b_strict:.6f}", f"{rss_lin_strict:.6e}"])
        w.writerow(["strict_loglog_powerlaw", f"{p_strict:.6f}", f"{np.exp(log_c_strict):.6f}", ""])
        w.writerow(["strict_constant", f"{const_strict:.6f}", "", f"{rss_const_strict:.6e}"])
        w.writerow(["eff_linear", f"{a_eff:.6f}", f"{b_eff:.6f}", f"{rss_lin_eff:.6e}"])
        w.writerow(["eff_constant", f"{const_eff:.6f}", "", f"{rss_const_eff:.6e}"])

    # Print summary table
    print()
    print("Summary table:")
    print(f"  {'k':>3} {'n':>5} {'mean_strict':>12} {'med_strict':>12} {'max_strict':>12} "
          f"{'mean_eff':>10} {'max_eff':>10} {'ent/log':>10}")
    for s in summary:
        print(f"  {s['k']:>3} {s['n_coprime']:>5} "
              f"{s['preimage_strict_mean']:>12.2f} {s['preimage_strict_median']:>12.0f} "
              f"{s['preimage_strict_max']:>12} "
              f"{s['preimage_eff_mean']:>10.2f} {s['preimage_eff_max']:>10} "
              f"{s['entropy_over_log_npre_median']:>10.4f}")

    print()
    print(f"Outputs in {OUTDIR}")


if __name__ == "__main__":
    main()
