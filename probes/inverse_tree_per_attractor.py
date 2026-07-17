"""
inverse_tree_per_attractor.py — test the spectral closure route for ⟨σ_S | j⟩.

Hypothesis (user pressure-test): per-attractor inverse-tree growth rate λ_j
closes ⟨σ_S | j⟩ via:
  ⟨σ_S | j⟩ ≈ log(N/m_j) / log(λ_j) + correction

If λ_j has structural closed form (e.g., from a per-attractor inverse-tree
residue-chain transition matrix), this bypasses the trajectory-measure
invariance requirement.

Method:
  For each m_j ∈ {5, 85, 341}, build inverse Syracuse tree by BFS:
    Start at m_j (depth 0). At each node n with n ≢ 0 mod 3, generate
    predecessors at v = v_min, v_min+2, ... where v_min = 1 if n ≡ 2 mod 3
    else 2. Predecessor m_pred = (n·2^v - 1)/3. Cap at m_pred ≤ N_max.
  Record count(d, j) = number of ancestors at depth d for attractor j.
  Compute ⟨σ_S | j⟩ = Σ d·count(d) / Σ count(d) (uniform weighting per ancestor).
  Compute λ_j = exp(slope of log count(d) vs d in growth regime).

Compare:
  Empirical ⟨σ_S | j⟩ via forward sim (Result 30) at same N
  Predicted ⟨σ_S | j⟩ via inverse-tree direct enumeration
  Spectral λ_j: does it close ⟨σ_S | j⟩ via μ_inv(j) = log λ_j?

Cap N_max = 10^7 (manageable; inverse tree of m_2 has ~5M nodes).
"""
import sys
import time
from pathlib import Path

import numpy as np
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2 * LOG2 - LOG3


def build_inverse_tree_count(m_j, N_max, max_depth=200):
    """BFS inverse Syracuse tree from m_j, capped at m_pred ≤ N_max.

    Returns:
        count_per_depth: dict {depth: int} of ancestor counts
        depth_logm_sum: dict {depth: float} of sum of log(m) at each depth
        depth_logm_sumsq: dict {depth: float} of sum of (log m)^2 at each depth
        total_count: total ancestors found
    """
    count_per_depth = defaultdict(int)
    logm_sum_per_depth = defaultdict(float)
    logm_sumsq_per_depth = defaultdict(float)

    # Start: depth 0 is m_j itself
    current_layer = [m_j]
    count_per_depth[0] = 1
    logm_sum_per_depth[0] = np.log(float(m_j))
    logm_sumsq_per_depth[0] = np.log(float(m_j))**2

    for depth in range(1, max_depth + 1):
        next_layer = []
        for n in current_layer:
            n_mod_3 = n % 3
            if n_mod_3 == 0:
                continue  # no Syracuse predecessors
            v = 1 if n_mod_3 == 2 else 2
            # Predecessors at v, v+2, v+4, ... while m_pred ≤ N_max
            n_2v = n << v  # n * 2^v
            while True:
                num = n_2v - 1
                if num <= 0:
                    break
                if num % 3 == 0:
                    m_pred = num // 3
                    if m_pred > N_max:
                        break
                    if m_pred > 1 and (m_pred & 1):  # odd, > 1
                        next_layer.append(m_pred)
                else:
                    # shouldn't happen if v_min was chosen correctly
                    pass
                v += 2
                n_2v <<= 2  # multiply by 4 for next v
                if n_2v > 3 * N_max + 1:
                    break
        if not next_layer:
            break
        count_per_depth[depth] = len(next_layer)
        logs = np.log(np.array(next_layer, dtype=np.float64))
        logm_sum_per_depth[depth] = float(logs.sum())
        logm_sumsq_per_depth[depth] = float((logs**2).sum())
        current_layer = next_layer

    total_count = sum(count_per_depth.values())
    return count_per_depth, logm_sum_per_depth, logm_sumsq_per_depth, total_count


def fit_growth_rate(count_per_depth, depth_min=2, depth_max=None):
    """Fit log count(d) = a + b·d for d in [depth_min, depth_max].
    Returns slope b = log(λ_j), λ_j = exp(b)."""
    depths = sorted(count_per_depth.keys())
    if depth_max is None:
        # Use the regime where count is growing (before saturation)
        depth_max = max(d for d in depths if count_per_depth[d] >= count_per_depth[max(depths)])
    pts = [(d, count_per_depth[d]) for d in depths if depth_min <= d <= depth_max and count_per_depth[d] > 0]
    if len(pts) < 3:
        return None
    d_arr = np.array([p[0] for p in pts])
    log_c = np.log([p[1] for p in pts])
    p = np.polyfit(d_arr, log_c, 1)
    slope, intercept = p
    pred = slope * d_arr + intercept
    R2 = 1 - np.sum((log_c - pred)**2) / np.sum((log_c - log_c.mean())**2)
    return {'slope': slope, 'lambda_j': np.exp(slope), 'R2': R2,
            'd_min': pts[0][0], 'd_max': pts[-1][0]}


def main():
    print("# Per-attractor inverse-tree spectral analysis")
    print("# Test: does ⟨σ_S | j⟩ ≈ log(N/m_j) / log(λ_j) close it?")
    print()

    # Use multiple N to check growth structure
    N_values = [10**6, 10**7]
    j_targets = [2, 4, 5, 7, 8]

    for N_max in N_values:
        print(f"## N_max = {N_max:,}")
        log_N = np.log(N_max)
        for jt in j_targets:
            m_j = (4**jt - 1) // 3
            log_mj = np.log(float(m_j))
            t0 = time.perf_counter()
            counts, logm_sum, logm_sumsq, total = build_inverse_tree_count(
                m_j, N_max, max_depth=300
            )
            elapsed = time.perf_counter() - t0
            if total < 100:
                print(f"  j={jt}, m_j={m_j}: only {total} ancestors found (skip)")
                continue
            # ⟨σ_S | j⟩ via direct enumeration (uniform weighting per ancestor)
            sigma_emp = sum(d * c for d, c in counts.items()) / total
            # ⟨log m | j⟩
            log_m_emp = sum(logm_sum[d] for d in logm_sum) / total
            # Variance
            log_m_sq = sum(logm_sumsq[d] for d in logm_sumsq) / total
            log_m_sd = np.sqrt(max(log_m_sq - log_m_emp**2, 0))
            # P(j) under uniform on [1, N]
            P_j = total / (N_max / 2)  # half are odd
            # W_j from this
            W_j_pred = sigma_emp - log_m_emp / LOG43 - 1 + log_mj / LOG43
            # Spectral fit
            growth_fit = fit_growth_rate(counts, depth_min=3)
            if growth_fit:
                # Spectral prediction: σ ≈ (log N - log m_j) / log(λ_j)  + correction
                # at saturation (when m_pred reaches N)
                spectral_sigma_pred = (log_N - log_mj) / growth_fit['slope'] if growth_fit['slope'] > 0 else None
                lam = growth_fit['lambda_j']
                R2 = growth_fit['R2']
                lam_str = f"λ_{jt}={lam:.5f}"
            else:
                spectral_sigma_pred = None
                lam_str = "fit failed"
            print(f"  j={jt} (m_j={m_j}, ancestors={total:,}, {elapsed:.1f}s):")
            print(f"     ⟨σ_S|j⟩_empirical_inverse_tree = {sigma_emp:.3f}")
            if spectral_sigma_pred is not None:
                print(f"     ⟨σ_S|j⟩_spectral = log(N/m_j)/log(λ_j) = {spectral_sigma_pred:.3f}  ({lam_str}, R²={R2:.4f})")
            print(f"     ⟨log m|j⟩ = {log_m_emp:.3f} (vs log N - 1 = {log_N - 1:.3f})")
            print(f"     SD log m | j = {log_m_sd:.3f}")
            print(f"     P(j) = {P_j:.4f}")
            print(f"     W_j_predicted = {W_j_pred:+.3f}")
        print()


if __name__ == "__main__":
    main()
