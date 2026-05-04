"""
Per-a magnitude pattern for primitive Fourier coefficients on Z/3^k Z.

Steps:
  1. Build R66 Markov chain on (Z/3^k Z)* coprime-to-3 residues at k=2..6
  2. Compute |mu_hat(a/3^k)|^2 = |sum_r pi_r exp(2*pi*i*a*r/3^k)|^2 for all primitive a
  3. Group by conjugate pairs (a, 3^k-a)
  4. Examine multiplicative structure: order in (Z/3^k Z)*, power of generator g=2
  5. Test parametric forms: |mu_hat|^2 vs (a mod divisor of phi(3^k))
  6. Verify averaging: sum_a |mu_hat(a/3^k)|^2 -> 7/15
  7. Characterize asymptotic distribution shape via histograms

Outputs:
  experiments_output/87_per_a_values.csv
  experiments_output/87_conjugate_pair_summary.csv
  experiments_output/87_multiplicative_structure.csv
  experiments_output/87_per_a_pattern_log.txt
"""
import sys
import io
import math
import cmath
from pathlib import Path
from collections import defaultdict

import numpy as np
import polars as pl

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

OUT = Path("C:/Collatz")
EXP_OUT = Path("C:/Collatz/experiments_output")

results_log = []
def log(s):
    print(s, flush=True)
    results_log.append(s)


# ============================================================
# R66 Markov chain on coprime-to-3 residues mod 3^k
# Transition: T(r) = (3r+1) · 2^(-v) mod 3^k, v ~ Geom(1/2)
# ============================================================

def build_markov_chain(k):
    N = 3**k
    M = 2 * 3**(k - 1)  # ord_{3^k}(2)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime_states = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    n_states = len(coprime_states)

    K = np.zeros((n_states, n_states), dtype=np.float64)
    Z_v = 1.0 - 2.0 ** (-M)
    for r in coprime_states:
        for r_v in range(1, M + 1):
            p = 2.0 ** (-r_v) / Z_v
            target = ((3 * r + 1) * powers_inv2[r_v - 1]) % N
            K[state_idx[r], state_idx[target]] += p
    return K, coprime_states


def stationary(K):
    eigvals, eigvecs = np.linalg.eig(K.T)
    idx = np.argmax(np.real(eigvals))
    pi = np.real(eigvecs[:, idx])
    if pi.sum() < 0: pi = -pi
    pi = pi / pi.sum()
    return pi


def fourier_at(a, k, pi, states):
    N = 3 ** k
    s = 0.0 + 0.0j
    for r, p in zip(states, pi):
        s += p * cmath.exp(2j * math.pi * a * r / N)
    return s


# ============================================================
# Multiplicative structure helpers
# ============================================================

def find_generator(k):
    """Find primitive root mod 3^k. For 3^k with k>=1, g=2 is a generator."""
    N = 3**k
    order_target = 2 * 3**(k - 1)
    for g in [2, 5, 11]:
        if math.gcd(g, N) != 1: continue
        # Verify order
        order = 1
        x = g
        while x != 1:
            x = (x * g) % N
            order += 1
            if order > order_target: break
        if order == order_target:
            return g
    raise ValueError(f"No generator found for 3^{k}")


def discrete_log(a, g, k):
    """Discrete log: find j such that g^j ≡ a mod 3^k."""
    N = 3 ** k
    order = 2 * 3 ** (k - 1)
    x = 1
    for j in range(order):
        if x == a: return j
        x = (x * g) % N
    return -1  # not found


def order_of(a, k):
    """Order of a in (Z/3^k Z)*."""
    N = 3 ** k
    order_target = 2 * 3 ** (k - 1)
    x = a
    o = 1
    while x != 1:
        x = (x * a) % N
        o += 1
        if o > order_target: return -1
    return o


# ============================================================
# Main
# ============================================================

def main():
    log("=" * 80)
    log("PER-a MAGNITUDE PATTERN for |mu_hat(a/3^k)|^2 (primitive Fourier coefficients)")
    log("=" * 80)

    rows_per_a = []
    rows_conj_pair = []
    rows_mult = []

    K_LEVELS = [2, 3, 4, 5, 6]

    for k in K_LEVELS:
        N = 3 ** k
        order_units = 2 * 3 ** (k - 1)
        log(f"\n=== k = {k}, 3^k = {N}, |units| = phi(3^k) = {order_units} ===")

        K_mat, states = build_markov_chain(k)
        pi = stationary(K_mat)
        log(f"  Markov chain: {len(states)} coprime states")
        log(f"  Stationary: max pi = {pi.max():.6f}, min pi = {pi.min():.6f}")

        # All primitive a (gcd(a, 3^k) = 1)
        primitives = [a for a in range(1, N) if math.gcd(a, N) == 1]

        # Generator and discrete logs
        try:
            g = find_generator(k)
            log(f"  Generator g = {g}, order = {order_units}")
        except Exception as e:
            log(f"  Generator search failed: {e}")
            g = None

        # Compute |mu_hat|^2 per primitive a
        per_a = {}
        for a in primitives:
            mh = fourier_at(a, k, pi, states)
            per_a[a] = abs(mh) ** 2

        # Sum and average
        S_k = sum(per_a.values())
        avg = S_k / len(primitives)
        max_a = max(per_a, key=per_a.get)
        min_a = min(per_a, key=per_a.get)
        log(f"\n  S_k = sum_a |mu_hat|^2 = {S_k:.6f}  (target: -> 7/15 = {7/15:.6f})")
        log(f"  avg = S_k / |units| = {avg:.6f}  (target: (7/30)·3^(-(k-1)) = {(7/30) / 3**(k-1):.6f})")
        log(f"  max at a = {max_a}: |mu_hat|^2 = {per_a[max_a]:.6f}")
        log(f"  min at a = {min_a}: |mu_hat|^2 = {per_a[min_a]:.6f}")
        log(f"  max/avg ratio = {per_a[max_a] / avg:.4f}, min/avg = {per_a[min_a] / avg:.4f}")

        # Group by conjugate pairs (a, 3^k - a)
        pairs = []
        seen = set()
        for a in primitives:
            if a in seen: continue
            a_conj = N - a
            if a == a_conj:
                pairs.append((a, a))
            else:
                pairs.append((a, a_conj))
                seen.add(a_conj)
            seen.add(a)
        # Sort pairs by |mu_hat|^2 descending
        pairs.sort(key=lambda p: -per_a[p[0]])
        log(f"\n  Conjugate pairs (sorted by |mu_hat|^2 desc):")
        log(f"    {'pair':>10}  {'|mu_hat|^2':>12}  {'a/(2*3^(k-1))':>15}  {'a*g (mod N)':>11}")

        for pi_idx, (a, a_conj) in enumerate(pairs):
            v_pair = per_a[a]
            log_a = discrete_log(a, g, k) if g else -1
            ratio_to_avg = v_pair / avg
            if pi_idx < 10:
                log(f"    ({a:>3},{a_conj:>3})  {v_pair:>12.6f}  {ratio_to_avg:>15.4f}  {log_a:>11}")
            rows_conj_pair.append({
                'k': k, 'a': a, 'a_conj': a_conj, 'mu_hat_sq': v_pair,
                'discrete_log_a': log_a, 'ratio_to_avg': ratio_to_avg,
            })

        if len(pairs) > 10:
            log(f"    ... ({len(pairs) - 10} more pairs)")

        # Per-a record
        for a in primitives:
            log_a = discrete_log(a, g, k) if g else -1
            ord_a = order_of(a, k)
            rows_per_a.append({
                'k': k, 'a': a, 'mu_hat_sq': per_a[a],
                'discrete_log': log_a, 'order_in_units': ord_a,
                'a_mod_3': a % 3,
            })

        # Multiplicative structure analysis
        # Group by discrete-log mod small divisors of order_units
        order_units = 2 * 3 ** (k - 1)
        for div in [2, 3]:
            if order_units % div != 0: continue
            log_groups = defaultdict(list)
            for a in primitives:
                la = discrete_log(a, g, k)
                log_groups[la % div].append(per_a[a])
            log(f"\n  Group by discrete_log mod {div}:")
            for grp, vals in sorted(log_groups.items()):
                log(f"    log mod {div} = {grp}: count={len(vals)}, mean={np.mean(vals):.6f}, max={max(vals):.6f}, min={min(vals):.6f}")
                rows_mult.append({
                    'k': k, 'div': div, 'log_mod': grp,
                    'count': len(vals), 'mean': float(np.mean(vals)),
                    'max': float(max(vals)), 'min': float(min(vals)),
                })

        # Distribution shape
        vals_arr = np.array(list(per_a.values()))
        normalized = vals_arr * (2 * 3 ** (k - 1)) / (7 / 15)
        log(f"\n  Normalized q_a = |mu_hat|^2 · 2·3^(k-1) / (7/15):")
        log(f"    mean = {normalized.mean():.4f} (target: 1.0)")
        log(f"    std  = {normalized.std():.4f}")
        log(f"    min  = {normalized.min():.4f}")
        log(f"    max  = {normalized.max():.4f}")
        log(f"    25th = {np.percentile(normalized, 25):.4f}")
        log(f"    50th = {np.percentile(normalized, 50):.4f}")
        log(f"    75th = {np.percentile(normalized, 75):.4f}")

    pl.DataFrame(rows_per_a).write_csv(EXP_OUT / "87_per_a_values.csv")
    pl.DataFrame(rows_conj_pair).write_csv(EXP_OUT / "87_conjugate_pair_summary.csv")
    pl.DataFrame(rows_mult).write_csv(EXP_OUT / "87_multiplicative_structure.csv")

    # ============================================================
    # Look for closed-form pattern: does |mu_hat(a/3^k)|^2 depend on a only via a mod ?
    # ============================================================
    log("\n" + "=" * 80)
    log("CLOSED-FORM SEARCH: does |mu_hat(a/3^k)|^2 depend on a only via a mod ?")
    log("=" * 80)

    log("\n  Test: at k=3, group all 18 primitive a by various a mod m. Constant within group?\n")

    k_test = 3
    N = 3 ** k_test
    primitives = [a for a in range(1, N) if math.gcd(a, N) == 1]
    K_mat, states = build_markov_chain(k_test)
    pi = stationary(K_mat)
    per_a = {a: abs(fourier_at(a, k_test, pi, states))**2 for a in primitives}

    g = find_generator(k_test)
    log(f"  k={k_test}, generator g={g}, |units|={2 * 3**(k_test-1)} = 18\n")

    log(f"  All 18 primitive a, |mu_hat|^2, discrete log mod 18, log mod 9, log mod 6, log mod 3, log mod 2:")
    log(f"  {'a':>3}  {'|mu_hat|^2':>12}  {'log_a':>5}  {'mod9':>5}  {'mod6':>5}  {'mod3':>5}  {'mod2':>5}")
    rows_test = []
    for a in primitives:
        v = per_a[a]
        la = discrete_log(a, g, k_test)
        log(f"  {a:>3}  {v:>12.6f}  {la:>5}  {la%9:>5}  {la%6:>5}  {la%3:>5}  {la%2:>5}")
        rows_test.append({'a': a, 'val': v, 'log': la})

    # Test grouping
    for mod in [2, 3, 6, 9, 18]:
        groups = defaultdict(list)
        for r in rows_test:
            groups[r['log'] % mod].append(r['val'])
        # Check if values within each group are constant (small std)
        max_std = 0
        for grp, vals in groups.items():
            std = np.std(vals)
            if std > max_std: max_std = std
        log(f"  Group by log mod {mod}: max within-group std = {max_std:.6f}, max value variation = {max_std/np.mean(list(per_a.values())):.4f} relative")

    # Test: maybe pairs of (log mod something) determine value?
    # Closer look at conjugate pairs and cube/square classes
    log(f"\n  Cube and square classes:")
    for a in primitives:
        la = discrete_log(a, g, k_test)
        is_cube = (la % 3 == 0)
        is_square = (la % 2 == 0)
        log(f"  a={a:>3}: |mu|^2={per_a[a]:.4f}  is_cube={is_cube} is_square={is_square}")

    # ============================================================
    # Asymptotic distribution shape analysis
    # ============================================================
    log("\n" + "=" * 80)
    log("ASYMPTOTIC DISTRIBUTION SHAPE")
    log("=" * 80)

    log("\n  Normalized values q_a = |mu_hat(a/3^k)|^2 · 2·3^(k-1) / (7/15)")
    log("  q has mean 1 by construction; examine its distribution shape\n")

    log(f"  {'k':>3}  {'#primitive':>10}  {'mean':>7}  {'std':>7}  {'min':>7}  {'max':>7}  {'25th':>7}  {'50th':>7}  {'75th':>7}")
    distribution_rows = []
    for k in K_LEVELS:
        N = 3 ** k
        primitives = [a for a in range(1, N) if math.gcd(a, N) == 1]
        K_mat, states = build_markov_chain(k)
        pi = stationary(K_mat)
        vals_arr = np.array([abs(fourier_at(a, k, pi, states))**2 for a in primitives])
        normalized = vals_arr * (2 * 3 ** (k - 1)) / (7 / 15)
        log(f"  {k:>3}  {len(primitives):>10}  {normalized.mean():>7.4f}  {normalized.std():>7.4f}  {normalized.min():>7.4f}  {normalized.max():>7.4f}  {np.percentile(normalized, 25):>7.4f}  {np.percentile(normalized, 50):>7.4f}  {np.percentile(normalized, 75):>7.4f}")
        distribution_rows.append({
            'k': k, 'n_primitive': len(primitives),
            'mean': float(normalized.mean()), 'std': float(normalized.std()),
            'min': float(normalized.min()), 'max': float(normalized.max()),
            'p25': float(np.percentile(normalized, 25)),
            'p50': float(np.percentile(normalized, 50)),
            'p75': float(np.percentile(normalized, 75)),
        })

    pl.DataFrame(distribution_rows).write_csv(EXP_OUT / "87_distribution_shape.csv")

    # Test for limiting distribution shape: bimodal? exponential? uniform?
    log("\n  Distribution shape diagnosis:")
    log(f"  Higher k normalized values:")
    k = K_LEVELS[-1]
    N = 3 ** k
    primitives = [a for a in range(1, N) if math.gcd(a, N) == 1]
    K_mat, states = build_markov_chain(k)
    pi = stationary(K_mat)
    vals_arr = np.array([abs(fourier_at(a, k, pi, states))**2 for a in primitives])
    normalized = vals_arr * (2 * 3 ** (k - 1)) / (7 / 15)

    # Histogram bins
    log(f"\n  Histogram at k={k} (n_primitive = {len(primitives)}):")
    bins = np.linspace(0, normalized.max(), 11)
    counts, edges = np.histogram(normalized, bins=bins)
    log(f"  {'bin lo':>7}  {'bin hi':>7}  {'count':>7}  {'frac':>7}")
    for i in range(len(counts)):
        log(f"  {edges[i]:>7.3f}  {edges[i+1]:>7.3f}  {counts[i]:>7}  {counts[i]/len(normalized):>7.4f}")

    # Test for exponential distribution: mean = 1, P(q>x) = exp(-x)
    log(f"\n  Test exponential dist: P(q > 1) should be ~0.368 = e^-1")
    p_above_1 = (normalized > 1).mean()
    log(f"  Empirical P(q > 1) = {p_above_1:.4f}")
    log(f"  Test exponential: P(q > 2) should be ~0.135 = e^-2")
    p_above_2 = (normalized > 2).mean()
    log(f"  Empirical P(q > 2) = {p_above_2:.4f}")

    # Verdict
    log("\n" + "=" * 80)
    log("VERDICT")
    log("=" * 80)
    log("\n  See per_a_magnitude_pattern.md for full analysis.")

    (EXP_OUT / "87_per_a_pattern_log.txt").write_text(
        "\n".join(results_log), encoding="utf-8")


if __name__ == "__main__":
    main()
