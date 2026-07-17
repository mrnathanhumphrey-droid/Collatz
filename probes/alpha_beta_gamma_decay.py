"""
alpha_beta_gamma_decay.py — Compute (α, β, γ) sub-cell mass-share decay structure
for level-lifting transitions k → k+1 in the Markov chain on (Z/3^k Z)*.

Each residue r mod 3^k (coprime to 3) lifts to 3 residues mod 3^(k+1):
  r, r + 3^k, r + 2·3^k

Mass shares: (α̃_r, β̃_r, γ̃_r) = (π_{k+1}[r], π_{k+1}[r+3^k], π_{k+1}[r+2·3^k]) / π_k[r]
Deviation from uniform: d_r = (α̃_r - 1/3, β̃_r - 1/3, γ̃_r - 1/3)

Test decay: max |d_r(k)| ~ (1/2)^k? Other rate?
"""
import sys
import math
import json
import os
from fractions import Fraction
from collections import defaultdict
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"


def build_markov_rational(k):
    """Build K_k as matrix of Fractions on coprime-to-3 residues mod 3^k."""
    N = 3**k
    M = 2 * 3**(k-1)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime_states = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    n = len(coprime_states)
    K = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    Z_v = Fraction(2**M - 1, 2**M)
    for r in coprime_states:
        for r_v in range(1, M + 1):
            p = Fraction(1, 2**r_v) / Z_v
            target = ((3 * r + 1) * powers_inv2[r_v - 1]) % N
            K[state_idx[r]][state_idx[target]] += p
    return K, coprime_states


def stationary_rational(K):
    n = len(K)
    A = [[K[j][i] - (Fraction(1) if i == j else Fraction(0)) for j in range(n)] for i in range(n)]
    A[n - 1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n - 1] = Fraction(1)
    for col in range(n):
        pivot = -1
        for row in range(col, n):
            if A[row][col] != 0:
                pivot = row
                break
        if pivot == -1:
            raise ValueError(f"Singular at col {col}")
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            b[col], b[pivot] = b[pivot], b[col]
        piv = A[col][col]
        for j in range(col, n):
            A[col][j] /= piv
        b[col] /= piv
        for row in range(n):
            if row != col and A[row][col] != 0:
                factor = A[row][col]
                for j in range(col, n):
                    A[row][j] -= factor * A[col][j]
                b[row] -= factor * b[col]
    return b


def stationary_numerical(k):
    """Faster stationary computation via numpy power iteration for higher k."""
    N = 3**k
    M = 2 * 3**(k-1)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime_states = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    n = len(coprime_states)
    K = np.zeros((n, n))
    Z_v = 1.0 - 2.0**(-M)
    for r in coprime_states:
        for r_v in range(1, M + 1):
            p = 2.0**(-r_v) / Z_v
            target = ((3 * r + 1) * powers_inv2[r_v - 1]) % N
            K[state_idx[r], state_idx[target]] += p

    # Power iteration on K^T (left eigvec)
    pi = np.ones(n) / n
    for _ in range(100):
        pi = pi @ K
        pi = pi / pi.sum()
    return pi, coprime_states


def compute_alpha_beta_gamma(pi_k, coprime_k, pi_kp1, coprime_kp1, k):
    """For each r in coprime_k, compute (α̃, β̃, γ̃) over 3 lifts to level k+1."""
    N_k = 3**k
    state_idx_k = {r: i for i, r in enumerate(coprime_k)}
    state_idx_kp1 = {r: i for i, r in enumerate(coprime_kp1)}

    results = []  # list of (r, α, β, γ, α̃, β̃, γ̃, d_max, d_norm)
    for r in coprime_k:
        # 3 lifts: r, r + 3^k, r + 2*3^k
        lifts = [r, r + N_k, r + 2 * N_k]
        masses = []
        for lift in lifts:
            if lift in state_idx_kp1:
                masses.append(pi_kp1[state_idx_kp1[lift]])
            else:
                masses.append(0.0)
        total = sum(masses)
        if total == 0:
            continue
        # Should equal pi_k[r] by consistency
        pi_k_r = pi_k[state_idx_k[r]]
        ratio = total / pi_k_r if pi_k_r != 0 else 1.0
        if abs(float(ratio) - 1.0) > 0.01:
            print(f"  WARN: consistency check fail at k={k}, r={r}: total={total} vs pi_k[r]={pi_k_r}")

        alpha_t, beta_t, gamma_t = (m / total for m in masses)
        d = (alpha_t - 1/3, beta_t - 1/3, gamma_t - 1/3)
        d_max = max(abs(x) for x in d)
        d_norm = math.sqrt(sum(x**2 for x in d))

        results.append({
            'r': r,
            'pi_k_r': float(pi_k_r),
            'masses': [float(m) for m in masses],
            'shares': [alpha_t, beta_t, gamma_t],
            'deviation': list(d),
            'd_max': d_max,
            'd_norm': d_norm,
        })
    return results


def main():
    print("# (α, β, γ) sub-cell mass-share decay across k=1→5")
    print()

    # Compute stationary distributions
    pis = {}
    for k in [1, 2, 3, 4, 5, 6]:
        if k <= 3:
            # Exact rational
            K, coprime = build_markov_rational(k)
            pi_q = stationary_rational(K)
            pi = [float(p) for p in pi_q]
            pis[k] = (pi, coprime, pi_q)
        else:
            # Numerical
            pi, coprime = stationary_numerical(k)
            pis[k] = (list(pi), coprime, None)
        print(f"  k={k}: {len(pis[k][1])} states, stationary computed")

    # For each k → k+1 transition, compute (α, β, γ) per r
    print()
    print("# Sub-cell decomposition per level-lifting:")
    summary_rows = []
    all_data = {}
    for k in [1, 2, 3, 4, 5]:
        pi_k, coprime_k, _ = pis[k]
        pi_kp1, coprime_kp1, _ = pis[k+1]
        results = compute_alpha_beta_gamma(pi_k, coprime_k, pi_kp1, coprime_kp1, k)
        if not results:
            continue
        d_maxes = [r['d_max'] for r in results]
        d_norms = [r['d_norm'] for r in results]
        # ψ = α² + β² + γ² for each r
        psis = [sum(s**2 for s in r['shares']) for r in results]
        psi_min = min(psis); psi_max = max(psis); psi_avg = sum(psis)/len(psis)
        summary_rows.append({
            'k': k,
            'n_residues': len(results),
            'd_max_max': max(d_maxes),
            'd_max_avg': sum(d_maxes)/len(d_maxes),
            'd_norm_max': max(d_norms),
            'd_norm_avg': sum(d_norms)/len(d_norms),
            'psi_max': psi_max,
            'psi_min': psi_min,
            'psi_avg': psi_avg,
            'psi_uniform': 1/3,
        })
        all_data[k] = results

    # Print summary
    print()
    print(f"  {'k':>3}  {'#r':>4}  {'max |d|':>10}  {'avg |d|':>10}  {'max d-norm':>12}  {'ψ avg':>8}  {'ψ - 1/3 avg':>12}")
    for s in summary_rows:
        print(f"  {s['k']:>3}  {s['n_residues']:>4}  {s['d_max_max']:>10.6f}  {s['d_max_avg']:>10.6f}  {s['d_norm_max']:>12.6f}  {s['psi_avg']:>8.6f}  {s['psi_avg']-1/3:>12.6f}")

    # Geometric decay test
    print()
    print(f"# Geometric decay test (max |d_r| ratio per level)")
    print(f"  {'k → k+1':>8}  {'max|d|_k':>10}  {'max|d|_{k+1}':>14}  {'ratio':>8}")
    for i in range(1, len(summary_rows)):
        prev = summary_rows[i-1]
        curr = summary_rows[i]
        ratio = curr['d_max_max'] / prev['d_max_max']
        print(f"  {prev['k']}→{curr['k']:<3}  {prev['d_max_max']:>10.6f}  {curr['d_max_max']:>14.6f}  {ratio:>8.4f}")

    print()
    print(f"# log(max|d|) vs k linear fit:")
    if len(summary_rows) >= 3:
        ks = [s['k'] for s in summary_rows]
        logs = [math.log(s['d_max_max']) for s in summary_rows]
        slope, intercept = np.polyfit(ks, logs, 1)
        rate = math.exp(slope)
        print(f"  slope = {slope:.4f}, rate = exp(slope) = {rate:.4f}")
        print(f"  Conjectured rate 1/2: {1/2:.4f}, slope = log(1/2) = {math.log(1/2):.4f}")
        print(f"  Match? {'YES' if abs(rate - 0.5) < 0.05 else 'NO, actual rate ' + f'{rate:.4f}'}")

    # ψ test (R70 mechanism: ψ → 1/3 means uniform split, gives S∞ = 7/15)
    print()
    print(f"# ψ_avg vs k convergence to 1/3:")
    print(f"  {'k':>3}  {'ψ_avg':>10}  {'ψ_avg - 1/3':>12}  {'ratio_(k+1)/k':>14}")
    prev_dev = None
    for s in summary_rows:
        dev = s['psi_avg'] - 1/3
        ratio = dev / prev_dev if prev_dev else float('nan')
        print(f"  {s['k']:>3}  {s['psi_avg']:>10.6f}  {dev:>12.6f}  {ratio if prev_dev else 'n/a':>14}")
        prev_dev = dev

    # ψ = 3/7 at k=1 should be exact
    print()
    if summary_rows:
        first = summary_rows[0]
        print(f"  ψ at k=1 (lifting from k=1 to k=2): {first['psi_avg']:.6f}")
        print(f"  Expected from R70: 3/7 = {3/7:.6f}")
        print(f"  Match: {'YES (exact)' if abs(first['psi_avg'] - 3/7) < 1e-10 else 'NO'}")

    # Save data
    out_csv = os.path.join(OUTDIR, "alpha_beta_gamma_values.csv")
    with open(out_csv, 'w') as f:
        f.write("k,r,pi_k_r,alpha,beta,gamma,alpha_share,beta_share,gamma_share,d_alpha,d_beta,d_gamma,d_max,d_norm,psi\n")
        for k, results in all_data.items():
            for res in results:
                shares = res['shares']
                d = res['deviation']
                psi = sum(s**2 for s in shares)
                f.write(f"{k},{res['r']},{res['pi_k_r']:.10f},{res['masses'][0]:.10f},{res['masses'][1]:.10f},{res['masses'][2]:.10f},"
                        f"{shares[0]:.10f},{shares[1]:.10f},{shares[2]:.10f},"
                        f"{d[0]:.10f},{d[1]:.10f},{d[2]:.10f},{res['d_max']:.10f},{res['d_norm']:.10f},{psi:.10f}\n")
    print(f"\n[save] {out_csv}")

    out_summary = os.path.join(OUTDIR, "deviation_decay.csv")
    with open(out_summary, 'w') as f:
        f.write("k,n_residues,d_max_max,d_max_avg,d_norm_max,d_norm_avg,psi_avg,psi_min,psi_max\n")
        for s in summary_rows:
            f.write(f"{s['k']},{s['n_residues']},{s['d_max_max']:.10f},{s['d_max_avg']:.10f},{s['d_norm_max']:.10f},{s['d_norm_avg']:.10f},{s['psi_avg']:.10f},{s['psi_min']:.10f},{s['psi_max']:.10f}\n")
    print(f"[save] {out_summary}")

    # Hierarchy: do (α, β, γ) at k cluster by mod-3^j class for j < k?
    print()
    print(f"# Sub-cell pattern: do similar (α̃, β̃, γ̃) cluster by lower-mod class?")
    for k in [2, 3, 4]:
        if k not in all_data:
            continue
        results = all_data[k]
        # Group by r mod 3 (most basic)
        by_class = defaultdict(list)
        for res in results:
            cls = res['r'] % 3
            by_class[cls].append(res)
        print(f"\n  k={k}, grouped by r mod 3:")
        for cls in sorted(by_class.keys()):
            items = by_class[cls]
            d_maxes = [it['d_max'] for it in items]
            psis = [sum(s**2 for s in it['shares']) for it in items]
            print(f"    class {cls} mod 3: n={len(items)}, max|d|={max(d_maxes):.4f}, avg|d|={sum(d_maxes)/len(d_maxes):.4f}, ψ_avg={sum(psis)/len(psis):.4f}")


if __name__ == "__main__":
    main()
