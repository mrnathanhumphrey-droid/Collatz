"""
Probe: B-smooth / B-rough partition of (Z/3^k)* and Plancherel mass localization.

For each r in (Z/3^k)*, take the smallest positive lift n = r and run the forward
Collatz orbit (n -> n/2 if even, 3n+1 if odd) until n=1 or max_steps, capped at
max_value. Track the maximum prime factor across all integers in the orbit.

Classify r as B-smooth if max_prime <= B (and orbit reached 1), else B-rough.

For each (k, B):
  X_k         = 3^k * sum_r pi_k(r)^2
  X_k_smooth  = 3^k * sum_{r smooth} pi_k(r)^2
  X_k_rough   = X_k - X_k_smooth
  count_share_smooth = |smooth| / |coprime|
  mass_share_smooth  = X_k_smooth / X_k
  mass_ratio_smooth  = mass_share_smooth / count_share_smooth
                     (= 1 means pi^2 mass is uniformly distributed across smoothness;
                      != 1 means localization)

Tested k = 5, 6, 7 with B in {7, 50, 100, 1000, 10000}.
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np
from sympy import primefactors

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = Path("C:/Collatz/probe_smoothness")
OUTDIR.mkdir(exist_ok=True)

S_INF = 7.0 / 15.0


def order_of_2_mod(N):
    o, p = 1, 2 % N
    while p != 1:
        p = (p * 2) % N
        o += 1
    return o


def build_K(k):
    """Unperturbed K_k builder (truncated-Geom(1/2) over v in {1..M=ord_{3^k}(2)})."""
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
    return K, coprime, idx, M


def power_iter(K, max_iter=50000, tol=1e-15):
    n = K.shape[0]
    pi = np.ones(n) / n
    KT = K.T.copy()
    diff = np.inf
    for it in range(max_iter):
        pi_new = KT @ pi
        s = pi_new.sum()
        pi_new /= s
        diff = float(np.linalg.norm(pi_new - pi, ord=np.inf))
        if diff < tol:
            return pi_new, it + 1, diff
        pi = pi_new
    return pi, max_iter, diff


def collatz_orbit_max_prime(start, max_steps=200, max_value=10 ** 12):
    """Run forward Collatz from `start` until n=1 or max_steps or n > max_value.
    Returns (max_prime, length, reached_1, diverged)."""
    n = start
    max_p = 1
    if n > 1:
        ps = primefactors(n)
        if ps:
            max_p = max(ps)
    length = 0
    reached_1 = (n == 1)
    diverged = False
    for step in range(max_steps):
        if n == 1:
            reached_1 = True
            break
        if n % 2 == 0:
            n = n // 2
        else:
            n = 3 * n + 1
            if n > max_value:
                diverged = True
                break
        length += 1
        if n > 1:
            ps = primefactors(n)
            if ps:
                p = max(ps)
                if p > max_p:
                    max_p = p
    if n == 1:
        reached_1 = True
    return max_p, length, reached_1, diverged


def classify_residues(coprime, max_steps=200):
    """For each r in coprime, run orbit from start=r and tabulate max prime."""
    out = []
    for r in coprime:
        max_p, length, reached_1, diverged = collatz_orbit_max_prime(r, max_steps=max_steps)
        out.append({
            "r": r, "max_prime": max_p, "length": length,
            "reached_1": reached_1, "diverged": diverged,
        })
    return out


def main():
    B_list = [7, 50, 100, 1000, 10000]
    k_list = [5, 6, 7]

    # Per-residue partition rows: (k, r, max_prime, length, reached_1, diverged, pi_value)
    partition_rows = []
    summary_rows = []

    print(f"{'k':>3} {'B':>6} {'n_smooth':>9} {'n_rough':>8} {'n_total':>8} "
          f"{'count_share_sm':>14} {'mass_share_sm':>13} {'mass_ratio':>10} "
          f"{'X_k':>9} {'X_smooth':>9} {'X_rough':>9}")
    print("-" * 110)

    for k in k_list:
        N = 3 ** k
        t0 = time.time()
        K, coprime, idx, M = build_K(k)
        pi, iters, residual = power_iter(K)
        X_k = N * float(np.sum(pi ** 2))
        elapsed_pi = time.time() - t0

        t0 = time.time()
        residue_data = classify_residues(coprime, max_steps=200)
        elapsed_orbit = time.time() - t0

        max_primes = [d["max_prime"] for d in residue_data]
        diverged = sum(1 for d in residue_data if d["diverged"])
        not_reached = sum(1 for d in residue_data if not d["reached_1"])
        max_orbit_len = max(d["length"] for d in residue_data)

        print(f"# k={k}: pi in {elapsed_pi:.1f}s  iters={iters}  "
              f"residual={residual:.1e}  X_k={X_k:.6f}  M=ord_{{3^k}}(2)={M}")
        print(f"# k={k}: orbits in {elapsed_orbit:.1f}s  "
              f"max_p range [{min(max_primes)}, {max(max_primes)}]  "
              f"diverged={diverged}  not_reached_1={not_reached}  "
              f"max_orbit_len={max_orbit_len}")

        for d in residue_data:
            partition_rows.append({
                "k": k, **d, "pi": float(pi[idx[d["r"]]]),
            })

        for B in B_list:
            smooth_set = set()
            rough_set = set()
            for d in residue_data:
                if d["diverged"] or not d["reached_1"]:
                    rough_set.add(d["r"])
                elif d["max_prime"] <= B:
                    smooth_set.add(d["r"])
                else:
                    rough_set.add(d["r"])
            n_smooth = len(smooth_set)
            n_rough = len(rough_set)
            n_total = len(coprime)
            X_smooth = N * float(sum(pi[idx[r]] ** 2 for r in smooth_set))
            X_rough = N * float(sum(pi[idx[r]] ** 2 for r in rough_set))

            count_share = n_smooth / n_total if n_total else 0.0
            mass_share = X_smooth / X_k if X_k > 0 else 0.0
            mass_ratio = mass_share / count_share if count_share > 0 else float('nan')

            print(f"{k:>3} {B:>6} {n_smooth:>9} {n_rough:>8} {n_total:>8} "
                  f"{count_share:>14.4f} {mass_share:>13.4f} {mass_ratio:>10.4f} "
                  f"{X_k:>9.4f} {X_smooth:>9.4f} {X_rough:>9.4f}")
            summary_rows.append({
                "k": k, "B": B,
                "n_smooth": n_smooth, "n_rough": n_rough, "n_total": n_total,
                "X_k": X_k, "X_smooth": X_smooth, "X_rough": X_rough,
                "count_share_smooth": count_share,
                "mass_share_smooth": mass_share,
                "mass_ratio_smooth": mass_ratio,
            })

    # Save partition CSV (per-residue per-k)
    with open(OUTDIR / "result_smooth_rough_partition.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "r", "max_prime", "orbit_length", "reached_1", "diverged", "pi"])
        for r in partition_rows:
            w.writerow([r["k"], r["r"], r["max_prime"], r["length"],
                        int(r["reached_1"]), int(r["diverged"]), f"{r['pi']:.10e}"])

    # Save summary CSV
    with open(OUTDIR / "result_S_k_conditional.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "B", "n_smooth", "n_rough", "n_total",
                    "X_k", "X_smooth", "X_rough",
                    "count_share_smooth", "mass_share_smooth", "mass_ratio_smooth"])
        for r in summary_rows:
            w.writerow([r["k"], r["B"], r["n_smooth"], r["n_rough"], r["n_total"],
                        f"{r['X_k']:.10e}", f"{r['X_smooth']:.10e}", f"{r['X_rough']:.10e}",
                        f"{r['count_share_smooth']:.6f}", f"{r['mass_share_smooth']:.6f}",
                        f"{r['mass_ratio_smooth']:.6f}"])

    # Cross-level eps_k attempt: for each B, compute S_k_smooth = X_k_smooth - X_{k-1}_smooth
    # using each level's own classification. Note: this is empirical, not a Plancherel
    # restriction; partitions at different k's are formed from different lifts.
    eps_rows = []
    print()
    print("Cross-level S_k_smooth = X_k_smooth - X_{k-1}_smooth (empirical, partitions independent across k):")
    print(f"{'k':>3} {'B':>6} {'X_k_smooth':>11} {'X_km1_smooth':>13} {'S_k_smooth':>11} "
          f"{'X_k_rough':>11} {'X_km1_rough':>12} {'S_k_rough':>11} "
          f"{'S_k_total':>11} {'eps_k':>11}")
    print("-" * 130)
    by_kB = {(r["k"], r["B"]): r for r in summary_rows}
    for B in B_list:
        for k in [6, 7]:
            if (k, B) in by_kB and (k - 1, B) in by_kB:
                Xk_s = by_kB[(k, B)]["X_smooth"]
                Xkm1_s = by_kB[(k - 1, B)]["X_smooth"]
                S_s = Xk_s - Xkm1_s
                Xk_r = by_kB[(k, B)]["X_rough"]
                Xkm1_r = by_kB[(k - 1, B)]["X_rough"]
                S_r = Xk_r - Xkm1_r
                S_total = S_s + S_r
                eps_k = S_total - S_INF
                print(f"{k:>3} {B:>6} {Xk_s:>11.6f} {Xkm1_s:>13.6f} {S_s:>11.6f} "
                      f"{Xk_r:>11.6f} {Xkm1_r:>12.6f} {S_r:>11.6f} "
                      f"{S_total:>11.6f} {eps_k:>+11.4e}")
                eps_rows.append({
                    "k": k, "B": B,
                    "S_k_smooth": S_s, "S_k_rough": S_r, "S_k_total": S_total,
                    "eps_k": eps_k,
                })

    with open(OUTDIR / "result_eps_conditional.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "B", "S_k_smooth", "S_k_rough", "S_k_total", "eps_k"])
        for r in eps_rows:
            w.writerow([r["k"], r["B"],
                        f"{r['S_k_smooth']:.10e}", f"{r['S_k_rough']:.10e}",
                        f"{r['S_k_total']:.10e}", f"{r['eps_k']:+.6e}"])

    print()
    print("Outputs:")
    print(f"  {OUTDIR / 'result_smooth_rough_partition.csv'}")
    print(f"  {OUTDIR / 'result_S_k_conditional.csv'}")
    print(f"  {OUTDIR / 'result_eps_conditional.csv'}")


if __name__ == "__main__":
    main()
