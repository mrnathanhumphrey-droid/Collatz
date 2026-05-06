"""
Probe: c=7/45 perturbation stability test (methodology choice (a) — rounded perturbation).

For each tested k, build perturbed Markov chain K_k(eps) where the integer multiplier 3 is
replaced by the real (3+eps):

    base_real = (3 + eps) * r + 1
    base_int  = round_half_up(base_real)        # nearest integer, ties -> up
    base_mod  = base_int mod 3^k
    For each v in {1..M} with truncated-Geom(1/2) weight 2^(-v) / (1 - 2^(-M)):
        target = (base_mod * inv2^v) mod 3^k
        K[r, target] += weight

Definitions (matching c_seven_forty_fifth.md / result_q_sweep_test_2):
    X_k(eps) := 3^k * sum_r pi_k(eps; r)^2
    S_k(eps) := X_k(eps) - X_{k-1}(eps)
    eps_k(perturb)  := S_k(eps) - 7/15

So at every eps we build K at level k AND k-1 to form S_k.

Leak handling:
- "Full leak" (base_int ≡ 0 mod 3): row would have all targets ≡ 0 mod 3,
  i.e. all transitions land on non-coprime states. Treated as self-loop.
- Per-v leak (target ≡ 0 mod 3): rare; row renormalized to sum to 1.

eps grid: {±0.05, ±0.02, ±0.01, ±0.005, ±0.002, ±0.001, 0}. Spans below
rounding threshold (no perturbation) through "many states perturbed".
For k=5 max coprime r=242; threshold for any perturbation |eps| >= 0.5/242 ~ 0.00207.
"""

import csv
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = Path("C:/Collatz/probe_perturbation")
OUTDIR.mkdir(exist_ok=True)

S_INF = 7.0 / 15.0


def round_half_up(x):
    return int(math.floor(x + 0.5))


def order_of_2_mod(N):
    o, p = 1, 2 % N
    while p != 1:
        p = (p * 2) % N
        o += 1
    return o


def build_K_eps(k, eps):
    """Return (K, n, n_perturbed, n_full_leak, n_per_v_leak, base_shift_dist, M)."""
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
    n_perturbed = 0
    n_full_leak = 0
    n_per_v_leak = 0
    base_shift_dist = {}

    for i, r in enumerate(coprime):
        base_real = (3.0 + eps) * r + 1.0
        base_int = round_half_up(base_real)
        shift = base_int - (3 * r + 1)
        base_shift_dist[shift] = base_shift_dist.get(shift, 0) + 1
        if shift != 0:
            n_perturbed += 1
        if base_int % 3 == 0:
            n_full_leak += 1
            continue
        base_mod = base_int % N
        for v in range(1, M + 1):
            target = (base_mod * inv2_pow[v]) % N
            if target % 3 == 0:
                n_per_v_leak += 1
                continue
            j = idx[target]
            K[i, j] += weights[v]

    row_sums = K.sum(axis=1)
    for i in range(n):
        if row_sums[i] > 0:
            K[i, :] /= row_sums[i]
        else:
            K[i, i] = 1.0

    return K, n, n_perturbed, n_full_leak, n_per_v_leak, base_shift_dist, M


def power_iter(K, max_iter=50000, tol=1e-15):
    n = K.shape[0]
    pi = np.ones(n) / n
    KT = K.T.copy()
    diff = np.inf
    for it in range(max_iter):
        pi_new = KT @ pi
        s = pi_new.sum()
        if s <= 0:
            return pi, it, np.inf
        pi_new /= s
        diff = float(np.linalg.norm(pi_new - pi, ord=np.inf))
        if diff < tol:
            return pi_new, it + 1, diff
        pi = pi_new
    return pi, max_iter, diff


def stationary_X(k, eps):
    """Build K_k(eps), find pi_k, return X_k = 3^k * ||pi_k||^2 plus diagnostics."""
    N = 3 ** k
    K, n, n_pert, n_leak, n_pv_leak, shifts, M = build_K_eps(k, eps)
    pi, iters, residual = power_iter(K)
    X_k = N * float(np.sum(pi ** 2))
    return {
        "X_k": X_k, "pi": pi, "K": K, "n": n,
        "n_perturbed": n_pert, "n_full_leak": n_leak, "n_per_v_leak": n_pv_leak,
        "shifts": shifts, "M": M, "iters": iters, "residual": residual,
    }


def main():
    eps_list = [-0.05, -0.02, -0.01, -0.005, -0.002, -0.001,
                0.0,
                +0.001, +0.002, +0.005, +0.01, +0.02, +0.05]

    k_target = 5
    print(f"Building K at level k={k_target} and k={k_target - 1} for each eps")
    print(f"S_k(eps) = X_k(eps) - X_{{k-1}}(eps); 7/15 = {S_INF:.10f}")
    print()

    # Build at both k=5 and k=4 for every eps
    print(f"{'k':>2} {'eps':>9} {'X_k':>14} {'iters':>6} {'residual':>10} "
          f"{'pert/n':>10} {'leak':>5} {'pv-leak':>7} {'M':>4} {'elapsed_s':>9}")
    print("-" * 105)

    rows_kk = []     # one row per (k, eps)
    for k in [k_target - 1, k_target]:
        N = 3 ** k
        for eps in eps_list:
            t0 = time.time()
            res = stationary_X(k, eps)
            elapsed = time.time() - t0
            print(f"{k:>2} {eps:>+9.4f} {res['X_k']:>14.10f} {res['iters']:>6} "
                  f"{res['residual']:>10.2e} {res['n_perturbed']:>3}/{res['n']:>3}  "
                  f"{res['n_full_leak']:>5} {res['n_per_v_leak']:>7} {res['M']:>4} {elapsed:>9.2f}")
            rows_kk.append({
                "k": k, "eps": eps, "X_k": res["X_k"], "iters": res["iters"],
                "residual": res["residual"], "n": res["n"], "n_perturbed": res["n_perturbed"],
                "n_full_leak": res["n_full_leak"], "n_per_v_leak": res["n_per_v_leak"],
                "M": res["M"], "shifts": res["shifts"], "elapsed_s": elapsed,
            })

    # Form S_k(eps) = X_k(eps) - X_{k-1}(eps)
    X_by = {(r["k"], r["eps"]): r["X_k"] for r in rows_kk}
    pert_by = {(r["k"], r["eps"]): r["n_perturbed"] for r in rows_kk}
    leak_by = {(r["k"], r["eps"]): r["n_full_leak"] for r in rows_kk}

    print()
    print(f"{'k':>2} {'eps':>9} {'X_{k-1}':>14} {'X_k':>14} {'S_k':>14} {'eps_k':>13} "
          f"{'pert(k)':>8} {'pert(k-1)':>10}")
    print("-" * 100)

    summary_rows = []
    for eps in eps_list:
        if (k_target, eps) in X_by and (k_target - 1, eps) in X_by:
            X_km1 = X_by[(k_target - 1, eps)]
            X_k = X_by[(k_target, eps)]
            S_k = X_k - X_km1
            eps_k = S_k - S_INF
            pert_k = pert_by[(k_target, eps)]
            pert_km1 = pert_by[(k_target - 1, eps)]
            print(f"{k_target:>2} {eps:>+9.4f} {X_km1:>14.10f} {X_k:>14.10f} "
                  f"{S_k:>14.10f} {eps_k:>+13.6e} {pert_k:>8} {pert_km1:>10}")
            summary_rows.append({
                "k": k_target, "eps": eps, "X_km1": X_km1, "X_k": X_k,
                "S_k": S_k, "eps_k": eps_k,
                "n_perturbed_k": pert_k, "n_perturbed_km1": pert_km1,
                "n_full_leak_k": leak_by[(k_target, eps)],
                "n_full_leak_km1": leak_by[(k_target - 1, eps)],
            })

    with open(OUTDIR / "result_perturbation_S_k_curve.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "eps", "X_km1", "X_k", "S_k", "eps_k",
                    "n_perturbed_k", "n_perturbed_km1",
                    "n_full_leak_k", "n_full_leak_km1"])
        for r in summary_rows:
            w.writerow([r["k"], r["eps"], f"{r['X_km1']:.15e}", f"{r['X_k']:.15e}",
                        f"{r['S_k']:.15e}", f"{r['eps_k']:+.6e}",
                        r["n_perturbed_k"], r["n_perturbed_km1"],
                        r["n_full_leak_k"], r["n_full_leak_km1"]])

    print()
    print("Base-integer shift histogram per eps at k=5 (shift = base_int - (3r+1)):")
    for r in rows_kk:
        if r["k"] == k_target:
            shifts_str = ", ".join(f"{s:+d}:{c}" for s, c in sorted(r["shifts"].items()))
            print(f"  eps={r['eps']:+.4f}  {{{shifts_str}}}")

    print()
    print("Numerical derivative dS_k/d(eps) at eps=0 (centered, forward, backward):")
    by_eps_S = {r["eps"]: r["S_k"] for r in summary_rows}
    deriv_rows = []
    for h in [0.001, 0.002, 0.005, 0.01, 0.02, 0.05]:
        if h in by_eps_S and -h in by_eps_S and 0.0 in by_eps_S:
            dSde_c = (by_eps_S[h] - by_eps_S[-h]) / (2 * h)
            dSde_f = (by_eps_S[h] - by_eps_S[0.0]) / h
            dSde_b = (by_eps_S[0.0] - by_eps_S[-h]) / h
            print(f"  k={k_target}  h={h:.4f}  centered={dSde_c:+.6e}  "
                  f"forward={dSde_f:+.6e}  backward={dSde_b:+.6e}")
            deriv_rows.append([k_target, h, dSde_c, dSde_f, dSde_b])

    with open(OUTDIR / "result_perturbation_derivative.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "h", "centered_dS_de", "forward_dS_de", "backward_dS_de"])
        for r in deriv_rows:
            w.writerow([r[0], r[1], f"{r[2]:+.6e}", f"{r[3]:+.6e}", f"{r[4]:+.6e}"])

    canonical_S_5 = 0.46551492
    eps_zero_S = next(r["S_k"] for r in summary_rows if r["eps"] == 0.0)
    print()
    print(f"Sanity check: at eps=0, S_5 (probe) = {eps_zero_S:.10f}")
    print(f"                       canonical = {canonical_S_5}")
    print(f"                       diff      = {eps_zero_S - canonical_S_5:+.4e}")

    print()
    print(f"Outputs:")
    print(f"  {OUTDIR / 'result_perturbation_S_k_curve.csv'}")
    print(f"  {OUTDIR / 'result_perturbation_derivative.csv'}")


if __name__ == "__main__":
    main()
