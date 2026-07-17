"""
c_tilde_q7_k4_probe.py
======================
Push q=7 to k=4 to test whether c~_7 ≈ 0.78 deviation from (q-3)/q = 0.571 is
finite-k transient or structural.

Estimate: q=7 k=4 has N = 6·7^3 = 2058 states, M = 3·7^3 = 1029. Roughly 5×
heavier than q=11 k=3 (N=1210) which took 4.57 hr. So ~20-25 hr in background.

Loads q=7 k=1..3 from cache; computes X_7^4 fresh; updates cache.
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8")

CACHE = r"C:\Collatz\experiments_output\result_q_sweep_test_2_cache.json"


def order_of_two(N):
    m = 1; v = 2 % N
    while v != 1:
        v = (v * 2) % N; m += 1
    return m


def build_markov_q(q: int, k: int):
    N = q ** k
    M = order_of_two(N)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime = [r for r in range(N) if r % q != 0]
    state_idx = {r: i for i, r in enumerate(coprime)}
    n = len(coprime)
    K = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    Z_v = Fraction(2 ** M - 1, 2 ** M)
    for r in coprime:
        for r_v in range(1, M + 1):
            p = Fraction(1, 2 ** r_v) / Z_v
            target = ((q * r + 1) * powers_inv2[r_v - 1]) % N
            K[state_idx[r]][state_idx[target]] += p
    return K, coprime, M


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
                pivot = row; break
        if pivot == -1: raise ValueError("singular")
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            b[col], b[pivot] = b[pivot], b[col]
        piv = A[col][col]
        for j in range(col, n): A[col][j] /= piv
        b[col] /= piv
        for row in range(n):
            if row != col and A[row][col] != 0:
                f = A[row][col]
                for j in range(col, n): A[row][j] -= f * A[col][j]
                b[row] -= f * b[col]
    return b


def main():
    print("=" * 72)
    print("q=7 k=4 probe: testing if c~_7 anomaly is finite-k transient")
    print("=" * 72)

    # Load cache
    with open(CACHE) as fh:
        d = json.load(fh)
    Xs = {}
    for key, val in d.items():
        q, k = map(int, key.split(","))
        if q == 7:
            Xs[k] = Fraction(int(val["X_num"]), int(val["X_den"]))

    print(f"  Loaded q=7 cache: k={sorted(Xs.keys())}")
    for k in sorted(Xs.keys()):
        print(f"    X_{k} = {float(Xs[k]):.10f}")
    print()

    # Compute X_7^4
    print("  Computing X_7^4 (N=2058 states, M=1029, ~20-25 hr expected)...", flush=True)
    t0 = time.time()
    K, coprime, M = build_markov_q(7, 4)
    t_build = time.time() - t0
    print(f"    K built in {t_build:.1f}s, N={len(coprime)}, M={M}", flush=True)

    t1 = time.time()
    pi = stationary_rational(K)
    t_solve = time.time() - t1
    X_4 = Fraction(7 ** 4) * sum(p * p for p in pi)
    print(f"    Stationary solved in {t_solve:.1f}s", flush=True)
    print(f"    X_4 = {float(X_4):.10f}", flush=True)

    # Save to cache
    d["7,4"] = {
        "X_num": str(X_4.numerator),
        "X_den": str(X_4.denominator),
        "n_states": len(coprime),
    }
    with open(CACHE, "w") as fh:
        json.dump(d, fh, indent=2)
    print(f"  [cache updated: {CACHE}]")
    print()

    Xs[4] = X_4
    Xs[0] = Fraction(1)

    # Compute S_k and c~_q sequence
    print("=" * 72)
    print("Updated c~_7 sequence")
    print("=" * 72)
    print(f"  {'k':>3}  {'S_k':>15}  {'c~_7 = S_k/(7/3)^k':>25}")
    for k in [1, 2, 3, 4]:
        Sk = Xs[k] - Xs[k-1]
        cq = float(Sk / (Fraction(7, 3) ** k))
        print(f"  {k:>3}  {float(Sk):>15.6f}  {cq:>25.6f}")
    print()

    # Compare to (q-3)/q = 4/7
    pred = 4.0 / 7
    cq_4 = float((Xs[4] - Xs[3]) / (Fraction(7, 3) ** 4))
    delta = cq_4 - pred
    print(f"  (q-3)/q = 4/7 = {pred:.6f}")
    print(f"  c~_7 at k=4 = {cq_4:.6f}")
    print(f"  delta(q=7) at k=4 = {delta:+.4f}")
    print()
    print(f"  Comparison to k=3 delta: {0.7826 - pred:+.4f}")
    if abs(delta) < 0.05:
        print("  -> Anomaly is FINITE-K transient. q=7 settling toward (q-3)/q.")
    elif abs(delta - (0.7826 - pred)) < 0.02:
        print("  -> Anomaly STABLE at k=4. q=7 has structural deviation; non-finite-k.")
    else:
        print(f"  -> Anomaly partially decreasing. delta moved from "
              f"{0.7826 - pred:+.4f} (k=3) to {delta:+.4f} (k=4). Needs k=5.")


if __name__ == "__main__":
    main()
