"""
c_tilde_q17_probe.py
====================
Compute c~_17 at k=2 to distinguish:
  - finite-k transient at q=7 (predicts c~_17 ≈ 14/17 = 0.8235)
  - structural non-prim-root correction (predicts c~_17 ~ 1.02 like q=7)
  - other (richer structure)

q=17 has ord(2 mod 17) = 8, NOT primitive root (q-1 = 16). Same class as q=7.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8")


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
    print("q=17 probe: c~_17 at k=2")
    print("=" * 72)
    print()

    # Sanity: ord(2 mod 17) and ord(2 mod 17^2)
    print(f"  ord(2 mod 17)   = {order_of_two(17)}     (expected 8, NOT prim root)")
    print(f"  ord(2 mod 17^2) = {order_of_two(17**2)}   (expected 136 = 8·17)")
    print()

    Xs = {0: Fraction(1)}

    for k in [1, 2]:
        t0 = time.time()
        K, coprime, M = build_markov_q(17, k)
        t_build = time.time() - t0
        n = len(coprime)
        t1 = time.time()
        pi = stationary_rational(K)
        t_solve = time.time() - t1
        Xk = Fraction(17 ** k) * sum(p * p for p in pi)
        Xs[k] = Xk
        print(f"  k={k}: build {t_build:.2f}s  solve {t_solve:.2f}s  N={n}  M={M}")
        print(f"    X_{k} = {Xk}")
        print(f"    X_{k} ~ {float(Xk):.10f}")
        print()

    S1 = Xs[1] - Xs[0]
    S2 = Xs[2] - Xs[1]
    print(f"  S_1^(17) = {S1}  ~ {float(S1):.10f}")
    print(f"  S_2^(17) = {S2}  ~ {float(S2):.10f}")
    print()

    c1 = float(S1 / (Fraction(17, 3) ** 1))
    c2 = float(S2 / (Fraction(17, 3) ** 2))
    print(f"  c~_17 sequence (k=1, 2): {c1:.6f}, {c2:.6f}")
    print(f"  S_{2}/S_{1} ratio = {float(S2/S1):.4f}  (expected 17/3 = {17/3:.4f})")
    print()

    # Predictions
    print("=" * 72)
    print("Decision rule")
    print("=" * 72)
    pred_main = (17 - 3) / 17
    pred_q7_like = pred_main + 0.21
    print(f"  (q-3)/q hypothesis: c~_17 ≈ {pred_main:.4f} = 14/17")
    print(f"  q=7-like deviation: c~_17 ≈ {pred_q7_like:.4f} (= 14/17 + 0.21)")
    print()
    print(f"  Empirical c~_17 (k=2) = {c2:.6f}")
    print(f"    deviation from 14/17 = {c2 - pred_main:+.4f}")
    print(f"    delta(q=17) = {c2 - pred_main:+.4f}")
    print(f"    compare delta(q=7) = +0.2112, delta(q=11) = +0.0015, delta(q=13) = +0.0006")
    print()

    if abs(c2 - pred_main) < 0.05:
        print("  -> CLOSE TO 14/17. q=7's deviation likely finite-k transient at q=7.")
        print("     (q-3)/q hypothesis SUPPORTED at q=17.")
    elif abs(c2 - pred_q7_like) < 0.05:
        print("  -> CLOSE TO q=7-like deviation. Non-prim-root pattern is structural.")
        print("     Need to characterize the non-prim-root correction.")
    else:
        print(f"  -> NEITHER. c~_17 = {c2:.4f} is between predictions.")
        print("     Suggests richer structure than current hypothesis allows.")


if __name__ == "__main__":
    main()
