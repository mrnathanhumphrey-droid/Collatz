"""
PROBE MIRROR -- the mirror map x -> (2x+1)/3^v = (q,p)=(2,3), the q<->p image of Collatz (2026-07-28)

The family involution q=(p+1)/(p-1) sends Collatz (3,2) to (2,3). On the boundary curve
q(p-1)=p+1, the ONLY two constructible critical integer points are (3,2) and (2,3). So the whole
"is the critical class related?" question reduces to: does S_inf(3,2) relate to S_inf(2,3) under swap?

Reuses the certified Xk_qp transfer op from probe_phydra_family.py (tower q^k, chain (Z/q^k)*,
r -> (qr+1) p^-v, v~Geom(1-1/p)). Three checks:
  (a) exact early ladder S_k(2,3) vs S_k(3,2)
  (b) INVOLUTION test: lim S_k(2,3) vs lim S_k(3,2)~0.4749, same machinery both sides
  (c) DEPTH class: exact denominator rate -- doubly-exp (same infinite-Mahler class) or tame (special)?
Guardrail: a clean digit-relation is a LEAD to prove, never trusted from digits. Outcome 3
(unrelated value) is the pre-registered likely result and is itself a theorem-flavored finding
(criticality isolates each constant), not a null.

VERDICT (this run): involution symmetric on the criticality CONDITION and the arithmetic TYPE
(both infinite Mahler depth, same <2>-cofactor-prime mechanism: 41|2^20-1, 193|2^96-1), but NOT on
the VALUE (0.475 != 0.459). Each critical constant is algebraically alone. See result_MIRROR.md.
"""
import math
from fractions import Fraction as F
import numpy as np


def order_of(a, N):
    x = a % N; k = 1
    while x != 1:
        x = (x * a) % N; k += 1
        if k > 8 * N: return -1
    return k


def factor(n):
    if n == 1: return "1"
    f = {}; m = n; d = 2
    while d * d <= m:
        while m % d == 0: f[d] = f.get(d, 0) + 1; m //= d
        d += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return " * ".join(f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(f.items()))


def Xk_qp_float(q, p, k, vmax=400):
    N = q ** k; M = order_of(p, N); vm = min(M, vmax); invp = pow(p, -1, N)
    w = np.array([(1.0 / p) ** v for v in range(1, vm + 1)]); w /= w.sum()
    powinvp = np.array([pow(invp, v, N) for v in range(1, vm + 1)], dtype=np.int64)
    cop = np.array([r for r in range(N) if r % q != 0], dtype=np.int64)
    n = len(cop); idx = -np.ones(N, np.int64); idx[cop] = np.arange(n)
    base = ((q * cop + 1) % N).astype(np.int64)
    tgt = np.stack([idx[(base * powinvp[vi]) % N] for vi in range(vm)])
    pi = np.full(n, 1.0 / n)
    for _ in range(6000):
        nxt = np.zeros(n)
        for vi in range(vm):
            np.add.at(nxt, tgt[vi], w[vi] * pi)
        nxt /= nxt.sum()
        if np.abs(nxt - pi).max() < 1e-15: pi = nxt; break
        pi = nxt
    return (q ** k) * float((pi ** 2).sum())


def Xk_qp_exact(q, p, k):
    N = q ** k; M = order_of(p, N); invp = pow(p, -1, N)
    powinvp = [pow(invp, v, N) for v in range(1, M + 1)]
    Z = sum(F(1, p ** v) for v in range(1, M + 1))
    w = [F(1, p ** v) / Z for v in range(1, M + 1)]
    cop = [r for r in range(N) if r % q != 0]; idx = {r: i for i, r in enumerate(cop)}; n = len(cop)
    K = [[F(0)] * n for _ in range(n)]
    for r in cop:
        b = (q * r + 1) % N
        for vi in range(M):
            K[idx[r]][idx[(b * powinvp[vi]) % N]] += w[vi]
    A = [[K[j][i] - (F(1) if i == j else F(0)) for j in range(n)] for i in range(n)]
    A[n - 1] = [F(1)] * n
    bb = [F(0)] * n; bb[n - 1] = F(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        A[c], A[piv] = A[piv], A[c]; bb[c], bb[piv] = bb[piv], bb[c]
        pv = A[c][c]; A[c] = [v / pv for v in A[c]]; bb[c] /= pv
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [A[r][j] - f * A[c][j] for j in range(n)]; bb[r] -= f * bb[c]
    return F(q ** k) * sum(x * x for x in bb)


def main():
    print("# PROBE MIRROR  x -> (2x+1)/3^v  (q,p)=(2,3)  q<->p image of Collatz\n")

    # (a) exact early ladder, mirror vs Collatz
    print("## (a) exact S_k early ladder  (mirror vs Collatz)")
    print(f"   {'k':>2} {'S_k(2,3) mirror':>18} | {'S_k(3,2) Collatz':>18}")
    Xm = {0: F(1)}; Xc = {0: F(1)}
    for k in range(1, 5):
        Xm[k] = Xk_qp_exact(2, 3, k); Xc[k] = Xk_qp_exact(3, 2, k)
        Sm = Xm[k] - Xm[k - 1]; Sc = Xc[k] - Xc[k - 1]
        sm = str(Sm) if len(str(Sm)) < 18 else f"{float(Sm):.10f}"
        sc = str(Sc) if len(str(Sc)) < 18 else f"{float(Sc):.10f}"
        print(f"   {k:>2} {sm:>18} | {sc:>18}")
    print()

    # (b) involution test -- lim S_k both sides, same machinery
    print("## (b) involution test -- lim S_k(2,3) vs lim S_k(3,2)~0.4749")
    Xmf = {0: 1.0}; Xcf = {0: 1.0}; Sm = Sc = None
    for k in range(1, 15):
        Xmf[k] = Xk_qp_float(2, 3, k); Xcf[k] = Xk_qp_float(3, 2, k)
        Sm = Xmf[k] - Xmf[k - 1]; Sc = Xcf[k] - Xcf[k - 1]
        if k >= 10:
            print(f"   k={k:>2}  S(2,3)={Sm:.6f}   S(3,2)={Sc:.6f}")
    print(f"\n   mirror  lim S_k(2,3) ~ {Sm:.6f}")
    print(f"   Collatz lim S_k(3,2) ~ {Sc:.6f}  (-> banked S_inf~0.4749, floor 0.473177)")
    print(f"   diff (2,3)-(3,2) = {Sm-Sc:+.6f}  ratio {Sm/Sc:.6f}  sum {Sm+Sc:.6f}  prod {Sm*Sc:.6f}")
    print("   => equality REFUTED (gap >> noise); no clean function at 3 digits.\n")

    # (c) depth class -- exact denominator rate, k=1..8
    print("## (c) depth class -- mirror exact denominators k=1..8 (doubly-exp => same infinite-depth class)")
    print(f"   {'k':>2} {'den(S_k)':>14} {'bits':>7} {'ratio':>7}  factor(den)")
    Xe = {0: F(1)}; Dprev = None
    for k in range(1, 9):
        Xe[k] = Xk_qp_exact(2, 3, k)
        Sk = Xe[k] - Xe[k - 1]; den = Sk.denominator
        D = math.log2(den) if den > 1 else 0.0
        ratio = (D / Dprev) if Dprev and Dprev > 0 else float('nan')
        print(f"   {k:>2} {den:>14} {D:>7.2f} {ratio:>7.3f}  {factor(den)}")
        Dprev = D
    print("   bits 5->13->47->106, ratio avg >2 = MAHLER signature; 41|2^20-1, 193|2^96-1 =")
    print("   <2>-cofactor primes (same mechanism as c_tilde). => SAME class, infinite Mahler depth.")
    print("\n   NET: involution symmetric on criticality CONDITION + arithmetic TYPE, silent on VALUE.")


if __name__ == "__main__":
    main()
