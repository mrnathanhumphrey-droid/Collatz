"""
PROBE P-HYDRA FAMILY -- the (3,p)-Hydra Plancherel constants: same class? (2026-07-28)

The genuine infinite-depth siblings of our Collatz S_inf vary the VALUATION PRIME p (Siegel p-Hydra:
x -> (3x+1)/p^v), NOT q (q=3 is the only contracting integer). Chain on (Z/3^k)*: r -> (3r+1) p^-v,
v ~ Geom(1-1/p). p=2 is Collatz; p=3 degenerates (v_3(3r+1)=0); p=5,7,11,... are the siblings.

Wilson's Q2 asks for a PSLQ independence screen on the LIMIT values S_inf^(3,p). But we have only ~3
digits of S_inf^(3,2) (infinite Mahler order => no acceleration; rate 0.867 => ~16 levels/digit), so
PSLQ-on-limits is NOT honestly runnable. Instead run the STRUCTURAL form of Q1/Q2, exact on finite-level
rationals: do the siblings share the infinite-depth signature (doubly-exp denominator rate)?

  A  build the family, gate p=2 vs known S_k (2/3, 10/21, ...); float convergence + reachable precision.
  B  exact S_k^(3,p) denominator rate: doubly-exp (ratio ~ related to p) => infinite-depth => SAME CLASS,
     vs single-exp/bounded => finite-depth => DIFFERENT class. (The MAHLER test, per sibling.)
  C  verdict on same-class + honest PSLQ feasibility.
"""
import os, sys, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np


def order_of(a, N):
    x = a % N; k = 1
    while x != 1:
        x = (x * a) % N; k += 1
        if k > 4 * N: return -1
    return k


# ---------- float power-iteration (fast, high k, for convergence/precision) ----------
def Xk_float(p, k, vmax=300):
    N = 3 ** k; M = order_of(p, N); vm = min(M, vmax)
    invp = pow(p, -1, N)
    w = np.array([(1.0 / p) ** v for v in range(1, vm + 1)]); w /= w.sum()
    powinvp = np.array([pow(invp, v, N) for v in range(1, vm + 1)], dtype=np.int64)
    cop = np.array([r for r in range(N) if r % 3 != 0], dtype=np.int64)
    n = len(cop); idx = -np.ones(N, np.int64); idx[cop] = np.arange(n)
    base = ((3 * cop + 1) % N).astype(np.int64)
    tgt = np.empty((vm, n), np.int64)
    for vi in range(vm):
        tgt[vi] = idx[(base * powinvp[vi]) % N]
    pi = np.full(n, 1.0 / n)
    for _ in range(4000):
        nxt = np.zeros(n)
        for vi in range(vm):
            np.add.at(nxt, tgt[vi], w[vi] * pi)
        nxt /= nxt.sum()
        if np.abs(nxt - pi).max() < 1e-15:
            pi = nxt; break
        pi = nxt
    return (3 ** k) * float((pi ** 2).sum())


# ---------- exact rational (few k, for the denominator-rate signature) ----------
def Xk_exact(p, k):
    N = 3 ** k; M = order_of(p, N); invp = pow(p, -1, N)
    powinvp = [pow(invp, v, N) for v in range(1, M + 1)]
    Z = sum(F(1, p ** v) for v in range(1, M + 1))
    w = [F(1, p ** v) / Z for v in range(1, M + 1)]
    cop = [r for r in range(N) if r % 3 != 0]; idx = {r: i for i, r in enumerate(cop)}; n = len(cop)
    # transition matrix K (row=from)
    K = [[F(0)] * n for _ in range(n)]
    for r in cop:
        b = (3 * r + 1) % N
        for vi in range(M):
            K[idx[r]][idx[(b * powinvp[vi]) % N]] += w[vi]
    # stationary via exact Gaussian (pi K = pi): solve (K^T - I) pi = 0, sum pi = 1
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
    pi = bb
    return F(3 ** k) * sum(x * x for x in pi)


def main():
    print("# PROBE P-HYDRA FAMILY -- (3,p)-Hydra Plancherel constants, same class?\n")
    primes = [2, 5, 7]

    # ---- A: build + gate + float convergence ----
    print("## A  family build, gate (p=2 vs known), float convergence + reachable precision")
    print("   Q3 LEAD: is S_{k+1}/S_k -> 3(p-1)/(p+1)?  (=1 at p=2 boundary; >1 diverge for p>2)")
    known = {1: F(2, 3), 2: F(10, 21)}
    for p in [2, 5, 7, 11]:
        Xs = {0: 1.0}
        for k in range(1, 11):
            Xs[k] = Xk_float(p, k)
        S = {k: Xs[k] - Xs[k - 1] for k in range(1, 11)}
        gate = ""
        if p == 2:
            g = max(abs(S[1] - float(known[1])), abs(S[2] - float(known[2])))
            gate = f"  [gate p=2 vs 2/3,10/21: {g:.1e}]"
        # Aitken accel on S_k (limit estimate) + rate
        ks = list(range(4, 11)); vals = [S[k] for k in ks]
        rate = (vals[-1] - vals[-2]) / (vals[-2] - vals[-3]) if abs(vals[-2] - vals[-3]) > 0 else float('nan')
        # Aitken
        s = vals
        ait = s[-1] - (s[-1] - s[-2]) ** 2 / ((s[-1] - s[-2]) - (s[-2] - s[-3])) if len(s) >= 3 else s[-1]
        Sinf_est = 2 * ait  # S_inf = 2 lim T; but here S_k IS the Plancherel increment; lim S_k -> ?
        pred = 3 * (p - 1) / (p + 1)
        print(f"   p={p:>2}: S_k k=7..10 = {', '.join(f'{S[k]:.5f}' for k in range(7,11))}{gate}")
        print(f"        S_{{k+1}}/S_k ~ {rate:.5f}   vs 3(p-1)/(p+1) = {pred:.5f}   |diff|={abs(rate-pred):.1e}")
    print()

    # ---- B: exact denominator-rate signature per sibling (the MAHLER test) ----
    print("## B  exact S_k^(3,p): denominator bits log_p(den) and ratio -> doubly-exp (infinite-depth)?")
    for p in primes:
        Xe = {0: F(1)}
        kmax = 5 if p == 2 else 4          # exact gets heavy; p=2 cheaper (smaller weights)
        t0 = time.time()
        Se = {}
        for k in range(1, kmax + 1):
            Xe[k] = Xk_exact(p, k)
            Se[k] = Xe[k] - Xe[k - 1]
        print(f"   p={p}  (exact k=1..{kmax}, {time.time()-t0:.0f}s):")
        print(f"     {'k':>2} {'S_k^(3,p)':>16} {'log_p(den) bits':>16} {'ratio':>7}")
        Dprev = None
        for k in range(1, kmax + 1):
            D = math.log(Se[k].denominator, p) if Se[k] != 0 else 0
            ratio = (D / Dprev) if Dprev and Dprev > 0 else float('nan')
            v = str(Se[k]) if len(str(Se[k])) < 20 else f"{float(Se[k]):.10f}"
            print(f"     {k:>2} {v:>16} {D:>16.2f} {ratio:>7.3f}")
            Dprev = D
        print(f"     [ratio -> 3 (>2) => infinite Mahler order (same class as Collatz); bounded => finite-depth]")
    print()

    # ---- D: unified 2-parameter law ratio(q,p) = q(p-1)/(p+1), mixed cases ----
    print("## D  UNIFIED LAW  S_{k+1}/S_k -> q(p-1)/(p+1)  (q/3 at p=2; 3(p-1)/(p+1) at q=3)")

    def Xk_qp(q, p, k, vmax=300):
        N = q ** k; M = order_of(p, N); vm = min(M, vmax); invp = pow(p, -1, N)
        w = np.array([(1.0 / p) ** v for v in range(1, vm + 1)]); w /= w.sum()
        powinvp = np.array([pow(invp, v, N) for v in range(1, vm + 1)], dtype=np.int64)
        cop = np.array([r for r in range(N) if r % q != 0], dtype=np.int64)
        n = len(cop); idx = -np.ones(N, np.int64); idx[cop] = np.arange(n)
        base = ((q * cop + 1) % N).astype(np.int64)
        tgt = np.stack([idx[(base * powinvp[vi]) % N] for vi in range(vm)])
        pi = np.full(n, 1.0 / n)
        for _ in range(4000):
            nxt = np.zeros(n)
            for vi in range(vm):
                np.add.at(nxt, tgt[vi], w[vi] * pi)
            nxt /= nxt.sum()
            if np.abs(nxt - pi).max() < 1e-15: pi = nxt; break
            pi = nxt
        return (q ** k) * float((pi ** 2).sum())

    print(f"   {'(q,p)':>8} {'meas S_{k+1}/S_k':>16} {'q(p-1)/(p+1)':>14} {'|diff|':>9}")
    for q, p in [(3, 2), (5, 2), (7, 2), (3, 5), (3, 7), (5, 7), (7, 5)]:
        kmax = 8 if q <= 5 else 7
        X = {0: 1.0}
        for k in range(1, kmax + 1): X[k] = Xk_qp(q, p, k)
        S = {k: X[k] - X[k - 1] for k in range(1, kmax + 1)}
        meas = S[kmax] / S[kmax - 1]
        pred = q * (p - 1) / (p + 1)
        print(f"   {f'({q},{p})':>8} {meas:>16.5f} {pred:>14.5f} {abs(meas-pred):>9.1e}")
    print("   => boundary (ratio=1) is the CURVE q(p-1)=p+1; Collatz (3,2) sits on it (3*1=3=p+1).\n")

    print("## C  VERDICT")
    print("   (read ratios above: all siblings ratio->3 => SAME CLASS, all infinite-depth => PSLQ-on-values")
    print("    infeasible AND unnecessary at structure level; independence is a Hardouin-criterion PEN result,")
    print("    not a digit hunt. If any sibling's ratio is bounded => it's a finite-depth 'nice' member.)")


if __name__ == "__main__":
    main()
