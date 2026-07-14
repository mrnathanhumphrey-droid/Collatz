"""
Probe 4 (PRE_REG_4_CTILDE_ORD2) — does the c~_q deviation track ord_q(2)?

c~_q(k) = S_k/(q/3)^k,  S_k = X_k - X_{k-1},  X_k = q^k * sum_i pi_i^2,
pi = stationary of the qx+1 Syracuse chain on coprime residues mod q^k.
delta_q = c~_q(2) - (q-3)/q.   Test whether delta_q decreases with ord_q(2).

q=7 (ord 3) is the unique small outlier (delta ~ +0.21). New primes at controlled
orders break the degeneracy: 31(5), 127(7), 73(9); baselines 5(4),11(10),13(12),17(8).

Exact (Fraction) stationary for q<=17 (validation); float sparse power-iteration
for large q (chain mixes fast, |lambda_2|~1e-4 per result_qspectrum). Cross-checked.
"""
import sys, numpy as np
import scipy.sparse as sp
from fractions import Fraction
sys.stdout.reconfigure(encoding="utf-8")

def order_of_two(N):
    m, v = 1, 2 % N
    while v != 1:
        v = (v*2) % N; m += 1
    return m

def X_exact(q, k):
    """Exact X_k via rational stationary (small q only)."""
    N = q**k; M = order_of_two(N); inv2 = pow(2, -1, N)
    powinv = [pow(inv2, v, N) for v in range(1, M+1)]
    coprime = [r for r in range(N) if r % q != 0]
    idx = {r: i for i, r in enumerate(coprime)}; n = len(coprime)
    Z = Fraction(2**M - 1, 2**M)
    K = [[Fraction(0)]*n for _ in range(n)]
    for r in coprime:
        for v in range(1, M+1):
            p = Fraction(1, 2**v)/Z
            t = ((q*r+1)*powinv[v-1]) % N
            K[idx[r]][idx[t]] += p
    # stationary via Gaussian elimination over Q
    A = [[K[j][i] - (Fraction(1) if i == j else 0) for j in range(n)] for i in range(n)]
    A[n-1] = [Fraction(1)]*n; b = [Fraction(0)]*n; b[n-1] = Fraction(1)
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        if piv != c: A[c], A[piv] = A[piv], A[c]; b[c], b[piv] = b[piv], b[c]
        d = A[c][c]
        A[c] = [x/d for x in A[c]]; b[c] /= d
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [A[r][j]-f*A[c][j] for j in range(n)]; b[r] -= f*b[c]
    return Fraction(q**k) * sum(p*p for p in b)

def X_float(q, k):
    """Float X_k via sparse power iteration."""
    N = q**k; M = order_of_two(N); inv2 = pow(2, -1, N)
    powinv = [pow(inv2, v, N) for v in range(1, M+1)]
    coprime = [r for r in range(N) if r % q != 0]
    idx = {r: i for i, r in enumerate(coprime)}; n = len(coprime)
    Z = (2**M - 1)/2**M
    rows, cols, vals = [], [], []
    for r in coprime:
        i = idx[r]
        for v in range(1, M+1):
            t = ((q*r+1)*powinv[v-1]) % N
            rows.append(i); cols.append(idx[t]); vals.append((0.5**v)/Z)
    K = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    Kt = K.T.tocsr()
    pi = np.full(n, 1.0/n)
    for _ in range(200):
        nxt = Kt.dot(pi); nxt /= nxt.sum()
        if np.abs(nxt - pi).sum() < 1e-15:
            pi = nxt; break
        pi = nxt
    return q**k * float(np.dot(pi, pi))

if __name__ == "__main__":
    QS = [5, 7, 11, 13, 17, 31, 73, 127]
    print("Probe 4: c~_q deviation vs ord_q(2)\n" + "="*72)
    # validation: float vs exact for small q
    print("Validation (float X_2 vs exact X_2):")
    for q in [5, 7, 11, 13]:
        xe, xf = float(X_exact(q, 2)), X_float(q, 2)
        print(f"  q={q}: exact={xe:.9f} float={xf:.9f}  {'OK' if abs(xe-xf)<1e-7 else 'MISMATCH'}")
    print()

    rows = []
    print(f"{'q':>4} {'ord2':>5} {'c~_q(2)':>10} {'(q-3)/q':>9} {'delta':>9}")
    for q in QS:
        X0, X1, X2 = 1.0, X_float(q, 1), X_float(q, 2)
        S2 = X2 - X1
        ctilde = S2 / (q/3)**2
        base = (q-3)/q
        d = ctilde - base
        ordq = order_of_two(q)
        rows.append((q, ordq, ctilde, base, d))
        print(f"{q:>4} {ordq:>5} {ctilde:>10.5f} {base:>9.5f} {d:>+9.5f}")

    print("\nSorted by ord_q(2):")
    print(f"{'ord2':>5} {'q':>4} {'delta':>10}")
    for q, o, c, b, d in sorted(rows, key=lambda x: x[1]):
        print(f"{o:>5} {q:>4} {d:>+10.5f}")

    # regressions
    ordv = np.array([r[1] for r in rows], float)
    qv   = np.array([r[0] for r in rows], float)
    dv   = np.array([r[4] for r in rows], float)
    print("\nRegressions of |delta| on candidate covariates (R^2):")
    for name, x in [("ord", ordv), ("1/ord", 1/ordv), ("ord/q", ordv/qv), ("1/(ord*q)", 1/(ordv*qv))]:
        A = np.vstack([x, np.ones_like(x)]).T
        coef, *_ = np.linalg.lstsq(A, np.abs(dv), rcond=None)
        pred = A @ coef
        ss_res = np.sum((np.abs(dv)-pred)**2); ss_tot = np.sum((np.abs(dv)-np.abs(dv).mean())**2)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else float('nan')
        print(f"  |delta| ~ {name:>10}: R^2={r2:.4f}  slope={coef[0]:+.5f}")
    print("\n(H_ORD fires if |delta| is monotone decreasing in ord_q(2) with the new primes ordered.)")
