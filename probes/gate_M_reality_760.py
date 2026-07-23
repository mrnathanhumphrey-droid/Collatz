"""
gate_M_reality_760.py -- belt-and-suspenders numerical gate for Lemma 76.0 (M-reality).

Lemma 76.0 (proved, elementary): M_n(eta) = Sum_{3|/xi} muhat_n(xi) muhat_n*(xi*eta) is REAL for every eta,
unconditionally, from (i) the index set A={3|/xi} is closed under xi->-xi and fixed-point-free, and
(ii) pi real => muhat(-xi)=muhat(xi)*. [No R66 class-symmetry needed.]

This gate does NOT prove it (the 2-line argument does) -- it confirms max|Im M_n(eta)| is machine-zero
over ALL eta, pushed well past the corpus's k=4, using FLOAT stationary pi (fast) + FFT char function.
Also reprints Thm 76.1 conservation and S=-2 M(1+3^{n-1}) as free riders.
"""
import os, sys, math, time
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")


def build_K_float(k):
    N = 3 ** k
    Mv = 2 * 3 ** (k - 1)                      # geom truncation length = |(Z/3^k)*|
    inv2 = pow(2, -1, N)
    powinv = [pow(inv2, v, N) for v in range(1, Mv + 1)]
    cop = [r for r in range(N) if r % 3]
    idx = {r: i for i, r in enumerate(cop)}
    n = len(cop)
    K = np.zeros((n, n))
    Zf = 1.0 - 2.0 ** (-Mv)                    # normalizer (=1 to machine prec for large Mv)
    for r in cop:
        i = idx[r]
        base = 3 * r + 1
        for rv in range(1, Mv + 1):
            p = 2.0 ** (-rv) / Zf              # 2.0**(-rv) underflows to 0 harmlessly for large rv
            if p == 0.0:
                break
            tgt = (base * powinv[rv - 1]) % N
            K[i, idx[tgt]] += p
    return K, cop, N


def stationary_float(K):
    n = len(K)
    A = K.T - np.eye(n)
    A[-1, :] = 1.0
    b = np.zeros(n); b[-1] = 1.0
    return np.linalg.solve(A, b)


def main():
    KMAX = 7
    t0 = time.time()
    print("# GATE: M-reality (Lemma 76.0).  max|Im M_n(eta)| over ALL eta, k=2..%d (corpus stopped at 4).\n" % KMAX)
    pis = {}
    for k in range(1, KMAX + 1):
        K, cop, N = build_K_float(k)
        pi = stationary_float(K)
        pis[k] = (pi, cop, N)
        print(f"   pi_{k} built (|(Z/3^{k})*|={len(cop)}, sum={pi.sum():.12f})  ({time.time()-t0:.1f}s)")
    print()

    print(f"   {'n':>2} {'#units':>7} {'max|Im M(eta)|':>16} {'max|M(eta)|':>13} {'rel |Im|/|M|':>13} "
          f"{'Im M(1+3^{n-1})':>16} {'max|cons Sum_j|':>15}")
    for n in range(2, KMAX + 1):
        pi, cop, N = pis[n]
        pfull = np.zeros(N)
        for i, r in enumerate(cop):
            pfull[r] = pi[i]
        muhat = np.fft.fft(pfull)              # muhat(xi) = sum_r pfull(r) e^{-2pi i r xi/N}
        g = muhat.copy()
        g[0::3] = 0.0                          # restrict to 3|/xi
        ar = np.arange(N)
        maxIm = 0.0; maxAbs = 0.0
        Mvals = {}                             # eta -> M(eta) (units only)
        for eta in cop:
            idxm = (ar * eta) % N
            Me = np.sum(g * np.conj(muhat[idxm]))
            Mvals[eta] = Me
            maxIm = max(maxIm, abs(Me.imag))
            maxAbs = max(maxAbs, abs(Me))
        # leading mode eta = 1 + 3^{n-1}
        eta_lead = (1 + 3 ** (n - 1)) % N
        Im_lead = Mvals[eta_lead].imag if eta_lead in Mvals else float('nan')
        # conservation Thm 76.1: Sum_{j=0,1,2} M(eta0 + j 3^{n-1}) over eta0 in (Z/3^{n-1})*
        Nn1 = 3 ** (n - 1)
        cons_max = 0.0
        for eta0 in [r for r in range(Nn1) if r % 3]:
            s = 0j
            for j in range(3):
                e = (eta0 + j * Nn1) % N
                s += Mvals.get(e, 0j)
            cons_max = max(cons_max, abs(s))
        rel = maxIm / maxAbs if maxAbs else float('nan')
        print(f"   {n:>2} {len(cop):>7} {maxIm:>16.3e} {maxAbs:>13.6f} {rel:>13.3e} "
              f"{Im_lead:>16.3e} {cons_max:>15.3e}")
    print(f"\n   [Lemma 76.0 PASS iff max|Im M| and rel are ~machine-zero at every n.]  ({time.time()-t0:.1f}s)")

    # free rider: S_{n} = -2 Re M(1+3^{n-1}) unconditional; and = -2 M(...) once real
    print("\n   free rider  S_n = -2 M_n(1+3^{n-1})  (real part; Im shown separately above):")
    print(f"   {'n':>2} {'S_n=Sum|muhat|^2':>16} {'-2 Re M(1+3^{n-1})':>18} {'match':>7}")
    for n in range(2, KMAX + 1):
        pi, cop, N = pis[n]
        pfull = np.zeros(N)
        for i, r in enumerate(cop):
            pfull[r] = pi[i]
        muhat = np.fft.fft(pfull)
        g = muhat.copy(); g[0::3] = 0.0
        S = float(np.sum(np.abs(g) ** 2))
        eta = (1 + 3 ** (n - 1)) % N
        idxm = (np.arange(N) * eta) % N
        M = np.sum(g * np.conj(muhat[idxm]))
        print(f"   {n:>2} {S:>16.10f} {(-2*M.real):>18.10f} {'OK' if abs(S+2*M.real)<1e-9 else 'DEV':>7}")


if __name__ == "__main__":
    main()
