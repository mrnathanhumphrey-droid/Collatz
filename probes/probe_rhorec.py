"""
PROBE RHOREC -- definitive: does Wilson's rho-hat recursion REPRODUCE fft(rho_n)? (2026-07-26)

In the codebase dlog convention: chi_a(Y) = e(2pi i a * d_n[(Y-1)//3] / 3^n). For Y=3X+1: (Y-1)//3 = X => chi_a(3X+1)
= e(2pi i a d_n[X]/3^n). And dlog(X_n) = dlog((3X+1)(-2)^-v) = d_n[X] - v*h (h=dlog(-2)=(3^n+1)/2). So
    rho-hat_n(a) = [Sum_{v>=1} 2^-v e(-2pi i a v h/3^n)] * [Sum_X nu_{n-1}(X) e(2pi i a d_n[X]/3^n)]
                 =  vfactor(a) * E_affine(a),   vfactor(a) = z/(2-z),  z=(-1)^a e^{-i pi a/3^n}.
(Wilson's D~=1/(2-z) is vfactor up to the unit phase z; the bound |vfactor|^2=1/(5-4(-1)^a cos) is unaffected.)
GATE: build g[s]=Sum_{X: d_n[X]=s} nu_{n-1}(X); E_affine=fft(g); RHS=vfactor*E_affine; compare to fft(rho_n).
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu


def dense_rho(nu_k, N, r):
    d = R10.dlog_table(r)
    mu = np.zeros(N)
    for X, w in nu_k.items():
        mu[(X - 1) // 3 % N] += float(w)
    rr = np.zeros(N); rr[np.fromiter((d[a] for a in range(N)), np.int64, N)] = mu
    return rr / rr.sum()


def main():
    t0 = time.time()
    print("# PROBE RHOREC -- does rho-hat_n = vfactor * E_affine reproduce fft(rho_n)?\n")
    nus = build_nu(0.5, 6)
    for n in range(2, 7):
        N = 3 ** n; Nm = 3 ** (n - 1)
        dn = R10.dlog_table(n)
        rho_n = dense_rho(nus[n], N, n)
        lhs = np.fft.fft(rho_n)                       # fft(rho_n)(a)
        # E_affine: g[s] = sum_{X in nu_{n-1}} w * [d_n[(X-1)//3 ... ] ] -- chi_a(3X+1)=e(a d_n[X]/N), X = pre-index at n-1
        tot = sum(float(w) for w in nus[n - 1].values())
        g = np.zeros(N)
        for X, w in nus[n - 1].items():
            Xp = (X - 1) // 3 % Nm                    # pre-dlog index of X at level n-1
            g[dn[Xp]] += float(w) / tot               # chi_a(3X+1) uses d_n at pre-index Xp
        E_aff = np.fft.fft(g)
        a = np.arange(N)
        z = ((-1.0) ** a) * np.exp(-1j * np.pi * a / N)
        vfac = z / (2 - z)
        rhs = vfac * E_aff
        # compare on primitive a (a=0 is the trivial mode)
        prim = (a % 3 != 0)
        rel = np.abs(rhs[prim] - lhs[prim]).mean() / (np.abs(lhs[prim]).mean() + 1e-30)
        # also try Xp = X directly (residue as index) in case convention differs
        g2 = np.zeros(N)
        for X, w in nus[n - 1].items():
            g2[dn[X % Nm]] += float(w) / tot
        rhs2 = vfac * np.fft.fft(g2)
        rel2 = np.abs(rhs2[prim] - lhs[prim]).mean() / (np.abs(lhs[prim]).mean() + 1e-30)
        print(f"   n={n}: rel(pre-index) = {rel:.3e} {'REPRODUCES' if rel<1e-9 else ('close' if rel<1e-3 else 'NO')}"
              f"  |  rel(residue-index) = {rel2:.3e} {'REPRODUCES' if rel2<1e-9 else ('close' if rel2<1e-3 else 'NO')}")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
