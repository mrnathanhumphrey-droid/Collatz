"""
PROBE P6E (Wilson) -- gate the collapsed cross-parity autocorrelation formula + read off the true constants (2026-07-26).

Wilson's explicit form, from nu_o = 1/2 nu_e(.+1) + 1/2 beta (beta=(m_1)_*nu), R_e := nu_e (corr) nu_e~:
    X(m) = 1/2[R_e(m+1)+R_e(m-1)] + 1/2[(nu_e corr beta)(m)+(beta corr nu_e)(m)]      (base-2 lag m)
         = R_e(m) + 1/2 D^2 R_e(m) + boundary
where X = cross-parity autocorrelation = nu_e (corr) nu_o~ + nu_o (corr) nu_e~. Base-4 lag k = base-2 lag 2k.

DELIVERABLES:
 (A) GATE X(m) == 1/2[R_e(m+1)+R_e(m-1)] + boundary  for ALL m, j=2..6. Machine precision (algebra from P6D collapse).
 (B) CONSTANTS (symbolic, Fractions): reindex Lambda = Sum_{k>=1} 4^-k A(k), A(k)=X(2k), collect coeff of R_e(n).
     Predict c_n = (5/4)2^-n (odd n>=3), 1/8 (n=1), 0 (even). = Wilson's constants. Compare to 2^-n (twice odd D~).
     => kernel = odd-part-of-D~ smeared by the (I + 1/2 D^2) flanking kernel; 5/4 = 1+1/4 (two flanks), 1/8 = n=1 defect.
 (C) D~ IDENTITY (Fractions): Sum_{n odd} 2^-n z^n == D~(z)-D~(-z), D~(u)=1/(2-u)=Sum_{v>=0} 2^-(v+1) u^v (a>=0 conv).
 (D) VALIDATION: in base-2, X(2k) (cross) vs S(2k):=R_ee+R_oo minus cross (same) reproduces P6B rule 3-nmid k=>cross.
     And ratio X(2) / A_j(1)_certified (bridge) to tie base-2 residue autocorr to the certified channel.

Reuses probe_p6d.build_base2 (certified base-2 one-step) + probe_p6b.shellA (certified channel). No new transport.
"""
import os, sys, time
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_p6d import build_base2
from probe_p1 import build_level
from probe_p6b import shellA


def corr(f, g):
    """circular cross-correlation corr[m] = sum_t f[t] g[(t-m) mod N]."""
    return np.fft.ifft(np.fft.fft(f) * np.conj(np.fft.fft(g))).real


def main():
    t0 = time.time()
    print("# PROBE P6E -- collapsed cross-parity autocorrelation: gate + constants\n")

    # ---------- (A) GATE the X(m) formula ----------
    print("## (A) GATE  X(m) == 1/2[R_e(m+1)+R_e(m-1)] + boundary   (all m, j=2..6)")
    for n in (2, 3, 4, 5, 6):
        S = build_base2(n); twoN = S['twoN']
        nu_e = S['R_e']; nu_o = S['R_o']; beta = 2.0 * S['B']       # beta = (m_1)_*nu = 2*B (B=1/2 push_1)
        Ree = corr(nu_e, nu_e)
        X = corr(nu_e, nu_o) + corr(nu_o, nu_e)
        bnd = 0.5 * (corr(nu_e, beta) + corr(beta, nu_e))
        pred = 0.5 * (np.roll(Ree, -1) + np.roll(Ree, 1)) + bnd      # roll(-1)[m]=Ree[m+1]; roll(1)[m]=Ree[m-1]
        gate = np.max(np.abs(X - pred))
        # also show the second-difference phrasing residual
        D2 = np.roll(Ree, -1) - 2 * Ree + np.roll(Ree, 1)
        pred2 = Ree + 0.5 * D2 + bnd
        gate2 = np.max(np.abs(X - pred2))
        print(f"   n={n} (twoN={twoN}):  max|X - [1/2(R_e(m+1)+R_e(m-1))+bnd]| = {gate:.2e}   "
              f"|X - [R_e+1/2 D^2 R_e+bnd]| = {gate2:.2e}   [{'EXACT' if max(gate,gate2)<1e-12 else 'approx'}]")
    print()

    # ---------- (B) CONSTANTS: symbolic reindex of Lambda = Sum 4^-k A(k), A(k)=X(2k) ----------
    print("## (B) CONSTANTS (symbolic): coeff of R_e(n) in Lambda = Sum_{k>=1} 4^-k * 1/2[R_e(2k+1)+R_e(2k-1)]")
    KMAX = 40
    c = {}
    for k in range(1, KMAX + 1):
        w = Fr(1, 4 ** k) * Fr(1, 2)
        for n in (2 * k + 1, 2 * k - 1):
            c[n] = c.get(n, Fr(0)) + w
    print("   n : c_n           (5/4)2^-n ?   2^-n (twice odd D~)   ratio c_n/2^-n")
    for n in (1, 3, 5, 7, 9):
        cn = c[n]; wilson = Fr(5, 4) * Fr(1, 2 ** n) if n >= 3 else Fr(1, 8)
        twice = Fr(1, 2 ** n)
        print(f"   {n} : {str(cn):>14}  {'MATCH' if cn == wilson else 'NO':>5}  "
              f"{str(twice):>10}   {str(cn / twice):>6}")
    print("   => odd n>=3: c_n=(5/4)2^-n (=1+1/4, the two flanks); n=1: 1/8 (upper flank only, k=0 absent).")
    print("      kernel = odd-part-of-D~ under (I + 1/2 D^2); NOT bare 2^-n. Wilson's constants confirmed exact.")
    print("   [CARE, unresolved: A(k)=X(2k) assumes cross for ALL k. Only 3-nmid k are cross; 3|k are SAME-parity")
    print("    (separate kernel). True Lambda = Sum_{3-nmid k}4^-k X(2k) + Sum_{3|k}4^-k A_same(k). See (D) for the split.]\n")

    # ---------- (C) D~ identity ----------
    print("## (C) D~ identity: Sum_{n odd}2^-n z^n == D~(z)-D~(-z), D~(u)=Sum_{v>=0}2^-(v+1)u^v = 1/(2-u)")
    NC = 12
    lhs = {n: Fr(1, 2 ** n) if n % 2 == 1 else Fr(0) for n in range(1, NC)}
    Dz = {v: Fr(1, 2 ** (v + 1)) for v in range(NC)}          # D~(z) coeffs
    rhs = {n: Dz.get(n, Fr(0)) - (Fr(-1) ** n) * Dz.get(n, Fr(0)) for n in range(1, NC)}
    ok = all(lhs[n] == rhs[n] for n in range(1, NC))
    print(f"   coefficient match to z^{NC-1}: {'IDENTITY HOLDS' if ok else 'FAILS'}  "
          f"(e.g. n=1:{lhs[1]}=={rhs[1]}, n=3:{lhs[3]}=={rhs[3]}, n=5:{lhs[5]}=={rhs[5]})\n")

    # ---------- (D) VALIDATION: base-2 parity rule + channel tie ----------
    print("## (D) base-2 X(2k): cross vs same reproduces P6B (3-nmid k => cross); and X(2)/A_j(1)_certified")
    for n in (3, 4, 5):
        S = build_base2(n); twoN = S['twoN']
        nu_e = S['R_e']; nu_o = S['R_o']
        Xc = corr(nu_e, nu_o) + corr(nu_o, nu_e)                  # cross-parity autocorr
        Ss = corr(nu_e, nu_e) + corr(nu_o, nu_o)                  # same-parity autocorr
        Ntot = twoN
        row = []
        for k in (1, 2, 3, 4, 5, 6):
            m = (2 * k) % Ntot
            verdict = 'CROSS' if abs(Xc[m]) > 10 * abs(Ss[m]) else ('SAME' if abs(Ss[m]) > 10 * abs(Xc[m]) else 'mix')
            row.append(f"k={k}:{verdict}")
        print(f"   n={n}: " + "  ".join(row) + f"   [3-nmid k => CROSS, 3|k(k=3,6) => SAME]")
        # channel tie: certified A_j(1) via bridge vs base-2 X(2)
        L = build_level(n)
        A1 = shellA(L, L['What'], 1)
        X2 = float(Xc[(2) % Ntot])
        print(f"        certified A_{n}(1)={A1:+.6f}   base-2 X(2)={X2:+.6f}   ratio X(2)/A1={X2/A1 if abs(A1)>1e-12 else float('nan'):+.4f}")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
