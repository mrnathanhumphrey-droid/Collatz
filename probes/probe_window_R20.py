"""
PROBE R20 -- THE THIN WINDOW. Reuses R7/R9/R10. Per-r shape statistics (no fit).

Key identity (Wilson, verified): w(k)=1/(4chi_k(4)-1)=Sum_{m>=1}4^{-m}e(-mk/3^r) => Lambda_r = Sum_{m>=1} 4^{-m} A_r(m),
  A_r(m)=gamma_r(tau_m)-gamma_{r-1}(tau_m)=C_{r+1}(m)/3 (REAL exact, all m). The R7 channel engine IS the
  deviation-field route. Crude |A_r(m)|<=S_r + weight => tail beyond m~r is <=(4/3)4^{-r}S_r (summable). So the
  needed control is on A_r(m) for m<~r only -- a logarithmically thin window.
R20-A window m=1..r.  R20-B m=9,27 settle R13-C.  R20-C additive argmax trajectory.  R20-D A_r(3^{r-1})=?=-S_r/2
(vacate). R20-E Lambda_r = Sum 4^{-m}A_r(m) == OffDiag/2 (engine==deviation-field).
"""
import os, sys, math, cmath, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import probe_engine_R7 as R7
import probe_gamma_R9 as R9
import probe_charledger_R10 as R10

v3 = R9.v3
_hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                    'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in _hist.items()}
S = {k: F(7, 15) + EPS[k] for k in EPS}


def Lunif(r):
    tr = F(3 ** r, 4 ** (3 ** r) - 1) - F(3 ** (r - 1), 4 ** (3 ** (r - 1)) - 1)
    return S[r] * tr / (2 * 3 ** (r - 1))


def A_r(mu, r, m):
    return R9.gamma(mu[r], r, R9.tau(m, r)) - R9.gamma(mu[r - 1], r - 1, R9.tau(m, r - 1))


def mu_hat(mu_r, r, xi):
    N = 3 ** r
    return sum(complex(p) * cmath.exp(2j * math.pi * (xi * a % N) / N) for a, p in mu_r.items())


def main():
    print("# PROBE R20 -- THE THIN WINDOW.\n")
    mu = {1: R7.mu1()}
    for k in range(2, 8):
        mu[k] = R7.build_mu(mu[k - 1], k)

    Lam = {r: (EPS[r + 1] - EPS[r]) / 2 for r in range(1, 8)}
    Lun = {r: Lunif(r) for r in range(1, 8)}
    b = {r: (Lam[r] - Lun[r]) / S[r] for r in range(2, 8)}

    # full A-spectra (exact) for norms
    print("(building full A-spectra for the stratum norms, r=3..7 ...)")
    Aspec = {r: {m: A_r(mu, r, m) for m in range(1, 3 ** r)} for r in range(3, 8)}
    normA = {r: sum(a * a for a in Aspec[r].values()) for r in Aspec}

    # ================= R20-A =================
    print("\n## R20-A  THE m<=r WINDOW (measurement, NO fit): all m=1..r, ratio to stratum-typical, weighted running sum")
    for r in range(3, 8):
        print(f"   r={r}  (Lambda_r={float(Lam[r]):+.6e}, b_r={float(b[r]):+.6e}, ||delta||^2_A={float(normA[r]):.4f})")
        print(f"     {'m':>3} {'v3':>3} {'A_r(m)':>13} {'|A|^2/typ':>10} {'4^-m A_r(m)':>14} {'running sum':>14}")
        run = F(0)
        for m in range(1, r + 1):
            A = Aspec[r][m]; j = v3(m); Nj = 2 * 3 ** (r - 1 - j)
            typ = (float(normA[r]) / r) / Nj
            term = F(1, 4 ** m) * A; run += term
            print(f"     {m:>3} {j:>3} {float(A):>+13.6f} {float(A)**2/typ:>10.3f} {float(term):>+14.6e} {float(run):>+14.6e}")
        print(f"     [running sum at m=r vs Lambda_r={float(Lam[r]):+.6e}: gap={float(run-Lam[r]):+.2e}]")
    print("   [Q: ratio O(1) across the whole window? running sum saturated by m~5?]\n")

    # ================= R20-B =================
    print("## R20-B  m=9 and m=27, SETTLED (measurement): is j>=2 exceptional or was R13-C reading an oscillation?")
    for mm in (9, 27):
        print(f"   --- m={mm} (v3={v3(mm)}) ---")
        print(f"     {'r':>2} {'gamma_r(tau_m)':>16} {'A_r(m)=succ.diff':>18} {'|A|^2':>12} {'|A|^2/typ':>10} {'sign':>5}")
        gprev = None
        for r in range(2, 8):
            g = R9.gamma(mu[r], r, R9.tau(mm, r))
            A = (g - gprev) if gprev is not None else None
            if A is None:
                print(f"     {r:>2} {float(g):>16.6f} {'--':>18} {'--':>12} {'--':>10} {'--':>5}")
            else:
                if mm >= 3 ** r:
                    typ = float('nan'); ratio = float('nan')
                else:
                    j = v3(mm); Nj = 2 * 3 ** (r - 1 - j)
                    typ = (float(normA[r]) / r) / Nj if r in normA else float('nan')
                    ratio = float(A) ** 2 / typ
                sg = '+' if A > 0 else '-'
                note = " (=-S_r/2 DEF)" if mm == 3 ** (r - 1) else ""
                print(f"     {r:>2} {float(g):>16.6f} {float(A):>+18.6f} {float(A)**2:>12.3e} {ratio:>10.3f} {sg:>5}{note}")
            gprev = g
    print("   [R13-C reported m=9 diffs +0.027,+0.073,+0.046 (r=5,6,7). Monotone growth => exceptional; turnover/sign => oscillation.]\n")

    # ================= R20-C =================
    print("## R20-C  ADDITIVE ARGMAX TRAJECTORY (measurement, NO fit; 'power of 2' is VACUOUS -- 2 is a primitive root)")
    print(f"   {'r':>2} {'argmax xi':>10} {'x=xi/N':>9} {'dist-to-0':>10} {'xi_r/xi_{r-1}':>13} {'(2/3)^r':>9}")
    prevxi = None
    for r in range(2, 8):
        N = 3 ** r
        xs = {xi: abs(mu_hat(mu[r], r, xi)) for xi in range(1, N) if xi % 3 != 0}
        xstar = max(xs, key=xs.get)
        d0 = min(xstar / N, 1 - xstar / N)
        rat = f"{xstar/prevxi:.4f}" if prevxi else "--"
        print(f"   {r:>2} {xstar:>10} {xstar/N:>9.5f} {d0:>10.5f} {rat:>13} {(2/3)**r:>9.5f}")
        prevxi = xstar
    print("   [Q: x->0 (near-trivial additive slow mode)? at what rate? is r=3 (19) a non-monotone excursion?]\n")

    # ================= R20-D =================
    print("## R20-D  THE 0.233 CHECK (ALGEBRAIC, forced): is A_r(3^{r-1}) == -S_r/2 exactly?  (R9-D precedent: vacate)")
    allvac = True
    for r in range(2, 8):
        m = 3 ** (r - 1)
        A = A_r(mu, r, m)
        eq = (A == -S[r] / 2)
        allvac = allvac and eq
        print(f"   r={r}: A_r(3^{{r-1}})={float(A):+.6f}  == -S_r/2={float(-S[r]/2):+.6f}? {eq}")
    print(f"   => {'VACATE: the coarsest non-DC mode is -S_r/2 (two conjugate members, DC-split) -- true and vacuous, ' if allvac else 'NOT definitional: '}"
          f"{'NOT a second constant (=-S_inf/2=-7/30). R19-D A-side argmax was this artifact.' if allvac else 'print both.'}\n")

    # ================= R20-E =================
    print("## R20-E  WEIGHT IDENTITY (ALGEBRAIC, forced): Lambda_r = Sum_{m=1..3^r} 4^-m A_r(m) == OffDiag_{r+1}/2")
    okE = True
    for r in range(2, 6):
        P = 3 ** r; geom = 1 - F(1, 4 ** P)
        tot = sum(F(1, 4 ** m) / geom * A_r(mu, r, m) for m in range(1, P + 1))
        L10 = R10.Lambda_r(mu, r)[0]
        off = (S[r + 1] - S[r]) / 2
        good = (tot == L10 == off)
        okE = okE and good
        print(f"   r={r}: Sum 4^-m A_r(m) == Lambda_r(R10) == OffDiag/2 ? {tot==L10} & {L10==off}  [{'OK' if good else 'DEV'}]")
    print(f"   => R20-E {'PASS -- deviation-field route == R7 channel engine (one computation, not two)' if okE else 'DEV'}")


if __name__ == "__main__":
    main()
