"""
PROBE R28 -- THE GAP-CORRELATION OPERATOR (R27-E). Certified renewal X'=1+3*2^{-v}X on 1+3Z/3^{r+1}.

Wilson's operator:
  nu_hat_r(xi) = e(xi/3^{r+1}) * E_v[nu_hat_{r-1}(xi 2^{-v})],  nu = law of X_r.
  R_r(d) := sum_{eta in Z/3^{r+1}} nu_hat_r(eta) * conj(nu_hat_r(eta*2^d)).
  Claimed: S_r = 3 * sum_d P(d) R_{r-1}(d),  P(0)=(1-lam)/(1+lam), P(+-d)=(1-lam)lam^d/(1+lam).
  d=0 term = 3P(0) R_{r-1}(0) = S_{r-1} at lam=1/2 (Diagonal Flatness). Odd-d vanish (conjugate-kill, ord_3(2)=2).

A: GATE -- verify S_r = 3 sum_d P(d) R_{r-1}(d) exact-ish, r=2..7; odd-d vanish; d=0 = S_{r-1}. (normalization probed.)
B: kappa measured -- m=+-1 twisted sum / untwisted, per r. Is kappa<1 every r, stabilizing?
C: channel mixing -- |m=+-2|,|+-3| relative to |m=+-1|, per r.
D: real+pair fit -- deflate Lambda by fixed real mode 0.5, fit the complex pair; reproduce signs -,-,+,+,+,-,+ & period.
"""
import os, sys, math, json, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np


def build_nu(lam, RMAX, tol=1e-18):
    """nu_r = law of X_r under X'=1+3*2^{-v}X (certified), as dense float arrays. nu_r on Z/3^{r+1}."""
    nus = {}
    nu = {1: 1.0}                                  # nu_0 = delta_1 on Z/3
    nus[0] = nu
    for r in range(1, RMAX + 1):
        M = 3 ** (r + 1); Mr = 3 ** r; inv2 = pow(2, -1, Mr)
        new = {}
        for Xp, wp in nu.items():
            u = inv2; v = 1
            while (1 - lam) * lam ** (v - 1) > tol:
                pv = (1 - lam) * lam ** (v - 1)
                t = (u * Xp) % Mr
                Xr = (1 + 3 * t) % M
                new[Xr] = new.get(Xr, 0.0) + pv * wp
                u = (u * inv2) % Mr; v += 1
        nu = new; nus[r] = nu
    return nus


def nu_hat(nu, M):
    dense = np.zeros(M)
    for X, w in nu.items():
        dense[X] = w
    return np.conj(np.fft.fft(dense))              # nu_hat(eta)=sum_X nu(X)e(+2pi i eta X/M)


def R_of_d(nh, M, d):
    idx = (np.arange(M) * pow(2, d, M)) % M
    return np.sum(nh * np.conj(nh[idx]))


def main():
    print("# PROBE R28 -- THE GAP-CORRELATION OPERATOR (R27-E).\n")
    hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                        'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
    EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in hist.items()}
    Sbank = {k: float(F(7, 15) + EPS[k]) for k in EPS}
    Lam = {r: (EPS[r + 1] - EPS[r]) / 2 for r in range(1, 8)}

    lam = 0.5
    P = lambda d: (1 - lam) * lam ** abs(d) / (1 + lam)
    nus = build_nu(lam, 7)
    NH = {r: nu_hat(nus[r], 3 ** (r + 1)) for r in range(0, 8)}

    # ================= R28-A =================
    # normalization: my R_{r-1}(0)=3 X_{r-1}, so the closing identity is X_r = sum_d P(d) R_{r-1}(d)
    # (= Wilson's S_r=3 sum P(d)R with his R = mine/3). Increment X_r-X_{r-1} = S_r.
    Xb = {0: 1.0}
    for k in range(1, 8):
        Xb[k] = Xb[k - 1] + Sbank[k]
    print("## R28-A  GATE: X_r = sum_d P(d) R_{r-1}(d)  (lam=1/2); odd-d vanish exact; d=0 diag = X_{r-1}; increments = S_r")
    print(f"   {'r':>2} {'sum P(d)R(d)':>14} {'X_r bank':>10} {'ratio':>9} {'P(0)R(0)=diag':>14} {'X_{r-1}':>9} {'incr':>9} {'S_r':>9} {'max|odd R|':>11}")
    okA = True
    prev = None
    for r in range(2, 8):
        M1 = 3 ** r; nh = NH[r - 1]
        R = {d: R_of_d(nh, M1, d) for d in range(-40, 41)}
        odd = max(abs(R[d]) for d in range(-39, 40, 2))
        tot = sum(P(d) * R[d].real for d in range(-40, 41))
        diag = P(0) * R[0].real
        incr = (tot - prev) if prev is not None else float('nan')
        ok = abs(tot / Xb[r] - 1) < 1e-6 and odd < 1e-10 and abs(diag - Xb[r - 1]) < 1e-6
        okA = okA and ok
        print(f"   {r:>2} {tot:>14.6f} {Xb[r]:>10.6f} {tot/Xb[r]:>9.6f} {diag:>14.6f} {Xb[r-1]:>9.6f} "
              f"{incr:>9.6f} {Sbank[r]:>9.6f} {odd:>11.2e}")
        prev = tot
    print(f"   => R28-A {'GATE PASS -- operator closes exact, odd-d killed exact, diagonal flatness; R27-E CERTIFIED' if okA else 'FAIL (#43)'}")
    r = 5; nh = NH[r - 1]; M1 = 3 ** r
    print(f"   R_{r-1}(d) d=-4..4 (odd killed): " + "  ".join(f"{d}:{R_of_d(nh,M1,d).real:+.4f}" for d in range(-4, 5)) + "\n")

    # ================= R28-B =================
    print("## R28-B  kappa MEASURED: |m=+-1 twisted channel| / |untwisted|, per r (lam=1/2 and subcritical)")
    for lm in (0.5, 0.55, 0.6):
        Pl = lambda d: (1 - lm) * lm ** abs(d) / (1 + lm)
        nusl = build_nu(lm, 7); NHl = {rr: nu_hat(nusl[rr], 3 ** (rr + 1)) for rr in range(0, 8)}
        row = []
        for r in range(3, 8):
            M1 = 3 ** r; nh = NHl[r - 1]
            R0 = R_of_d(nh, M1, 0).real
            R2 = R_of_d(nh, M1, 2).real                # gap d=2 = m=+1 channel
            kap = abs(R2) / abs(R0) if R0 else float('nan')
            row.append(f"r{r}:{kap:.4f}")
        print(f"   lam={lm}: |R(2)/R(0)| per r = " + "  ".join(row) + f"   [2lam^2={2*lm*lm:.4f}]")
    print("   [Q: kappa<1 every r? stabilizing? (this is kappa as observable, not a fitted residual)]\n")

    # ================= R28-C =================
    print("## R28-C  CHANNEL MIXING: |R(4)|,|R(6)| (m=+-2,+-3) relative to |R(2)| (m=+-1), lam=1/2, per r")
    for r in range(3, 8):
        M1 = 3 ** r; nh = NH[r - 1]
        R2 = abs(R_of_d(nh, M1, 2)); R4 = abs(R_of_d(nh, M1, 4)); R6 = abs(R_of_d(nh, M1, 6))
        print(f"   r={r}: |R4/R2|={R4/R2:.4f}  |R6/R2|={R6/R2:.4f}")
    print("   [small => bound |lam2|/rho<2lam^2 nearly rigorous; large => channel-mixing gap, measured.]\n")

    # ================= R28-D =================
    print("## R28-D  REAL+PAIR FIT: deflate Lambda by fixed real mode 0.5, fit complex pair; signs & period")
    L = {r: float(Lam[r]) for r in range(1, 8)}
    mu = {r: L[r + 1] - 0.5 * L[r] for r in range(1, 7)}          # (E-0.5)Lambda removes the 0.5 mode
    print("   mu_r = Lambda_{r+1}-0.5 Lambda_r : " + "  ".join(f"{mu[r]:+.3e}" for r in range(1, 7)))
    # 2-term recurrence on mu (the complex pair), robust solves at r=2,3,4
    for i in (2, 3, 4):
        det = mu[i] * mu[i] - mu[i - 1] * mu[i + 1]
        if abs(det) < 1e-20:
            print(f"   solve{{{i},{i+1}}}: near-singular"); continue
        a = (mu[i + 1] * mu[i] - mu[i + 2] * mu[i - 1]) / det
        b = (mu[i] * mu[i + 2] - mu[i + 1] * mu[i + 1]) / det
        disc = a * a + 4 * b
        if disc < 0:
            z = complex(a / 2, math.sqrt(-disc) / 2)
            per = 2 * math.pi / abs(cmath.phase(z))
            print(f"   solve{{{i},{i+1}}}: a={a:+.4f} b={b:+.4f} -> |z|={abs(z):.4f} period={per:.3f} (complex pair)")
        else:
            r1 = (a + math.sqrt(disc)) / 2; r2 = (a - math.sqrt(disc)) / 2
            print(f"   solve{{{i},{i+1}}}: a={a:+.4f} b={b:+.4f} -> real roots {r1:+.4f},{r2:+.4f}")
    # full 3-mode reconstruction: roots 0.5, z, zbar from best solve; propagate & check signs
    print("   [if mu gives a clean complex pair, period resolves & 3-mode {0.5,z,zbar} should reproduce signs -,-,+,+,+,-,+]")


if __name__ == "__main__":
    main()
