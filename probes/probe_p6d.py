"""
PROBE P6D (Wilson) -- the cross-term COLLAPSE in the certified base-2 coordinate, boundary term carried (2026-07-26).

Wilson certified the reindex: 2 is a primitive root mod 3^{n+1}, so dlog_2 : units -> Z/(2*3^n) is a bijection with
x2 = +1, and CRT Z/(2*3^n) = Z/2 (the <4>-coset = branch parity = x mod 3) x Z/3^n (base-4 dlog). In THIS coordinate
x2^-1 is a plain shift by -1. P6B's roll-scan failed only because it lived in the FOLDED base-4 coordinate where
x2^-1 flips the coset (not a roll). Here: no fold.

The collapse is an EXACT algebraic identity (not an approximate roll). Since x2^-1 . push_a = push_{a+1} exactly
(one more division by 2), and nu_o = Sum_{a odd>=1} 2^-a push_a, nu_e = Sum_{a even>=2} 2^-a push_a (a=0 excluded):
    1/2 (x2^-1)_* nu_e = 1/2 Sum_{a even>=2} 2^-a push_{a+1} = Sum_{b odd>=3} 2^-b push_b = nu_o - 2^-1 push_1(nu).
=>  nu_o = 1/2 (x2^-1)_* nu_e  +  2^-1 (m_1)_* nu       [Wilson's boundary term, EXACT, carried explicitly].

BUILD: certified SINGLEREC one-step in residue domain, indexed by base-2 dlog (unfolded). For coprime x mod 3^n,
Y=3x+1 mod q=3^{n+1} (a unit; 2 primitive => dlog_2 defined). x' = Y*2^-a mod q; a even -> R_e, a odd -> R_o,
a==1 -> boundary B. All exact.

GATE-1 (reindex/stationarity): reduce R=R_e+R_o mod 3^n == the stationary measure nuW (fixed point).
GATE-2 (the collapse, machine precision): R_o[t] == 1/2 R_e[(t+1) mod 2*3^n] + B[t].  Report WITH and WITHOUT B.

Reuses stationary_trunc (certified forward measure). No new transport -- one certified Syracuse step, reindexed.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import gcd
import numpy as np
from probe_27_high_k_rho_q5 import stationary_trunc


def build_base2(n, Amax=80):
    q = 3 ** (n + 1); Nn = 3 ** n; twoN = 2 * Nn
    # base-2 dlog table over units mod q (2 is a primitive root => enumerates all 2*3^n units)
    D2 = {}; g = 1
    for t in range(twoN):
        D2[g] = t; g = (g * 2) % q
    assert len(D2) == twoN, f"dlog_2 not a bijection: {len(D2)} != {twoN}"
    inv2 = pow(2, -1, q)

    piW, _ = stationary_trunc(3, n)
    cp = np.array([r for r in range(Nn) if gcd(r, 3) == 1], dtype=np.int64)
    nu = np.asarray(piW, float); nu = nu / nu.sum()               # stationary measure over coprime mod 3^n
    nu_full = np.zeros(Nn); nu_full[cp] = nu

    R_e = np.zeros(twoN); R_o = np.zeros(twoN); B = np.zeros(twoN)  # boundary = a==1 (m_1) piece
    for x, w in zip(cp.tolist(), nu.tolist()):
        Y = (3 * x + 1) % q
        p = 1
        for a in range(1, Amax):
            p = (p * inv2) % q
            t = D2[(Y * p) % q]
            m = w * (0.5 ** a)
            if a % 2 == 0:
                R_e[t] += m
            else:
                R_o[t] += m
                if a == 1:
                    B[t] += m
    return dict(q=q, Nn=Nn, twoN=twoN, D2=D2, nu_full=nu_full, R_e=R_e, R_o=R_o, B=B)


def main():
    t0 = time.time()
    print("# PROBE P6D -- cross-term collapse in base-2 coordinate, boundary carried\n")
    print("Identity to verify (exact):  nu_o = 1/2 (x2^-1)_* nu_e + 2^-1 (m_1)_* nu\n")

    for n in (2, 3, 4, 5, 6):
        S = build_base2(n); twoN = S['twoN']; Nn = S['Nn']
        R_e, R_o, B = S['R_e'], S['R_o'], S['B']

        # GATE-1: R reduces mod 3^n to the stationary measure (fixed point of the certified step)
        R = R_e + R_o
        # invert dlog_2: position t -> residue (2^t mod q) mod 3^n
        red = np.zeros(Nn)
        g = 1; q = S['q']
        for t in range(twoN):
            red[g % Nn] += R[t]; g = (g * 2) % q
        gate1 = np.max(np.abs(red - S['nu_full']))

        # mass split (should be P(a even)=1/3, P(a odd)=2/3 -- the seed)
        me, mo = R_e.sum(), R_o.sum()

        # GATE-2: the collapse.  x2^-1 pushforward = shift dlog by -1 => (x2^-1 nu_e)[t] = R_e[t+1]
        shift_Re = np.roll(R_e, -1)                # shift_Re[t] = R_e[(t+1) mod twoN]
        pred_noB = 0.5 * shift_Re
        pred_B = 0.5 * shift_Re + B
        scale = np.max(np.abs(R_o)) + 1e-30
        res_noB = np.max(np.abs(R_o - pred_noB)) / scale
        res_B = np.max(np.abs(R_o - pred_B)) / scale
        # is the residual-without-B exactly B ?  (i.e. does dropping the boundary miss precisely 2^-1 push_1)
        res_isB = np.max(np.abs((R_o - pred_noB) - B)) / scale

        print(f"## n={n} (twoN=2*3^{n}={twoN})")
        print(f"   GATE-1 fixed-point  max|reduce(R)-nu| = {gate1:.2e}   "
              f"[mass a-even={me:.5f} (1/3?), a-odd={mo:.5f} (2/3?)]")
        print(f"   COLLAPSE  |R_o - 1/2 shift(R_e) - B| = {res_B:.2e}   "
              f"[{'EXACT' if res_B < 1e-12 else 'approx'}]")
        print(f"   without boundary B:      residual = {res_noB:.2e}   "
              f"(dropped term == B? {res_isB:.2e})")
        print()

    print("# interpretation: EXACT collapse => nu_o carries NO information beyond nu_e + the a=1 boundary.")
    print("# d_1 = A_j(1) = <cross-parity> then reduces to a functional of the single sub-measure rho_e (+boundary).")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
