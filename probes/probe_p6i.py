"""
PROBE P6I (Wilson) -- the deparitied residual rate: does Lambda_i settle, and to what S_inf? (2026-07-26)

Wilson found the factor 4: R_e-hat = C_rho-hat/(17-8cos phi), Re w = (4cos-1)/(17-8cos), denominators cancel, so
the (4,-1) kernel IS Re w's NUMERATOR and nu_e absorbs its denominator. T_i = S_{i+1}/2 exactly (P6H). The ONE open
number: the decay rate of the residual Lambda_i = T_i - T_{i-1}. Raw endpoint (0.977 over i=8..12) ignores the 2-cycle
wobble -- HARD RULE: deparity / two-step before quoting any rate.

To deparity honestly needs more than 6 points. So: VECTORIZE build_base2 (same certified SINGLEREC one-step, numpy
bincount instead of a Python dict loop) -- GATE bit-for-bit vs the original -- then push i to ~15-16 and deparity.

Reuses probe_27.stationary_trunc + gates vs probe_p6d.build_base2. No new transport (refactor, bit-for-bit gated).
"""
import os, sys, time
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_27_high_k_rho_q5 import stationary_trunc
from probe_p6d import build_base2 as build_base2_slow


def build_base2_fast(n, Amax=80):
    q = 3 ** (n + 1); Nn = 3 ** n; twoN = 2 * Nn
    D2 = np.full(q, -1, dtype=np.int64)
    g = 1
    for t in range(twoN):                      # base-2 dlog table (inherently sequential)
        D2[g] = t; g = (g * 2) % q
    inv2 = pow(2, -1, q)
    piW, _ = stationary_trunc(3, n)
    r = np.arange(Nn); cp = r[r % 3 != 0]
    nu = np.asarray(piW, float); nu = nu / nu.sum()
    Y = (3 * cp + 1) % q
    R_e = np.zeros(twoN); R_o = np.zeros(twoN); B = np.zeros(twoN)
    p = 1
    for a in range(1, Amax):
        p = (p * inv2) % q
        t = D2[(Y * p) % q]
        contrib = np.bincount(t, weights=nu * (0.5 ** a), minlength=twoN)
        if a % 2 == 0:
            R_e += contrib
        else:
            R_o += contrib
            if a == 1:
                B += contrib
    return dict(q=q, Nn=Nn, twoN=twoN, R_e=R_e, R_o=R_o, B=B)


def autocorr(f):
    F = np.fft.fft(f); return np.fft.ifft(F * np.conj(F)).real


def main():
    t0 = time.time()
    print("# PROBE P6I -- deparitied residual rate\n")

    # ---- GATE: fast == slow, bit-for-bit ----
    print("## GATE build_base2_fast == build_base2 (bit-for-bit, i=2..7)")
    for i in range(2, 8):
        Sf = build_base2_fast(i); Ss = build_base2_slow(i)
        e = max(np.max(np.abs(Sf['R_e'] - Ss['R_e'])), np.max(np.abs(Sf['R_o'] - Ss['R_o'])))
        print(f"   i={i}: max|fast-slow| = {e:.1e}")
    print()

    # ---- push deep ----
    IMAX = 15
    R0 = {}; R2 = {}
    for i in range(1, IMAX + 1):
        ti = time.time()
        S = build_base2_fast(i); Re = autocorr(S['R_e'])
        R0[i] = float(Re[0]); R2[i] = float(Re[2])
        if i >= 11:
            print(f"   [build i={i}: {time.time()-ti:.1f}s]")
    T = {i: 3 ** i * (4 * R2[i] - R0[i]) for i in range(1, IMAX + 1)}
    T[0] = 1.0 / 3
    Lam = {i: T[i] - T[i - 1] for i in range(1, IMAX + 1)}

    print(f"\n## the residual Lambda_i = T_i - T_(i-1), T_i = S_(i+1)/2  (i=1..{IMAX})")
    print(f"   {'i':>2} {'T_i':>12} {'Lambda_i':>13} {'T_i-7/30':>12}")
    for i in range(1, IMAX + 1):
        print(f"   {i:>2} {T[i]:>12.8f} {Lam[i]:>+13.8f} {T[i]-7/30:>+12.8f}")

    # ---- DEPARITY: two-step rate (removes the 2-cycle wobble) + even/odd subsequences ----
    print("\n## DEPARITIED rate of the residual (HARD RULE: two-step, not raw consecutive)")
    print(f"   {'i':>2} {'Lam_i':>12} {'raw Lam_i/Lam_(i-1)':>19} {'two-step (Lam_i/Lam_(i-2))^(1/2)':>32}")
    for i in range(8, IMAX + 1):
        raw = Lam[i] / Lam[i - 1] if abs(Lam[i - 1]) > 1e-18 else float('nan')
        two = (Lam[i] / Lam[i - 2]) ** 0.5 if Lam[i - 2] > 0 and Lam[i] > 0 else float('nan')
        print(f"   {i:>2} {Lam[i]:>+12.8f} {raw:>19.4f} {two:>32.5f}")
    # even/odd i subsequence geometric rates (deparitied by construction)
    ev = [i for i in range(8, IMAX + 1) if i % 2 == 0]
    od = [i for i in range(8, IMAX + 1) if i % 2 == 1]
    def subrate(idx):
        v = np.array([Lam[i] for i in idx])
        if len(v) < 2 or np.any(v <= 0):
            return float('nan'), v
        # geometric rate per single level from the even/odd chain (step 2 -> sqrt)
        rr = (v[1:] / v[:-1]) ** 0.5
        return rr, v
    re, ve = subrate(ev); ro, vo = subrate(od)
    print(f"   even-i {ev}: Lam={np.array([Lam[i] for i in ev])}")
    print(f"       per-level rate (sqrt of 2-step): {re}")
    print(f"   odd-i  {od}: Lam={np.array([Lam[i] for i in od])}")
    print(f"       per-level rate (sqrt of 2-step): {ro}")

    # ---- implied S_inf under a geometric tail at rate rho (deparitied), from i=IMAX ----
    print("\n## implied S_inf = 2*T_inf, T_inf = T_IMAX + Sum_{j>IMAX} Lambda_j, geometric tail Lambda~Lam_IMAX*rho^(j-IMAX)")
    for rho in (0.90, 0.93, 0.95, 0.977, 1.0):
        if rho < 1:
            tail = Lam[IMAX] * rho / (1 - rho)
            Sinf = 2 * (T[IMAX] + tail)
            print(f"   rho={rho:.3f}: tail Sum={tail:+.5f}  S_inf={Sinf:.5f}  (7/15={7/15:.5f}, 0.477 ref)")
        else:
            print(f"   rho=1.000: tail DIVERGES -> S_inf unbounded (residual must decay)")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
