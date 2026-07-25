"""
PROBE RATIO-2 -- is the oscillating component of Re dhat_r(1) SUBDOMINANT (dying vs the dominant mode)?

Wilson's spec (2026-07-25). RATIO settled that rho_r = d1_{r+1}/d1_r oscillates. RATIO-2 asks the only version
that touches the sign: does the oscillation DECAY relative to the dominant geometric mode?

Primary observable (fit-free): D_r := ln rho_r - ln rho_{r-1} = ln d1_{r+1} - 2 ln d1_r + ln d1_{r-1}
  = 2nd difference of ln d1_r; annihilates the dominant geometric rate EXACTLY. D_r IS the oscillating component.

Two independent computations of d1_r = Re dhat_r(1):
  (P) profile path (validated G0): rho=dlog pushforward of nu; prof=|FFT rho|^2; delta=prof/S-1/M on prim;
      d_n = FFT(delta)[n].real.  Gives d1,d2,d3.
  (F) EXACT closed form (derived): d1 = [2C(1) - (C(N/3-1)+C(N/3+1))] / [2(C(0)-C(N/3))], N=3^r,
      C(d)=sum_s rho(s)rho((s+d)%N).  (sum_{k in prim} cos(2pi k a/N)=N[N|a]-(N/3)[(N/3)|a] kills all but 5 lags.)
  (F) is available in EXACT rationals (reference) and float (cross-check). P and F must agree.

Run order (Wilson): R2-A precision-floor GATE first; only if it passes, R2-B/C/D.
Predictions (pre-registered): subdominant -> D_14=-0.0052 in [-0.0104,-0.0026], signs(D12..D15)=+,+,-,+, dip r=14.
                              not-subdominant -> |D_14| >= 0.0363.
"""
import os, sys, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu

# ---------------- dlog as an int array (one multiplicative pass; avoids the dict + comprehension) --------
def dlog_array(r):
    Mp = 3 ** (r + 1); N = 3 ** r
    arr = np.empty(N, dtype=np.int64); el = 1
    for t in range(N):
        arr[((el - 1) // 3) % N] = t
        el = (el * 4) % Mp
    return arr

def rho_dense(nu, j):
    N = 3 ** j
    mu = np.zeros(N)
    for X, w in nu.items():
        mu[(X - 1) // 3 % N] += float(w)
    arr = dlog_array(j)
    rho = np.zeros(N)
    np.add.at(rho, arr, mu)
    del arr, mu
    return rho, N

# ---------------- (P) profile path: d1,d2,d3 = Re dhat(1,2,3) ----------------
def d_hat_123(rho, N):
    prof = np.abs(np.fft.fft(rho)) ** 2
    k = np.arange(N)
    primmask = (k % 3 != 0)                     # k=0 excluded (3|0); k in 1..N-1, 3 nmid k
    M = int(primmask.sum())
    S = float(prof[primmask].sum())
    dd = np.zeros(N)
    dd[primmask] = prof[primmask] / S - 1.0 / M
    del prof
    H = np.fft.fft(dd)
    return H[1].real, H[2].real, H[3].real

# ---------------- (F) closed form, float ----------------
def d1_formula_float(rho, N):
    n3 = N // 3
    C = lambda dl: float(np.dot(rho, np.roll(rho, -dl)))
    return (2 * C(1) - (C(n3 - 1) + C(n3 + 1))) / (2 * (C(0) - C(n3)))

# ---------------- (F) closed form, EXACT (Fraction) ----------------
def build_nu_exact(RMAX, vmax=120):
    nus = {0: {1: F(1)}}; nu = nus[0]
    for r in range(1, RMAX + 1):
        M = 3 ** (r + 1); Mr = 3 ** r; inv2 = pow(2, -1, Mr)
        new = {}
        for Xp, wp in nu.items():
            u = inv2
            for v in range(1, vmax + 1):
                Xr = (1 + 3 * ((u * Xp) % Mr)) % M
                new[Xr] = new.get(Xr, F(0)) + F(1, 2 ** v) * wp
                u = (u * inv2) % Mr
        nu = new; nus[r] = nu
    return nus

def rho_exact_dict(nu, j):
    N = 3 ** j; mu = {}
    for X, w in nu.items():
        a = (X - 1) // 3 % N; mu[a] = mu.get(a, F(0)) + w
    d = R10.dlog_table(j); rho = {}
    for a, m in mu.items():
        t = d[a]; rho[t] = rho.get(t, F(0)) + m
    return rho, N

def d1_exact(rho, N):
    n3 = N // 3
    def C(dl):
        s = F(0)
        for k, val in rho.items():
            w = rho.get((k + dl) % N)
            if w is not None:
                s += val * w
        return s
    return (2 * C(1) - (C(n3 - 1) + C(n3 + 1))) / (2 * (C(0) - C(n3)))


def main():
    t0 = time.time()
    print("# PROBE RATIO-2 -- is the oscillating component of Re dhat_r(1) subdominant?\n")

    # ===================== R2-A  PRECISION FLOOR (gate) =====================
    print("## R2-A  PRECISION FLOOR (gate; signal |D_14|~0.0052, need floor < 0.0026)\n")
    nprod = build_nu(0.5, 13)

    print("### A0  profile-path (P) vs closed-form (F), float, r=2..12  [both must agree]")
    d1P = {}; d1F = {}; d2 = {}; d3 = {}
    for r in range(2, 13):
        rho, N = rho_dense(nprod[r], r)
        a1, a2, a3 = d_hat_123(rho, N)
        d1P[r] = a1; d2[r] = a2; d3[r] = a3
        d1F[r] = d1_formula_float(rho, N)
        del rho
    worst = max(abs(d1P[r] - d1F[r]) / abs(d1P[r]) for r in range(2, 13))
    print(f"   max rel|P - F| over r=2..12 = {worst:.2e}  [{'OK both agree' if worst < 1e-9 else 'FAIL'}]")
    print("   d1_r (P): " + " ".join(f"{d1P[r]:+.4e}" for r in range(2, 13)) + "\n")

    print("### A1  exact-rational reference (r<=7):  |float(P) - exact| / |exact|")
    nex = build_nu_exact(7); sig = {}
    for r in range(2, 8):
        rho_e, N = rho_exact_dict(nex[r], r)
        de = d1_exact(rho_e, N)
        sig[r] = abs(d1P[r] - float(de)) / abs(float(de))
        print(f"   r={r}: exact={float(de):+.12e}  float(P)={d1P[r]:+.12e}  rel={sig[r]:.2e}")
    print(f"   sigma_rel r=2..7: " + " ".join(f"{sig[r]:.1e}" for r in range(2, 8)) + "\n")

    print("### A2  working-end floor (r=12,13)")
    ntight = build_nu(0.5, 13, tol=1e-26)
    for r in (12, 13):
        rp, N = rho_dense(nprod[r], r); dp = d_hat_123(rp, N)[0]
        rt, _ = rho_dense(ntight[r], r); dt = d_hat_123(rt, N)[0]
        fm = d1_formula_float(rp, N)
        del rp, rt
        print(f"   r={r}: tol(1e-18 vs 1e-26) rel={abs(dp-dt)/abs(dp):.2e} ; P-vs-F rel={abs(dp-fm)/abs(dp):.2e}")
    rp, N = rho_dense(nprod[13], 13); dp = d_hat_123(rp, N)[0]
    rt, _ = rho_dense(ntight[13], 13); dt = d_hat_123(rt, N)[0]
    sig_work = max(abs(dp - dt) / abs(dp), abs(dp - d1_formula_float(rp, N)) / abs(dp), sig[7])
    del rp, rt
    dD = 4 * sig_work
    print(f"\n   working sigma_rel = {sig_work:.2e}  =>  floor delta_D ~ 4*sigma = {dD:.2e}  vs half-signal 0.0026")
    if dD >= 0.0026:
        print(f"\n   *** GATE FAIL: floor {dD:.2e} >= 0.0026. INCONCLUSIVE, stopping before R2-B. ***")
        print(f"# ({time.time()-t0:.1f}s)"); return
    print(f"   *** GATE PASS: floor {dD:.2e} << 0.0026. proceeding. ***\n")

    # ===================== R2-B  EXTEND THE LADDER =====================
    print("## R2-B  extend d1_r (and d2,d3) to r=13,14,15,(16)")
    d1 = dict(d1P)
    del nprod, ntight
    RTOP = 15
    try:
        ndeep = build_nu(0.5, 16); RTOP = 16
        print(f"   build_nu reached r=16 ({time.time()-t0:.1f}s)")
    except MemoryError:
        ndeep = build_nu(0.5, 15); RTOP = 15
        print(f"   build_nu r=16 OOM -> r=15 ({time.time()-t0:.1f}s)")
    for r in range(13, RTOP + 1):
        try:
            rho, N = rho_dense(ndeep[r], r)
            a1, a2, a3 = d_hat_123(rho, N)
            fm = d1_formula_float(rho, N); del rho
            d1[r] = a1; d2[r] = a2; d3[r] = a3
            print(f"   r={r}: d1={a1:+.12e}  (P-vs-F rel {abs(a1-fm)/abs(a1):.1e})  d2={a2:+.4e}  d3={a3:+.4e}  ({time.time()-t0:.1f}s)")
        except MemoryError:
            print(f"   r={r}: OOM -> ladder stops at r={r-1}"); RTOP = r - 1; break
    print("\n### full d1_r table (auditable)")
    for r in sorted(d1):
        print(f"   r={r:>2}  d1={d1[r]:+.12e}   d2={d2.get(r, float('nan')):+.6e}   d3={d3.get(r, float('nan')):+.6e}")
    print()

    # ===================== R2-C  DIP TEST =====================
    print("## R2-C  rho_r, D_r, dip test")
    rr = {r: d1[r + 1] / d1[r] for r in sorted(d1) if r + 1 in d1}
    D = {r: math.log(rr[r]) - math.log(rr[r - 1]) for r in sorted(rr) if r - 1 in rr}
    dr = {r: rr[r] - rr[r - 1] for r in sorted(rr) if r - 1 in rr}
    print(f"   {'r':>3} {'rho_r':>10} {'D_r':>11} {'d_rho':>11}")
    for r in sorted(rr):
        print(f"   {r:>3} {rr[r]:>10.4f} {(f'{D[r]:+.4f}' if r in D else ''):>11} {(f'{dr[r]:+.4f}' if r in dr else ''):>11}")
    print("\n   PRED: subdominant D_14=-0.0052 [-0.0104,-0.0026]; not-sub |D_14|>=0.0363; signs(D12..15)=+,+,-,+; dip r=14")
    for r in (12, 13, 14, 15):
        if r in D:
            print(f"      D_{r} = {D[r]:+.5f}   drho_{r} = {dr[r]:+.5f}")
    if 14 in D and 10 in D and D[14] / D[10] > 0:
        mu = (D[14] / D[10]) ** 0.25; rho1 = rr.get(11, rr[max(rr)])
        print(f"\n   mu = (D_14/D_10)^(1/4) = {mu:.4f} ;  rho_c = mu*rho1 = {mu*rho1:.4f}  (rho1~{rho1:.3f})")
    if 14 in D:
        d14 = D[14]
        v = ("SUBDOMINANT/DYING (in band)" if -0.0104 <= d14 <= -0.0026 else
             "NOT SUBDOMINANT (|D14|>=0.0363) -> deep tail oscillates, 7/15 live" if abs(d14) >= 0.0363 else
             "over-damped/near-floor -- report mu,rho_c, no verdict" if abs(d14) < 0.0026 else
             "INTERMEDIATE (decaying but slower than predicted) -- report mu,rho_c, no verdict")
        print(f"\n   READ D_14={d14:+.5f}: {v}")

    # ===================== R2-D  SHAPE =====================
    print("\n## R2-D  per-period mean of D (is rho converging?)")
    pm = lambda lo, hi: (lambda vs: (sum(vs) / len(vs) if vs else float('nan'), len(vs)))([D[r] for r in range(lo, hi + 1) if r in D])
    m1, n1 = pm(7, 10); m2, n2 = pm(11, 14)
    print(f"   mean D r=7..10  = {m1:+.4f} (n={n1})  [banked +0.0444]")
    print(f"   mean D r=11..14 = {m2:+.4f} (n={n2})")
    if not math.isnan(m2):
        print("   => " + ("trend SHRINKING: 3:1 asymmetry = trend+oscillation, Aitken 0.908 meaningful."
                          if abs(m2) < 0.5 * abs(m1) else
                          "trend FLAT ~0.044: rho NOT converging, 0.908 premature, tail factor 10.75 understated."))
    print(f"\n# ({time.time()-t0:.1f}s)  RTOP={RTOP}")


if __name__ == "__main__":
    main()
