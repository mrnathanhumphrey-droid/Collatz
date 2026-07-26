"""
PROBE CHANNEL_ID -- verification of the pen identity:  d1_r = A_r(1)/A_r(0)   (r >= 2).

Derivation (pen, 2026-07-25): the three target lags {1, m-1, m+1} (m=3^{r-1}) are the SAME base lag 1 in the
three fiber channels (m-1 = 1+2m mod 3m, C even). For fixed s, {s+1+dm : d=0,1,2} is exactly the full fiber over
(s+1 mod m) [no carries: the 3 points differ by m, 3m=0], and the tower gives sum_j rho_r(t+jm) = rho_{r-1}(t)
EXACTLY (nu_r mod 3^r = nu_{r-1}; dlog commutes with reduction). Hence
    sum_d C^{(r)}(k+dm) = C^{(r-1)}(k)   for any k, so
    Num = 2C(1)-C(m-1)-C(m+1) = 3 p_r(1) - p_{r-1}(1),   Den = 2[C(0)-C(m)] = 3 p_r(0) - p_{r-1}(0),
    d1_r = A_r(1)/A_r(0),  A_r(m') = gamma_r(m') - gamma_{r-1}(m'),  gamma_r = 3^r p_r  (banked charledger objects).
Corollary (same proof): Re dhat_r(n) = A_r(n)/A_r(0) for all 3 nmid n -- the WHOLE mode ladder = channel ledger.

Equivalences: d1>0 <=> A_r(1)>0 (m=1 channel grows; MOON monotone+) <=> q_r(1) = p_r(1)/p_{r-1}(1) > 1/3
(conditional lifting of the x4-collision beats neutral; = S2's measured survival ~0.3335, m=1, thru r~38).

VERIFY (falsifiers, in order):
 V1 exact (Fractions, r=2..5): tower fold rho_r -> rho_{r-1} EXACT equality; identity vs banked exact d1.
 V2 float (r=2..16): identity vs banked d1 ladder (build_nu r<=11; scratchpad rho_12..16).
 V3 corollary n=2 vs banked d2 (r=12..16).
 V4 the Lambda weld: Lam_r = sum_{m>=1} 4^{-m} A_r(m) vs banked LAM_NU (r=12..16) -- ties to the eps ladder.
 V5 the q-ladder: q_r(1), 3q_r-1 for r=2..16 (the exact-side version of S2's survival plateau).
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu
from probe_ratio2 import build_nu_exact

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"
D1 = {2: 5.714285714286e-2, 3: 1.982809084493e-2, 4: 8.643329308003e-3, 5: 8.039124362558e-3,
      6: 7.742271269236e-3, 7: 5.786625297669e-3, 8: 4.665965059268e-3, 9: 3.894809120578e-3,
      10: 3.600723714723e-3, 11: 3.211120395286e-3, 12: 2.963565168845e-3, 13: 2.696227117865e-3,
      14: 2.440864723362e-3, 15: 2.187132833772e-3, 16: 1.939224678822e-3}
D2 = {12: -1.851543e-4, 13: +7.733991e-5, 14: -3.950496e-5, 15: +1.516829e-4, 16: +1.568780e-4}
LAM_NU = {12: 3.3677e-4, 13: 3.1971e-4, 14: 2.8672e-4, 15: 2.6193e-4, 16: 2.3426e-4}
D1EX = {2: F(2, 35), 3: None, 4: None, 5: None}   # r=2 exact known 4/70; others compared as floats


def rho_exact_norm(nu, r):
    N = 3 ** r; d = R10.dlog_table(r)
    rho = {}
    for X, w in nu.items():
        s = d[(X - 1) // 3 % N]
        rho[s] = rho.get(s, F(0)) + w
    tot = sum(rho.values())
    return {s: w / tot for s, w in rho.items()}, N


def C_ex(rho, N, k):
    return sum(w * rho.get((s + k) % N, F(0)) for s, w in rho.items())


def main():
    t0 = time.time()
    print("# PROBE CHANNEL_ID -- verify d1_r = A_r(1)/A_r(0)\n")

    # ---------------- V1 exact ----------------
    print("## V1  EXACT (Fractions), r=2..5: tower fold + identity")
    nex = build_nu_exact(5)
    rex = {}
    for r in range(1, 6):
        rex[r], _ = rho_exact_norm(nex[r], r)
    for r in range(2, 6):
        N = 3 ** r; m = N // 3
        fold = {}
        for s, w in rex[r].items():
            fold[s % m] = fold.get(s % m, F(0)) + w
        ok = all(fold.get(t, F(0)) == rex[r - 1].get(t, F(0)) for t in set(fold) | set(rex[r - 1]))
        # identity: d1 = (3 p_r(1) - p_{r-1}(1)) / (3 p_r(0) - p_{r-1}(0))
        pr1, pr0 = C_ex(rex[r], N, 1), C_ex(rex[r], N, 0)
        pm1, pm0 = C_ex(rex[r - 1], N // 3, 1), C_ex(rex[r - 1], N // 3, 0)
        d1id = (3 * pr1 - pm1) / (3 * pr0 - pm0)
        # direct 5-lag d1 for cross-check
        d1direct = (2 * C_ex(rex[r], N, 1) - C_ex(rex[r], N, m - 1) - C_ex(rex[r], N, m + 1)) / (2 * (pr0 - C_ex(rex[r], N, m)))
        same = (d1id == d1direct)
        rel = abs(float(d1id) - D1[r]) / D1[r]
        extra = f"  == 2/35 exact: {d1id == F(2,35)}" if r == 2 else ""
        print(f"   r={r}: tower fold EXACT: {ok} | identity == 5-lag (exact Fraction): {same} | vs banked d1 rel {rel:.1e}{extra}")
    print()

    # ---------------- float rho, r=2..16 ----------------
    nus = build_nu(0.5, 11)
    rho = {}
    for r in range(1, 12):
        N = 3 ** r
        mu = np.zeros(N)
        for X, w in nus[r].items():
            mu[(X - 1) // 3 % N] += float(w)
        d = R10.dlog_table(r)
        rr = np.zeros(N); rr[np.fromiter((d[a] for a in range(N)), np.int64, N)] = mu
        rho[r] = rr / rr.sum()
    del nus
    for r in range(12, 17):
        rr = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy"))
        rho[r] = rr / rr.sum()
    print(f"  (rho 1..16 ready, {time.time()-t0:.1f}s)\n")

    MM = 13
    p = {0: {k: 1.0 for k in range(MM + 1)}}          # level 0: trivial group, nu_0=delta, p=1 for all lags
    for r in range(1, 17):
        p[r] = {k: float(np.dot(rho[r], np.roll(rho[r], -k))) for k in range(MM + 1)}

    # ---------------- V2 identity vs banked d1 ----------------
    print("## V2  identity d1 = (3p_r(1)-p_{r-1}(1))/(3p_r(0)-p_{r-1}(0)) vs banked, r=2..16")
    worst = 0
    for r in range(2, 17):
        d1id = (3 * p[r][1] - p[r - 1][1]) / (3 * p[r][0] - p[r - 1][0])
        rel = abs(d1id - D1[r]) / D1[r]; worst = max(worst, rel)
        print(f"   r={r:>2}: {d1id:+.9e} vs {D1[r]:+.9e}  rel {rel:.1e}")
    print(f"   worst rel = {worst:.2e}  [{'IDENTITY VERIFIED' if worst < 1e-6 else 'FAIL'}]\n")

    # ---------------- V3 corollary n=2 ----------------
    print("## V3  corollary Re dhat(2) = A_r(2)/A_r(0) vs banked d2, r=12..16")
    for r in range(12, 17):
        d2id = (3 * p[r][2] - p[r - 1][2]) / (3 * p[r][0] - p[r - 1][0])
        print(f"   r={r}: {d2id:+.6e} vs banked {D2[r]:+.6e}  rel {abs(d2id-D2[r])/abs(D2[r]):.1e}")
    print()

    # ---------------- V4 Lambda weld ----------------
    print("## V4  Lambda weld: sum_{m=1..%d} 4^-m A_r(m) vs banked LAM_NU (nu-route Lambda)" % MM)
    for r in range(12, 17):
        A = {m: (3 ** r) * p[r][m] - (3 ** (r - 1)) * p[r - 1][m] for m in range(MM + 1)}
        lam = sum(4.0 ** -m * A[m] for m in range(1, MM + 1))
        print(f"   r={r}: {lam:+.5e} vs banked {LAM_NU[r]:+.5e}  rel {abs(lam-LAM_NU[r])/LAM_NU[r]:.1e}")
    print()

    # ---------------- V5 q-ladder ----------------
    print("## V5  the q-ladder (exact-side version of S2's survival plateau): q_r(1)=p_r(1)/p_{r-1}(1)")
    print(f"   {'r':>2} {'gamma_r(1)':>11} {'A_r(1)':>11} {'q_r(1)':>9} {'3q-1':>10} | {'A_r(0)':>9} {'X_r=gam(0)':>11}")
    for r in range(1, 17):
        g1 = (3 ** r) * p[r][1]; g1m = (3 ** (r - 1)) * p[r - 1][1] if r >= 1 else float('nan')
        A1 = g1 - g1m
        q = p[r][1] / p[r - 1][1]
        g0 = (3 ** r) * p[r][0]; A0 = g0 - (3 ** (r - 1)) * p[r - 1][0]
        print(f"   {r:>2} {g1:>11.6f} {A1:>+11.2e} {q:>9.6f} {3*q-1:>+10.2e} | {A0:>9.5f} {g0:>11.4f}")
    print("   [gamma_0(1)=1, gamma_1(1)=2/3 (hand-exact: nu_1=(0,1/3,2/3) in dlog): the channel FALLS first,")
    print("    then A_r(1)>0 for all r>=2 = d1>0. q_16 ~ 0.33375; S2 measured q~0.3335 (m=1) thru r~38.]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
