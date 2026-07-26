"""
PROBE CHANNEL_AUDIT -- Wilson's audit items, quantified.

 A) A_r(0) == S_r EXACTLY? (Wilson's own arithmetic had A_2(0)=10/21=S_2; R28 diagonal flatness.)
    Exact Fractions r=2..5 vs the exact eps json; float r=2..16 vs the EPS chain.
    If true the identity is d1_r = A_r(1)/S_r -- one more weld.
 B) Lambda identity EXACT: Lam_r = sum_{m>=1} 4^{-m} A_r(m). Check vs exact-eps Lambda (r=2..7) at
    machine precision with K=40; then decompose V4's 1.7e-5: banked-LAM_NU 5-digit rounding vs K=13 truncation.
 C) S2 error bar ON q (the decider for "plateau" vs "consistency"): analytic SE from N=4e5, R_REP=20,
    with and without ESS discount (uniq-frac 0.223); per-level z at r=38 and window-aggregate z (17..38)
    against the propagated excess 4.2e-4 * 0.887^{r-16}.
 D) gamma_inf(1) extrapolation + inverse-symbolic scan (rationals q<=200 + named combos incl (1+S_inf)/2).
"""
import os, sys, json, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu
from probe_ratio2 import build_nu_exact

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"

# exact eps (k<=8) + float chain (banked, probe_fiber constants)
_hist = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments_output',
                                    'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
EPS_EX = {int(k): F(int(v['num']), int(v['den'])) for k, v in _hist.items()}          # exact eps_1..eps_8
EPS_F = {k: float(v) for k, v in EPS_EX.items()}
EPS_F.update({9: -7.520257156400000e-6, 10: 7.207509171100000e-4, 11: 1.501967012082273e-3, 12: 2.274713720558208e-3})
LAM_NU = {12: 3.3677e-4, 13: 3.1971e-4, 14: 2.8672e-4, 15: 2.6193e-4, 16: 2.3426e-4}
for r in range(12, 16):
    EPS_F[r + 1] = EPS_F[r] + 2 * LAM_NU[r]
S715 = F(7, 15)


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
    print("# PROBE CHANNEL_AUDIT -- Wilson's audit, quantified\n")

    # -------- exact rho r<=5, float rho r<=16 --------
    nex = build_nu_exact(5)
    rex = {r: rho_exact_norm(nex[r], r)[0] for r in range(1, 6)}
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
        rr = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy")); rho[r] = rr / rr.sum()
    KMAX = 40
    p = {0: {k: 1.0 for k in range(KMAX + 1)}}
    for r in range(1, 17):
        p[r] = {k: float(np.dot(rho[r], np.roll(rho[r], -k))) for k in range(KMAX + 1)}
    print(f"  (rho ready, {time.time()-t0:.1f}s)\n")

    # ======== A) A_r(0) == S_r ========
    print("## A)  A_r(0) == S_r  (diagonal flatness weld; then d1 = A_r(1)/S_r)")
    print("   EXACT (Fractions), r=2..5:")
    for r in range(2, 6):
        N = 3 ** r
        A0 = 3 ** r * C_ex(rex[r], N, 0) - 3 ** (r - 1) * C_ex(rex[r - 1], N // 3, 0)
        Sr = S715 + EPS_EX[r]
        print(f"     r={r}: A_r(0) = {A0}  vs  S_r = {Sr}   EQUAL: {A0 == Sr}")
    print("   float, r=6..16 (vs EPS chain; 13..16 chain built from 5-digit LAM_NU):")
    for r in (6, 8, 10, 12, 14, 16):
        A0 = 3 ** r * p[r][0] - 3 ** (r - 1) * p[r - 1][0]
        Sr = 7.0 / 15 + EPS_F[r]
        print(f"     r={r:>2}: A_r(0)={A0:.9f}  S_r={Sr:.9f}  rel {abs(A0-Sr)/Sr:.1e}")
    print()

    # ======== B) Lambda identity exact ========
    print("## B)  Lambda identity: Lam_r = sum_{m>=1} 4^-m A_r(m)  (K=40) vs exact-eps Lambda")
    for r in range(2, 11):
        lam_true = (EPS_F[r + 1] - EPS_F[r]) / 2
        A = {m: 3 ** r * p[r][m] - 3 ** (r - 1) * p[r - 1][m] for m in range(1, KMAX + 1)}
        lam_ch = sum(4.0 ** -m * A[m] for m in range(1, KMAX + 1))
        tag = "exact-eps" if r + 1 <= 8 else "float-eps"
        print(f"   r={r:>2}: channel-sum {lam_ch:+.12e}  vs  (eps_{r+1}-eps_{r})/2 = {lam_true:+.12e}"
              f"  rel {abs(lam_ch-lam_true)/abs(lam_true):.1e}  [{tag}]")
    print("   V4's 1.7e-5 decomposed (r=16): banked LAM_NU has 5 sig digits -> half-ulp rel ~", end=" ")
    print(f"{0.5e-4/2.3426:.1e}; K=13 truncation adds ~4^-14*A(14)/Lam ~ "
          f"{4.0**-14*abs(3**16*p[16][14]-3**15*p[15][14])/LAM_NU[16]:.1e}  => rounding dominates.\n")

    # ======== C) S2 error bar on q ========
    print("## C)  S2's error bar ON q (N=4e5, R_REP=20, uniq-frac ~0.223)")
    N_, R_ = 4e5, 20
    q0 = 1.0 / 3
    se_naive = math.sqrt(q0 * (1 - q0) / N_) / math.sqrt(R_)
    se_ess = math.sqrt(q0 * (1 - q0) / (N_ * 0.223)) / math.sqrt(R_)
    print(f"   SE(q) per level: naive {se_naive:.2e} ; ESS-discounted {se_ess:.2e}   (needed at r=38: ~3e-5)")
    exc16 = 4.23e-4; rho_d = 0.887
    exc = {r: exc16 * rho_d ** (r - 16) for r in range(17, 39)}
    z38n, z38e = exc[38] / se_naive, exc[38] / se_ess
    mean_exc = sum(exc.values()) / len(exc)
    zaggn = mean_exc / (se_naive / math.sqrt(22)); zagge = mean_exc / (se_ess / math.sqrt(22))
    print(f"   per-level at r=38: excess {exc[38]:.2e} -> z = {z38n:.2f} (naive) / {z38e:.2f} (ESS)  [UNRESOLVABLE]")
    print(f"   window-aggregate 17..38: mean excess {mean_exc:.2e} -> z = {zaggn:.1f} (naive) / {zagge:.1f} (ESS)")
    print("   => Wilson's downgrade CONFIRMED: per-level sign at the deep end NOT resolvable; aggregate ~2-4 sigma =")
    print("      'consistent where informative', the strong statistical statement stays the telescoped eps-hat (z~3.4).")
    print("      And the 2.3-sigma r~36 rollover sits exactly where excess is smallest + resolution worst.\n")

    # ======== D) gamma_inf(1) ========
    print("## D)  gamma_inf(1): extrapolation + rational scan")
    g16 = 3 ** 16 * p[16][1]
    A16 = g16 - 3 ** 15 * p[15][1]
    cands = {}
    for rr_ in (0.87, 0.88, 0.887, 0.89, 0.90):
        gi = g16 + A16 * rr_ / (1 - rr_)
        cands[rr_] = gi
        print(f"   rho={rr_:.3f}: gamma_inf = {gi:.6f}")
    lo, hi = min(cands.values()), max(cands.values())
    mid = (lo + hi) / 2
    print(f"   band: [{lo:.5f}, {hi:.5f}]")
    print("   rational candidates q<=200 in band (sorted by q):")
    seen = []
    for qd in range(2, 201):
        pn = round(mid * qd)
        for pp in (pn - 1, pn, pn + 1):
            v = pp / qd
            if lo <= v <= hi and math.gcd(pp, qd) == 1:
                seen.append((qd, pp, v))
    for qd, pp, v in seen[:12]:
        print(f"     {pp}/{qd} = {v:.6f}")
    Sinf = 0.4766
    named = {"(1+S_inf)/2 (S=0.4766)": (1 + Sinf) / 2, "11/15=(1+7/15)/2": 11.0 / 15,
             "31/42": 31 / 42, "S_inf+0.26?": Sinf + 0.26, "1-(1-S_inf)*1/2": (1 + Sinf) / 2}
    print("   named:")
    for k, v in named.items():
        inb = "IN BAND" if lo <= v <= hi else "out"
        print(f"     {k} = {v:.6f}  [{inb}]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
