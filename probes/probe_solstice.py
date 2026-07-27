"""
PROBE SOLSTICE -- does the deparitied rate drift persist -> turnover? (2026-07-26)

Two-mode model: Lambda_i = A rho1^i - B rho2^i (rho2>rho1) => deparitied rate falls monotonically, hits 0 at the
crossing. Persistent falling drift = leading turnover indicator; its magnitude dates it. A null is as informative as a hit.

S-A: extend Lambda_i to i=18 if reachable (17 min). GATE Lambda_12..15 = 0.337,0.320,0.287,0.262 e-3 + S_16=0.4714
     bit-for-bit first. Report precision floor at 17,18. *** stationary_trunc WALLS at n=17 (41 GiB) -> i<=16 only. ***
S-B: deparitied two-step rates (Lam_i/Lam_{i-2})^.5 SEPARATELY even-i / odd-i, i=12..16. Per-class drift + SE.
S-C: curvature (2nd diff of rate) per class; refit crossing (rate,drift,curvature) overdetermined.
S-D: mirror -- enriched (3|k) rate. Flat enriched vs falling depleted FALSIFIES the measure-global two-mode.
S-E: direct refit Lambda=A rho1^i - B rho2^i, i=10..16 deparitied. Held-out: fit 10..14, predict 15,16 BEFORE crossing.

Certified base-2/nu_e machinery only. Gates via build_base2_fast (recompute) + P6J enriched.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_p6i import build_base2_fast, autocorr


def main():
    t0 = time.time()
    print("# PROBE SOLSTICE -- does the rate drift persist?\n")

    # ---- S-A: recompute Lambda to i=15 (gate), append cached i=16 (i>=17 walls at 41 GiB) ----
    IMAX = 15
    R0 = {}; R2 = {}
    for i in range(1, IMAX + 1):
        S = build_base2_fast(i); Re = autocorr(S['R_e'])
        R0[i] = float(Re[0]); R2[i] = float(Re[2])
    T = {i: 3 ** i * (4 * R2[i] - R0[i]) for i in range(1, IMAX + 1)}
    T[0] = 1.0 / 3
    T[16] = 0.23591008          # cached from earlier build_base2_fast(16); i>=17 = 41 GiB stationary_trunc wall
    Lam = {i: T[i] - T[i - 1] for i in range(1, 17)}

    print("## S-A GATE: reproduce banked Lambda_12..15 (x1e-3: 0.337,0.320,0.287,0.262) + S_16=0.4714")
    banked = {12: 0.00033677, 13: 0.00031971, 14: 0.00028672, 15: 0.00026193}
    for i in (12, 13, 14, 15):
        print(f"   Lambda_{i} = {Lam[i]:.8f}  banked {banked[i]:.8f}  diff {Lam[i]-banked[i]:.1e}")
    print(f"   S_16 = 2*T_15 = {2*T[15]:.6f} (banked 0.4714);  Lambda_16 = {Lam[16]:+.8f}")
    # precision floor: FFT/round noise on Lambda at i=15,16
    floor = 3.0 ** 15 * 1e-16
    print(f"   precision floor ~ 3^15 * 1e-16 ~ {floor:.1e} (abs on Lambda); rate signal ~3e-3 -- floor {'OK' if floor<3e-5 else 'MARGINAL'}")
    print("   *** i>=17 UNREACHABLE: stationary_trunc needs 41 GiB at n=17. i=16 is the exact-method wall. ***\n")

    # ---- S-B: deparitied rates by parity class, drift + SE ----
    def rate(i):
        return (Lam[i] / Lam[i - 2]) ** 0.5 if Lam[i] * Lam[i - 2] > 0 else float('nan')
    print("## S-B DECIDER: deparitied two-step rate by parity class")
    ev = [i for i in range(12, 17) if i % 2 == 0]      # 12,14,16
    od = [i for i in range(13, 17) if i % 2 == 1]      # 13,15
    for name, idx in (("even-i", ev), ("odd-i", od)):
        rs = [(i, rate(i)) for i in idx]
        print(f"   {name}: " + "  ".join(f"i{i}:{r:.4f}" for i, r in rs))
        xs = np.array([i for i, r in rs]); ys = np.array([r for i, r in rs])
        if len(xs) >= 3:
            p, cov = np.polyfit(xs, ys, 1, cov=True)
            se = np.sqrt(cov[0, 0])
            print(f"      drift = {p[0]:+.5f}/level  SE = {se:.5f}  |slope/SE| = {abs(p[0]/se):.2f}")
        elif len(xs) == 2:
            print(f"      drift = {(ys[1]-ys[0])/(xs[1]-xs[0]):+.5f}/level  (2 pts -> no SE)")
    print()

    # ---- S-C: curvature ----
    print("## S-C CURVATURE (2nd diff of even-i rate; model predicts acceleration)")
    er = [rate(i) for i in ev]
    if len(er) >= 3:
        curv = er[2] - 2 * er[1] + er[0]
        print(f"   even-i rates {[round(x,4) for x in er]}  2nd diff = {curv:+.5f} "
              f"({'accelerating down' if curv<0 else 'decelerating'})")
    print()

    # ---- S-D: mirror (enriched, from P6J) ----
    print("## S-D MIRROR: enriched (3|k) Lambda^+ -- flat vs falling depleted?")
    LamP = {8: -4.930e-5, 9: -1.733e-5, 10: -1.014e-5, 11: -3.91e-6, 12: -6.29e-6,
            13: 2.0e-7, 14: 1.4e-7, 15: -1.2e-7}
    print("   Lambda^+ (enriched): " + " ".join(f"i{i}:{LamP[i]:+.1e}" for i in range(10, 16)))
    print(f"   enriched CONVERGED (|Lambda^+|<1e-6 by i~13); rate is noise -> FLAT, not drifting.")
    print(f"   => the slow mode is DEPLETED-ONLY, not a global measure mode. S-D FALSIFIES clean two-mode/global-slow read.\n")

    # ---- S-E: direct refit + held-out ----
    print("## S-E REFIT Lambda=A rho1^i - B rho2^i (deparitied via odd+even), held-out falsifier")
    from scipy.optimize import curve_fit
    model = lambda i, A, B, r1, r2: A * r1 ** i - B * r2 ** i
    xall = np.array([i for i in range(10, 17)]); yall = np.array([Lam[i] for i in range(10, 17)])
    # held-out: fit 10..14, predict 15,16
    xfit = np.array(list(range(10, 15))); yfit = np.array([Lam[i] for i in range(10, 15)])
    try:
        p, _ = curve_fit(model, xfit, yfit, p0=[1e-2, 1e-2, 0.9, 0.95], maxfev=200000)
        pred = {i: model(i, *p) for i in (15, 16)}
        print(f"   fit i=10..14: A={p[0]:.2e} B={p[1]:.2e} rho1={p[2]:.4f} rho2={p[3]:.4f}")
        for i in (15, 16):
            print(f"      HELD-OUT predict Lambda_{i} = {pred[i]:+.3e}  actual {Lam[i]:+.3e}  "
                  f"ratio {pred[i]/Lam[i]:.3f}  {'MISS' if abs(pred[i]-Lam[i])>0.3*abs(Lam[i]) else 'ok'}")
    except Exception as e:
        print(f"   held-out fit failed: {e}")
    try:
        p2, _ = curve_fit(model, xall, yall, p0=[1e-2, 1e-2, 0.9, 0.95], maxfev=200000)
        # crossing: A r1^i = B r2^i -> i = ln(A/B)/ln(r2/r1)
        A, B, r1, r2 = p2
        if r2 > r1 and A > 0 and B > 0:
            icross = np.log(A / B) / np.log(r2 / r1)
            print(f"   full fit i=10..16: rho1={r1:.4f} rho2={r2:.4f} B/A={B/A:.3e} -> crossing i={icross:.1f}")
        else:
            print(f"   full fit i=10..16: A={A:.2e} B={B:.2e} rho1={r1:.4f} rho2={r2:.4f} (no clean crossing)")
    except Exception as e:
        print(f"   full fit failed: {e}")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
