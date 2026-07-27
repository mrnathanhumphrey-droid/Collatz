"""
PROBE P6B (Wilson) -- resolve the parity mechanism + test the cross-term collapse (2026-07-26).

(1) MECHANISM (odd/even-lag corollary): decompose A_j(m) by parity for m=1..6. Which m are cross-parity (same=0)
    vs same-parity (cross=0)? Settles the load-bearing sentence from DATA (Claude asserted odd-lag=cross; Wilson's
    dlog arith gives same-parity). No assertion -- measure.
(2) COLLAPSE (Wilson pen): nu_o = 1/2 (x2^-1)_* nu_e exactly (odd branch = even, halved+translated). Test at the
    dlog-profile level: is rho_o[s] = c * rho_e[(s-delta)%N] for some (c,delta)? Report c, delta, residual. If exact,
    the m=1 cross-term collapses to a SINGLE sub-measure (rho_e) autocorrelation at ONE lag (ratio -2 = gen of <2>/<4>).
(3) CONVENTION: stationary_trunc uses a>=1 (v>=1). Wilson's caveat: a>=1 leaves a boundary term 1/2 (m_1)_* nu vs a>=0.

Reuses probe_p1.build_level/bridge + probe_p6.partial_pihat. Per level. No new transport.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_p1 import build_level, bridge
from probe_p6 import partial_pihat


def shellA(L, What_vec, m):
    rhat = bridge(L, What_vec); Nn = L['Nn']; a = np.arange(Nn); prim = a % 3 != 0
    P = np.abs(rhat[prim]) ** 2; ap = a[prim]
    return float(np.sum(P * np.cos(2 * np.pi * ap * m / Nn)))


def rho_from(L, What_vec):
    """dlog profile whose conj(fft) = bridge(What_vec)."""
    return np.fft.ifft(np.conj(bridge(L, What_vec))).real


def main():
    t0 = time.time()
    print("# PROBE P6B -- parity mechanism + cross-term collapse\n")

    # ---------- (1) MECHANISM: parity decomposition of A_j(m) ----------
    print("## (1) MECHANISM: A_j(m) = A[even] + A[odd] + cross, for lags m=1..6")
    print("   [cross==all & same==0 => CROSS-parity lag ; cross==0 => SAME-parity lag]")
    for j in (4, 6):
        L = build_level(j); N = 3 ** j; W = L['What']
        pe = partial_pihat(W, N, lambda a: a % 2 == 0)
        po = partial_pihat(W, N, lambda a: a % 2 == 1)
        print(f"   -- j={j} --")
        print(f"   {'m':>2} {'A_all':>10} {'A_even':>10} {'A_odd':>10} {'cross':>10} {'verdict':>10}")
        for m in range(1, 7):
            Aa = shellA(L, W, m); Ae = shellA(L, pe, m); Ao = shellA(L, po, m)
            cr = Aa - Ae - Ao
            same = abs(Ae) + abs(Ao)
            verd = 'CROSS' if abs(cr) > 10 * (same + 1e-15) else ('SAME' if same > 10 * abs(cr) else 'mixed')
            print(f"   {m:>2} {Aa:>+10.3e} {Ae:>+10.3e} {Ao:>+10.3e} {cr:>+10.3e} {verd:>10}")
    print()

    # ---------- (2) COLLAPSE: rho_o = c * shift_delta rho_e ? ----------
    print("## (2) COLLAPSE: is rho_o[s] = c * rho_e[(s-delta)%N] ?  (Wilson: c=1/2, delta=dlog(2^-1))")
    for j in (4, 6, 8):
        L = build_level(j); N = 3 ** j; W = L['What']
        pe = partial_pihat(W, N, lambda a: a % 2 == 0)
        po = partial_pihat(W, N, lambda a: a % 2 == 1)
        rho_e = rho_from(L, pe); rho_o = rho_from(L, po); rho_full = rho_from(L, W)
        # gate: rho_e+rho_o == rho_full
        gate = np.max(np.abs(rho_e + rho_o - rho_full))
        # find best shift delta: maximize <rho_o, roll(rho_e,delta)>; then c = that / ||rho_e||^2
        ne = float(np.dot(rho_e, rho_e))
        corr = np.array([float(np.dot(rho_o, np.roll(rho_e, d))) for d in range(N)]) if N <= 20000 else None
        if corr is not None:
            delta = int(np.argmax(np.abs(corr)))
            c = corr[delta] / ne
            resid = np.max(np.abs(rho_o - c * np.roll(rho_e, delta))) / (np.max(np.abs(rho_o)) + 1e-30)
            # is delta = 2^-1 mod N ?
            inv2 = pow(2, -1, N)
            print(f"   j={j}: rho_e+rho_o==rho_full {gate:.1e} | best shift delta={delta} "
                  f"(2^-1 mod N={inv2}, match {'YES' if delta==inv2 else 'no'}), c={c:.5f} (1/2?), "
                  f"residual {resid:.2e} [{'EXACT' if resid<1e-9 else 'approx'}]")
        else:
            print(f"   j={j}: N={N} too large for full shift scan; gate rho_e+rho_o {gate:.1e}")
    print()

    # ---------- (3) the collapsed form: A_j(1) vs 1/2 autocorr(rho_e) at the shifted lag ----------
    print("## (3) A_j(1) vs the collapsed single-measure form  1/2 <rho_e, shift_{1-delta} rho_e> * 3^j")
    for j in (4, 6, 8):
        L = build_level(j); N = 3 ** j; W = L['What']
        pe = partial_pihat(W, N, lambda a: a % 2 == 0)
        rho_e = rho_from(L, pe)
        A1 = shellA(L, W, 1)
        inv2 = pow(2, -1, N)
        # cross-parity lag-1 collapses to lag (1 - dlog(2^-1)); Wilson: 1-h == h == 2^-1
        lag = inv2
        collapsed = 3.0 ** j * float(np.dot(rho_e, np.roll(rho_e, -lag)))   # <rho_e, shift_lag rho_e>*3^j
        print(f"   j={j}: A_j(1)={A1:+.5f}  1/2*collapsed(lag=2^-1)={0.5*collapsed:+.5f}  "
              f"ratio {0.5*collapsed/A1 if abs(A1)>1e-9 else float('nan'):.3f}")
    print("   [convention: stationary_trunc uses a>=1 (v>=1); Wilson's exact collapse assumes a>=0, leftover 1/2(m_1)nu.]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
