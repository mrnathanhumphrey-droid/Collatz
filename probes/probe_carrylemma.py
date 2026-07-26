"""
PROBE CARRYLEMMA -- gate Wilson's carry identity + attenuation mechanism (2026-07-25).

Setup (level pair M=3^r): x = class mod M (integer rep in [0,M)), c(x) = floor(4x/M) in {0,1,2,3} (DF Prop 5.1
carry), L(x) = 4x mod M. New digits D, D' at weight M; lift event D' = c + D (mod 3) [since 4 = 1 mod 3].
  pi_x(d) = nu_hi(x + d M)/nu_low(x),  dpi = pi - 1/3,  (T_c f)(d) = f(d+c).
LEMMA:  P(lift|x) - 1/3 = <dpi_x, T_{c(x)} dpi_{L(x)}>,   q_r(1) - 1/3 = E_mu <dpi_X, T_c dpi_{4X}>,
        mu(x) prop nu_low(x) nu_low(L(x))   [survival conditioning FULLY ABSORBED into the mu-weighting].
MECHANISM: u=<dpi_x,dpi_L>, v=T1-rot, w=T2-rot; u+v+w=0 identically (dpi mean-zero on Z/3). Carry cells
{0,1,2,3} -> rotations {0,1,2,0}: rotation 0 gets TWO cells => under equal cells + independence,
E = u/2 - u/4 = u/4: ATTENUATION BY EXACTLY 4, SIGN PRESERVED. Exact form (no independence):
  q-1/3 = E[(1_{c in {0,3}} - 1_{c=2}) u] + E[(1_{c=1} - 1_{c=2}) v].

GATES:
 G0 lemma identity: q_r(1)-1/3 == E_mu<dpi,T_c dpi>  (exact Fractions r=4,5; float r=4..16).
 G1 u+v+w = 0 per x (machine zero).
 G2 two-term decomposition == q-1/3 exactly.
 P1 (Wilson pred 1): E_mu[u] vs 4(q-1/3)  (ratio ~1 <=> c independent of u; at r=16 predicted E[u]~1.7e-3).
 P2 (Wilson pred 2): cell measures P0..P3 under mu (p1~p2 => v-term drops). Plus E[u|c=j] and the (c,u) coupling.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
from probe_gapop_R28 import build_nu
from probe_ratio2 import build_nu_exact

RTOP = 16


def dense(nu, mod):
    a = np.zeros(mod)
    for X, w in nu.items():
        a[X % mod] += float(w)
    return a / a.sum()


def gates_float(nu_hi_d, r):
    M = 3 ** r
    hi = nu_hi_d                                    # dense over 3^{r+1}, normalized
    low = hi[:M] + hi[M:2 * M] + hi[2 * M:3 * M]
    xs = np.nonzero(low > 0)[0]
    c = (4 * xs) // M
    L = (4 * xs) % M
    wgt = low[xs] * low[L]
    m = wgt > 0
    xs, c, L, wgt = xs[m], c[m], L[m], wgt[m]
    Z = wgt.sum()
    # dpi vectors
    dpx = np.stack([hi[xs + d * M] / low[xs] for d in range(3)], 1) - 1.0 / 3
    dpL = np.stack([hi[L + d * M] / low[L] for d in range(3)], 1) - 1.0 / 3
    u = (dpx * dpL).sum(1)
    v = (dpx * np.roll(dpL, -1, axis=1)).sum(1)
    wr = (dpx * np.roll(dpL, -2, axis=1)).sum(1)
    zero = np.max(np.abs(u + v + wr))
    # coherence with actual rotation c mod 3
    cm = c % 3
    t = np.where(cm == 0, u, np.where(cm == 1, v, wr))
    lem = float((wgt * t).sum() / Z)
    # direct q
    Mp = 3 * M
    idx4 = (4 * np.arange(Mp)) % Mp
    p_hi = float(np.dot(hi, hi[idx4]))
    p_low = float(Z)
    q = p_hi / p_low
    # decomposition + predictions
    Eu = float((wgt * u).sum() / Z)
    Ev = float((wgt * v).sum() / Z)
    T1 = float((wgt * (((c == 0) | (c == 3)).astype(float) - (c == 2)) * u).sum() / Z)
    T2 = float((wgt * (((c == 1)).astype(float) - (c == 2)) * v).sum() / Z)
    P = [float(wgt[c == j].sum() / Z) for j in range(4)]
    Euc = [float((wgt[c == j] * u[c == j]).sum() / wgt[c == j].sum()) if (c == j).any() else float('nan')
           for j in range(4)]
    return dict(q=q, lem=lem, zero=zero, Eu=Eu, Ev=Ev, T1=T1, T2=T2, P=P, Euc=Euc)


def gates_exact(nex, r):
    M = 3 ** r
    hi = {}
    for X, w in nex[r].items():
        hi[X % (3 * M)] = hi.get(X % (3 * M), F(0)) + w
    tot = sum(hi.values())
    hi = {k: w / tot for k, w in hi.items()}
    low = {}
    for X, w in hi.items():
        low[X % M] = low.get(X % M, F(0)) + w
    Z = F(0); acc = F(0); accu = F(0); T1 = F(0); T2 = F(0); Pc = [F(0)] * 4
    for x, lw in low.items():
        L = (4 * x) % M; c = (4 * x) // M
        lwL = low.get(L, F(0))
        w0 = lw * lwL
        if w0 == 0:
            continue
        dpx = [hi.get(x + d * M, F(0)) / lw - F(1, 3) for d in range(3)]
        dpL = [hi.get(L + d * M, F(0)) / lwL - F(1, 3) for d in range(3)]
        u = sum(dpx[d] * dpL[d] for d in range(3))
        v = sum(dpx[d] * dpL[(d + 1) % 3] for d in range(3))
        wr = sum(dpx[d] * dpL[(d + 2) % 3] for d in range(3))
        assert u + v + wr == 0
        t = (u, v, wr)[c % 3]
        Z += w0; acc += w0 * t; accu += w0 * u; Pc[c] += w0
        T1 += w0 * ((1 if c in (0, 3) else 0) - (1 if c == 2 else 0)) * u
        T2 += w0 * ((1 if c == 1 else 0) - (1 if c == 2 else 0)) * v
    # exact q
    p_hi = F(0)
    for X, w in hi.items():
        p_hi += w * hi.get((4 * X) % (3 * M), F(0))
    q = p_hi / Z
    return q - F(1, 3), acc / Z, accu / Z, (T1 + T2) / Z


def main():
    t0 = time.time()
    print("# PROBE CARRYLEMMA -- gate the carry identity + x4 attenuation\n")

    print("## EXACT (Fractions) r=4,5: lemma identity q-1/3 == E_mu<dpi,T_c dpi> ; u+v+w=0 asserted per x")
    nex = build_nu_exact(5)
    for r in (4, 5):
        d_q, d_lem, d_u, d_T = gates_exact(nex, r)
        print(f"   r={r}: q-1/3 = {float(d_q):+.9e} | lemma = {float(d_lem):+.9e} | EQUAL as Fractions: {d_q == d_lem}"
              f" | T1+T2 == q-1/3: {d_T == d_q} | E[u] = {float(d_u):+.6e}")
    print()

    print(f"## FLOAT r=4..{RTOP}  (build_nu to {RTOP} ... ~6-10 min)")
    nus = build_nu(0.5, RTOP)
    print(f"   built ({time.time()-t0:.1f}s)\n")
    print(f"   {'r':>2} {'q-1/3':>11} {'lemma':>11} {'rel':>8} {'E[u]':>11} {'E[u]/4(q-1/3)':>13} "
          f"{'T1':>11} {'T2':>11} {'uvw0':>8}")
    rows = {}
    for r in range(4, RTOP + 1):
        hi = dense(nus[r], 3 ** (r + 1))
        g = gates_float(hi, r)
        del hi
        ex = g['q'] - 1.0 / 3
        rel = abs(g['lem'] - ex) / abs(ex)
        ratio = g['Eu'] / (4 * ex)
        rows[r] = g
        print(f"   {r:>2} {ex:>+11.4e} {g['lem']:>+11.4e} {rel:>8.1e} {g['Eu']:>+11.4e} {ratio:>13.4f} "
              f"{g['T1']:>+11.4e} {g['T2']:>+11.4e} {g['zero']:>8.1e}")
        if r in (8, 12, 16):
            P = g['P']; Euc = g['Euc']
            print(f"        cells P0..P3 = {P[0]:.4f} {P[1]:.4f} {P[2]:.4f} {P[3]:.4f}   "
                  f"E[u|c] = {Euc[0]:+.3e} {Euc[1]:+.3e} {Euc[2]:+.3e} {Euc[3]:+.3e}")
    print()
    g16 = rows[RTOP]
    ex16 = g16['q'] - 1.0 / 3
    print("## READS")
    print(f"   P1: E[u] at r={RTOP} = {g16['Eu']:+.4e} vs predicted ~4x(q-1/3) = {4*ex16:+.4e} "
          f"(ratio {g16['Eu']/(4*ex16):.3f}; ~1 => c independent of u; else (c,u) coupling live = delta-2 answer)")
    print(f"   P2: cells P1 vs P2: {g16['P'][1]:.5f} vs {g16['P'][2]:.5f} "
          f"({'~equal => v-term drops' if abs(g16['P'][1]-g16['P'][2]) < 0.01 else 'SKEWED => v-term live'}); "
          f"T2/(q-1/3) = {g16['T2']/ex16:+.3f}")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
