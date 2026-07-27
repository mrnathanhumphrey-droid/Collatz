"""
PROBE P1 (Wilson) -- amplitude-vs-phase GATE. Result: the bridge is BLIND to the channel's support (2026-07-26).

Refactor BRIDGE2's EXACT transfer, GATE bit-for-bit, then test whether arg(pi-hat) can reach gamma. The gate
DISCRIMINATED: the additive->multiplicative bridge only covers the PRIMITIVE (3-nmid a) part of rho-hat; the channel
gamma lives ENTIRELY in the IMPRIMITIVE (3|a, coarse) part (P4), where the Gauss sum tau_a = 0. So arg(pi-hat) -- which
the bridge maps only to 3-nmid a -- cannot affect gamma. gamma is carried by the frozen coarse marginals (tower-fixed),
= AMPLITUDE, and NOT the fine |pi-hat| work.

GATE (Wilson's bit-for-bit condition): factored bridge == BRIDGE2 direct on 3-nmid a (BRIDGE2's actual validated
domain -- it only tested a in {1,2,4,5,7}). Report the 3|a breakdown + the gamma carrier support to prove blindness.

n=7. Reuses BRIDGE2's chi/DL/tau/offset verbatim + stationary_trunc. Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import gcd
import numpy as np
from probe_27_high_k_rho_q5 import stationary_trunc


def Mclass(k):
    return 5.0 / 3.0 if k % 3 == 0 else 2.0 / 3.0


def build_level(n):
    q = 3 ** (n + 1); Nn = 3 ** n
    DL = {}; g = 1
    for s in range(Nn):
        DL[g] = s; g = (g * 4) % q
    ts = np.array([t for t in range(q) if t % 3 != 0], dtype=np.int64)
    tp = np.where(ts % 3 == 2, (-ts) % q, ts)
    DLt = np.array([DL[int(x)] for x in tp], dtype=np.int64)
    eq = np.exp(2j * np.pi * ts / q)
    piW, _ = stationary_trunc(3, n)
    cpW = np.array([r for r in range(Nn) if gcd(r, 3) == 1], dtype=np.int64)
    nuW = np.asarray(piW, float); nuW = nuW / nuW.sum()
    Y = (3 * cpW + 1) % q
    DLfY = np.array([DL[int(np.where(y % 3 == 2, (-y) % q, y))] for y in Y], dtype=np.int64)
    rho = np.bincount(DLfY, weights=nuW, minlength=Nn)           # certified dlog profile
    dense = np.zeros(Nn, dtype=complex); dense[cpW] = nuW
    What = np.conj(np.fft.fft(dense))
    M = np.exp(-2j * np.pi * np.outer(np.arange(Nn), DLt) / Nn)
    tau = M @ eq
    return dict(q=q, Nn=Nn, ts=ts, eq=eq, What=What, M=M, tau=tau, rho=rho)


def bridge(L, What_vec):
    muY = np.exp(2j * np.pi * L['ts'] / L['q']) * What_vec[L['ts'] % L['Nn']]
    return (L['M'] @ muY) / L['tau']


def main():
    t0 = time.time()
    n = 7
    print(f"# PROBE P1 -- amplitude-vs-phase GATE: the bridge is blind to the channel (n={n})\n")
    L = build_level(n); Nn = L['Nn']
    a = np.arange(Nn)
    rhat0 = bridge(L, L['What'])                                 # factored bridge
    rd = np.conj(np.fft.fft(L['rho']))                           # BRIDGE2 direct (= P4 rho-hat)
    diff = np.abs(rhat0 - rd)

    # --- GATE on BRIDGE2's actual domain (3-nmid a) ---
    prim = (a % 3 != 0)
    imprim = (a % 3 == 0) & (a > 0)
    gate_rel = diff[prim].max() / np.abs(rd[prim]).max()
    print(f"## GATE (BRIDGE2 domain = 3-nmid a): max rel = {gate_rel:.2e}  "
          f"[{'BIT-FOR-BIT OK' if gate_rel < 1e-10 else 'FAIL'}]  ({prim.sum()} primitive chars)")
    print(f"   |tau_a|: 3-nmid a -> {np.abs(L['tau'][prim]).mean():.2f} (= sqrt(q)={np.sqrt(L['q']):.2f}); "
          f"3|a,a>0 -> {np.abs(L['tau'][imprim]).mean():.2e} (IMPRIMITIVE, ~0)")
    print(f"   => bridge FAILS on 3|a (max rel {diff[imprim].max()/np.abs(rd[imprim]).max():.2f}); "
          f"dominant carrier a=N/3=3^{n-1}: bridge |rho-hat|^2={abs(rhat0[3**(n-1)])**2:.4f} "
          f"vs true {abs(rd[3**(n-1)])**2:.4f}, tau={abs(L['tau'][3**(n-1)]):.1e}\n")

    # --- gamma carrier support: 3-nmid a vs 3|a ---
    P = np.abs(rd) ** 2
    def gcontrib(k, mask):
        return float(np.sum(P[mask] * np.cos(2 * np.pi * a[mask] * k / Nn)))
    print("## gamma carrier support: contribution to gamma_n(k)-1 from 3-nmid a (bridge domain) vs 3|a (channel)")
    print(f"   {'k':>2} {'gamma-1':>9} {'from 3-nmid a':>14} {'from 3|a':>10} {'3|a share':>10}")
    for k in (1, 2, 3, 6):
        g1 = gcontrib(k, a > 0)
        gp = gcontrib(k, prim); gi = gcontrib(k, imprim)
        print(f"   {k:>2} {g1:>+9.4f} {gp:>+14.5f} {gi:>+10.4f} {100*gi/g1 if abs(g1)>1e-9 else 0:>9.2f}%")
    print("   [3|a carries ~all of gamma; the bridge (arg pi-hat operator) only reaches 3-nmid a => cannot affect gamma.]")

    # --- the cascade, split by class, all carried by 3|a ---
    print("\n## cascade c_k = gamma_n(k) - M_class(k)  (the channel-distinguishing part; all in 3|a)")
    print(f"   {'k':>2} {'3|k':>4} {'M':>6} {'gamma':>8} {'c_k':>9}")
    for k in range(1, 8):
        g = 1.0 + gcontrib(k, a > 0)
        print(f"   {k:>2} {'yes' if k%3==0 else 'no':>4} {Mclass(k):>6.3f} {g:>8.4f} {g-Mclass(k):>+9.4f}")

    print("\n## VERDICT: arg(pi-hat) cannot carry the channel.")
    print("   gamma = Sum|rho-hat|^2 e(ak/N) (Wiener-Khinchin, amplitude-only in rho-hat). Its support is 3|a (coarse,")
    print("   = frozen mod-3^j marginals, tower-fixed). The bridge maps arg(pi-hat) ONLY to 3-nmid a, where tau_a=sqrt(q);")
    print("   on 3|a it has tau_a=0 and is BLIND. So the channel is AMPLITUDE-carried by the frozen coarse marginals,")
    print("   and NOT by the fine |pi-hat| work (sup/tail, all 3-nmid a) -- the seam is a THIRD object. Scramble is moot:")
    print("   the operator that would carry phase is blind to gamma's support.")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
