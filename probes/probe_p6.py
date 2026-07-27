"""
PROBE P6 (Wilson) -- the 2-adic ORIGIN of the phase (per level, 2026-07-26).

The cascade is arg(pi-hat); arg(pi-hat) = arg integral over Z_2; which part of the DOMAIN makes it?
Conditioning on first L branch vars (a_1..a_L) = L-fold unrolled SINGLEREC (certified). Each cylinder contributes:
weight 2^-(a_1+..+a_L), an explicit phase, and pi-hat_j at 3^L xi 2^-(sum) mod 3^j. DO NOT build new decomp -- unroll.

A GATE: depth-L cylinder sum == pi-hat_j(xi), L=1,2,3, sample xi, levels j<=8. Machine precision or STOP.
B COHERENCE: per xi, coherence = |Sum contrib|/Sum|contrib| (1=aligned,0=cancelling); which a carries max modulus.
C INFLUENCE: recompute shell A_j(1) with cylinders (a_1) restricted to subsets -> which 2-adic regions make the +shell.
D SIGN: sign(A_j(1)) per subset. If one family carries +, complement -, that family = mechanism.
PRE-REG: branch vars iid => no distinguished region a priori (homogeneity). ONLY break = two-primes: 2^-a mod 3^j
depends on a mod ord_{3^j}(2)=a mod 2*3^{j-1}. Structure should appear on a mod 2*3^{j-1}, NOT a. Report both.
FAILURE MODE (real answer): uniformly low diffuse coherence across every cylinder/depth = phase has NO local source.

Reuses probe_p1.build_level/bridge (factored+gated) + fwd_hat. Per level. No new transport.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_p1 import build_level, bridge


def partial_pihat(What, N, cond, Amax=64):
    """L=1 cylinder reconstruction of pi-hat restricted to a satisfying cond(a)."""
    xi = np.arange(N); inv2 = pow(2, -1, N)
    out = np.zeros(N, dtype=complex); p = 1
    for a in range(1, Amax):
        p = (p * inv2) % N
        if cond(a):
            w = (xi * p) % N
            out += (0.5 ** a) * np.exp(2j * np.pi * w / N) * What[(3 * w) % N]
    return out


def cyl_sum(What, N, xi, L, Amax=52):
    """depth-L unrolled SINGLEREC at a single xi (recursive)."""
    inv2 = pow(2, -1, N); tot = [0j]
    def rec(cx, depth, wt, ph):
        if depth == L:
            tot[0] += wt * np.exp(2j * np.pi * (ph % N) / N) * What[cx % N]; return
        p = 1
        for a in range(1, Amax):
            p = (p * inv2) % N
            w = (cx * p) % N
            rec((3 * w) % N, depth + 1, wt * 0.5 ** a, ph + w)
    rec(xi % N, 0, 1.0, 0)
    return tot[0]


def shellA1(L, What_vec):
    rhat = bridge(L, What_vec); Nn = L['Nn']; a = np.arange(Nn); prim = a % 3 != 0
    P = np.abs(rhat[prim]) ** 2; ap = a[prim]
    return float(np.sum(P * np.cos(2 * np.pi * ap / Nn)))


def main():
    t0 = time.time()
    print("# PROBE P6 -- 2-adic origin of the phase (per level)\n")

    # ---------- A GATE ----------
    print("## P6-A GATE: depth-L cylinder sum == pi-hat_j(xi)  (L=1,2,3; 15 sample primitive xi)")
    rng = np.random.default_rng(0)
    for j in (4, 6, 8):
        L = build_level(j); N = 3 ** j; W = L['What']
        prim = np.array([x for x in range(1, N) if x % 3 != 0])
        samp = rng.choice(prim, size=min(15, len(prim)), replace=False)
        row = []
        for Ldep in (1, 2, 3):
            err = max(abs(cyl_sum(W, N, int(x), Ldep) - W[x]) for x in samp)
            row.append(f"L={Ldep}:{err:.1e}")
        print(f"   j={j}: " + "  ".join(row))
    print("   [must be ~machine; else unrolling mis-indexed.]\n")

    # ---------- B COHERENCE ----------
    print("## P6-B COHERENCE (L=1): |Sum_a term_a|/Sum_a|term_a| over primitive xi; argmax-a")
    for j in (4, 6, 8):
        L = build_level(j); N = 3 ** j; W = L['What']
        prim = np.array([x for x in range(1, N) if x % 3 != 0])
        inv2 = pow(2, -1, N)
        num = np.abs(W[prim])
        den = np.zeros(len(prim)); p = 1
        maxmod = np.zeros(len(prim)); argmax_a = np.zeros(len(prim), dtype=int)
        for a in range(1, 60):
            p = (p * inv2) % N
            w = (prim * p) % N
            m = (0.5 ** a) * np.abs(W[(3 * w) % N])
            den += m
            upd = m > maxmod; maxmod = np.where(upd, m, maxmod); argmax_a = np.where(upd, a, argmax_a)
        coh = num / (den + 1e-30)
        print(f"   j={j}: coherence mean={coh.mean():.4f} median={np.median(coh):.4f} "
              f"min={coh.min():.4f} max={coh.max():.4f} | argmax-a: median={int(np.median(argmax_a))} "
              f"(a=1 share {np.mean(argmax_a==1)*100:.0f}%)")
    print("   [uniformly low+diffuse => no local source (real answer). high/variable => structured.]\n")

    # ---------- C/D INFLUENCE + SIGN ----------
    print("## P6-C/D INFLUENCE: shell A_j(1) with a_1 restricted to subsets  (full A_1 gated in P1LVL)")
    subsets = [
        ("all", lambda a: True),
        ("a=1", lambda a: a == 1),
        ("a<=3", lambda a: a <= 3),
        ("a>=4", lambda a: a >= 4),
        ("a even", lambda a: a % 2 == 0),
        ("a odd", lambda a: a % 2 == 1),
        ("a=1 mod3", lambda a: a % 3 == 1),
        ("a=2 mod3", lambda a: a % 3 == 2),
        ("a=0 mod3", lambda a: a % 3 == 0),
    ]
    for j in (4, 6, 8):
        L = build_level(j); N = 3 ** j; W = L['What']
        print(f"   -- j={j} (ord_3^j(2)=2*3^{j-1}={2*3**(j-1)}) --")
        for name, cond in subsets:
            pw = partial_pihat(W, N, cond)
            A = shellA1(L, pw)
            print(f"      {name:>10}: A_j(1) = {A:>+10.5f}  sign {'+' if A>0 else '-'}")

    # ---------- two-primes test: A_j(1) by a mod ord (where ord < Amax, i.e. j<=4) ----------
    print("\n## P6 two-primes test: A_j(1) restricted to a == r mod ord  (j=4, ord=54; the a-mod-ord structure)")
    j = 4; L = build_level(j); N = 3 ** j; W = L['What']; ord_ = 2 * 3 ** (j - 1)
    vals = []
    for r in range(ord_):
        A = shellA1(L, partial_pihat(W, N, lambda a, rr=r, o=ord_: a % o == rr))
        vals.append((r, A))
    vals.sort(key=lambda t: -abs(t[1]))
    print("   top |A| by residue r mod 54: " + ", ".join(f"r={r}:{A:+.4f}" for r, A in vals[:8]))
    allA = np.array([A for _, A in vals])
    print(f"   spread over 54 residues: mean {allA.mean():+.5f}, std {allA.std():.5f}, "
          f"pos {np.mean(allA>0)*100:.0f}% / neg {np.mean(allA<0)*100:.0f}%")
    print("   [pre-reg: structure on a-mod-ord (two primes meet) vs null on a. flag which.]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
