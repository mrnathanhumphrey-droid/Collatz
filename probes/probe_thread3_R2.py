"""
PROBE R2 -- THE AMPLITUDE (the 7/10 gate; thread 3 session two). Cheap, exact.
A_k := 3^k * P(pair agrees mod 3^k) = 3^k * sum_{m>=k} g_m,  g_m = a_m - a_{m-1}/3,  a_m = ||pi_m||^2.
Recursion S_k = A_k - A_{k+1}/3 = 3^k g_k = c_k (R1's welded shell sequence). Since c_k -> 7/15,
A_inf = (3/2)(7/15) = 7/10 (algebraic).
R2-A: A_k -> 7/10 from the frozen chain; consistency S_k reproduces welded 2/3, 10/21.
R2-B: CLOSED FORM of the c0-mode overlap from Real-T1's ell_0 = R_0(e_rho)/R_0(0) on the zero-carry sector,
      against the independent-pair state v_indep = w (x) w.  <ell_0|v_indep> = Sum_e R_0(e)^2 / R_0(0).
      PRE-REG: = 7/10 exactly. Report the exact rational either way.
R2-C: q=5,7 A_k -> 0 geometrically (the gap as vanishing agreement amplitude).
Guards: no fitting; exact rationals where feasible; deviations as deviations.
"""
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from probe_phase2a_q2b_q6 import build_M_gen, subgroup

SEVEN_10 = Fr(7, 10); SEVEN_15 = Fr(7, 15)


def full_M(q, L, lam=0.5):
    qL = q ** L; sub = subgroup(2, qL); ordn = len(sub)
    raw = [lam ** d for d in range(1, ordn + 1)]
    M, idx, n = build_M_gen(q, L, 2, raw)
    return M.tocsr(), idx


def a_sequence(M, idx, K):
    """a_k = 1^T M^k v0 (float, frozen). a_0 = 1."""
    n = M.shape[0]; v = np.zeros(n); v[idx[(1, 1, 0)]] = 1.0
    a = [1.0]
    for _ in range(K):
        v = M.dot(v); a.append(v.sum())
    return a                                                    # a[m] = a_m, m=0..K


def exact_a(q, L, K, lam=Fr(1, 2)):
    """Exact a_m = ||pi_m||^2 for m<=L (frozen chain exact regime)."""
    qL = q ** L; sub = subgroup(2, qL); ordn = len(sub); inv = pow(2, -1, qL)
    raw = [lam ** d for d in range(1, ordn + 1)]; Z = sum(raw); w = [r / Z for r in raw]
    mult = [(pow(inv, d, qL), w[d - 1]) for d in range(1, ordn + 1)]
    v = defaultdict(Fr); v[(1, 1, 0)] = Fr(1); a = [Fr(1)]
    for _ in range(K):
        nv = defaultdict(Fr)
        for (a_, b_, g), val in v.items():
            for (ga, wa) in mult:
                ap = (a_ * ga) % qL
                for (gb, wb) in mult:
                    bp = (b_ * gb) % qL
                    T = (ap - bp) % qL
                    if (g + T) % q == 0:
                        gp = ((g + T) // q) % qL
                        nv[(ap, bp, gp)] += val * wa * wb
        v = nv; a.append(sum(v.values()))
    return a                                                    # a[m], m=0..K exact


# ============================ R2-A ============================
def R2A():
    print(f"\n{'='*88}\n## R2-A  A_k = 3^k * P(agree mod 3^k) -> 7/10  (frozen chain; consistency with welded 2/3,10/21)")
    M, idx = full_M(3, 3)
    a = a_sequence(M, idx, K=45)
    g = [a[m] - a[m-1] / 3 for m in range(1, len(a))]           # g[m-1] = g_m, m=1..
    # A_k = 3^k * sum_{m>=k} g_m  (tail: g_m ~ (7/15) 3^{-m})
    A = {}
    for k in range(1, 11):
        tail = sum(g[m-1] for m in range(k, len(a)))            # partial sum up to m=44
        A[k] = (3.0 ** k) * tail
    print(f"     k :      A_k            S_k=A_k-A_(k+1)/3   [S_k should track c_k=welded shell]")
    for k in range(1, 10):
        Sk = A[k] - A[k+1] / 3
        print(f"     {k:2d}: {A[k]:+.8f}        {Sk:+.8f}")
    print(f"   >> A_10 = {A[10]:.8f}   (7/10 = 0.70000000; frozen L=3 -> ~0.66 by the same -2.5% shell truncation as R1)")
    # exact consistency anchors: S_1, S_2 from exact a
    ae = exact_a(3, 3, 3)
    ge = [ae[m] - ae[m-1]/3 for m in range(1, 4)]               # g_1,g_2,g_3 exact
    S1 = 3 * ge[0]; S2 = 9 * ge[1]
    print(f"   >> EXACT consistency: S_1 = 3 g_1 = {S1} ({'==2/3' if S1==Fr(2,3) else 'DEV'}),  "
          f"S_2 = 9 g_2 = {S2} ({'==10/21' if S2==Fr(10,21) else 'DEV'})  "
          f"[S_k = 3^k g_k = c_k, welded]")
    print(f"   >> ALGEBRAIC: A_inf = (3/2) S_inf = (3/2)(7/15) = {Fr(3,2)*SEVEN_15} "
          f"({'== 7/10 EXACT' if Fr(3,2)*SEVEN_15==SEVEN_10 else 'DEV'})")


# ============================ R2-B ============================
def Rk_autocorr(wf, D, k=0):
    """R_k(e) = sum_delta w_delta w_{delta+e} omega^{delta}, k=0 -> real autocorrelation."""
    if k == 0:
        return [sum(wf[j] * wf[(j + e) % D] for j in range(D)) for e in range(D)]
    raise NotImplementedError


def R2B():
    print(f"\n{'='*88}\n## R2-B  CLOSED FORM: c0-mode overlap <ell_0|v_indep> = Sum_e R_0(e)^2 / R_0(0)  (PRE-REG 7/10)")
    for L in (2, 3):
        qL = 3 ** L; sub = subgroup(2, qL); D = len(sub)
        raw = [Fr(1, 2) ** d for d in range(1, D + 1)]; Z = sum(raw); wf = [r / Z for r in raw]
        R0 = Rk_autocorr(wf, D, 0)
        R0_0 = R0[0]                                            # = Sum w^2 = c_0 (finite-D)
        num = sum(r * r for r in R0)                            # Sum_e R_0(e)^2
        overlap = num / R0_0
        print(f"   L={L} (D={D}): R_0(0)=Sum w^2 = {R0_0} ({float(R0_0):.6f});  "
              f"Sum_e R_0(e)^2 = {num} ({float(num):.6f})")
        print(f"        <ell_0|v_indep> = Sum_e R_0(e)^2 / R_0(0) = {overlap} = {float(overlap):.8f}   "
              f"{'== 7/10 EXACT (CROWN)' if overlap==SEVEN_10 else f'!= 7/10 (dev {float(overlap-SEVEN_10):+.4f}); reported as the exact rational it IS'}")


# ============================ R2-C ============================
def R2C():
    print(f"\n{'='*88}\n## R2-C  q=5,7 CONTRAST: A_k -> 0 geometrically (gap = vanishing agreement amplitude)")
    for q, L in [(5, 2), (7, 2)]:
        M, idx = full_M(q, L)
        a = a_sequence(M, idx, K=40)
        g = [a[m] - a[m-1] / 3 for m in range(1, len(a))]
        A = [(3.0 ** k) * sum(g[m-1] for m in range(k, len(a))) for k in range(1, 11)]
        print(f"   q={q} (L={L}): A_k, k=1..10:")
        print("     " + "  ".join(f"{v:.5f}" for v in A))
        rats = [A[i+1]/A[i] if A[i] != 0 else float('nan') for i in range(len(A)-1)]
        print(f"     A_(k+1)/A_k: " + "  ".join(f"{r:+.4f}" for r in rats) +
              f"   (pre-reg -> 0 geometrically; 3/q={3/q:.4f})")


def main():
    print("# PROBE R2 -- THE AMPLITUDE (7/10 gate). Cheap, exact. No fitting.")
    R2A()
    R2B()
    R2C()


if __name__ == "__main__":
    main()
