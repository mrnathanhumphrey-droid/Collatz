"""
GATE for Wilson's THEOREM (Real-T1): the closed-form left eigenvectors of the real
q=3 pair operator are twisted autocorrelations of the halving weights.

M = build_M_gen(3, L, 2, [lam^d]), lam=1/2. D = 2*3^(L-1). Folded weights w_delta ~ 2^-delta.
omega = exp(2i pi k / D).  Twisted autocorrelation  R_k(e) = sum_delta w_delta w_{delta+e} omega^delta.
Claim: for each k in Z/D,
   ell_k(a,b,gamma) = omega^{-e_a} * R_k(e_rho)/R_k(0) * [gamma==0]
is an exact LEFT eigenvector, eigenvalue c_k = R_k(0) = sum_delta w_delta^2 omega^delta.
   e_a   = dlog_2(a)              (a = 2^{e_a} mod 3^L)
   e_rho = dlog_2(b * a^{-1})     (the class coordinate e = log2(b a^-1))

GATE: residual ||M^T ell_k - c_k ell_k||_inf / ||ell_k||_inf.  Pre-registered EXACT <= 1e-12,
for ALL D members, at L=3 (18 members).  Direct methods only (sparse mat-vec; no ARPACK).
Lock convention on L=2 (reproduce Wilson's 6/6) first.
"""
import numpy as np

from probe_phase2a_q2b_q6 import build_M_gen, subgroup

Q = 3


def setup(L, lam=0.5):
    qL = Q ** L
    sub = subgroup(2 % qL, qL)
    D = len(sub)
    dlog = {}
    x = 1 % qL
    for e in range(D):
        dlog[x] = e
        x = (x * 2) % qL
    raw = [lam ** d for d in range(1, D + 1)]
    wf = np.array(raw) / sum(raw)              # wf[j] = weight of phase 2^-(j+1), j=0..D-1
    M, idx, n = build_M_gen(Q, L, 2, raw)
    states = [None] * n
    for s, i in idx.items():
        states[i] = s
    return qL, sub, D, wf, M, idx, n, states, dlog


def Rk_array(wf, D, k):
    om = np.exp(2j * np.pi * k / D)
    powers = om ** (np.arange(D) + 1)          # omega^delta, delta=1..D  (matches c_k in probe_H)
    R = np.zeros(D, complex)
    for e in range(D):
        R[e] = np.sum(wf * np.roll(wf, -e) * powers)   # sum_j wf[j] wf[(j+e)%D] omega^{j+1}
    return R                                    # R[0] = c_k


def build_ellk(k, L, qL, D, wf, n, states, dlog):
    om = np.exp(2j * np.pi * k / D)
    R = Rk_array(wf, D, k)
    R0 = R[0]
    ell = np.zeros(n, complex)
    for i, (a, b, g) in enumerate(states):
        if g != 0:
            continue
        e_a = dlog[a % qL]
        e_rho = dlog[(b * pow(a, -1, qL)) % qL]
        ell[i] = om ** (-e_a) * R[e_rho] / R0
    return ell, R0


def gate(L):
    qL, sub, D, wf, M, idx, n, states, dlog = setup(L)
    MT = M.T.tocsr()
    print(f"## L={L}  D={D}  dim={n}   (18-member gate)" if D == 18 else f"## L={L}  D={D}  dim={n}")
    worst = 0.0
    for k in range(D):
        ell, ck = build_ellk(k, L, qL, D, wf, n, states, dlog)
        resid = MT.dot(ell) - ck * ell
        rel = np.max(np.abs(resid)) / (np.max(np.abs(ell)) + 1e-300)
        worst = max(worst, rel)
        tag = "PASS" if rel <= 1e-12 else "FAIL"
        print(f"   k={k:<2} c_k=R_k(0)={ck.real:+.6f}{ck.imag:+.6f}j  ||M^T l - c l||/||l|| = {rel:.2e}  [{tag}]")
    print(f"   => worst residual over all {D} members: {worst:.2e}  "
          f"[{'ALL PASS (<=1e-12)' if worst <= 1e-12 else 'FAIL'}]")
    return worst


def main():
    print("# GATE -- Wilson's Real-T1: ell_k = omega^{-e_a} R_k(e_rho)/R_k(0) [gamma=0]  exact left eigvec\n")
    w2 = gate(2)     # lock convention (Wilson: 6/6)
    print()
    w3 = gate(3)     # the requested 18-member gate
    print()
    print(f"# VERDICT: L=2 {'6/6 PASS' if w2<=1e-12 else 'FAIL'} (convention locked); "
          f"L=3 {'18/18 PASS' if w3<=1e-12 else 'FAIL'} (worst {w3:.1e}).")


if __name__ == "__main__":
    main()
