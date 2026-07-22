"""
PROBE R26 -- THE DIRECT EIGENVALUE (lam2 via Prony / matrix-pencil). Reuses the R25 deep renewal builder.

Wilson's redirect + prediction:
 - R25-C reading REVERSED: clean points |lam2|/rho = 0.69(eps=.1), 0.57(eps=.05) DECREASE => moving AWAY from 1.
   The "->1" was unresolved small-eps artifacts read as a trend (#32/#40 failure mode). Corrected.
 - Derived: for p_v=(1-lam)lam^{v-1}, pair-gap P(d)=(1-lam)lam^|d|/(1+lam); leading eig rho=3P(0)=3(1-lam)/(1+lam)
   (=confirmed decay rate); subdominant = gap-+-2 channel => |lam2|/rho = 2*lam^2 (~5% below => ~0.95*2lam^2).
   At criticality 2lam^2 = 1/2: GAP SURVIVES, C continuous at 1/2, route closes. Also P(+-2)/P(0)=1/2 at lam=1/2.

The exact second-moment operator is on growing spaces (R16 crux), so we read lam2 as the SUBDOMINANT eigenvalue of
the finite COMPANION matrix of the linear recurrence that S_r (subcritical) / eps_r=S_r-7/15 (critical) satisfies --
a well-posed finite eigenproblem, no r->inf extrapolation, valid AT lam=1/2.

 E1 (gate): |lam2|/rho must reproduce 0.69 (eps=.1), 0.57 (eps=.05). Miss => method wrong, stop.
 E2: |lam2|/rho ~ 0.95*2lam^2 across eps, -> ~1/2 as eps->0.
 E3: run AT lam=1/2 (rho=1); read |lam2| directly, compare 0.493/0.503 and 1/2.
 E4: arg(lam2) at lam=1/2 -> period 2pi/|arg| (compare the period-9 chased since R13).
"""
import os, sys, math, json, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np


def build_mu_qf(arr, k, q, lam, tol=1e-18):
    M = q ** k; inv2 = pow(2, -1, M)
    a_idx = np.nonzero(arr)[0]; a_val = arr[a_idx]
    base = (1 + q * a_idx) % M
    mu = np.zeros(M); u = inv2; v = 1
    while (1 - lam) * lam ** (v - 1) > tol:
        wv = (1 - lam) * lam ** (v - 1)
        mu += np.bincount((u * base) % M, weights=wv * a_val, minlength=M)
        u = (u * inv2) % M; v += 1
    return mu


def shells(q, lam, RMAX):
    arr = np.array([1.0]); Y = {0: 1.0}
    for k in range(1, RMAX + 1):
        arr = build_mu_qf(arr, k, q, lam)
        Y[k] = q ** k * float(np.sum(arr * arr))
    return {r: Y[r] - Y[r - 1] for r in range(1, RMAX + 1)}


def prony_roots(x, M):
    """Companion-matrix (Prony) eigenvalues of a signal ~ sum of M exponentials. Returns roots sorted by |.| desc."""
    x = np.asarray(x, float); N = len(x)
    A = np.array([[x[n - k] for k in range(1, M + 1)] for n in range(M, N)])
    b = np.array([x[n] for n in range(M, N)])
    a, *_ = np.linalg.lstsq(A, b, rcond=None)
    roots = np.roots(np.concatenate([[1.0], -a]))
    return sorted(roots, key=lambda z: -abs(z))


def main():
    print("# PROBE R26 -- THE DIRECT EIGENVALUE (lam2 via Prony/companion matrix).\n")
    q = 3
    hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                        'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
    EPS_exact = {int(k): float(v['num']) / float(v['den']) for k, v in hist.items()}   # eps_k = S_k - 7/15, exact k<=8

    # ================= E3 / E4 : AT CRITICALITY lam=1/2 =================
    print("## E3/E4  CRITICAL lam=1/2 (rho=1): read |lam2| and arg(lam2) directly from eps_r = S_r - 7/15")
    RC = 16
    Sc = shells(q, 0.5, RC)
    eps = {r: Sc[r] - 7 / 15 for r in range(1, RC + 1)}
    # cross-check float vs exact eps (k<=8)
    dev = max(abs(eps[k] - EPS_exact[k]) for k in range(1, 9))
    print(f"   float eps_r vs exact ledger (k<=8): max dev = {dev:.2e}  [{'OK' if dev<1e-6 else 'FLOAT ISSUE'}]")
    print(f"   eps_r (r=1..{RC}): " + "  ".join(f"{eps[r]:+.3e}" for r in range(1, RC + 1)))
    for M in (2, 3, 4):
        roots = prony_roots([eps[r] for r in range(3, RC + 1)], M)   # skip r=1,2 transient
        lead = roots[0]
        per = (2 * math.pi / abs(cmath.phase(lead))) if abs(cmath.phase(lead)) > 1e-9 else float('inf')
        print(f"   M={M}: dominant lam2 = {lead.real:+.5f}{lead.imag:+.5f}i  |lam2|={abs(lead):.5f}  "
              f"arg={math.degrees(cmath.phase(lead)):+.2f} deg -> period={per:.3f}")
        if M == 3:
            print(f"        (all roots: " + ", ".join(f"{r.real:+.3f}{r.imag:+.3f}i" for r in roots) + ")")
    print(f"   [PRE-REG: |lam2| vs R18-A exact ratio 0.493/0.503 and 1/2; period vs 9.  gap survives iff |lam2|<1.]\n")

    # ================= E1 / E2 : SUBCRITICAL =================
    print("## E1/E2  SUBCRITICAL: leading root = rho (gate), subdominant = lam2; |lam2|/rho vs 0.95*2lam^2, R25-C")
    print(f"   {'eps':>6} {'lam':>6} {'rho_pred':>9} {'rho(Prony)':>11} {'|lam2|':>8} {'|lam2|/rho':>10} "
          f"{'0.95*2lam^2':>11} {'2lam^2':>8} {'R25-C':>7}")
    R25C = {0.1: 0.69, 0.05: 0.57}
    RS = 14
    for eps_ in (0.1, 0.05, 0.02, 0.01, 0.005):
        lam = 0.5 + eps_; rho = 3 * (1 - lam) / (1 + lam)
        S = shells(q, lam, RS)
        roots = prony_roots([S[r] for r in range(3, RS + 1)], 3)     # skip early transient
        # leading root ~ rho (real positive); subdominant = next by magnitude
        rho_p = roots[0].real
        lam2 = roots[1]
        ratio = abs(lam2) / rho_p
        r25 = f"{R25C[eps_]:.2f}" if eps_ in R25C else "-"
        print(f"   {eps_:>6} {lam:>6.3f} {rho:>9.5f} {rho_p:>11.5f} {abs(lam2):>8.5f} {ratio:>10.5f} "
              f"{0.95*2*lam*lam:>11.5f} {2*lam*lam:>8.5f} {r25:>7}")
    print("   [E1 gate: |lam2|/rho ~ 0.69, 0.57 at eps=.1,.05. E2: -> ~1/2 as eps->0 (via 0.95*2lam^2).")
    print("    Prony rho(leading) must also match rho_pred=3(1-lam)/(1+lam).]")


if __name__ == "__main__":
    main()
