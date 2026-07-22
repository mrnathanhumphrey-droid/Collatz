"""
PROBE R25 -- THE SPECTRAL GAP (subcritical). Reuses the renewal builder, now DEEP (r->14).

Wilson's redirect:
 (A) R24 def audit: amplitude/(eps*X_inf) must = 8/(3+2eps) exactly [since 1-rho=8eps/(3+2eps)]. R24 was off because
     amplitude used Y_inf=1+Sum S_r while eps*X_inf used X_inf=Sum S_r -- the X0=1 term. Recompute consistently.
 (B) DROP X_inf. Use C(lam) := lim_r S_r/rho^r  (per-r, no tail, no truncation). rho=3(1-lam)/(1+lam)=q*Sum p_v^2
     is the exact LEADING eigenvalue (collision rate); C is its spectral-projection weight. rho(1/2)=1 => C(1/2)=S_inf,
     so the theorem is the BOUNDARY VALUE C(1/2)=7/15 of a function analytic on (1/2,1).
 (C) THE GATEKEEPER: C extends continuously to lam=1/2 iff the SPECTRAL GAP survives, iff |lam2|/rho stays < 1 as
     eps->0. Measure |lam2|/rho from S_r/rho^r's deviation from its plateau. (The critical period-9 envelope~0.98 is
     what a subdominant complex pair near the circle looks like -- the live risk.)
 (D) Richardson in eps on the PLATEAU C(lam), not on X_inf.

Deep build: p_v=(1-lam)lam^{v-1} is negligible past v~70, so truncate the v-loop (exact to float) -> reach r=14.
"""
import os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np


def build_mu_qf(arr, k, q, lam, tol=1e-18):
    """Direct renewal in FLOAT (numpy): iterate v=1,2,... (p_v=(1-lam)lam^{v-1}), truncate when negligible."""
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
    S = {r: Y[r] - Y[r - 1] for r in range(1, RMAX + 1)}
    return S, Y


def aitken_last(seq):
    """Aitken delta^2 from the last triple; returns estimate."""
    a, b, c = seq[-3], seq[-2], seq[-1]
    d2 = c - 2 * b + a
    return (a - (b - a) ** 2 / d2) if d2 != 0 else c


def main():
    print("# PROBE R25 -- THE SPECTRAL GAP (subcritical), deep build r->14.\n")
    q = 3; RMAX = 14
    EPS = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002]

    # ---------- R25-A: definition audit ----------
    print("## R25-A  DEFINITION AUDIT: 1-rho = 8eps/(3+2eps); amplitude/(eps*X_inf) must = 8/(3+2eps)")
    print(f"   {'eps':>6} {'rho=3(1-l)/(1+l)':>16} {'1-rho':>10} {'8eps/(3+2eps)':>14} {'match':>7}")
    for eps in EPS:
        lam = 0.5 + eps; rho = 3 * (1 - lam) / (1 + lam)
        chk = 8 * eps / (3 + 2 * eps)
        print(f"   {eps:>6} {rho:>16.8f} {1-rho:>10.6f} {chk:>14.6f} {'OK' if abs(1-rho-chk)<1e-12 else 'DEV':>7}")
    print("   (R24 discrepancy CAUSE: amplitude used Y_inf=1+Sum S_r, eps*X_inf used X_inf=Sum S_r; the X0=1 term.")
    print("    Fix: amplitude := (1-rho)*X_inf; then amplitude/(eps*X_inf)=(1-rho)/eps=8/(3+2eps) exactly. But NOTE")
    print("    (1-rho)X_inf = C*rho + subdominant corrections, NOT clean C -- which is why we DROP X_inf for (B).)\n")

    data = {}
    for eps in EPS:
        lam = 0.5 + eps; rho = 3 * (1 - lam) / (1 + lam)
        S, Y = shells(q, lam, RMAX)
        p = {r: S[r] / rho ** r for r in range(1, RMAX + 1)}      # S_r/rho^r -> C(lam)
        data[eps] = (lam, rho, S, Y, p)

    # ---------- R25-B: the plateau ----------
    print("## R25-B  THE PLATEAU  C(lam)=lim_r S_r/rho^r  (per-r, NO tail; theorem = C(1/2)=7/15)")
    Cvals = {}
    for eps in EPS:
        lam, rho, S, Y, p = data[eps]
        C_est = aitken_last([p[r] for r in range(1, RMAX + 1)])
        Cvals[eps] = C_est
        tail = "  ".join(f"{p[r]:.5f}" for r in range(RMAX - 5, RMAX + 1))
        print(f"   eps={eps:<5} rho={rho:.5f}: p_r (r={RMAX-5}..{RMAX}) = {tail}   Aitken C={C_est:.6f}")
    print("   [does p_r flatten (plateau exists = gap survives) and C(lam) -> 7/15=0.466667 as eps->0?]\n")

    # ---------- R25-C: THE GATEKEEPER |lam2|/rho ----------
    print("## R25-C  THE GATEKEEPER: |lam2|/rho vs eps  (subdominant/dominant eigenvalue ratio)")
    print("   from d_r = p_{r+1}-p_r ~ C2(lam2/rho)^r*(lam2/rho-1): ratio |d_{r+1}/d_r| -> |lam2|/rho")
    print(f"   {'eps':>6} {'rho':>8} {'d_r ratios (last 6)':>44} {'|lam2|/rho (geom, last 5)':>26} {'osc?':>5}")
    gate = {}
    for eps in EPS:
        lam, rho, S, Y, p = data[eps]
        d = {r: p[r + 1] - p[r] for r in range(1, RMAX)}
        rr = [r for r in range(1, RMAX)]
        ratios = [d[rr[i + 1]] / d[rr[i]] if d[rr[i]] != 0 else float('nan') for i in range(len(rr) - 1)]
        # geometric rate of |d| over last 5 steps
        k = 5
        rate = (abs(d[RMAX - 1]) / abs(d[RMAX - 1 - k])) ** (1.0 / k) if abs(d[RMAX - 1 - k]) > 0 else float('nan')
        osc = any(ratios[i] < 0 for i in range(max(0, len(ratios) - 6), len(ratios)))
        gate[eps] = rate
        rshow = "  ".join(f"{x:+.3f}" for x in ratios[-6:])
        print(f"   {eps:>6} {rho:>8.5f} {rshow:>44} {rate:>26.5f} {'YES' if osc else 'no':>5}")
    print("   [GATEKEEPER: |lam2|/rho bounded < 1 as eps->0 => gap survives, C continuous at 1/2, ROUTE CLOSES.")
    print("    -> 1 => gap shuts at the point of interest; osc (sign flips) => complex pair = the period-9 risk.]\n")

    # ---------- R25-D: Richardson in eps on the plateau ----------
    print("## R25-D  RICHARDSON IN eps on the PLATEAU C(lam) (not on X_inf)")
    es = EPS[:]
    print(f"   {'eps':>6} {'C(lam)=plateau':>15}")
    for eps in es:
        print(f"   {eps:>6} {Cvals[eps]:>15.6f}")
    # linear-in-eps Richardson on the two smallest eps
    e1, e2 = es[-1], es[-2]
    C0_lin = Cvals[e1] + (Cvals[e1] - Cvals[e2]) / (e2 - e1) * e1
    # Aitken on the C(eps) sequence (last triple)
    C0_ait = aitken_last([Cvals[e] for e in es])
    print(f"   linear-in-eps extrap (eps={e1},{e2}) -> C(1/2) = {C0_lin:.6f}")
    print(f"   Aitken on C(eps) sequence        -> C(1/2) = {C0_ait:.6f}   [target 7/15 = {7/15:.6f}]")
    print("   [CAVEAT: valid ONLY if (C) shows the gap survives; if |lam2|/rho->1 the plateau C(lam) is ill-defined near 1/2.]")


if __name__ == "__main__":
    main()
