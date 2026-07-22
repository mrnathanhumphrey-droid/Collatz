"""
PROBE R24 -- THE SUBCRITICAL SCALING LAW. Reuses R7/R23 builder. Hold q=3, step off criticality to lam=1/2+eps.

Rationale (Wilson): S_inf=2 Sum_chi |nu_hat|^2 w(chi) is CONDITIONALLY convergent at criticality (terms don't
shrink; only cancellation converges) => no finite prefix determines it, extrapolators impose a model (the 0.72-0.81
scatter = model misspecification). Subcritically (lam>1/2) the SAME sum is POSITIVE-TERM, ABSOLUTELY convergent,
clean geometric decay S_r ~ C(lam) rho^r, rho=3(1-lam)/(1+lam)<1. Float-safe, no cancellation.

  X_inf(lam) := Sum_{r>=1} S_r = (lim_r 3^r ||mu_r||^2) - 1.   1-rho ~ (8/3)eps.
  PRE-REG: eps * X_inf(1/2+eps) -> 7/40 = 0.175   (equivalently amplitude C(lam)=(1-rho)X_inf -> 7/15).
Two INDEPENDENT extrapolators must AGREE (that agreement is the trust, unlike the critical scatter):
  (i) Aitken delta^2 on the partial sums Y_r=3^r||mu_r||^2;
  (ii) exact-rho geometric tail: Y_inf = Y_R + S_R * rho/(1-rho).
Support = 2*3^{r-1} (lam-independent), so build caps at r~10; the extrapolation does the rest -- LEGITIMATE here
because the decay is a known single-mode geometric, NOT the critical oscillation.
Trap avoided (Wilson): extrapolate X_inf, NOT f(tau_m) -- subcritically X_inf is finite so f has NO identity
singularity; the singularity is born exactly at lam=1/2.
"""
import os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np


def ord2(M):
    o, x = 1, 2 % M
    while x != 1:
        x = (x * 2) % M; o += 1
    return o


def build_mu_qf(arr, k, q, lam):
    """Generalized renewal in FLOAT via numpy (bincount). arr: dense length q^{k-1}."""
    M = q ** k; inv2 = pow(2, -1, M); ordv = ord2(M)
    a_idx = np.nonzero(arr)[0]; a_val = arr[a_idx]
    base = np.array([(1 + q * int(a)) % M for a in a_idx], dtype=np.int64)
    mu = np.zeros(M); u = inv2; dnm = 1 - lam ** ordv
    for j in range(1, ordv + 1):
        wv = (1 - lam) * lam ** (j - 1) / dnm
        idx = (u * base) % M
        mu += np.bincount(idx, weights=wv * a_val, minlength=M)
        u = (u * inv2) % M
    return mu


def aitken(seq):
    L = []
    for i in range(len(seq) - 2):
        d1 = seq[i + 1] - seq[i]; d2 = seq[i + 2] - 2 * seq[i + 1] + seq[i]
        L.append(seq[i] - d1 * d1 / d2 if d2 != 0 else float('nan'))
    return L


def main():
    print("# PROBE R24 -- THE SUBCRITICAL SCALING LAW. q=3, lam=1/2+eps.\n")
    q = 3; RMAX = 10
    EPS = [0.1, 0.05, 0.02, 0.01]
    summary = []
    for eps in EPS:
        lam = 0.5 + eps
        rho = 3 * (1 - lam) / (1 + lam)          # exact predicted decay rate
        arr = np.array([1.0])
        Y = {0: 1.0}
        for k in range(1, RMAX + 1):
            arr = build_mu_qf(arr, k, q, lam)
            Y[k] = q ** k * float(np.sum(arr * arr))
        S = {r: Y[r] - Y[r - 1] for r in range(1, RMAX + 1)}
        print(f"## eps={eps}  lam={lam:.3f}  rho_pred=3(1-lam)/(1+lam)={rho:.6f}  (1-rho={1-rho:.5f}, (8/3)eps={8/3*eps:.5f})")
        print(f"   {'r':>2} {'Y_r=3^r||mu||^2':>16} {'S_r':>12} {'S_r/S_{r-1}':>12}")
        for r in range(1, RMAX + 1):
            rat = S[r] / S[r - 1] if r > 1 and S[r - 1] != 0 else float('nan')
            print(f"   {r:>2} {Y[r]:>16.8f} {S[r]:>12.6e} {rat:>12.6f}")
        # extrapolator (i): Aitken on Y_r
        La = [x for x in aitken([Y[r] for r in range(1, RMAX + 1)]) if not math.isnan(x)]
        Y_aitken = La[-1] if La else float('nan')
        spread = (max(La[-4:]) - min(La[-4:])) if len(La) >= 2 else float('nan')
        # extrapolator (ii): exact-rho geometric tail from r=RMAX
        Y_tail = Y[RMAX] + S[RMAX] * rho / (1 - rho)
        # measured-rho tail (cross-check)
        rho_meas = S[RMAX] / S[RMAX - 1]
        Y_tailm = Y[RMAX] + S[RMAX] * rho_meas / (1 - rho_meas)
        agree = abs(Y_aitken - Y_tail)
        Xinf = Y_tail - 1                          # Wilson's X_inf = Sum_r S_r
        print(f"   extrap Y_inf: Aitken={Y_aitken:.6f} (spread {spread:.1e}); rho-tail(exact)={Y_tail:.6f}; "
              f"rho-tail(meas rho={rho_meas:.4f})={Y_tailm:.6f}  | agree |A-tail|={agree:.1e}")
        print(f"   X_inf=Sum_r S_r = Y_inf-1 = {Xinf:.6f}")
        print(f"   >>> eps*X_inf = {eps*Xinf:.6f}   [PRE-REG -> 7/40 = 0.175000]   "
              f"amplitude (1-rho)*Y_inf = {(1-rho)*Y_tail:.6f} [-> 7/15=0.466667]\n")
        summary.append((eps, lam, rho, Y_tail, Xinf, eps * Xinf, (1 - rho) * Y_tail, spread, agree))

    print("## R24 SCALING-LAW SUMMARY (the test: eps*X_inf -> 7/40 as eps->0)")
    print(f"   {'eps':>6} {'rho':>8} {'X_inf':>10} {'eps*X_inf':>11} {'(1-rho)Yinf':>12} {'Aitken spread':>14} {'estim agree':>12}")
    for eps, lam, rho, Yt, Xinf, eXi, amp, sp, ag in summary:
        print(f"   {eps:>6} {rho:>8.5f} {Xinf:>10.5f} {eXi:>11.6f} {amp:>12.6f} {sp:>14.1e} {ag:>12.1e}")
    print("   [PRE-REG: eps*X_inf column -> 0.175000; amplitude column -> 0.466667. Extrapolators must AGREE")
    print("    (small spread/agree = Wilson's 'calm subcritical' promise verified; large = model still misspecified).]")
    print("   [CAVEAT: test is 'does C(lam)->7/15', theorem approached from OUTSIDE; needs continuity of C at criticality (unproved).]")


if __name__ == "__main__":
    main()
