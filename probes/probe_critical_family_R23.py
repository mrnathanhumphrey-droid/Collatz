"""
PROBE R23 -- THE CRITICAL FAMILY + f(tau) EXTRAPOLATION. Reuses R7/R9. Tests Wilson's decisive q-sweep and the
skeptic's convergence check.

Generalized renewal at (q, lambda), R7-VALIDATED build form (byte-gate vs banked q=3 confirms it):
  mu_k on Z/q^k, index r = 2^{-v}(1+q a) mod q^k, a~mu_{k-1}, v~Geom(lambda), p_v=(1-lam)lam^{v-1}.
  Class weight for v==j mod ord_{q^k}(2): (1-lam)lam^{j-1}/(1-lam^ord).
  [The certified R13-D orbit form is X'=1+3*2^{-v}X on 1+3Z_3 (the audit's X=2^{-a}(1+3X') is a WRONG rearrangement,
   not used); this code uses the byte-gated R7 builder, so the q=3 sanity gate is the guarantor.]
Shell: X_k = q^k ||mu_k||^2 (=gamma_k(0)); S_k = X_k - X_{k-1}; S_inf = lim S_k.
Criticality lam_c(q)=(q-1)/(q+1). Conjecture 2 (critical family):
  S_inf(q) = M4/M3 = (3q^2+1)/(2q(q^2+1)).

THE DECISIVE TEST (Wilson): S_inf(q)~3/(2q) is SMOOTH but every finite S_r inherits ord_q(2) which swings wildly
  (ord=2,4,3,10,12 at q=3,5,7,11,13; S_1=0.667,0.492,1.459,... non-monotone). If S_k -> smooth prediction despite
  ord swinging 3->12, Mersenne structure is real; if S_k tracks ord_q(2), Conjecture 2 is DEAD.
  PRE-REG: q=5->19/65=0.292308, q=7->37/175=0.211429, q=11->91/671=0.135618, q=13->127/1105=0.114932.
  (lam-power denominators explode => EXACT infeasible for q>=5; float build, exact only for the q=3 byte-gate.)
Functional-equation move (reviewer #2) is RETIRED: Wilson proved rho->rho' is not well-defined (nu and nu_g have
  identical ratio laws but different transports) -- no self-map exists. Not pursued here.

R23-A critical family (forced): q=3 byte-gate + q=5,7,11,13 float -> smooth-vs-ord decision.
R23-B f(tau_1),f(tau_2) exact r<=7 + float push r->10 + Aitken; report extrapolated value AND spread across orders.
R23-C rational hunt (GATED on B locking: denom<200 + must survive; else premature).
"""
import os, sys, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_engine_R7 as R7
import probe_gamma_R9 as R9


def build_mu_q_lam(q, mu_prev, k, lam):
    """Exact generalized renewal (Fraction lam)."""
    M = q ** k; inv2 = pow(2, -1, M)
    ordv, x = 1, 2 % M
    while x != 1:
        x = (x * 2) % M; ordv += 1
    denom = 1 - lam ** ordv
    mu = {}
    for j in range(1, ordv + 1):
        wv = (1 - lam) * lam ** (j - 1) / denom
        u = pow(inv2, j, M)
        for a, pa in mu_prev.items():
            r = (u * (1 + q * a)) % M
            mu[r] = mu.get(r, F(0)) + wv * pa
    return mu


def shells_exact(q, lam, K):
    mu = {0: {0: F(1)}}; X = {0: F(1)}; S = {}
    for k in range(1, K + 1):
        mu[k] = build_mu_q_lam(q, mu[k - 1], k, lam)
        X[k] = q ** k * sum(p * p for p in mu[k].values())
        S[k] = X[k] - X[k - 1]
    return S, mu


def ord2(M):
    o, x = 1, 2 % M
    while x != 1:
        x = (x * 2) % M; o += 1
    return o


def build_mu_qf(arr, k, q, lam):
    """Generalized renewal in FLOAT via numpy dense arrays. arr: dense length q^{k-1}."""
    M = q ** k; inv2 = pow(2, -1, M); ordv = ord2(M)
    a_idx = np.nonzero(arr)[0]; a_val = arr[a_idx]
    base = (1 + q * a_idx.astype(object)) % M          # object to avoid int64 overflow for big q^k
    mu = np.zeros(M); u = inv2; denom = 1 - lam ** ordv
    for j in range(1, ordv + 1):
        wv = (1 - lam) * lam ** (j - 1) / denom
        idx = np.array((u * base) % M, dtype=np.int64)
        np.add.at(mu, idx, wv * a_val)
        u = (u * inv2) % M
    return mu


def shells_float(q, lam, K):
    arr = np.array([1.0]); X = {0: 1.0}; S = {}
    for k in range(1, K + 1):
        arr = build_mu_qf(arr, k, q, lam)
        X[k] = q ** k * float(np.sum(arr * arr)); S[k] = X[k] - X[k - 1]
    return S


def Mk(q, lam, k):
    return (1 - lam) ** k / (1 - lam ** k)


def gamma_float(mu_arr, r, t):
    M = 3 ** r
    a = np.nonzero(mu_arr)[0]
    partner = (a + t * (1 + 3 * a)) % M
    return (3 ** r) * float(np.sum(mu_arr[a] * mu_arr[partner]))


def aitken(seq):
    """Aitken delta^2 extrapolation; returns list of L estimates."""
    L = []
    for i in range(len(seq) - 2):
        d1 = seq[i + 1] - seq[i]; d2 = seq[i + 2] - 2 * seq[i + 1] + seq[i]
        L.append(seq[i] - d1 * d1 / d2 if d2 != 0 else float('nan'))
    return L


def main():
    print("# PROBE R23 -- THE CRITICAL FAMILY + f(tau) EXTRAPOLATION.\n")

    # ================= R23-A =================
    print("## R23-A  CRITICAL FAMILY (forced): S_inf(q)=(3q^2+1)/(2q(q^2+1)); lam_c(q)=(q-1)/(q+1)")
    banked = dict(R7.S)
    # SANITY GATE: q=3 lam=1/2 EXACT byte-exact vs banked (guarantor of the generalized code)
    S3e, _ = shells_exact(3, F(1, 2), 5)
    okgate = all(S3e[k] == banked[k] for k in range(1, 6))
    print(f"   SANITY GATE q=3,lam=1/2 EXACT: S_k==banked k=1..5? {okgate}  (pred 28/60=7/15={float(F(7,15)):.6f})")
    # float builder cross-check at q=3
    S3f = shells_float(3, 0.5, 5)
    print(f"   float builder q=3: S_5(float)={S3f[5]:.6f} vs exact {float(S3e[5]):.6f}  "
          f"[{'float-builder OK' if abs(S3f[5]-float(S3e[5]))<1e-9 else 'FLOAT BUG'}]")
    # Wilson's hand values S_1(q): q=5->416/845=0.49231, q=7->1998/1369=1.45946
    QS = ((5, 2 / 3, 6, 19 / 65, 4), (7, 3 / 4, 5, 37 / 175, 3), (11, 5 / 6, 4, 91 / 671, 10), (13, 6 / 7, 4, 127 / 1105, 12))
    print(f"\n   {'q':>3} {'lam_c':>6} {'ord_q(2)':>8} {'S_1':>9} {'S_kmax':>9} {'target S_inf':>13} {'Aitken S_inf':>13} {'hit?':>6}")
    for q, lam, K, tgt, ordq in QS:
        Sq = shells_float(q, lam, K)
        seq = [Sq[k] for k in range(1, K + 1)]
        Ls = [x for x in aitken(seq) if not math.isnan(x)]
        ext = Ls[-1] if Ls else float('nan')
        hit = abs(ext - tgt) < 5e-3 and abs(Sq[K] - tgt) < abs(Sq[1] - tgt)
        print(f"   {q:>3} {lam:>6.4f} {ordq:>8} {Sq[1]:>9.5f} {Sq[K]:>9.5f} {tgt:>13.6f} {ext:>13.6f} {'YES' if hit else 'no':>6}")
        eps = "  ".join(f"S{k}-tgt={Sq[k]-tgt:+.2e}" for k in range(1, K + 1))
        print(f"       {eps}")
    print("   [DECISION: S_k tracks SMOOTH prediction (eps shrinks, Aitken->target) despite ord swinging 3..12 = Conj-2 alive;")
    print("    S_k tracks ord_q(2) (S_1 non-monotone persists) = Conj-2 DEAD. q=3 byte-gate makes the 4 numbers meaningful.]\n")

    # ================= R23-B =================
    print("## R23-B  f(tau_1),f(tau_2): exact r<=7 + FLOAT push + Aitken (skeptic's convergence check)")
    mu = {1: R7.mu1()}
    for k in range(2, 8):
        mu[k] = R7.build_mu(mu[k - 1], k)
    g1 = {}; g2 = {}
    for r in range(1, 8):
        g1[r] = float(R9.gamma(mu[r], r, R9.tau(1, r)))
        g2[r] = float(R9.gamma(mu[r], r, R9.tau(2, r)))
    # float push r=8..RMAX
    arr = np.zeros(3 ** 7)
    for a, p in mu[7].items():
        arr[a] = float(p)
    RMAX = 10
    for r in range(8, RMAX + 1):
        arr = build_mu_qf(arr, r, 3, 0.5)
        g1[r] = gamma_float(arr, r, R9.tau(1, r))
        g2[r] = gamma_float(arr, r, R9.tau(2, r))
        print(f"   [built r={r} float, supp={int(np.count_nonzero(arr))}]")
    print(f"   {'r':>2} {'gamma_r(tau_1)':>16} {'gamma_r(tau_2)':>16}")
    for r in range(1, RMAX + 1):
        tag = "" if r <= 7 else "  (float)"
        print(f"   {r:>2} {g1[r]:>16.8f} {g2[r]:>16.8f}{tag}")
    for lab, g in (("f(tau_1)", g1), ("f(tau_2)", g2)):
        seq = [g[r] for r in range(3, RMAX + 1)]      # skip early transient
        L = aitken(seq)
        Lc = [x for x in L if not math.isnan(x)]
        spread = (max(Lc) - min(Lc)) if Lc else float('nan')
        print(f"   {lab}: raw r={RMAX}={g[RMAX]:.6f}; Aitken L (r>=3): {[round(x,5) for x in L]}")
        print(f"        extrapolated ~ {Lc[-1]:.6f} (last), spread across orders = {spread:.2e}  "
              f"[{'LOCKED' if spread<1e-3 else 'NOT locked (limits may not be limits)'}]")
    print()

    # ================= R23-C =================
    print("## R23-C  RATIONAL HUNT (GATED on B locking; denom<200; must survive extension)")
    def hunt(x, maxden=200):
        fr = F(x).limit_denominator(maxden)
        return fr, abs(float(fr) - x)
    L1 = [x for x in aitken([g1[r] for r in range(3, RMAX + 1)]) if not math.isnan(x)]
    L2 = [x for x in aitken([g2[r] for r in range(3, RMAX + 1)]) if not math.isnan(x)]
    sp1 = (max(L1) - min(L1)) if L1 else 9; sp2 = (max(L2) - min(L2)) if L2 else 9
    for lab, Lc, sp, susp in (("f(tau_1)", L1, sp1, None), ("f(tau_2)", L2, sp2, F(10, 21))):
        if sp >= 1e-3:
            print(f"   {lab}: extrapolation spread {sp:.1e} >= 1e-3 -> NOT locked -> rational hunt PREMATURE (per pre-reg).")
        else:
            fr, err = hunt(Lc[-1])
            print(f"   {lab}: locked at {Lc[-1]:.6f}; best rational (den<200) = {fr} (err {err:.1e})"
                  + (f"; suspect {susp}={float(susp):.6f} gap {abs(float(susp)-Lc[-1]):.1e}" if susp else ""))
        if susp:
            print(f"        NOTE (Wilson): f(tau_2)=10/21 is SUSPECT -- level-2 quantity as a limit (R9-D/R20-D vacated this shape); test not believe.")


if __name__ == "__main__":
    main()
