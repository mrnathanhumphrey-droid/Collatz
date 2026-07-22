"""
PROBE R30 -- THE DIGIT PROBABILITY. Recast the theorem as conditional-digit equidistribution.

Wilson's reframe: the requirement is NOT sqrt-cancellation. The theorem needs sum_r |A_r(m)| < infty, so
|A_r(m)| = O(r^{-1-delta}) suffices -- and it is observed GEOMETRIC. The object:

  A_r(m) = sum_{ord chi = 3^r} |nu_hat(chi)|^2 * chi(4^{-m})     (m-th Fourier coeff of spectral measure sigma_r)

Digit form. With p_r := Pr[X' ≡ 4^{-m} X mod 3^{r+1}], gamma_r := 3^r p_r, q_r(m) := p_r/p_{r-1}:
    A_r(m) = gamma_{r-1}(tau_m) * (3 q_r(m) - 1),      so   A_r(m) -> 0  <=>  q_r(m) -> 1/3.
The next 3-adic digit of the ratio becomes asymptotically uniform, conditional on agreement so far.

TWO INDEPENDENT PATHS welded here (never before directly, for m>=1):
  collision side:  gamma_r(tau_m) = 3^r Pr[X'≡4^{-m}X mod 3^{r+1}]   (probe_gamma_R9.gamma; collision-partner)
  character side:  A_r(m) = sum_u g_r(u) c_{3^r}((u-m) mod 3^r)       (probe_charledger_R10; dlog + Ramanujan)
  identity:        A_r(m) = gamma_r - gamma_{r-1} = gamma_{r-1}(3q_r-1).   [m=0 => S_r, the R10-B/R28 weld]

A: measure q_r(m), 3q_r-1, gamma_{r-1}(3q_r-1) vs banked A_r(m). GATE: last two agree exactly (the identity).
B: rate (3q_{r+1}-1)/(3q_r-1) per m -- same across m? ~|lam2|~1/2? (measurement, NO fit)
C: full next-digit distribution (3 values) conditional on agreement -- where does the excess over uniform sit?
D: is q_r sensitive to WHICH digits matched or only HOW MANY? (depth-Markov => scalar recursion => provable gap)
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_engine_R7 as R7
from probe_gamma_R9 import gamma, tau, v3
from probe_charledger_R10 import dlog_table, autocorr_dlog
from probe_gapop_R28 import build_nu

MS = [1, 2, 3, 4]                              # ratio indices tau_m ~ 4^{-m}
R_EXACT = 7                                     # exact mu block (char side + exact gate); float beyond
R_FLOAT = 10


def A_char(g, r, m):
    """Character-side A_r(m) = sum_u g_r(u) c_{3^r}((u-m) mod 3^r).  Exact rational. (m=0 => S_r.)
       g = cached autocorr_dlog(mu_r, r)."""
    N = 3 ** r
    return sum(g[u] * R7.cram(r, (u - m) % N) for u in range(N))


def p_from_nu(dense, M, m):
    """float p_r = Pr[X' ≡ 4^{-m} X mod M], M=3^{r+1}, from dense law of X_r."""
    fac = pow(pow(4, -1, M), m, M)             # 4^{-m} mod M
    idx = (np.arange(M) * fac) % M
    return float(np.sum(dense * dense[idx]))


def digit_dist_exact(mu_r, r, m):
    """Exact next-digit distribution conditional on agreement mod 3^r.
       X=1+3a (a in Z/3^r) -> X mod 3^{r+1}; over iid pair, digit j=((X'-4^{-m}X)/3^r) mod 3 in {0,1,2}.
       Returns (agree_mass, [d0,d1,d2]) as Fractions."""
    M1, M2 = 3 ** r, 3 ** (r + 1)
    fac = pow(pow(4, -1, M2), m, M2)           # 4^{-m} mod 3^{r+1}
    massXp = {}                                 # residue (X' mod M2) -> mass
    for a, pa in mu_r.items():
        res = (1 + 3 * (a % M1)) % M2
        massXp[res] = massXp.get(res, F(0)) + pa
    dd = [F(0), F(0), F(0)]
    for a, pa in mu_r.items():
        Y = (fac * ((1 + 3 * (a % M1)) % M2)) % M2
        for j in range(3):
            tgt = (Y + j * M1) % M2
            m2 = massXp.get(tgt)
            if m2:
                dd[j] += pa * m2
    agree = dd[0] + dd[1] + dd[2]
    return agree, dd


def q_by_prefix(mu_r, r, m, modk):
    """q_r conditioned on the prefix X mod 3^modk (WHICH digits). Returns {c: (q_class, agree_mass)}."""
    M1, M2 = 3 ** r, 3 ** (r + 1)
    Mk = 3 ** modk
    fac = pow(pow(4, -1, M2), m, M2)
    massXp = {}
    for a, pa in mu_r.items():
        res = (1 + 3 * (a % M1)) % M2
        massXp[res] = massXp.get(res, F(0)) + pa
    cls = {}                                    # c -> [d0_mass, agree_mass]
    for a, pa in mu_r.items():
        X = (1 + 3 * (a % M1)) % M2
        c = X % Mk
        Y = (fac * X) % M2
        d0 = massXp.get(Y % M2, F(0))           # j=0 (match)
        ag = d0 + massXp.get((Y + M1) % M2, F(0)) + massXp.get((Y + 2 * M1) % M2, F(0))
        e = cls.setdefault(c, [F(0), F(0)])
        e[0] += pa * d0; e[1] += pa * ag
    return {c: (float(v[0] / v[1]) if v[1] else float('nan'), v[1]) for c, v in cls.items()}


def main():
    print("# PROBE R30 -- THE DIGIT PROBABILITY.  q_r(m) -> 1/3  <=>  A_r(m) -> 0.\n")
    t0 = time.time()
    mu = {1: R7.mu1()}
    for k in range(2, R_EXACT + 1):
        mu[k] = R7.build_mu(mu[k - 1], k)
        print(f"   built mu_{k}: |supp|={len(mu[k])}  ({time.time()-t0:.1f}s)")
    print()

    # exact gamma_r(tau_m) and character-side A_r(m)  [cache autocorr_dlog once per level]
    gam = {m: {r: gamma(mu[r], r, tau(m, r)) for r in range(1, R_EXACT + 1)} for m in MS}
    gcache = {}
    for r in range(1, R_EXACT + 1):
        gcache[r] = autocorr_dlog(mu[r], r)
        print(f"   autocorr g_{r} cached  ({time.time()-t0:.1f}s)")
    Ach = {m: {r: A_char(gcache[r], r, m) for r in range(1, R_EXACT + 1)} for m in MS}
    print()

    # float extension via nu
    nus = build_nu(0.5, R_FLOAT)
    dense = {r: np.zeros(3 ** (r + 1)) for r in range(1, R_FLOAT + 1)}
    for r in range(1, R_FLOAT + 1):
        for X, w in nus[r].items():
            dense[r][X] = w
    gamf = {m: {r: 3 ** r * p_from_nu(dense[r], 3 ** (r + 1), m) for r in range(1, R_FLOAT + 1)} for m in MS}

    # ============ R30-A ============
    print("## R30-A  MEASURE q_r(m) + IDENTITY GATE  [gamma_{r-1}(3q_r-1)  vs  banked A_r(m) (char side)]")
    print("   exact block r=2..%d; PRE-REG: the two rightmost columns agree EXACTLY (a miss = derivation wrong).\n" % R_EXACT)
    okA = True
    for m in MS:
        print(f"   --- m={m}  (tau_{m}, v3={v3(m)};  f(4^-{m}) ~ lim gamma_r) ---")
        print(f"   {'r':>2} {'gamma_r':>10} {'q_r=p/p':>9} {'3q_r-1':>10} {'g_{r-1}(3q-1)':>13} {'A_r(m) char':>13} {'gate':>6}")
        for r in range(2, R_EXACT + 1):
            gr, grm1 = gam[m][r], gam[m][r - 1]
            q = gr / (3 * grm1)                          # p_r/p_{r-1}
            three_q = 3 * q - F(1)
            digit_side = grm1 * three_q                  # = gamma_r - gamma_{r-1}
            achar = Ach[m][r]
            good = (digit_side == achar)
            okA = okA and good
            print(f"   {r:>2} {float(gr):>10.6f} {float(q):>9.6f} {float(three_q):>+10.6f} "
                  f"{float(digit_side):>+13.6f} {float(achar):>+13.6f} {'PASS' if good else 'FAIL':>6}")
            if not good:
                print(f"      MISMATCH r={r}: digit_side={digit_side}\n                A_char   ={achar}")
        print()
    print(f"   => R30-A {'GATE PASS -- identity A_r(m)=gamma_{r-1}(3q_r-1) certified exactly, m=1..4, r=2..%d' % R_EXACT if okA else 'FAIL -- derivation wrong (mismatch above)'}\n")

    # convergence gamma_r -> f, and q_r -> 1/3, into the float range
    print("   gamma_r(tau_m) -> f(4^-m)  and  q_r -> 1/3 ?  (exact r<=%d, float beyond; Wilson anchor A_r(1)~0.014,0.006 @ r~8,9)" % R_EXACT)
    for m in MS:
        cells = []
        for r in range(2, R_FLOAT + 1):
            g = gamf[m][r]; gm1 = gamf[m][r - 1]
            q = g / (3 * gm1) if gm1 else float('nan')
            cells.append(f"r{r}:g={g:.4f},q={q:.4f}")
        print(f"   m={m}: " + "  ".join(cells))
    print()

    # ============ R30-B ============
    print("## R30-B  RATE + m-INDEPENDENCE  ratio (3q_{r+1}-1)/(3q_r-1) per m  (measurement, NO fit)")
    print("   [same across m => one mixing chain; value ~|lam2|~0.5 ?]  (exact where available, float beyond)")
    def three_q(g, m, r):
        return 3 * (g[m][r] / (3 * g[m][r - 1])) - (F(1) if isinstance(g[m][r], F) else 1.0)
    for m in MS:
        cells = []
        for r in range(3, R_FLOAT):
            # use exact in exact block, else float
            if r + 1 <= R_EXACT:
                num = float(three_q(gam, m, r + 1)); den = float(three_q(gam, m, r))
            else:
                num = float(three_q(gamf, m, r + 1)); den = float(three_q(gamf, m, r))
            cells.append(f"r{r}->{r+1}:{num/den:+.4f}" if den else f"r{r}:na")
        print(f"   m={m}: " + "  ".join(cells))
    # also A_{r+1}/A_r directly (char side, exact)
    print("   [ref] A_{r+1}(m)/A_r(m) char-side exact:")
    for m in MS:
        cells = [f"r{r}->{r+1}:{float(Ach[m][r+1]/Ach[m][r]):+.4f}" for r in range(1, R_EXACT) if Ach[m][r] != 0]
        print(f"   m={m}: " + "  ".join(cells))
    print()

    # ============ R30-C ============
    print("## R30-C  NEXT-DIGIT DISTRIBUTION (3 values) conditional on agreement mod 3^r  [uniform=(1/3,1/3,1/3)]")
    print("   digit j=((X'-4^{-m}X)/3^r) mod 3; j=0 is the match (prob q_r). Where does the excess sit?")
    for m in MS:
        print(f"   --- m={m} ---")
        for r in range(2, min(R_EXACT, 6) + 1):
            agree, dd = digit_dist_exact(mu[r], r, m)
            probs = [float(x / agree) for x in dd]
            skew = [p - 1 / 3 for p in probs]
            print(f"   r={r}: P(d=0,1,2)=({probs[0]:.5f},{probs[1]:.5f},{probs[2]:.5f})  "
                  f"excess=({skew[0]:+.5f},{skew[1]:+.5f},{skew[2]:+.5f})  q_r=P(d=0)={probs[0]:.5f}")
        print()

    # ============ R30-D ============
    print("## R30-D  CONDITIONING CHECK: does q_r depend on WHICH digits matched, or only HOW MANY (depth)?")
    print("   q_r conditioned on prefix X mod 3^k (k=1: X mod3 trivial; k=2: X mod9 in {1,4,7}; k=3: X mod27).")
    print("   [q_r constant across prefix classes => depth-Markov => scalar recursion => provable gap.]")
    for m in (1, 2):
        print(f"   --- m={m} ---")
        for r in range(3, min(R_EXACT, 6) + 1):
            for modk in (2, 3):
                byc = q_by_prefix(mu[r], r, m, modk)
                items = sorted(byc.items())
                qs = [q for _, (q, _) in items if q == q]  # drop nan
                spread = (max(qs) - min(qs)) if qs else float('nan')
                shown = "  ".join(f"X%{3**modk}={c}:q={q:.5f}" for c, (q, w) in items[:6])
                print(f"   r={r} mod3^{modk}: {shown}   [spread={spread:.2e}]")
        print()
    print(f"   (total {time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
