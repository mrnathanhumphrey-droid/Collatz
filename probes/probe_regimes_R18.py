"""
PROBE R18 -- REGIMES AND ROUGHNESS. Reuses R7/R10. Settles the rate-conflation (A) + the roughness/branch structure.

A: exact ratio table d_k, Lambda_r, b_r (=bulk/S_r) -- is the exact-regime |ratio| near sqrt(1/3)=0.5774?
   Separately-labelled column: the float-era rho~0.984 fit (phase_routeB_prime_eps_fit.json) -- a DIFFERENT object
   (2-mode envelope decay of the SIGNED d_k over k=7..13, sign-change near k=9), not a term ratio.
B: roughness of delta -- |nu_hat(chi_k)|^2 binned by v3(k) (R8 W_j strata), EXACT via
   Sum_{v3(k)=j} |nu_hat|^2 = 3^j Sum_u g_r(u) c_{3^{r-j}}(u mod 3^{r-j})  (Ramanujan). Low-order fraction shrink?
C: branch factorization T = U+ D+ + U- D- : v-parity split of the R16-A-certified renewal (ord=2*3^{k-1} EVEN =>
   parity constant within each residue class mod ord). Check sum reproduces mu_k EXACT; DC split 1/3,2/3;
   circle-avg weights Sum p_v^2 = 1/15 (even) + 4/15 (odd) = 1/3.
D: max-coefficient profile -- max_{3|/xi}|mu_hat_r(xi)| (ADDITIVE) vs typical sqrt(S_r/(2*3^{r-1})). Spike or not?
E: R85 rung-1 feasibility -- one line.
"""
import os, sys, math, cmath, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import probe_engine_R7 as R7
import probe_charledger_R10 as R10

# --- exact epsilon table (d_k = S_k - 7/15), certified through k=8 ---
_hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                    'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in _hist.items()}
S = {k: F(7, 15) + EPS[k] for k in EPS}
KMAX = max(EPS)                              # 8

SQ13 = 1 / math.sqrt(3)                      # 0.57735 -- the sqrt(1/3) reference


def Lunif(r):
    tr = F(3 ** r, 4 ** (3 ** r) - 1) - F(3 ** (r - 1), 4 ** (3 ** (r - 1)) - 1)
    return S[r] * tr / (2 * 3 ** (r - 1))


def mu_hat(mu_r, r, xi):
    N = 3 ** r
    return sum(complex(p) * cmath.exp(2j * math.pi * (xi * a % N) / N) for a, p in mu_r.items())


def build_parity(mu_prev, k):
    """Renewal split by v-parity. ord=2*3^{k-1} is EVEN => within a residue class j mod ord, all v share j's parity.
       Returns (mu_even, mu_odd, mass_even, mass_odd)."""
    M = 3 ** k
    inv2 = pow(2, -1, M)
    ordv, x = 1, 2 % M
    while x != 1:
        x = (x * 2) % M; ordv += 1
    denom = 1 - F(1, 2 ** ordv)
    me, mo = {}, {}
    mass_e = mass_o = F(0)
    for j in range(1, ordv + 1):
        wv = F(1, 2 ** j) / denom
        u = pow(inv2, j, M)
        tgt = me if j % 2 == 0 else mo
        if j % 2 == 0:
            mass_e += wv
        else:
            mass_o += wv
        for a, pa in mu_prev.items():
            rr = (u * (1 + 3 * a)) % M
            tgt[rr] = tgt.get(rr, F(0)) + wv * pa
    return me, mo, mass_e, mass_o


def main():
    print("# PROBE R18 -- REGIMES AND ROUGHNESS.\n")
    mu = {1: R7.mu1()}
    for k in range(2, 8):
        mu[k] = R7.build_mu(mu[k - 1], k)

    # exact Lambda_r = (eps_{r+1}-eps_r)/2, bulk_r = Lambda_r - Lambda_r^unif, b_r = bulk_r/S_r
    Lam = {r: (EPS[r + 1] - EPS[r]) / 2 for r in range(1, KMAX)}         # r=1..7
    Lun = {r: Lunif(r) for r in range(1, KMAX)}
    bulk = {r: Lam[r] - Lun[r] for r in range(1, KMAX)}                  # bulk_1=0
    b = {r: bulk[r] / S[r] for r in range(1, KMAX)}

    # ================= R18-A =================
    print("## R18-A  EXACT RATIO TABLE (measurement, NO fit; settles the rho<->sqrt(1/3) conflation)")
    print(f"   sqrt(1/3) = {SQ13:.5f}  <- the generic sqrt-cancellation term-ratio reference")
    print(f"   {'k':>2} {'d_k=eps_k (exact)':>26} {'float':>13} {'|d_k/d_{k-1}|':>13} {'sign':>5}")
    for k in range(1, KMAX + 1):
        rr = abs(float(EPS[k] / EPS[k - 1])) if k > 1 else float('nan')
        sg = '+' if EPS[k] > 0 else '-'
        print(f"   {k:>2} {str(EPS[k])[:26]:>26} {float(EPS[k]):>+13.6e} {rr:>13.5f} {sg:>5}")
    print(f"   {'r':>2} {'Lambda_r (exact=(e_{r+1}-e_r)/2)':>34} {'float':>13} {'|L_r/L_{r-1}|':>13} {'sign':>5}")
    for r in range(1, KMAX):
        rr = abs(float(Lam[r] / Lam[r - 1])) if r > 1 else float('nan')
        sg = '+' if Lam[r] > 0 else '-'
        print(f"   {r:>2} {str(Lam[r])[:34]:>34} {float(Lam[r]):>+13.6e} {rr:>13.5f} {sg:>5}")
    print(f"   {'r':>2} {'bulk_r=S_r*b_r (exact)':>26} {'b_r=bulk/S_r':>14} {'|bk_r/bk_{r-1}|':>15} {'|b_r/b_{r-1}|':>13} {'sign':>5}")
    for r in range(2, KMAX):
        rbk = abs(float(bulk[r] / bulk[r - 1])) if r > 2 else float('nan')
        rb = abs(float(b[r] / b[r - 1])) if r > 2 else float('nan')
        sg = '+' if bulk[r] > 0 else '-'
        print(f"   {r:>2} {float(bulk[r]):>+26.6e} {float(b[r]):>+14.6e} {rbk:>15.5f} {rb:>13.5f} {sg:>5}")
    # separately-labelled float-era column
    ff = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                     'phase_routeB_prime_eps_fit.json')))
    fit = ff['k=7..13_2mode']
    print("   --- SEPARATELY LABELLED: float-era rho fit (phase_routeB_prime_eps_fit.json, k=7..13, 2-mode) ---")
    print(f"   fit is on the SIGNED d_k over k=7..13 (envelope decay of an OSCILLATION, sign-change ~k=9), NOT a term ratio:")
    print(f"     rho(envelope)={fit['rho']:.4f}  period={fit['period']:.2f}  theta={fit['theta']:.4f}  loss={fit['loss']:.1e}")
    print(f"     d_k actual k=7..13: {[f'{x:+.2e}' for x in fit['actual']]}")
    print(f"   [Q: is exact-regime |ratio| near {SQ13:.3f}? is the float rho the same object? -- table only.]\n")

    # ================= R18-B =================
    print("## R18-B  ROUGHNESS OF delta (EXACT): |nu_hat(chi_k)|^2 binned by v3(k) (R8 W_j strata)")
    print("   Sum_{v3(k)=j} |nu_hat|^2 = Sum_u g_r(u) c_{3^{r-j}}(u mod 3^{r-j})  [Ramanujan, exact rational]")
    print("   j=0 = FINEST (order 3^r chars, roughest);  j=r-1 = COARSEST (order-3 chars).")
    for r in range(4, 8):
        N = 3 ** r
        g = R10.autocorr_dlog(mu[r], r)
        binm = {}
        for j in range(0, r):
            Nj = 3 ** (r - j)
            # Sum_{v3(k)=j, k!=0} zeta^{ku} = c_{3^{r-j}}(u mod 3^{r-j}) : k=3^j m runs m over ONE period (no 3^j).
            binm[j] = sum(g[u] * R7.cram(r - j, u % Nj) for u in range(N))
        tot = sum(binm.values())
        # exact cross-check: Sum_{k!=0}|nu_hat|^2 = N*g[0]-1
        chk = N * g[0] - 1
        fine = binm[0] / tot                                     # v3=0 fraction (finest)
        loworder = sum(binm[j] for j in range(r // 2, r)) / tot  # coarse half (low order)
        print(f"   r={r}: ||delta||^2={float(tot):.6f} (chk N*g0-1={float(chk):.6f} {'OK' if tot==chk else 'DEV'})  "
              f"frac[v3=0,finest]={float(fine):.4f}  frac[coarse half]={float(loworder):.4f}")
        row = "  ".join(f"j={j}:{float(binm[j]/tot):.4f}" for j in range(0, r))
        print(f"        per-stratum frac(v3=j): {row}")
    print("   [Q: does the fine (v3=0) fraction GROW with r (rough field) and coarse/low-order fraction shrink?]\n")

    # ================= R18-C =================
    print("## R18-C  BRANCH FACTORIZATION (forced): T = U+D+ + U-D-  (v-parity split of R16-A renewal)")
    okC = True
    for k in range(2, 6):
        me, mo, mass_e, mass_o = build_parity(mu[k - 1], k)
        # reconstruct: mu_even + mu_odd == mu_k (R7.build_mu)
        rec = dict(me)
        for a, p in mo.items():
            rec[a] = rec.get(a, F(0)) + p
        # compare to canonical mu[k]
        keys = set(rec) | set(mu[k])
        same = all(rec.get(a, F(0)) == mu[k].get(a, F(0)) for a in keys)
        okC = okC and same and (mass_e == F(1, 3)) and (mass_o == F(2, 3))
        print(f"   k={k}: mu_even+mu_odd == mu_k? {same}   DC split: mass(+even)={mass_e} (=1/3? {mass_e==F(1,3)}) "
              f"mass(-odd)={mass_o} (=2/3? {mass_o==F(2,3)})")
    # circle-average branch weights <|D_+-|^2> = Sum_{v even/odd} p_v^2, p_v=2^{-v} (closed form, exact)
    we = F(1, 16) / (1 - F(1, 16))      # Sum_{v even>=2} 4^{-v} = (1/16)/(15/16) = 1/15
    wo = F(1, 4) / (1 - F(1, 16))       # Sum_{v odd>=1}  4^{-v} = (1/4)/(15/16)  = 4/15
    print(f"   circle-avg weights <|D|^2>=Sum p_v^2 (closed form): even={we} (=1/15? {we==F(1,15)}) "
          f"odd={wo} (=4/15? {wo==F(4,15)})  sum={we+wo} (=1/3? {we+wo==F(1,3)})")
    print(f"   => R18-C {'PASS -- branch form CERTIFIED (exact renewal split, weights exact)' if okC else 'FAIL (#39)'}\n")

    # ================= R18-D =================
    print("## R18-D  MAX-COEFFICIENT PROFILE (measurement; the Prop-1.17 quantity)")
    print(f"   {'r':>2} {'max|mu_hat(xi)|_3|/xi':>22} {'typical=sqrt(S_r/2/3^{r-1})':>28} {'ratio max/typ':>14} {'||delta||^2_add':>16}")
    for r in range(2, 8):
        N = 3 ** r
        prim = [xi for xi in range(1, N) if xi % 3 != 0]
        mags = [abs(mu_hat(mu[r], r, xi)) for xi in prim]
        mx = max(mags)
        l2 = sum(m * m for m in mags)                    # additive primitive Plancherel
        typ = math.sqrt(float(S[r]) / (2 * 3 ** (r - 1)))
        print(f"   {r:>2} {mx:>22.6f} {typ:>28.6f} {mx/typ:>14.4f} {l2:>16.6f}")
    print("   [Q: does max stay within a constant factor of typical, or is there a growing spike?]\n")

    # ================= R18-E =================
    print("## R18-E  R85 RUNG-1 FEASIBILITY (one line)")
    print("   Exact-rational: DEAD at r=8 -- Lambda_8=(eps_9-eps_8)/2 needs eps_9 which is FLOAT-ONLY (k>=9 wall),")
    print("   AND mu_8 autocorr is ~19M-Fraction (supp 4374). Float: CHEAP (Bluestein/FFT autocorr O(N log N), N=6561).")
    print("   Verdict: exact gate DEAD (eps_9 wall), float measurement CHEAP -- only a numeric rung-1, no exact term.")


if __name__ == "__main__":
    main()
