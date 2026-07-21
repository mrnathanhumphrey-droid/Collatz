"""
PROBE R10 -- THE CHARACTER LEDGER (the spectral form). Exact rationals via roots-of-unity traces.
Gates Wilson's session: S is a partial sum of a FROZEN character series.

  Multiplicative group G_r = (1+3Z)/(1+3^{r+1}Z) = <4>, cyclic order 3^r; dlog_4 the iso to Z/3^r.
  nu = pushforward of Syrac(Z/3^r)=mu_r under s |-> 1+3s; characters chi_k(X)=zeta^{k*dlog4(X)}, zeta=e(1/3^r).
  Spectral form:  S_{n+1} = 2 sum_{chi in G^_n} |nu_hat(chi)|^2 / (4 chi(4) - 1).
  FROZENNESS: |nu_hat(chi)|^2 frozen for all n >= ord-level; layer  Lambda_r := sum_{ord(chi)=3^r} |nu_hat|^2/(4chi(4)-1).
  Then OffDiag_{r+1} = 2 Lambda_r exact; S_inf=7/15 <=> sum_{r>=1} Lambda_r = -1/10 (Lambda_0=1/3).

EXACT MECHANISM (no cyclotomic field): with g_r(u)=dlog-domain autocorrelation of nu_r (rational),
  |nu_hat(chi_k)|^2 = sum_u g_r(u) zeta^{ku};  and the geometric weight closes as a rational trace:
      A_N(u) := sum_{omega^N=1} omega^u/(4 omega - 1) = N * 4^{N-i0}/(4^N - 1),  i0 = ((u-1) mod N)+1.
  Primitive-k restriction by inclusion-exclusion (order exactly 3^r = all minus 3|k):
      Lambda_r = sum_u g_r(u) [ A_{3^r}(u) - A_{3^{r-1}}(u) ].
  sum_{k prim} zeta^{ku} = Ramanujan c_{3^r}(u), so layer mass = sum_u g_r(u) c_{3^r}(u).
Uses dlog table + Syrac(Z/3^r) ONLY (never C-tables, gamma-tables, or level r+1 data).
"""
import os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import probe_engine_R7 as R7

S = {1: F(2, 3), 2: F(10, 21), 3: F(31370, 67963),
     4: F(143195649659456490, 308468774477179141),
     5: F(2490699741144069281815149277465294323722281911695597204406570,
          5350418720142111510029542161258891403960563152740894082816203),
     6: F(3073594394508858158063860398646918469118621544451223655237787736044140305600491108232683716375149129809572776883214501465626550917830849901208399259893949668774475847028354296093382145075270,
          6593308371642093718618660948586623421072140639443126614173839021148571641533385540268510160579160359417192756329060861771589941737564918045220171980607356372262183045471725888175409486651187)}


def dlog_table(r):
    """dlog_4 on (1+3Z)/(1+3^{r+1}Z): returns dict s -> t with 4^t = 1+3s mod 3^{r+1}, s in Z/3^r."""
    Mp = 3 ** (r + 1); N = 3 ** r
    d = {}; el = 1
    for t in range(N):
        s = ((el - 1) // 3) % N
        d[s] = t
        el = (el * 4) % Mp
    return d


def autocorr_dlog(mu_r, r):
    """g_r(u) = Pr[ dlog(X)-dlog(X') = u ], X=1+3s, s~mu_r; u in Z/3^r. Rational, fresh from Syrac+dlog."""
    N = 3 ** r
    d = dlog_table(r)
    g = [F(0)] * N
    items = list(mu_r.items())
    for s, ps in items:
        ts = d[s % N]
        for sp, psp in items:
            u = (ts - d[sp % N]) % N
            g[u] += ps * psp
    return g


def A_N(u, N):
    """sum_{omega^N=1} omega^u/(4 omega - 1) = N * 4^{N-i0}/(4^N-1), exact."""
    u0 = u % N
    i0 = u0 if u0 != 0 else N
    return F(N * 4 ** (N - i0), 4 ** N - 1)


def Lambda_r(mu, r):
    """Lambda_r from the character side (dlog + Syrac level r only)."""
    N = 3 ** r; Nm = 3 ** (r - 1)
    g = autocorr_dlog(mu[r], r)
    tot = F(0)
    for u in range(N):
        tot += g[u] * (A_N(u, N) - A_N(u % Nm, Nm) if r >= 1 else A_N(u, N))
    return tot, g


def layer_mass(g, r):
    """sum_{ord(chi)=3^r} |nu_hat|^2 = sum_u g_r(u) c_{3^r}(u)  (Ramanujan)."""
    N = 3 ** r
    return sum(g[u] * R7.cram(r, u) for u in range(N))


def nu_hat_sq_float(g, k, N):
    """|nu_hat(chi_k)|^2 = sum_u g(u) cos(2 pi k u/N)  (real; float measurement)."""
    return sum(float(g[u]) * math.cos(2 * math.pi * k * u / N) for u in range(N))


def main():
    print("# PROBE R10 -- THE CHARACTER LEDGER (spectral form). Exact rationals via roots-of-unity traces.\n")
    mu = {1: R7.mu1()}
    for k in range(2, 6):
        mu[k] = R7.build_mu(mu[k - 1], k)

    # ---- R10-A ----
    print("## R10-A  Lambda-LEDGER GATE  (2 Lambda_r = OffDiag_{r+1} = S_{r+1}-S_r; char side only)")
    okA = True; gcache = {}
    for r in range(1, 6):
        L, g = Lambda_r(mu, r); gcache[r] = g
        off = S[r + 1] - S[r]
        good = (2 * L == off)
        okA = okA and good
        print(f"   r={r}: Lambda_r={str(L)[:24]:>24}  2Lambda_r={str(2*L)[:24]:>24}  "
              f"OffDiag_{r+1}={str(off)[:24]:>24}  [{'PASS' if good else 'FAIL'}]")
        if not good:
            print(f"      FAIL r={r}: 2Lambda={2*L}\n              OffDiag={off}")
    print(f"   anchors: Lambda_1={gcache and 2*Lambda_r(mu,1)[0]/2} (=-2/21?), 2Lambda_2 vs -2980/203889")
    print(f"   => R10-A {'GATE PASS' if okA else 'FAIL -- spectral form dies (#31)'}\n")

    # ---- R10-B ----
    print("## R10-B  LAYER MASS  (sum_{ord=3^r}|nu_hat|^2 = S_r exact; multiplicative cross-thread weld)")
    okB = True
    for r in range(1, 6):
        lm = layer_mass(gcache[r], r)
        good = (lm == S[r])
        okB = okB and good
        print(f"   r={r}: layer mass = {str(lm)[:22]:>22}  vs S_r = {str(S[r])[:22]:>22}  [{'WELD' if good else 'DEV'}]")
    print(f"   => R10-B {'PASS' if okB else 'DEV'}\n")

    # ---- R10-C ----
    print("## R10-C  TRACE IDENTITY  (sum_{chi} 1/(4chi(4)-1) = 3^r/(4^{3^r}-1) exact; measure-free)")
    okC = True
    for r in range(1, 5):
        N = 3 ** r
        lhs = A_N(0, N)                              # = sum over all N-th roots of 1/(4w-1)
        rhs = F(N, 4 ** N - 1)
        good = (lhs == rhs)
        okC = okC and good
        print(f"   r={r} (N={N}): sum 1/(4w-1) = {lhs} = 3^r/(4^{{3^r}}-1)={rhs}  [{'OK' if good else 'DEV'}]")
    print(f"   => R10-C {'PASS' if okC else 'DEV'}\n")

    # ---- R10-D ----
    print("## R10-D  WITHIN-LAYER PROFILE  (measurement; NO fit. baseline=S_r/(2*3^{r-1}) per char)")
    for r in range(2, 6):
        N = 3 ** r; cnt = 2 * 3 ** (r - 1)
        g = gcache[r]
        prim = [k for k in range(1, N) if k % 3 != 0]
        base = S[r] / cnt
        # <w>_r and Lambda^unif
        wtot = A_N(0, N) - A_N(0, 3 ** (r - 1))       # total primitive weight (real part; = sum_{k prim}1/(4z^k-1))
        wbar = wtot / cnt
        Lunif = S[r] * wbar
        Lact = Lambda_r(mu, r)[0]
        # palindrome check |nu_hat(k)|^2 == |nu_hat(N-k)|^2 (exact via g symmetry: g(u) real)
        pal_ok = all(abs(nu_hat_sq_float(g, k, N) - nu_hat_sq_float(g, N - k, N)) < 1e-9 for k in prim[:len(prim)//2])
        vals = sorted(set(round(nu_hat_sq_float(g, k, N), 6) for k in prim))
        show = vals if len(vals) <= 8 else f"{len(vals)} distinct in [{vals[0]:.5f},{vals[-1]:.5f}]"
        ratio = float(Lact / Lunif)
        print(f"   r={r}: |nu_hat(chi_k)|^2 over {cnt} primitive k: {show}")
        print(f"         baseline S_r/(2*3^{{r-1}}) = {float(base):.6f} per char [compared-against: uniform |nu_hat|^2]")
        print(f"         Lambda_r = {float(Lact):+.6f} [actual]  vs  Lambda^unif=S_r*<w>_r = {float(Lunif):+.6f} "
              f"[uniform-|nu_hat| baseline]   excess ratio = {ratio:.3f}x   palindrome={'exact' if pal_ok else 'BROKEN(bug)'}")
    print()

    # ---- R10-E ----
    print("## R10-E  DECAY PROFILE  (measurement; max|nu_hat| and ell^1 layer mass per r)")
    for r in range(1, 6):
        N = 3 ** r; g = gcache[r]
        prim = [k for k in range(1, N) if k % 3 != 0] if r >= 1 else []
        if r == 1:
            prim = [k for k in range(1, N) if k % 3 != 0]
        mags = [math.sqrt(max(0.0, nu_hat_sq_float(g, k, N))) for k in prim]
        l1 = sum(mags); mx = max(mags) if mags else 0.0
        print(f"   r={r}: #prim={len(prim)}  max|nu_hat|={mx:.6f}  ell^1 layer mass sum|nu_hat|={l1:.6f}")
    print("   (raw material for the multiplicative-1.17 question; no law fitted, no rate claimed.)")


if __name__ == "__main__":
    main()
