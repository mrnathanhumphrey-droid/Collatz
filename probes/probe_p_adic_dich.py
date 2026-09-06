"""
PROBE p-ADIC-DICH (Wilson spec, 2026-09-06) -- does the enrichment dichotomy gamma_inf(k)>1 <=> p|k
lift off p=3?  CHANNELDICH proved (p=3): gamma_inf(k) ENRICHED (>1 vs white) iff 3|k, DEPLETED off it.
Engine = fiber-tower telescoping Sum_j rho_r(t+j p^{r-1}) = rho_{r-1}(t).  Wilson argues it survives changing
the prime.  THIS tests p=5,7 for the p-adic affine renewal  a_r = q^{-v}(1 + p a_{r-1}) mod p^r, v~Geom(1/2).

FAITHFULNESS (matched to the p=3 machinery, copied structure, generalized p,q):
 * chain build = build_nu_exact (probe_ratio2) with 3->p, 2->q.  Tracks X=1+p*a; a_r=q^{-v}(1+p a_{r-1}) mod p^r.
 * dlog coord = dlog_table (probe_charledger_R10): "dlog_4 on 1+3Z".  Generator 4 = 1+3 = 1+p (principal-unit
   generator; at p=3 it coincides with q^2=4, a p=3 accident).  GENERALIZE generator -> (1+p).  (q^2=4 is NOT a
   principal unit mod 5 or 7, so 1+p is the only faithful generalization of the coordinate.)
 * gamma_r(k) = p^r C_r(k), C_r(k)=<rho_r, shift_k rho_r>.  White baseline uniform rho => C=1/p^r => gamma=1.

GATES (kill-first): (0) white baseline gamma==1 ; (1) p=3 reproduction: my builder == build_nu_exact AND the
3|k dichotomy reproduced.  Only then p=5 ; if p=5 FAILS CLEAN, STOP before p=7.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
from probe_ratio2 import build_nu_exact          # p=3 reference exact builder (for the faithfulness assert)

VMAX = 48                                          # tail 2^-48 ~ 3.6e-15 (normalization cancels it in the fold)
KMAX = 14


# ------------------------------------------------------------------ generalized exact chain build
def build_chain_exact(p, q, RMAX, const=1, shift=None, vmax=VMAX):
    """a_r = q^{-v} ( const + shift * a_{r-1} ) mod p^r,  v~Geom(1/2).  Track X = const + shift*a on Z/p^{r+1}.
       Defaults const=1, shift=p  ==  build_nu_exact generalized (the primary p-adic affine renewal)."""
    if shift is None:
        shift = p
    nus = {0: {const % (p): F(1)}}                  # a_0 = 0  -> X_0 = const
    nu = nus[0]
    for r in range(1, RMAX + 1):
        M = p ** (r + 1); Mr = p ** r; invq = pow(q, -1, Mr)
        new = {}
        for Xp, wp in nu.items():
            u = invq
            for v in range(1, vmax + 1):
                Xr = (const + shift * ((u * Xp) % Mr)) % M
                new[Xr] = new.get(Xr, F(0)) + F(1, 2 ** v) * wp
                u = (u * invq) % Mr
        nu = new; nus[r] = nu
    return nus


# ------------------------------------------------------------------ generalized dlog coordinate (generator 1+p)
def dlog_table_p(p, r):
    """dlog_{1+p} on the principal units (1+pZ)/(1+p^{r+1}Z): dict a -> t with (1+p)^t = 1+p*a mod p^{r+1}."""
    Mp = p ** (r + 1); N = p ** r; g = 1 + p
    d = {}; el = 1
    for t in range(N):
        d[((el - 1) // p) % N] = t
        el = (el * g) % Mp
    return d


def rho_exact(nu, p, r, const=1):
    """dlog-pushforward of the chain measure, EXACT rationals, normalized. a = (X-const)//p."""
    N = p ** r; d = dlog_table_p(p, r)
    rho = {}
    for X, w in nu.items():
        a = (X - const) // p % N
        s = d[a]
        rho[s] = rho.get(s, F(0)) + w
    tot = sum(rho.values())
    return {s: w / tot for s, w in rho.items()}, N


def rho_dense(nu, p, r, const=1):
    N = p ** r; d = dlog_table_p(p, r)
    rr = np.zeros(N)
    for X, w in nu.items():
        a = (X - const) // p % N
        rr[d[a]] += float(w)
    return rr / rr.sum(), N


def gammas(rr, kmax=KMAX):
    N = len(rr)
    return {k: N * float(np.dot(rr, np.roll(rr, -k))) for k in range(kmax + 1)}


def tower_fold_exact(nu_r, nu_rm1, p, r, const=1):
    """dlog-coordinate tower check: fold rho_r by (s mod p^{r-1}) and compare EXACTLY to rho_{r-1}."""
    rho_r, _ = rho_exact(nu_r, p, r, const)
    rho_rm1, _ = rho_exact(nu_rm1, p, r - 1, const)
    m = p ** (r - 1)
    fold = {}
    for s, w in rho_r.items():
        fold[s % m] = fold.get(s % m, F(0)) + w
    keys = set(fold) | set(rho_rm1)
    return all(fold.get(t, F(0)) == rho_rm1.get(t, F(0)) for t in keys)


def marginal_consistent(nu_r, nu_rm1, p, r, const=1):
    """RAW (coordinate-free) tower check: reduce a-values of nu_r mod p^{r-1}, compare law to nu_{r-1}'s a-values.
       Used for the hypothesis-boundary variants where the dlog coordinate meaning changes."""
    def amarg(nu, rr, mod):
        out = {}
        for X, w in nu.items():
            a = (X - const) // p % (p ** rr) % mod
            out[a] = out.get(a, F(0)) + w
        tot = sum(out.values())
        return {a: w / tot for a, w in out.items()}
    lo = amarg(nu_r, r, p ** (r - 1))
    hi = amarg(nu_rm1, r - 1, p ** (r - 1))
    keys = set(lo) | set(hi)
    return all(lo.get(t, F(0)) == hi.get(t, F(0)) for t in keys)


def dichotomy_table(p, q, RMAX, const=1, shift=None, label=""):
    """build exact chain, gamma_r(k) for r=2..RMAX, extrapolate gamma_inf(k), print 2x2.  Returns (verdict, table)."""
    nus = build_chain_exact(p, q, RMAX, const=const, shift=shift)
    gam = {}
    for r in range(1, RMAX + 1):
        rr, _ = rho_dense(nus[r], p, r, const)
        gam[r] = gammas(rr)
    gam[0] = {k: 1.0 for k in range(KMAX + 1)}
    A = {r: {k: gam[r][k] - gam[r - 1][k] for k in range(KMAX + 1)} for r in range(1, RMAX + 1)}

    def ginf(k):
        # geometric-tail extrapolation if the last A-ratios are stable in (-1,1); else deepest gamma_r
        ratios = [A[r][k] / A[r - 1][k] for r in range(RMAX - 1, RMAX + 1)
                  if r - 1 >= 2 and abs(A[r - 1][k]) > 1e-15]
        if ratios:
            rho = float(np.median(ratios))
            if 0 < rho < 1:
                return gam[RMAX][k] + A[RMAX][k] * rho / (1 - rho), rho
        return gam[RMAX][k], float('nan')

    print(f"\n## {label}  p={p}, q={q}, exact r=2..{RMAX}   (dichotomy: gamma_inf(k)>1 <=> {p}|k ?)")
    print(f"   {'k':>2} {'p|k':>4} " + " ".join(f"g_{r}".rjust(9) for r in range(2, RMAX + 1))
          + f" {'g_inf':>9} {'rho':>6} {'sign':>9} {'pred':>9} {'':>6}")
    tab = {}
    ok = True
    for k in range(1, KMAX + 1):
        gi, rho = ginf(k)
        pk = (k % p == 0)
        enr = "ENRICHED" if gi > 1 else "depleted"
        pred = "enriched" if pk else "depleted"
        match = (gi > 1) == pk
        ok = ok and match
        tab[k] = (pk, gi, gi > 1, match)
        row = " ".join(f"{gam[r][k]:9.5f}" for r in range(2, RMAX + 1))
        print(f"   {k:>2} {'YES' if pk else '  .':>4} {row} {gi:9.5f} {rho:6.3f} {enr:>9} {pred:>9} "
              f"{'' if match else 'MISS':>6}")
    # 2x2
    n_pk_enr = sum(1 for k in tab if tab[k][0] and tab[k][2])
    n_pk_dep = sum(1 for k in tab if tab[k][0] and not tab[k][2])
    n_np_enr = sum(1 for k in tab if not tab[k][0] and tab[k][2])
    n_np_dep = sum(1 for k in tab if not tab[k][0] and not tab[k][2])
    print(f"   2x2 [rows p|k / p-nmid ; cols enriched / depleted]:")
    print(f"        p|k   : enriched {n_pk_enr:2d}   depleted {n_pk_dep:2d}")
    print(f"        p-nmid: enriched {n_np_enr:2d}   depleted {n_np_dep:2d}")
    clean_lift = (n_pk_dep == 0 and n_np_enr == 0 and n_pk_enr > 0)
    clean_fail = (n_pk_dep > 0 or n_np_enr > 0)
    verdict = "LIFTS" if (ok and clean_lift) else ("FAILS-CLEAN" if clean_fail and not clean_lift else "PARTIAL")
    print(f"   => {label} VERDICT: {verdict}   (all-k match={ok})")
    return verdict, tab, gam, nus


def main():
    t0 = time.time()
    print("# PROBE p-ADIC-DICH -- does gamma_inf(k)>1 <=> p|k lift off p=3?\n")

    # ================= GATE 0: white baseline =================
    print("## GATE 0 -- white baseline: uniform rho => gamma(k)==1 exactly (normalization sanity)")
    okw = True
    for p, r in [(3, 4), (5, 3), (7, 3)]:
        N = p ** r
        g = gammas(np.ones(N) / N)
        bad = max(abs(g[k] - 1.0) for k in range(1, KMAX + 1))
        okw = okw and bad < 1e-12
        print(f"   p={p} r={r}: max|gamma_white(k)-1| = {bad:.2e}  {'OK' if bad<1e-12 else 'FAIL'}")
    if not okw:
        print("   *** GATE 0 FAIL: normalization/coordinate wrong. STOP. ***"); return
    print("   GATE 0 PASS.\n")

    # ================= GATE 1: p=3 faithfulness + reproduction =================
    print("## GATE 1 -- p=3 faithfulness (my builder == build_nu_exact) + reproduce 3|k dichotomy")
    mine = build_chain_exact(3, 2, 5)
    ref = build_nu_exact(5, vmax=VMAX)             # match vmax so the truncated weights are comparable
    ident = all(mine[r] == ref[r] for r in range(0, 6))
    print(f"   build_chain_exact(3,2) == build_nu_exact (r=0..5), EXACT dict equality: {ident}")
    if not ident:
        # localize
        for r in range(6):
            if mine[r] != ref[r]:
                print(f"     first mismatch at r={r}: |mine|={len(mine[r])} |ref|={len(ref[r])}")
                break
        print("   *** GATE 1 FAIL: builder not faithful. STOP. ***"); return
    v3, tab3, gam3, nus3 = dichotomy_table(3, 2, 5, label="GATE1 p=3 repro")
    if v3 != "LIFTS":
        print("   *** GATE 1 FAIL: p=3 dichotomy not reproduced. coordinate/gamma mismatch. STOP. ***"); return
    print(f"   GATE 1 PASS: p=3 faithful + dichotomy reproduced ({time.time()-t0:.1f}s).\n")

    # ================= PRIMARY: p=5 (kill-first) =================
    v5, tab5, gam5, nus5 = dichotomy_table(5, 2, 5, label="PRIMARY p=5")
    print(f"   ({time.time()-t0:.1f}s)")
    if v5 == "FAILS-CLEAN":
        print("\n## KILL-FIRST: p=5 FAILS CLEAN -> dichotomy is 3-specific. STOPPING before p=7 (per spec).")
        _secondary(5, nus5, gam5, tab5, only_tower=True, t0=t0)
        print(f"\n# DONE ({time.time()-t0:.1f}s)  VERDICT p=5 = FAILS-CLEAN")
        return

    # ================= PRIMARY: p=7 =================
    v7, tab7, gam7, nus7 = dichotomy_table(7, 2, 5, label="PRIMARY p=7")
    print(f"   ({time.time()-t0:.1f}s)")

    # ================= SECONDARY (only if primary lifts) =================
    lifts = (v5 == "LIFTS" and v7 == "LIFTS")
    print(f"\n## PRIMARY VERDICTS: p=5 = {v5} ; p=7 = {v7}  ->  {'BOTH LIFT' if lifts else 'NOT both LIFT'}")
    _secondary(5, nus5, gam5, tab5, only_tower=not lifts, t0=t0)
    _secondary(7, nus7, gam7, tab7, only_tower=not lifts, t0=t0)

    print(f"\n# DONE ({time.time()-t0:.1f}s)  VERDICTS: p=5={v5}, p=7={v7}")


def _secondary(p, nus, gam, tab, only_tower, t0):
    print(f"\n## SECONDARY p={p}")
    # (a) TOWER CHECK direct, exact Fractions, r=3,4
    print(f"   (a) TOWER (fiber telescoping) EXACT check, dlog coords, r=3,4:")
    for r in (3, 4):
        ok = tower_fold_exact(nus[r], nus[r - 1], p, r)
        print(f"       Sum_j rho_{r}(t+j*{p}^{r-1}) == rho_{r-1}(t) exactly: {ok}")
    if only_tower:
        return
    # (c) enrichment magnitude at p|k multiples (does it grow with fiber depth like p=3 did?)
    print(f"   (c) enrichment magnitude at p|k (grows at deeper fiber multiples?):")
    for k in range(1, KMAX + 1):
        if k % p == 0:
            print(f"       gamma_inf({k}) = {tab[k][1]:.5f}  ({'ENRICHED' if tab[k][2] else 'depleted'})")
    # (b) HYPOTHESIS-BOUNDARY: is the p-divisible affine shift load-bearing?
    print(f"   (b) HYPOTHESIS-BOUNDARY (marginal-consistency / tower survival), r=3,4:")
    #   (b1) unit affine constant u!=1 (keep p-divisible shift)
    u = 2 if p != 2 else 3
    nb1 = build_chain_exact(p, 2, 4, const=u, shift=p)
    for r in (3, 4):
        ok = marginal_consistent(nb1[r], nb1[r - 1], p, r, const=u)
        print(f"       [const u={u}, shift={p}]  nu_{r} mod {p}^{r-1} == nu_{r-1}: {ok}")
    #   (b2) NON-p-divisible shift (full-range), keep const=1  -- predicted to BREAK
    nb2 = build_chain_exact(p, 2, 4, const=1, shift=1)
    for r in (3, 4):
        ok = marginal_consistent(nb2[r], nb2[r - 1], p, r, const=1)
        print(f"       [const 1, shift=1 (NON-p-div)]  nu_{r} mod {p}^{r-1} == nu_{r-1}: {ok}  (expect False=BREAK)")


if __name__ == "__main__":
    main()
