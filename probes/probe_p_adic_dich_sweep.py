"""
p-ADIC-DICH-SWEEP (Wilson spec, 2026-09-06) -- establish the 3-way dissociation across MANY primes.
Reuses validated probe_p_adic_dich machinery (build_chain_exact, dlog_table_p, rho_exact, rho_dense,
tower_fold_exact, gammas).  Sweep p in {5,7,11,13,17,19,23,29,31}, q=2.  Report each claim as a COUNT /9:
  (1) ENGINE  Sum_j rho_r(t+j p^{r-1})=rho_{r-1} exact   (2) FORWARD p|k=>gamma>1 + gamma(p)<gamma(2p)
  (3) CONVERSE  is p-nmid k depleted? -> the enriched-non-multiple LEAK SET per p (empty=holds).
Then HUNT the 9 leak-sets for a clean description (<q>-orbits / residues / QR / size ~ (p-1)/ord_p(q)) or
honest null.  Secondary: p=11,13 with q=3,5 (engine/forward q-robust? leak-set change with q?).

Compute discipline [[no_heavy_compute]]: EXACT (Fraction) build only shallow (engine fold, p^r<=~30k keys);
FLOAT chain for the deeper gamma tables (adaptive depth, p^RG<=300k).  vmax small (fold is vmax-independent;
gamma sign robust to 2^-vmax).  KMAX=64 so k=p and k=2p are in range for every prime up to 31.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_p_adic_dich as P

KMAX = 64
VMF = 40          # float-build vmax (tail 2^-40 ~ 9e-13)
VME = 20          # exact-build vmax (fold cancels it exactly)
PRIMES = [5, 7, 11, 13, 17, 19, 23, 29, 31]


def build_chain_float(p, q, RMAX, const=1, shift=None, vmax=VMF):
    if shift is None:
        shift = p
    nus = {0: {const % p: 1.0}}; nu = nus[0]
    for r in range(1, RMAX + 1):
        M = p ** (r + 1); Mr = p ** r; invq = pow(q, -1, Mr); new = {}
        for Xp, wp in nu.items():
            u = invq
            for v in range(1, vmax + 1):
                Xr = (const + shift * ((u * Xp) % Mr)) % M
                new[Xr] = new.get(Xr, 0.0) + (0.5 ** v) * wp
                u = (u * invq) % Mr
        nu = new; nus[r] = nu
    return nus


def rg_for(p):                                   # deepest float-gamma level: p^RG <= 300k, clamped [3,5]
    r = 3
    while p ** (r + 1) <= 300000 and r < 5:
        r += 1
    return r


def engine_pair(p):                              # exact-fold levels: (3,4) if cheap else (2,3)
    return (3, 4) if p ** 4 <= 60000 else (2, 3)


def orbits_mod(p, q):                            # <q>-orbits on (Z/p)^*
    seen = set(); orbs = []
    for c in range(1, p):
        if c in seen:
            continue
        o = []; x = c
        while x not in seen:
            seen.add(x); o.append(x); x = (x * q) % p
        orbs.append(sorted(o))
    return orbs


def qr_mod(p):
    return sorted({(x * x) % p for x in range(1, p)})


def analyze_prime(p, q, verbose=True):
    # white baseline gate
    N3 = p ** 3
    gw = P.gammas(np.ones(N3) / N3, KMAX)
    if max(abs(gw[k] - 1.0) for k in range(1, KMAX + 1)) >= 1e-11:
        return {'p': p, 'q': q, 'white_ok': False}
    RG = rg_for(p)
    # ---- float gamma tables r=2..RG ----
    nusf = build_chain_float(p, q, RG)
    gam = {0: {k: 1.0 for k in range(KMAX + 1)}}
    for r in range(2, RG + 1):
        rr, _ = P.rho_dense(nusf[r], p, r)
        gam[r] = P.gammas(rr, KMAX)
    del nusf
    gtop = gam[RG]; gprev = gam[RG - 1] if (RG - 1) in gam else gam[RG]
    sign = {k: (gtop[k] > 1) for k in range(1, KMAX + 1)}
    stable = {k: ((gtop[k] > 1) == (gprev[k] > 1)) for k in range(1, KMAX + 1)}
    # ---- ENGINE: exact tower fold ----
    ra, rb = engine_pair(p)
    nuse = P.build_chain_exact(p, q, rb, vmax=VME)
    eng_a = P.tower_fold_exact(nuse[ra], nuse[ra - 1], p, ra)
    eng_b = P.tower_fold_exact(nuse[rb], nuse[rb - 1], p, rb)
    del nuse
    engine_ok = eng_a and eng_b
    # ---- FORWARD: all in-range multiples enriched + fiber-depth gamma(p)<gamma(2p) ----
    mults = [k for k in range(1, KMAX + 1) if k % p == 0]
    mult_enr = [k for k in mults if sign[k]]
    fiber_ok = (2 * p <= KMAX) and (gtop[p] < gtop[2 * p])
    forward_ok = (len(mult_enr) == len(mults)) and fiber_ok
    # ---- CONVERSE: leak set (enriched non-multiples) ----
    leak = [k for k in range(1, KMAX + 1) if k % p != 0 and sign[k]]
    converse_holds = (len(leak) == 0)
    res = {
        'p': p, 'q': q, 'white_ok': True, 'RG': RG, 'engine_pair': (ra, rb),
        'engine_ok': engine_ok, 'forward_ok': forward_ok, 'fiber_ok': fiber_ok,
        'mults': mults, 'mult_enr': mult_enr, 'g_p': gtop[p] if p <= KMAX else None,
        'g_2p': gtop[2 * p] if 2 * p <= KMAX else None,
        'leak': leak, 'converse_holds': converse_holds,
        'leak_resid': sorted({k % p for k in leak}),
        'unstable_signs': [k for k in range(1, KMAX + 1) if not stable[k]],
        'gtop': gtop,
    }
    if verbose:
        print(f"  p={p:>2} q={q} RG={RG} eng{(ra,rb)}: ENGINE={'OK' if engine_ok else 'FAIL'} "
              f"FORWARD={'OK' if forward_ok else 'FAIL'}(mult {len(mult_enr)}/{len(mults)}, "
              f"fiber g{p}={gtop[p]:.3f}<g{2*p}={gtop[2*p]:.3f}={fiber_ok}) "
              f"CONVERSE {'HOLDS' if converse_holds else 'FAILS'} leak={leak}")
    return res


def main():
    t0 = time.time()
    print("# p-ADIC-DICH-SWEEP -- 3-way dissociation across 9 primes (q=2)\n")
    print("## PRIMARY SWEEP  p in", PRIMES, " q=2\n")
    R = [analyze_prime(p, 2) for p in PRIMES]
    R = [r for r in R if r.get('white_ok')]
    n = len(R)
    eng = sum(r['engine_ok'] for r in R)
    fwd = sum(r['forward_ok'] for r in R)
    convf = sum(not r['converse_holds'] for r in R)
    print(f"\n## COUNTS over {n} primes")
    print(f"   (1) ENGINE holds        : {eng}/{n}")
    print(f"   (2) FORWARD holds        : {fwd}/{n}   (p|k enriched AND gamma(p)<gamma(2p))")
    print(f"   (3) CONVERSE FAILS       : {convf}/{n}   (nonempty leak set)")
    diss = (eng == n and fwd == n and convf == n)
    print(f"   => DISSOCIATION ESTABLISHED: {diss}  (needs engine {n}/{n} + forward {n}/{n} + converse-fails {n}/{n})")

    print(f"\n## LEAK SETS (enriched non-multiples of p, k=1..{KMAX})")
    for r in R:
        oq = orbits_mod(r['p'], r['q'])
        qr = qr_mod(r['p'])
        # is leak-residue-set a union of <q>-orbits?
        lr = set(r['leak_resid'])
        orb_union = any(set(o) == lr for o in oq) or (lr and lr == set().union(*[set(o) for o in oq if set(o) <= lr]))
        orb_closed = lr == set().union(*([set(o) for o in oq if set(o) & lr] or [set()]))
        print(f"   p={r['p']:>2}: |leak|={len(r['leak']):>2}  leak={r['leak']}")
        print(f"         resid mod p = {r['leak_resid']}  | QR mod p = {qr} | <2>-orbits = {oq}")
        print(f"         leak-resid is union of <2>-orbits: {orb_closed}  | unstable-sign k: {r['unstable_signs']}")

    print(f"\n## CHARACTERIZATION SUMMARY")
    print(f"   {'p':>3} {'ord_p(2)':>8} {'(p-1)/ord':>9} {'|leakres|':>9} {'|leak|':>7} {'leakres=<2>-orbits?':>20}")
    for r in R:
        p = r['p']; oq = orbits_mod(p, 2)
        op = len(next(o for o in oq if 1 in o))            # ord_p(2) = size of orbit of 1... actually orbit of any; use of 1
        # ord_p(2):
        x = 2 % p; o = 1
        while x != 1:
            x = (x * 2) % p; o += 1
        lr = set(r['leak_resid'])
        orb_closed = lr == set().union(*([set(orb) for orb in oq if set(orb) & lr] or [set()]))
        print(f"   {p:>3} {o:>8} {(p-1)//o if (p-1)%o==0 else (p-1)/o:>9} {len(r['leak_resid']):>9} {len(r['leak']):>7} {str(orb_closed):>20}")

    # ---- SECONDARY: q-variation on p=11,13 ----
    print(f"\n## SECONDARY q-variation (p=11,13 ; q=3,5) -- engine/forward q-robust? leak-set vs q?")
    for p in (11, 13):
        for q in (3, 5):
            analyze_prime(p, q)

    print(f"\n# DONE ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
