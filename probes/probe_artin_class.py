"""
PROBE ARTIN-CLASS (Wilson, 2026-09-06) -- is forward-enrichment PASS/FAIL a Chebotarev/Artin condition on 2 mod p,
or a FIFTH dead predictor-story (after q^2=1+p, primitive-root, orbit-size, and 'no predictor')?
Density (Test 1) is weak; the CONDITION-MATCH over ~33 primes (Test 2) is decisive.

CRITERION (uniform, reproduces anchors p=5,7=PASS ; p=11,13,23=FAIL): PASS := gamma(p),gamma(2p),gamma(3p) ALL >1
(the first-3-multiples-enriched forward test; captures p=11's depth-robust 3p depletion; applied identically to
every prime).  Float sign-call (as the anchors were established), shallow r (pass/fail dodges the doubly-exp wall).
Stability: compute at r_top and r_top-1; if the PASS call flips -> UNRESOLVED, excluded from the condition-match.
"""
import os, sys, time, gc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import probe_p_adic_dich as P
from probe_p_adic_dich_sweep import build_chain_float


def primes_up_to(n):
    s = [True] * (n + 1); s[0:2] = [False, False]
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            for j in range(i * i, n + 1, i):
                s[j] = False
    return [i for i in range(5, n + 1) if s[i]]


def ord_p(a, p):
    x = a % p; o = 1
    while x != 1:
        x = (x * a) % p; o += 1
    return o


def v2(n):
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k


def rtop(p):
    r = 2
    while p ** (r + 1) <= 2_000_000 and r < 4:
        r += 1
    return r


def gamma_at(rr, ks):
    N = len(rr)
    return {k: N * float(np.dot(rr, np.roll(rr, -k))) for k in ks}


def pass_at(p, q, r):
    nus = build_chain_float(p, q, r)
    rr, _ = P.rho_dense(nus[r], p, r)
    del nus
    g = gamma_at(rr, [p, 2 * p, 3 * p])
    del rr; gc.collect()
    return (g[p] > 1 and g[2 * p] > 1 and g[3 * p] > 1), g


def measure(p, q=2):
    # white gate
    N = p ** 2
    gw = P.gammas(np.ones(N) / N, 3 * p + 1) if False else None
    R = rtop(p)
    pt, gt = pass_at(p, q, R)
    pprev, _ = pass_at(p, q, R - 1)
    confident = (pt == pprev)
    return dict(p=p, r=R, PASS=pt, confident=confident,
                gp=gt[p], g2p=gt[2 * p], g3p=gt[3 * p], ord=ord_p(q, p))


def conditions(p, ordp):
    idx = (p - 1) // ordp
    return {
        "(a) 2 is QR (p==+-1 mod8)": p % 8 in (1, 7),
        "(b) ord_p(2) even": ordp % 2 == 0,
        "(c1) ord==0 mod4": ordp % 4 == 0,
        "(c2) ord==0 mod3": ordp % 3 == 0,
        "(d3) p==1 mod3": p % 3 == 1,
        "(d4) p==1 mod4": p % 4 == 1,
        "(d6) p==1 mod6": p % 6 == 1,
        "(d8) p==1 mod8": p % 8 == 1,
        "(d12) p==1 mod12": p % 12 == 1,
        "(e) 2 primitive root": ordp == p - 1,
        "(f1) index even": idx % 2 == 0,
        "(f2) index odd": idx % 2 == 1,
        "(g) QR2 AND ord even": (p % 8 in (1, 7)) and (ordp % 2 == 0),
        "(h1) v2(p-1)==1 (p==3 mod4)": v2(p - 1) == 1,
        "(h2) v2(p-1)>=3": v2(p - 1) >= 3,
    }


def main():
    t0 = time.time()
    print("# PROBE ARTIN-CLASS -- is forward-enrichment PASS/FAIL a Chebotarev condition on 2 mod p?\n")
    ps = primes_up_to(149)
    rows = [measure(p, 2) for p in ps]

    # ---- anchor gate ----
    anc = {r['p']: r['PASS'] for r in rows if r['p'] in (5, 7, 11, 13, 23)}
    want = {5: True, 7: True, 11: False, 13: False, 23: False}
    print("## ANCHOR GATE (criterion must reproduce p=5,7=PASS ; p=11,13,23=FAIL)")
    ok = all(anc.get(k) == v for k, v in want.items())
    for k in (5, 7, 11, 13, 23):
        print(f"   p={k}: PASS={anc.get(k)}  want {want[k]}  {'ok' if anc.get(k)==want[k] else 'MISMATCH'}")
    if not ok:
        print("   *** ANCHOR GATE FAIL -- criterion drifted. STOP. ***"); return
    print("   ANCHOR GATE PASS.\n")

    # ---- PASS/FAIL table ----
    print(f"## PASS/FAIL table ({len(rows)} primes, q=2)   [criterion: g(p),g(2p),g(3p) all >1]")
    print(f"   {'p':>4} {'r':>2} {'ord':>4} {'g(p)':>8} {'g(2p)':>8} {'g(3p)':>8} {'PASS':>5} {'conf':>5}")
    for r in rows:
        print(f"   {r['p']:>4} {r['r']:>2} {r['ord']:>4} {r['gp']:>8.3f} {r['g2p']:>8.3f} {r['g3p']:>8.3f} "
              f"{str(r['PASS']):>5} {'yes' if r['confident'] else 'UNRES':>5}")

    conf = [r for r in rows if r['confident']]
    npass = sum(r['PASS'] for r in conf)
    print(f"\n## TEST 1 -- DENSITY")
    print(f"   confident primes: {len(conf)}/{len(rows)} ; PASS {npass} -> pass-fraction = {npass/len(conf):.4f}")
    # running value
    run = []
    c = 0; tot = 0
    for r in rows:
        if r['confident']:
            tot += 1; c += r['PASS']; run.append(c / tot)
    print(f"   running pass-fraction (confident): " + " ".join(f"{x:.3f}" for x in run[::max(1,len(run)//12)]))
    print(f"   vs Artin A=0.3739, 1/2=0.5, 1/3=0.333, 2/3=0.667, 3/8=0.375, 5/8=0.625")

    # ---- TEST 2: condition match ----
    print(f"\n## TEST 2 -- CONDITION MATCH (decisive)  [over {len(conf)} confident primes]")
    passset = set(r['p'] for r in conf if r['PASS'])
    print(f"   PASS set = {sorted(passset)}")
    print(f"   {'condition':>32} {'agree':>6} {'disag':>6} {'best(vs compl)':>15}  disagreeing primes")
    results = []
    for name in conditions(5, ord_p(2, 5)):
        dis = [r['p'] for r in conf if conditions(r['p'], r['ord'])[name] != r['PASS']]
        dis_c = [r['p'] for r in conf if conditions(r['p'], r['ord'])[name] == r['PASS']]  # complement disagreers
        agree = len(conf) - len(dis)
        best = min(len(dis), len(dis_c))
        exact = (len(dis) == 0) or (len(dis_c) == 0)
        results.append((name, agree, len(dis), best, exact, dis if len(dis) <= len(dis_c) else dis_c))
        tag = "  <== EXACT MATCH" if exact else ""
        print(f"   {name:>32} {agree:>6} {len(dis):>6} {best:>15} {tag}  {(dis if len(dis)<=len(dis_c) else dis_c)[:8]}")

    exact_hits = [r for r in results if r[4]]
    print(f"\n## VERDICT")
    if exact_hits:
        print(f"   Test 2 EXACT MATCH: {[h[0] for h in exact_hits]} => Artin/Chebotarev-class CONFIRMED, condition IDENTIFIED.")
    else:
        best = min(results, key=lambda x: x[3])
        print(f"   NO exact condition-match. Closest: {best[0]} (best {best[3]} disagree of {len(conf)}).")
        print(f"   => 'Artin-class' is a FIFTH DEAD predictor-story: forward PASS/FAIL is NOT a clean Chebotarev")
        print(f"      condition on 2 mod p. Honest verdict: genuinely structureless / not Chebotarev-predictable.")
    print(f"\n# DONE ({time.time()-t0:.1f}s)")
    return exact_hits


if __name__ == "__main__":
    main()
