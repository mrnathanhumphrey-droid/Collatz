"""
p-ADIC-DICH-ORBIT (Wilson, 2026-09-06) -- is FORWARD pass/fail predicted by the <2>-ORBIT SIZE ord_p(2)?
From the strat sweep: forward OK iff ord_p(2) small (<=8), FAIL iff large (>=10), 0 exceptions in 16 primes.
Wilson's "fails when a specific orbit is filled" = the <2>-orbit is LARGE.  This VERIFIES it and runs the
DECOUPLING test that separates ord from prime-size:
  p=127 ord=7  (LARGE prime, tiny orbit)  -> predict forward OK   (the killer test)
  p= 73 ord=9  (threshold)                -> discriminator
  p= 71 ord=35 (large orbit)              -> predict forward FAIL  (matched shallow depth vs 127)
Reuses probe_p_adic_dich_sweep machinery.  Per-prime KMAX=2p+2 so k=p,2p are in range for every prime.
"""
import os, sys, time, gc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import probe_p_adic_dich as P
from probe_p_adic_dich_sweep import build_chain_float

BASE = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]
DECOUPLE = [71, 73, 127]          # large primes: 71 ord35, 73 ord9, 127 ord7


def ord_p(a, p):
    x = a % p; o = 1
    while x != 1:
        x = (x * a) % p; o += 1
    return o


def rg_for(p):                    # deepest float level p^RG<=300k, clamp [2,5]
    r = 2
    while p ** (r + 1) <= 300000 and r < 5:
        r += 1
    return r


def gamma_at(rr, ks):
    N = len(rr)
    return {k: N * float(np.dot(rr, np.roll(rr, -k))) for k in ks}


def forward_core(p, q=2, rforce=None):
    RG = rforce if rforce else rg_for(p)
    nus = build_chain_float(p, q, RG)
    rr, _ = P.rho_dense(nus[RG], p, RG)
    del nus
    ks = [p, 2 * p] + [j * p for j in (3, 4, 5) if j * p <= 2 * p + 2 * p]  # p,2p,(3p..5p if small)
    ks = sorted(set([k for k in ks if k < len(rr)]))
    g = gamma_at(rr, ks)
    del rr; gc.collect()
    gp, g2p = g[p], g[2 * p]
    core = (gp > 1) and (g2p > 1) and (gp < g2p)
    mult_j = [j for j in (1, 2, 3, 4, 5) if j * p in g]
    mult_enr = [j for j in mult_j if g[j * p] > 1]
    return dict(p=p, RG=RG, gp=gp, g2p=g2p, core=core,
                mult_j=mult_j, mult_enr=mult_enr, g=g)


def main():
    t0 = time.time()
    print("# p-ADIC-DICH-ORBIT -- forward pass/fail vs <2>-orbit size ord_p(2)?  (q=2)\n")
    rows = []
    for p in BASE + DECOUPLE:
        o = ord_p(2, p)
        r = forward_core(p)
        rows.append(dict(ord=o, **r, decouple=(p in DECOUPLE)))
    rows.sort(key=lambda x: x['ord'])

    print(f"## sorted by ord_p(2)   (forward_core = k=p enr AND k=2p enr AND gamma(p)<gamma(2p))")
    print(f"   {'p':>4} {'ord':>4} {'RG':>3} {'g(p)':>9} {'g(2p)':>9} {'forward':>8} {'mult_enr/j':>12} {'':>6}")
    for x in rows:
        tag = " <-DECOUPLE" if x['decouple'] else ""
        print(f"   {x['p']:>4} {x['ord']:>4} {x['RG']:>3} {x['gp']:>9.3f} {x['g2p']:>9.3f} "
              f"{'OK' if x['core'] else 'FAIL':>8} {str(len(x['mult_enr']))+'/'+str(len(x['mult_j'])):>12}{tag}")

    # threshold search
    print(f"\n## THRESHOLD: is there T with  forward-OK <=> ord < T ?")
    ok_ords = sorted(x['ord'] for x in rows if x['core'])
    no_ords = sorted(x['ord'] for x in rows if not x['core'])
    print(f"   ord of forward-OK primes  : {ok_ords}")
    print(f"   ord of forward-FAIL primes: {no_ords}")
    if ok_ords and no_ords:
        clean = max(ok_ords) < min(no_ords)
        print(f"   clean separation: {clean}  (max OK ord {max(ok_ords)} {'<' if clean else '>='} min FAIL ord {min(no_ords)})")
        if clean:
            print(f"   => THRESHOLD in ({max(ok_ords)}, {min(no_ords)}]  -- forward OK IFF ord_p(2) <= {max(ok_ords)}")

    # decoupling verdict
    print(f"\n## DECOUPLING VERDICT (ord vs prime-size)")
    for x in rows:
        if x['decouple']:
            pred = "OK (small orbit)" if x['ord'] <= 9 else "FAIL (large orbit)"
            got = "OK" if x['core'] else "FAIL"
            hit = "CONFIRMS ord-not-size" if ((x['ord'] <= 9) == x['core']) else "BREAKS the ord rule"
            print(f"   p={x['p']} (LARGE, ord={x['ord']}, RG={x['RG']}): predict {pred:>18} -> got {got:>4}  [{hit}]")

    # ---- the CLEAN decoupling: matched size (~100) AND matched depth r=3, varying ord ----
    print(f"\n## CLEAN DECOUPLING at MATCHED r=3 (size all ~100, only ord varies): 127(ord7) vs 73(ord9) vs 71(ord35)")
    for p in (127, 73, 71):
        o = ord_p(2, p)
        r = forward_core(p, rforce=3)
        pred = "OK (small orbit)" if o <= 9 else "FAIL (large orbit)"
        print(f"   p={p} ord={o:>2} r=3: g(p)={r['gp']:.3f} g(2p)={r['g2p']:.3f} -> "
              f"{'OK' if r['core'] else 'FAIL':>4}  (mult_enr {len(r['mult_enr'])}/{len(r['mult_j'])})  predict {pred}")
    print(f"\n# DONE ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
