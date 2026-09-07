"""
p-ADIC-DICH-STRAT (Wilson #1, 2026-09-06) -- is FORWARD pass/fail predicted by (p-1)/ord_p(2)
(i.e. by whether 2 is a primitive root / the number of <2>-cosets)?  n=9 showed a 7/9 tendency
(2 primitive root => forward tends to FAIL).  Test it on a bigger stratified sample: primes 5..61.
Reuses probe_p_adic_dich_sweep.analyze_prime (validated).  Adaptive depth, q=2.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import probe_p_adic_dich_sweep as SW

PRIMES = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]


def ord_p(a, p):
    x = a % p; o = 1
    while x != 1:
        x = (x * a) % p; o += 1
    return o


def main():
    t0 = time.time()
    print("# p-ADIC-DICH-STRAT -- forward pass/fail vs (p-1)/ord_p(2) over primes 5..61 (q=2)\n")
    rows = []
    for p in PRIMES:
        r = SW.analyze_prime(p, 2, verbose=False)
        if not r.get('white_ok'):
            print(f"  p={p}: white gate FAIL, skipped"); continue
        o = ord_p(2, p)
        cosets = (p - 1) // o
        primroot = (cosets == 1)
        gp = r['gtop'][p] if p <= SW.KMAX else None
        g2p = r['gtop'][2 * p] if 2 * p <= SW.KMAX else None
        kp_enr = (gp is not None and gp > 1)
        rows.append(dict(p=p, RG=r['RG'], o=o, cosets=cosets, primroot=primroot,
                         forward=r['forward_ok'], kp_enr=kp_enr, gp=gp, g2p=g2p,
                         nmult=len(r['mults']), nmult_enr=len(r['mult_enr']), leak_n=len(r['leak'])))
        print(f"  p={p:>2} RG={r['RG']} ord2={o:>2} cosets={cosets} primroot2={str(primroot):>5}  "
              f"FORWARD={'OK ' if r['forward_ok'] else 'FAIL'} (k=p {'enr' if kp_enr else 'DEPL'} g={gp:.3f}, "
              f"mult {len(r['mult_enr'])}/{len(r['mults'])})")

    print(f"\n## CONTINGENCY: forward-holds vs (2 is primitive root)")
    pr_ok = sum(1 for x in rows if x['primroot'] and x['forward'])
    pr_no = sum(1 for x in rows if x['primroot'] and not x['forward'])
    np_ok = sum(1 for x in rows if not x['primroot'] and x['forward'])
    np_no = sum(1 for x in rows if not x['primroot'] and not x['forward'])
    npr = pr_ok + pr_no; nnp = np_ok + np_no
    print(f"   {'':>18}{'forward OK':>11}{'forward FAIL':>13}")
    print(f"   2 primitive root {pr_ok:>11}{pr_no:>13}   (n={npr})")
    print(f"   2 NOT prim root  {np_ok:>11}{np_no:>13}   (n={nnp})")
    print(f"   tendency 'primroot => forward fails': {pr_no}/{npr} ; 'non-primroot => forward ok': {np_ok}/{nnp}")

    print(f"\n## same, using the crisper single indicator  k=p ENRICHED (gamma(p)>1)")
    a = sum(1 for x in rows if x['primroot'] and x['kp_enr'])
    b = sum(1 for x in rows if x['primroot'] and not x['kp_enr'])
    c = sum(1 for x in rows if not x['primroot'] and x['kp_enr'])
    d = sum(1 for x in rows if not x['primroot'] and not x['kp_enr'])
    print(f"   {'':>18}{'k=p enr':>9}{'k=p depl':>10}")
    print(f"   2 primitive root {a:>9}{b:>10}")
    print(f"   2 NOT prim root  {c:>9}{d:>10}")

    print(f"\n## exceptions to 'non-primroot => forward OK' and 'primroot => forward FAIL'")
    for x in rows:
        exc = (x['primroot'] and x['forward']) or (not x['primroot'] and not x['forward'])
        if exc:
            print(f"   p={x['p']:>2} primroot2={x['primroot']} forward={x['forward']} "
                  f"(k=p g={x['gp']:.3f}, cosets={x['cosets']}) <-- exception")

    print(f"\n## does forward track cosets=(p-1)/ord more finely than the binary primroot?")
    from collections import defaultdict
    byc = defaultdict(lambda: [0, 0])
    for x in rows:
        byc[x['cosets']][0 if x['forward'] else 1] += 1
    for cst in sorted(byc):
        ok, no = byc[cst]
        print(f"   cosets={cst:>2}: forward OK {ok}, FAIL {no}")

    print(f"\n# DONE ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
