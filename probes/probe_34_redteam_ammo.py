"""
PROBE 34 -- red-team ammo for L3 (the 3 falsifiable checks in PHASE3_L3_REDTEAM.md sec.7).

CHECK 1 (Attack 2 -- does the gap survive L->inf?): extend direct rho_k for q=5 to k=9.
   If r_5 stays ~0.62 (not drifting toward 1), the gap is robust in the true (L->inf) limit.
CHECK 2 (Attack 3 -- is the Krylov/tower restriction legit?): track the tower mode's |A|
   across L=1,2 for q=5,7. If |A|~0 at every L (not growing), the tower is genuinely outside
   Krylov(v0) and r_q is clean. (q=3's near-1 mode is r_3 itself, huge |A| -- not a tower.)
CHECK 3 (Attack 1/6 -- the s>=2 regime): tabulate s=v_q(2^d-1), d=ord_q(2) for small primes.
   s>=2 <=> q Wieferich (q^2 | 2^{q-1}-1); smallest is 1093 -- UNCOMPUTABLE. Confirms our 6
   tested primes are ALL s=1, so the s>=2 boundary regime is entirely UNTESTED numerically.

NOT AT STAKE: R10-R33. This is diagnostic ammo, not a new claim.
"""
import numpy as np
from probe_6_conservation_generalize import order_of_two
from probe_30_rq_pin_prony import cross_lean
from probe_25_transfer_operator_Aprime import build_M
from probe_32_offdiag_spectrum import spectrum

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def check1():
    log("## CHECK 1 (Attack 2) -- does r_5 drift toward 1? extend direct rho_k to k=9")
    cr = {}
    for k in range(4, 10):
        vmax = 48 if k >= 9 else 60
        n = (5 - 1) * 5 ** (k - 1)
        try:
            c, nn = cross_lean(5, k, vmax)
        except MemoryError:
            log(f"   k={k}: MemoryError (n={n}) -- SKIP"); break
        cr[k] = c
        log(f"   cross({k}) = {c:.12f}  (n={nn})")
    ck = {k: cr[k] - cr.get(k - 1, 0.0) for k in cr if k - 1 in cr or k == min(cr)}
    ks = sorted(cr)
    rho = []
    for i in range(len(ks) - 1):
        a, b = ks[i], ks[i + 1]
        if b == a + 1 and a in ck and b in ck and abs(ck[a]) > 1e-11:
            rho.append((b, ck[b] / ck[a]))
    log(f"   rho_k (=c_(k+1)/c_k): {[('k%d:%.4f' % (k-1, r)) for k, r in rho]}")
    vals = [r for _, r in rho]
    if vals:
        log(f"   last 3 rho: {['%.4f' % v for v in vals[-3:]]}  -> r_5 settling ~{np.mean(vals[-3:]):.3f}")
        log(f"   VERDICT: {'STABLE (no drift toward 1) -> gap robust' if max(vals[-3:]) < 0.75 else 'DRIFTING -- Attack 2 has teeth'}")
    log("")


def check2():
    log("## CHECK 2 (Attack 3) -- tower mode |A| vs L for q=5,7 (should stay ~0)")
    for q in [5, 7]:
        log(f"   q={q}:")
        for L in [1, 2]:
            M, idx, n = build_M(q, L)
            dense = n <= 400
            modes = spectrum(M, idx[(1, 1, 0)], dense)
            l1 = max(m[0].real for m in modes)
            # tower = modes with |mu| in [0.7,1) that are NOT the Perron (mu~1) and NOT r_q
            # (r_q = largest |mu|<0.7-ish with real amplitude). Report max |A| among |mu| in [0.7,0.999)
            tower = [(abs(lam / l1), abs(A)) for lam, A in modes
                     if 0.7 < abs(lam / l1) < 0.999 and abs(lam / l1 - 1) >= 2e-3]
            tmax = max((a for _, a in tower), default=0.0)
            # r_q mode amplitude (largest |mu|<0.999 with |A|>1e-6 excluding Perron)
            carr = sorted([(abs(lam / l1), abs(A)) for lam, A in modes if abs(A) > 1e-6
                           and abs(lam / l1 - 1) >= 2e-3 and abs(lam / l1) < 0.999],
                          key=lambda t: -t[0])
            rq = carr[0] if carr else (float('nan'), 0.0)
            log(f"      L={L} (dim={n}): tower max|A|={tmax:.2e}   r_q mode |mu|={rq[0]:.4f} |A|={rq[1]:.3e}")
        log("")
    log("   VERDICT: tower |A| ~0 at every L (not growing) => tower outside Krylov(v0),")
    log("            r_q clean. If tower |A| grows with L, Attack 3 has teeth.")
    log("")


def check3():
    log("## CHECK 3 (Attack 1/6) -- s=v_q(2^d-1) spectrum: where does s>=2 first occur?")
    log(f"   {'q':>5} {'d=ord':>6} {'s=v_q(2^d-1)':>13} {'2 prim?':>8} {'tested?':>8}")
    tested = {3, 5, 7, 11, 13, 19, 29}
    primes = [p for p in range(3, 260) if all(p % j for j in range(2, int(p ** 0.5) + 1))]
    smax = 1
    first_s2 = None
    for q in primes:
        d = order_of_two(q)
        val = 2 ** d - 1
        s = 0
        while val % q == 0:
            s += 1; val //= q
        prim = "Y" if d == q - 1 else "n"
        tag = "*tested" if q in tested else ""
        if s >= 2 and first_s2 is None:
            first_s2 = q
        if q in tested or s >= 2 or q < 20:
            log(f"   {q:>5} {d:>6} {s:>13} {prim:>8} {tag:>8}")
    log("")
    log(f"   ALL tested primes {sorted(tested)} have s=1 (confirmed above).")
    log(f"   First prime with s>=2 in [3,260): {first_s2 if first_s2 else 'NONE'}")
    log("   (s>=2 <=> q^2 | 2^d-1 <=> Wieferich; smallest = 1093, then 3511 -- both far beyond")
    log("    compute reach: q=1093 needs phi(1093^k) states, k=2 already ~1.2M, k>=3 ~1e9.)")
    log("   => VERDICT: the s>=2 regime is ENTIRELY UNTESTED. Attack 1/6 lives at Wieferich")
    log("      primes we cannot reach; L3's behavior there must be settled by the MECHANISM,")
    log("      not computation. This is a genuine hole flagged, not closed.")
    log("")


def main():
    log("# PROBE 34 -- L3 red-team ammo (3 falsifiable checks)")
    log("")
    check1()
    check2()
    check3()
    with open("result_34_redteam_ammo_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
