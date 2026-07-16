"""
PROBE 31 -- test the PRIMITIVE-ROOT law: for primes where 2 is a primitive root (H=<2>
full), is r_q ~ 3/q, and does r_q*q/3 oscillate around 1?

PRE-REGISTRATION.
------------------------------------------------------------------
LEAD (from the audit): r_q hugs 3/q for primitive-root q {3,5,11,13} but q=7 (2 NOT
primitive, d=3) sits below. And r_q*q/3 looks like it alternates around 1 (5:1.03,
7:0.91,...). This probe EXTENDS the primitive-root set {5,11,13,19,29} to SEE the pattern.

HONEST CAVEAT (stated before running): r_q is only SOLID at q=5,7 (settled rho over many
k). For q>=11 the direct method is oscillation-limited and the matrix pencil is order-
unstable (R30). So q>=11 estimates are FUZZY (report order-2 AND order-3 pencil + rho-tail
as a range). NO closed form is committed; we are LOOKING, not fitting.

GATE: cross(k) matches known exact/prior values at low k (<1e-7). Reuses validated tools.

BUDGET: n*vmax <= 1.5e8. q=19 k<=5, q=29 k<=4. SAID where it stops.

NOT AT STAKE: R10-R30.
"""
import numpy as np
from probe_6_conservation_generalize import order_of_two
from probe_30_rq_pin_prony import cross_lean, matrix_pencil

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def is_primitive_root_2(q):
    return order_of_two(q) == q - 1


def main():
    log("# PROBE 31 -- primitive-root law: r_q vs 3/q for primes where 2 is primitive")
    log("")
    # solid known values (from R27/R30)
    KNOWN = {3: 1.0, 5: 0.62, 7: 0.39}
    KMAX = {11: 6, 13: 6, 19: 5, 29: 4}
    VMAX = {11: 52, 13: 44, 19: 48, 29: 44}
    NCAP = 1.5e8
    rq = dict(KNOWN)
    for q in [11, 13, 19, 29]:
        d = order_of_two(q)
        prim = "primitive" if d == q - 1 else f"NOT primitive (d={d})"
        log(f"## q={q}  ({prim})")
        cr = {}
        for k in range(2, KMAX[q] + 1):
            n = (q - 1) * q ** (k - 1)
            if n * VMAX[q] > NCAP:
                log(f"   k={k}: n*vmax={n*VMAX[q]:.1e} > cap -- STOP (SAID)")
                break
            try:
                c, nn = cross_lean(q, k, VMAX[q])
            except MemoryError:
                log(f"   k={k}: MemoryError -- SKIP (SAID)")
                break
            cr[k] = c
            log(f"   cross({k})={c:.12f} (n={nn})")
        ck = [cr[k] - cr.get(k - 1, 0.0) for k in sorted(cr)]
        rhos = [ck[i + 1] / ck[i] for i in range(len(ck) - 1) if abs(ck[i]) > 1e-13]
        log(f"   rho_k: {['%.4f' % r for r in rhos]}")
        z2 = matrix_pencil(ck, 2)
        z3 = matrix_pencil(ck, 3) if len(ck) >= 4 else []
        d2 = abs(z2[0]) if z2 else float('nan')
        d3 = abs(z3[0]) if z3 else float('nan')
        log(f"   pencil dom |z|: order-2={d2:.4f}  order-3={d3:.4f}   (range = fuzz)")
        # estimate: average of available estimators, flag range
        ests = [x for x in [d2, d3] + (rhos[-2:] if rhos else []) if x == x and 0 < x < 1]
        rq[q] = float(np.median(ests)) if ests else float('nan')
        log(f"   => r_{q} ~ {rq[q]:.3f}  (median of estimators; FUZZY for q>=11)")
        log("")

    log("## PRIMITIVE-ROOT TABLE: r_q vs 3/q  (* = solid; others fuzzy)")
    log(f"   {'q':>4} {'2 prim?':>8} {'r_q':>8} {'3/q':>8} {'r_q*q/3':>9} {'solid?':>7}")
    for q in [3, 5, 7, 11, 13, 19, 29]:
        if q not in rq:
            continue
        prim = "Y" if is_primitive_root_2(q) else "n"
        ratio = rq[q] * q / 3
        solid = "*" if q in (3, 5, 7) else ""
        log(f"   {q:>4} {prim:>8} {rq[q]:>8.3f} {3/q:>8.3f} {ratio:>9.3f} {solid:>7}")
    log("")
    log("## READ: among PRIMITIVE q (Y), does r_q*q/3 stay near 1 (with a wobble)?")
    log("   The non-primitive q=7 (n) is the control -- it should sit off the line.")
    log("   FUZZY q>=11: treat r_q*q/3 as +/-0.1. Pattern only suggestive until pinned.")
    flush()


def flush():
    with open("result_31_primitive_root_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
