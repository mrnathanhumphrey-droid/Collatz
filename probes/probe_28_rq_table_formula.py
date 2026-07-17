"""
PROBE 28 -- ROUTE 1 (option A): extend the r_q table (q=11,13) and hunt the closed form.

PRE-REGISTRATION (written before running).
------------------------------------------------------------------
R27 pinned r_3=1, r_5~0.61, r_7~0.39 (high-k rho_k, v-truncated stationary, gate-validated).
3/q brackets them but isn't exact. Two competing hypotheses for the closed form:
  H_RAT  : r_q = 3/q                      -> r_11=0.2727, r_13=0.2308
  H_EXP  : r_q = a^{q-3} (ln r linear in q). FIT on {5,7}, PREDICT {11,13} (out-of-sample).
           From r_5,r_7: ln r_5=-0.494, ln r_7=-0.942 -> slope -0.224/unit q, a=exp(-0.224)
           -> predicts r_11=exp(-0.224*8)=~0.167, r_13=~0.107.
These differ ~1.6x at q=11 -> r_11 alone discriminates. (Neither is committed to win --
priors 0-for-8; both are laid out to be tested, per feedback_r2_cannot_discriminate.)

METHOD: cross_trunc (R27) -> rho_k = c_{k+1}/c_k. r_q = robust tail estimate (mean of the
last few rho_k, which oscillate around the limit per R26's +/- mode pair).

GATE: cross(k) vs R22/R23 exact for low k (<1e-9). Reuses R27's validated cross_trunc.

DECISION (pre-committed, structural): the winning form must EXTRAPOLATE to q=11 AND q=13
within the oscillation band; a form that fits {5,7} but misses {11,13} is REJECTED. If
NEITHER fits, report "no elementary closed form found; r_q is the operator's subdominant
eigenvalue ratio" -- an honest negative, not a forced fit.

BUDGET: n*vmax<=3e7. q=11 k<=5 (9.4M), q=13 k<=5 (22M), q=5 k<=8, q=7 k<=7. SAID per q.

NOT AT STAKE: R10-R27, R5, R6, R7, R12, THEOREM_C_745.
"""
import numpy as np
from probe_27_high_k_rho_q5 import cross_trunc

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


REF = {(5, 2): 0.038875214298, (5, 3): 0.059623775122, (7, 2): 0.047025340793,
       (7, 3): 0.068042081671, (11, 2): 0.000113752827, (11, 3): 0.000134296401,
       (13, 2): 0.0, (13, 3): 0.0}


def rho_seq(q, kmax, ncap=3e7):
    cr = {}
    for k in range(2, kmax + 1):
        n = (q - 1) * q ** (k - 1)
        if n * 64 > ncap:
            log(f"   q={q} k={k}: n*vmax={n*64:.1e} > cap -- STOP (SAID)")
            break
        c, nn = cross_trunc(q, k)
        cr[k] = c
        g = ""
        if (q, k) in REF and REF[(q, k)] > 0:
            rel = abs(c - REF[(q, k)]) / abs(REF[(q, k)])
            g = f"  gate |rel|={rel:.1e}{' OK' if rel < 1e-8 else ' CHK'}"
        log(f"   cross({k})={c:.12f} (n={nn}){g}")
    ck = {k: cr[k] - cr.get(k - 1, 0.0) for k in cr}
    rho = {}
    for k in sorted(ck):
        if k + 1 in ck and abs(ck[k]) > 1e-11 and abs(ck[k + 1]) > 1e-11:
            rho[k] = ck[k + 1] / ck[k]
    return rho


def main():
    log("# PROBE 28 -- extend r_q table (q=11,13) + hunt closed form")
    log("# H_RAT: r_q=3/q  vs  H_EXP: r_q=a^{q-3} (fit {5,7}, predict {11,13} OOS)")
    log("")

    KMAX = {5: 8, 7: 7, 11: 5, 13: 5}
    rq = {}
    for q in [5, 7, 11, 13]:
        log(f"## q={q}")
        rho = rho_seq(q, KMAX[q])
        ks = sorted(rho)
        log(f"   rho_k: {['%.5f' % rho[k] for k in ks]}")
        # robust tail estimate: mean of last min(3, all) rho
        tail = [rho[k] for k in ks[-3:]] if ks else []
        est = float(np.mean(tail)) if tail else float('nan')
        band = (max(tail) - min(tail)) if len(tail) > 1 else float('nan')
        rq[q] = est
        log(f"   => r_{q} ~ {est:.4f}  (tail spread +/-{band/2 if band==band else float('nan'):.4f})")
        log("")

    rq[3] = 1.0
    log("## r_q TABLE")
    log(f"   {'q':>4} {'r_q (measured)':>16} {'3/q':>10} {'a^(q-3) OOS':>14}")
    # fit H_EXP on q=5,7
    import math
    x = [5, 7]; y = [math.log(rq[5]), math.log(rq[7])]
    b = (y[1] - y[0]) / (x[1] - x[0]); a0 = y[0] - b * x[0]
    def hexp(q): return math.exp(a0 + b * q)
    log(f"   (H_EXP fit on q=5,7: r_q = exp({a0:.4f} + {b:.4f} q); a=exp(b)={math.exp(b):.4f})")
    for q in [3, 5, 7, 11, 13]:
        oos = hexp(q) if q in (11, 13) else (hexp(q))
        tag = " <-OOS" if q in (11, 13) else ""
        log(f"   {q:>4} {rq[q]:>16.4f} {3/q:>10.4f} {hexp(q):>14.4f}{tag}")
    log("")
    log("## VERDICT")
    for q in [11, 13]:
        e_rat = abs(rq[q] - 3 / q)
        e_exp = abs(rq[q] - hexp(q))
        log(f"   q={q}: |r-3/q|={e_rat:.4f}   |r-a^(q-3)|={e_exp:.4f}   "
            f"-> {'3/q closer' if e_rat < e_exp else 'exp closer'}")
    log("")
    log("   READ: whichever form matches the OOS points q=11,13 within the oscillation band")
    log("   wins; if neither, r_q has no elementary closed form (= operator eigenvalue ratio).")
    flush()


def flush():
    with open("result_28_rq_formula_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
