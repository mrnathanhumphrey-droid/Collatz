"""
PROBE 27 -- ROUTE 1 (option 3): high-k exact rho_k for q=5,7 via v-TRUNCATED stationary.
Settle r_q: does rho_k -> ~0.60 (R26 modal, 3/q) or keep climbing?

PRE-REGISTRATION (written before running).
------------------------------------------------------------------
R22 stopped at k=5 (q=5) because the full stationary needs n*M nonzeros with M=ord_2(q^k)
(=12500 at q=5 k=6 -> 156M, heavy). BUT the transition weight is 2^{-v}; v>64 contributes
<2^{-64}~5e-20. TRUNCATE the v-sum at vmax=64 -> matrix is n*vmax (~800k at q=5 k=6),
trivial, and the stationary pi is identical to ~5e-20. This lets us reach k=8 (q=5),
k=7 (q=7) and read rho_k directly -- the ground-truth test of R26's modal r_5~0.603.

GATE (validates the truncation): v-truncated cross(k) must equal R22/R23 exact cross(k)
for k<=5 (q=5), k<=3 (q=7) to <1e-10. If it deviates, vmax is too small -> raise it.

  cross(k) = ||pi_k||^2 / P2^k - 1 - ratio_within(q,k),  pi from v-truncated power iteration.
  c_k = cross(k)-cross(k-1);  rho_k = c_{k+1}/c_k.

PRE-COMMITTED (structure only, value NOT predicted -- priors 0-for-8):
  If R26 is right, rho_k(q=5) SETTLES near 0.60 (possibly oscillating, since a mu=-0.60
  mode accompanies +0.603). If instead rho_k keeps CLIMBING toward 1, R26's amplitude read
  was wrong and r_5 is near 1 (R25). This probe DISCRIMINATES.
  Numerical floor: rho unreliable once |c_k| < 1e-10 (flagged). For r~0.6, c_k~0.6^k stays
  >1e-10 to k~13, well beyond reach -> clean.

BUDGET: n*vmax. q=5: k<=8 (n=312500, 20M) ; q=7: k<=7 (n=705894, 45M). Cap n*vmax<=5e7.
  No heavy compute; no Lambda. SAID per q where it stops.

NOT AT STAKE: R10-R26, R5's rate, R6, R7, R12, THEOREM_C_745.
"""
import sys
from math import gcd
import numpy as np
import scipy.sparse as sp

from probe_6_conservation_generalize import order_of_two
from probe_8_selfsimilar_overlap import sum_p2_exact
from probe_15_tower_k_count import ratio_within

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def stationary_trunc(q, k, vmax=64):
    """Stationary pi with the v-sum truncated at vmax (weights 2^{-v}, v>vmax negligible)."""
    N = q ** k
    inv2 = pow(2, -1, N)
    M = order_of_two(N)
    vm = min(vmax, M)
    cp = np.array([r for r in range(N) if gcd(r, q) == 1], dtype=np.int64)
    n = len(cp)
    Z = 1.0 - 2.0 ** (-M)                    # = (2^M-1)/2^M, no overflow for large M
    inv_idx = np.full(N, -1, dtype=np.int64)
    inv_idx[cp] = np.arange(n)
    base_t = (q * cp + 1) % N
    rows_l, cols_l, vals_l = [], [], []
    inv2v = 1
    src = np.arange(n)
    for v in range(1, vm + 1):
        inv2v = (inv2v * inv2) % N
        t = (base_t * inv2v) % N
        rows_l.append(src)
        cols_l.append(inv_idx[t])
        vals_l.append(np.full(n, (0.5 ** v) / Z))
    K = sp.csr_matrix((np.concatenate(vals_l),
                       (np.concatenate(rows_l), np.concatenate(cols_l))), shape=(n, n))
    Kt = K.T.tocsr()
    pi = np.full(n, 1.0 / n)
    for _ in range(4000):
        nxt = Kt.dot(pi)
        s = nxt.sum()
        if s == 0:
            break
        nxt /= s
        if np.abs(nxt - pi).sum() < 1e-15:
            pi = nxt
            break
        pi = nxt
    return pi, n


def cross_trunc(q, k, vmax=64):
    pi, n = stationary_trunc(q, k, vmax)
    M = order_of_two(q ** k)
    nrm = float(np.dot(pi, pi))
    diag = float(sum_p2_exact(M) ** k)
    return nrm / diag - 1.0 - ratio_within(q, k), n


def main():
    log("# PROBE 27 -- high-k rho_k via v-truncated stationary. Settle r_q (q=5,7).")
    log("# Pre-reg: GATE (trunc==exact k<=5) / discriminate rho->0.60 vs ->1")
    log("")

    # R22/R23 exact-ish references for the gate
    REF = {(5, 2): 0.038875214298, (5, 3): 0.059623775122, (5, 4): 0.070157562299,
           (5, 5): 0.076734752,    (7, 2): 0.047025340793, (7, 3): 0.068042081671}

    KMAX = {5: 8, 7: 7}
    NCAP = 5e7
    for q in [5, 7]:
        log(f"## q={q}")
        cr = {}
        for k in range(2, KMAX[q] + 1):
            n = (q - 1) * q ** (k - 1)
            if n * 64 > NCAP:
                log(f"   k={k}: n*vmax={n*64:.1e} > cap -- STOP (SAID)")
                break
            c, nn = cross_trunc(q, k)
            cr[k] = c
            gate = ""
            if (q, k) in REF:
                rel = abs(c - REF[(q, k)]) / abs(REF[(q, k)])
                gate = f"  gate vs exact: |rel|={rel:.2e} {'OK' if rel < 1e-9 else ('ok~' if rel<1e-6 else 'FAIL')}"
            log(f"   cross({k}) = {c:.12f}   (n={nn}){gate}")
        # increments and rho
        ck = {k: cr[k] - cr.get(k - 1, 0.0) for k in cr if k >= 2}
        ck[min(cr)] = cr[min(cr)]  # cross(1)=0
        log("")
        log(f"   {'k':>3} {'c_k':>16} {'rho_k=c_{k+1}/c_k':>20}")
        for k in sorted(ck):
            if k + 1 in ck and abs(ck[k]) > 1e-10:
                r = ck[k + 1] / ck[k]
                flag = "" if abs(ck[k + 1]) > 1e-10 else "  (floor)"
                log(f"   {k:>3} {ck[k]:>16.6e} {r:>20.8f}{flag}")
            else:
                log(f"   {k:>3} {ck[k]:>16.6e} {'--':>20}")
        log("")

    log("## READ: if rho_k(q=5) settles near 0.60 -> R26 modal confirmed, 3/5 vindicated.")
    log("   if rho_k climbs toward 1 -> R25 (r near 1) and R26 modal was wrong.")
    log("   oscillation around 0.60 is EXPECTED (R26: mu=+0.603 with a mu=-0.60 partner).")
    flush()


def flush():
    with open("result_27_high_k_rho_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
