"""
PROBE 30 -- pin r_11, r_13 via matrix-pencil (Prony) on exact c_k + a k=6 compute extension.

PRE-REGISTRATION.
------------------------------------------------------------------
R28 failed to pin r_11,r_13 because rho_k=c_{k+1}/c_k oscillates (R26: a +r/-r complex-mode
pair) and only 3 transient points were available. FIX (two levers):
  (1) MATRIX PENCIL: c_k = sum_i beta_i z_i^k (complex geometric modes). Fit the modes
      directly from exact c_k; r_q = |z_dominant| -- robust to oscillation, uses existing data.
  (2) COMPUTE: extend q=11 to k=6 (~2.5 GB, v-truncated int32) and attempt q=13 k=6
      (guarded MemoryError) for one more c-value each.

GATE: cross(k) matches R22/R23/R27 exact at low k (<1e-8). Reuses the validated v-truncation.

DECISION: report r_q from BOTH rho-tail and matrix-pencil; they should agree for q=5,7
(where rho settled). For q=11,13 the pencil is the estimate. Then re-test candidates
(3/q, exp) on the pinned {r_5,r_7,r_11,r_13}. NO closed form is committed (priors 0-for-8);
report which if any survives, else "no elementary form" stands.

NOT AT STAKE: R10-R29.
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


def stat_trunc_lean(q, k, vmax=52):
    """Memory-lean v-truncated stationary: int32 indices, build K^T directly."""
    N = q ** k
    inv2 = pow(2, -1, N)
    M = order_of_two(N)
    vm = min(vmax, M)
    cp = np.array([r for r in range(N) if gcd(r, q) == 1], dtype=np.int64)
    n = len(cp)
    Z = 1.0 - 2.0 ** (-M)
    inv_idx = np.full(N, -1, dtype=np.int32)
    inv_idx[cp] = np.arange(n, dtype=np.int32)
    base_t = (q * cp + 1) % N
    src = np.arange(n, dtype=np.int32)
    rows, cols, vals = [], [], []      # rows=dst (Kt row), cols=src
    inv2v = 1
    for v in range(1, vm + 1):
        inv2v = (inv2v * inv2) % N
        t = (base_t * inv2v) % N
        rows.append(inv_idx[t])
        cols.append(src)
        vals.append(np.full(n, (0.5 ** v) / Z, dtype=np.float64))
    Kt = sp.csr_matrix((np.concatenate(vals),
                        (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
    del rows, cols, vals
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


def cross_lean(q, k, vmax=52):
    pi, n = stat_trunc_lean(q, k, vmax)
    M = order_of_two(q ** k)
    nrm = float(np.dot(pi, pi))
    diag = float(sum_p2_exact(M) ** k)
    return nrm / diag - 1.0 - ratio_within(q, k), n


def matrix_pencil(c, order):
    """Prony via linear predictor: c_n = sum_{i=1}^order a_i c_{n-i}; return roots by |.|."""
    N = len(c)
    if N < order + 1:
        return []
    rows = [[c[n - i] for i in range(1, order + 1)] for n in range(order, N)]
    rhs = [c[n] for n in range(order, N)]
    a, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
    comp = np.zeros((order, order))
    comp[0, :] = a
    for i in range(1, order):
        comp[i, i - 1] = 1.0
    z = np.linalg.eigvals(comp)
    return sorted(z, key=lambda t: -abs(t))


REF = {(5, 2): 0.038875214298, (7, 2): 0.047025340793,
       (11, 2): 0.000113752827, (13, 3): 0.000196127139}


def main():
    log("# PROBE 30 -- pin r_11,r_13 via matrix-pencil + k=6 extension")
    log("")
    KMAX = {5: 8, 7: 7, 11: 6, 13: 6}
    VMAX = {5: 64, 7: 64, 11: 52, 13: 44}
    NCAP = {5: 5e7, 7: 5e7, 11: 1.2e8, 13: 2.6e8}
    rq = {3: 1.0}
    for q in [5, 7, 11, 13]:
        log(f"## q={q}")
        cr = {}
        for k in range(2, KMAX[q] + 1):
            n = (q - 1) * q ** (k - 1)
            if n * VMAX[q] > NCAP[q]:
                log(f"   k={k}: n*vmax={n*VMAX[q]:.1e} > cap -- STOP (SAID)")
                break
            try:
                c, nn = cross_lean(q, k, VMAX[q])
            except MemoryError:
                log(f"   k={k}: MemoryError (n={n}) -- SKIP, needs Lambda (SAID)")
                break
            cr[k] = c
            g = ""
            if (q, k) in REF:
                rel = abs(c - REF[(q, k)]) / abs(REF[(q, k)])
                g = f"  gate|rel|={rel:.1e}{' OK' if rel<1e-7 else ' CHK'}"
            log(f"   cross({k})={c:.12f} (n={nn}){g}")
        ck = [cr[k] - cr.get(k - 1, 0.0) for k in sorted(cr)]
        ks = sorted(cr)
        # rho tail
        rhos = [ck[i + 1] / ck[i] for i in range(len(ck) - 1) if abs(ck[i]) > 1e-12]
        log(f"   rho_k: {['%.5f' % r for r in rhos]}")
        # matrix pencil on c_k (try orders 2 and 3)
        log(f"   c_k: {['%.4e' % v for v in ck]}")
        for order in [2, 3]:
            if len(ck) >= 2 * order + 1 or len(ck) >= order + 2:
                z = matrix_pencil(ck, order)
                zs = "  ".join(f"{complex(zz):.4f}" for zz in z)
                dom = abs(z[0]) if z else float('nan')
                log(f"   pencil order-{order}: modes |z|={['%.4f'%abs(zz) for zz in z]}  dom={dom:.4f}")
        # best estimate: prefer order-2 dominant |z| (the physical +/- pair modulus)
        z2 = matrix_pencil(ck, 2)
        rq[q] = abs(z2[0]) if z2 else (rhos[-1] if rhos else float('nan'))
        log(f"   => r_{q} ~ {rq[q]:.4f}  (pencil order-2 dominant)")
        log("")

    # ---- candidate test on pinned table ----
    log("## r_q TABLE (pinned) + candidate closed forms")
    import math
    log(f"   {'q':>4} {'r_q':>10} {'3/q':>10} {'1/(q-3)?':>12} {'2/(q-1)':>10}")
    for q in [3, 5, 7, 11, 13]:
        v = rq.get(q, float('nan'))
        alt1 = 1.0 / (q - 3) if q != 3 else float('inf')
        alt2 = 2.0 / (q - 1)
        log(f"   {q:>4} {v:>10.4f} {3/q:>10.4f} {alt1:>12.4f} {alt2:>10.4f}")
    log("")
    # fit exp on 5,7 predict 11,13 (OOS)
    if all(k in rq for k in (5, 7, 11, 13)):
        b = (math.log(rq[7]) - math.log(rq[5])) / 2
        a0 = math.log(rq[5]) - b * 5
        log(f"   H_EXP (fit 5,7): r_q=exp({a0:.3f}{b:+.3f}q)")
        for q in (11, 13):
            log(f"      q={q}: exp pred {math.exp(a0+b*q):.4f}  vs measured {rq[q]:.4f}  "
                f"vs 3/q {3/q:.4f}")
    log("")
    log("## READ: pencil r_q robust to oscillation. If a candidate hits all of {5,7,11,13}")
    log("   within ~0.02, it's the law; else 'no elementary closed form' stands (R28).")
    flush()


def flush():
    with open("result_30_rq_pin_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
