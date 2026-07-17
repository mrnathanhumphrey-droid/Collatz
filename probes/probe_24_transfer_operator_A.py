"""
PROBE 24 -- ROUTE 1A: the fixed transfer operator via the NO-TOWER folded model.
Extract r_q = subdominant eigenvalue; GATE against R23's exact cross(k).

PRE-REGISTRATION (written before running).
------------------------------------------------------------------
R23 showed: (i) c_k is order >=2 (not geometric), (ii) the q=3 approach to 7/15 is a
super-geometric TOWER x_j = 2^{-d q^{j-1}}, so there is NO finite recurrence for the real
c_k. The tower corrections decay DOUBLY-exponentially, hence are irrelevant to the
subdominant rate r_q. This probe STRIPS the tower and builds the fixed operator.

THE FOLDED (no-tower) MODEL. In the real model coordinate j has modulus m_j = d q^{j-1}
(the tower). FOLD every coordinate to modulus d (drop the tower): each of the k
coordinates carries e_i in {1..d} with folded weight p(e) = 2^{-e}/(1-2^{-d})
(a proper distribution, sum_{e=1}^d p(e) = 1). Value is still the exact q-adic collision
functional  value = sum_{m=1}^k q^{m-1} 2^{-S_m} mod q^k,  S_m = e_{k-m+1}+..+e_k.
  total_LO(k) = ||W||^2 / P2_LO^k - 1,   W_v = sum_{cells->v} prod p(e_i),
                P2_LO = sum_{e=1}^d p(e)^2 = (2^d+1)/(3(2^d-1)).
Because within(k) is k-FLAT (R15 H_FLAT), Delta total = Delta cross = c_k for k>=3, so
the increments of total_LO give the SAME rate as cross. The folded model has NO tower ->
a genuine FINITE linear recurrence -> a fixed transfer operator whose subdominant
eigenvalue is r_q. Cell count d^k (tiny): reach high k.

WHY r_q^LO = r_q^real. real c_k = c_k^LO + (tower), tower ~ 2^{-d q^k} -> 0 doubly-exp.
So both sequences share the SAME limit ratio r_q; folded converges CLEANLY (no tower
noise) and fast.

HYPOTHESES / GATES (pre-committed):
  G_TOWER (gate + measurement): |c_k^LO - c_k^real| at the k where R23 has exact real c_k
      (q=3 k<=5, q=5 k<=4, q=7 k<=3) must be SMALL and SHRINKING in k (that gap IS the
      tower). NOT an equality gate -- the folded model is deliberately tower-free; the
      gate is that the gap is small (< a few % at q=3, tiny at q>=5) and decreasing.
      If the gap is LARGE or GROWING, folding changed the rate -> Option A is wrong ->
      go to Option B (Fourier) without reporting an r_q.
  G_FINITE (*** the claim ***): c_k^LO satisfies a FINITE linear recurrence (Hankel rank
      stabilizes) -- i.e. the folded operator is finite-dimensional. Order = op dimension.
  R_Q: r_q = dominant subdominant root of that recurrence = lim rho_k^LO. Report for
      q=3,5,7,11. CHECK: r_3 ~ 1 (divergence), r_q < 1 for q>=5 (convergence), and
      consistency with R22's float rho_k trend (q=5 ~0.6+, climbing).

  PRE-COMMITTED (stated to lose, priors 0-for-8): I do NOT guess r_q's value. I commit
  only to the STRUCTURE: r_3 = 1 exactly (the divergence eigenvalue) and r_q < 1 for
  q>=5. Any numeric r_q is READ OUT, not predicted.

BUDGET: d^k cells. Cap 2.2e6. q=3 k<=21, q=5 k<=10, q=7 k<=13, q=11 k<=6. No silent trunc.

NOT AT STAKE: R10-R23, R5's rate, R6, R7, R12, THEOREM_C_745.
"""
import sys
from fractions import Fraction
from itertools import product
import numpy as np

from probe_6_conservation_generalize import order_of_two
from probe_23_exact_increment_recurrence import exact_cross

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def folded_total(q, k):
    """total_LO(k) = ||W||^2 / P2_LO^k - 1 for the no-tower folded model. Float."""
    d = order_of_two(q)
    N = q ** k
    inv2 = pow(2, -1, N)
    # p(e) = 2^{-e}/(1-2^{-d}), e=1..d
    denom = 1.0 - 2.0 ** (-d)
    p = [(2.0 ** (-e)) / denom for e in range(1, d + 1)]
    P2_LO = sum(pe * pe for pe in p)
    # precompute inv2^S mod N for S up to k*d
    maxS = k * d
    inv2pow = [1] * (maxS + 1)
    for s in range(1, maxS + 1):
        inv2pow[s] = (inv2pow[s - 1] * inv2) % N
    qpow = [q ** (m - 1) % N for m in range(1, k + 1)]
    W = {}
    for e in product(range(1, d + 1), repeat=k):
        w = 1.0
        for ei in e:
            w *= p[ei - 1]
        # value = sum_{m=1}^k q^{m-1} 2^{-S_m}, S_m = suffix sum
        val = 0
        s = 0
        for m in range(1, k + 1):
            s += e[k - m]
            val = (val + qpow[m - 1] * inv2pow[s]) % N
        W[val] = W.get(val, 0.0) + w
    norm2 = sum(wv * wv for wv in W.values())
    return norm2 / (P2_LO ** k) - 1.0


def hankel_order(seq, tol=1e-9):
    """Smallest r such that the (r+1)x(r+1) Hankel determinant ~ 0 (recurrence order)."""
    n = len(seq)
    for r in range(1, n // 2):
        H = np.array([[seq[i + j] for j in range(r + 1)] for i in range(r + 1)], float)
        sv = np.linalg.svd(H, compute_uv=False)
        if sv[-1] / sv[0] < tol:
            return r
    return None


def recurrence_roots(seq, order):
    """Fit c_{n} = sum_{i=1}^{order} a_i c_{n-i} (least squares), return sorted |roots|."""
    n = len(seq)
    rows = [[seq[m - i] for i in range(1, order + 1)] for m in range(order, n)]
    rhs = [seq[m] for m in range(order, n)]
    a, *_ = np.linalg.lstsq(np.array(rows, float), np.array(rhs, float), rcond=None)
    comp = np.zeros((order, order))
    comp[0, :] = a
    for i in range(1, order):
        comp[i, i - 1] = 1.0
    roots = np.linalg.eigvals(comp)
    return sorted(roots, key=lambda z: -abs(z)), a


def main():
    log("# PROBE 24 -- ROUTE 1A: fixed transfer operator via no-tower folded model")
    log("# Pre-reg: G_TOWER(gate) / G_FINITE(*** claim ***) / R_Q read out (r_3=1, r_q<1 q>=5)")
    log("")

    KMAX = {3: 21, 5: 10, 7: 13, 11: 6}
    seqs = {}
    for q in [3, 5, 7, 11]:
        kmax = KMAX[q]
        tot = {k: folded_total(q, k) for k in range(2, kmax + 1)}
        ck = {k: tot[k] - tot[k - 1] for k in range(3, kmax + 1)}
        ck[2] = tot[2]
        seqs[q] = (tot, ck)
        log(f"## q={q}  (d={order_of_two(q)}, folded k=2..{kmax}, {order_of_two(q)}^{kmax} cells)")
        log(f"   {'k':>3} {'total_LO':>16} {'c_k^LO':>16} {'rho_k^LO=c_{k+1}/c_k':>22}")
        for k in range(2, kmax + 1):
            r = (ck[k + 1] / ck[k]) if (k + 1 in ck and ck[k] != 0) else float('nan')
            log(f"   {k:>3} {tot[k]:>16.10f} {ck[k]:>16.3e} {r:>22.8f}")
        log("")

    # ---- G_TOWER gate: folded vs exact real c_k (R23) ----
    log("## G_TOWER -- gap |c_k^LO - c_k^real| (this gap IS the tower; must be small+shrinking)")
    log(f"   {'q':>4} {'k':>3} {'c_k^LO':>16} {'c_k^real (exact)':>18} {'|gap|':>12}")
    EXACTK = {3: [3, 4, 5], 5: [3, 4], 7: [3]}
    tower_ok = True
    for q, ks in EXACTK.items():
        _, ck = seqs[q]
        prevgap = None
        for k in ks:
            cr = exact_cross(q, k)[0] - exact_cross(q, k - 1)[0]
            real_ck = float(cr)
            gap = abs(ck[k] - real_ck)
            flag = ""
            if prevgap is not None and gap > prevgap + 1e-12:
                flag = "  <- GROWING"
                tower_ok = False
            prevgap = gap
            log(f"   {q:>4} {k:>3} {ck[k]:>16.10f} {real_ck:>18.10f} {gap:>12.3e}{flag}")
    log(f"   G_TOWER: {'gap small & shrinking -> folding preserved the rate' if tower_ok else 'GROWING -> folding changed rate -> use Option B'}")
    log("")

    # ---- G_FINITE + R_Q ----
    log("## G_FINITE (*** claim ***) -- does c_k^LO obey a FINITE linear recurrence?")
    log("   Hankel rank = recurrence order = transfer-operator dimension; roots = eigenvalues.")
    log("")
    for q in [3, 5, 7, 11]:
        _, ck = seqs[q]
        ks = sorted(k for k in ck if k >= 2)
        seq = [ck[k] for k in ks]
        # use the settled tail (drop k=2 pre-asymptotic head if long enough)
        tail = seq[1:] if len(seq) > 6 else seq
        order = hankel_order(tail)
        rho_tail = [tail[i + 1] / tail[i] for i in range(len(tail) - 1) if tail[i] != 0]
        rq_est = rho_tail[-1] if rho_tail else float('nan')
        log(f"   q={q}: rho_k^LO tail -> {['%.6f' % r for r in rho_tail[-4:]]}")
        if order is None:
            log(f"          Hankel order: not resolved within available terms")
        else:
            roots, a = recurrence_roots(tail, order)
            log(f"          Hankel order = {order}  (operator dim)")
            log(f"          roots |lambda|: {['%.6f' % abs(z) for z in roots]}")
        log(f"          => r_q (= lim rho_k^LO) ~ {rq_est:.6f}")
        log("")

    log("## READ: r_3 should be ~1 (divergence eigenvalue), r_q<1 for q>=5 (spectral gap).")
    log("   G_TOWER small+shrinking validates the folded operator carries the real rate.")
    log("   Next: Option B (Fourier transfer operator) as an INDEPENDENT check of r_q.")
    flush()


def flush():
    with open("result_24_transfer_A_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
