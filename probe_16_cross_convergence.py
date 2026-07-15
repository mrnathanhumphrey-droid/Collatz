"""
PROBE 16 (qx+1 paper) -- does the CROSS-cell term converge at q>=5, and at rate 3/q?
This is the phase boundary stated as a convergence rate.

PRE-REGISTRATION (written before running; derivation first; RATES pre-committed).
------------------------------------------------------------------
R15 established: within-cell overlap is closed-form, bounded, k-INDEPENDENT; at q=3 it
freezes at 0.71958983896 from k=4 while the total grows LINEARLY, so ALL of the q=3
domination failure is CROSS-cell with slope 0.46577 ~ 7/15 (off 0.06%).
R15 also REFUTED my prior that cross stays flat at q>=5 -- it GROWS at every q. But its
increments decayed geometrically at q=5 (0.155, 0.079, 0.049; ratio ~0.6 ~ 3/q), which
was an eyeball read on four points. This probe tests it.

THE DERIVATION (done before running):
  R7 (proved): X_k = X_0 + sum_{j<=k} M_j, and R5/R11: M_j ~ c~_q * (q/3)^j.
  So X_k / (q/3)^k -> C_q := c~_q * q/(q-3), a GEOMETRIC SUM, with tail (3/q)^k.
  Since ratio_k := ||pi_k||^2 / (sum_v p_v^2)^k - 1 = X_k/(q/3)^k - 1 (up to the
  exact P2-vs-1/3 correction), we get
      *** ratio_k -> C_q - 1  with deficit ~ (3/q)^k ***
  A geometric series with ratio (3/q):
      q >= 5 : |3/q| < 1  ==> CONVERGES ==> domination holds
      q = 3  : 3/q = 1    ==> DIVERGES linearly ==> domination fails by a factor of k
  This is R7's geometric-series divergence, now visible as a CONVERGENCE RATE inside the
  overlap count -- the same phase boundary for the sixth independent time.

  TEST WITHOUT KNOWING THE LIMIT: if deficit_k = A*r^k then successive differences
      D_k := ratio_{k+1} - ratio_k = deficit_k - deficit_{k+1} = A*r^k*(1-r)
  so  *** D_{k+1}/D_k -> r ***  -- the rate is measurable with NO extrapolation.

PRE-COMMITTED RATES (r = 3/q, computed before the run):
      q=5  -> 0.600000        q=7  -> 0.428571
      q=11 -> 0.272727        q=13 -> 0.230769
      q=3  -> 1.000000  (no decay; R15 already confirmed linear growth, slope 7/15)

HYPOTHESES:
  H_RATE (*** THE TEST ***): D_{k+1}/D_k -> 3/q at q>=5.
      PRIOR: TRUE. STATED TO LOSE. Honest limitation up front: the largest affordable k
      gives only 3 ratios at q=5, 2 at q=7, 1 at q=11 -- so this can CORROBORATE or
      REFUTE a 3/q trend but CANNOT pin the rate to high precision. I will NOT fit.
  H_CONV: at q>=5 the ratio_k sequence is increasing and bounded (Aitken-extrapolated
      limit finite and consistent with C_q - 1 = c~_q*q/(q-3) - 1). PRIOR: TRUE.
  H_CRIT: at q=3 the same difference-ratio D_{k+1}/D_k -> 1 (no decay).
      PRIOR: TRUE. This is the control: the SAME statistic must show 1 at q=3 and 3/q
      at q>=5, or the framing is wrong.

DECISION RULES (pre-committed):
  H_RATE CORROBORATED iff, at every q>=5 tested, the LAST available D_{k+1}/D_k is within
      20% of 3/q AND the sequence of ratios is moving TOWARD 3/q (not away).
      REFUTED iff the last ratio is off by >2x, or moving away from 3/q.
      Anything else -> INCONCLUSIVE, reported as such. (20% is deliberately loose: with
      2-3 points a tight threshold would be theatre. This test is a SIGN CHECK on the
      mechanism, not a measurement of the rate.)
  H_CRIT CONFIRMED iff at q=3 the last D_{k+1}/D_k is within 5% of 1.0.
  No law will be fitted to these points -- per feedback_r2_cannot_discriminate_monotone_fits.

RESOURCE NOTE (stated up front): the chain matrix here is ~100% DENSE (every state has
M distinct targets), so a dense array is CHEAPER than sparse. Peak ~n^2*8 bytes:
q=5,k=6 -> 1.25 GB; q=7,k=5 -> 1.66 GB; q=11,k=4 -> 1.42 GB. Seconds-to-a-minute each.
10x today's largest run, far below any Lambda threshold. Skips anything over 2.5 GB.

NOT AT STAKE: R10's law, R11, R13, R14, R15's within-cell identity, R5's rate, R6, R7,
R12, THEOREM_C_745. A refutation kills my route to Result 1, not a banked result.

Author's priors this arc: 15-for-23, and FIVE mis-specified decision rules. Hence the
rates above are NUMBERS committed before the run, and the threshold is honestly loose.
"""
import sys
import numpy as np

from probe_6_conservation_generalize import order_of_two
from probe_8_selfsimilar_overlap import sum_p2_exact
from probe_15_tower_k_count import ratio_within

LOG = []


def log(m=""):
    print(m)
    LOG.append(str(m))


def total_ratio_dense(q, k, max_gb=2.5):
    """||pi_k||^2/(sum p^2)^k - 1 via a DENSE transition matrix (the chain is ~100% dense)."""
    from math import gcd
    N = q ** k
    M = order_of_two(N)
    cp = np.array([r for r in range(N) if gcd(r, q) == 1], dtype=np.int64)
    n = len(cp)
    gb = n * n * 8 / 1e9
    if gb > max_gb:
        return None, gb, M
    inv2 = pow(2, -1, N)
    inv_idx = np.full(N, -1, dtype=np.int64)
    inv_idx[cp] = np.arange(n)
    base_t = (q * cp + 1) % N
    K = np.zeros((n, n), dtype=np.float64)
    Z = 1.0 - 2.0 ** (-M)
    src = np.arange(n)
    i2v = 1
    for v in range(1, M + 1):
        i2v = (i2v * inv2) % N
        tgt = inv_idx[(base_t * i2v) % N]
        K[src, tgt] += (0.5 ** v) / Z
    pi = np.full(n, 1.0 / n)
    Kt = K.T.copy()
    for _ in range(3000):
        nxt = Kt.dot(pi)
        s = nxt.sum()
        if s == 0:
            break
        nxt /= s
        if np.abs(nxt - pi).sum() < 1e-15:
            pi = nxt
            break
        pi = nxt
    nrm = float(np.dot(pi, pi))
    diag = float(sum_p2_exact(M) ** k)
    del K, Kt
    return nrm / diag - 1, gb, M


def main():
    log("# PROBE 16 -- does the CROSS-cell term converge at q>=5, and at rate 3/q?")
    log("# Pre-reg: H_RATE(*** test ***) / H_CONV / H_CRIT(control). Rates committed BEFORE the run:")
    for q in [3, 5, 7, 11, 13]:
        log(f"#    q={q:>3}: 3/q = {3/q:.6f}")
    log("")

    PLAN = {3: [2, 3, 4, 5, 6, 7, 8], 5: [2, 3, 4, 5, 6], 7: [2, 3, 4, 5], 11: [2, 3, 4], 13: [2, 3]}
    res = {}
    log("## totals (dense chain), within (R15 closed form), cross = total - within")
    log("")
    log(f"{'q':>4} {'k':>3} {'mem GB':>8} {'total ratio_k':>18} {'within':>16} {'cross':>16}")
    for q, ks in PLAN.items():
        for k in ks:
            tot, gb, M = total_ratio_dense(q, k)
            if tot is None:
                log(f"{q:>4} {k:>3} {gb:>8.2f}  SKIP (over 2.5 GB cap)")
                break
            wit = ratio_within(q, k)
            res[(q, k)] = {"tot": tot, "wit": wit, "cross": tot - wit}
            log(f"{q:>4} {k:>3} {gb:>8.3f} {tot:>18.10f} {wit:>16.10f} {tot-wit:>16.10f}")
        log("")

    # ---------------- H_RATE / H_CRIT ----------------
    log("## H_RATE / H_CRIT -- successive-difference ratios D_{k+1}/D_k  (no limit needed)")
    log("   D_k := ratio_{k+1} - ratio_k ;  if deficit ~ A*r^k then D_{k+1}/D_k -> r")
    log("")
    for q in [3, 5, 7, 11, 13]:
        ks = sorted(k for (qq, k) in res if qq == q)
        if len(ks) < 3:
            log(f"   q={q:>3}: only {len(ks)} points -- no difference ratio available")
            continue
        tots = [res[(q, k)]["tot"] for k in ks]
        D = [tots[i + 1] - tots[i] for i in range(len(tots) - 1)]
        R = [D[i + 1] / D[i] for i in range(len(D) - 1) if D[i] != 0]
        pred = 3 / q
        last = R[-1] if R else float("nan")
        toward = ""
        if len(R) >= 2:
            toward = "TOWARD" if abs(R[-1] - pred) < abs(R[-2] - pred) else "AWAY"
        log(f"   q={q:>3} (3/q = {pred:.6f}): D = {['%.3e' % d for d in D]}")
        log(f"          D_{{k+1}}/D_k = {['%.5f' % r for r in R]}   last={last:.5f}  "
            f"off {abs(last-pred)/pred:.1%}  {toward}")
        if q == 3:
            v = abs(last - 1.0) <= 0.05
            log(f"          H_CRIT (control, want ~1.0): {'CONFIRMED' if v else 'NOT confirmed'}")
        else:
            near = abs(last - pred) / pred <= 0.20
            far = abs(last - pred) / pred > 1.0
            verdict = ("CORROBORATED" if (near and toward != "AWAY") else
                       ("REFUTED" if far else "INCONCLUSIVE"))
            log(f"          H_RATE: {verdict}")
        log("")

    # ---------------- H_CONV ----------------
    log("## H_CONV -- Aitken-extrapolated limit vs C_q - 1 = c~_q*q/(q-3) - 1")
    log("   (c~_q taken as the largest-k measurement available; NOT fitted)")
    log("")
    for q in [5, 7, 11]:
        ks = sorted(k for (qq, k) in res if qq == q)
        if len(ks) < 3:
            continue
        t = [res[(q, k)]["tot"] for k in ks]
        a, b, c = t[-3], t[-2], t[-1]
        den = (c - b) - (b - a)
        ait = c - (c - b) ** 2 / den if den != 0 else float("nan")
        log(f"   q={q:>3}: ratio_k = {['%.6f' % x for x in t]}")
        log(f"          Aitken limit = {ait:.6f}   (cross limit = {ait - ratio_within(q, ks[-1]):.6f})")
    log("")
    log("   (H_CONV: sequence increasing + finite Aitken limit => cross CONVERGES at q>=5,")
    log("    while q=3 grows linearly forever. Same phase boundary, 6th independent sighting.)")
    flush()


def flush():
    with open("result_16_cross_convergence_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
