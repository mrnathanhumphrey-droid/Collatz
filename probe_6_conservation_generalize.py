"""
PROBE 6 (qx+1 paper) -- scope the domination identity: does R76 generalize to (Z/q^k)*?

PRE-REGISTRATION (written before running; falsifier first).
------------------------------------------------------------------
Target: R5 flagged "uniform diagonal-self-overlap domination on (Z/q^k)*, an
R76-style conservation identity generalized from q=3" as the ONE line left to
upgrade S_k^(q) ~ (q/3)^k from mechanism to theorem.

Reading result_76_conservation_law.md, R76 is TWO theorems, and they have very
different generalization prospects. This probe separates them.

  Thm 76.1 (CONSERVATION):  sum_{j=0..2} M_{n+1}(eta_0 + j*3^n) = 0.
    Proof = complete-character-sum vanishing. The j-sum is
      sum_{j=0}^{q-1} e^{2pi i r xi j / q},
    which is 0 unless q | r*xi. With gcd(r,q)=1 (chain support) this is q ∤ xi.
    NOTHING in that argument uses q=3.
    => H_CONS: conservation generalizes VERBATIM to every odd q.

  Thm 76.3 (LEADING-MODE):  S_{n+1} = -2 * M_{n+1}(1 + 3^n).
    Proof needs Lemma 76.2: among the 3 lifts of eta_0, exactly ONE is
    self-inverse and the other TWO are mutual inverses. Combined with
    M(eta)=M(eta^-1) real, conservation becomes S + 2M = 0 -> S = -2M.
    For general odd prime q there are q lifts of eta_0=1: lift 1 is
    self-inverse, and -1 is NOT a lift of 1 (for q^k > 2), so the other q-1
    lifts form (q-1)/2 mutual-inverse PAIRS.
      q=3 -> (q-1)/2 = 1 pair  -> conservation = 1 equation, 1 unknown -> SOLVED.
      q>=5 -> (q-1)/2 >= 2 pairs -> conservation = 1 equation, >=2 unknowns.
    => H_LEAD_BREAKS: 76.3 does NOT port; q=3 is special precisely because
       (q-1)/2 = 1. This is a REAL obstruction, not a missing line.

  BUT the underdetermination dissolves if the pair-moments coincide:
    => H_EQUAL (the falsifier / the hopeful branch): M_{k+1}(1 + j*q^k) is
       INDEPENDENT of j in {1..q-1}. Then conservation gives
         S_{k+1} = -(q-1) * M_{k+1}(1 + q^k),
       which at q=3 reproduces 76.3's "-2" as "-(q-1)". If H_EQUAL holds, the
       leading-mode identity generalizes after all and the "-2" was never a 2.

Falsifier-first discipline: H_EQUAL is the claim that would SAVE the port, so it
is the one tested hardest and first. If H_EQUAL is false, report the obstruction.
Author's prior is stated to lose: priors are 0-for-5 this arc.

VALIDATION GATE (must pass or the run is void):
  q=3,k=2: S_2 = 10/21 = 0.4761904762 and M_2(4) = -0.2380952381 (R76 sec.3 table).

DECISION RULES (pre-committed):
  H_CONS   CONFIRMED iff max |sum_j M| < 1e-10 for every tested (q,k,eta_0).
  H_EQUAL  CONFIRMED iff spread over j of M_{k+1}(1+j q^k) < 1e-10, all tested q.
           REFUTED if any spread > 1e-8.
  Anything between 1e-10 and 1e-8 -> INCONCLUSIVE, report as such, no verdict.

Not at stake: THEOREM_C_745, Th 78.1-78.3, R81b, eps_k. This probe touches only
the standalone qx+1 paper's remaining rigor step.

Reuses: probe_5_universal_rate.X_gen chain construction (general q),
        bilinear_pair_operator.compute_M_n structure (generalized off q=3).
"""
import sys
from math import gcd
import numpy as np

LOG = []


def log(m=""):
    print(m)
    LOG.append(str(m))


def order_of_two(N):
    """Multiplicative order of 2 mod N."""
    o, x = 1, 2 % N
    while x != 1:
        x = (x * 2) % N
        o += 1
        if o > N:
            raise RuntimeError("2 not invertible")
    return o


def stationary(q, k):
    """Stationary pi of the qx+1 Syracuse chain on residues coprime to q mod q^k.
    Mirrors probe_5_universal_rate.X_gen exactly (float power iteration)."""
    N = q ** k
    inv2 = pow(2, -1, N)
    M = order_of_two(N)
    cp = np.array([r for r in range(N) if gcd(r, q) == 1], dtype=np.int64)
    n = len(cp)
    Z = (2 ** M - 1) / 2 ** M
    inv_idx = np.full(N, -1, dtype=np.int64)
    inv_idx[cp] = np.arange(n)
    base_t = (q * cp + 1) % N
    rows_l, cols_l, vals_l = [], [], []
    inv2v = 1
    src = np.arange(n)
    for v in range(1, M + 1):
        inv2v = (inv2v * inv2) % N
        t = (base_t * inv2v) % N
        rows_l.append(src)
        cols_l.append(inv_idx[t])
        vals_l.append(np.full(n, (0.5 ** v) / Z))
    import scipy.sparse as sp
    K = sp.csr_matrix(
        (np.concatenate(vals_l), (np.concatenate(rows_l), np.concatenate(cols_l))),
        shape=(n, n),
    )
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
    return pi, cp, N


def M_all(pi, cp, N):
    """M(eta) = sum_{xi coprime} mu_hat(xi) * conj(mu_hat(xi*eta)) for eta coprime.
    mu_hat(xi) = sum_r pi(r) e^{-2pi i r xi/N} = DFT of zero-padded pi."""
    dense = np.zeros(N, dtype=np.float64)
    dense[cp] = pi
    mh = np.fft.fft(dense)  # mh[xi] = sum_r dense[r] e^{-2pi i r xi/N}
    out = {}
    mh_cp = mh[cp]
    for eta in cp:
        idx = (cp * int(eta)) % N
        out[int(eta)] = complex(np.sum(mh_cp * np.conj(mh[idx])))
    return out, mh


def main():
    log("# PROBE 6 -- does R76 conservation/leading-mode generalize to (Z/q^k)*?")
    log("# Pre-registered: H_CONS (verbatim port) / H_LEAD_BREAKS / H_EQUAL (falsifier)")
    log("")

    # ---------------- VALIDATION GATE ----------------
    log("## VALIDATION GATE -- q=3,k=2 vs R76 sec.3 table")
    pi, cp, N = stationary(3, 2)
    Mv, _ = M_all(pi, cp, N)
    S2 = Mv[1].real
    M4 = Mv[4].real
    ok1 = abs(S2 - 10 / 21) < 1e-9
    ok2 = abs(M4 - (-0.2380952381)) < 1e-9
    log(f"   M_2(1)=S_2 : got {S2:.10f}  target {10/21:.10f}  {'OK' if ok1 else 'FAIL'}")
    log(f"   M_2(1+3)=M_2(4): got {M4:.10f}  target -0.2380952381  {'OK' if ok2 else 'FAIL'}")
    log(f"   S_2 = -2*M_2(4)? {S2:.10f} vs {-2*M4:.10f}")
    if not (ok1 and ok2):
        log("   GATE: FAIL -> run VOID, no verdict.")
        flush()
        sys.exit(1)
    log("   GATE: PASS")
    log("")

    # ---------------- H_CONS ----------------
    log("## H_CONS -- conservation  sum_{j=0..q-1} M_{k+1}(eta_0 + j*q^k) = 0 ?")
    log("   (eta_0 ranges over (Z/q^k)*; lifts into (Z/q^{k+1})*)")
    cons_worst = 0.0
    for q, k in [(3, 1), (3, 2), (5, 1), (5, 2), (7, 1), (7, 2), (11, 1), (9, 1), (15, 1)]:
        kk = k + 1
        try:
            pi2, cp2, N2 = stationary(q, kk)
        except Exception as e:
            log(f"   q={q} k+1={kk}: skip ({e})")
            continue
        M2, _ = M_all(pi2, cp2, N2)
        base = q ** k
        worst = 0.0
        for eta0 in range(base):
            if gcd(eta0, q) != 1:
                continue
            tot = 0j
            for j in range(q):
                lift = (eta0 + j * base) % N2
                if lift in M2:
                    tot += M2[lift]
            worst = max(worst, abs(tot))
        cons_worst = max(cons_worst, worst)
        log(f"   q={q:2d} k={k} -> k+1={kk} (N={N2:5d}): max|sum_j M| = {worst:.3e}")
    log(f"   => worst over all: {cons_worst:.3e}")
    log(f"   H_CONS: {'CONFIRMED (verbatim port)' if cons_worst < 1e-10 else 'REFUTED/INCONCLUSIVE'}")
    log("")

    # ---------------- lift pairing structure ----------------
    log("## Lift structure at eta_0=1: how many self-inverse / how many pairs?")
    for q in [3, 5, 7, 11]:
        k = 2
        N2 = q ** k
        base = q ** (k - 1)
        lifts = [(1 + j * base) % N2 for j in range(q)]
        selfinv = [h for h in lifts if (h * h) % N2 == 1]
        log(f"   q={q:2d}: {q} lifts of 1 mod {N2}; self-inverse={selfinv}; "
            f"=> (q-1)/2 = {(q-1)//2} mutual-inverse pair(s)")
    log("   (q=3 is the ONLY odd q with exactly ONE pair -> 1 eqn / 1 unknown.)")
    log("")

    # ---------------- H_EQUAL (the falsifier / hopeful branch) ----------------
    log("## H_EQUAL -- is M_{k+1}(1 + j*q^k) independent of j in {1..q-1}?")
    log("   If YES: S_{k+1} = -(q-1)*M_{k+1}(1+q^k)  [76.3's '-2' is really -(q-1)] -> port SAVED.")
    log("   If NO : conservation is 1 eqn in (q-1)/2 >= 2 unknowns -> port BLOCKED.")
    verdicts = {}
    for q, k in [(3, 1), (3, 2), (5, 1), (5, 2), (7, 1), (7, 2), (11, 1)]:
        kk = k + 1
        try:
            pi2, cp2, N2 = stationary(q, kk)
        except Exception as e:
            log(f"   q={q}: skip ({e})")
            continue
        M2, _ = M_all(pi2, cp2, N2)
        base = q ** k
        S = M2[1].real
        vals = []
        for j in range(1, q):
            lift = (1 + j * base) % N2
            vals.append(M2[lift].real if lift in M2 else float("nan"))
        spread = max(vals) - min(vals)
        pred = -(q - 1) * vals[0]
        verdicts[(q, k)] = spread
        log(f"   q={q:2d} k={k}: M(1+j*q^k) over j=1..{q-1} =")
        log(f"          {['%+.8f' % v for v in vals]}")
        log(f"          spread={spread:.3e}   S={S:+.8f}   -(q-1)*M(1+q^k)={pred:+.8f}"
            f"   match={'YES' if abs(S-pred) < 1e-9 else 'NO'}")
    log("")
    worst_spread = max(v for v in verdicts.values())
    log(f"   worst spread over all tested (q,k): {worst_spread:.3e}")
    if worst_spread < 1e-10:
        log("   H_EQUAL: CONFIRMED -> leading-mode identity GENERALIZES as S=-(q-1)*M.")
    elif worst_spread > 1e-8:
        log("   H_EQUAL: REFUTED -> q=3 special; 76.3 does NOT port. Obstruction is REAL.")
    else:
        log("   H_EQUAL: INCONCLUSIVE (between thresholds) -- no verdict per pre-reg.")
    flush()


def flush():
    with open("result_6_conservation_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
