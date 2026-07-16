"""
PROBE 25 -- ROUTE 1A': the WEIGHTED CASCADE transfer operator M on state (a,b,gamma).
r_q = lambda_2(M)/lambda_1(M). Gate: sum(M^k v0) == exact ||pi_k||^2 for k<=L.

PRE-REGISTRATION (written before running).
------------------------------------------------------------------
R24 killed the tower-free fold (collisions ARE the tower). The CORRECT fixed operator
keeps the shift as a weighted transition and carries the carry. Build it and read r_q as
an eigenvalue ratio.

THE OPERATOR (depth-L truncation of the carry). Two independent addresses; collision to
depth k <=> sum_{m=1}^k q^{m-1} T_m ≡ 0 mod q^k, T_m = 2^{-S_m} - 2^{-S'_m}, processed
digit by digit: condition_m = (gamma_{m-1}+T_m ≡ 0 mod q), carry gamma_m=(gamma_{m-1}+T_m)/q.
STATE s=(a,b,gamma): a=2^{-S} mod q^L and b=2^{-S'} mod q^L (both in H=<2> mod q^L),
gamma in Z/q^L. ONE transfer step = prepend one coordinate to each address:
  pick multipliers g_a=2^{-Delta}, g_b=2^{-Delta'} in H with weights 2^{-Delta}/Z (folded
  geometric over the period ord_2(q^L)); a'=a*g_a, b'=b*g_b; T=a'-b' mod q^L;
  KEEP the branch iff (gamma+T) ≡ 0 mod q; new gamma'=((gamma+T)//q) mod q^L.
M is the resulting sub-stochastic nonnegative matrix. Start v0 = indicator(a=1,b=1,gamma=0).
  sum(M^k v0) = P(collide to depth k) = ||pi_k||^2   [EXACT for k<=L; approx beyond].

WHY r_q = lambda_2/lambda_1. ||pi_k||^2 = sum(M^k v0) ~ sum_i A_i lambda_i^k. The diagonal
diag = P2^k ~ lambda_1^k (top eigenvalue governs the ||pi_k||^2 decay, ~1/3 for all q by
R5). Then cross(k) = ||pi_k||^2/P2^k - 1 - within ~ (A_1-1-within) + A_2 (lambda_2/lambda_1)^k.
So rho_k = c_{k+1}/c_k -> lambda_2/lambda_1 =: r_q.
  => r_3 = 1  <=>  lambda_1 DEGENERATE (Jordan block -> linear growth, divergence).
  => r_q < 1  <=>  SPECTRAL GAP of M  <=>  convergence. THIS is the Phase-3 theorem.

HYPOTHESES / GATES (pre-committed):
  G_GATE (*** the hard gate ***): sum(M_L^k v0) == float ||pi_k||^2 (from stationary) to
      <1e-9 rel for ALL k<=L, every q. If FALSE the operator is built wrong -> STOP & debug
      (do NOT report any r_q). This is exact by construction if the build is correct.
  G_LAMBDA1: lambda_1(M) ~ 1/3 for every q (R5: ||pi_k||^2 ~ 3^{-k}).
  R_Q (read out, NOT predicted -- priors 0-for-8): r_q = |lambda_2|/lambda_1. Check
      r_3 -> 1 (degenerate top), r_q < 1 for q>=5, and CONVERGENCE in L (L=1,2,3),
      and consistency with R22's float rho_k (q=5 ~0.6+ climbing, q=7 ~0.36-0.45, q=11 ~0.18).

BUDGET: dim = (d q^{L-1})^2 * q^L = d^2 q^{3L-2}. Dense eig for small, sparse eigs(k=6) else.
  L=1: 12/80/63/1100.  L=2: 324/10000/21609/1.46M.  L=3(q=3 only): 8748.
  Cap dim 2e5 -> q=11 only L=1; q=5,7 L<=2; q=3 L<=3. No silent trunc (SAID per q).

NOT AT STAKE: R10-R24, R5's rate, R6, R7, R12, THEOREM_C_745.
"""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from probe_6_conservation_generalize import stationary, order_of_two

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def two_subgroup(qL):
    s, x = [], 1
    while True:
        s.append(x)
        x = (x * 2) % qL
        if x == 1:
            break
    return s


def build_M(q, L):
    qL = q ** L
    inv2 = pow(2, -1, qL)
    sub = two_subgroup(qL)
    ordL = len(sub)                       # d*q^{L-1}
    Z = 1.0 - 2.0 ** (-ordL)
    mult = [(pow(inv2, delta, qL), (2.0 ** (-delta)) / Z) for delta in range(1, ordL + 1)]
    states = [(a, b, g) for a in sub for b in sub for g in range(qL)]
    idx = {s: i for i, s in enumerate(states)}
    n = len(states)
    rows, cols, vals = [], [], []
    for (a, b, g) in states:
        i = idx[(a, b, g)]
        for (ga, wa) in mult:
            ap = (a * ga) % qL
            for (gb, wb) in mult:
                bp = (b * gb) % qL
                T = (ap - bp) % qL
                if (g + T) % q == 0:
                    gp = ((g + T) // q) % qL
                    j = idx[(ap, bp, gp)]
                    rows.append(j); cols.append(i); vals.append(wa * wb)
    M = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    return M, idx, n


def collision_mass(M, idx, k):
    n = M.shape[0]
    v = np.zeros(n)
    v[idx[(1, 1, 0)]] = 1.0
    out = []
    for _ in range(k):
        v = M.dot(v)
        out.append(v.sum())
    return out


def top_eigs(M, howmany=6):
    n = M.shape[0]
    if n <= 400:
        w = np.linalg.eigvals(M.toarray())
    else:
        try:
            w = spla.eigs(M, k=min(howmany, n - 2), which='LM', return_eigenvectors=False)
        except Exception:
            w = np.linalg.eigvals(M.toarray())
    return sorted(w, key=lambda z: -abs(z))


def main():
    log("# PROBE 25 -- ROUTE 1A': weighted cascade operator M; r_q = lambda_2/lambda_1")
    log("# Pre-reg: G_GATE(*** hard ***) / G_LAMBDA1(~1/3) / R_Q read out (r_3=1 degenerate, r_q<1 q>=5)")
    log("")

    PLAN = {3: [1, 2, 3], 5: [1, 2], 7: [1, 2], 11: [1]}
    DIMCAP = 200_000
    results = {}

    for q in [3, 5, 7, 11]:
        d = order_of_two(q)
        # exact-ish ||pi_k||^2 reference via stationary (float, R17 residual <=2.2e-16).
        # ONLY need k <= max(L)+1 for the gate; guard by state count so q=11 never blows up.
        ref = {}
        kmax_ref = max(PLAN[q]) + 1
        for k in range(1, kmax_ref + 1):
            if (q - 1) * q ** (k - 1) > 200_000:
                break
            try:
                pi, cp, _ = stationary(q, k)
                ref[k] = float(np.dot(pi, pi))
            except Exception:
                break
        for L in PLAN[q]:
            dim = (d * q ** (L - 1)) ** 2 * q ** L
            if dim > DIMCAP:
                log(f"## q={q} L={L}: dim={dim} > {DIMCAP} -- SKIPPED (SAID, not silent)")
                continue
            M, idx, n = build_M(q, L)
            mass = collision_mass(M, idx, min(7, 6))
            # gate: mass[k-1] vs ref[k] for k<=L
            log(f"## q={q}  L={L}  (dim={n})")
            log(f"   {'k':>3} {'sum(M^k v0)':>18} {'||pi_k||^2 (ref)':>18} {'|rel|':>10} {'exact k<=L?':>12}")
            gate_ok = True
            for k in range(1, len(mass) + 1):
                if k in ref:
                    rel = abs(mass[k - 1] - ref[k]) / abs(ref[k])
                    exk = "yes" if k <= L else "approx"
                    if k <= L and rel >= 1e-9:
                        gate_ok = False
                    log(f"   {k:>3} {mass[k-1]:>18.12f} {ref[k]:>18.12f} {rel:>10.2e} {exk:>12}")
            log(f"   G_GATE (k<=L): {'PASS' if gate_ok else '*** FAIL -> operator built wrong, no r_q reported ***'}")
            if not gate_ok:
                log("")
                continue
            eigs = top_eigs(M)
            l1 = eigs[0]
            l2 = eigs[1] if len(eigs) > 1 else float('nan')
            rq = abs(l2) / abs(l1)
            results[(q, L)] = (abs(l1), abs(l2), rq)
            log(f"   lambda_1 = {abs(l1):.8f}  (expect ~1/3={1/3:.6f})   "
                f"lambda_2 = {complex(l2):.6f}   |l2| = {abs(l2):.8f}")
            log(f"   => r_q = |lambda_2|/lambda_1 = {rq:.8f}")
            log(f"      top |eigs|: {['%.6f' % abs(z) for z in eigs[:6]]}")
            log("")

    # ---- convergence + phase-boundary read ----
    log("## R_Q -- convergence in L and phase boundary")
    log(f"   {'q':>4} " + " ".join(f"{'L='+str(L):>12}" for L in [1, 2, 3]) + f"{'R22 rho trend':>16}")
    trend = {3: "->1", 5: "~0.6+ rising", 7: "~0.36-0.45", 11: "~0.18"}
    for q in [3, 5, 7, 11]:
        cells = []
        for L in [1, 2, 3]:
            if (q, L) in results:
                cells.append(f"{results[(q,L)][2]:>12.6f}")
            else:
                cells.append(f"{'-':>12}")
        log(f"   {q:>4} " + " ".join(cells) + f"{trend[q]:>16}")
    log("")
    log("## READ: r_3 -> 1 means lambda_1 degenerate (Jordan -> linear divergence).")
    log("   r_q<1 for q>=5 = SPECTRAL GAP of M = the Phase-3 convergence, as an eigenvalue.")
    log("   Convergence of r_q(L) in L validates the depth truncation. G_GATE proves M is")
    log("   the right operator (reproduces ||pi_k||^2 exactly to depth L).")
    flush()


def flush():
    with open("result_25_transfer_Aprime_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
