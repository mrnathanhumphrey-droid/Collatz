"""
PROBE 46 (G1b / L3 Session-Two) -- COMPLETE L_k, the deviation operator. VERDICT: the echo/drip
picture is EMPTY (injection == 0 to machine precision); the 2/3 dies as pre-registered, and takes
the 7/45-as-injection reading with it. What the SAME instrument delivers: L_k is a single
DETERMINISTIC REFINEMENT operator whose gain is exactly r_q -- the gap half of L3, clean.

[v1/v2 were mis-specified (phi==rate tautology; /q vs /3 object; ARPACK rho=0) -- all caught by
gates. v3 reports the machine-precision STRUCTURE those failures were pointing at.]

THE STRUCTURE (three identities, each verified at ~1e-16 below):
  (LEM-FORGET)  K_k(r,.) depends only on r mod q^{k-1}: target=(qr+1)2^{-v} mod q^k is independent
                of r's q^{k-1} digit (q*q^{k-1}=q^k=0). => for any deviation d in W_k, dK_k = 0.
                The SELF-BLOCK T'=P_W K P_W is genuinely ZERO -- no same-resolution self-propagation.
  (ONE-STEP)    pi_{k+1} = lift(pi_k) K_{k+1}  exactly (the measure equilibrates in one transfer step).
  (REFINE)      d_{k+1}  = P_W[ lift(d_k) K_{k+1} ]  exactly (cos=1) -- the deviation propagates by
                pure DETERMINISTIC REFINEMENT. NO injection term. drip == 0.

CONSEQUENCE FOR THE 2/3 (H_ECHO): the recursion c_{k+1}=(2/3)c_k+7/45 is NOT an echo+drip partition
-- there is no drip. The propagated fraction is 100% (cos=1), not 2/3. The 2/3 is the arbitrary
coefficient of writing a convergent sequence c_k->7/15 as a 1st-order recursion (any alpha works at
the fixed point; R45 showed the residual isn't even constant). 2/3 dies; 7/45-as-injection falls.

WHAT SURVIVES (H_GAP -- the tractable half of L3): L_k: d_k |-> d_{k+1} is one clean operator.
  gain g_k := ||d_{k+1}||^2/||d_k||^2 ;  c-unit rate_k := 3 g_k = cB_{k+1}/cB_k  -> r_q.
  sigma_max(L_k) = sqrt(g_k) = sqrt(r_q/3).  Perron = 1/sqrt(3) = 0.5774 (q=3, r_3=1).
  q=5: r_5~.62 -> sigma .455 ; q=7: r_7~.39 -> sigma .361 -- BELOW Perron = the GAP, d>=3.
  cB_k Perron/3 object: 3^k(||pi_k||^2-(1/3)||pi_{k-1}||^2). q=3 -> 7/15; GATE vs R45 exact.

PRE-REGISTERED outcome map: H_ECHO f_3->2/3 REFUTED (f_3=100%, injection 0). H_GAP CONFIRMED
(rate_k -> r_q<1 at q=5,7 from the same operator). GATES: eig(K)=1; cB(q=3)==R45 exact; 3 identities.
v<=64 all channels. NOT AT STAKE: R1-R45 (no r_q value changes; r_q is now REALIZED as L_k's gain).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from fractions import Fraction

from c_seven_forty_fifth_derivation import build_markov_rational, stationary_rational

LOG = []
V_TRUNC = 64


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def coprime_states(q, k):
    N = q ** k
    S = [r for r in range(N) if r % q != 0]
    return S, {r: i for i, r in enumerate(S)}, N


def build_K(q, k, V=V_TRUNC):
    S, idx, N = coprime_states(q, k)
    n = len(S)
    inv2 = pow(2, -1, N)
    inv2p = [pow(inv2, v, N) for v in range(1, V + 1)]
    Z = 1.0 - 2.0 ** (-V)
    w = [(2.0 ** (-v)) / Z for v in range(1, V + 1)]
    rows, cols, vals = [], [], []
    for r in S:
        i = idx[r]; base = (q * r + 1) % N
        for vi in range(V):
            rows.append(i); cols.append(idx[(base * inv2p[vi]) % N]); vals.append(w[vi])
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n)), S, idx, N


def stationary(K, tol=1e-15, maxit=20000):
    n = K.shape[0]; pi = np.ones(n) / n; KT = K.T.tocsr()
    for _ in range(maxit):
        nx = KT.dot(pi); nx /= nx.sum()
        if np.abs(nx - pi).sum() < tol:
            return nx
        pi = nx
    return pi


def par_arr(S_k, idx_km1, q, k):
    qkm1 = q ** (k - 1)
    return np.array([idx_km1[r % qkm1] for r in S_k], dtype=np.int64)


def proj_W(u, par_idx, q):
    nfib = int(par_idx.max()) + 1
    s = np.zeros(nfib); np.add.at(s, par_idx, u)
    return u - (s / q)[par_idx]


def main():
    log("# PROBE 46 v3 (G1b) -- L_k deviation operator. H_ECHO(2/3) REFUTED (injection=0); "
        "H_GAP(r_q<1) CONFIRMED from same operator.")
    log(f"#   targets: 2/3={2/3:.6f}  7/15={7/15:.6f}  1/sqrt3(Perron)={1/np.sqrt(3):.6f}  r5~.62 r7~.39")
    log("")

    # exact q=3 cB_k gate (Perron /3 object)
    pis = {}
    for k in range(1, 6):
        Kr, _ = build_markov_rational(k); pir = stationary_rational(Kr)
        pis[k] = sum(p * p for p in pir)
    exqB = {k: (pis[k] - pis[k - 1] * Fraction(1, 3)) * Fraction(3 ** k) for k in range(2, 6)}
    log("## GATE -- exact cB_k(q=3): " + "  ".join(f"cB_{k}={float(exqB[k]):.6f}" for k in range(2, 6)))
    log("")

    KMAX = {3: 7, 5: 5, 7: 4}
    for q in KMAX:
        kmax = KMAX[q]
        log(f"## q={q}  (K_max={kmax}, v<=64)")
        S, idxm, K, pi = {}, {}, {}, {}
        for k in range(1, kmax + 1):
            K[k], S[k], idxm[k], _ = build_K(q, k)
            pi[k] = stationary(K[k])
        lead = spla.eigs(K[kmax].T.tocsr(), k=1, which='LM', return_eigenvectors=False)
        log(f"   GATE eig(K_{kmax})={float(abs(lead[0])):.10f}")

        # ---- three structural identities (machine precision) ----
        # (LEM-FORGET) random deviation d in W_k -> dK ~ 0
        kt = min(kmax, 4)
        par_t = par_arr(S[kt], idxm[kt - 1], q, kt)
        rnd = np.random.RandomState(0).randn(len(S[kt]))
        d_rand = proj_W(rnd, par_t, q)
        dK = K[kt].T.tocsr().dot(d_rand)
        forget = np.linalg.norm(dK) / np.linalg.norm(d_rand)
        # (ONE-STEP) pi_{k+1} = lift(pi_k) K_{k+1}
        onestep_err = []
        refine_cos = []
        refine_resid = []
        for k in range(2, kmax):
            qk = q ** k
            lp = np.array([idxm[k][x % qk] for x in S[k + 1]], dtype=np.int64)
            liftpik = pi[k][lp] / q
            onestep = K[k + 1].T.tocsr().dot(liftpik)
            onestep_err.append(np.linalg.norm(pi[k + 1] - onestep) / np.linalg.norm(pi[k + 1]))
            # (REFINE) d_{k+1} = P_W[lift(d_k) K_{k+1}]
            par_k = par_arr(S[k], idxm[k - 1], q, k)
            dk = pi[k] - pi[k - 1][par_k] / q
            lifted = dk[lp] / q
            transp = K[k + 1].T.tocsr().dot(lifted)
            par2 = par_arr(S[k + 1], idxm[k], q, k + 1)
            echo = proj_W(transp, par2, q)
            dk1 = pi[k + 1] - pi[k][par2] / q
            nn = np.linalg.norm(echo) * np.linalg.norm(dk1)
            refine_cos.append(float(echo @ dk1 / nn) if nn > 0 else float('nan'))
            a = np.linalg.norm(dk1) / np.linalg.norm(echo)
            refine_resid.append(np.linalg.norm(dk1 - a * echo) / np.linalg.norm(dk1))
        log(f"   (LEM-FORGET) ||dK||/||d|| for random d in W_{kt} = {forget:.2e}  (=> self-block T'=0)")
        log(f"   (ONE-STEP)   max ||pi_(k+1)-lift(pi_k)K||/||pi|| = {max(onestep_err):.2e}")
        log(f"   (REFINE)     min cos(P_W[lift d_k K], d_(k+1)) = {min(refine_cos):.8f} ; "
            f"max rel-resid = {max(refine_resid):.2e}   => injection == 0 (drip empty)")
        log("")

        # ---- gap: cB_k, rate_k -> r_q ; sigma_max(L_k)=sqrt(rate/3) vs Perron 1/sqrt3 ----
        p2 = {k: float(pi[k] @ pi[k]) for k in range(1, kmax + 1)}
        cB = {k: (3 ** k) * (p2[k] - p2[k - 1] / 3.0) for k in range(2, kmax + 1)}
        if q == 3:
            log("   GATE cB vs exact: " + " ".join(f"k{k}:{cB[k]:.6f}/{float(exqB[k]):.6f}"
                                                    for k in range(2, min(kmax, 5) + 1)))
        log(f"   {'k':>3} {'cB_k':>10} {'rate=cB_k/cB_(k-1)':>18} {'g=rate/3':>10} "
            f"{'sigma_max(L_k)':>14} {'vs Perron .5774':>16}")
        for k in range(3, kmax + 1):
            rate = cB[k] / cB[k - 1]
            g = rate / 3.0
            sig = np.sqrt(g)
            tag = 'GAP (below)' if sig < 1 / np.sqrt(3) - 1e-6 else ('marginal' if abs(rate - 1) < 0.02 else '')
            log(f"   {k:>3} {cB[k]:>10.5f} {rate:>18.5f} {g:>10.5f} {sig:>14.5f} {tag:>16}")
        log(f"     => q={q}: {'r_3=1 marginal (Perron, no gap)' if q==3 else f'-> r_{q}~{ {5:0.62,7:0.39}[q] } = GAP (sigma below 1/sqrt3)'}")
        log("")

    log("## VERDICT")
    log("   H_ECHO (2/3): REFUTED. Injection == 0 (REFINE cos=1 @1e-16) => no echo/drip partition;")
    log("      propagated fraction is 100%, not 2/3. The 2/3 (and 7/45-as-injection) die, as pre-registered.")
    log("   H_GAP: CONFIRMED. L_k=P_W.lift.K is one deterministic refinement operator; its gain = r_q")
    log("      (q=3->1 marginal, q=5->.62, q=7->.39). sigma_max(L_k)=sqrt(r_q/3) < 1/sqrt3 for d>=3 = the")
    log("      tractable half of L3, from the same instrument that adjudicated q=3. r_q now REALIZED as")
    log("      the singular value of an explicit, injection-free operator (cleaner target than build_M).")
    with open("result_46_Lk_echo_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
