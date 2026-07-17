"""
GATES G0b + G0c (L3 Phase 0 closure) -- weld the object across coordinates by MEASUREMENT.

L3_DEFINITIONS v2: one r_q, three coordinates.
  (M) build_M pair operator: r_q = |lambda_2|/lambda_1  (modal).
  (c) c_k = 3^k ||pi_k - lift(pi_{k-1})||^2 = 3^k(||pi_k||^2 - (1/q)||pi_{k-1}||^2);  r_q = approach
      rate of c_k -> c_inf.  Dictionary: c_k = (3/q)^k S_0(k), S_0(k)=X_k-X_{k-1}, X_k=q^k||pi_k||^2.
  (A) A(z)=Sum S_0(i) z^i -- shares the dictionary with (c); not separately fired (per spec: G0c's
      real content is (M) vs (c)).

G0b: c_k approach rate at HIGH k (v<=64, R27 method). PRE-REG -> 0.62 (q=5), 0.39 (q=7); q=3 -> 1.
     approach rate rho_c(k) = (c_{k+1}-c_k)/(c_k-c_{k-1}) -> r_q  (c_k = c_inf + A r_q^k + ...).
G0c: EQUIVALENCE. (M) modal r_q vs (c) approach rate, at q=5 AND 7 (STANDING TRAP: q=3 alone is
     blind since 1/q=1/3). Declared-same banned; measured-same within error bars closes Phase 0.

Dictionary consistency (exact-algebra gate): c_k ==? (3/q)^k (X_k - X_{k-1}).
NOT AT STAKE: R1-R46. This MEASURES whether the three coordinates are one r_q.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from probe_25_transfer_operator_Aprime import build_M, top_eigs
from probe_6_conservation_generalize import order_of_two

V = 64
LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def K_sparse(q, k, Vt=V):
    N = q ** k
    S = [r for r in range(N) if r % q != 0]
    idx = {r: i for i, r in enumerate(S)}
    n = len(S)
    inv2 = pow(2, -1, N)
    inv2p = [pow(inv2, v, N) for v in range(1, Vt + 1)]
    Z = 1.0 - 2.0 ** (-Vt)
    w = [(2.0 ** (-v)) / Z for v in range(1, Vt + 1)]
    rows, cols, vals = [], [], []
    for r in S:
        i = idx[r]; base = (q * r + 1) % N
        for vi in range(Vt):
            rows.append(i); cols.append(idx[(base * inv2p[vi]) % N]); vals.append(w[vi])
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n)), N


def stationary_sparse(K, tol=1e-14, maxit=100000):
    n = K.shape[0]; pi = np.ones(n) / n; KT = K.T.tocsr()
    for _ in range(maxit):
        nx = KT.dot(pi); nx /= nx.sum()
        if np.abs(nx - pi).sum() < tol:
            return nx
        pi = nx
    return pi


def main():
    log("# GATES G0b + G0c -- weld r_q across coordinates (M),(c). v<=64.")
    log(f"#   PRE-REG (c) approach rate -> 0.62 (q=5), 0.39 (q=7), 1.0 (q=3).")
    log("")

    KMAX = {3: 8, 5: 7, 7: 6}
    rM_by_q = {}
    rc_by_q = {}
    for q in KMAX:
        kmax = KMAX[q]
        log(f"## q={q}  (d=ord_q(2)={order_of_two(q)})")

        # ---- (c) coordinate: high-k c_k and approach rate ----
        p2 = {}
        for k in range(1, kmax + 1):
            K, N = K_sparse(q, k)
            pi = stationary_sparse(K)
            p2[k] = float(pi @ pi)
        X = {k: (q ** k) * p2[k] for k in range(1, kmax + 1)}
        S0 = {k: X[k] - X[k - 1] for k in range(2, kmax + 1)}
        ck = {k: (3 ** k) * (p2[k] - p2[k - 1] / q) for k in range(2, kmax + 1)}
        # dictionary gate
        dict_ok = all(abs(ck[k] - (3.0 / q) ** k * S0[k]) < 1e-9 * max(1, abs(ck[k]))
                      for k in range(2, kmax + 1))
        log(f"   dictionary c_k==(3/q)^k S_0(k): {'OK' if dict_ok else 'FAIL'}")
        log(f"   c_k: " + " ".join(f"{ck[k]:.5f}" for k in range(2, kmax + 1)))
        # approach rate rho_c(k) = (c_{k+1}-c_k)/(c_k-c_{k-1})
        rho = {}
        for k in range(3, kmax):
            dnum = ck[k + 1] - ck[k]; dden = ck[k] - ck[k - 1]
            rho[k] = dnum / dden if abs(dden) > 1e-18 else float('nan')
        log(f"   approach rho_c(k)=(dc_{{k+1}})/(dc_k): " + " ".join(f"k{k}:{rho[k]:+.4f}" for k in rho))
        # settled estimate = last value(s)
        rc = rho[max(rho)] if rho else float('nan')
        rc_settled = np.mean([rho[k] for k in sorted(rho)[-2:]]) if len(rho) >= 2 else rc
        rc_by_q[q] = rc_settled
        log(f"   => (c) approach rate settled ~ {rc_settled:+.4f}  (last {rho[max(rho)]:+.4f})")

        # ---- (M) coordinate: build_M modal r_q ----
        mvals = []
        for L in [1, 2]:
            dim = (order_of_two(q) * q ** (L - 1)) ** 2 * q ** L
            if dim > 60000:
                log(f"   (M) L={L}: dim={dim} > 60000 skipped")
                continue
            M, idx, n = build_M(q, L)
            eigs = top_eigs(M)
            rM = abs(eigs[1]) / abs(eigs[0])
            mvals.append((L, rM))
            log(f"   (M) build_M L={L} (dim={n}): |l1|={abs(eigs[0]):.5f} |l2|={abs(eigs[1]):.5f}  r_M={rM:.4f}")
        rM_by_q[q] = mvals[-1][1] if mvals else float('nan')
        log("")

    # ---- G0c verdict ----
    log("## G0c EQUIVALENCE -- (M) modal vs (c) approach rate, q=5 AND 7 (q=3 = boundary, blind-trap):")
    log(f"   {'q':>4} {'(M) r_M':>10} {'(c) approach':>13} {'|diff|':>8} {'agree?':>8}")
    ok = True
    for q in [3, 5, 7]:
        rM = rM_by_q.get(q, float('nan')); rc = rc_by_q.get(q, float('nan'))
        diff = abs(rM - rc)
        # error bar: build_M L=2 modal is coarse; c-rate has small-k transient. Bar ~0.08.
        agree = diff < 0.08
        if q in (5, 7) and not agree:
            ok = False
        log(f"   {q:>4} {rM:>10.4f} {rc:>13.4f} {diff:>8.4f} {('yes' if agree else 'NO'):>8}")
    log("")
    log(f"## G0b: (c) approach rate vs pre-reg 0.62(q5)/0.39(q7): "
        f"q5={rc_by_q.get(5):.3f} q7={rc_by_q.get(7):.3f}")
    log(f"## G0c: (M) vs (c) measured-same at q=5 AND 7: {'PASS -- Phase 0 object welded' if ok else 'NEEDS MORE k / walk-back'}")
    with open("result_G0bc_equivalence_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
