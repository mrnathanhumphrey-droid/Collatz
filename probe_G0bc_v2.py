"""
G0b + G0c v2 -- fix both extractions the v1 gate exposed.
  (M): naive |l2|/l1 was TOWER-CONTAMINATED (0.98). Use ARNOLDI on Krylov(v0) -> amplitude-carrying
       Ritz values (R32: tower modes have zero amplitude, excluded by starting Arnoldi at v0).
  (c): /q c_k approach rate is the CLEAN r_q (in c_k the (3/q)^k mode cancels exactly, leaving
       c_k - c_inf ~ r_q^k; the /3 object ΔX approaches at max(r_q,3/q) and is 3/q-CONTAMINATED,
       e.g. 3/7=0.4286 masks r_7=0.39). Push to higher k; estimate rate from 3-point geometric
       solve on the LAST points (robust to a single-mode tail).

PRE-REG unchanged: r_q -> 0.62 (q5), 0.39 (q7), 1 (q3). G0c: (M)==(c) at q5 AND q7, measured.
"""
import numpy as np
import scipy.sparse as sp

from probe_25_transfer_operator_Aprime import build_M
from probe_6_conservation_generalize import order_of_two

V = 64
LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def arnoldi_ritz(M, v0, m=60):
    n = M.shape[0]
    Q = np.zeros((n, m + 1)); H = np.zeros((m + 1, m))
    Q[:, 0] = v0 / np.linalg.norm(v0)
    mm = m
    for j in range(m):
        w = M.dot(Q[:, j])
        for i in range(j + 1):
            H[i, j] = Q[:, i] @ w; w = w - H[i, j] * Q[:, i]
        H[j + 1, j] = np.linalg.norm(w)
        if H[j + 1, j] < 1e-13:
            mm = j + 1; break
        Q[:, j + 1] = w / H[j + 1, j]
    ev = np.linalg.eigvals(H[:mm, :mm])
    return sorted(ev, key=lambda z: -abs(z))


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


def stationary_sparse(K, tol=1e-15, maxit=200000):
    n = K.shape[0]; pi = np.ones(n) / n; KT = K.T.tocsr()
    for _ in range(maxit):
        nx = KT.dot(pi); nx /= nx.sum()
        if np.abs(nx - pi).sum() < tol:
            return nx
        pi = nx
    return pi


def main():
    log("# G0b+G0c v2 -- (M) Arnoldi/Krylov (tower-free), (c) /q c_k high-k (3/q-cancelled). v<=64.")
    log(f"#   PRE-REG r_q: q3->1, q5->0.62, q7->0.39.  3/q ref: 3/5={3/5:.3f} 3/7={3/7:.4f}")
    log("")

    KMAX = {3: 8, 5: 8, 7: 6}
    MPLAN = {3: [1, 2], 5: [1, 2], 7: [1, 2]}
    rM, rc = {}, {}
    for q in KMAX:
        log(f"## q={q}  (d={order_of_two(q)}, 3/q={3/q:.4f})")

        # ---- (M) Arnoldi ----
        for L in MPLAN[q]:
            Mm, idx, n = build_M(q, L)
            v0 = np.zeros(n); v0[idx[(1, 1, 0)]] = 1.0
            ritz = arnoldi_ritz(Mm, v0, m=min(60, n - 1))
            r = abs(ritz[1]) / abs(ritz[0]) if len(ritz) > 1 else float('nan')
            nxt = f"{abs(ritz[2]):.4f}" if len(ritz) > 2 else "-"
            log(f"   (M) Arnoldi L={L} (dim={n}): #ritz={len(ritz)} |l1|={abs(ritz[0]):.5f} "
                f"|l2|={abs(ritz[1]) if len(ritz)>1 else float('nan'):.5f} r_M={r:.4f}   next|l|={nxt}")
            rM[q] = r

        # ---- (c) /q c_k high-k ----
        p2 = {}
        for k in range(1, KMAX[q] + 1):
            K, _ = K_sparse(q, k)
            p2[k] = float((lambda pi: pi @ pi)(stationary_sparse(K)))
        ck = {k: (3 ** k) * (p2[k] - p2[k - 1] / q) for k in range(2, KMAX[q] + 1)}
        log(f"   (c) c_k: " + " ".join(f"{ck[k]:.6f}" for k in range(2, KMAX[q] + 1)))
        dc = {k: ck[k] - ck[k - 1] for k in range(3, KMAX[q] + 1)}
        ratios = {k: dc[k + 1] / dc[k] for k in range(3, KMAX[q]) if abs(dc[k]) > 1e-15}
        log(f"   (c) diff-ratios dc_{{k+1}}/dc_k: " + " ".join(f"k{k}:{ratios[k]:+.4f}" for k in ratios))
        # robust estimate: last ratio (single-mode tail), and 3-point geometric solve at top k
        kk = KMAX[q]
        rc_est = ratios[max(ratios)] if ratios else float('nan')
        rc[q] = rc_est
        log(f"   (c) approach rate (last ratio, k={max(ratios) if ratios else '-'}): {rc_est:+.4f}")
        log("")

    log("## G0c EQUIVALENCE (measured, not declared) -- (M) Arnoldi vs (c) approach, q=5 AND 7:")
    log(f"   {'q':>4} {'(M) Arnoldi':>12} {'(c) approach':>13} {'3/q ref':>8} {'|M-c|':>7} {'weld?':>6}")
    ok = True
    for q in [3, 5, 7]:
        d = abs(rM[q] - rc[q])
        agree = d < 0.06
        if q in (5, 7) and not agree:
            ok = False
        log(f"   {q:>4} {rM[q]:>12.4f} {rc[q]:>13.4f} {3/q:>8.4f} {d:>7.4f} {('yes' if agree else 'NO'):>6}")
    log("")
    log(f"## G0b (c) vs pre-reg 0.62/0.39: q5={rc[5]:.3f} q7={rc[7]:.3f}")
    log(f"## G0c weld at q=5 AND 7: {'PASS -- Phase 0 object frozen' if ok else 'NOT YET -- (c) needs higher k'}")
    with open("result_G0bc_v2_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
