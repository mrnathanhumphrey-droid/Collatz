"""
GATE G0 (L3 Phase 0) -- FRESH implementation of L3_DEFINITIONS.md sec 1-7, importing NOTHING from
prior probes. Proves the WRITTEN object == the MEASURED object, and reports whether the operator
norm sigma_max(L_k) equals the stationary gain gamma_k (the Phase-2a fork, previewed).

Coded literally from the page:
  sec1 U_k = coprime-to-q residues mod q^k
  sec2 K_k row-stochastic, branch r->(qr+1)2^{-v} mod q^k, weight 2^{-v}/Z_V, V=64; measures act LEFT
  sec3 pi_k stationary (pi K = pi)
  sec4 rho_k reduction, lift_k uniform spread /q
  sec5 W_{k+1} fiberwise mean-zero; P_W subtract fiber mean
  sec6 d_k = pi_k - lift(pi_{k-1})
  sec7 L_k = P_W . (.K_{k+1}) . lift_k ;  sigma_max = ||L_k|| (full SVD) ;  gamma_k = ||L_k d_k||/||d_k||

PASS BAR (pre-registered): (i) gamma_k reproduces R46 sqrt(r_q/3) [q3->0.5774, q5->~.45, q7->~.36];
  (ii) q=3 -> 1/sqrt3 to 1e-10; (iii) REFINE ||L_k d_k - d_{k+1}||/||d_{k+1}|| < 1e-12.
  REPORT: sigma_max vs gamma  (== => operator-norm target valid; > => bound the dominant mode).
"""
import numpy as np

V = 64


def U(q, k):
    N = q ** k
    S = [r for r in range(N) if r % q != 0]
    return S, {r: i for i, r in enumerate(S)}, N


def K_dense(q, k):
    S, idx, N = U(q, k)
    n = len(S)
    inv2 = pow(2, -1, N)
    inv2p = [pow(inv2, v, N) for v in range(1, V + 1)]
    Z = 1.0 - 2.0 ** (-V)
    M = np.zeros((n, n))
    for r in S:
        i = idx[r]; base = (q * r + 1) % N
        for v in range(1, V + 1):
            M[i, idx[(base * inv2p[v - 1]) % N]] += (2.0 ** (-v)) / Z
    return M, S, idx, N


def stationary(K):
    n = K.shape[0]
    pi = np.ones(n) / n
    for _ in range(20000):
        nx = pi @ K
        nx /= nx.sum()
        if np.abs(nx - pi).sum() < 1e-15:
            return nx
        pi = nx
    return pi


def fiber_ids(S, q, k):
    """parent index (r mod q^{k-1}) as a contiguous fiber id per state of U_k."""
    qkm1 = q ** (k - 1)
    par = np.array([r % qkm1 for r in S])
    uniq = {p: i for i, p in enumerate(sorted(set(par)))}
    return np.array([uniq[p] for p in par], dtype=np.int64)


def proj_W(u, fib, q):
    nf = int(fib.max()) + 1
    s = np.zeros(nf); np.add.at(s, fib, u)
    return u - (s / q)[fib]


def helmert(q):
    """orthonormal (q x (q-1)) matrix with mean-zero columns (contrasts)."""
    A = np.zeros((q, q - 1))
    for j in range(1, q):
        col = np.zeros(q)
        col[:j] = 1.0
        col[j] = -j
        A[:, j - 1] = col / np.linalg.norm(col)
    return A


def W_basis(S, idx, q, k):
    """orthonormal columns spanning W_k = fiberwise-mean-zero on U_k. Fibers = level-(k-1) parents."""
    fib = fiber_ids(S, q, k)
    H = helmert(q)
    nf = int(fib.max()) + 1
    n = len(S)
    cols = []
    # children of each fiber, in stable order
    members = [[] for _ in range(nf)]
    for i, f in enumerate(fib):
        members[f].append(i)
    for f in range(nf):
        ch = members[f]                     # length q
        for c in range(q - 1):
            v = np.zeros(n)
            for a, gi in enumerate(ch):
                v[gi] = H[a, c]
            cols.append(v)
    B = np.array(cols).T                     # n x dimW
    return B, fib


def main():
    print("# GATE G0 -- fresh code of L3_DEFINITIONS sec1-7. sigma_max(L_k) & gamma_k, q=3,5,7.")
    print(f"#   1/sqrt3 = {1/np.sqrt(3):.10f}   (q=3 target).  r5~.62->{np.sqrt(0.62/3):.4f}  r7~.39->{np.sqrt(0.39/3):.4f}")
    print()
    KMAX = {3: 6, 5: 4, 7: 3}
    for q in KMAX:
        kmax = KMAX[q]
        print(f"## q={q}")
        S, idx, Kd, pi = {}, {}, {}, {}
        for k in range(1, kmax + 2):                      # need level kmax+1 for L_{kmax}
            Kd[k], S[k], idx[k], _ = K_dense(q, k)
            pi[k] = stationary(Kd[k])
        print(f"   {'k':>3} {'dimW_k':>7} {'sigma_max(L_k)':>15} {'gamma_k(stat)':>14} "
              f"{'sig==gam?':>10} {'REFINE err':>11}")
        for k in range(2, kmax + 1):
            # bases
            Bk, fibk = W_basis(S[k], idx[k], q, k)
            Bk1, fibk1 = W_basis(S[k + 1], idx[k + 1], q, k + 1)
            # lift map U_k -> U_{k+1}: (lift u)(r') = u(r' mod q^k)/q
            qk = q ** k
            liftpar = np.array([idx[k][r % qk] for r in S[k + 1]], dtype=np.int64)
            KT = Kd[k + 1].T
            # build M = B_{k+1}^T @ K_{k+1}^T @ lift @ B_k   (dimW_{k+1} x dimW_k)
            dw = Bk.shape[1]
            Mcols = np.empty((Bk1.shape[1], dw))
            for j in range(dw):
                u = Bk[:, j]
                lifted = u[liftpar] / q
                transp = KT @ lifted
                Mcols[:, j] = Bk1.T @ transp
            sv = np.linalg.svd(Mcols, compute_uv=False)
            sigma = sv[0]
            # stationary gain gamma_k and REFINE
            par_k = fiber_ids(S[k], q, k)
            dk = pi[k] - pi[k - 1][np.array([idx[k - 1][r % (q ** (k - 1))] for r in S[k]])] / q
            lifted_d = dk[liftpar] / q
            Ldk = proj_W(KT @ lifted_d, fibk1, q)
            gamma = np.linalg.norm(Ldk) / np.linalg.norm(dk)
            dk1 = pi[k + 1] - pi[k][np.array([idx[k][r % qk] for r in S[k + 1]])] / q
            refine = np.linalg.norm(Ldk - dk1) / np.linalg.norm(dk1)
            eq = "yes" if abs(sigma - gamma) < 1e-9 else f"NO d={sigma-gamma:+.2e}"
            print(f"   {k:>3} {dw:>7} {sigma:>15.10f} {gamma:>14.10f} {eq:>10} {refine:>11.2e}")
        print()
    print("## READ: (i) gamma_k vs R46 sqrt(r_q/3); (ii) q=3 -> 1/sqrt3; (iii) REFINE<1e-12; ")
    print("##       sigma_max==gamma? -> operator-norm target valid; sigma>gamma -> dominant-mode target.")


if __name__ == "__main__":
    main()
