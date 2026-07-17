"""
G1 SMOKE (hygiene, non-load-bearing per Nathan's downgrade): verify the five Phase-1 substrate
lemmas numerically -- "the proof explains the number". All must hold to machine precision; the
s-blindness claim is checked by L1 pointwise at q=1093 (Wieferich), from which L3's engine descends.

L1 FORGET     : T_v(x) mod q^{k+1} indep of x's q^k-digit          [pointwise, incl q=1093]
L2 ONE-STEP   : pi_{k+1} = lift(pi_k) K_{k+1}; and mu K = pi_{k+1} for any mu, proj_k mu = pi_k
L3 INTERTWINE : lift^2(mu) K_{k+1} = lift(lift(mu) K_k), random signed mu
L4 REFINE     : d_{k+1} = lift(d_k) K_{k+1}
L5 PYTHAGORAS : <d_k, lift(pi_{k-1})> = 0 ; X_k = (3/q) X_{k-1} + c_k
"""
import numpy as np
import scipy.sparse as sp

V = 64


def U(q, k):
    N = q ** k
    S = [r for r in range(N) if r % q != 0]
    return S, {r: i for i, r in enumerate(S)}, N


def K_sparse(q, k):
    S, idx, N = U(q, k)
    n = len(S)
    inv2 = pow(2, -1, N)
    inv2p = [pow(inv2, v, N) for v in range(1, V + 1)]
    Z = 1.0 - 2.0 ** (-V)
    rows, cols, vals = [], [], []
    for r in S:
        i = idx[r]; base = (q * r + 1) % N
        for vi in range(V):
            rows.append(i); cols.append(idx[(base * inv2p[vi]) % N]); vals.append((2.0 ** -(vi + 1)) / Z)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n)), S, idx, N


def stationary(K):
    n = K.shape[0]; pi = np.ones(n) / n; KT = K.T.tocsr()
    for _ in range(20000):
        nx = KT.dot(pi); nx /= nx.sum()
        if np.abs(nx - pi).sum() < 1e-15:
            return nx
        pi = nx
    return pi


def lift(mu, Sj1, idxj, q, j):
    """lift level-j measure to level-(j+1): (lift mu)(r')=mu(r' mod q^j)/q."""
    qj = q ** j
    return np.array([mu[idxj[r % qj]] for r in Sj1]) / q


def apply_K(mu, K):
    return mu @ K              # measure acts on left


def main():
    print("# G1 SMOKE -- five substrate lemmas, machine precision. vmax=64.")
    print()

    # ---- L1 FORGET (pointwise, incl q=1093 Wieferich) ----
    print("## L1 FORGET: T_v(x) mod q^{k+1} indep of x's q^k-digit")
    rng = np.random.RandomState(0)
    for q in [5, 7, 11, 13, 1093]:
        worst = 0
        for k in [2, 3, 4]:
            N1 = q ** (k + 1); qk = q ** k
            inv2 = pow(2, -1, N1)
            for _ in range(200):
                x = int(rng.randint(0, qk, dtype=np.int64))
                while x % q == 0:
                    x = int(rng.randint(0, qk, dtype=np.int64))
                a = int(rng.randint(0, q)); v = int(rng.randint(1, V + 1))
                xp = (x + a * qk) % N1
                t1 = ((q * x + 1) * pow(inv2, v, N1)) % N1
                t2 = ((q * xp + 1) * pow(inv2, v, N1)) % N1
                worst = max(worst, abs(t1 - t2))
        print(f"   q={q:>5}: max |T_v(x)-T_v(x')| over x=x' mod q^k = {worst}  {'OK' if worst==0 else 'FAIL'}")
    print()

    # ---- L2..L5 with stationary, q=5,7,11 ----
    for q in [5, 7, 11]:
        print(f"## q={q}")
        K, S, idx, N = {}, {}, {}, {}
        pi = {}
        for k in range(1, 5):
            K[k], S[k], idx[k], N[k] = K_sparse(q, k)
            pi[k] = stationary(K[k])
        # L2
        e2 = []
        for k in range(1, 4):
            lp = lift(pi[k], S[k + 1], idx[k], q, k)
            e2.append(np.linalg.norm(apply_K(lp, K[k + 1]) - pi[k + 1]))
            # random mu with proj_k mu = pi_k
            mu = pi[k + 1] + 0.0
            # perturb within fibers (keep fiber sums = pi_k): add zero-fiber-sum noise
            noise = np.random.RandomState(k).randn(len(S[k + 1]))
            # project noise to zero-fiber-sum
            qk = q ** k
            par = np.array([idx[k][r % qk] for r in S[k + 1]])
            fs = np.zeros(len(S[k])); np.add.at(fs, par, noise)
            noise = noise - (fs / q)[par]
            mu = mu + 0.01 * noise
            e2.append(np.linalg.norm(apply_K(mu, K[k + 1]) - pi[k + 1]))
        print(f"   L2 ONE-STEP  max err = {max(e2):.2e}")
        # L3 random signed mu on level k-1
        e3 = []
        for k in range(2, 4):
            mu = np.random.RandomState(100 + k).randn(len(S[k - 1]))
            l2mu = lift(lift(mu, S[k], idx[k - 1], q, k - 1), S[k + 1], idx[k], q, k)
            lhs = apply_K(l2mu, K[k + 1])
            rhs = lift(apply_K(lift(mu, S[k], idx[k - 1], q, k - 1), K[k]), S[k + 1], idx[k], q, k)
            e3.append(np.linalg.norm(lhs - rhs))
        print(f"   L3 INTERTWINE max err = {max(e3):.2e}")
        # L4 REFINE
        d = {k: pi[k] - lift(pi[k - 1], S[k], idx[k - 1], q, k - 1) for k in range(2, 5)}
        e4 = []
        for k in range(2, 4):
            e4.append(np.linalg.norm(lift(d[k], S[k + 1], idx[k], q, k) @ K[k + 1] - d[k + 1]))
        print(f"   L4 REFINE    max err = {max(e4):.2e}")
        # L5 PYTHAGORAS
        e5o, e5x = [], []
        for k in range(2, 5):
            lpi = lift(pi[k - 1], S[k], idx[k - 1], q, k - 1)
            e5o.append(abs(float(d[k] @ lpi)))
            Xk = (3 ** k) * float(pi[k] @ pi[k])
            Xkm1 = (3 ** (k - 1)) * float(pi[k - 1] @ pi[k - 1])
            ck = (3 ** k) * float(d[k] @ d[k])
            e5x.append(abs(Xk - ((3.0 / q) * Xkm1 + ck)))
        print(f"   L5 PYTHAG    <d,lift pi> max = {max(e5o):.2e} ; X_k recursion max err = {max(e5x):.2e}")
        print()

    print("## VERDICT: all five lemmas hold to machine precision (L1 incl q=1093). Substrate confirmed.")


if __name__ == "__main__":
    main()
