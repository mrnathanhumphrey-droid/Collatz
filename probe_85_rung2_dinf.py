"""
Probe 85 rung 2 — the CONVERGENT generalized observable D_inf(p) = lim (3/p)^k S_k.

Axes A/B showed S_k(p) diverges as (p/3)^k off p=3, i.e. ||d_k||^2 ~ 3^{-k} for ALL p.
So the p-independent convergent observable is
    D_inf(p) := lim_k (3/p)^k * S_k(p),   D_inf(3) = 7/15  (anchor).
Question: does the cyclotomic numerator N_p = Phi_p(2) = 2^p - 1 sit cleanly in D_inf(p)?
Prelim (coarse) data hinted back-out D(p)=(2^p-1)/D_inf = 15, 63 = 2^4-1, 2^6-1 -> 2^{p+1}-1,
but p=7 (low k) broke it. Need clean high-k data -> use an EXACT Plancherel identity that
kills the dense-matrix memory wall, and sparse stationary, to push k high over many primes.

Exact identity (no dense mu-hat matrix):
   S_k = sum_{p nmid xi} |mu_hat_k(xi)|^2
       = N * ||pi||^2  -  p^{k-1} * sum_{c in Z/p^{k-1}} ( mass of pi on {r : r = c mod p^{k-1}} )^2
   [from sum_{p nmid xi} e(-(r-r') xi / N) = N*[r=r'] - p^{k-1}*[r = r' mod p^{k-1}] ].
Base-2 halving throughout (Collatz halving); map r->(p r+1) 2^{-v}. Not at stake: THEOREM_C_745.
"""
import sys
import numpy as np
from scipy import sparse
sys.stdout.reconfigure(encoding="utf-8")

LOG = []
def log(m=""):
    print(m); LOG.append(str(m))

def build_sparse(p, k, Vtrunc=64):
    N = p ** k
    inv2 = pow(2, -1, N)
    pw = [pow(inv2, v, N) for v in range(1, Vtrunc + 1)]
    w = np.array([2.0 ** (-v) for v in range(1, Vtrunc + 1)]); w /= w.sum()
    states = [r for r in range(N) if r % p != 0]
    idx = {r: i for i, r in enumerate(states)}
    n = len(states)
    rows = np.empty(n * Vtrunc, dtype=np.int64)
    cols = np.empty(n * Vtrunc, dtype=np.int64)
    vals = np.empty(n * Vtrunc, dtype=np.float64)
    t = 0
    for r in states:
        i = idx[r]; base = (p * r + 1) % N
        for v in range(Vtrunc):
            rows[t] = i; cols[t] = idx[(base * pw[v]) % N]; vals[t] = w[v]; t += 1
    K = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
    return K, np.array(states, dtype=np.int64)

def stationary_sparse(K, iters=6000, tol=1e-15):
    n = K.shape[0]
    pi = np.full(n, 1.0 / n)
    Kt = K.T.tocsr()
    for _ in range(iters):
        nx = Kt @ pi; nx /= nx.sum()
        if np.abs(nx - pi).sum() < tol:
            return nx
        pi = nx
    return pi

def S_k_exact(p, k):
    """Exact-identity S_k (no dense DFT)."""
    N = p ** k
    K, states = build_sparse(p, k)
    pi = stationary_sparse(K)
    term1 = N * float(np.sum(pi ** 2))
    # group pi mass by residue mod p^{k-1}
    Nm1 = p ** (k - 1)
    classes = states % Nm1
    mass = np.zeros(Nm1, dtype=np.float64)
    np.add.at(mass, classes, pi)
    term2 = Nm1 * float(np.sum(mass ** 2))
    return term1 - term2

def aitken(s2, s1, s0):
    d = (s2 - s1) - (s1 - s0)
    return s2 if abs(d) < 1e-18 else s2 - (s2 - s1) ** 2 / d

def main():
    log("# PROBE 85 rung 2 — D_inf(p) = lim (3/p)^k S_k ; is the numerator cyclotomic 2^p-1?")
    log("")
    plan = {3: [6, 7, 8], 5: [4, 5, 6], 7: [3, 4, 5], 11: [3, 4], 13: [3, 4]}
    res = {}
    for p in (3, 5, 7, 11, 13):
        ks = plan[p]
        Sk = {}
        for k in ks:
            Sk[k] = S_k_exact(p, k)
        Dk = {k: (3.0 / p) ** k * Sk[k] for k in ks}
        if len(ks) >= 3:
            Dinf = aitken(Dk[ks[2]], Dk[ks[1]], Dk[ks[0]])
        else:
            Dinf = Dk[ks[-1]]   # 2 points: take highest k (already flat)
        Np = 2 ** p - 1
        res[p] = (Sk, Dk, Dinf, Np)
        log(f"## p={p}   N_p=2^{p}-1={Np}")
        log(f"   S_k:      " + "  ".join(f"S_{k}={Sk[k]:.6g}" for k in ks))
        log(f"   (3/p)^k S_k: " + "  ".join(f"D_{k}={Dk[k]:.8f}" for k in ks))
        log(f"   D_inf = {Dinf:.8f}")
        if abs(Dinf) > 1e-9:
            D = Np / Dinf
            log(f"   back-out D(p) = N_p/D_inf = {D:.5f}   "
                f"[2^(p+1)-1={2**(p+1)-1}, p(p+2)={p*(p+2)}, p^2+p+... ]")
        log("")
    log(f"ANCHOR p=3: D_inf={res[3][2]:.8f} vs 7/15={7/15:.8f} (resid {res[3][2]-7/15:+.1e})")
    log("")
    log("## READOUT — does D(p)=N_p/D_inf follow one clean structural law across all p?")
    log("   If D(p)=2^{p+1}-1 holds (15,63,255,...) -> the '7' AND the denominator are")
    log("   both CYCLOTOMIC (2-power), a genuine structural identity, not p=3 coincidence.")
    for p in (3, 5, 7, 11, 13):
        Dinf = res[p][2]; Np = res[p][3]
        if abs(Dinf) > 1e-9:
            D = Np / Dinf
            log(f"   p={p:2d}: D={D:9.4f}   2^(p+1)-1={2**(p+1)-1:6d}   ratio {D/(2**(p+1)-1):.4f}")
    with open(r"C:\Collatz\result_85_rung2_dinf_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
    log("")
    log("[wrote] result_85_rung2_dinf_log.txt")

if __name__ == "__main__":
    main()
