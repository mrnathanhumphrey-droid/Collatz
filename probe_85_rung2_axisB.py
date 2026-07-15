"""
Probe 85 rung 2, AXIS B — vary the MULTIPLIER p (3x+1 -> p*x+1), keep base-2 halving.

Axis A (vary halving base m) broke S_k convergence (S_k diverges for m!=2). N(2-omega)
mixes TWO bases: the halving 2 and the cube-root omega (from the multiplier 3). Axis A
moved the halving and wrecked the 2-adic convergence. Axis B moves the OTHER base:
  map: r -> (p*r+1) * 2^{-v}  on (Z/p^k)*,  v ~ Geom(1/2)  [base-2 halving PRESERVED]
Now omega -> primitive p-th root; the Eisenstein/cyclotomic norm generalizes to
  N_p := prod_{primitive p-th roots w} (2 - w) = Phi_p(2) = 2^p - 1     (=7 at p=3).
Base-2 halving is untouched, so S_k(p) should still converge (rate 1/p).

TEST: does S_inf(p) isolate the cyclotomic numerator 2^p-1?  Back out D(p)=N_p/S_inf and
see if it is structural. p=3 MUST give 7/15 (anchor). p in {3,5,7}, primes != 2.
Not at stake: THEOREM_C_745 (this asks WHERE its 7 comes from).
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

LOG = []
def log(m=""):
    print(m); LOG.append(str(m))

def build_chain(p, k, Vtrunc=64):
    """r -> (p*r+1)*2^{-v} on (Z/p^k)*, weight P(v) prop 2^{-v}. Float transition matrix."""
    N = p ** k
    inv2 = pow(2, -1, N)
    pow2 = [pow(inv2, v, N) for v in range(1, Vtrunc + 1)]     # 2^{-v} mod N
    w = np.array([2.0 ** (-v) for v in range(1, Vtrunc + 1)])
    w /= w.sum()
    states = [r for r in range(N) if r % p != 0]
    idx = {r: i for i, r in enumerate(states)}
    n = len(states)
    K = np.zeros((n, n), dtype=float)
    for r in states:
        i = idx[r]
        base = (p * r + 1) % N
        for v in range(Vtrunc):
            K[i, idx[(base * pow2[v]) % N]] += w[v]
    return K, states

def stationary(K, iters=4000, tol=1e-15):
    n = K.shape[0]
    pi = np.full(n, 1.0 / n)
    for _ in range(iters):
        nx = pi @ K; nx /= nx.sum()
        if np.abs(nx - pi).sum() < tol:
            return nx
        pi = nx
    return pi

def S_k(p, k):
    """S_k(p) = sum_{xi: p nmid xi} |mu_hat_k(xi)|^2."""
    N = p ** k
    K, states = build_chain(p, k)
    pi = stationary(K)
    r = np.array(states, dtype=float)
    xis = np.array([x for x in range(N) if x % p != 0], dtype=float)
    E = np.exp((-2j * np.pi / N) * np.outer(xis, r))
    mu = E @ pi
    return float(np.sum(np.abs(mu) ** 2))

def aitken(s2, s1, s0):
    d = (s2 - s1) - (s1 - s0)
    return s2 if abs(d) < 1e-18 else s2 - (s2 - s1) ** 2 / d

def main():
    log("# PROBE 85 rung 2 AXIS B — vary multiplier p (p*x+1), base-2 halving fixed")
    log("# cyclotomic numerator N_p = Phi_p(2) = 2^p - 1  (=7 at p=3, the '7' in 7/45)")
    log("")
    # per-p k-lists sized to keep the dense mu-hat matrix feasible (n = p^{k-1}(p-1))
    plan = {3: [4, 5, 6], 5: [3, 4, 5], 7: [2, 3, 4]}
    res = {}
    for p in (3, 5, 7):
        Np = 2 ** p - 1
        ks = plan[p]
        Sk = {k: S_k(p, k) for k in ks}
        Sinf = aitken(Sk[ks[2]], Sk[ks[1]], Sk[ks[0]])
        res[p] = (Sk, Sinf, Np)
        log(f"## p={p}   N_p=Phi_p(2)=2^{p}-1={Np}")
        log(f"   S_k: " + "  ".join(f"S_{k}={Sk[k]:.8f}" for k in ks))
        log(f"   S_inf (Aitken) = {Sinf:.8f}")
        if Sinf and abs(Sinf) > 1e-9:
            D = Np / Sinf
            log(f"   back-out D(p) = N_p / S_inf = {D:.6f}")
            # candidate structural D forms to eyeball
            log(f"      compare: p*(p+2)={p*(p+2)}, p^2+? , 3*(...)  |  p=3 gives D=15=3*5")
        log("")
    # anchor
    S3 = res[3][1]
    log(f"ANCHOR p=3: S_inf={S3:.8f} vs 7/15={7/15:.8f} (resid {S3-7/15:+.2e}); "
        f"c=S_inf/3={S3/3:.8f} vs 7/45={7/45:.8f}")
    log("")
    log("## READOUT")
    log("   Does S_inf(p) converge (base-2 halving preserved)?  Does N_p=2^p-1 sit cleanly")
    log("   in the numerator (D(p) structural)?  p=3 -> 7/15 is the anchor; p=5,7 decide")
    log("   whether the '7' is the cyclotomic Phi_p(2) or a p=3 coincidence.")
    with open(r"C:\Collatz\result_85_rung2B_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
    log("")
    log("[wrote] result_85_rung2B_log.txt")

if __name__ == "__main__":
    main()
