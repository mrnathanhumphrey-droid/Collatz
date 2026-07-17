"""
Probe 85 rung 2 — is the '7' in 7/45 the Eisenstein norm N(m-omega)?  DECISIVE via
the base-m generalization (vary the hidden parameter, the arc's method).

7/45 = (1/3)*S_inf, S_inf = 7/15 (THEOREM_C_745). The Geom(1/2) halving has PGF
G(z)=z/(2-z); generalize the halving base 2 -> m (coprime to 3): G_m(z)=z/(m-z),
Eisenstein norm N(m-omega) = (m-omega)(m-omega_bar) = m^2+m+1  (=7 at m=2).

At m=2 the denominator 5=1+4 is AMBIGUOUS: 4 = class-mass ratio (2/3)^2/(1/3)^2
(R77, m-INDEPENDENT) OR 4 = m^2 (PGF, m-DEPENDENT). These diverge for m!=2:
    E (Eisenstein) : S_inf(m) = (m^2+m+1) / (3*(m^2+1))     [num N(m-omega), den 5->1+m^2]
    C (class-ratio): S_inf(m) = (m^2+m+1) / 15              [num N(m-omega), den fixed 15]
    N (neither)    : the m=2 factorization was coincidence.

Chain (base m): stationary pi_k of  r -> (3r+1)*m^{-v},  v ~ P(v) prop m^{-v}, on
(Z/3^k)*. S_k(m) = sum_{3 nmid xi} |mu_hat_k(xi)|^2.  m=2 MUST give 7/15 (anchor).
Float power-iteration stationary (Phase-0 lesson); Aitken-extrapolate S_inf from k=3,4,5.
Not at stake: THEOREM_C_745 itself (this asks WHERE its 7 and 5 come from, not whether c=7/45).
"""
import sys, cmath, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

LOG = []
def log(m=""):
    print(m); LOG.append(str(m))

def build_chain_base_m(k, m, Vtrunc=400):
    """Transition matrix (float) of r->(3r+1)*m^{-v} on (Z/3^k)*, weight P(v) prop m^{-v}."""
    N = 3 ** k
    invm = pow(m, -1, N)
    powm = [pow(invm, v, N) for v in range(1, Vtrunc + 1)]       # m^{-v} mod N
    w = np.array([m ** (-v) for v in range(1, Vtrunc + 1)], dtype=float)
    w /= w.sum()                                                  # normalized Geom(1/m)-tail
    states = [r for r in range(N) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(states)}
    n = len(states)
    K = np.zeros((n, n), dtype=float)
    for r in states:
        i = idx[r]
        base = (3 * r + 1) % N
        for v in range(Vtrunc):
            tgt = (base * powm[v]) % N
            K[i, idx[tgt]] += w[v]
    return K, states

def stationary_float(K, iters=2000, tol=1e-15):
    n = K.shape[0]
    pi = np.full(n, 1.0 / n)
    for _ in range(iters):
        nx = pi @ K
        nx /= nx.sum()
        if np.abs(nx - pi).sum() < tol:
            return nx
        pi = nx
    return pi

def S_k(k, m):
    """S_k(m) = sum_{xi in Z/3^k, 3 nmid xi} |mu_hat_k(xi)|^2, mu_hat_k(xi)=sum_r pi(r)e(-r xi/3^k)."""
    N = 3 ** k
    K, states = build_chain_base_m(k, m)
    pi = stationary_float(K)
    r = np.array(states, dtype=float)
    xis = np.array([x for x in range(N) if x % 3 != 0], dtype=float)
    # mu_hat[xi] = sum_r pi[r] exp(-2pi i r xi / N)   (vectorized over xi)
    ang = -2.0 * np.pi / N
    # matrix e^{i ang r xi}: shape (len(xis), n)
    E = np.exp(1j * ang * np.outer(xis, r))
    mu = E @ pi
    return float(np.sum(np.abs(mu) ** 2))

def aitken(s2, s1, s0):
    """Aitken Delta^2 on three consecutive S_k (s0<s1<s2 in k): accelerate to S_inf."""
    denom = (s2 - s1) - (s1 - s0)
    if abs(denom) < 1e-18:
        return s2
    return s2 - (s2 - s1) ** 2 / denom

def main():
    log("# PROBE 85 rung 2 — Eisenstein norm test for the '7' (and '5') in 7/45")
    log("# base-m Syracuse chain; S_inf(m) vs E=(m^2+m+1)/(3(m^2+1)) and C=(m^2+m+1)/15")
    log("")
    M_LIST = [2, 4, 5, 7]
    K_LIST = [3, 4, 5]
    rows = {}
    for m in M_LIST:
        Sk = {}
        for k in [2] + K_LIST:
            Sk[k] = S_k(k, m)
        Sinf = aitken(Sk[5], Sk[4], Sk[3])
        rows[m] = (Sk, Sinf)
        Nm = m * m + m + 1
        pred_E = Nm / (3 * (m * m + 1))
        pred_C = Nm / 15
        log(f"## m={m}   N(m-w)=m^2+m+1={Nm}")
        log(f"   S_k: " + "  ".join(f"S_{k}={Sk[k]:.8f}" for k in [2] + K_LIST))
        log(f"   S_inf (Aitken k=3,4,5) = {Sinf:.8f}")
        log(f"     E  (m^2+m+1)/(3(m^2+1)) = {pred_E:.8f}   resid {Sinf-pred_E:+.2e}")
        log(f"     C  (m^2+m+1)/15         = {pred_C:.8f}   resid {Sinf-pred_C:+.2e}")
        # also: back out the denominator D from S_inf = Nm / D  -> D = Nm/S_inf
        D = Nm / Sinf if Sinf else float('nan')
        log(f"     back-out D = N(m-w)/S_inf = {D:.5f}   (E predicts 3(m^2+1)={3*(m*m+1)}, "
            f"C predicts 15)")
        log("")
    # anchor check
    S2 = rows[2][1]
    log(f"ANCHOR m=2: S_inf={S2:.8f} vs 7/15={7/15:.8f}  (resid {S2-7/15:+.2e}); "
        f"c=S_inf/3={S2/3:.8f} vs 7/45={7/45:.8f}")
    log("")
    # verdict: which form does the data pick across m!=2?
    log("## VERDICT")
    e_ok = c_ok = True
    for m in M_LIST:
        Sinf = rows[m][1]; Nm = m*m+m+1
        e_ok = e_ok and abs(Sinf - Nm/(3*(m*m+1))) < 3e-4
        c_ok = c_ok and abs(Sinf - Nm/15) < 3e-4
    if e_ok and not c_ok:
        log("   -> EISENSTEIN (E): S_inf(m)=(m^2+m+1)/(3(m^2+1)). BOTH the 7=N(m-w) AND the")
        log("      5=1+m^2 are PGF/Eisenstein; the m=2 '5=1+4=class-ratio' reading is the")
        log("      coincidence. Your N(2-w)/(3^2(1+4)) decomposition is STRUCTURAL.")
    elif c_ok and not e_ok:
        log("   -> CLASS-RATIO (C): numerator is N(m-w)=m^2+m+1 (the 7 IS Eisenstein) but the")
        log("      denominator 5 is a FIXED class-mass constant, not 1+m^2. Half the conjecture.")
    elif e_ok and c_ok:
        log("   -> UNDECIDED at these m (E and C too close); need larger m spread.")
    else:
        log("   -> NEITHER: the m=2 factorization 7/(9*5) does not generalize; likely")
        log("      numerical coincidence. Report straight.")
    with open(r"C:\Collatz\result_85_rung2_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
    log("")
    log("[wrote] result_85_rung2_log.txt")

if __name__ == "__main__":
    main()
