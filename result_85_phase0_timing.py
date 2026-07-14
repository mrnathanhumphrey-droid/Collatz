"""
Probe 85 Phase 0 — MEASURED dim-scaling of the DWM moment (no target moments computed).
Times the existing G1 (3-alternating) computation at n=3,4,5 to fit the real dim-exponent,
then extrapolates to n=6 at V_MAX=16. Phase-0 feasibility only.
"""
import sys, time, cmath
import numpy as np
sys.path.insert(0, r'C:\Collatz')
from bilinear_pair_operator import build_markov_rational, stationary_rational

def compute_G1(N_LEVEL, V_MAX):
    N = 3**N_LEVEL
    TPI = 2j*cmath.pi/N
    K, coprime = build_markov_rational(N_LEVEL)
    pi_q = stationary_rational(K)
    dim = len(coprime)
    idx = {r: i for i, r in enumerate(coprime)}
    pi_f = np.array([float(p) for p in pi_q])
    inv2 = pow(2, -1, N)
    powinv = [pow(inv2, v, N) for v in range(0, 4*V_MAX+2)]

    def Mt(v, vp, j, b):
        M = np.zeros((dim, dim), complex)
        if v == vp: return M
        xj = (pow(3, 2*j-2, N)*pow(inv2, b, N)) % N
        pd = (powinv[v]-powinv[vp]) % N
        et = v+vp
        for i, xi in enumerate(coprime):
            t = idx.get((xi*powinv[et]) % N, -1)
            if t < 0: continue
            M[i, t] += cmath.exp(-TPI*xi*xj*pd)
        return M

    def Off_mean(j, b):
        M = np.zeros((dim, dim), complex)
        for v in range(1, V_MAX+1):
            for vp in range(1, V_MAX+1):
                if v == vp: continue
                M += 2.0**(-v-vp)*Mt(v, vp, j, b)
        return M

    t0 = time.time()
    Off1 = Off_mean(1, 0)
    se = 0j; tr = 0j; nit = 0
    for v1 in range(1, V_MAX+1):
        for vp1 in range(1, V_MAX+1):
            if v1 == vp1: continue
            w1 = 2.0**(-v1-vp1); b1 = v1+vp1
            X1 = Mt(v1, vp1, 1, 0) - Off1
            Off2 = Off_mean(2, b1)
            for v2 in range(1, V_MAX+1):
                for vp2 in range(1, V_MAX+1):
                    if v2 == vp2: continue
                    w2 = 2.0**(-v2-vp2)
                    X2 = Mt(v2, vp2, 2, b1) - Off2
                    P = X1 @ X2 @ X1
                    se += w1*w2*complex(P.sum())
                    tr += w1*w2*complex(np.einsum('i,ii->', pi_f, P))
                    nit += 1
    return time.time()-t0, nit, dim, se.real

if __name__ == "__main__":
    print("Probe 85 Phase 0 — measured dim-scaling (G1, 3-alternating)\n"+"="*60)
    rows = []
    for n in (3, 4, 5):
        t, nit, dim, val = compute_G1(n, 8)   # V_MAX=8 for the scaling fit
        rows.append((n, dim, t));
        print(f"  n={n} dim={dim:>4} V_MAX=8: {t:7.2f}s  ({nit} iters, G1_se={val:.4e})")
    # dim exponent fit
    import math
    dims = np.array([r[1] for r in rows], float); ts = np.array([r[2] for r in rows], float)
    a, logC = np.polyfit(np.log(dims), np.log(ts), 1)
    print(f"\n  fit: time ~ dim^{a:.2f}  (dim^3 = compute-bound matmul; dim<3 = overhead-bound)")
    # V_MAX check at n=3
    t8,_,_,_ = compute_G1(3, 8); t16,_,_,_ = compute_G1(3, 16)
    print(f"  V_MAX check n=3: t(8)={t8:.3f}s t(16)={t16:.3f}s ratio={t16/t8:.2f} (V_MAX^4 -> 16.0)")
    # extrapolate to n=6, V_MAX=16
    C = math.exp(logC)
    for n_t, dim_t in [(5, 162), (6, 486)]:
        t_v16 = C * dim_t**a * (16/8)**4
        print(f"  EXTRAPOLATED n={n_t} (dim={dim_t}) V_MAX=16: {t_v16:8.1f}s = {t_v16/3600:.2f} h  (G1, 2 reductions)")
    print("\n  full bridge = ~2x (G1+G2) x 2 (extra reductions delta_1/vac_pi) ~ 4x the above.")
    print("  memory: matrices are dim^2*16B -> 3.8 MB at n=6. COMPUTE-bound, not memory-bound.")
