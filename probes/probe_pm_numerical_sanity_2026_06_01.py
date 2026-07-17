"""
Compute p_m via numerical offset_distribution for m=0,1,2,3,4
Verify against Fraction-based exact p_0, p_1, p_2 from earlier.

This pinpoints whether:
  - The numerical method is correct (and the Fraction p_2 was truncation-noisy)
  - OR there's an error in my numerical p_3
"""
from __future__ import annotations
import sys, gc, time, json
from mpmath import mp, mpf, mpc
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 30
q = 17

def legendre(x, qq):
    x %= qq
    return 0 if x == 0 else (1 if pow(x, (qq-1)//2, qq) == 1 else -1)

chi_L_table = [legendre(x, q) for x in range(q)]
omega_re = {1: 1, 2: -1, 4: 1, 8: -1, 9: -1, 13: 1, 15: -1, 16: 1}
omega_im = {3: 1, 5: 1, 6: -1, 7: -1, 10: -1, 11: -1, 12: 1, 14: 1}

def offset_distribution_mp(q, n, A_MAX):
    N = q ** n
    inv2 = pow(2, -1, N)
    P_U = [mpf(0)] * N
    p_val = inv2
    half = mpf(1) / mpf(2)
    cur_w = half
    for a in range(1, A_MAX + 1):
        P_U[p_val] += cur_w
        p_val = (p_val * inv2) % N
        cur_w *= half
    total = sum(P_U[i] for i in range(N))
    P_U = [x / total for x in P_U]
    u_sup = [i for i in range(N) if P_U[i] != 0]
    u_wt  = [P_U[i] for i in u_sup]
    P_S = P_U
    for j in range(n - 1, 0, -1):
        P_V = [mpf(0)] * N
        for i in range(N):
            P_V[(1 + q * i) % N] += P_S[i]
        P_new = [mpf(0)] * N
        for u, w in zip(u_sup, u_wt):
            for i in range(N):
                P_new[(u * i) % N] += w * P_V[i]
        P_S = P_new
        del P_V
        gc.collect()
    return P_S

def compute_pm(m, A_MAX=80):
    """Compute p_m(σ) and return (c_m, |<omega, p_m>|^2, p_m_vec)."""
    n = m + 1
    N = q ** n
    qm = q ** m
    t0 = time.time()
    P_X = offset_distribution_mp(q, n, A_MAX)
    p_sigma = [mpf(0)] * q
    for j in range(1, q):
        d = j * qm
        Pd = mpf(0)
        for y in range(N):
            Pd += P_X[y] * P_X[(y - d) % N]
        p_sigma[j] = Pd
    total = sum(p_sigma)
    p_m = [v / total for v in p_sigma]
    c_m = sum(mpf(chi_L_table[s]) * p_m[s] for s in range(q))
    Re = sum(omega_re[s] * p_m[s] for s in range(1, q) if s in omega_re)
    Im = sum(omega_im[s] * p_m[s] for s in range(1, q) if s in omega_im)
    mag_sq = Re * Re + Im * Im
    t1 = time.time()
    return p_m, c_m, mag_sq, Re, Im, t1 - t0

# Reference values
c0_ref = mpf(19) / mpf(127)
c1_ref = mpf("0.15317823005468567112594736436994376136486793738682")
c2_ref = mpf("0.15324792077908743677169673975108527199790")
c3_ref = mpf("0.15300533169151402670180823886191235877480")
c4_ref = mpf("0.15298870904682116522680466011500498985497")
c5_ref = mpf("0.15298899941352190530916684536500000000000")
c_inf  = mpf("0.15298912060588517527891674877413229926086222622334")

# Fraction-based exact p_2 (from earlier JSON)
with open("C:/Collatz/pm_distributions_2026_05_31.json") as f:
    data_exact = json.load(f)
from fractions import Fraction
p_2_exact = {int(s): Fraction(*v) for s, v in data_exact["p_2_rational"].items()}

print(f"=== Numerical p_m via offset_distribution at dps=30, A_MAX=80 ===\n")

for m in [0, 1, 2, 3, 4]:
    print(f"--- m = {m} ---")
    p_m, c_m, mag_sq, Re, Im, dt = compute_pm(m, A_MAX=80 if m <= 3 else 100)
    ref = [c0_ref, c1_ref, c2_ref, c3_ref, c4_ref][m]
    print(f"  c({m}) numerical = {c_m}")
    print(f"  c({m}) reference = {ref}")
    print(f"  c({m}) diff       = {float(c_m - ref):+.3e}")
    print(f"  |<omega, p_{m}>|^2 = {mag_sq}")
    print(f"  |<omega>|^2 - c_inf = {float(mag_sq - c_inf):+.3e}")
    print(f"  |<omega>|^2 - c({m})  = {float(mag_sq - c_m):+.3e}")
    print(f"  Time: {dt:.1f}s")
    # Compare p_2 with Fraction-based
    if m == 2:
        print(f"  p_2(σ) comparison (numerical vs Fraction-exact):")
        for sig in range(1, q):
            p_frac = float(p_2_exact[sig].numerator) / float(p_2_exact[sig].denominator)
            diff = float(p_m[sig]) - p_frac
            print(f"    σ={sig:2d}: num={float(p_m[sig]):.10f}, frac={p_frac:.10f}, diff={diff:+.2e}")
    print()

print("=== Summary ===")
print("Look for: does |<omega, p_m>|^2 converge to c_inf as m increases?")
print("If YES → conjecture supported")
print("If NO  → conjecture refuted; depth-1 agreement was coincidence")
