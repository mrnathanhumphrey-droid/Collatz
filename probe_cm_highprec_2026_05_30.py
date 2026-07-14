"""
probe_cm_highprec_2026_05_30.py

Compute c(2), c(3) at mpmath dps=40 via direct offset_distribution.
Then use Δ_1 (EXACT from c(1) rational), Δ_2, Δ_3 at high precision,
plus Δ_4, Δ_5 at float64, to solve damped oscillation recurrence
and extract c_∞ at the highest precision the data allows.

PSLQ on T1 = c_∞ - c(1), T2 = Δ_∞, T3 = c_∞ against expanded L-value basis.
"""
from __future__ import annotations
import sys, gc, time
from mpmath import mp, mpf, mpc, log, sqrt, pi, exp, digamma, pslq, identify
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 40
q = 17

def legendre(x, qq):
    x %= qq
    return 0 if x == 0 else (1 if pow(x, (qq-1)//2, qq) == 1 else -1)

def offset_distribution_mp(q, n, A_MAX):
    """Compute P_X (Syracuse-distributed q-adic at depth n) as list of mpf, length q^n."""
    N = q ** n
    inv2 = pow(2, -1, N)
    # Initial P_U: single-digit Syracuse distribution
    P_U = [mpf(0)] * N
    p_val = inv2
    cur_w = mpf("0.5")
    half = mpf("0.5")
    for a in range(1, A_MAX + 1):
        P_U[p_val] += cur_w
        p_val = (p_val * inv2) % N
        cur_w *= half
    # normalize
    total = mpf(0)
    for x in P_U: total += x
    P_U = [x / total for x in P_U]

    # sparse u_sup
    u_sup = []
    u_wt = []
    for i, w in enumerate(P_U):
        if w != 0:
            u_sup.append(i)
            u_wt.append(w)
    print(f"    |u_sup| = {len(u_sup)}", flush=True)

    P_S = list(P_U)
    arange = list(range(N))

    for j in range(n - 1, 0, -1):
        t_j = time.time()
        # P_V[(1 + q*i) mod N] = P_S[i]
        P_V = [mpf(0)] * N
        for i in range(N):
            v_idx = (1 + q * i) % N
            P_V[v_idx] += P_S[i]
        # P_new[(u*i) mod N] += w * P_V[i]
        P_new = [mpf(0)] * N
        for u, w in zip(u_sup, u_wt):
            for i in range(N):
                idx = (u * i) % N
                P_new[idx] += w * P_V[i]
        P_S = P_new
        del P_V
        gc.collect()
        print(f"    j={j} done in {time.time()-t_j:.1f}s", flush=True)
    return P_S

def c_m_at_n(n, m, A_MAX=150):
    """Compute c(m) at FFT level n at current mpf precision."""
    N = q ** n
    print(f"  computing offset dist at n={n}, N={N}, A_MAX={A_MAX}...", flush=True)
    t0 = time.time()
    P_X = offset_distribution_mp(q, n, A_MAX)
    print(f"  P_X done in {time.time()-t0:.1f}s", flush=True)

    # P_D(d) for d = j*q^m, j=1..q-1 (only the q-1 values we need)
    qm = q ** m
    chi_table = [legendre(x, q) for x in range(q)]

    t1 = time.time()
    num = mpf(0); den = mpf(0)
    for j in range(1, q):
        d = j * qm
        Pd = mpf(0)
        for y in range(N):
            yp = (y - d) % N
            Pd += P_X[y] * P_X[yp]
        num += Pd * mpf(chi_table[j])
        den += Pd
    c_val = num / den
    print(f"  P_D extraction in {time.time()-t1:.1f}s", flush=True)
    return c_val

# === Compute c(2) at n=3 dps=40 ===
print(f"=== Computing c(2) at dps={mp.dps} ===")
t_total = time.time()
c2 = c_m_at_n(n=3, m=2, A_MAX=150)
print(f"\nc(2) = {c2}")
print(f"(time: {time.time()-t_total:.1f}s)\n")

# === Compute c(3) at n=4 dps=40 ===
print(f"=== Computing c(3) at dps={mp.dps} ===")
t_total = time.time()
c3 = c_m_at_n(n=4, m=3, A_MAX=150)
print(f"\nc(3) = {c3}")
print(f"(time: {time.time()-t_total:.1f}s)\n")

# === c(1) EXACT from rational ===
from fractions import Fraction
c1_exact = Fraction(265011804960406635465672455997699, 1730087916969634762193659498034425)
c0_exact = Fraction(19, 127)
c1 = mpf(c1_exact.numerator) / mpf(c1_exact.denominator)
c0 = mpf(c0_exact.numerator) / mpf(c0_exact.denominator)

# === Δ_m ===
delta1 = c1 - c0
delta2 = c2 - c0
delta3 = c3 - c0

print(f"=== Δ_m at high precision ===")
print(f"  Δ_1 = {delta1}")
print(f"  Δ_2 = {delta2}")
print(f"  Δ_3 = {delta3}")

# Use existing float64 c(4), c(5) at limited precision
c4_f64 = mpf("0.1529887090468210")
c5_f64 = mpf("0.1529889994135218")
delta4 = c4_f64 - c0
delta5 = c5_f64 - c0
print(f"  Δ_4 = {delta4}  (float64 ~ 12 digits)")
print(f"  Δ_5 = {delta5}  (float64 ~ 12 digits)")

# === Solve damped oscillation recurrence ===
# u Δ_{m+1} - v Δ_m + K = Δ_{m+2}, m=1, 2, 3
# Use Δ_1, Δ_2, Δ_3 at high precision, Δ_4, Δ_5 at float64.
# Equation A (m=1): u Δ_2 - v Δ_1 + K = Δ_3  (all high precision)
# Equation B (m=2): u Δ_3 - v Δ_2 + K = Δ_4  (mixed)
# Equation C (m=3): u Δ_4 - v Δ_3 + K = Δ_5  (mixed)
#
# Eq A alone gives K = Δ_3 - u Δ_2 + v Δ_1 in terms of (u, v).
# Plug into B, C to get 2 equations in (u, v).

# Linear system in (u, v, K):
import mpmath as mpmath
# Matrix M and rhs:
M = mpmath.matrix([
    [delta2, -delta1, mpf(1)],
    [delta3, -delta2, mpf(1)],
    [delta4, -delta3, mpf(1)],
])
rhs = mpmath.matrix([delta3, delta4, delta5])
sol = mpmath.lu_solve(M, rhs)
u, v, K = sol[0], sol[1], sol[2]

print(f"\n=== Damped osc recurrence fit ===")
print(f"  u = {u}")
print(f"  v = {v}")
print(f"  K = {K}")

# Derived
rho = sqrt(v) if v > 0 else mpf("NaN")
if v > 0:
    cos_th = u / (2 * rho)
    if abs(cos_th) <= 1:
        theta = mpmath.acos(cos_th)
        print(f"  ρ = {rho}")
        print(f"  θ = {theta}")
A = K / (1 - u + v)
print(f"  A = {A}")
print(f"\nc_∞ = c(0) + A = {c0 + A}")

# === Targets ===
T1 = c0 + A - c1   # c_∞ - c(1)
T2 = A             # c_∞ - 19/127 = Δ_∞
T3 = c0 + A        # c_∞

print(f"\n=== Targets ===")
print(f"  T1 = c_∞ - c(1)   = {T1}")
print(f"  T2 = Δ_∞ = c_∞-c(0) = {T2}")
print(f"  T3 = c_∞          = {T3}")

# === Build L-value candidate basis ===
g_17 = 3
g_5 = 2
g_13 = 2

def build_chi(q, g, r, phi_q):
    table = [mpc(0)] * q
    x = 1
    for k in range(phi_q):
        table[x] = exp(2 * pi * mpc(0, 1) * r * k / phi_q)
        x = (x * g) % q
    return table

def L1_chi(chi_table, q):
    total = mpc(0)
    for a in range(1, q):
        total += chi_table[a] * digamma(mpf(a)/mpf(q))
    return -total / mpf(q)

print(f"\n=== Computing L(1, chi) values ===")
chi_4_17 = build_chi(17, g_17, 4, 16)
chi_2_17 = build_chi(17, g_17, 8, 16)
chi_4_5  = build_chi(5, g_5, 1, 4)
chi_4_13 = build_chi(13, g_13, 3, 12)

L_4_17 = L1_chi(chi_4_17, 17)
L_2_17 = L1_chi(chi_2_17, 17)
L_4_5  = L1_chi(chi_4_5, 5)
L_4_13 = L1_chi(chi_4_13, 13)

print(f"  L(1, chi_4_17) = {L_4_17}")
print(f"  L(1, chi_2_17) = {L_2_17}")
print(f"  L(1, chi_4_5)  = {L_4_5}")
print(f"  L(1, chi_4_13) = {L_4_13}")

# === Candidate vector ===
basis = {
    '1': mpf(1),
    'Re L(1,chi_4_17)': mpf(L_4_17.real),
    'Im L(1,chi_4_17)': mpf(L_4_17.imag),
    'L(1,chi_2_17)/2': mpf(L_2_17.real) / 2,  # standard convention
    'Re L(1,chi_4_5)': mpf(L_4_5.real),
    'Im L(1,chi_4_5)': mpf(L_4_5.imag),
    'Re L(1,chi_4_13)': mpf(L_4_13.real),
    'Im L(1,chi_4_13)': mpf(L_4_13.imag),
    'log(4+√17)/√17': log(mpf(4)+sqrt(mpf(17)))/sqrt(mpf(17)),
    '1/√17': mpf(1)/sqrt(mpf(17)),
    '1/√5': mpf(1)/sqrt(mpf(5)),
    'π/√17': pi/sqrt(mpf(17)),
    'log(2)': log(mpf(2)),
    'log(17)': log(mpf(17)),
}

# === PSLQ ===
print(f"\n=== PSLQ search at higher precision ===")
# Determine effective precision
print(f"  c_∞ raw precision estimate: ~{mp.dps - 8} digits (after 3x3 amplification)")
# Use tol matching estimated precision
tol = mpf(10) ** (-(mp.dps - 10))
print(f"  Using tol = 10^-{mp.dps - 10}, maxcoeff = 10^8")
maxcoeff = 10**8

cand_names = list(basis.keys())
cand_vals = list(basis.values())
for target_name, target_val in [('T1 = c_∞ - c(1)', T1), ('T2 = Δ_∞', T2), ('T3 = c_∞', T3)]:
    print(f"\n--- {target_name} = {target_val} ---")
    vec = [target_val] + cand_vals
    rel = pslq(vec, tol=tol, maxcoeff=maxcoeff)
    if rel is None:
        print(f"  No relation found at tol={float(tol):.2e}")
    else:
        terms = []
        if rel[0] != 0:
            terms.append(f"({rel[0]})·{target_name}")
        for c, name in zip(rel[1:], cand_names):
            if c != 0:
                terms.append(f"({c:+d})·{name}")
        print(f"  Relation: {' + '.join(terms)}")
        residual = sum(c * v for c, v in zip(rel, vec))
        print(f"  Residual: {residual}")
