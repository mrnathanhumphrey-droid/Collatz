"""
Path (a.1): re-derive L exactly.

Compute the full distribution p_m(σ) := P(σ_m = σ | v_q(D) = m) at depths m=0,1,2
using the SAME nested character-sum machinery that computed c(0), c(1), c(2) exactly.

Replace χ_2(shift) accumulation with indicator 1[shift = σ] accumulation per σ.

Output: 17-vector p_m(σ) for σ ∈ {1..16} (and σ=0 must have weight 0 due to v_q(D)=m
conditioning).

Then:
  - Verify c(m) = ⟨χ_2, p_m⟩ matches reference c(m) values
  - Fit p_m = p_∞ + B(σ)·z^m·v(σ) + B(σ)·z̄^m·v̄(σ) using known z ≈ 0.034 + 0.068i
  - Solve for p_∞ via linear system over m=0,1,2
  - PSLQ each component of p_∞ against algebraic / cyclotomic basis
"""
from __future__ import annotations
import sys, time
from fractions import Fraction
sys.stdout.reconfigure(encoding="utf-8")

q = 17
ord_2 = 8
ord_2_qsq = ord_2 * q          # 136
ord_2_q3  = ord_2 * q * q      # 2312
qsq = q * q
q3  = q ** 3

# Power tables
inv2_q   = pow(2, -1, q)
inv2_qsq = pow(2, -1, qsq)
inv2_q3  = pow(2, -1, q3)
pow_inv2_q   = [pow(inv2_q, a, q)   for a in range(ord_2)]
pow_inv2_qsq = [pow(inv2_qsq, a, qsq) for a in range(ord_2_qsq)]
pow_inv2_q3  = [pow(inv2_q3, a, q3)   for a in range(ord_2_q3)]

def chi_2(x):
    x %= q
    if x == 0: return 0
    return 1 if pow(x, (q-1)//2, q) == 1 else -1

def W_8(r):
    if r == 0:
        return Fraction(1, 2**ord_2 - 1)
    return Fraction(2**(ord_2 - r), 2**ord_2 - 1)

def W_136(r):
    if r == 0:
        return Fraction(1, 2**ord_2_qsq - 1)
    return Fraction(2**(ord_2_qsq - r), 2**ord_2_qsq - 1)

# === Depth-0 M[s, σ] matrix ===
# M[s, σ] = Σ_(ar,br) W_8(ar) W_8(br) · 1[(s + 2^(-ar) - 2^(-br)) mod q = σ AND != 0]
# T[s] = Σ_σ M[s, σ] (total non-zero mass at shift s)
print("Building depth-0 M matrix...")
M = [[Fraction(0)]*q for _ in range(q)]
T_depth0 = [Fraction(0)] * q
for s in range(q):
    for ar in range(ord_2):
        for br in range(ord_2):
            val = (s + pow_inv2_q[ar] - pow_inv2_q[br]) % q
            w = W_8(ar) * W_8(br)
            if val != 0:
                M[s][val] += w
                T_depth0[s] += w

# Recover N[s] = Σ_σ χ_2(σ) M[s, σ]
N_depth0 = [sum(chi_2(sig) * M[s][sig] for sig in range(q)) for s in range(q)]
print(f"  c(0) check: N[0]/T[0] = {N_depth0[0]/T_depth0[0]}")
print(f"  Expected: 19/127 = {Fraction(19, 127)}")
print(f"  Match: {N_depth0[0]/T_depth0[0] == Fraction(19, 127)}")

# p_0(σ) = M[0, σ] / T[0]
p_0 = [M[0][sig] / T_depth0[0] for sig in range(q)]
print(f"\n=== p_0(σ) ===")
for sig in range(q):
    if sig == 0:
        continue
    print(f"  p_0({sig:2d}) = {p_0[sig]}  ≈ {float(p_0[sig]):.6f}")

# Verify sum = 1
print(f"  Sum = {sum(p_0)}")
# Verify c(0) = <chi_2, p_0>
c0_check = sum(chi_2(sig) * p_0[sig] for sig in range(q))
print(f"  <chi_2, p_0> = {c0_check} ≈ {float(c0_check):.10f}")
print(f"  Reference c(0) = {Fraction(19, 127)} ≈ {19/127:.10f}")

# === Depth-1: p_1(σ) ===
print("\n=== Depth-1 ===")
# p_1(σ) = (1/Z_1) Σ_κ w_κ · M[2κ mod q, σ]
# Z_1 = Σ_κ w_κ · T[2κ mod q]
# w_κ = 1/3 if κ=0 else 1/(3·2^(8|κ|))
K_inner = 12
num_p1 = [Fraction(0)] * q  # vector of σ
den_p1 = Fraction(0)
for k in range(-K_inner, K_inner+1):
    wk = Fraction(1, 3) if k == 0 else Fraction(1, 3 * 2**(ord_2 * abs(k)))
    sk = (2 * k) % q
    for sig in range(q):
        num_p1[sig] += wk * M[sk][sig]
    den_p1 += wk * T_depth0[sk]
p_1 = [num_p1[sig] / den_p1 for sig in range(q)]
print(f"  Sum p_1 = {sum(p_1)}")
c1_check = sum(chi_2(sig) * p_1[sig] for sig in range(q))
print(f"  <chi_2, p_1> = {float(c1_check):.15f}")
ref_c1_num = 265011804960406635465672455997699
ref_c1_den = 1730087916969634762193659498034425
print(f"  Reference c(1) = {ref_c1_num/ref_c1_den:.15f}")
print(f"  Match diff: {float(c1_check) - ref_c1_num/ref_c1_den:.2e}")

# Show p_1 components
print(f"\n  p_1(σ) for σ in (Z/17)*:")
for sig in range(1, q):
    print(f"    p_1({sig:2d}) = {float(p_1[sig]):.8f}   chi_2={chi_2(sig):+d}")

# === Depth-2: p_2(σ) ===
print("\n=== Depth-2 ===")
K_outer = 4
A_MAX_a1 = 25  # tail 2^(-50), enough for our precision

# Case A: κ_1 = 0 → contribution (1/3) · p_1 / 1 (mass)
# Case B: κ_1 ≠ 0, level-2 compensates
print(f"  K_outer={K_outer}, A_MAX_a1={A_MAX_a1}")
case_A_num = [Fraction(1, 3) * num_p1[sig] for sig in range(q)]
case_A_den = Fraction(1, 3) * den_p1

case_B_num = [Fraction(0) for sig in range(q)]
case_B_den = Fraction(0)

t0 = time.time()
total_configs = 0
for kappa_1 in range(-K_outer, K_outer + 1):
    if kappa_1 == 0:
        continue
    t_k = time.time()
    if kappa_1 > 0:
        two_neg_8k_q3 = pow_inv2_q3[8 * kappa_1 % ord_2_q3]
    else:
        two_neg_8k_q3 = pow(2, 8 * abs(kappa_1), q3)
    expr_q3 = (two_neg_8k_q3 - 1) % q3
    assert expr_q3 % q == 0

    for a_1 in range(1, A_MAX_a1 + 1):
        w_kappa_a1 = Fraction(1, 2**(2*a_1 + ord_2*abs(kappa_1)))
        pow_inv2_a1_q3 = pow_inv2_q3[a_1 % ord_2_q3]
        delta_1_q3 = (pow_inv2_a1_q3 * expr_q3) % q3
        if delta_1_q3 % q != 0:
            continue
        delta_1_div_q = delta_1_q3 // q
        delta_1_1 = delta_1_div_q % q
        delta_1_2 = (delta_1_div_q // q) % q
        s_constraint = (-delta_1_1) % q

        for r_a in range(ord_2):
            for r_b in range(ord_2):
                if (pow_inv2_q[r_a] - pow_inv2_q[r_b]) % q != s_constraint:
                    continue
                lifts_A = list(range(r_a if r_a > 0 else ord_2, ord_2_qsq + 1, ord_2))
                lifts_B = list(range(r_b if r_b > 0 else ord_2, ord_2_qsq + 1, ord_2))
                for A_2 in lifts_A:
                    for B_2 in lifts_B:
                        pow_A2_qsq = pow_inv2_qsq[A_2 % ord_2_qsq]
                        pow_B2_qsq = pow_inv2_qsq[B_2 % ord_2_qsq]
                        delta_2_qsq = (pow_A2_qsq - pow_B2_qsq) % qsq
                        delta_2_0 = delta_2_qsq % q
                        if delta_2_0 != s_constraint:
                            continue
                        delta_2_div_q = (delta_2_qsq - delta_2_0) // q
                        delta_2_1 = delta_2_div_q % q
                        shift = (delta_1_2 + delta_2_1) % q
                        weight = w_kappa_a1 * W_136(A_2 % ord_2_qsq) * W_136(B_2 % ord_2_qsq)
                        for sig in range(q):
                            case_B_num[sig] += weight * M[shift][sig]
                        case_B_den += weight * T_depth0[shift]
                        total_configs += 1
    print(f"    κ_1={kappa_1}: configs total {total_configs}, time {time.time()-t_k:.1f}s")

print(f"  Total configs: {total_configs}, time {time.time()-t0:.1f}s")

num_p2 = [case_A_num[sig] + case_B_num[sig] for sig in range(q)]
den_p2 = case_A_den + case_B_den
p_2 = [num_p2[sig] / den_p2 for sig in range(q)]
c2_check = sum(chi_2(sig) * p_2[sig] for sig in range(q))
print(f"  <chi_2, p_2> = {float(c2_check):.15f}")
print(f"  Reference c(2) = 0.15324792077909")
print(f"  Match diff: {float(c2_check) - 0.15324792077908744:.2e}")

print(f"\n  p_2(σ) for σ in (Z/17)*:")
for sig in range(1, q):
    print(f"    p_2({sig:2d}) = {float(p_2[sig]):.8f}   chi_2={chi_2(sig):+d}")

# === Save p_0, p_1, p_2 to file for downstream PSLQ ===
import json
out = {
    "q": q,
    "p_0": {str(sig): [int(p_0[sig].numerator), int(p_0[sig].denominator)] for sig in range(q)},
    "p_1_approx_K12_K0_4": {str(sig): float(p_1[sig]) for sig in range(q)},
    "p_2_approx_K0_4_A25": {str(sig): float(p_2[sig]) for sig in range(q)},
    "p_1_rational": {str(sig): [int(p_1[sig].numerator), int(p_1[sig].denominator)] for sig in range(q)},
    "p_2_rational": {str(sig): [int(p_2[sig].numerator), int(p_2[sig].denominator)] for sig in range(q)},
}
with open("C:/Collatz/pm_distributions_2026_05_31.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nSaved to C:/Collatz/pm_distributions_2026_05_31.json")

# === Spectral fit: p_m(σ) = p_∞(σ) + ρ^m·[B_c(σ)·cos(mθ+φ_σ)] (real form) ===
# With c(m) damping: ρ ≈ 0.076, θ ≈ arctan(2) ≈ 1.1071
# 3 unknowns per σ: p_∞(σ), B_cos(σ) = B(σ)cos(φ_σ), B_sin(σ) = B(σ)sin(φ_σ)
# p_m(σ) = p_∞(σ) + ρ^m·[B_cos(σ)·cos(mθ) - B_sin(σ)·sin(mθ)]
import numpy as np
rho = 0.076
theta = np.arctan(2)
print(f"\n=== Spectral fit at ρ={rho}, θ={theta:.5f} ===")
print("σ | p_inf(σ) extracted | B_cos | B_sin | (c_inf contribution)")
p_inf_est = []
for sig in range(q):
    if sig == 0:
        p_inf_est.append(0.0)
        continue
    pm = [float(p_0[sig]), float(p_1[sig]), float(p_2[sig])]
    # 3x3 system: pm[m] = p_inf + rho^m (B_cos cos(m theta) - B_sin sin(m theta))
    M_fit = np.array([
        [1.0, 1.0,                     0.0                  ],
        [1.0, rho*np.cos(theta),       -rho*np.sin(theta)   ],
        [1.0, rho**2*np.cos(2*theta),  -rho**2*np.sin(2*theta)],
    ])
    sol = np.linalg.solve(M_fit, pm)
    p_inf, B_cos, B_sin = sol
    p_inf_est.append(p_inf)
    print(f"  {sig:2d} | {p_inf:+.8f} | {B_cos:+.6e} | {B_sin:+.6e}")

# Check: c_inf_predicted = Σ chi_2(σ) p_inf(σ)
c_inf_pred = sum(chi_2(sig) * p_inf_est[sig] for sig in range(q))
print(f"\nc_inf predicted by fit: {c_inf_pred:.10f}")
print(f"Reference c_inf:          0.15298912060588...")
print(f"Diff:                    {c_inf_pred - 0.15298912060588517:.2e}")

# QR/NQR sums
QR_list = [1, 2, 4, 8, 16, 15, 13, 9]
NQR_list = [3, 5, 6, 7, 10, 11, 12, 14]
qr_sum = sum(p_inf_est[s] for s in QR_list)
nqr_sum = sum(p_inf_est[s] for s in NQR_list)
print(f"\nQR  sum: {qr_sum:.10f}")
print(f"NQR sum: {nqr_sum:.10f}")
print(f"Diff (QR-NQR) = c_inf: {qr_sum - nqr_sum:.10f}")
