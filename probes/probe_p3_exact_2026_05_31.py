"""
Compute exact p_3(σ) distribution at depth 3.

Extends probe_pm_distribution depth-2 machinery to depth 3.

Conditioning: v_q(D) = 3, i.e., d_0 = d_1 = d_2 = 0, d_3 ≠ 0.
σ_3 = d_3 = δ_1^(3) + δ_2^(2) + δ_3^(1) + δ_4^(0) mod q.

Cases:
  A. κ_1 = 0 (level-1 fully dominant): contribution = (1/3) × p_2 of iid copy
     - δ_1 = 0 globally, so all level-1 terms vanish
     - The conditioning becomes depth-2 starting from level-2
  B. κ_1 ≠ 0 (level-1 sub-dominant): nested cascade with level-2, level-3 alignment

Output: p_3(σ) for σ ∈ (Z/17)*, then |<omega, p_3>|^2 exact rational.
"""
from __future__ import annotations
import sys, time
from fractions import Fraction
sys.stdout.reconfigure(encoding="utf-8")

q = 17
ord_2 = 8
ord_2_q = ord_2 * q              # 136
ord_2_qsq = ord_2 * q * q        # 2312
ord_2_q3 = ord_2 * q * q * q     # 39304
ord_2_q4 = ord_2 * q ** 4        # 668168
qsq = q * q
q3 = q ** 3
q4 = q ** 4

# Precompute power tables
inv2_q   = pow(2, -1, q)
inv2_q3  = pow(2, -1, q3)
inv2_q4  = pow(2, -1, q4)
inv2_qsq = pow(2, -1, qsq)
pow_inv2_q   = [pow(inv2_q, a, q)   for a in range(ord_2)]
pow_inv2_q3  = [pow(inv2_q3, a, q3) for a in range(ord_2_q3)]
pow_inv2_q4  = [pow(inv2_q4, a, q4) for a in range(ord_2_q4)]
pow_inv2_qsq = [pow(inv2_qsq, a, qsq) for a in range(ord_2_qsq)]

def chi_L(x):
    x %= q
    if x == 0: return 0
    return 1 if pow(x, (q-1)//2, q) == 1 else -1

def W_8(r):
    if r == 0:
        return Fraction(1, 2**ord_2 - 1)
    return Fraction(2**(ord_2 - r), 2**ord_2 - 1)

def W_qsq(r):
    """Weight for residue r mod 8q^2 = 2312."""
    if r == 0:
        return Fraction(1, 2**ord_2_qsq - 1)
    return Fraction(2**(ord_2_qsq - r), 2**ord_2_qsq - 1)

def W_q(r):
    """Weight for residue r mod 8q = 136."""
    if r == 0:
        return Fraction(1, 2**ord_2_q - 1)
    return Fraction(2**(ord_2_q - r), 2**ord_2_q - 1)

t_start = time.time()

# === Depth-0 M[s, σ] matrix ===
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

# === Depth-1 numerator/denominator p_1 ===
# p_1(σ) = (1/Z_1) Σ_κ w_κ M[2κ mod q, σ]
print("Building depth-1 p_1...")
K_inner = 12
num_p1 = [Fraction(0)] * q
den_p1 = Fraction(0)
for k in range(-K_inner, K_inner+1):
    wk = Fraction(1, 3) if k == 0 else Fraction(1, 3 * 2**(ord_2 * abs(k)))
    sk = (2 * k) % q
    for sig in range(q):
        num_p1[sig] += wk * M[sk][sig]
    den_p1 += wk * T_depth0[sk]
print(f"  depth-1 done, time={time.time()-t_start:.1f}s")

# === Depth-2 ===
print("Building depth-2 (case A + B)...")
K_outer = 4
A_MAX_a1 = 25

# Case A (κ_1=0 → contribution to depth-2 case): (1/3) * (depth-1 result on iid copy)
# Need numerator/denominator at the depth-2 level
caseA_p2_num = [Fraction(1, 3) * num_p1[sig] for sig in range(q)]
caseA_p2_den = Fraction(1, 3) * den_p1

# Case B (κ_1 ≠ 0): nested level-2 alignment
caseB_p2_num = [Fraction(0)] * q
caseB_p2_den = Fraction(0)
configs_B_depth2 = 0

for kappa_1 in range(-K_outer, K_outer + 1):
    if kappa_1 == 0:
        continue
    if kappa_1 > 0:
        two_neg_8k_q3 = pow_inv2_q3[(ord_2 * kappa_1) % ord_2_q3]
    else:
        two_neg_8k_q3 = pow(2, ord_2 * abs(kappa_1), q3)
    expr_q3 = (two_neg_8k_q3 - 1) % q3

    for a_1 in range(1, A_MAX_a1 + 1):
        w_kappa_a1 = Fraction(1, 2**(2*a_1 + ord_2 * abs(kappa_1)))
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
                        shift_depth2 = (delta_1_2 + delta_2_1) % q
                        weight = w_kappa_a1 * W_qsq(A_2 % ord_2_qsq) * W_qsq(B_2 % ord_2_qsq)
                        for sig in range(q):
                            caseB_p2_num[sig] += weight * M[shift_depth2][sig]
                        caseB_p2_den += weight * T_depth0[shift_depth2]
                        configs_B_depth2 += 1

print(f"  depth-2 done, configs_B={configs_B_depth2}, time={time.time()-t_start:.1f}s")

# Combine depth-2
num_p2 = [caseA_p2_num[sig] + caseB_p2_num[sig] for sig in range(q)]
den_p2 = caseA_p2_den + caseB_p2_den

# Verify c(2) and |<omega, p_2>|^2
p_2 = [num_p2[sig] / den_p2 for sig in range(q)]
c_2_compute = sum(Fraction(chi_L(sig)) * p_2[sig] for sig in range(q))
print(f"  c(2) computed = {float(c_2_compute):.20f}, expected ≈ 0.15310841")

# === DEPTH 3 ===
print("\nBuilding depth-3 (case A + B)...")

# Case A (κ_1=0): contribution = (1/3) * (depth-2 result on iid copy)
caseA_p3_num = [Fraction(1, 3) * num_p2[sig] for sig in range(q)]
caseA_p3_den = Fraction(1, 3) * den_p2

# Case B (κ_1 ≠ 0): nested level-2 + level-3 alignment
# Now we also need δ_1^(3), so δ_1 at q^4 precision
caseB_p3_num = [Fraction(0)] * q
caseB_p3_den = Fraction(0)
configs_B_depth3 = 0

print(f"  K_outer={K_outer}, A_MAX_a1={A_MAX_a1}")

for kappa_1 in range(-K_outer, K_outer + 1):
    if kappa_1 == 0:
        continue
    t_k = time.time()

    # Need delta_1 at q^4 precision now (for δ_1^(3))
    if kappa_1 > 0:
        two_neg_8k_q4 = pow_inv2_q4[(ord_2 * kappa_1) % ord_2_q4]
    else:
        two_neg_8k_q4 = pow(2, ord_2 * abs(kappa_1), q4)
    expr_q4 = (two_neg_8k_q4 - 1) % q4

    for a_1 in range(1, A_MAX_a1 + 1):
        w_kappa_a1 = Fraction(1, 2**(2*a_1 + ord_2 * abs(kappa_1)))
        pow_inv2_a1_q4 = pow_inv2_q4[a_1 % ord_2_q4]
        delta_1_q4 = (pow_inv2_a1_q4 * expr_q4) % q4
        if delta_1_q4 % q != 0:
            continue
        delta_1_div_q = delta_1_q4 // q
        delta_1_1 = delta_1_div_q % q
        delta_1_2 = (delta_1_div_q // q) % q
        delta_1_3 = ((delta_1_div_q // q) // q) % q
        s_constraint_lvl2 = (-delta_1_1) % q  # δ_2^(0) constraint

        for r_a2 in range(ord_2):
            for r_b2 in range(ord_2):
                if (pow_inv2_q[r_a2] - pow_inv2_q[r_b2]) % q != s_constraint_lvl2:
                    continue
                lifts_A2 = list(range(r_a2 if r_a2 > 0 else ord_2, ord_2_qsq + 1, ord_2))
                lifts_B2 = list(range(r_b2 if r_b2 > 0 else ord_2, ord_2_qsq + 1, ord_2))
                for A_2 in lifts_A2:
                    for B_2 in lifts_B2:
                        pow_A2_qsq = pow_inv2_qsq[A_2 % ord_2_qsq]
                        pow_B2_qsq = pow_inv2_qsq[B_2 % ord_2_qsq]
                        delta_2_qsq = (pow_A2_qsq - pow_B2_qsq) % qsq
                        delta_2_0 = delta_2_qsq % q
                        if delta_2_0 != s_constraint_lvl2:
                            continue
                        delta_2_div_q = (delta_2_qsq - delta_2_0) // q
                        delta_2_1 = delta_2_div_q % q
                        delta_2_2 = (delta_2_div_q // q) % q

                        # Level-3 constraint: δ_3^(0) = -(δ_1^(2) + δ_2^(1)) mod q
                        s_constraint_lvl3 = (-(delta_1_2 + delta_2_1)) % q

                        weight_lvl12 = w_kappa_a1 * W_qsq(A_2 % ord_2_qsq) * W_qsq(B_2 % ord_2_qsq)

                        for r_a3 in range(ord_2):
                            for r_b3 in range(ord_2):
                                if (pow_inv2_q[r_a3] - pow_inv2_q[r_b3]) % q != s_constraint_lvl3:
                                    continue
                                lifts_A3 = list(range(r_a3 if r_a3 > 0 else ord_2, ord_2_q + 1, ord_2))
                                lifts_B3 = list(range(r_b3 if r_b3 > 0 else ord_2, ord_2_q + 1, ord_2))
                                for A_3 in lifts_A3:
                                    for B_3 in lifts_B3:
                                        # δ_3 at q² precision
                                        pow_A3_q = pow_inv2_qsq[A_3 % ord_2_qsq]
                                        pow_B3_q = pow_inv2_qsq[B_3 % ord_2_qsq]
                                        delta_3_qsq = (pow_A3_q - pow_B3_q) % qsq
                                        delta_3_0 = delta_3_qsq % q
                                        if delta_3_0 != s_constraint_lvl3:
                                            continue
                                        delta_3_1 = ((delta_3_qsq - delta_3_0) // q) % q

                                        # σ_3 = δ_1^(3) + δ_2^(2) + δ_3^(1) + δ_4^(0) mod q
                                        # δ_4^(0) comes from (A_4, B_4) iid Geom(1/2) mod 8
                                        # = same structure as depth-0 with shift = (δ_1_3 + δ_2_2 + δ_3_1)
                                        shift_depth3 = (delta_1_3 + delta_2_2 + delta_3_1) % q

                                        weight_full = weight_lvl12 * W_q(A_3 % ord_2_q) * W_q(B_3 % ord_2_q)

                                        # Accumulate into M[shift_depth3][σ]
                                        for sig in range(q):
                                            caseB_p3_num[sig] += weight_full * M[shift_depth3][sig]
                                        caseB_p3_den += weight_full * T_depth0[shift_depth3]
                                        configs_B_depth3 += 1
    print(f"    κ_1={kappa_1} done in {time.time()-t_k:.1f}s (total elapsed {time.time()-t_start:.1f}s), configs_B={configs_B_depth3}")

print(f"\nTotal configs at depth 3 case B: {configs_B_depth3}")
print(f"Total compute time: {time.time()-t_start:.1f}s")

# Combine
num_p3 = [caseA_p3_num[sig] + caseB_p3_num[sig] for sig in range(q)]
den_p3 = caseA_p3_den + caseB_p3_den

p_3 = [num_p3[sig] / den_p3 for sig in range(q)]

# Compute c(3) and |<omega, p_3>|^2
c_3_compute = sum(Fraction(chi_L(sig)) * p_3[sig] for sig in range(q))
print(f"\nc(3) computed = {float(c_3_compute):.30f}")
print(f"c(3) expected ≈ 0.15300533169151")

# omega = chi_4
omega_re = {1: 1, 2: -1, 4: 1, 8: -1, 9: -1, 13: 1, 15: -1, 16: 1}
omega_im = {3: 1, 5: 1, 6: -1, 7: -1, 10: -1, 11: -1, 12: 1, 14: 1}

Re = Fraction(0)
Im = Fraction(0)
for sig in range(1, q):
    if sig in omega_re:
        Re += omega_re[sig] * p_3[sig]
    if sig in omega_im:
        Im += omega_im[sig] * p_3[sig]
mag_sq = Re * Re + Im * Im

print(f"\n=== |<omega, p_3>|^2 exact ===")
print(f"  Re<omega, p_3> = {float(Re):+.20f}")
print(f"  Im<omega, p_3> = {float(Im):+.20f}")
print(f"  |<omega, p_3>|^2 = {float(mag_sq):.30f}")

c_inf_str = "0.15298912060588517527891674877413229926086222622334"
from mpmath import mp, mpf
mp.dps = 60
c_inf_ref = mpf(c_inf_str)
mag_sq_mp = mpf(mag_sq.numerator) / mpf(mag_sq.denominator)
diff = mag_sq_mp - c_inf_ref
print(f"\n  c_inf reference = {c_inf_str}")
print(f"  |<omega, p_3>|^2 - c_inf = {float(diff):+.6e}")

# Compare to predicted (30× decay from 2.27e-6)
predicted_if_exact = 2.27e-6 / 30
print(f"  Predicted if 30× decay continues: ≈ {predicted_if_exact:+.1e}")
print(f"  Plateau if coincidence:           ≈ +2.27e-6")

# Save
import json
out = {
    "q": q,
    "p_3": {str(sig): [int(p_3[sig].numerator), int(p_3[sig].denominator)] for sig in range(q)},
    "c_3": [int(c_3_compute.numerator), int(c_3_compute.denominator)],
    "omega_inner_sq": [int(mag_sq.numerator), int(mag_sq.denominator)],
    "Re_omega": [int(Re.numerator), int(Re.denominator)],
    "Im_omega": [int(Im.numerator), int(Im.denominator)],
    "diff_from_c_inf": float(diff),
}
with open("C:/Collatz/p3_distribution_2026_05_31.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nSaved to C:/Collatz/p3_distribution_2026_05_31.json")
