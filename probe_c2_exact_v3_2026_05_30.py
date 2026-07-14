"""
probe_c2_exact_v3_2026_05_30.py

FIX: a_1 marginal conditional on κ_1 ≠ 0 is Geom(4), giving absolute weight
     2^(-(2a_1 + 8|κ_1|)) for each (κ_1, a_1) config.

Iterates a_1 ∈ {1..A_MAX_a1} with explicit Geom(4) weight rather than W_136.
"""
from fractions import Fraction
import sys, time
sys.stdout.reconfigure(encoding="utf-8")

q = 17
ord_2 = 8
ord_2_qsq = 8 * q       # 136
ord_2_q3 = 8 * q * q    # 2312
qsq = q * q             # 289
q3 = q ** 3             # 4913

# Power tables
inv2_q = pow(2, -1, q)
pow_inv2_q = [pow(inv2_q, a, q) for a in range(ord_2)]

inv2_qsq = pow(2, -1, qsq)
pow_inv2_qsq = [pow(inv2_qsq, a, qsq) for a in range(ord_2_qsq)]

inv2_q3 = pow(2, -1, q3)
pow_inv2_q3 = [pow(inv2_q3, a, q3) for a in range(ord_2_q3)]

def chi_2(x):
    x %= q
    if x == 0: return 0
    return 1 if pow(x, (q-1)//2, q) == 1 else -1

def W_8(r):
    if r == 0:
        return Fraction(1, 2**ord_2 - 1)
    return Fraction(2**(ord_2 - r), 2**ord_2 - 1)

def W_136(r):
    """Geom(2) marginal at mod 136. For (A_2, B_2) iid Geom(2)."""
    if r == 0:
        return Fraction(1, 2**ord_2_qsq - 1)
    return Fraction(2**(ord_2_qsq - r), 2**ord_2_qsq - 1)

# Depth-0 N(s), T(s)
N = {}; T = {}
for s in range(q):
    Ns = Fraction(0); Ts = Fraction(0)
    for ar in range(ord_2):
        for br in range(ord_2):
            val = (s + pow_inv2_q[ar] - pow_inv2_q[br]) % q
            w = W_8(ar) * W_8(br)
            if val != 0:
                Ts += w
                Ns += w * chi_2(val)
    N[s] = Ns
    T[s] = Ts

c0 = N[0] / T[0]
print(f"c(0) = {c0} = {float(c0):.15f}")

# c(1)
K_inner = 12
num_c1 = Fraction(0); den_c1 = Fraction(0)
for k in range(-K_inner, K_inner+1):
    wk = Fraction(1, 3) if k == 0 else Fraction(1, 3 * 2**(ord_2 * abs(k)))
    sk = (2 * k) % q
    num_c1 += wk * N[sk]
    den_c1 += wk * T[sk]
c1_val = num_c1 / den_c1
print(f"c(1) = {float(c1_val):.15f}")

# Case A: κ_1 = 0 → contribution = (1/3) · c(1) num/den
case_A_num = Fraction(1, 3) * num_c1
case_A_den = Fraction(1, 3) * den_c1

# Case B: κ_1 ≠ 0, level-2 compensates, ABSOLUTE WEIGHT for each config:
#   P = 2^(-(2a_1 + 8|κ_1|)) · W_136(A_2 mod 136) · W_136(B_2 mod 136)
K_outer = 4
A_MAX_a1 = 30  # Geom(4) tail: 2^(-2·30) = 2^(-60), enough precision

print(f"\nCase B with K_outer={K_outer}, A_MAX_a1={A_MAX_a1}...")
case_B_num = Fraction(0)
case_B_den = Fraction(0)

t0 = time.time()
total_configs = 0
for kappa_1 in range(-K_outer, K_outer + 1):
    if kappa_1 == 0:
        continue
    t_k = time.time()

    # Precompute (2^(-8κ_1) - 1) mod q³
    if kappa_1 > 0:
        two_neg_8k_q3 = pow_inv2_q3[8 * kappa_1 % ord_2_q3]
    else:
        two_neg_8k_q3 = pow(2, 8 * abs(kappa_1), q3)
    expr_q3 = (two_neg_8k_q3 - 1) % q3
    assert expr_q3 % q == 0, f"v_q(2^(-8κ)-1) < 1 for κ = {kappa_1}"

    for a_1 in range(1, A_MAX_a1 + 1):
        # Geom(4)-like absolute weight: 2^(-(2a_1 + 8|κ_1|))
        w_kappa_a1 = Fraction(1, 2**(2*a_1 + ord_2*abs(kappa_1)))

        # δ_1 at full q³-precision (use a_1 mod 2312, but a_1 ≤ 30 < 2312)
        pow_inv2_a1_q3 = pow_inv2_q3[a_1 % ord_2_q3]
        delta_1_q3 = (pow_inv2_a1_q3 * expr_q3) % q3
        # Verify v_q ≥ 1
        if delta_1_q3 % q != 0:
            continue
        delta_1_div_q = delta_1_q3 // q
        delta_1_1 = delta_1_div_q % q
        delta_1_2 = (delta_1_div_q // q) % q

        # Constraint on (A_2, B_2): δ_2^(0) ≡ -δ_1^(1) mod q
        s_constraint = (-delta_1_1) % q

        # Iterate (r_a, r_b) ∈ {0..7}² satisfying constraint
        for r_a in range(ord_2):
            for r_b in range(ord_2):
                if (pow_inv2_q[r_a] - pow_inv2_q[r_b]) % q != s_constraint:
                    continue
                # Lifts to mod 136
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
                        case_B_num += weight * N[shift]
                        case_B_den += weight * T[shift]
                        total_configs += 1
    print(f"    κ_1={kappa_1} done in {time.time()-t_k:.1f}s, configs={total_configs}")

print(f"\nTotal configs: {total_configs}, time {time.time()-t0:.1f}s")

# Combine
num_c2 = case_A_num + case_B_num
den_c2 = case_A_den + case_B_den
c2 = num_c2 / den_c2

print(f"\n=== Result ===")
print(f"case_A_num = {float(case_A_num):.10e}, case_A_den = {float(case_A_den):.10e}")
print(f"case_B_num = {float(case_B_num):.10e}, case_B_den = {float(case_B_den):.10e}")
print(f"\nc(2) computed = {float(c2):.15f}")
print(f"c(2) measured = 0.15324792077909  (FFT)")
print(f"diff          = {float(c2) - 0.15324792077908744:.2e}")
print(f"Δ_2 computed  = {float(c2 - c0):.15e}")
print(f"Δ_2 expected  = {0.15324792077908744 - 19/127:.15e}")
