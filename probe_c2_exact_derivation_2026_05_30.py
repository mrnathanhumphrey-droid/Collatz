"""
probe_c2_exact_derivation_2026_05_30.py

Explicit derivation of c(2) via the boundary integral framework.

Structure:
  c(2) = Case_A + Case_B (both as num/den with depth-0 N(s), T(s))

  Case A (κ_1 = 0): level-1 dominant, level-2 has its own sub-dom k_2 ∈ Z.
    Identical to c(1) machinery but with (A_2, B_2) playing the level-1 role
    and (A_3, B_3) playing the level-2 (depth-0) role.
    Contribution = (1/3) × (c(1)-style numerator and denominator).

  Case B (κ_1 ≠ 0): level-1 sub-dom, level-2 compensates.
    δ_1^(1) = 2κ_1 · 2^(-a_1) mod q forces δ_2^(0) = -2κ_1·2^(-a_1) mod q.
    Lift (A_2, B_2) to mod 136 (ord_{q²}(2)). For each lift, compute δ_2^(1)
    at q²-precision. Then σ_2 = δ_1^(2) + δ_2^(1) + δ_3^(0) mod q.

Compare result to FFT-measured c(2) = 0.15324792078...
"""
from fractions import Fraction
import sys
sys.stdout.reconfigure(encoding="utf-8")

q = 17
ord_2 = 8         # ord_q(2)
ord_2_sq = 136    # ord_{q²}(2) = 8 × 17
qsq = q * q       # 289

# pow_inv2 mod q
inv2_q = pow(2, -1, q)
pow_inv2_q = [pow(inv2_q, a, q) for a in range(ord_2)]
print("pow_inv2 mod 17:", pow_inv2_q)

# pow_inv2 mod q²
inv2_qsq = pow(2, -1, qsq)
pow_inv2_qsq = [pow(inv2_qsq, a, qsq) for a in range(ord_2_sq)]

# Geom-conditional weight mod 8
def W_8(r):
    if r == 0:
        return Fraction(1, 2**ord_2 - 1)
    return Fraction(2**(ord_2 - r), 2**ord_2 - 1)

# Geom-conditional weight mod 136 (full, summed over Geom)
def W_136(r):
    if r == 0:
        return Fraction(1, 2**ord_2_sq - 1)
    return Fraction(2**(ord_2_sq - r), 2**ord_2_sq - 1)

def chi_2(x):
    x %= q
    if x == 0: return 0
    return 1 if pow(x, (q-1)//2, q) == 1 else -1

# === Depth-0 N(s), T(s) ===
N = {}
T = {}
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
print(f"\nc(0) = {c0} = {float(c0):.15f}  (expected 19/127)")

# === c(1) via the c(1) script formula (validation) ===
K = 12
num_c1 = Fraction(0); den_c1 = Fraction(0)
for k in range(-K, K+1):
    wk = Fraction(1, 3) if k == 0 else Fraction(1, 3 * 2**(ord_2 * abs(k)))
    sk = (2 * k) % q
    num_c1 += wk * N[sk]
    den_c1 += wk * T[sk]
c1 = num_c1 / den_c1
print(f"c(1) = {float(c1):.15f}  (expected 0.153178230055)")
print(f"      diff = {float(c1) - 0.153178230055:.2e}")

# === c(2) Case A (κ_1 = 0) ===
# Case A: level-1 dominant, level-2 has its own sub-dom k_2 ∈ Z.
# σ_2 = δ_2^(1)·κ_2 + δ_3^(0) (analogous to c(1)'s σ_1).
# Contribution to c(2) num/den = (1/3) × (c(1) num and den), with the (1/3)
# factor from P(κ_1 = 0 | level-1 alignment).

# Actually the structure of Case A is identical to c(1) with (A_2, B_2) -> (A_1, B_1) and
# (A_3, B_3) -> (A_2, B_2). So Case A num/den contribution to c(2) is just c(1) num/den
# weighted by (1/3) for the level-1 dominant prob.

# But Case B also contributes. Total c(2) num = Σ contributions, den = Σ contributions.

# === c(2) Case B (κ_1 ≠ 0, level-2 compensates) ===
# For each κ_1 ∈ {±1, ±2, ..., ±K_outer} and each a_1 mod 8 ∈ {1..8}:
#   - δ_1^(1) = 2·κ_1·pow_inv2_q[a_1 mod 8] mod q
#   - δ_1^(2) = 3·κ_1·pow_inv2_q[a_1 mod 8] mod q
#     (from (2^(-8κ_1) - 1)/q² ≡ 3κ_1 mod q via LTE; needs verification)
#   - Constraint: δ_2^(0) ≡ -δ_1^(1) mod q
#   - For each (A_2 mod 136, B_2 mod 136) with pow_inv2_qsq[A_2] - pow_inv2_qsq[B_2] ≡ -δ_1^(1) mod q:
#     - Compute δ_2 = (pow_inv2_qsq[A_2] - pow_inv2_qsq[B_2]) mod qsq
#     - δ_2^(0) = δ_2 mod q (= -δ_1^(1) by constraint)
#     - δ_2^(1) = ((δ_2 - δ_2^(0)) // q) mod q  [careful with negative integers, mod qsq first]
#     - σ_2 shift = (δ_1^(2) + δ_2^(1)) mod q
#     - weight = (κ_1 Geom-weight) × W_136(a_1) × W_136(A_2 mod 136) × W_136(B_2 mod 136)
#       But a_1 is part of the κ_1 sub-dom weight; need to integrate.

# Let me carefully derive the weight structure.
#
# Marginal of (A_1, B_1, κ_1) conditional on level-1 alignment:
#   A_1 ∈ {a + 8 max(κ, 0)}, B_1 ∈ {a + 8 max(-κ, 0)}, a ≥ 1.
#   Marginal: P(a_1 = a, κ_1 = κ) = 2^(-(2a + 8|κ|))  (for κ ≠ 0)
#                                = 2^(-2a)               (for κ = 0)
#   Total mass: Σ_a 2^(-2a) + 2·Σ_{κ ≥ 1} Σ_a 2^(-(2a+8κ)) = (1/3)·(1 + 2·2^(-8)/(1-2^(-8)))
#
# For Case B with κ_1 ≠ 0 (κ given):
#   P(a_1 = a, κ_1 = κ | level-1 align) ∝ 2^(-(2a + 8|κ|))
#
# Conditional weight as a function of (a_1 mod 8, κ_1):
#   marginal over a_1 in residue class r ∈ {1..8}:
#   P(a_1 ≡ r mod 8, κ_1 = κ) = (Σ_{a ≡ r mod 8, a≥1} 2^(-2a)) · 2^(-8|κ|)

# Σ_{a ≡ r mod 8, a≥1} 2^(-2a) for r ∈ {1..7}: a ∈ {r, r+8, r+16, ...}.
# = Σ_{j≥0} 2^(-2(r + 8j)) = 2^(-2r) · 1/(1 - 2^(-16))
# For r = 8: a ∈ {8, 16, 24, ...}, sum = 2^(-16) · 1/(1 - 2^(-16))

# Let's denote U(r) = Σ_{a ≡ r mod 8, a≥1} 2^(-2a). Then:
# U(r) = 2^(-2r) · 2^16 / (2^16 - 1) for r ∈ {1..7}
# U(0) = U(8) = 2^(-16) · 2^16 / (2^16 - 1) = 1/(2^16 - 1)

# Wait U(8) means a ≡ 0 mod 8: a ∈ {8, 16, 24, ...}. First a = 8. So U(0) ≡ U(8) since we identify a mod 8 = 0 with a ≡ 8 mod 8 by convention. Use 0 here.

def U(r):
    """Σ_{a ≡ r mod 8, a≥1} 2^(-2a) as Fraction."""
    if r == 0:  # a ∈ {8, 16, ...}
        return Fraction(1, 2**16 - 1)
    return Fraction(2**(16 - 2*r), 2**16 - 1)

# Sanity: sum over r should be 1/3
total_U = sum(U(r) for r in range(ord_2))
print(f"\nSanity: Σ_r U(r) = {total_U} (expected 1/3 = {Fraction(1,3)})")
assert total_U == Fraction(1, 3), f"U normalization failed: {total_U}"

# Now the Case B contribution:
# c(2) Case B num = Σ_{κ_1 ≠ 0, a_1 res, A_2 mod 136, B_2 mod 136 | constraint}
#                   2·U(a_1 res) · 2^(-8|κ_1|) · W_136(A_2) · W_136(B_2) · N[shift]
# (factor 2 from ±κ_1 symmetry, sum |κ_1| ≥ 1)
# Similarly for den.

# Wait: for each κ_1 sign separately, we have its own (a_1, A_2, B_2) configurations.
# Let's iterate κ_1 ∈ {±1, ±2, ..., ±K_outer} explicitly.

K_outer = 4  # truncate |κ_1| ≤ K_outer; weight 2^(-8|κ|) drops fast

print(f"\nComputing Case B (κ_1 ≠ 0, level-2 compensates)...")
case_B_num = Fraction(0)
case_B_den = Fraction(0)

# For each κ_1 in nonzero range:
for kappa_1 in range(-K_outer, K_outer + 1):
    if kappa_1 == 0:
        continue
    print(f"  κ_1 = {kappa_1}...", flush=True)
    # κ_1 weight (Geom-like)
    w_kappa = Fraction(1, 2**(ord_2 * abs(kappa_1)))
    # For each a_1 mod 8 residue:
    for a_1_res in range(ord_2):
        u_weight = U(a_1_res)
        # Use a_1_res = a_1 mod 8 → 2^(-a_1) mod q = pow_inv2_q[a_1_res]
        pow_inv2_a1 = pow_inv2_q[a_1_res]
        # δ_1^(1) = 2·κ_1·pow_inv2_a1 mod q
        delta_1_1 = (2 * kappa_1 * pow_inv2_a1) % q
        # Now constraint on (A_2, B_2): δ_2^(0) ≡ -δ_1^(1) mod q
        s_constraint = (-delta_1_1) % q
        # Now compute δ_1^(2) - need to carefully derive.
        # 2^(-8κ_1) - 1 mod qsq. By LTE-like expansion:
        # 2^(-8) mod qsq, then (2^(-8))^κ_1 = 2^(-8κ_1) mod qsq.
        pow_neg8 = pow_inv2_qsq[8 % ord_2_sq]  # 2^(-8) mod qsq
        if kappa_1 > 0:
            pow_neg8_kappa = pow(pow_neg8, kappa_1, qsq)
        else:
            # 2^(8|κ|) - need positive direction. 2^(-8κ) for κ < 0 = 2^(8|κ|).
            pow_neg8_kappa = pow(pow(2, 8, qsq), -kappa_1, qsq)
        delta_1_full = (pow_inv2_qsq[a_1_res] * (pow_neg8_kappa - 1)) % qsq
        # Actually δ_1 = 2^(-A_1) - 2^(-B_1) at full precision.
        # For κ_1 = +1, A_1 = a_1+8: 2^(-A_1) - 2^(-B_1) = 2^(-a_1-8) - 2^(-a_1) = 2^(-a_1)(2^(-8) - 1).
        # For κ_1 = +1: δ_1 = pow_inv2_qsq[a_1_res] · (pow_neg8 - 1) mod qsq
        #   But "pow_inv2_qsq[a_1_res]" treats a_1 = a_1_res as a SPECIFIC value, not averaged.
        # We need 2^(-a_1) mod qsq, which depends on a_1 mod 136 NOT just mod 8.

        # Hmm. The a_1 mod 8 determines 2^(-a_1) mod q (period 8). For mod q², period is 136.
        # So a_1 mod 136 matters at q² level.

        # This breaks our approach. We need to lift a_1 to mod 136 too. Let me handle this.
        pass

# Let me redesign: lift a_1 to mod 136.
print("\n[Reset] Lifting a_1 to mod 136 (the proper q²-precision)...")

case_B_num = Fraction(0)
case_B_den = Fraction(0)

for kappa_1 in range(-K_outer, K_outer + 1):
    if kappa_1 == 0:
        continue
    print(f"  κ_1 = {kappa_1}...", flush=True)
    # κ_1 marginal weight (relative, Geom-like)
    w_kappa = Fraction(1, 2**(ord_2 * abs(kappa_1)))

    # For each a_1 mod 136:
    for a_1_full in range(1, ord_2_sq + 1):  # a_1 ∈ {1, ..., 136}, then periodic
        pow_inv2_a1_qsq = pow_inv2_qsq[a_1_full % ord_2_sq]
        pow_inv2_a1_q = pow_inv2_a1_qsq % q
        # weight on a_1: 2^(-a_1) summed over lifts → use W_136(a_1 mod 136) but since we're enumerating mod 136 explicitly, the weight is 2^(-a_1) for THIS lift, normalized somehow.
        #
        # Actually: marginal P(a_1 = a) = 2^(-a) for a ≥ 1. Integrating into mod-136 classes:
        # For each r ∈ {1..136}: P(a_1 ≡ r mod 136) = 2^(-r) / (1 - 2^(-136)) = W_136(r mod 136).
        #
        # And conditional on (κ_1 ≠ 0, a_1 = a): same.
        #
        # So weight for (a_1 mod 136 = r) configuration:
        w_a1 = W_136(a_1_full % ord_2_sq)

        # Compute δ_1 at qsq-precision for this (κ_1, a_1) config
        if kappa_1 > 0:
            # A_1 = a_1 + 8κ_1, B_1 = a_1
            pow_A1 = pow_inv2_qsq[(a_1_full + 8 * kappa_1) % ord_2_sq]
            pow_B1 = pow_inv2_a1_qsq
        else:
            # A_1 = a_1, B_1 = a_1 + 8|κ_1|
            pow_A1 = pow_inv2_a1_qsq
            pow_B1 = pow_inv2_qsq[(a_1_full + 8 * abs(kappa_1)) % ord_2_sq]
        delta_1_full = (pow_A1 - pow_B1) % qsq  # in {0, ..., qsq-1}
        # δ_1^(0) and δ_1^(1) from delta_1_full
        # δ_1^(0) = delta_1_full mod q (should be 0 by alignment, but let's verify)
        delta_1_0 = delta_1_full % q
        if delta_1_0 != 0:
            print(f"  *** v_q(δ_1) < 1 for (κ_1={kappa_1}, a_1={a_1_full}) - ERROR? ***")
            continue
        delta_1_div_q = delta_1_full // q  # this is δ_1 / q in Z (since v_q ≥ 1)
        delta_1_1 = delta_1_div_q % q  # first non-trivial q-digit
        delta_1_2 = (delta_1_div_q // q) % q  # next q-digit (but this needs qsq-precision in delta_1_div_q, which requires v_q ≥ 2 of δ_1 OR we need more precision)
        #
        # Hmm, delta_1_full is mod qsq, so delta_1_div_q is mod q (after dividing by q).
        # delta_1_2 = (delta_1_div_q // q) % q is NOT well-defined from mod-q² δ_1 alone.
        # Need δ_1 mod q³ to get δ_1^(2).
        #
        # But (2^(-A_1) - 2^(-B_1)) has v_q = 1 (k_1 = ±1 coprime to q), so we have:
        # δ_1 = δ_1^(1)·q + δ_1^(2)·q² + δ_1^(3)·q³ + ... (no δ_1^(0) since v_q ≥ 1)
        #
        # δ_1^(2) requires δ_1 mod q³. We need pow_inv2 mod q³ = mod 4913.

# OK we need to lift further. Let me redo with q³-precision.
print("\n[Reset 2] Need q³-precision for δ_1^(2). Lifting pow_inv2 to mod q³.")
print("This requires ord_{q³}(2) = 8·17² = 2312 entries.")
print("Memory/time: ~2312 entries, feasible.")
