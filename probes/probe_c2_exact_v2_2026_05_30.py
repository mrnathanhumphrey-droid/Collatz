"""
probe_c2_exact_v2_2026_05_30.py

Clean c(2) derivation. Two cases:

  Case A: κ_1 = 0, level-2 has own sub-dom κ_2. Structurally identical to c(1).
          Relative contribution to c(2) num/den: same as c(1) num/den (since the
          mass is "absorbed" — Case A configs are 1/3 of the total mass, but Case B
          configs ALSO carry 1/3 mass scale, so the relative factors cancel in the
          ratio. Actually it's cleaner to define both with absolute weights.)

  Case B: κ_1 ≠ 0, level-2 (A_2, B_2) constrained to compensate δ_1^(1) =
          2κ_1·2^(-a_1). Need δ_2^(1) at q²-precision and δ_1^(2) at q³-precision.

We work in absolute weights and form c(2) = num/den.
"""
from fractions import Fraction
import sys, time
sys.stdout.reconfigure(encoding="utf-8")

q = 17
ord_2 = 8         # ord_q(2)
ord_2_qsq = 8 * q       # = 136
ord_2_q3  = 8 * q * q   # = 2312
qsq = q * q       # 289
q3  = q ** 3      # 4913

# === Precomputed power tables ===
# 2^(-a) mod q for a in 0..ord_2-1
inv2_q = pow(2, -1, q)
pow_inv2_q = [pow(inv2_q, a, q) for a in range(ord_2)]

# 2^(-a) mod q² for a in 0..ord_2_qsq-1
inv2_qsq = pow(2, -1, qsq)
pow_inv2_qsq = [pow(inv2_qsq, a, qsq) for a in range(ord_2_qsq)]

# 2^(-a) mod q³ for a in 0..ord_2_q3-1 (needed for δ_1^(2))
inv2_q3 = pow(2, -1, q3)
pow_inv2_q3 = []
val = 1
for a in range(ord_2_q3):
    pow_inv2_q3.append(val)
    val = (val * inv2_q3) % q3

print(f"Power tables built (ord_q={ord_2}, ord_q²={ord_2_qsq}, ord_q³={ord_2_q3})")

# === Helpers ===
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

def W_2312(r):
    if r == 0:
        return Fraction(1, 2**ord_2_q3 - 1)
    return Fraction(2**(ord_2_q3 - r), 2**ord_2_q3 - 1)

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
print(f"c(0) = {c0} = {float(c0):.15f}  (expected 19/127 ≈ 0.14961)")

# === c(1) absolute ===
# In c(1) script weights: w_0 = 1/3, w_κ = 1/(3·2^(8|κ|)) for κ ≠ 0.
# This corresponds to absolute P(κ_1 = κ | unconstrained Geom × Geom).
K_inner = 12
num_c1 = Fraction(0); den_c1 = Fraction(0)
for k in range(-K_inner, K_inner+1):
    wk = Fraction(1, 3) if k == 0 else Fraction(1, 3 * 2**(ord_2 * abs(k)))
    sk = (2 * k) % q
    num_c1 += wk * N[sk]
    den_c1 += wk * T[sk]
c1_val = num_c1 / den_c1
print(f"c(1) = {float(c1_val):.15f}  (expected 0.15317823)")

# === Case A: κ_1 = 0, level-2 sub-dom (κ_2 ∈ Z) ===
# σ_2 with κ_1=0 = depth-2 sub-dom κ_2 shift + (A_3, B_3) depth-0 character.
# Structurally identical to c(1) but at one level deeper.
# Conditional probability of κ_1 = 0 given level-1 alignment: this is the
# same (1/3) marginal that's already in num_c1, den_c1. So Case A is:
case_A_num = Fraction(1, 3) * num_c1
case_A_den = Fraction(1, 3) * den_c1

# Wait — let me think about this more carefully. In c(1):
#   c(1) = E[χ_2(σ_1) | v_q(D) ≥ 1]
# The conditioning v_q ≥ 1 picks (A_1, B_1) such that 8 | A_1 - B_1.
# Within this, κ_1 = 0 (A_1 = B_1) has prob (1/3)/(sum κ) = (1/3) / (256/255 · 1/3) ≈ 0.996.
#
# For c(2) = E[χ_2(σ_2) | v_q(D) = 2], we condition further. Case A is
# (κ_1 = 0) ∧ (v_q(D) = 2 via level-2 deepening). The level-2 structure
# is same as c(1)'s level-1.
#
# So Case A contribution to c(2) num = P(κ_1 = 0 | v_q = 2) · c(1)_num_relative
#                                    = (...) · num_c1.
#
# But Case A's "(A_3, B_3)" plays the role of c(1)'s "(A_2, B_2)" in N(s)/T(s).
# So the structural ratio Case A_num/Case A_den = c(1)_num/c(1)_den = c(1).
#
# Then Case A contribution = (weight) · c(1).

# Hmm — to get the absolute weights right, let me just compute Case A num/den
# absolutely from level-1 + level-2 marginals.
#
# Case A absolute: P(κ_1 = 0, level-2 has κ_2) · contribution
#   = (1/3) · num_c1  for num
#   = (1/3) · den_c1  for den
# where num_c1, den_c1 are themselves built from level-2 sub-dom weights × depth-0 character moments.

# Wait — c1_val = num_c1 / den_c1 is c(1). If Case A contributes (1/3)·num_c1 and (1/3)·den_c1,
# then (Case A only) gives c(2)_Case_A_only = num_c1/den_c1 = c(1).

# That makes sense because Case A IS c(1) by direct analogy.

# So Case A contribution to c(2) ratio = c(1) weighted by P(Case A given v_q = 2).
# And Case B contribution = (different shift sum) weighted by P(Case B given v_q = 2).
# Total: c(2) = [P(A)·c(1) + P(B)·(case B avg)] (a weighted average).

# In ratio form:
#   num_c2 = case_A_num + case_B_num
#   den_c2 = case_A_den + case_B_den
#   c(2) = num_c2 / den_c2

# === Case B: κ_1 ≠ 0, (A_2, B_2) compensates ===
K_outer = 3  # truncate |κ_1| ≤ K_outer; 2^(-8K_outer) ~ 10^(-7) for K=3, so ~25-digit truncation error

print(f"\nComputing Case B with K_outer={K_outer}...")
case_B_num = Fraction(0)
case_B_den = Fraction(0)

t0 = time.time()
total_configs = 0
for kappa_1 in range(-K_outer, K_outer + 1):
    if kappa_1 == 0:
        continue
    t_k = time.time()
    # Geom-weight of κ_1 = (1/3) · 2^(-8|κ_1|)
    w_kappa = Fraction(1, 3) * Fraction(1, 2**(ord_2 * abs(kappa_1)))

    # Loop a_1 mod 136 (full q² precision for δ_1^(1) since we need q² for δ_1^(2))
    # Actually δ_1^(2) needs q³ precision in δ_1, which needs a_1 mod 2312.
    # Hmm let me reconsider.

    # δ_1 = 2^(-A_1) - 2^(-B_1). A_1 = a_1 + 8 max(κ_1,0), B_1 = a_1 + 8 max(-κ_1,0).
    # For κ_1 = +1, A_1 = a_1+8, B_1 = a_1. δ_1 = 2^(-a_1) (2^(-8) - 1).
    # 2^(-a_1) mod q³ depends on a_1 mod 2312.
    # (2^(-8) - 1) mod q³ is a specific number (depending on 2^(-8) mod q³).

    # δ_1 mod q³ = 2^(-a_1) (2^(-8) - 1) mod q³.
    # δ_1^(1) = (δ_1 // q) mod q = ((2^(-a_1) (2^(-8) - 1)) // q) mod q.
    #   This only depends on 2^(-a_1) mod q²·17 = q² level... wait depends on a_1 mod 136.
    # Let me just compute:
    #   δ_1 mod q² depends on a_1 mod 136 (and κ_1).
    #   δ_1 / q mod q (the first non-zero q-digit) depends on a_1 mod 136.
    #   δ_1 mod q³ depends on a_1 mod 2312.
    #   δ_1^(2) depends on a_1 mod 2312.

    # So a_1 must be lifted to mod 2312 for full q³ precision of δ_1^(2).
    # That's 2312 values per κ_1. Then for each, A_2 mod 136 needed... 2312·136² = 43M per κ_1.
    # Too slow at 6 κ_1 values = 256M configs at Fraction = hours.

    # Truncate: a_1 lift to mod 136 only (ignoring q³ correction to δ_1^(2)).
    # Then δ_1^(2) approximation: derived from δ_1 mod q² extrapolated.

    # Actually we want δ_1^(2) mod q. If we only have δ_1 mod q², we miss the q² coefficient.
    # Better approach: derive δ_1^(2) as a closed form via (2^(-8κ) - 1)/q² mod q.

    # (2^(-8κ_1) - 1)/q² mod q can be derived analytically.
    # 2^(-8κ_1) - 1 has v_q = 1 (LTE). So (2^(-8κ_1) - 1)/q has v_q = 0.
    # (2^(-8κ_1) - 1)/q mod q² can be computed from 2^(-8κ_1) mod q³.
    # Specifically, we need 2^(-8κ_1) mod q³ for (2^(-8κ_1) - 1)/q² mod q.

    # Let me compute this directly.
    # 2^(-8) mod q³: we have pow_inv2_q3[8] = 2^(-8) mod q³.
    # 2^(-8κ_1) mod q³ = (2^(-8))^κ_1 mod q³ for positive κ_1, or 2^(8|κ_1|) for negative.
    if kappa_1 > 0:
        two_neg_8k_q3 = pow(pow_inv2_q3[8 % ord_2_q3], kappa_1, q3)
    else:
        two_neg_8k_q3 = pow(pow(2, 8, q3), abs(kappa_1), q3)
    # (2^(-8κ_1) - 1) mod q³ (interpret -1 mod q³, so this can be in {0..q³-1})
    expr_q3 = (two_neg_8k_q3 - 1) % q3
    # Verify v_q ≥ 1: expr_q3 % q should = 0
    if expr_q3 % q != 0:
        print(f"  ERROR: 2^(-8κ_1) - 1 not ≡ 0 mod q for κ_1 = {kappa_1}")
        continue
    # (2^(-8κ_1) - 1) / q mod q²: well-defined since v_q ≥ 1
    expr_div_q = expr_q3 // q  # mod q²
    delta_1_1_per_a1 = expr_div_q % q  # contributes 2^(-a_1)·this to δ_1^(1)
    expr_div_qsq = expr_div_q // q  # mod q (this is the q^2 coefficient)
    delta_1_2_per_a1 = expr_div_qsq % q  # contributes 2^(-a_1)·this to δ_1^(2)

    # For each a_1 mod 136 (the relevant precision for δ_2^(1) at q² is mod 136):
    for a_1 in range(1, ord_2_qsq + 1):  # a_1 = 1, 2, ..., 136 (representing residue class mod 136)
        # Actually a_1_res mod 136 = a_1 % 136 (which for a_1 in {1..136} just = a_1, with 136 ≡ 0).
        a_1_mod_qsq = a_1 % ord_2_qsq
        a_1_mod_q = a_1 % ord_2  # mod 8
        # Geom weight on a_1 residue mod 136: W_136(a_1_mod_qsq).
        w_a1 = W_136(a_1_mod_qsq)

        pow_inv2_a1_qsq = pow_inv2_qsq[a_1_mod_qsq]
        pow_inv2_a1_q = pow_inv2_a1_qsq % q

        # δ_1^(1) = pow_inv2_a1_q · delta_1_1_per_a1 mod q
        delta_1_1 = (pow_inv2_a1_q * delta_1_1_per_a1) % q
        # δ_1^(2) = pow_inv2_a1_q · delta_1_2_per_a1 mod q (approximation - exact only if δ_1 doesn't have cross terms)
        # Actually δ_1 = 2^(-a_1) · (2^(-8κ_1) - 1). So δ_1 / q = 2^(-a_1) · expr_div_q.
        # δ_1 mod q^3 = 2^(-a_1) · (2^(-8κ_1) - 1) mod q^3 = pow_inv2_a1_q3 · expr_q3 mod q^3.
        # need pow_inv2_a1 mod q^3 not just q^2. So a_1 mod 2312, not 136. Approximation:
        # For a_1 mod 136 we have 2^(-a_1) mod q^2. The q^2 coefficient of 2^(-a_1) requires q^3 lift.
        #
        # Let me compute δ_1^(2) approximately as (delta_1_per_a1 · 2^(-a_1)_q^2 piece) mod q.
        # Actually: δ_1 = 2^(-a_1) · expr. δ_1 / q = 2^(-a_1) · (expr/q).
        # δ_1^(2) = (δ_1/q^2) mod q = (2^(-a_1) · (expr/q^2)) mod q ... if expr/q^2 has integer part.
        #
        # Simpler: compute δ_1 mod q^3 directly.
        # Need 2^(-a_1) mod q^3, which requires a_1 mod 2312.
        # For a_1 ∈ {1..136} = first period of mod 136, the lift to mod 2312 is just a_1 itself.
        # For a_1 > 136 it cycles. We're iterating a_1 in {1..136}, so use a_1 directly.
        pow_inv2_a1_q3 = pow_inv2_q3[a_1 % ord_2_q3]  # a_1 < ord_2_q3 = 2312 here
        delta_1_q3 = (pow_inv2_a1_q3 * expr_q3) % q3
        # δ_1^(0) should be 0
        assert delta_1_q3 % q == 0, f"v_q(δ_1) < 1 for a_1={a_1}, κ_1={kappa_1}"
        delta_1_div_q_actual = delta_1_q3 // q  # mod q²
        delta_1_1_actual = delta_1_div_q_actual % q
        delta_1_2_actual = (delta_1_div_q_actual // q) % q
        # These should match the "per_a1" values multiplied by 2^(-a_1) mod q.
        # delta_1_1_actual should = (pow_inv2_a1_q · delta_1_1_per_a1) % q
        # delta_1_2_actual is the true value (using full q³ precision).
        delta_1_1 = delta_1_1_actual
        delta_1_2 = delta_1_2_actual

        # Constraint on (A_2, B_2): δ_2^(0) ≡ -δ_1^(1) mod q
        s_constraint = (-delta_1_1) % q

        # Now sum over (A_2 mod 136, B_2 mod 136) with 2^(-A_2) - 2^(-B_2) ≡ s_constraint mod q.
        # The constraint is on mod 8, but we need mod 136 for δ_2^(1).
        # For each (r_a, r_b) ∈ {0..7}² with pow_inv2_q[r_a] - pow_inv2_q[r_b] ≡ s_constraint mod q:
        #   For each lift j_a ∈ {0..16} (A_2 = r_a + 8·j_a, j_a ≥ 1 if r_a = 0 else j_a ≥ 0; also A_2 ≤ 136 in first period):
        #     For each lift j_b ∈ {0..16}:
        #       Compute δ_2^(1) and accumulate.
        # But careful: we want a representative period, so A_2 ∈ {1, ..., 136}. With r_a = 0, A_2 ∈ {8, 16, ..., 136}.
        #              With r_a ≠ 0, A_2 ∈ {r_a, r_a+8, ..., r_a+128}. So 17 lifts each.

        for r_a in range(ord_2):
            for r_b in range(ord_2):
                if (pow_inv2_q[r_a] - pow_inv2_q[r_b]) % q != s_constraint:
                    continue
                # For each lift to mod 136
                # A_2 ∈ {r_a, r_a+8, r_a+16, ..., r_a+128} if r_a > 0, else {8, 16, ..., 136}
                lifts_A = list(range(r_a if r_a > 0 else ord_2, ord_2_qsq + 1, ord_2))
                lifts_B = list(range(r_b if r_b > 0 else ord_2, ord_2_qsq + 1, ord_2))
                for A_2 in lifts_A:
                    for B_2 in lifts_B:
                        pow_A2_qsq = pow_inv2_qsq[A_2 % ord_2_qsq]
                        pow_B2_qsq = pow_inv2_qsq[B_2 % ord_2_qsq]
                        delta_2_qsq = (pow_A2_qsq - pow_B2_qsq) % qsq
                        # δ_2^(0)
                        delta_2_0 = delta_2_qsq % q
                        if delta_2_0 != s_constraint:
                            # Lift didn't preserve the constraint at q² level (shouldn't happen but safety)
                            continue
                        # δ_2^(1) = ((δ_2 - δ_2^(0))/q) mod q
                        delta_2_div_q = (delta_2_qsq - delta_2_0) // q
                        delta_2_1 = delta_2_div_q % q

                        # σ_2 shift = (δ_1^(2) + δ_2^(1)) mod q
                        shift = (delta_1_2 + delta_2_1) % q

                        # Weight: w_kappa × w_a1 × W_136(A_2) × W_136(B_2)
                        weight = w_kappa * w_a1 * W_136(A_2 % ord_2_qsq) * W_136(B_2 % ord_2_qsq)

                        # Accumulate with N(shift), T(shift)
                        case_B_num += weight * N[shift]
                        case_B_den += weight * T[shift]
                        total_configs += 1
    print(f"    κ_1={kappa_1} done in {time.time()-t_k:.1f}s, total configs={total_configs}")

print(f"\nTotal Case B configs: {total_configs}")
print(f"Case B compute time: {time.time()-t0:.1f}s")

# === Combine ===
num_c2 = case_A_num + case_B_num
den_c2 = case_A_den + case_B_den
c2 = num_c2 / den_c2

print(f"\n=== Result ===")
print(f"Case A num = {float(case_A_num):.10e}")
print(f"Case A den = {float(case_A_den):.10e}")
print(f"Case B num = {float(case_B_num):.10e}")
print(f"Case B den = {float(case_B_den):.10e}")
print(f"\nc(2) computed = {float(c2):.15f}")
print(f"c(2) measured = 0.15324792077909  (FFT, mpmath dps=40)")
print(f"diff = {float(c2) - 0.15324792077908744:.2e}")

# Save the exact rational
print(f"\nExact c(2) numerator   = {c2.numerator}")
print(f"Exact c(2) denominator = {c2.denominator}")
print(f"Num digits: {len(str(c2.numerator))}, Den digits: {len(str(c2.denominator))}")
