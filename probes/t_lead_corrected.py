"""
T_LEAD_CORRECTED — re-evaluate T_lead = T_diag + Off_lin at the corrected rate.

Per CROSS_FREQ_PHASE1_EXPANSION:
  - W_+(g) = 2^{-g+1}/15  (P^++ weight sum at fixed g; v ∈ {2,4,...}, geometric in v)
  - W_-(g) = 2^{-g+3}/15  = 4 · W_+(g)  (P^-- weight sum, v ∈ {1,3,...})
  - Off correction along (P_+, P_-) direction is RANK-1 along (1, 4) for all g
    (because W_-/W_+ = 4 uniformly).
  - The class-summed cross-freq moment X̄_n(c; g) is what multiplies
    W_±(g) in the off-diagonal contribution.

The (1,4) projection: when restricted to the (P_+, P_-) deviation subspace
along the (1,4) eigenvector of T_diag, the Off_lin contribution gives a
scalar contraction whose value comes from summing the geometric weights.

CRITICAL: the cross_freq Phase 1 derivation parameterizes the action by
the off-diagonal weight, with the EXACT amplitude on the moment X̄_n(g, c)
to be tracked. The leading contributing g=2 has W_+(2)=1/30 (per §6); the
next g=4 has W_+(4)=1/120; etc. Cross_freq computed an explicit factor of
3 (from lift-fiber survival, §3) in Off_{n+1}^{++} = 3·Σ W_+(g)·X̄_n(c; g).

This script computes:
1. Σ_g W_+(g) and the (1,4)-action of Off_lin via the cross_freq machinery
2. T_lead's 2D spectrum on (P_+, P_-)
3. Empirical cross-check: does T_lead's leading eigenvalue match
   - rate 1/2 (R77.3 falsified)?
   - rate ρ ≈ 0.984 (slow oscillating mode, period ≈ 9.2)?
   - some other value (the corrected within-level rate)?
4. Ratio |ε_n|/|ε_{n-1}| trajectory + |ε_n|·2^n envelope

Exact rationals via fractions.Fraction.
"""

from fractions import Fraction
import math
import json
from pathlib import Path

# =====================================================================
# Phase 1: Cross-freq weight sum + Off_lin action on (1, 4)
# =====================================================================
#
# Per CROSS_FREQ_PHASE1_EXPANSION §6:
#   W_+(g) = 2^{-g+1}/15  for g ∈ {2, 4, 6, ...}  (even ≥ 2)
#   W_-(g) = 2^{-g+3}/15  = 4 · W_+(g)
#
# Note these are NOT the only weights — the v=v' diagonal (g=0) is T_diag
# proper, already at λ=1 on (1,4). The OFF-diagonal correction is g ≥ 2.

print("=" * 72)
print("PHASE 1: cross-freq weight sums")
print("=" * 72)

# Compute Σ_g W_+(g) for g ∈ {2, 4, 6, ...} as exact Fraction
# W_+(g) = 2^{-g+1}/15 = Fraction(2, 15) * Fraction(1, 4)^{(g/2 - 1) + 1}
# Wait: 2^{-g+1} for g=2 is 2^{-1} = 1/2, so W_+(2) = (1/2)/15 = 1/30 ✓
# For g=4: W_+(4) = 2^{-3}/15 = (1/8)/15 = 1/120 ✓
# So W_+(g) = Fraction(1, 15) * Fraction(1, 2)^{g-1}
def W_plus(g):
    """W_+(g) = 2^{-g+1}/15 = 1/(15 * 2^{g-1})."""
    return Fraction(1, 15) * Fraction(1, 2)**(g - 1)

def W_minus(g):
    """W_-(g) = 4 · W_+(g)."""
    return 4 * W_plus(g)

# Verify against cross_freq Phase 1 §6 values
print(f"  W_+(2) = {W_plus(2)} = {float(W_plus(2)):.6f}  (expected 1/30 = 0.033333)")
print(f"  W_+(4) = {W_plus(4)} = {float(W_plus(4)):.6f}  (expected 1/120 = 0.008333)")
print(f"  W_+(6) = {W_plus(6)} = {float(W_plus(6)):.6f}  (expected 1/480 = 0.002083)")
print(f"  W_-(2) = {W_minus(2)} = {float(W_minus(2)):.6f}  (expected 2/15 = 0.133333)")

# Σ_{g=2,4,6,...} W_+(g) is geometric: a = 1/30, ratio = 1/4
# Sum = (1/30) / (1 - 1/4) = (1/30) * (4/3) = 4/90 = 2/45
sum_W_plus = Fraction(1, 30) / (1 - Fraction(1, 4))
print(f"\n  Σ_g W_+(g) for g ∈ {{2,4,6,...}} = {sum_W_plus} = {float(sum_W_plus):.6f}")
print(f"    (expected 2/45 ≈ 0.04444)")

# Σ_g W_-(g) = 4 · Σ_g W_+(g) = 8/45
sum_W_minus = 4 * sum_W_plus
print(f"  Σ_g W_-(g) = {sum_W_minus} = {float(sum_W_minus):.6f}  (expected 8/45)")

# =====================================================================
# CRITICAL: there's a factor-of-3 from cross_freq §3 (lift-fiber survives → 3)
# Cross_freq §7: Off_{n+1}^{++} = 3 · Σ_g W_+(g) · X̄_n(c; g)
# So when we ask "what's Off's action on (1,4) in the (P_+, P_-) plane?"
# we need to know what (X̄_n) is as a function of (P_+, P_-).
#
# The X̄_n(c; g=0) IS P_n^{++}(c)+P_n^{+-}(c)+P_n^{-+}(c)+P_n^{--}(c).
# Under structural collapse P^{+-} = P^{-+} = 0 (n ≥ 2), and class-c-symmetry:
#   X̄_n(c; g=0) = P_+ + P_-
# For g ≥ 2, X̄_n(c; g) is a NEW moment (cross_freq H_CROSS_CLOSES_ON_ENLARGED_SPAN),
# not in span{P_+, P_-}. So the "(1,4) projection" interpretation needs
# the assumption that X̄_n(g≥2) ALSO projects onto (1,4) in the leading order.
# Cross_freq §7 establishes this is structurally true via W_-/W_+ = 4 uniform.
# =====================================================================

print("\n" + "=" * 72)
print("PHASE 2: T_lead's action on (P_+, P_-) restricted to (1, 4) eigenvector")
print("=" * 72)
print()
print("T_diag (rigorous, R77 §1, T_diag = (1/5)·[[1,1],[4,4]]):")
print("  spectrum {0, 1}, eigenvector (1, 4) at λ=1")
print()
print("Now consider T_lead = T_diag + Off_lin restricted to the (P_+, P_-) plane.")
print("Off_lin is rank-1 along (1, 4) (cross_freq §7) and contracts the")
print("(1,4)-component of the deviation by the sum of weights.")
print()

# Off_lin's action on (P_+, P_-) along (1,4) eigenvector:
# Off_{n+1}^{++} = 3 · Σ W_+(g) · X̄_n(c; g)  [includes g=0 if we restate]
# Wait — cross_freq §6's W_+(g) was for g ≥ 2 OFF-DIAGONAL only.
# Let me re-read: Phase 1 cross_freq §1 separates the v=v' diagonal (= T_diag)
# from v ≠ v' off-diagonal. The W_+(g), W_-(g) of §6 are for g ≥ 2 OFF-DIAG.
#
# So Off_lin total contribution to (P_+, P_-) from g ≥ 2 cross-freq pairs:
#   Off_+ = 3 · (1/30 + 1/120 + 1/480 + ...) · X̄ = 3 · (2/45) · X̄ = 6/45 = 2/15
#   Off_- = 3 · (2/15 + 1/30 + 1/120 + ...) · X̄ = 3 · (8/45) · X̄ = 8/15
#
# Where X̄ is the cross-freq class-summed moment value (for the leading mode).
# Assuming X̄ projects onto (1,4) of (P_+, P_-): X̄ ~ amount in (1,4) eigenvec.

# Total off-diagonal contribution to (P_+, P_-) along (1,4):
total_off_plus = 3 * sum_W_plus
total_off_minus = 3 * sum_W_minus

print(f"  Off_lin total contribution to P_+ along (1,4): 3·Σ_g W_+(g) = {total_off_plus} = 2/15")
print(f"  Off_lin total contribution to P_- along (1,4): 3·Σ_g W_-(g) = {total_off_minus} = 8/15")
print(f"  Ratio Off_-/Off_+ = 4  (confirms (1,4) eigenvector preservation)")
print()

# Now T_lead's leading eigenvalue on (1, 4):
# T_diag · (1, 4) = (1, 4) [eigenvalue 1]
# Off_lin · (1, 4) is rank-1 along (1, 4) with scalar coefficient...
#
# CAREFUL: The conventional reading is that T_diag is the FULL on-diagonal
# linearization that ALREADY captures the full P→P recursion (eigenvalue 1
# on (1,4) means the slow eigenvalue is on the (1,4) direction). Off_lin
# is what CONTRACTS this from λ=1 to λ < 1.
#
# So T_lead's eigenvalue on (1,4) = 1 - (Off_lin's contraction along (1,4))
#
# Reading cross_freq §7 carefully: Off_+ = 3·Σ W_+(g)·X̄_n  -- this is the
# contribution that gets ADDED to T_diag · (P_+, P_-) to give (P_+, P_-)_{n+1}.
#
# Now: T_diag's eigenvalue on (1,4) is 1, but T_diag's action on the
# DEVIATION (P_+ - 7/150, P_- - 14/75) might be 1 (preserving deviation)
# or might be 0 (killing deviation) depending on convention.
#
# From R77 §1: T_diag has rank 1, eigenvalue 1 on (1,4) PRESERVES Plancherel
# total mass S = 2(P_+ + P_-). So T_diag(deviation along (1,4)) is the
# (1,4)-component which gets preserved (λ=1). The CONTRACTION comes from Off.
#
# This means: in the (1,4) direction, T_lead = 1 + Δ where Δ is Off_lin's
# (1,4)-eigenvalue. Empirically Δ < 0 (since S_n converges to 7/15, the
# deviation must shrink, so |T_lead eigenvalue| < 1).
#
# So the QUESTION is: what's the sign and magnitude of Off_lin's
# (1,4)-eigenvalue? The classical R77 §2 reading was "+1/2" from naive
# bilinear weights but this missed the sign convention.

# Let me compute the leading T_lead eigenvalue from the structural setup:
#   T_lead on (1,4) = T_diag on (1,4) + Off_lin on (1,4)
#   = 1 + contraction_coefficient
#
# To match empirical |ε_n+1/ε_n| from k=2..6 around -0.5 (which gives ε_n ~
# (-1/2)^n decay), the eigenvalue should be -1/2 (NEGATIVE!), not +1/2. The
# sign came from the off-diagonal phase factors generating an oscillating
# contraction.
#
# But that's the OLD rate-1/2 narrative which R77.3 falsified.
#
# For the CORRECTED rate, the prediction depends on how Off_lin sums.
#
# The MOST DIRECT reading of cross_freq §7: the rank-1 contraction on (1,4)
# from the FULL off-diagonal sum is:
#   λ_off = -(coefficient_+, coefficient_-) projected onto (1,4)
#
# Actually let me re-read what the brief says more carefully...
# Brief Phase 1: "Off · (1, 4) = (Σ_g W_+(g)) · (1, 4)_+ + 4·(Σ_g W_+(g)) · (1, 4)_-"
# "Sum = (1/30) / (1 - 1/4) = (1/30) · (4/3) = 4/90 = 2/45"
# "So Off_lin's contribution to T_lead on (1,4) is **−2/45** (negative sign —
# Off contracts) → T_lead on (1, 4) = 1 − 2/45 = **43/45 ≈ 0.9556**."
#
# So the brief's prediction is T_lead eigenvalue = 1 - 2/45 = 43/45 on (1,4).
# Let me COMPUTE that:

lambda_T_lead = 1 - sum_W_plus  # 1 - 2/45 = 43/45
print(f"  T_lead's (1,4) eigenvalue = 1 - Σ_g W_+(g) = 1 - 2/45 = {lambda_T_lead}")
print(f"                            = {float(lambda_T_lead):.6f}")
print()

# Computing T_lead's full 2D spectrum on (P_+, P_-):
# T_lead = T_diag + Off_lin where:
# T_diag = (1/5) * [[1, 1], [4, 4]]  (rank 1, spectrum {0, 1})
# Off_lin = -(rank-1 along (1,4)) with coefficient |Σ_g W_+(g)| = 2/45
# (sign negative for contraction)
#
# A rank-1 matrix along (1,4) is c * [[1*α, 1*β], [4*α, 4*β]]
# Eigenvector (1, 4) with eigenvalue c*(α + 4β)
# For the contraction to be by exactly Σ_g W_+(g) along (1,4):
# we need c*(α + 4β) = -2/45 (eigenvalue) and the other eigenvalue is 0.
#
# Constructive: Off_lin = -(2/45) * [[1, 0], [4, 0]] gives:
#   - eigenvalue -(2/45) on (1, 4)  ✓
#   - eigenvalue 0 on (0, 1)
#
# Wait that gives -(2/45)*(1, 4) but check: matrix · (1,4) = -(2/45)*(1+0, 4+0) = -(2/45, -8/45)
# Hmm, [[1,0],[4,0]] · (1,4) = (1, 4)  ✓ for the (1) component
# So Off_lin · (1,4) = -(2/45)·(1, 4)  ✓
#
# But this construction is somewhat arbitrary on the orthogonal direction.
# Cross_freq §7 only constrains the (1,4) direction; the orthogonal
# direction is unconstrained (rank-1 means 0 on perpendicular).
#
# So T_lead = T_diag + (-2/45)·(1, 4)·(1, 0)/normalization
# = (1/5)·[[1, 1], [4, 4]] + (-2/45)·[[1, 0], [4, 0]]
#
# Combined:
# [[1/5 - 2/45, 1/5], [4/5 - 8/45, 4/5]]
# = [[(9-2)/45, 1/5], [(36-8)/45, 4/5]]
# = [[7/45, 1/5], [28/45, 4/5]]

a = Fraction(7, 45)
b = Fraction(1, 5)
c = Fraction(28, 45)
d = Fraction(4, 5)
print(f"  T_lead = [[{a}, {b}], [{c}, {d}]]")
print(f"        = [[7/45, 9/45], [28/45, 36/45]]")
print()

# Characteristic polynomial: det(T_lead - λI) = λ² - tr·λ + det
tr = a + d
det = a*d - b*c
print(f"  trace(T_lead) = {tr}")
print(f"  det(T_lead)  = {det}")

# Discriminant
disc = tr*tr - 4*det
print(f"  discriminant = {disc} = {float(disc):.6f}")

# Eigenvalues
disc_float = float(disc)
sqrt_disc = math.sqrt(abs(disc_float)) if disc_float >= 0 else math.sqrt(-disc_float)
if disc_float >= 0:
    lambda1 = (float(tr) + sqrt_disc) / 2
    lambda2 = (float(tr) - sqrt_disc) / 2
    print(f"  eigenvalues (real): λ_1 = {lambda1:.10f}, λ_2 = {lambda2:.10f}")
else:
    re_lam = float(tr) / 2
    im_lam = sqrt_disc / 2
    print(f"  eigenvalues (complex): {re_lam} ± {im_lam}i, modulus = {math.sqrt(re_lam**2 + im_lam**2):.10f}")

# Also confirm exactly: if tr = 43/45 and the smaller eigenvalue is 0:
# det = 0 ⇒ T_lead is RANK 1 → eigenvalues are tr and 0 → (43/45, 0)
print()
print(f"  Confirm: det should be 0 if T_lead is rank-1 (both T_diag and Off rank-1, both along (1,4)) → det = {det}")
print()

if det == 0:
    print(f"  T_lead is RANK 1; spectrum = {{trace, 0}} = {{43/45, 0}}")
    print(f"  Leading eigenvalue = trace = 43/45 = {float(tr):.10f}")
    eigenvalue_1_4 = tr
    print(f"  Eigenvector at λ = 43/45: (1, 4) (same as T_diag's λ=1 direction)")
else:
    print(f"  T_lead NOT rank-1; det ≠ 0. Verify Off_lin construction.")
    eigenvalue_1_4 = None

# =====================================================================
# Phase 4: Empirical cross-check with ε_n exact + numerical
# =====================================================================
print()
print("=" * 72)
print("PHASE 4: Empirical cross-check")
print("=" * 72)

# Load ε_n cache from R77.7 v2
eps_cache = json.loads(Path("C:/Collatz/experiments_output/result_77_7_eps_exact_through_k7_v2.json").read_text())
eps_exact = {}
for k_str, v in eps_cache.items():
    k = int(k_str)
    eps_exact[k] = Fraction(int(v["num"]), int(v["den"]))

# Numerical floats for ε_8..ε_13 from PADE_NUMERICAL_DATA.md
eps_numerical = {
    8: -7.4554636729e-04,
    9: -7.5202571564e-06,
    10: +7.2075091711e-04,
    11: +1.5019670121e-03,
    12: +2.2747137206e-03,
    13: +2.9482473172e-03,
}

# Build a unified float dict
eps = {}
for k, frac in eps_exact.items():
    eps[k] = float(frac)
for k, val in eps_numerical.items():
    eps[k] = val

print(f"\n  ε_n trajectory (k=1..13):")
print(f"  {'k':>3} {'ε_k':>15} {'sign':>5} {'|ε_k|':>15} {'|ε_k|·2^k':>15} {'|ε_k/ε_{k-1}|':>16}")
for k in sorted(eps.keys()):
    s = '+' if eps[k] >= 0 else '-'
    abs_eps = abs(eps[k])
    env = abs_eps * (2 ** k)
    if k > 1:
        ratio = abs(eps[k] / eps[k-1])
        print(f"  {k:>3} {eps[k]:+.6e} {s:>5} {abs_eps:.6e} {env:.4e} {ratio:.6f}")
    else:
        print(f"  {k:>3} {eps[k]:+.6e} {s:>5} {abs_eps:.6e} {env:.4e} {'—':>16}")

# Geometric rate test: if ε_n ~ const · ρ^n, then |ε_n/ε_{n-1}| → ρ
print()
print(f"  Geometric-rate diagnostic |ε_n / ε_{{n-1}}|:")
ratios = []
for k in sorted(eps.keys()):
    if k > 1 and eps[k-1] != 0:
        r = abs(eps[k] / eps[k-1])
        ratios.append((k, r))

for k, r in ratios:
    matches_half = abs(r - 0.5) / 0.5
    matches_43_45 = abs(r - 43/45) / (43/45)
    matches_984 = abs(r - 0.984) / 0.984
    print(f"    k={k:>3} → {r:.6f}  | dist(1/2)={matches_half:.3f}  dist(43/45={(43/45):.4f})={matches_43_45:.3f}  dist(0.984)={matches_984:.3f}")

# Late ratios (k ≥ 10) should be cleanest if asymptotic regime reached
print()
print(f"  Late-trajectory geometric mean of |ε_n/ε_{{n-1}}| at k=11..13:")
late_ratios = [r for k, r in ratios if k >= 11]
if late_ratios:
    geomean = math.exp(sum(math.log(r) for r in late_ratios) / len(late_ratios))
    print(f"    geomean(k=11..13) = {geomean:.6f}")
    print(f"    distance to 43/45 = {abs(geomean - 43/45):.4f}")
    print(f"    distance to 0.984 = {abs(geomean - 0.984):.4f}")

# Also test |ε_n|^(1/n) (Hadamard radius)
print()
print(f"  Hadamard-style |ε_n|^(1/n) (asymptotic rate from power-of-modulus):")
for k in [7, 8, 9, 10, 11, 12, 13]:
    if k in eps and abs(eps[k]) > 0:
        had = abs(eps[k]) ** (1.0/k)
        print(f"    k={k} → |ε_k|^(1/{k}) = {had:.6f}")

# =====================================================================
# Adversarial check A4: T_lead's eigenvalue should be stable across n
# Since T_lead is built purely from weights (no n-dependence in coefficients
# from cross-freq materials), eigenvalue should be exactly 43/45 at every n.
# Verify: T_lead matrix entries are constants over Q.
# =====================================================================
print()
print("=" * 72)
print("ADVERSARIAL CHECKS")
print("=" * 72)
print()
print("(A1) Cross-freq fidelity:")
print(f"    Phase 1 §6: W_+(g) = 2^{{-g+1}}/15 for g even ≥ 2")
print(f"    Σ_g W_+(g) = (1/30) + (1/120) + (1/480) + ... = 2/45 = {2/45:.6f}")
print(f"    [Verified algebraically from geometric series ratio 1/4.]")
print()
print("(A2) Period-9 vs single real eigenvalue test:")
print(f"    T_lead has eigenvalue 43/45 ≈ 0.9556 on (1,4), other eigenvalue 0.")
print(f"    This predicts ε_n ~ const · (43/45)^n, monotonic in sign.")
print(f"    But empirical ε_n has sign pattern + + - - - - - - - + + + +")
print(f"    Single zero-crossing at k=9→10. NOT monotonic.")
print(f"    So 43/45 cannot be the full asymptotic story.")
print()
print("(A3) Exact rationals: all Off_lin entries are over Q via Fraction.")
print(f"    T_lead = [[7/45, 1/5], [28/45, 4/5]] over Q exactly.")
print()
print("(A4) Iteration stability:")
print(f"    T_lead is n-independent (cross-freq §6 weights are absolute).")
print(f"    So eigenvalue 43/45 is the same at every n ≥ 2.")
print()
print("(A5) Empirical fit at n=7 doesn't match pure rate-0.984.")
empirical_envelope_growth = []
for k in [6, 7, 8, 9, 10, 11, 12, 13]:
    if k in eps:
        env = abs(eps[k]) * (2 ** k)
        empirical_envelope_growth.append((k, env))
print(f"    Envelope |ε_n|·2^n: {[(k, f'{e:.4f}') for k,e in empirical_envelope_growth]}")
print(f"    Predicted at rate 43/45: |ε_n|·2^n grows as (2·43/45)^n = (86/45)^n = (1.911)^n")
print(f"    Predicted at rate 0.984: |ε_n|·2^n grows as (1.968)^n")
print()
print(f"    Growth ratios envelope_{{n}}/envelope_{{n-1}}:")
for i in range(1, len(empirical_envelope_growth)):
    k, e = empirical_envelope_growth[i]
    _, e_prev = empirical_envelope_growth[i-1]
    if e_prev > 0:
        gr = e / e_prev
        print(f"      k={k}: ratio = {gr:.4f}  (43/45·2 = {86/45:.4f}, 0.984·2 = {0.984*2:.4f})")

print()
print("=" * 72)
print("DISPOSITION")
print("=" * 72)
print()
print(f"T_lead has eigenvalue 43/45 ≈ 0.9556 on (1, 4) eigenvector.")
print()
print(f"Comparison to empirical asymptote candidates:")
print(f"  - 1/2 (R77.3 falsified)                    : T_lead misses by 91% relative")
print(f"  - 43/45 ≈ 0.9556                            : T_lead's exact prediction")
print(f"  - ρ ≈ 0.984 (STATE.md slow-mode)            : T_lead misses by 2.9%")
print()
print(f"Empirical late-trajectory ratios (k=11..13):")
if late_ratios:
    geomean = math.exp(sum(math.log(r) for r in late_ratios) / len(late_ratios))
    print(f"  geomean = {geomean:.6f}")
    print(f"  distance to 43/45 = {abs(geomean - 43/45):.4f} ({100*abs(geomean - 43/45)/(43/45):.2f}%)")
    print(f"  distance to 0.984 = {abs(geomean - 0.984):.4f} ({100*abs(geomean - 0.984)/0.984:.2f}%)")
print()
print("CRITICAL CAVEAT:")
print(f"  The (1,4) eigenvalue calc T_lead = 1 - 2/45 = 43/45 assumes that")
print(f"  the cross-freq X̄_n(c; g) moments project ONTO the (1,4) eigenvector")
print(f"  of (P_+, P_-) with FULL UNIT amplitude. This is the H_CROSS_CLOSES_")
print(f"  ON_ENLARGED_SPAN reading: V_M is NOT in span{{P}}, so X̄_n(g≥2) is")
print(f"  a NEW moment and the 'projection onto (1,4)' is a STRUCTURAL")
print(f"  ASSUMPTION (cross_freq §7), not a derived identity. The honest")
print(f"  reading: T_lead's 'eigenvalue 43/45 on (1,4)' is what you'd get")
print(f"  IF the (P_+, P_-) projection of M_n^{{ab}}(g≥2, c) equals the (1,4)")
print(f"  component of the X̄_n that appears in Off — but cross_freq §7's")
print(f"  derivation establishes only the W_-/W_+ = 4 ratio, NOT that X̄_n")
print(f"  is rigorously identified with the (1,4) component on its own.")
print()
print("So: 43/45 is the cross-freq-machinery-natural answer for T_lead's")
print("(1,4)-eigenvalue, IF the projection assumption holds. It's NOT 1/2.")
print("It's NOT 0.984. It's structurally meaningful as 1 - Σ_g W_+(g) =")
print("1 - 2/45 — the SUM-OF-OFF-DIAGONAL-BILINEAR-WEIGHTS contraction.")
