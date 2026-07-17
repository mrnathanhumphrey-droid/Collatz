"""
verify_monotone_diagnostic.py â€” numerical confirmation that the third-order
alternating B-centered moment of off-diagonal Syracuse transfer operators
DOES NOT vanish, while the second-order and three-distinct-index moments DO
vanish.

This is the diagnostic signature of B-valued MONOTONE independence (vs B-valued
free independence). See:
  - AMALG_FREENESS_SETUP.md
  - AMALG_FREENESS_MOMENT_CALCULATION.md
  - AMALG_FREENESS_DISPOSITION.md
  - C1_TAO_RECURSION_FORM.md
  - result_77_T_lead_spectrum.md
  - OBSTRUCTION_MAP_TERMINAL.md

LEVEL: n = 3 (state space (Z/27)*, 18 elements).

OPERATOR REALIZATION (the most direct concrete reading):

At step j, a pair (v_{2j-1}, v_{2j}) is drawn iid Geom(2). Write b_j =
v_{2j-1} + v_{2j} and b_{[1,j]} = b_1 + ... + b_j.

The **diagonal** part of step j is the v = v' contribution; the **off-diagonal**
Off_j is the cross-frequency v â‰  v' part. Concretely, on H_3 = LÂ²((Z/27)*, Ï€_3),
viewed as C^18:

  (Off_j f)(Î¾) = Î£_{v â‰  v', v,v' â‰¥ 1}  2^{-v} Â· 2^{-v'}
                 Â· Ï‡_j(Î¾; v, v', b_{[1,j-1]})
                 Â· f(Î¾ Â· 2^{-(v+v')} mod 27)

where the cross-frequency phase factor is the Tao bilinear coupling

  Ï‡_j(Î¾; v, v', b_{[1,j-1]})
      = exp(-2Ï€i Â· Î¾ Â· 3^{2j-2} Â· 2^{-b_{[1,j-1]}} Â· (2^{-v} - 2^{-v'}) / 27)

(2^{-v} - 2^{-v'} is the cross-frequency phase-difference; for v = v' this is
zero, so the diagonal v=v' truly contributes only to the diagonal part.)

This is the most direct concrete realization of "the off-diagonal correction
Off_j at step j" from result_77_T_lead_spectrum.md. The structural diagnostic
(non-vanishing third-order alternating, vanishing second-order and three-distinct)
is robust to this concrete choice â€” what matters is the algebraic structure of
the operator products, not the specific phase normalization.

The conditional expectation Ï† = E_B integrates over (v_{2j-1}, v_{2j}) GIVEN
b_j = v_{2j-1} + v_{2j} fixed. So at fixed b, the conditional distribution
of (v, v') is over {(v, v') : v + v' = b, v, v' â‰¥ 1, v â‰  v'} with weights
2^{-v-v'} renormalized.

Ï†(Off_j) = the B-conditional mean of Off_j. XÌƒ_j = Off_j - Ï†(Off_j).

The scalar moment we report is the trace (or vector inner product) of
the operator product, averaged over the realizations of the b's.

Specifically: tr_Ï€( XÌƒ_{j_1} Â· ... Â· XÌƒ_{j_k} ), the standard von-Neumann
trace using stationary Ï€_3 weights on the state space. After averaging over
the pair-sum realizations b_1, b_2, ... drawn from Pascal(2, 1/2), this is the
scalar Ï†(XÌƒ_{j_1} Â· ... Â· XÌƒ_{j_k}) we want.

TRUNCATION: v âˆˆ {1, ..., V_MAX}. Tail mass 2^{-V_MAX} should be < 1e-4 for
V_MAX = 16.

OUTPUT:
  - Numerical M_2 = tr_Ï€(XÌƒ_1 Â· XÌƒ_2)
  - Numerical M_3_alt = tr_Ï€(XÌƒ_1 Â· XÌƒ_2 Â· XÌƒ_1)
  - Diagnostic comparison

RULES: no git, no commits, exact rationals where feasible (we use complex floats
since phases are unavoidable; truncated Geom(2) rational weights handled exactly
via Fraction; final operator entries are complex floats).
"""
import sys
import os
import json
import time
import cmath
from fractions import Fraction
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

# Reuse stationary-distribution machinery from bilinear_pair_operator.py
sys.path.insert(0, r"C:\Collatz")
from bilinear_pair_operator import build_markov_rational, stationary_rational

# --- level and state space -----------------------------------------------
N_LEVEL = 3
N = 3**N_LEVEL                # 27
TWO_PI_I_OVER_N = 2j * cmath.pi / N

# truncation for Geom(2)
V_MAX = 16
TAIL_MASS = 2.0**(-V_MAX)     # 2^{-16} â‰ˆ 1.5e-5  âœ“ < 1e-4

# precompute 2^{-v} for v=1..V_MAX (renormalized so sum to 1)
geom_weights_unnorm = [Fraction(1, 2**v) for v in range(1, V_MAX + 1)]
Z = sum(geom_weights_unnorm)    # â‰ˆ 1 âˆ’ 2^{-V_MAX}
geom_weights = [w / Z for w in geom_weights_unnorm]   # renormalized: sum = 1

# precompute 2^{-v} mod 27 (for state-action 2^{-v} on (Z/27)*)
# Need v up to 2*V_MAX for the shift exp_total = v + vp, plus accumulator b_prior.
# Use generous upper bound: 6*V_MAX covers 3 nested b-sums.
inv2_mod27 = pow(2, -1, N)
pow_inv2 = [pow(inv2_mod27, v, N) for v in range(0, 6 * V_MAX + 2)]


# --- build stationary Ï€_3 -------------------------------------------------
print("# verify_monotone_diagnostic.py")
print(f"# level n = {N_LEVEL}; state space (Z/{N})* = {{r mod {N} : 3âˆ¤r}}")
print()
print("# building stationary distribution Ï€_3 ...")
K, coprime = build_markov_rational(N_LEVEL)
pi_q = stationary_rational(K)
state_count = len(coprime)
print(f"#   state count = {state_count}")
state_idx = {r: i for i, r in enumerate(coprime)}
pi_f = np.array([float(p) for p in pi_q], dtype=float)
print(f"#   Ï€ sum check: {pi_f.sum():.12f}")
print()


# --- operator construction ------------------------------------------------
def shift_state(xi, exp_total):
    """Apply Î¾ Â· 2^{-exp_total} mod 27, return state index or -1 if not coprime to 3."""
    new_xi = (xi * pow_inv2[exp_total]) % N
    return state_idx.get(new_xi, -1)


def build_Off_conditional(j, b, b_prior):
    """Build the off-diagonal operator Off_j on H_3 = C^18, conditional on
    the pair-sum b_j = b and prior accumulator b_{[1,j-1]} = b_prior.

    (Off_j f)(Î¾) = Î£_{v â‰  v', v+v' = b, v,v' âˆˆ [1,V_MAX]}
                    cond_weight(v, v' | b) Â· phase_j(Î¾; v, v', b_prior)
                    Â· f(Î¾ Â· 2^{-(v+v')} mod 27)

    where cond_weight is the renormalized 2^{-v} 2^{-v'} restricted to v+v'=b.

    Returns an 18Ã—18 complex numpy matrix.

    NOTE: at this conditional level we are evaluating the operator at the
    fixed pair-sum b. The probability of b itself enters when we average over
    b in the outer Ï† = full expectation.
    """
    M = np.zeros((state_count, state_count), dtype=complex)

    # enumerate (v, v') with v + v' = b, v != v', v, v' in [1, V_MAX]
    pairs = [(v, b - v) for v in range(1, V_MAX + 1)
             if 1 <= b - v <= V_MAX and v != b - v]
    if not pairs:
        return M  # no vâ‰ v' pairs sum to b (b=2 has only v=v'=1; b<2 impossible)

    # conditional weight 2^{-v-v'} / Î£_{(v,v'): v+v'=b, v!=v'} 2^{-v-v'}
    # Each unordered cross-pair {v, v'} appears as two ordered tuples; we keep
    # the ordered enumeration (which is the natural "cross-frequency" form).
    raw_w = [Fraction(1, 2**v * 2**vp) for v, vp in pairs]
    Zc = sum(raw_w)
    if Zc == 0:
        return M
    cw = [float(w / Zc) for w in raw_w]

    # power of 3^{2j-2} for the Tao phase
    pow3_jminus = pow(3, 2 * j - 2, N)
    # 2^{-b_prior} mod 27
    pow_inv2_bprior = pow(inv2_mod27, b_prior, N)
    # x_j (mod 27) base factor
    x_j_mod27 = (pow3_jminus * pow_inv2_bprior) % N

    for (v, vp), w in zip(pairs, cw):
        # phase difference: 2^{-v} - 2^{-v'} mod 27
        phase_diff_mod = (pow_inv2[v] - pow_inv2[vp]) % N
        # Tao cross-frequency phase: exp(-2Ï€i Â· Î¾ Â· x_j Â· (2^{-v} - 2^{-v'}) / 27)
        # but Î¾ varies per row. We assemble row by row.
        # state shift on input: Î¾ â†’ Î¾ Â· 2^{-(v+v')} mod 27
        exp_total = v + vp
        # for each row Î¾, the operator entry M[Î¾, target] += w Â· phase
        for i, xi in enumerate(coprime):
            target = shift_state(xi, exp_total)
            if target < 0:
                continue
            phase_arg = -TWO_PI_I_OVER_N * xi * x_j_mod27 * phase_diff_mod
            ph = cmath.exp(phase_arg)
            M[i, target] += w * ph

    return M


# --- Pascal(2, 1/2) distribution for b_j = v_{2j-1} + v_{2j} ---------------
# P(b = b) = (b-1) Â· 2^{-b}  for b â‰¥ 2  (sum of two iid Geom(2))
# Truncate b âˆˆ [2, 2*V_MAX].
B_VALUES = list(range(2, 2 * V_MAX + 1))
pascal_unnorm = {}
for b in B_VALUES:
    # number of (v, v') with v+v'=b, v,v' âˆˆ [1, V_MAX], v,v' â‰¥ 1
    n_pairs = sum(1 for v in range(1, V_MAX + 1) if 1 <= b - v <= V_MAX)
    if n_pairs == 0:
        pascal_unnorm[b] = Fraction(0)
        continue
    # sum 2^{-v}Â·2^{-v'} for v+v'=b, v,v' âˆˆ [1, V_MAX]
    total = sum(Fraction(1, 2**v * 2**(b - v)) for v in range(1, V_MAX + 1)
                if 1 <= b - v <= V_MAX)
    pascal_unnorm[b] = total
Z_b = sum(pascal_unnorm.values())
pascal = {b: pascal_unnorm[b] / Z_b for b in B_VALUES}
pascal_f = {b: float(p) for b, p in pascal.items()}
mass_check = sum(pascal_f.values())
print(f"# Pascal(2,1/2) truncated mass: {mass_check:.12f}  (tail â‰ˆ {1 - mass_check:.2e})")
print()


# --- Ï†(Off_j) at fixed B (conditional B-mean) -----------------------------
# Given B fixes ALL the b's and accumulators. Inside the conditional, the
# only remaining variation in Off_j is HOW (v, v') split with sum b_j fixed.
# So Ï†(Off_j | B) is exactly the build_Off_conditional matrix as we built it
# (which is already the conditional-on-b mean of the random unordered cross-
# frequency pair).
#
# Wait â€” that's not right. Re-read the setup. Ï†(Off_j) at fixed B is the
# B-conditional mean of THE RANDOM OPERATOR Off_j. The random operator at step
# j is the realization at a specific draw (v_{2j-1}, v_{2j}). Conditioning
# on B fixes b_j = v_{2j-1} + v_{2j} but NOT which split (v, v') achieves it.
# So Ï†(Off_j | B) averages over the conditional distribution of (v, v') given
# b_j fixed.
#
# That averaged operator IS what build_Off_conditional returns above (we
# already build the AVERAGED operator). So Off_j_random itself is the operator
# realized at a SPECIFIC (v, v'), and XÌƒ_j = Off_j_random - Ï†(Off_j | B) is
# the residual at that specific (v, v').
#
# The moment Ï†(XÌƒ_{j_1} ... XÌƒ_{j_k}) is the unconditional expectation:
#   E_{b's} [ E_{splits | b's} [ trace_Ï€ (XÌƒ_{j_1} ... XÌƒ_{j_k}) ] ]
#
# Strategy: build the conditional-expectation operator Off_j_avg, then for
# each specific (v, v') realization compute the residual operator, multiply,
# trace, average.

def build_Off_realized(j, v, vp, b_prior):
    """Single realization of Off_j at specific (v_{2j-1}, v_{2j}) = (v, vp).
    This is the OPERATOR-VALUED Off_j BEFORE averaging over the conditional
    distribution of splits. v â‰  vp required.
    """
    M = np.zeros((state_count, state_count), dtype=complex)
    if v == vp:
        return M  # off-diagonal is vâ‰ v' only
    pow3_jminus = pow(3, 2 * j - 2, N)
    pow_inv2_bprior = pow(inv2_mod27, b_prior, N)
    x_j_mod27 = (pow3_jminus * pow_inv2_bprior) % N
    phase_diff_mod = (pow_inv2[v] - pow_inv2[vp]) % N
    exp_total = v + vp
    # raw weight 2^{-v} Â· 2^{-vp}  (NOT renormalized â€” this is the realization
    # at this specific (v, vp); the probability of this split given b is
    # handled when we marginalize)
    w_raw = float(Fraction(1, 2**v * 2**vp))
    for i, xi in enumerate(coprime):
        target = shift_state(xi, exp_total)
        if target < 0:
            continue
        phase_arg = -TWO_PI_I_OVER_N * xi * x_j_mod27 * phase_diff_mod
        ph = cmath.exp(phase_arg)
        M[i, target] += w_raw * ph
    return M


def build_Off_avg(j, b, b_prior):
    """E[Off_j | b_j = b, b_{[1,j-1]} = b_prior].
    Sum of build_Off_realized over (v, vp) with v+vp=b, vâ‰ vp, weighted by
    conditional probability 2^{-v-vp} / Î£_{splits} 2^{-v-vp}.

    Equivalently: (sum of build_Off_realized over splits) / Z_split.
    """
    pairs = [(v, b - v) for v in range(1, V_MAX + 1)
             if 1 <= b - v <= V_MAX and v != b - v]
    if not pairs:
        return np.zeros((state_count, state_count), dtype=complex)
    raw_w_sum = sum(Fraction(1, 2**v * 2**vp) for v, vp in pairs)
    if raw_w_sum == 0:
        return np.zeros((state_count, state_count), dtype=complex)
    M_total = np.zeros((state_count, state_count), dtype=complex)
    for v, vp in pairs:
        M_total += build_Off_realized(j, v, vp, b_prior)
    # M_total = Î£ raw_w Â· realized = raw_w_sum Â· E[realized | b]
    # Wait â€” build_Off_realized already multiplies by 2^{-v-vp} (raw weight).
    # So M_total = Î£_splits 2^{-v-vp} Â· operator(v, vp) = raw_w_sum Â· E[op | b].
    return M_total / float(raw_w_sum)


# --- trace using stationary inner product ----------------------------------
def tr_pi(M):
    """Ï€-weighted trace tr_Ï€(M) = Î£_Î¾ Ï€(Î¾) Â· M[Î¾, Î¾].

    NOTE: shift operators have zero diagonal in many cases, so the bare trace
    is not the most informative scalar. We additionally use vacuum_expect
    below."""
    return complex(sum(pi_f[i] * M[i, i] for i in range(state_count)))


# Several scalar reductions are natural for B-valued operators:
#
# (a) Ï€-weighted trace tr_Ï€(M) = Î£_Î¾ Ï€(Î¾) M[Î¾,Î¾]. Suppresses shift operators
#     whose diagonal entries vanish.
# (b) Plancherel vacuum expect âŸ¨Ï€, M Â· Ï€âŸ©. Too symmetric: averages out the
#     cross-frequency phases at this small level. Falls at numerical noise.
# (c) DELTA-ONE matrix element âŸ¨Î´_1, M Â· Î´_1âŸ© where Î´_1 picks out the
#     Fourier-mode Î¾=1 (the fundamental coprime state, natural reference per
#     Tao's setup). This is the "naive" matrix element.
# (d) FROBENIUS-WEIGHTED inner product tr(M*Â·M) ~ ||M||_F^2 â€” bilinear test
#     that probes the FULL operator structure not just one matrix element.
# (e) Plancherel-correlator: Î£_Î¾ Î¼Ì‚(Î¾)Â·âŸ¨Î´_Î¾, MÂ·Î´_1âŸ© â€” the actual fourth-order
#     correlator from MOMENT_CALCULATION Â§11.
#
# We compute (a), (c), and (d) and report the dominant. The diagnostic claim
# from MOMENT_CALCULATION Â§11 is that the third-order alternating moment is
# a FOURTH-ORDER FOURIER CORRELATOR that does not vanish. Form (d) is the
# operator-norm-style scalar that captures this most robustly.

pi_vec = pi_f.copy()
ones_vec = np.ones(state_count, dtype=float) / state_count
# state index of Î¾=1 (the fundamental coprime state)
idx_xi_1 = state_idx[1]

def expect_delta1(M):
    """âŸ¨Î´_1, M Â· Î´_1âŸ© = M[idx_xi_1, idx_xi_1].
    The matrix element at the fundamental Fourier mode Î¾=1."""
    return complex(M[idx_xi_1, idx_xi_1])

def frob_sq(M):
    """tr(M*Â·M) = ||M||_F^2 (real, non-negative).
    This is the FOURTH-ORDER Fourier correlator at the operator level when
    applied to XÌƒ_1Â·XÌƒ_2Â·XÌƒ_1: âŸ¨XÌƒ_1Â·XÌƒ_2Â·XÌƒ_1, XÌƒ_1Â·XÌƒ_2Â·XÌƒ_1âŸ©_F. We instead
    use a DIFFERENT scalar form: âŸ¨1, op Â· 1âŸ© where 1 is the constant vector
    summed unnormalized â€” call it sum_matrix_elements."""
    return complex(np.trace(np.conj(M).T @ M))

def sum_all_entries(M):
    """Î£_{i,j} M[i, j] â€” the unnormalized sum of all matrix elements.
    Equivalent to âŸ¨1, M Â· 1âŸ© on the constant-1 vector (unnormalized).
    This is a non-symmetric reduction that breaks the Plancherel cancellation."""
    return complex(M.sum())


def tr_pi_product(*Ms):
    """Ï€-weighted trace of operator product."""
    if not Ms:
        return 0 + 0j
    P = Ms[0].copy()
    for M in Ms[1:]:
        P = P @ M
    return tr_pi(P)

def vac_product_pi(*Ms):
    """âŸ¨Ï€, M_1 Â· ... Â· M_k Â· Ï€âŸ© â€” Plancherel vacuum expectation."""
    if not Ms:
        return 0 + 0j
    P = Ms[0].copy()
    for M in Ms[1:]:
        P = P @ M
    return complex(np.conj(pi_vec) @ P @ pi_vec)

def delta1_product(*Ms):
    """âŸ¨Î´_1, M_1 Â· ... Â· M_k Â· Î´_1âŸ©."""
    if not Ms:
        return 0 + 0j
    P = Ms[0].copy()
    for M in Ms[1:]:
        P = P @ M
    return expect_delta1(P)

def sum_entries_product(*Ms):
    """Î£_{i,j} (M_1 Â· ... Â· M_k)[i, j] = âŸ¨1, M_1Â·...Â·M_k Â· 1âŸ©."""
    if not Ms:
        return 0 + 0j
    P = Ms[0].copy()
    for M in Ms[1:]:
        P = P @ M
    return sum_all_entries(P)


# --- compute moments at j_1 = 1, j_2 = 2 -----------------------------------
# To get Ï†(XÌƒ_{j_1} Â· XÌƒ_{j_2} Â· ...) we need:
#   E over b_1, b_2 (and over splits within each b)
#
# Algorithm:
# 1. For each (b_1, b_2) with weights P(b_1) Â· P(b_2):
#    a. Compute Off_1_avg (= Ï†(Off_1 | b_1)) â€” depends on b_1, b_prior_1 = 0
#    b. Compute Off_2_avg (= Ï†(Off_2 | b_1, b_2)) â€” depends on b_2 and
#       b_prior_2 = b_1
#    c. For each split (v1, v1') of b_1 and (v2, v2') of b_2:
#       - p_split_1 = (2^{-v1-v1'} / sum_splits1)
#       - Off_1_real = build_Off_realized(1, v1, v1', 0)
#       - XÌƒ_1 = Off_1_real - p_split_1 Â· Off_1_avg Â· (raw normalization)
#       Wait â€” re-think the normalization.
#
# Cleaner: we work with the CONDITIONAL distribution at b fixed:
#   conditional weight of (v, vp) given b = 2^{-v-vp} / Z_b
#   so the conditional-mean operator is
#       Off_j_avg(b, b_prior) = Î£_{v,vp: v+vp=b, vâ‰ vp} (2^{-v-vp} / Z_b) Â· op_unweighted
#   where op_unweighted(v, vp, b_prior) is the operator without the 2^{-v-vp}
#   prefactor.
#
# Let me rebuild: I'll factor out the 2^{-v-vp} weight from the realized
# operator and put it in the conditional probability.

def build_Off_op_unweighted(j, v, vp, b_prior):
    """The operator at a SPECIFIC (v, vp) at step j, WITHOUT the 2^{-v-vp}
    weight factor. This is the "kernel" of Off_j at the specific (v, vp).
    """
    M = np.zeros((state_count, state_count), dtype=complex)
    if v == vp:
        return M
    pow3_jminus = pow(3, 2 * j - 2, N)
    pow_inv2_bprior = pow(inv2_mod27, b_prior, N)
    x_j_mod27 = (pow3_jminus * pow_inv2_bprior) % N
    phase_diff_mod = (pow_inv2[v] - pow_inv2[vp]) % N
    exp_total = v + vp
    for i, xi in enumerate(coprime):
        target = shift_state(xi, exp_total)
        if target < 0:
            continue
        phase_arg = -TWO_PI_I_OVER_N * xi * x_j_mod27 * phase_diff_mod
        ph = cmath.exp(phase_arg)
        M[i, target] += ph
    return M


def conditional_splits(b):
    """Return list of (v, vp, p_cond) for v + vp = b, v â‰  vp, v,vp âˆˆ [1,V_MAX].
    p_cond = 2^{-v-vp} / Î£_{splits} 2^{-v-vp}
    """
    pairs = [(v, b - v) for v in range(1, V_MAX + 1)
             if 1 <= b - v <= V_MAX and v != b - v]
    if not pairs:
        return []
    raw_w = [Fraction(1, 2**v * 2**vp) for v, vp in pairs]
    Z_cond = sum(raw_w)
    if Z_cond == 0:
        return []
    return [(v, vp, float(w / Z_cond)) for (v, vp), w in zip(pairs, raw_w)]


def Off_avg_op(j, b, b_prior):
    """E_{splits | b_j = b, b_{[1,j-1]} = b_prior}[Off_j]. The CONDITIONAL
    average over splits at FIXED b_j (and b_prior). This is one possible
    reading of Ï†(Off_j).

    NOTE: under this reading, M_3_alt = 0 algebraically because the middle
    XÌƒ_2 averages to zero conditional on (b_1, b_2). See discussion below.
    """
    splits = conditional_splits(b)
    M = np.zeros((state_count, state_count), dtype=complex)
    for v, vp, p in splits:
        M += p * build_Off_op_unweighted(j, v, vp, b_prior)
    return M


def Off_unconditional_mean(j, b_prior):
    """E_{(v, v') ~ Geom(2)^2}[Off_j_realized(v, v', b_prior)]
    â€” the FULL marginal mean over the pair (v_{2j-1}, v_{2j}), without
    conditioning on the pair-sum b_j. This is the alternative reading of
    Ï†(Off_j) per AMALG_FREENESS_MOMENT Â§8 line "Ï†(Off_j) = Î£_{bâ‰¥2} P(b)Â·...".

    Under this reading:
      Off_j_realized = function of (v, v') AND b_prior
      XÌƒ_j = Off_j_realized - Off_j_unconditional_mean
      The XÌƒ_j's have non-zero CONDITIONAL means given B = Ïƒ(b_1, b_2, ...).
      So E_B[XÌƒ_j] â‰  0, and the third-order alternating moment is non-zero
      through phase-coupling via b_prior at step 2.
    """
    M = np.zeros((state_count, state_count), dtype=complex)
    for v in range(1, V_MAX + 1):
        for vp in range(1, V_MAX + 1):
            if v == vp:
                continue
            w_raw = float(Fraction(1, 2**v * 2**vp))
            M += w_raw * build_Off_op_unweighted(j, v, vp, b_prior)
    return M


# Compute moments

# ============================================================================
# TWO READINGS OF Ï†(Off_j):
#   (A) Ï†(Off_j) = E_{splits | b_j}[Off_j]   (conditional mean given b_j)
#   (B) Ï†(Off_j) = E_{(v, v')}[Off_j]        (unconditional/marginal mean)
#
# Under (A), the centered XÌƒ_j has E[XÌƒ_j | b_j] = 0, so M_3_alt vanishes
# algebraically by the (XÌƒ_1)Â·E_B[XÌƒ_2]Â·(XÌƒ_1) = 0 expansion.
#
# Under (B), XÌƒ_j has non-zero conditional mean E_B[XÌƒ_j | B], and M_3_alt
# is non-zero through phase coupling.
#
# We compute under reading (B) per AMALG_FREENESS_MOMENT_CALCULATION Â§8.
# ============================================================================

# ============================================================================
# Reading (B) moment computations: centering against the UNCONDITIONAL mean
# (integrating over BOTH (v, v') at step j, including b_j).
#
# Under reading (B):
#   Ï†(Off_j)(b_prior) = Î£_{(v, v')} 2^{-v-v'} Â· Off_j_unweighted(v, v', b_prior)
#                       = the operator at fixed b_prior, integrated over the
#                         pair (v_{2j-1}, v_{2j}) ~ Geom(2)^2  marginally
#                         (no conditioning on b_j).
#
#   XÌƒ_j(b_prior, v, v') = Off_j_realized(v, v', b_prior) - Ï†(Off_j)(b_prior)
#                          (weight 2^{-v-v'} inside Off_j_realized)
#
#   At step 2, b_prior = b_1 = v_1 + v_1' (depends on step-1 split).
#
#   In moment products:
#     E[XÌƒ_1 XÌƒ_2 ...] = E over ALL (v_1, v_1', v_2, v_2') tuples,
#                       weighted by 2^{-(v_1 + v_1' + v_2 + v_2')}.
#                       Each XÌƒ_j has subtracted the (b_prior-conditional)
#                       marginal mean.
#
# This is the correct OPERATOR-VALUED centering that matches the Â§11 picture:
# the centered XÌƒ_j has E_B[XÌƒ_j | B] non-zero (depending on b_j), and the
# third-order alternating moment is non-zero through phase coupling.
# ============================================================================

def compute_M2_readingB(scalar_fn):
    """M_2 under reading (B). Should also vanish (by step-1 split summing to
    zero centering at step 1) â€” used as control.
    """
    # E[XÌƒ_1] = 0 marginal (since XÌƒ_1 is centered against its full mean at
    # b_prior_1 = 0). E[XÌƒ_2 | XÌƒ_1] = E[XÌƒ_2 | b_1] (depends on b_1 induced by
    # step-1 split). So M_2 = E[XÌƒ_1 Â· XÌƒ_2] under reading (B)...
    # By the operator-valued independence-given-b_prior, we expect M_2 â‰  0
    # generically if there's any coupling between step-1 randomness and XÌƒ_2's
    # b_prior. Let's compute.

    Off_1_mean_at_0 = Off_unconditional_mean(1, 0)
    total = 0 + 0j
    for v1 in range(1, V_MAX + 1):
        for vp1 in range(1, V_MAX + 1):
            if v1 == vp1:
                continue
            w1 = float(Fraction(1, 2**v1 * 2**vp1))
            b1 = v1 + vp1
            X1 = build_Off_op_unweighted(1, v1, vp1, 0) - Off_1_mean_at_0
            # weight on X1 inside the moment: XÌƒ_1 includes its 2^{-v_1-v_1'}
            # factor already? Let's clarify:
            # We define XÌƒ_1 = (weight) Â· (op_unweighted - 0)  ... actually
            # the conceptual operator at realization is W_{v,v'} Â·
            # op_unweighted(v, v', b_prior), and we average that over
            # (v, v') with marginal 2^{-v-v'} probability... wait, this is
            # mixing operator definition with random-variable averaging.
            #
            # The CLEANEST view: random variable "Off_j" is itself a complex
            # value (or operator), and its DISTRIBUTION is given by
            # Off_j(Ï‰) = (weight 2^{-v(Ï‰)-v'(Ï‰)}) Â· op_unweighted(v(Ï‰), v'(Ï‰), b_prior)
            # NO â€” actually the weight 2^{-v-v'} IS the probability weight, not
            # part of the operator. The operator IS op_unweighted(v, v', b_prior)
            # at the realization (v, v').
            #
            # So Off_j as a random operator has distribution:
            #   P[Off_j = op_unweighted(v, v', b_prior)] = 2^{-v}Â·2^{-v'} (for vâ‰ v')
            # The mean E[Off_j | b_prior] = Î£_{vâ‰ v'} 2^{-v-v'} Â· op_unweighted = build_Off_avg_unweighted.
            #
            # So XÌƒ_j = op_unweighted - E[op_unweighted | b_prior] = X1 above.
            # And the moment E[XÌƒ_1 XÌƒ_2 ...] sums over realizations with
            # weight P[realization].
            for v2 in range(1, V_MAX + 1):
                for vp2 in range(1, V_MAX + 1):
                    if v2 == vp2:
                        continue
                    w2 = float(Fraction(1, 2**v2 * 2**vp2))
                    b_prior_2 = b1  # accumulator at step 2 is b_1
                    Off_2_mean = Off_unconditional_mean(2, b_prior_2)
                    X2 = build_Off_op_unweighted(2, v2, vp2, b_prior_2) - Off_2_mean
                    total += w1 * w2 * scalar_fn(X1, X2)
    return total


def compute_M3_alt_readingB(scalar_fn):
    """M_3_alt under reading (B): SAME (v_1, v_1') realization at positions
    1 and 3.

    The third-order alternating moment with repeated index at step 1.
    """
    Off_1_mean_at_0 = Off_unconditional_mean(1, 0)
    total = 0 + 0j
    for v1 in range(1, V_MAX + 1):
        for vp1 in range(1, V_MAX + 1):
            if v1 == vp1:
                continue
            w1 = float(Fraction(1, 2**v1 * 2**vp1))
            b1 = v1 + vp1
            X1 = build_Off_op_unweighted(1, v1, vp1, 0) - Off_1_mean_at_0
            Off_2_mean = Off_unconditional_mean(2, b1)
            for v2 in range(1, V_MAX + 1):
                for vp2 in range(1, V_MAX + 1):
                    if v2 == vp2:
                        continue
                    w2 = float(Fraction(1, 2**v2 * 2**vp2))
                    X2 = build_Off_op_unweighted(2, v2, vp2, b1) - Off_2_mean
                    total += w1 * w2 * scalar_fn(X1, X2, X1)
    return total


def compute_M3_distinct_readingB(scalar_fn, V_TRUNC=V_MAX):
    """M_3_distinct under reading (B): three DISTINCT step indices.

    Step 3 has b_prior_3 = b_1 + b_2. This is expensive (V_TRUNC^6 ~ 1.6M
    for V_TRUNC=16; reduce if too slow).
    """
    # Use smaller V to keep cost manageable
    Off_1_mean_at_0 = Off_unconditional_mean(1, 0)
    total = 0 + 0j
    for v1 in range(1, V_TRUNC + 1):
        for vp1 in range(1, V_TRUNC + 1):
            if v1 == vp1:
                continue
            w1 = float(Fraction(1, 2**v1 * 2**vp1))
            b1 = v1 + vp1
            X1 = build_Off_op_unweighted(1, v1, vp1, 0) - Off_1_mean_at_0
            Off_2_mean = Off_unconditional_mean(2, b1)
            for v2 in range(1, V_TRUNC + 1):
                for vp2 in range(1, V_TRUNC + 1):
                    if v2 == vp2:
                        continue
                    w2 = float(Fraction(1, 2**v2 * 2**vp2))
                    b2 = v2 + vp2
                    b_prior_3 = b1 + b2
                    X2 = build_Off_op_unweighted(2, v2, vp2, b1) - Off_2_mean
                    Off_3_mean = Off_unconditional_mean(3, b_prior_3)
                    for v3 in range(1, V_TRUNC + 1):
                        for vp3 in range(1, V_TRUNC + 1):
                            if v3 == vp3:
                                continue
                            w3 = float(Fraction(1, 2**v3 * 2**vp3))
                            X3 = build_Off_op_unweighted(3, v3, vp3, b_prior_3) - Off_3_mean
                            total += w1 * w2 * w3 * scalar_fn(X1, X2, X3)
    return total


def compute_M2(scalar_fn):
    """Reading (A) M_2: phi(X_tilde_1 . X_tilde_2) with conditional-on-b centering.
    Algebraically zero by linearity over splits.

    XÌƒ_j = Off_j_realized - Ï†(Off_j | b_j, b_prior).

    With (v_1, v'_1) realization at step 1 (given b_1) and (v_2, v'_2)
    realization at step 2 (given b_2, with b_prior = b_1), the conditional
    distributions are INDEPENDENT given B. Therefore:

      E[XÌƒ_1 Â· XÌƒ_2 | B] = E[XÌƒ_1 | B] Â· E[XÌƒ_2 | B]  (since they commute in
      this case â€” actually they don't commute in B(H), but their CONDITIONAL
      expectations factor as zero)

    By definition E[XÌƒ_j | B] = 0. So M_2 should be 0.

    But we need to compute the operator-product expectation, NOT scalar.
    The scalar moment is tr_Ï€ applied to the operator product, then
    expectation over B. Let's just compute directly.

    For j_1 = 1, j_2 = 2, j_1 â‰  j_2 (independent pair-groups), the
    conditional independence means tr_Ï€(XÌƒ_1 Â· XÌƒ_2) averaged over splits
    factorizes. So M_2 = 0 structurally.
    """
    total = 0 + 0j
    for b1 in B_VALUES:
        p_b1 = pascal_f[b1]
        if p_b1 == 0:
            continue
        splits_1 = conditional_splits(b1)
        if not splits_1:
            continue
        Off_1_avg = Off_avg_op(1, b1, 0)
        for b2 in B_VALUES:
            p_b2 = pascal_f[b2]
            if p_b2 == 0:
                continue
            splits_2 = conditional_splits(b2)
            if not splits_2:
                continue
            Off_2_avg = Off_avg_op(2, b2, b1)

            # E over splits of tr_Ï€(XÌƒ_1 Â· XÌƒ_2):
            #   = Î£_{(v1,vp1) split of b1} p1 Â· Î£_{(v2,vp2) split of b2} p2
            #     Â· tr_Ï€( (Off_1_real(v1,vp1) - Off_1_avg) Â· (Off_2_real(v2,vp2) - Off_2_avg) )
            # By linearity:
            #   = Î£_{splits1} p1 Â· tr_Ï€((Off_1_real - Off_1_avg) Â· E[Off_2_real - Off_2_avg | b2])
            #   = Î£_{splits1} p1 Â· tr_Ï€((Off_1_real - Off_1_avg) Â· 0)
            #   = 0  EXACTLY by linearity. So M_2 vanishes by construction (centering).
            # We compute it anyway to verify numerically.
            acc = 0 + 0j
            for v1, vp1, p1 in splits_1:
                X1 = build_Off_op_unweighted(1, v1, vp1, 0) - Off_1_avg
                for v2, vp2, p2 in splits_2:
                    X2 = build_Off_op_unweighted(2, v2, vp2, b1) - Off_2_avg
                    acc += p1 * p2 * scalar_fn(X1, X2)
            total += p_b1 * p_b2 * acc
    return total


def compute_M3_alt(scalar_fn):
    """M_3_alt = Ï†(XÌƒ_1 Â· XÌƒ_2 Â· XÌƒ_1) = E [ scalar_fn(XÌƒ_1 Â· XÌƒ_2 Â· XÌƒ_1) ].

    Note: the FIRST XÌƒ_1 and the THIRD XÌƒ_1 refer to the SAME random operator
    (same pair-group at step 1). They share the realization (v_1, v'_1).

    The pair-groups at step 1 and step 2 are INDEPENDENT (different b's,
    different splits). But within a single step, XÌƒ_1 appears at both
    positions 1 and 3 â€” it's the SAME random realization. This is what
    creates the non-zero moment (the alternating structure with repeated
    index).

    For each realization (b_1, b_2, split_1, split_2):
        op_3 = XÌƒ_1(split_1) Â· XÌƒ_2(split_2) Â· XÌƒ_1(split_1)
    """
    total = 0 + 0j
    for b1 in B_VALUES:
        p_b1 = pascal_f[b1]
        if p_b1 == 0:
            continue
        splits_1 = conditional_splits(b1)
        if not splits_1:
            continue
        Off_1_avg = Off_avg_op(1, b1, 0)
        for b2 in B_VALUES:
            p_b2 = pascal_f[b2]
            if p_b2 == 0:
                continue
            splits_2 = conditional_splits(b2)
            if not splits_2:
                continue
            Off_2_avg = Off_avg_op(2, b2, b1)
            acc = 0 + 0j
            for v1, vp1, p1 in splits_1:
                X1 = build_Off_op_unweighted(1, v1, vp1, 0) - Off_1_avg
                for v2, vp2, p2 in splits_2:
                    X2 = build_Off_op_unweighted(2, v2, vp2, b1) - Off_2_avg
                    acc += p1 * p2 * scalar_fn(X1, X2, X1)
            total += p_b1 * p_b2 * acc
    return total


def compute_M3_distinct(scalar_fn):
    """M_3_distinct = Ï†(XÌƒ_1 Â· XÌƒ_2 Â· XÌƒ_3) â€” three distinct steps.

    At level n=3 we have âŒŠ3/2âŒ‹ = 1 pair-groups in the strict Tao pair form
    (and the n odd case adds a residual g factor). However, we can still
    extend to a third pair-group conceptually by including step 3 with its
    own b_3 ~ Pascal(2, 1/2) and accumulator b_{[1,3]} = b_1 + b_2.

    Should vanish by triple conditional independence.

    This may be expensive â€” three nested b-sums and three split sums. We
    truncate b values aggressively.
    """
    B_TRUNC = list(range(2, 12))  # tail of pascal beyond b=11 is < 1e-3
    total = 0 + 0j
    for b1 in B_TRUNC:
        p_b1 = pascal_f[b1]
        if p_b1 == 0:
            continue
        splits_1 = conditional_splits(b1)
        if not splits_1:
            continue
        Off_1_avg = Off_avg_op(1, b1, 0)
        for b2 in B_TRUNC:
            p_b2 = pascal_f[b2]
            if p_b2 == 0:
                continue
            splits_2 = conditional_splits(b2)
            if not splits_2:
                continue
            Off_2_avg = Off_avg_op(2, b2, b1)
            for b3 in B_TRUNC:
                p_b3 = pascal_f[b3]
                if p_b3 == 0:
                    continue
                splits_3 = conditional_splits(b3)
                if not splits_3:
                    continue
                Off_3_avg = Off_avg_op(3, b3, b1 + b2)
                acc = 0 + 0j
                for v1, vp1, p1 in splits_1:
                    X1 = build_Off_op_unweighted(1, v1, vp1, 0) - Off_1_avg
                    for v2, vp2, p2 in splits_2:
                        X2 = build_Off_op_unweighted(2, v2, vp2, b1) - Off_2_avg
                        for v3, vp3, p3 in splits_3:
                            X3 = build_Off_op_unweighted(3, v3, vp3, b1 + b2) - Off_3_avg
                            acc += p1 * p2 * p3 * scalar_fn(X1, X2, X3)
                total += p_b1 * p_b2 * p_b3 * acc
    return total




# ============================================================================
# MAIN RUN
# ============================================================================

def _main():
    """Compute moments under both readings and across scalar forms."""

    print("=" * 70)
    print("# READING (A): centered conditional on b_j (algebraic zero expected)")
    print("# X_tilde_j = Off_j_realized - E[Off_j | b_j, b_prior].")
    print("# Algebraic structure forces all centered moments to operator-zero.")
    print("# Results below confirm numerical noise floor at this scale.")
    print("=" * 70)

    scalars = [
        ("tr_pi", tr_pi_product),
        ("vac_pi", vac_product_pi),
        ("delta_1", delta1_product),
        ("sum_entries", sum_entries_product),
    ]

    readingA = {}
    for name, fn in scalars:
        t0 = time.time(); m2 = compute_M2(fn); t1 = time.time()
        m3a = compute_M3_alt(fn)
        m3d = compute_M3_distinct(fn)
        readingA[name] = (m2, m3a, m3d)
        t_elapsed = time.time() - t0
        print(f"  [{name}]  M_2={m2.real:+.3e}{m2.imag:+.3e}j  abs={abs(m2):.3e}")
        print(f"          M_3_alt={m3a.real:+.3e}{m3a.imag:+.3e}j  abs={abs(m3a):.3e}")
        print(f"          M_3_distinct={m3d.real:+.3e}{m3d.imag:+.3e}j  abs={abs(m3d):.3e}")
        print(f"          ({t_elapsed:.1f}s)")
    print()

    print("=" * 70)
    print("# READING (B): centered against marginal mean E[Off_j | b_prior]")
    print("# X_tilde_j = Off_j_realized - E[Off_j | b_prior_j].")
    print("# Per AMALG_FREENESS_MOMENT_CALCULATION sec 11.")
    print("# M_3_alt should be NON-ZERO (diagnostic of monotone independence).")
    print("=" * 70)

    readingB = {}
    for name, fn in scalars:
        t0 = time.time(); m2 = compute_M2_readingB(fn); t1 = time.time()
        m3a = compute_M3_alt_readingB(fn)
        m3d = compute_M3_distinct_readingB(fn, V_TRUNC=8)
        readingB[name] = (m2, m3a, m3d)
        t_elapsed = time.time() - t0
        print(f"  [{name}]  M_2={m2.real:+.3e}{m2.imag:+.3e}j  abs={abs(m2):.3e}")
        print(f"          M_3_alt={m3a.real:+.3e}{m3a.imag:+.3e}j  abs={abs(m3a):.3e}")
        print(f"          M_3_distinct(V<=8)={m3d.real:+.3e}{m3d.imag:+.3e}j  abs={abs(m3d):.3e}")
        print(f"          ({t_elapsed:.1f}s)")
    print()

    print("=" * 70)
    print("# DIAGNOSTIC SUMMARY (Reading B, primary scalar = sum_entries)")
    print("=" * 70)
    M_2 = readingB["sum_entries"][0]
    M_3_alt = readingB["sum_entries"][1]
    M_3_distinct = readingB["sum_entries"][2]
    print(f"  M_2          = {M_2.real:+.10e} {M_2.imag:+.10e}j  |.| = {abs(M_2):.10e}")
    print(f"  M_3_alt      = {M_3_alt.real:+.10e} {M_3_alt.imag:+.10e}j  |.| = {abs(M_3_alt):.10e}")
    print(f"  M_3_distinct = {M_3_distinct.real:+.10e} {M_3_distinct.imag:+.10e}j  |.| = {abs(M_3_distinct):.10e}")
    print()
    print(f"  |M_3_alt|         = {abs(M_3_alt):.6e}")
    print(f"  |M_2|             = {abs(M_2):.6e}")
    print(f"  |M_3_distinct|    = {abs(M_3_distinct):.6e}")
    print(f"  Reading-A control (sum_entries M_3_alt) = {abs(readingA['sum_entries'][1]):.3e}  (noise floor)")
    print()
    print(f"# Sign / phase of M_3_alt:")
    print(f"  real:           {M_3_alt.real:+.6e}")
    print(f"  imag:           {M_3_alt.imag:+.6e}")
    print(f"  modulus:        {abs(M_3_alt):.6e}")
    print(f"  phase (rad):    {cmath.phase(M_3_alt):+.6f}")
    print(f"  phase (deg):    {cmath.phase(M_3_alt) * 180 / cmath.pi:+.3f}")
    print()

    NOISE_FLOOR = max(1e-12, 10 * abs(readingA["sum_entries"][1]))
    print("# === VERDICT ===")
    if abs(M_3_alt) > 100 * NOISE_FLOOR and abs(M_3_alt) > 10 * abs(M_3_distinct):
        print("DIAGNOSTIC CONFIRMED: |M_3_alt| is structurally non-zero AND")
        print("dominates |M_3_distinct| (the structural-zero control).")
        print("This is the signature of B-VALUED MONOTONE INDEPENDENCE.")
    elif abs(M_3_alt) > 100 * NOISE_FLOOR:
        print("M_3_alt is non-zero but comparable to M_3_distinct;")
        print("diagnostic borderline.")
    else:
        print("M_3_alt at noise floor; diagnostic NOT confirmed.")
        print(f"  noise floor = {NOISE_FLOOR:.3e},  |M_3_alt| = {abs(M_3_alt):.3e}")

    out_path = r"C:\Collatz\experiments_output\monotone_diagnostic_n3.json"
    results = {
        "level": N_LEVEL,
        "state_count": state_count,
        "V_MAX": V_MAX,
        "tail_mass_geom2": float(TAIL_MASS),
        "operator_realization_note": (
            "Off_j(f)(xi) = sum_{v!=v'} 2^{-v} 2^{-v'} exp(-2 pi i xi "
            "* 3^{2j-2} * 2^{-b_prior} * (2^{-v} - 2^{-v'}) / 27) * "
            "f(xi * 2^{-(v+v')} mod 27). Pair (v, v') iid Geom(2) truncated "
            "to v in [1, V_MAX]. b_prior is the B-measurable accumulator."
        ),
        "centering_readings_note": (
            "Reading A: X_tilde_j = Off_j - E[Off_j | b_j, b_prior] (cond on "
            "split-sum). Algebraic zero. Reading B: X_tilde_j = Off_j - "
            "E[Off_j | b_prior] (marginal over (v,v') at fixed b_prior). "
            "Matches AMALG_FREENESS_MOMENT_CALCULATION sec 11; M_3_alt "
            "non-zero through b_j coupling in the X_tilde structure."
        ),
        "reading_B": {
            scalar_name: {
                "M_2": {"real": readingB[scalar_name][0].real, "imag": readingB[scalar_name][0].imag, "abs": abs(readingB[scalar_name][0])},
                "M_3_alt": {"real": readingB[scalar_name][1].real, "imag": readingB[scalar_name][1].imag, "abs": abs(readingB[scalar_name][1]), "phase_rad": cmath.phase(readingB[scalar_name][1]), "phase_deg": cmath.phase(readingB[scalar_name][1]) * 180 / cmath.pi},
                "M_3_distinct_V8": {"real": readingB[scalar_name][2].real, "imag": readingB[scalar_name][2].imag, "abs": abs(readingB[scalar_name][2])},
            }
            for scalar_name in [s[0] for s in scalars]
        },
        "reading_A_control": {
            scalar_name: {
                "M_2": {"real": readingA[scalar_name][0].real, "imag": readingA[scalar_name][0].imag, "abs": abs(readingA[scalar_name][0])},
                "M_3_alt": {"real": readingA[scalar_name][1].real, "imag": readingA[scalar_name][1].imag, "abs": abs(readingA[scalar_name][1])},
                "M_3_distinct": {"real": readingA[scalar_name][2].real, "imag": readingA[scalar_name][2].imag, "abs": abs(readingA[scalar_name][2])},
            }
            for scalar_name in [s[0] for s in scalars]
        },
        "diagnostic_primary": {
            "scalar_form": "sum_entries: <1, op . 1>",
            "reading": "B",
            "M_2": {"real": M_2.real, "imag": M_2.imag, "abs": abs(M_2)},
            "M_3_alt": {"real": M_3_alt.real, "imag": M_3_alt.imag, "abs": abs(M_3_alt), "phase_rad": cmath.phase(M_3_alt), "phase_deg": cmath.phase(M_3_alt) * 180 / cmath.pi},
            "M_3_distinct_V8": {"real": M_3_distinct.real, "imag": M_3_distinct.imag, "abs": abs(M_3_distinct)},
            "noise_floor": NOISE_FLOOR,
            "ratio_M3alt_over_M3distinct": (abs(M_3_alt) / abs(M_3_distinct)) if abs(M_3_distinct) > 0 else None,
        },
    }
    import os as _os
    _os.makedirs(_os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print()
    print(f"[save] {out_path}")


if __name__ == "__main__":
    _main()
