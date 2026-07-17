"""
verify_n4_alternating.py — D1 numerical confirmation.

Compute the n = 4 alternating moment

    M_4_alt := E_B( X_tilde_{j_1} . X_tilde_{j_2} . X_tilde_{j_1} . X_tilde_{j_2} )

at level n = 3, state space (Z/27)*, with j_1 = 1, j_2 = 2, under Reading (B)
centering (X_tilde_j = Off_j - E[Off_j | b_prior]) — the same reading that
gave |M_3_alt| = 0.10783 in Task 1.

H1' (HS 2014 Defn 2.2 verbatim, single-phi reading) predicts the peak-rule
substitution at position 2 OR position 4 forces RHS = 0; whereas the structural
argument in H1_PRIME_LOW_ORDER_CHECKS.md sec 4 predicts LHS != 0 generically.
D1 measures the gap.

CONTROLS:
  - M_4_distinct: four distinct indices j_1=1, j_2=2, j_3=3, j_4=4 (b_prior_4 = b_1+b_2+b_3)
    should be ~ 0 by full distinct-index independence.
  - M_2 distinct: already 1.076e-7 (sum_entries) from Task 1 — recompute here.
  - Reading-A control on M_4_alt (algebraic zero by linearity over splits).

FUBINI DECOMPOSITION (sec 4 of H1_PRIME_LOW_ORDER_CHECKS):
  Compute the conditional inner factor at fixed step-1 split (v_1, v_1'):
      F(v_1, v_1') := E_{(v_2, v_2')} [
                        X_tilde_2(v_2, v_2'; b_prior=b_1) .
                        X_tilde_1(v_1, v_1'; b_prior=0) .
                        X_tilde_2(v_2, v_2'; b_prior=b_1) ]
  This is a B-valued operator at step 1; the moment then collapses to
      M_4_alt = E_{(v_1, v_1')} [ X_tilde_1 . F(v_1, v_1') ]_scalar
  We verify the Fubini inner factor is non-zero (i.e. there's something
  positive to integrate against the outer X_tilde_1).

OUTPUT: experiments_output/n4_alternating_diagnostic.json
"""
import sys
import os
import json
import time
import cmath
from fractions import Fraction
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

sys.path.insert(0, r"C:\Collatz")
from bilinear_pair_operator import build_markov_rational, stationary_rational

# Level / state space (match Task 1)
N_LEVEL = 3
N = 3**N_LEVEL
TWO_PI_I_OVER_N = 2j * cmath.pi / N
V_MAX = 16

# Pre-computed primitives
geom_weights_unnorm = [Fraction(1, 2**v) for v in range(1, V_MAX + 1)]
Z_geom = sum(geom_weights_unnorm)
geom_weights = [w / Z_geom for w in geom_weights_unnorm]
inv2_mod27 = pow(2, -1, N)
pow_inv2 = [pow(inv2_mod27, v, N) for v in range(0, 8 * V_MAX + 2)]

print("# verify_n4_alternating.py (D1)")
print(f"# level n = {N_LEVEL}; state space (Z/{N})*")
print()
print("# building stationary distribution pi_3 ...")
K, coprime = build_markov_rational(N_LEVEL)
pi_q = stationary_rational(K)
state_count = len(coprime)
print(f"#   state count = {state_count}")
state_idx = {r: i for i, r in enumerate(coprime)}
pi_f = np.array([float(p) for p in pi_q], dtype=float)
print(f"#   pi sum check: {pi_f.sum():.12f}")
print()


def shift_state(xi, exp_total):
    new_xi = (xi * pow_inv2[exp_total]) % N
    return state_idx.get(new_xi, -1)


def build_Off_op_unweighted(j, v, vp, b_prior):
    """Operator at specific (v, vp), step j, prior accumulator b_prior — WITHOUT
    the 2^{-v-vp} weight. The weight enters multiplicatively when summing
    realizations.
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


def Off_unconditional_mean(j, b_prior):
    """E_{(v,v')}[Off_j realized] — full unconditional pair-marginal mean."""
    M = np.zeros((state_count, state_count), dtype=complex)
    for v in range(1, V_MAX + 1):
        for vp in range(1, V_MAX + 1):
            if v == vp:
                continue
            w_raw = float(Fraction(1, 2**v * 2**vp))
            M += w_raw * build_Off_op_unweighted(j, v, vp, b_prior)
    return M


def conditional_splits(b):
    """Reading-A conditional distribution of (v,v') given b_j = b."""
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
    """Reading-A conditional mean operator at fixed b_j = b, b_prior."""
    splits = conditional_splits(b)
    M = np.zeros((state_count, state_count), dtype=complex)
    for v, vp, p in splits:
        M += p * build_Off_op_unweighted(j, v, vp, b_prior)
    return M


# Pascal(2, 1/2) distribution for b_j (only needed for Reading-A control)
B_VALUES = list(range(2, 2 * V_MAX + 1))
pascal_unnorm = {}
for b in B_VALUES:
    n_pairs = sum(1 for v in range(1, V_MAX + 1) if 1 <= b - v <= V_MAX)
    if n_pairs == 0:
        pascal_unnorm[b] = Fraction(0)
        continue
    total = sum(Fraction(1, 2**v * 2**(b - v)) for v in range(1, V_MAX + 1)
                if 1 <= b - v <= V_MAX)
    pascal_unnorm[b] = total
Z_b = sum(pascal_unnorm.values())
pascal = {b: pascal_unnorm[b] / Z_b for b in B_VALUES}
pascal_f = {b: float(p) for b, p in pascal.items()}


# --- scalar reductions (match Task 1) -------------------------------------
pi_vec = pi_f.copy()
idx_xi_1 = state_idx[1]


def tr_pi(M):
    return complex(sum(pi_f[i] * M[i, i] for i in range(state_count)))


def expect_delta1(M):
    return complex(M[idx_xi_1, idx_xi_1])


def sum_all_entries(M):
    return complex(M.sum())


def _prod(Ms):
    P = Ms[0].copy()
    for M in Ms[1:]:
        P = P @ M
    return P


def tr_pi_product(*Ms):
    return tr_pi(_prod(Ms))


def vac_product_pi(*Ms):
    P = _prod(Ms)
    return complex(np.conj(pi_vec) @ P @ pi_vec)


def delta1_product(*Ms):
    return expect_delta1(_prod(Ms))


def sum_entries_product(*Ms):
    return sum_all_entries(_prod(Ms))


SCALARS = [
    ("tr_pi", tr_pi_product),
    ("vac_pi", vac_product_pi),
    ("delta_1", delta1_product),
    ("sum_entries", sum_entries_product),
]


# ==========================================================================
# Reading-B moment computations at n = 4
# ==========================================================================

def compute_M4_alt_readingB(scalar_fn, V_TRUNC=V_MAX):
    """M_4_alt := E[ X1 . X2 . X1 . X2 ] with same (v_1,v_1') at positions 1, 3
    and same (v_2,v_2') at positions 2, 4 (Reading B).

    Outer sum: over (v_1, v_1', v_2, v_2') with weight 2^{-v_1-v_1'-v_2-v_2'}.

    Cost: O(V_TRUNC^4 * 18^3) for the product. V_TRUNC=16 -> 65536 quadruples,
    each ~3 mat-muls of 18x18 -> ~ minutes. Doable.
    """
    Off_1_mean_0 = Off_unconditional_mean(1, 0)
    total = 0 + 0j
    for v1 in range(1, V_TRUNC + 1):
        for vp1 in range(1, V_TRUNC + 1):
            if v1 == vp1:
                continue
            w1 = float(Fraction(1, 2**v1 * 2**vp1))
            b1 = v1 + vp1
            X1 = build_Off_op_unweighted(1, v1, vp1, 0) - Off_1_mean_0
            Off_2_mean = Off_unconditional_mean(2, b1)
            for v2 in range(1, V_TRUNC + 1):
                for vp2 in range(1, V_TRUNC + 1):
                    if v2 == vp2:
                        continue
                    w2 = float(Fraction(1, 2**v2 * 2**vp2))
                    X2 = build_Off_op_unweighted(2, v2, vp2, b1) - Off_2_mean
                    total += w1 * w2 * scalar_fn(X1, X2, X1, X2)
    return total


def compute_M2_readingB(scalar_fn, V_TRUNC=V_MAX):
    """M_2 = E[ X1 . X2 ]. Reading B. (Recompute for cross-check.)"""
    Off_1_mean_0 = Off_unconditional_mean(1, 0)
    total = 0 + 0j
    for v1 in range(1, V_TRUNC + 1):
        for vp1 in range(1, V_TRUNC + 1):
            if v1 == vp1:
                continue
            w1 = float(Fraction(1, 2**v1 * 2**vp1))
            b1 = v1 + vp1
            X1 = build_Off_op_unweighted(1, v1, vp1, 0) - Off_1_mean_0
            Off_2_mean = Off_unconditional_mean(2, b1)
            for v2 in range(1, V_TRUNC + 1):
                for vp2 in range(1, V_TRUNC + 1):
                    if v2 == vp2:
                        continue
                    w2 = float(Fraction(1, 2**v2 * 2**vp2))
                    X2 = build_Off_op_unweighted(2, v2, vp2, b1) - Off_2_mean
                    total += w1 * w2 * scalar_fn(X1, X2)
    return total


def compute_M4_distinct_readingB(scalar_fn, V_TRUNC=6):
    """M_4_distinct := E[ X1 . X2 . X3 . X4 ] with j_1=1, j_2=2, j_3=3, j_4=4
    all distinct. Reading B. Should be ~ 0 by full distinct-index independence.

    Cost: V_TRUNC^8. V_TRUNC=6 -> 6^8 ~ 1.68M, each with 4 mat-muls and 4 Off
    means. Cap at V_TRUNC=6 for tractability. Per Geom(2): tail at v=6 is
    2^{-6} = 1.5e-2; total truncated mass for one (v,v') is ~ (1-2^{-6})^2 ~ 0.97.
    """
    Off_1_mean_0 = Off_unconditional_mean(1, 0)
    total = 0 + 0j
    n_quad = 0
    for v1 in range(1, V_TRUNC + 1):
        for vp1 in range(1, V_TRUNC + 1):
            if v1 == vp1:
                continue
            w1 = float(Fraction(1, 2**v1 * 2**vp1))
            b1 = v1 + vp1
            X1 = build_Off_op_unweighted(1, v1, vp1, 0) - Off_1_mean_0
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
                            b3 = v3 + vp3
                            b_prior_4 = b1 + b2 + b3
                            X3 = build_Off_op_unweighted(3, v3, vp3, b_prior_3) - Off_3_mean
                            Off_4_mean = Off_unconditional_mean(4, b_prior_4)
                            for v4 in range(1, V_TRUNC + 1):
                                for vp4 in range(1, V_TRUNC + 1):
                                    if v4 == vp4:
                                        continue
                                    w4 = float(Fraction(1, 2**v4 * 2**vp4))
                                    X4 = build_Off_op_unweighted(4, v4, vp4, b_prior_4) - Off_4_mean
                                    total += w1 * w2 * w3 * w4 * scalar_fn(X1, X2, X3, X4)
                                    n_quad += 1
    return total, n_quad


# ==========================================================================
# Reading-A control on M_4_alt (algebraic zero expected)
# ==========================================================================

def compute_M4_alt_readingA(scalar_fn):
    """Reading-A centering: E[X1 . X2 . X1 . X2] where X_j = Off_j - Off_avg_op(j, b_j, b_prior).
    Within a fixed (b_1, b_2), splits are independent across steps; the inner
    E[X_j | b_j, b_prior] = 0 forces the second X_2 factor to integrate to zero,
    so the moment must vanish algebraically. Numerical check.
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
                    acc += p1 * p2 * scalar_fn(X1, X2, X1, X2)
            total += p_b1 * p_b2 * acc
    return total


# ==========================================================================
# Peak-rule prediction (under H1' assumed-true)
# ==========================================================================

def compute_peak_rule_prediction(scalar_fn, peak_position=2, V_TRUNC=V_MAX):
    """H1' under Defn 2.2: at peak position i, substitute X_{j_i} -> phi(X_{j_i}).
    Under Reading-B, phi(X_tilde_j) = E[X_tilde_j] (marginal expectation).

    By centering: E[X_tilde_j] = 0 (the unconditional mean is what we subtracted).
    Thus the peak substitution at position 2 (or 4) gives X1 . 0 . X1 . X2 or
    X1 . X2 . X1 . 0, both with scalar 0. RHS = 0 in either case.

    This function returns 0 by construction — it's the H1'-predicted RHS, just
    for documentation symmetry.
    """
    return 0 + 0j


# ==========================================================================
# Fubini inner factor (sec 4 of H1_PRIME_LOW_ORDER_CHECKS)
# ==========================================================================

def compute_fubini_inner_factor(v1, vp1, scalar_fn, V_TRUNC=V_MAX):
    """At fixed step-1 realization (v_1, v_1') with b_prior_1 = 0, compute
    the inner B-valued operator (scalar reduction):

        F_scalar(v_1, v_1') = E_{(v_2, v_2')} [
                                scalar_fn( X_tilde_2 . X_tilde_1 . X_tilde_2 ) ]

    where X_tilde_1 is fixed at this (v_1, v_1') realization and X_tilde_2 is
    averaged over its pair (with b_prior_2 = b_1 = v_1 + vp_1).

    The Fubini argument in H1_PRIME_LOW_ORDER_CHECKS sec 4 says: if this inner
    factor is generically positive across (v_1, v_1'), then integrating against
    X_tilde_1 (which has marginal mean zero but variance > 0) gives a non-zero
    outer integral — confirming the structural mechanism for the M_4_alt
    non-vanishing.
    """
    Off_1_mean_0 = Off_unconditional_mean(1, 0)
    b1 = v1 + vp1
    X1 = build_Off_op_unweighted(1, v1, vp1, 0) - Off_1_mean_0
    Off_2_mean = Off_unconditional_mean(2, b1)
    acc = 0 + 0j
    for v2 in range(1, V_TRUNC + 1):
        for vp2 in range(1, V_TRUNC + 1):
            if v2 == vp2:
                continue
            w2 = float(Fraction(1, 2**v2 * 2**vp2))
            X2 = build_Off_op_unweighted(2, v2, vp2, b1) - Off_2_mean
            acc += w2 * scalar_fn(X2, X1, X2)
    return acc


# ==========================================================================
# MAIN
# ==========================================================================

def _main():
    print("=" * 70)
    print("# READING (B) — primary computations at n = 4")
    print(f"# V_TRUNC for M_4_alt = {V_MAX} (full)")
    print("# V_TRUNC for M_4_distinct = 6 (V^8 cost)")
    print("=" * 70)

    readingB = {}
    for name, fn in SCALARS:
        t0 = time.time()
        m2 = compute_M2_readingB(fn, V_TRUNC=V_MAX)
        t_m2 = time.time() - t0
        t0 = time.time()
        m4a = compute_M4_alt_readingB(fn, V_TRUNC=V_MAX)
        t_m4a = time.time() - t0
        readingB[name] = {"M_2": m2, "M_4_alt": m4a}
        print(f"  [{name}]")
        print(f"    M_2     = {m2.real:+.6e}{m2.imag:+.6e}j  |.| = {abs(m2):.6e}  ({t_m2:.1f}s)")
        print(f"    M_4_alt = {m4a.real:+.6e}{m4a.imag:+.6e}j  |.| = {abs(m4a):.6e}  ({t_m4a:.1f}s)")

    print()
    print("=" * 70)
    print("# READING (B) — M_4_distinct control (4 distinct indices, V_TRUNC=6)")
    print("=" * 70)
    distinct_results = {}
    for name, fn in SCALARS:
        t0 = time.time()
        m4d, n_quad = compute_M4_distinct_readingB(fn, V_TRUNC=6)
        t_m4d = time.time() - t0
        distinct_results[name] = m4d
        print(f"  [{name}] M_4_distinct = {m4d.real:+.6e}{m4d.imag:+.6e}j  |.| = {abs(m4d):.6e}  "
              f"(n_quad={n_quad}, {t_m4d:.1f}s)")

    print()
    print("=" * 70)
    print("# READING (A) — control on M_4_alt (algebraic zero expected)")
    print("=" * 70)
    readingA_control = {}
    for name, fn in SCALARS:
        t0 = time.time()
        m4a_A = compute_M4_alt_readingA(fn)
        t_a = time.time() - t0
        readingA_control[name] = m4a_A
        print(f"  [{name}] M_4_alt (Reading A) = {m4a_A.real:+.6e}{m4a_A.imag:+.6e}j  "
              f"|.| = {abs(m4a_A):.6e}  ({t_a:.1f}s)")

    print()
    print("=" * 70)
    print("# Fubini inner factor F(v_1, v_1') across a grid (Reading B, sum_entries)")
    print("# F = sum_{(v_2,v_2')} 2^{-v_2-v_2'} . scalar( X_tilde_2 . X_tilde_1 . X_tilde_2 )")
    print("=" * 70)
    fubini_results = {}
    fubini_grid_pairs = [(1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2),
                         (1, 4), (4, 1), (2, 4), (4, 2), (1, 5), (5, 1)]
    for scalar_name, scalar_fn in SCALARS:
        per_scalar = []
        for (v1, vp1) in fubini_grid_pairs:
            F_val = compute_fubini_inner_factor(v1, vp1, scalar_fn, V_TRUNC=V_MAX)
            per_scalar.append({"v1": v1, "vp1": vp1, "F_real": F_val.real,
                               "F_imag": F_val.imag, "F_abs": abs(F_val)})
        fubini_results[scalar_name] = per_scalar
        # Print compact for sum_entries
        if scalar_name == "sum_entries":
            print(f"  [{scalar_name}]")
            for entry in per_scalar:
                print(f"    (v1={entry['v1']:2d}, vp1={entry['vp1']:2d}) "
                      f"F = {entry['F_real']:+.4e}{entry['F_imag']:+.4e}j  "
                      f"|F| = {entry['F_abs']:.4e}")

    print()
    print("=" * 70)
    print("# DIAGNOSTIC SUMMARY")
    print("=" * 70)
    print()
    print("  Reading-B M_4_alt across scalar reductions:")
    for name, _ in SCALARS:
        m4a = readingB[name]["M_4_alt"]
        print(f"    {name:12s}: |M_4_alt| = {abs(m4a):.6e}")
    print()
    print("  Should-vanish controls:")
    print(f"    Reading-B M_2 (sum_entries)         = {abs(readingB['sum_entries']['M_2']):.6e}")
    print(f"    Reading-B M_4_distinct (sum_entries, V=6) = {abs(distinct_results['sum_entries']):.6e}")
    print(f"    Reading-A M_4_alt control (sum_entries)   = {abs(readingA_control['sum_entries']):.6e}")
    print()
    print("  H1' violation strength = |LHS - RHS_peak_rule| = |M_4_alt - 0| = |M_4_alt|")
    for name, _ in SCALARS:
        m4a = readingB[name]["M_4_alt"]
        m4d = distinct_results[name]
        ratio = abs(m4a) / abs(m4d) if abs(m4d) > 0 else float('inf')
        print(f"    {name:12s}: H1' gap = {abs(m4a):.6e},  ratio M4alt/M4dist = {ratio:.2e}")
    print()

    # Verdict
    primary = readingB["sum_entries"]["M_4_alt"]
    distinct_primary = distinct_results["sum_entries"]
    reading_A_floor = abs(readingA_control["sum_entries"])
    noise_floor = max(1e-12, 10 * reading_A_floor)
    print("  === VERDICT ===")
    if abs(primary) > 100 * noise_floor and abs(primary) > 10 * abs(distinct_primary):
        verdict = "H1_PRIME_VIOLATION_CONFIRMED_NUMERICALLY"
        print(f"  {verdict}")
        print(f"  |M_4_alt|             = {abs(primary):.6e}")
        print(f"  |M_4_distinct|        = {abs(distinct_primary):.6e}")
        print(f"  Reading-A noise floor = {reading_A_floor:.6e}")
        print(f"  Gap to RHS_peak_rule  = {abs(primary):.6e} (= |LHS|, since RHS=0)")
    elif abs(primary) > 100 * noise_floor:
        verdict = "BORDERLINE — M4_alt non-zero but comparable to distinct control"
        print(f"  {verdict}")
    else:
        verdict = "NOT_CONFIRMED — M4_alt at noise floor"
        print(f"  {verdict}")

    # save
    out_path = r"C:\Collatz\experiments_output\n4_alternating_diagnostic.json"
    results = {
        "task_id": "D1",
        "level": N_LEVEL,
        "state_count": state_count,
        "V_MAX": V_MAX,
        "V_TRUNC_M4_alt": V_MAX,
        "V_TRUNC_M4_distinct": 6,
        "j_1": 1,
        "j_2": 2,
        "operator_realization_note": (
            "Off_j(f)(xi) = sum_{v!=v'} 2^{-v} 2^{-v'} exp(-2 pi i xi * 3^{2j-2} "
            "* 2^{-b_prior} * (2^{-v} - 2^{-v'}) / 27) * f(xi * 2^{-(v+v')} mod 27). "
            "Same realization as verify_monotone_diagnostic.py."
        ),
        "centering_reading": (
            "Reading B (per AMALG_FREENESS_MOMENT_CALCULATION sec 11): "
            "X_tilde_j = Off_j - E[Off_j | b_prior]. Marginal centering over (v, v')."
        ),
        "scalar_reductions": [s[0] for s in SCALARS],
        "reading_B": {
            name: {
                "M_2": {"real": readingB[name]["M_2"].real,
                        "imag": readingB[name]["M_2"].imag,
                        "abs": abs(readingB[name]["M_2"])},
                "M_4_alt": {"real": readingB[name]["M_4_alt"].real,
                            "imag": readingB[name]["M_4_alt"].imag,
                            "abs": abs(readingB[name]["M_4_alt"]),
                            "phase_rad": cmath.phase(readingB[name]["M_4_alt"]),
                            "phase_deg": cmath.phase(readingB[name]["M_4_alt"]) * 180 / cmath.pi},
            }
            for name, _ in SCALARS
        },
        "reading_B_M_4_distinct_control": {
            name: {
                "real": distinct_results[name].real,
                "imag": distinct_results[name].imag,
                "abs": abs(distinct_results[name]),
            }
            for name, _ in SCALARS
        },
        "reading_A_M_4_alt_control": {
            name: {
                "real": readingA_control[name].real,
                "imag": readingA_control[name].imag,
                "abs": abs(readingA_control[name]),
            }
            for name, _ in SCALARS
        },
        "peak_rule_RHS_prediction": {
            "value": 0.0,
            "rationale": (
                "Under H1' Defn 2.2 with phi(X_tilde_j) = E[X_tilde_j] = 0 "
                "(marginal centering), peak-rule substitution at position 2 "
                "OR position 4 gives RHS = phi(X_1 . 0 . X_1 . X_2) = 0 or "
                "phi(X_1 . X_2 . X_1 . 0) = 0. Either way RHS = 0."
            ),
        },
        "H1_prime_violation_strength": {
            name: {
                "LHS": readingB[name]["M_4_alt"].real,
                "RHS_peak_rule": 0.0,
                "abs_gap": abs(readingB[name]["M_4_alt"]),
            }
            for name, _ in SCALARS
        },
        "fubini_inner_factor": fubini_results,
        "verdict": verdict,
        "diagnostic_primary": {
            "scalar_form": "sum_entries",
            "reading": "B",
            "M_4_alt": {"real": primary.real, "imag": primary.imag,
                        "abs": abs(primary),
                        "phase_rad": cmath.phase(primary),
                        "phase_deg": cmath.phase(primary) * 180 / cmath.pi},
            "M_4_distinct_V6": {"real": distinct_primary.real,
                                "imag": distinct_primary.imag,
                                "abs": abs(distinct_primary)},
            "noise_floor": noise_floor,
            "ratio_M4alt_over_M4distinct": (abs(primary) / abs(distinct_primary))
                if abs(distinct_primary) > 0 else None,
        },
    }
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print()
    print(f"[save] {out_path}")


if __name__ == "__main__":
    _main()
