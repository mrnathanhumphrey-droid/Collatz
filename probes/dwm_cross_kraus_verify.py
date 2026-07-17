"""
dwm_cross_kraus_verify.py — DWM-MP-G1 full closure.

Build Syracuse's Off_j operator EXPLICITLY as a sum of cross-Kraus M̃_{v,v'}
DWM building blocks, then compute ϕ(X̃_1·X̃_2·X̃_1) and verify it reproduces
Syracuse's measured value 0.108 (sum_entries scalar, Reading B marginal centering).

Cross-Kraus operator at step j conditional on pair-sum b and prior accumulator b_prior:

  M̃_{v,v'}^{(j, b_prior)}(b) · f(ξ) = √(cond_w_b(v,v')) · phase_cross_{v,v'}(ξ; b_prior, j)
                                       · f(ξ · 2^{-(v+v')} mod 27)

where cond_w_b(v,v') = 2^{-v-v'} / Σ_{(v'',v'''): v''+v'''=b, v''≠v'''} 2^{-v''-v'''},
and phase_cross_{v,v'}(ξ) = exp(-2πi · ξ · 3^{2j-2} · 2^{-b_prior} · (2^{-v} - 2^{-v'}) / 27).

Off_j(b, b_prior) = Σ_{(v,v'): v+v'=b, v≠v'} M̃_{v,v'}^{(j, b_prior)}(b)

After averaging over b distribution Pascal(2, 1/2) and b_prior:
  Off_j = E_b[Off_j(b, b_prior)] for fixed b_prior; X̃_j = Off_j - E_B[Off_j]
"""
import sys
import os
import cmath
import numpy as np
import json

sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, r'C:\Collatz')

from bilinear_pair_operator import build_markov_rational, stationary_rational

N_LEVEL = 3
N = 3**N_LEVEL
TWO_PI_I_OVER_N = 2j * cmath.pi / N
V_MAX = 16

print(f"# DWM cross-Kraus verification — DWM-MP-G1 full closure attempt")
print(f"# level n={N_LEVEL}, (Z/{N})*, V_MAX={V_MAX}")
print()

K_kernel, coprime = build_markov_rational(N_LEVEL)
pi_q = stationary_rational(K_kernel)
state_count = len(coprime)
state_idx = {r: i for i, r in enumerate(coprime)}
pi_f = np.array([float(p) for p in pi_q], dtype=float)
inv2_mod27 = pow(2, -1, N)
pow_inv2 = [pow(inv2_mod27, v, N) for v in range(0, 6 * V_MAX + 2)]


def shift_state(xi, exp_total):
    new_xi = (xi * pow_inv2[exp_total]) % N
    return state_idx.get(new_xi, -1)


# Pascal(2, 1/2) on b = v + v'
B_VALUES = list(range(2, 2 * V_MAX + 1))
pascal_unnorm = []
for b in B_VALUES:
    pairs_count_raw = sum(2.0**(-v - (b - v)) for v in range(1, V_MAX + 1)
                          if 1 <= b - v <= V_MAX)
    pascal_unnorm.append(pairs_count_raw)
Z_pascal = sum(pascal_unnorm)
pascal_w = [p / Z_pascal for p in pascal_unnorm]


def build_off_conditional_cross_kraus(j, b, b_prior):
    """Off_j(b, b_prior) = Σ_{v≠v', v+v'=b} cond_w · phase · shift via cross-Kraus."""
    # Enumerate (v, v') with v + v' = b, v != v', v, v' in [1, V_MAX]
    pairs = [(v, b - v) for v in range(1, V_MAX + 1)
             if 1 <= b - v <= V_MAX and v != b - v]
    if not pairs:
        return np.zeros((state_count, state_count), dtype=complex)

    # Conditional weights: 2^{-v-v'} normalized over cross-pairs summing to b
    raw_w = [2.0**(-v - vp) for v, vp in pairs]
    Zc = sum(raw_w)
    if Zc == 0:
        return np.zeros((state_count, state_count), dtype=complex)
    cw = [w / Zc for w in raw_w]

    M = np.zeros((state_count, state_count), dtype=complex)

    pow3_jm = pow(3, 2*j - 2, N)
    pow_inv2_bprior = pow(inv2_mod27, b_prior, N)
    x_j_mod = (pow3_jm * pow_inv2_bprior) % N

    for (v, vp), w in zip(pairs, cw):
        phase_diff_mod = (pow_inv2[v] - pow_inv2[vp]) % N
        exp_total = v + vp
        for i, xi in enumerate(coprime):
            target = shift_state(xi, exp_total)
            if target < 0:
                continue
            phase_arg = -TWO_PI_I_OVER_N * xi * x_j_mod * phase_diff_mod
            ph = cmath.exp(phase_arg)
            M[i, target] += w * ph
    return M


def center_marginal_B(Off_op):
    """Reading B marginal centering: X̃ = Off - projection onto stationary.
    Implementation matches verify_monotone_diagnostic.py."""
    proj = np.outer(np.ones(state_count), pi_f) @ Off_op
    return Off_op - proj


# SHARED-RANDOMNESS structure: X̃_1(b_1) appears TWICE with the SAME b_1
# This is the correct interpretation: a single random operator X̃_j(b_j) per step j,
# used at every occurrence in the monomial (the DWM "single fixed Kraus per step"
# structure that distinguishes Syracuse from iid-copies frameworks).
print("# Compute ϕ(X̃_1(b_1) · X̃_2(b_1, b_2) · X̃_1(b_1)) — SHARED b_1 in both X̃_1 positions")
print()
M_3_alt = 0.0
M_3_alt_trpi = 0.0
for b1, pw1 in zip(B_VALUES, pascal_w):
    if pw1 == 0:
        continue
    # X̃_1 at b_1 specifically — SAME operator used twice in the triple
    Off_1_b1 = build_off_conditional_cross_kraus(1, b1, 0)
    X_tilde_1_b1 = center_marginal_B(Off_1_b1)

    # X̃_2 conditional on b_prior_2 = b_1, averaged over b_2
    for b2, pw2 in zip(B_VALUES, pascal_w):
        if pw2 == 0:
            continue
        Off_2_b = build_off_conditional_cross_kraus(2, b2, b1)
        X_tilde_2_b1_b2 = center_marginal_B(Off_2_b)

        triple = X_tilde_1_b1 @ X_tilde_2_b1_b2 @ X_tilde_1_b1
        sum_ent = np.ones(state_count) @ triple @ np.ones(state_count)
        trpi = np.einsum('i,ii->', pi_f, triple)
        M_3_alt += pw1 * pw2 * sum_ent.real
        M_3_alt_trpi += pw1 * pw2 * trpi.real

print(f"  M_3_alt (sum_entries) cross-Kraus = {M_3_alt:.6e}")
print(f"  Syracuse direct measurement:        +1.078e-1")
print(f"  ratio                                = {M_3_alt / 0.1078:.6f}")
print()
print(f"  M_3_alt (tr_π) cross-Kraus  = {M_3_alt_trpi:.6e}")
print(f"  Syracuse direct (tr_π):       +5.36e-2 (approx, from Task 1)")
print()

# Save
out = {
    'description': 'DWM-MP-G1 cross-Kraus verification of 0.108',
    'level': N_LEVEL,
    'V_MAX': V_MAX,
    'M_3_alt_cross_kraus_sum_entries': M_3_alt,
    'M_3_alt_cross_kraus_tr_pi': M_3_alt_trpi,
    'syracuse_direct_sum_entries': 0.1078,
    'ratio_to_syracuse': M_3_alt / 0.1078,
}
out_path = r'C:\Collatz\experiments_output\dwm_cross_kraus_verify.json'
with open(out_path, 'w') as f:
    json.dump(out, f, indent=2, default=float)
print(f"[save] {out_path}")
