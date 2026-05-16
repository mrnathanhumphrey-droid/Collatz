"""
dwm_kraus_match_syracuse.py — DWM-MP-G1 exact-match closure.

Replicate Syracuse's verify_monotone_diagnostic.compute_M3_alt_readingB structure
exactly, but PHRASE the inner operator construction in DWM cross-Kraus building-
block language. Result should reproduce Syracuse's 0.1078 to floating-point
precision (since it IS the same operator-algebraic moment, just expressed in
DWM language).

The DWM cross-Kraus building block at SPECIFIC (v, v'):

  M̃_{v,v'}^{(j, b_prior)} = √(2^{-v-v'}) · phase_cross_{v,v'}(j, b_prior) · σ_{-(v+v')}

This is the "unweighted operator" times the √(weight) — i.e., the Stinespring
amplitude. The full operator is then weighted by 2^{-v-v'} (the Geom(2)² weight).

Off_j(v, v', b_prior) = "operator at this realization" = M̃_{v,v'} (single-realization)

E[Off_j | b_prior] = Σ_{(v,v'): v≠v'} 2^{-v-v'} · M̃_{v,v'}/√weight
                   = Σ 2^{-v-v'} · (phase · shift)
                   (= Off_unconditional_mean in Syracuse's notation)

X̃_j(v, v', b_prior) = M̃_{v,v'}/√weight · 2^{-v-v'}/2^{-v-v'} - E[Off_j | b_prior]
                     = (Off_op_unweighted at realization) - (Off_unconditional_mean)

Moment ϕ(X̃_1 · X̃_2 · X̃_1) = E over (v_1, v_1', v_2, v_2') of weights ·
                              scalar(X̃_1(v_1,v_1') · X̃_2(v_2,v_2',b_prior=v_1+v_1') · X̃_1(v_1,v_1'))

This is exactly Syracuse's structure, just re-expressed via DWM building blocks.
"""
import sys
import os
import cmath
from fractions import Fraction
import numpy as np
import json

sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, r'C:\Collatz')

from bilinear_pair_operator import build_markov_rational, stationary_rational

N_LEVEL = 3
N = 3**N_LEVEL
TWO_PI_I_OVER_N = 2j * cmath.pi / N
V_MAX = 16

print(f"# DWM cross-Kraus EXACT-MATCH attempt to Syracuse 0.1078")
print(f"# level n={N_LEVEL}, (Z/{N})*, V_MAX={V_MAX}")
print()

K_kernel, coprime = build_markov_rational(N_LEVEL)
pi_q = stationary_rational(K_kernel)
state_count = len(coprime)
state_idx = {r: i for i, r in enumerate(coprime)}
pi_f = np.array([float(p) for p in pi_q], dtype=float)
inv2_mod27 = pow(2, -1, N)
pow_inv2 = [pow(inv2_mod27, v, N) for v in range(0, 4 * V_MAX + 2)]


def shift_state(xi, exp_total):
    new_xi = (xi * pow_inv2[exp_total]) % N
    return state_idx.get(new_xi, -1)


def build_M_tilde_unweighted(v, vp, j, b_prior):
    """DWM cross-Kraus building block WITHOUT the 2^{-v-v'} weight prefactor.

    This matches Syracuse's build_Off_op_unweighted exactly: the operator at
    specific (v, v') realization, weight 1 (raw operator); the 2^{-v-v'} is
    applied OUTSIDE when integrating over realizations.

    Returns 18×18 complex matrix.
    """
    M = np.zeros((state_count, state_count), dtype=complex)
    if v == vp:
        return M
    pow3_jm = pow(3, 2*j - 2, N)
    pow_inv2_bprior = pow(inv2_mod27, b_prior, N)
    x_j_mod = (pow3_jm * pow_inv2_bprior) % N
    phase_diff_mod = (pow_inv2[v] - pow_inv2[vp]) % N
    exp_total = v + vp
    for i, xi in enumerate(coprime):
        target = shift_state(xi, exp_total)
        if target < 0:
            continue
        phase_arg = -TWO_PI_I_OVER_N * xi * x_j_mod * phase_diff_mod
        ph = cmath.exp(phase_arg)
        M[i, target] += ph
    return M


def Off_unconditional_mean(j, b_prior):
    """E_{(v, v') ~ Geom(2)^2 truncated, v≠v'}[Off_j_unweighted(v, v', b_prior)] · 2^{-v-v'}.

    Matches Syracuse's Off_unconditional_mean: sums 2^{-v-v'} · op_unweighted
    over (v, v') with v ≠ v', v, v' ∈ [1, V_MAX].
    """
    M = np.zeros((state_count, state_count), dtype=complex)
    for v in range(1, V_MAX + 1):
        for vp in range(1, V_MAX + 1):
            if v == vp:
                continue
            w_raw = 2.0**(-v - vp)
            M += w_raw * build_M_tilde_unweighted(v, vp, j, b_prior)
    return M


def sum_entries_scalar(*Ms):
    if not Ms:
        return 0 + 0j
    P = Ms[0].copy()
    for M in Ms[1:]:
        P = P @ M
    return complex(P.sum())


def tr_pi_scalar(*Ms):
    if not Ms:
        return 0 + 0j
    P = Ms[0].copy()
    for M in Ms[1:]:
        P = P @ M
    return complex(np.einsum('i,ii->', pi_f, P))


# Compute M_3_alt with Syracuse's exact structure ----------------------
print("# Computing M_3_alt under Syracuse's exact structure (raw realizations + Geom(2)² weights)")
print()

Off_1_mean_at_0 = Off_unconditional_mean(1, 0)

M_3_alt_sum_entries = 0 + 0j
M_3_alt_tr_pi = 0 + 0j

import time
t0 = time.time()
n_iter = 0
for v1 in range(1, V_MAX + 1):
    for vp1 in range(1, V_MAX + 1):
        if v1 == vp1:
            continue
        w1 = 2.0**(-v1 - vp1)
        b1 = v1 + vp1
        X1 = build_M_tilde_unweighted(v1, vp1, 1, 0) - Off_1_mean_at_0

        Off_2_mean_at_b1 = Off_unconditional_mean(2, b1)
        for v2 in range(1, V_MAX + 1):
            for vp2 in range(1, V_MAX + 1):
                if v2 == vp2:
                    continue
                w2 = 2.0**(-v2 - vp2)
                X2 = build_M_tilde_unweighted(v2, vp2, 2, b1) - Off_2_mean_at_b1
                contribution = w1 * w2 * sum_entries_scalar(X1, X2, X1)
                M_3_alt_sum_entries += contribution
                contribution_trpi = w1 * w2 * tr_pi_scalar(X1, X2, X1)
                M_3_alt_tr_pi += contribution_trpi
                n_iter += 1

t1 = time.time()
print(f"  iterations: {n_iter}, elapsed {t1-t0:.1f}s")
print()
print(f"  M_3_alt (sum_entries scalar) = {M_3_alt_sum_entries.real:+.6e} + {M_3_alt_sum_entries.imag:+.2e}i")
print(f"  Syracuse Task 1 reported:      +1.0783×10⁻¹")
print(f"  ratio = {M_3_alt_sum_entries.real / 0.10783:.6f}")
print()
print(f"  M_3_alt (tr_π scalar) = {M_3_alt_tr_pi.real:+.6e} + {M_3_alt_tr_pi.imag:+.2e}i")
print(f"  Syracuse Task 1 reported (tr_π): ~+5.36×10⁻²")
print()

# Save
out = {
    'description': 'DWM cross-Kraus moment matched to Syracuse exact structure',
    'level': N_LEVEL,
    'V_MAX': V_MAX,
    'n_iter': n_iter,
    'elapsed_sec': t1 - t0,
    'M_3_alt_sum_entries_real': M_3_alt_sum_entries.real,
    'M_3_alt_sum_entries_imag': M_3_alt_sum_entries.imag,
    'M_3_alt_tr_pi_real': M_3_alt_tr_pi.real,
    'M_3_alt_tr_pi_imag': M_3_alt_tr_pi.imag,
    'syracuse_target_sum_entries': 0.10783,
    'syracuse_target_tr_pi': 0.0536,
    'ratio_sum_entries': M_3_alt_sum_entries.real / 0.10783,
}
out_path = r'C:\Collatz\experiments_output\dwm_kraus_match_syracuse.json'
with open(out_path, 'w') as f:
    json.dump(out, f, indent=2, default=float)
print(f"[save] {out_path}")
