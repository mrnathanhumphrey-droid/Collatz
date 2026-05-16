"""
dwm_kraus_match_g2.py — DWM-MP-G2 closure via DWM cross-Kraus form.

Compute ϕ(X̃_1·X̃_2·X̃_1·X̃_2) at n=3 using Syracuse's exact realization-weight
structure (raw (v, v') iteration + Geom(2)² weights), reusing the building
blocks from dwm_kraus_match_syracuse.py.

The (j_1, j_2, j_1, j_2) pattern has SHARED (v_1, v_1') across the two X̃_1
positions and SHARED (v_2, v_2') across the two X̃_2 positions. So we iterate
over V_MAX² × V_MAX² = V_MAX^4 (same as M_3_alt cost, ~65k iterations).
"""
import sys
import cmath
import numpy as np
import json
import time

sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, r'C:\Collatz')

from bilinear_pair_operator import build_markov_rational, stationary_rational

N_LEVEL = 3
N = 3**N_LEVEL
TWO_PI_I_OVER_N = 2j * cmath.pi / N
V_MAX = 16

print(f"# DWM-MP-G2: ϕ(X̃_1·X̃_2·X̃_1·X̃_2) via DWM cross-Kraus form")
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
    M = np.zeros((state_count, state_count), dtype=complex)
    for v in range(1, V_MAX + 1):
        for vp in range(1, V_MAX + 1):
            if v == vp:
                continue
            w_raw = 2.0**(-v - vp)
            M += w_raw * build_M_tilde_unweighted(v, vp, j, b_prior)
    return M


print("# Computing M_4_alt = ϕ(X̃_1·X̃_2·X̃_1·X̃_2) — shared (v_1,v_1'), shared (v_2,v_2')")
print()

Off_1_mean_at_0 = Off_unconditional_mean(1, 0)

M_4_alt_sum_entries = 0 + 0j
M_4_alt_tr_pi = 0 + 0j
M_4_alt_delta1 = 0 + 0j
M_4_alt_vacpi = 0 + 0j

idx_xi_1 = state_idx[1]

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

                # 4-product: X1 · X2 · X1 · X2
                P12 = X1 @ X2
                P1212 = P12 @ P12  # = X1 X2 X1 X2

                # scalar reductions
                M_4_alt_sum_entries += w1 * w2 * P1212.sum()
                M_4_alt_tr_pi += w1 * w2 * np.einsum('i,ii->', pi_f, P1212)
                M_4_alt_delta1 += w1 * w2 * P1212[idx_xi_1, idx_xi_1]
                M_4_alt_vacpi += w1 * w2 * (pi_f @ P1212 @ pi_f)
                n_iter += 1

t1 = time.time()
print(f"  iterations: {n_iter}, elapsed {t1-t0:.1f}s")
print()
print(f"  M_4_alt (sum_entries)   = {M_4_alt_sum_entries.real:+.6e} + {M_4_alt_sum_entries.imag:+.2e}i")
print(f"  Syracuse Task 1:           +6.089×10⁻¹")
print(f"  ratio = {M_4_alt_sum_entries.real / 0.6089:.6f}")
print()
print(f"  M_4_alt (tr_π)          = {M_4_alt_tr_pi.real:+.6e}")
print(f"  M_4_alt (delta_1)       = {M_4_alt_delta1.real:+.6e}")
print(f"    Syracuse delta_1:        +5.742×10⁻²")
print(f"  M_4_alt (vac_π)         = {M_4_alt_vacpi.real:+.6e}")
print(f"    Syracuse vac_π:          +4.775×10⁻³")
print()

# tr_pi reference comparison
print(f"  Syracuse tr_pi reported: +5.357×10⁻² (Task 1)")
print(f"  ratio_trpi = {M_4_alt_tr_pi.real / 0.05357:.6f}")
print()

out = {
    'description': 'DWM-MP-G2: M_4_alt via DWM cross-Kraus, exact-match attempt to Syracuse',
    'level': N_LEVEL,
    'V_MAX': V_MAX,
    'n_iter': n_iter,
    'elapsed_sec': t1 - t0,
    'M_4_alt_sum_entries_real': M_4_alt_sum_entries.real,
    'M_4_alt_tr_pi_real': M_4_alt_tr_pi.real,
    'M_4_alt_delta1_real': M_4_alt_delta1.real,
    'M_4_alt_vacpi_real': M_4_alt_vacpi.real,
    'syracuse_targets': {
        'sum_entries': 0.6089,
        'tr_pi': 0.05357,
        'delta_1': 0.05742,
        'vac_pi': 0.004775,
    },
}
import os
out_path = r'C:\Collatz\experiments_output\dwm_kraus_match_g2.json'
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w') as f:
    json.dump(out, f, indent=2, default=float)
print(f"[save] {out_path}")
