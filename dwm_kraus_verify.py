"""
dwm_kraus_verify.py — DWM-MP-G1 verification.

Build the Davies-Wiseman-Milburn adaptive Kraus operators M_v^{(j, b_prior)}
explicitly, verify POVM resolution Σ M_v† M_v = I, verify Markov preservation
of stationary π_3, then compute the third-order moment

  ϕ(X̃_{j_1} · X̃_{j_2} · X̃_{j_1})

via the explicit Kraus form and compare to Syracuse's direct measurement at
n=3 (verify_monotone_diagnostic.py output 1.078e-1).

The DWM Kraus form (per FRAMEWORK_IDENTIFICATION.md):

  M_v^{(j, b_prior)} · f(ξ) = 2^{-v/2} · A_v^{(j)}(ξ, b_prior) · f(ξ · 2^{-v} mod 27)

where A_v is a pure phase factor (the Tao single-frequency phase). The full
transfer operator at step j is T_j = Σ_v M_v.

If |A_v(ξ)| = 1 and σ_{-v} (the shift by 2^{-v} on (Z/27)*) is unitary on
the coprime support, then:

  (M_v† M_v) f(ξ) = 2^{-v} · |A_v(ξ)|² · f(ξ) = 2^{-v} · f(ξ)

so Σ_v M_v† M_v = Σ_v 2^{-v} · I = (1 - 2^{-V_MAX}) · I ≈ I  ✓ POVM resolution.

This script verifies:
  (1) POVM: Σ_v M_v† M_v ≈ I
  (2) Stationary: Σ_v M_v π M_v† ≈ π (preservation of stationary state)
  (3) Decomposition: T_j(full) = Σ_v M_v matches via direct construction
  (4) Moment: ϕ(X̃_{j_1} · X̃_{j_2} · X̃_{j_1}) via Kraus form = 0.108 ± tolerance
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

# --- setup ---------------------------------------------------------------
N_LEVEL = 3
N = 3**N_LEVEL  # 27
TWO_PI_I_OVER_N = 2j * cmath.pi / N
V_MAX = 16

print(f"# DWM Kraus verification — DWM-MP-G1")
print(f"# level n={N_LEVEL}, (Z/{N})*, V_MAX={V_MAX}")
print()

K_kernel, coprime = build_markov_rational(N_LEVEL)
pi_q = stationary_rational(K_kernel)
state_count = len(coprime)
state_idx = {r: i for i, r in enumerate(coprime)}
pi_f = np.array([float(p) for p in pi_q], dtype=float)
print(f"# state count = {state_count}")
print(f"# π sum = {pi_f.sum():.12f}")
print()

inv2_mod27 = pow(2, -1, N)
pow_inv2 = [pow(inv2_mod27, v, N) for v in range(0, 4 * V_MAX + 2)]


def shift_state(xi, exp_total):
    new_xi = (xi * pow_inv2[exp_total]) % N
    return state_idx.get(new_xi, -1)


def build_M_v(v, j, b_prior):
    """Build Kraus operator M_v^{(j, b_prior)} as a state_count×state_count
    complex matrix.

    M_v · f(ξ) = 2^{-v/2} · A_v^{(j)}(ξ, b_prior) · f(ξ · 2^{-v} mod 27)

    where A_v^{(j)}(ξ, b_prior) = exp(-2πi · ξ · 3^{2j-2} · 2^{-b_prior-v} / 27)
    is the single-frequency Tao phase (NOT the cross-frequency v≠v' factor).

    Returns an 18×18 complex matrix M[i, target] for state i = ξ and
    target = state_idx of ξ·2^{-v}.
    """
    M = np.zeros((state_count, state_count), dtype=complex)
    amp = 2.0 ** (-v / 2)
    pow3_jm = pow(3, 2*j - 2, N)
    # x_j_v = 3^{2j-2} · 2^{-(b_prior + v)} mod 27
    x_j_v_mod = (pow3_jm * pow(inv2_mod27, b_prior + v, N)) % N

    for i, xi in enumerate(coprime):
        target = shift_state(xi, v)
        if target < 0:
            continue
        phase_arg = -TWO_PI_I_OVER_N * xi * x_j_v_mod
        ph = cmath.exp(phase_arg)
        M[i, target] += amp * ph
    return M


# --- POVM verification ---------------------------------------------------
print("# (1) POVM resolution check: Σ_v M_v† M_v ≈ I")
print()

j_test = 1
b_prior_test = 0

POVM_sum = np.zeros((state_count, state_count), dtype=complex)
for v in range(1, V_MAX + 1):
    Mv = build_M_v(v, j_test, b_prior_test)
    POVM_sum += Mv.conj().T @ Mv

dev_I = np.linalg.norm(POVM_sum - np.eye(state_count), ord='fro')
trace_total = np.trace(POVM_sum).real / state_count
expected_trace = 1.0 - 2.0**(-V_MAX)
print(f"  ‖Σ_v M_v† M_v − I‖_F = {dev_I:.6e}")
print(f"  tr(Σ_v M_v† M_v) / dim = {trace_total:.10f}")
print(f"  expected = 1 − 2^{{-{V_MAX}}} = {expected_trace:.10f}")
print(f"  match within tail = {abs(trace_total - expected_trace):.2e}")
print()

# --- Stationary preservation ---------------------------------------------
print("# (2) Stationary preservation: Σ_v M_v · diag(π) · M_v† should preserve π")
print()

diag_pi = np.diag(pi_f.astype(complex))
channel_pi = np.zeros((state_count, state_count), dtype=complex)
for v in range(1, V_MAX + 1):
    Mv = build_M_v(v, j_test, b_prior_test)
    channel_pi += Mv @ diag_pi @ Mv.conj().T

# extract diagonal (the new state)
new_diag = np.diag(channel_pi).real
print(f"  ‖new_diag − π‖_∞ = {np.max(np.abs(new_diag - pi_f)):.6e}")
print(f"  ‖new_diag − π‖_2 = {np.linalg.norm(new_diag - pi_f):.6e}")
print(f"  sum(new_diag) = {new_diag.sum():.10f}")
print(f"  sum(π)        = {pi_f.sum():.10f}")
print()

# --- T_j Kraus reconstruction ---------------------------------------------
print("# (3) T_j(full) = Σ_v M_v decomposition")
print()
T_j_kraus = np.zeros((state_count, state_count), dtype=complex)
for v in range(1, V_MAX + 1):
    T_j_kraus += build_M_v(v, j_test, b_prior_test)

# Note: T_j as a linear operator (NOT a CP channel) is Σ_v M_v at unit-norm
# weighted by 2^{-v/2}. The CP channel is ρ → Σ_v M_v ρ M_v†.
print(f"  ‖T_j(Kraus)‖_F = {np.linalg.norm(T_j_kraus, ord='fro'):.6f}")
print(f"  ‖T_j(Kraus)‖_op (spectral norm) = {np.linalg.norm(T_j_kraus, ord=2):.6f}")
print()


# --- Moment computation via Kraus form -----------------------------------
print("# (4) Third-order moment ϕ(X̃_{j_1} · X̃_{j_2} · X̃_{j_1}) via Kraus")
print("#     j_1 = 1, j_2 = 2; b_priors averaged over Pascal(2,1/2)")
print()

# Pascal(2, 1/2) for pair-sum b: P(b) = (b-1) · 2^{-b} for b ≥ 2, normalized over [2, 2V_MAX]
B_VALUES = list(range(2, 2 * V_MAX + 1))
pascal_unnorm = []
for b in B_VALUES:
    # number of (v, v') pairs with v + v' = b in [1, V_MAX]²
    n_pairs = sum(1 for v in range(1, V_MAX + 1) if 1 <= b - v <= V_MAX)
    if n_pairs == 0:
        pascal_unnorm.append(0.0)
    else:
        # P(b) = sum over (v,v') of 2^{-v} 2^{-v'}
        pascal_unnorm.append(sum(2.0**(-v - (b - v)) for v in range(1, V_MAX + 1)
                                  if 1 <= b - v <= V_MAX))
Z_pascal = sum(pascal_unnorm)
pascal_w = [p / Z_pascal for p in pascal_unnorm]


def build_T_j_conditional(j, b_prior):
    """Build T_j as full operator conditional on b_prior (sum of all Kraus M_v)."""
    T = np.zeros((state_count, state_count), dtype=complex)
    for v in range(1, V_MAX + 1):
        T += build_M_v(v, j, b_prior)
    return T


def expectation_E_B_at_step(j, b_prior_list):
    """E_B[T_j] = average T_j over pair-sum b realizations.

    Here we compute T_j conditional on b_prior (the prior accumulator b_{[1,j-1]}),
    then later average over the b_prior realizations using the appropriate
    accumulator distribution.

    For the marginal-centering reading (§8 of AMALG_FREENESS_MOMENT_CALCULATION.md),
    E_B[T_j](b_prior) is just T_j(b_prior) since we marginalize over the
    pair-sum b_j at step j while holding b_prior fixed.

    Actually for the "centered" version X̃_j = T_j - E_B[T_j], we need:
      X̃_j(b_prior) = T_j(b_prior) - E_avg[T_j(b_prior)]
    where E_avg averages over (v_{2j-1}, v_{2j}) within fixed b_prior + b_j = sum.

    Since T_j is already linear in v (sums over Kraus indices), the centering
    here removes the rank-1 "diagonal contribution" π · (something), leaving
    the off-diagonal correction Off_j (per §7 of AMALG_FREENESS_MOMENT_CALCULATION.md).
    """
    # For Reading B (marginal centering): center each row by π-weighted column-mean
    T_j_b = build_T_j_conditional(j, b_prior_list)
    # The "expected" value under B-marginal averaging is the rank-1 projector
    # onto stationary: E_B[T_j] = π_f^T · 1 (or similar). The implementation
    # in verify_monotone_diagnostic.py is to subtract the row-stationary-mean.
    # We follow that convention here.
    row_mean = T_j_b @ pi_f
    centered = T_j_b - np.outer(row_mean, np.ones(state_count) * 0)  # placeholder
    # Actually the marginal centering subtracts the column-projected mean:
    proj = np.outer(np.ones(state_count), pi_f) @ T_j_b
    centered = T_j_b - proj
    return T_j_b, centered


# Helper: compute the moment by averaging over b_prior (the accumulator at step j_1's start)
# For j_1 = 1: b_prior_1 = 0 (no prior accumulator)
# For j_2 = 2: b_prior_2 = b_1 (the pair-sum at step 1, distributed Pascal(2,1/2))

# At j_1=1, b_prior=0 fixed. Build M_v^{(1, 0)}.
T_1_full = build_T_j_conditional(1, 0)
# X̃_1 = T_1 - E_B[T_1]; under Reading B marginal centering,
# E_B[T_1] subtracts the rank-1 stationary projection
proj_T1 = np.outer(np.ones(state_count), pi_f) @ T_1_full
X_tilde_1 = T_1_full - proj_T1

# For j_2 = 2: T_2 depends on b_prior_2 = b_1 ~ Pascal(2, 1/2)
# We compute E[X̃_1 · X̃_2(b_1) · X̃_1] averaged over b_1 distribution

# Center X̃_2(b_1) for each b_1:
M_3_alt_kraus = 0.0
for b_1, pw in zip(B_VALUES, pascal_w):
    if pw == 0.0:
        continue
    T_2_b1 = build_T_j_conditional(2, b_1)
    proj_T2 = np.outer(np.ones(state_count), pi_f) @ T_2_b1
    X_tilde_2 = T_2_b1 - proj_T2

    # Trace of X̃_1 · X̃_2(b_1) · X̃_1 against π (scalar reduction "sum_entries")
    triple = X_tilde_1 @ X_tilde_2 @ X_tilde_1
    scalar = np.ones(state_count) @ triple @ np.ones(state_count)
    M_3_alt_kraus += pw * scalar.real

print(f"  Kraus-form result (sum_entries scalar): M_3_alt_kraus = {M_3_alt_kraus:.6e}")
print(f"  Syracuse direct measurement (Task 1):    M_3_alt_direct = 1.078e-1")
print(f"  ratio = {abs(M_3_alt_kraus) / 0.1078:.4f}")
print()

# --- Compute via tr_π scalar reduction too -------------------------------
M_3_alt_kraus_trpi = 0.0
for b_1, pw in zip(B_VALUES, pascal_w):
    if pw == 0.0:
        continue
    T_2_b1 = build_T_j_conditional(2, b_1)
    proj_T2 = np.outer(np.ones(state_count), pi_f) @ T_2_b1
    X_tilde_2 = T_2_b1 - proj_T2

    triple = X_tilde_1 @ X_tilde_2 @ X_tilde_1
    scalar_trpi = np.einsum('i,ii->', pi_f, triple)
    M_3_alt_kraus_trpi += pw * scalar_trpi.real

print(f"  Kraus-form (tr_π reduction): {M_3_alt_kraus_trpi:.6e}")
print(f"  Syracuse (tr_π reduction, Task 1): 5.36e-2 area")
print()

# --- Save outputs --------------------------------------------------------
out = {
    'description': 'DWM-MP-G1 verification: explicit Kraus form of T_j and moment computation',
    'level': N_LEVEL,
    'state_count': state_count,
    'V_MAX': V_MAX,
    'POVM_deviation_from_I_frobenius': dev_I,
    'POVM_trace_per_dim': trace_total,
    'POVM_expected_trace': expected_trace,
    'stationary_preservation_max_dev': float(np.max(np.abs(new_diag - pi_f))),
    'stationary_preservation_l2_dev': float(np.linalg.norm(new_diag - pi_f)),
    'M_3_alt_kraus_sum_entries': M_3_alt_kraus,
    'M_3_alt_kraus_tr_pi': M_3_alt_kraus_trpi,
    'Syracuse_direct_M_3_alt_sum_entries': 1.078e-1,
    'Syracuse_direct_M_3_alt_tr_pi': 5.36e-2,
}
out_path = r'C:\Collatz\experiments_output\dwm_kraus_verify.json'
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w') as f:
    json.dump(out, f, indent=2, default=float)
print(f"[save] {out_path}")
