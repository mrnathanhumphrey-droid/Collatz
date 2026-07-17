"""
d2_extract_phi_X_squared.py — D2 follow-up.

Compute the single-index second moments

    phi(X_tilde_j^2)  for j = j_1, j_2

at level n = 3, state space (Z/27)*, V_MAX=16, under Reading B
(marginal centering X_tilde_j = Off_j - E[Off_j | b_prior]) — the same
operator realization and reading used in
  C:/Collatz/verify_n4_alternating.py
  C:/Collatz/verify_monotone_diagnostic.py

These second moments are needed to apply the BMT / bigraph predicted
moment formulas to Syracuse rows (b), (d), (f):
  Row (f) under TENSOR-digraph BMT or bigraph with tensor edge (j_1,j_2):
      phi(a_1 a_2 a_3 a_4) = phi(X_tilde_{j_1}^2) . phi(X_tilde_{j_2}^2)

For j_1 we evaluate at b_prior_1 = 0 (root of accumulator chain) and
average over the within-pair split (v_1, v_1').

For j_2 we marginalize:
  - over the j_2 within-pair split (v_2, v_2'),
  - over the j_2 b_prior = b_1 = v_1+v_1' (since X_tilde_2 carries a
    b_prior-dependent phase, its squared expectation depends on b_prior;
    the natural single-index phi(X_tilde_2^2) for moment-formula comparison
    is the b_prior-averaged value).

We report under all 4 scalar reductions (tr_pi, vac_pi, delta_1, sum_entries).

OUTPUT: C:/Collatz/experiments_output/d2_phi_X_squared.json
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

N_LEVEL = 3
N = 3 ** N_LEVEL
TWO_PI_I_OVER_N = 2j * cmath.pi / N
V_MAX = 16

geom_weights_unnorm = [Fraction(1, 2 ** v) for v in range(1, V_MAX + 1)]
Z_geom = sum(geom_weights_unnorm)
geom_weights = [w / Z_geom for w in geom_weights_unnorm]
inv2_mod27 = pow(2, -1, N)
pow_inv2 = [pow(inv2_mod27, v, N) for v in range(0, 8 * V_MAX + 2)]

print("# d2_extract_phi_X_squared.py")
print(f"# level n = {N_LEVEL}; state space (Z/{N})*; V_MAX = {V_MAX}")
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
    """Off_j realization at split (v,vp), prior accumulator b_prior; WITHOUT 2^{-v-vp}."""
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
    M = np.zeros((state_count, state_count), dtype=complex)
    for v in range(1, V_MAX + 1):
        for vp in range(1, V_MAX + 1):
            if v == vp:
                continue
            w_raw = float(Fraction(1, 2 ** v * 2 ** vp))
            M += w_raw * build_Off_op_unweighted(j, v, vp, b_prior)
    return M


# scalar reductions
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
    return complex(np.conj(pi_f) @ P @ pi_f)


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


def compute_phi_X1_squared(scalar_fn, V_TRUNC=V_MAX):
    """phi(X_tilde_1^2) = E_{(v_1,v_1')} [ scalar( X_tilde_1 . X_tilde_1 ) ]

    X_tilde_1 = Off_1(v_1, v_1'; b_prior=0) - E_{(v,v')}[Off_1; b_prior=0]
    Outer sum: over (v_1, v_1') with weight 2^{-v_1-v_1'}. SAME pair at both positions.
    """
    Off_1_mean_0 = Off_unconditional_mean(1, 0)
    total = 0 + 0j
    for v1 in range(1, V_TRUNC + 1):
        for vp1 in range(1, V_TRUNC + 1):
            if v1 == vp1:
                continue
            w1 = float(Fraction(1, 2 ** v1 * 2 ** vp1))
            X1 = build_Off_op_unweighted(1, v1, vp1, 0) - Off_1_mean_0
            total += w1 * scalar_fn(X1, X1)
    return total


def compute_phi_X2_squared(scalar_fn, V_TRUNC=V_MAX):
    """phi(X_tilde_2^2) = E_{(v_2,v_2'; b_prior=b_1)} [ scalar( X_tilde_2 . X_tilde_2 ) ]
    averaged over b_1 = v_1+v_1' under (v_1, v_1') marginal (this is the
    "global" single-index phi(X_2^2) consistent with Reading B over the
    full prior accumulator).
    """
    total = 0 + 0j
    for v1 in range(1, V_TRUNC + 1):
        for vp1 in range(1, V_TRUNC + 1):
            if v1 == vp1:
                continue
            w_outer = float(Fraction(1, 2 ** v1 * 2 ** vp1))
            b1 = v1 + vp1
            Off_2_mean = Off_unconditional_mean(2, b1)
            inner = 0 + 0j
            for v2 in range(1, V_TRUNC + 1):
                for vp2 in range(1, V_TRUNC + 1):
                    if v2 == vp2:
                        continue
                    w2 = float(Fraction(1, 2 ** v2 * 2 ** vp2))
                    X2 = build_Off_op_unweighted(2, v2, vp2, b1) - Off_2_mean
                    inner += w2 * scalar_fn(X2, X2)
            total += w_outer * inner
    return total


def _main():
    print("=" * 70)
    print("# Computing phi(X_tilde_1^2) and phi(X_tilde_2^2) under Reading B")
    print(f"# V_TRUNC = {V_MAX}")
    print("=" * 70)

    results = {}
    for name, fn in SCALARS:
        t0 = time.time()
        phi_X1_sq = compute_phi_X1_squared(fn, V_TRUNC=V_MAX)
        t1 = time.time() - t0

        t0 = time.time()
        phi_X2_sq = compute_phi_X2_squared(fn, V_TRUNC=V_MAX)
        t2 = time.time() - t0

        product = phi_X1_sq * phi_X2_sq
        results[name] = {
            "phi_X1_squared": phi_X1_sq,
            "phi_X2_squared": phi_X2_sq,
            "product": product,
        }
        print(f"  [{name}]")
        print(f"    phi(X1^2) = {phi_X1_sq.real:+.6e}{phi_X1_sq.imag:+.6e}j  "
              f"|.| = {abs(phi_X1_sq):.6e}  ({t1:.1f}s)")
        print(f"    phi(X2^2) = {phi_X2_sq.real:+.6e}{phi_X2_sq.imag:+.6e}j  "
              f"|.| = {abs(phi_X2_sq):.6e}  ({t2:.1f}s)")
        print(f"    product   = {product.real:+.6e}{product.imag:+.6e}j  "
              f"|.| = {abs(product):.6e}")
        print()

    # Comparison to Syracuse measured row (f)
    print("=" * 70)
    print("# Comparison: BMT-tensor-digraph / bigraph(tensor edge) prediction for row (f)")
    print("# Predicted: phi(X1^2) . phi(X2^2)")
    print("# Measured (D1): M_4_alt = phi(X_tilde_1 . X_tilde_2 . X_tilde_1 . X_tilde_2)")
    print("=" * 70)

    # Load D1 measured values
    D1_M4_alt = {
        "tr_pi": 5.357224817678e-02,
        "vac_pi": 4.775478985215e-03,
        "delta_1": 5.742025710034e-02,
        "sum_entries": 6.088793223229e-01,
    }
    for name, _ in SCALARS:
        predicted = abs(results[name]["product"])
        measured = D1_M4_alt[name]
        ratio = predicted / measured if measured > 0 else float("inf")
        gap = abs(predicted - measured)
        print(f"  [{name}]")
        print(f"    predicted (phi(X1^2).phi(X2^2)) = {predicted:.6e}")
        print(f"    measured  (M_4_alt) ............ = {measured:.6e}")
        print(f"    ratio predicted/measured ....... = {ratio:.4f}")
        print(f"    abs gap ........................ = {gap:.6e}")

    out = {
        "task_id": "D2_phi_X_squared",
        "level": N_LEVEL,
        "state_count": state_count,
        "V_MAX": V_MAX,
        "centering_reading": (
            "Reading B (marginal centering): X_tilde_j = Off_j - E[Off_j | b_prior]. "
            "Same realization as verify_n4_alternating.py."
        ),
        "results_per_scalar_reduction": {
            name: {
                "phi_X1_squared": {
                    "real": results[name]["phi_X1_squared"].real,
                    "imag": results[name]["phi_X1_squared"].imag,
                    "abs": abs(results[name]["phi_X1_squared"]),
                },
                "phi_X2_squared": {
                    "real": results[name]["phi_X2_squared"].real,
                    "imag": results[name]["phi_X2_squared"].imag,
                    "abs": abs(results[name]["phi_X2_squared"]),
                },
                "product_phi_X1sq_phi_X2sq": {
                    "real": results[name]["product"].real,
                    "imag": results[name]["product"].imag,
                    "abs": abs(results[name]["product"]),
                },
                "D1_measured_M_4_alt": D1_M4_alt[name],
                "ratio_predicted_over_measured": (
                    abs(results[name]["product"]) / D1_M4_alt[name]
                    if D1_M4_alt[name] > 0 else None
                ),
            }
            for name, _ in SCALARS
        },
        "note_on_BMT_monotone_digraph": (
            "Under monotone digraph (j_2, j_1) in E since j_1<j_2: kerG[i=(j_1,j_2,j_1,j_2)]: "
            "1~3 since (i_2, i_1)=(j_2,j_1) IS edge; 2~4 requires (i_3, i_2)=(j_1,j_2) edge, "
            "but (j_1,j_2) NOT in E (would need j_2<j_1). So kerG[i] = {{1,3},{2},{4}} "
            "and BMT predicts phi(X1^2).phi(X2).phi(X4) = phi(X1^2).0.0 = 0. "
            "NOT phi(X1^2).phi(X2^2) as the auditor claimed."
        ),
        "note_on_BMT_tensor_digraph": (
            "Under tensor (complete) digraph: kerG[i] = ker[i] = {{1,3},{2,4}}. "
            "BMT predicts phi(X1.X3).phi(X2.X4) = phi(X1^2).phi(X2^2). Numerically tested above."
        ),
        "note_on_row_b": (
            "Row (b) m=2, i=(j_1,j_2): kerG[i]={{1},{2}} under both monotone and tensor digraph. "
            "BMT predicts phi(a_1).phi(a_2) = 0 . 0 = 0. Matches Syracuse ~0."
        ),
        "note_on_row_d": (
            "Row (d) m=3, i=(j_1,j_2,j_1): "
            "Monotone digraph: kerG[i]={{1,3},{2}} (1~3 via edge (j_2,j_1) in E). "
            "BMT predicts phi(X1.X3).phi(X2) = phi(X1^2).0 = 0. MISMATCH with measured 0.108. "
            "Tensor digraph: kerG[i]=ker[i]={{1,3},{2}}, same prediction = 0. MISMATCH."
        ),
    }
    out_path = r"C:\Collatz\experiments_output\d2_phi_X_squared.json"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print()
    print(f"[save] {out_path}")


if __name__ == "__main__":
    _main()
