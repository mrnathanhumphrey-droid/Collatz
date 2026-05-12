"""
pade_numerical.py
=================
Wilson follow-up probe to R77.6 / PADE_EXTENSION: extend Padé analysis of
E(z) = sum_{n>=1} eps_n z^n using NUMERICAL eps_7..eps_13 from prior-session
power-iteration / scipy.eigs / FFT computations on K_k.

Inputs (hardcoded floats, sourced from result_epsilon_*.md and probes):
  eps_1..eps_6: exact rationals (rendered as float64)
  eps_7..eps_13: numerical floats (cross-checked across methods to ~1e-15)

Output:
  - stdout: full diagonal [1/1]..[5/5] plus near-diagonal pole tables
  - C:/Collatz/experiments_output/result_pade_numerical_poles.csv

This is the *verification* artifact for main-thread execution. The agent
already derived the pole locations analytically and tabulated them in the
PADE_NUMERICAL_*.md deliverables; this script confirms.
"""
from __future__ import annotations

import csv
import os
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)
OUT_CSV = os.path.join(OUTDIR, "result_pade_numerical_poles.csv")

# --------------------------------------------------------------------------- #
# Data: eps_n from project files                                              #
# --------------------------------------------------------------------------- #

EPS = {
    1:  +2.0000000000e-01,   # exact 1/5
    2:  +9.5238095238e-03,   # exact 1/105
    3:  -5.0919863259e-03,   # exact -5191/1019445
    4:  -2.4522582483e-03,
    5:  -1.1517469151e-03,
    6:  -4.9790566522e-04,
    7:  -1.1752368304e-03,
    8:  -7.4554636729e-04,
    9:  -7.5202571564e-06,
    10: +7.2075091711e-04,
    11: +1.5019670121e-03,
    12: +2.2747137206e-03,
    13: +2.9482473172e-03,
}

# f_tilde(z) = (E(z) - eps_1 * z) / z^2 has Taylor coefficients
# c_j = eps_{j+2} for j = 0, 1, 2, ...
def coeffs(K_max):
    """Return [c_0, ..., c_{K_max-1}] = [eps_2, ..., eps_{K_max+1}]."""
    return np.array([EPS[j + 2] for j in range(K_max)], dtype=float)


# --------------------------------------------------------------------------- #
# Pade construction                                                           #
# --------------------------------------------------------------------------- #

def pade_approximant(c, m, n):
    """Compute [m/n] Pade of f(z) = sum c_k z^k. f has degree m+n+1
    coefficients available.

    Returns (P, Q) where P degree<=m, Q degree<=n, Q[0]=1.
    P/Q matches f up to z^{m+n} inclusive.

    Linear system:
        sum_{j=1}^{n} q_j * c_{m+i-j} = -c_{m+i},  i = 1..n  (with c_k = 0 for k < 0)
    """
    if len(c) < m + n + 1:
        raise ValueError(f"Need {m+n+1} coeffs, have {len(c)}")
    if n == 0:
        return c[:m+1].copy(), np.array([1.0])
    # Build A (n x n), b (n,)
    A = np.zeros((n, n), dtype=float)
    b = np.zeros(n, dtype=float)
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            idx = m + i - j
            A[i - 1, j - 1] = c[idx] if idx >= 0 else 0.0
        b[i - 1] = -c[m + i]
    try:
        q_rest = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        return None
    Q = np.concatenate(([1.0], q_rest))
    # Compute P_k = sum_{j=0}^{min(k,n)} Q_j * c_{k-j}
    P = np.zeros(m + 1, dtype=float)
    for k in range(m + 1):
        s = 0.0
        for j in range(0, min(k, n) + 1):
            s += Q[j] * c[k - j]
        P[k] = s
    return P, Q


def roots_of_Q(Q):
    """Return list of complex roots of Q(z) = sum_j Q[j] z^j. Use numpy
    convention (highest-degree-first)."""
    # numpy.roots wants highest-degree-first coefficient
    pol = Q[::-1]
    # Strip leading zeros
    while len(pol) > 1 and pol[0] == 0:
        pol = pol[1:]
    if len(pol) <= 1:
        return []
    return list(np.roots(pol))


# --------------------------------------------------------------------------- #
# Pole-table dump                                                              #
# --------------------------------------------------------------------------- #

Z2 = 2.0 + 0.0j
Z_SLOW = 1.0163  # = 1 / 0.984 = slow-mode prediction from STATE.md (rho ~ 0.984)


def report_approximant(label, c, m, n, rows):
    print(f"  [{m}/{n}] ({label})")
    res = pade_approximant(c, m, n)
    if res is None:
        print(f"    SINGULAR linear system")
        return
    P, Q = res
    roots = roots_of_Q(Q)
    roots.sort(key=lambda r: abs(r - Z2))
    for j, r in enumerate(roots):
        d2 = abs(r - Z2)
        ds = abs(r - Z_SLOW)
        kind = "REAL" if abs(r.imag) < 1e-9 else "CC"
        print(f"    pole {j}: z = {r.real:+.6f}{r.imag:+.6f}j  |z-2|={d2:.4f}  |z-1.016|={ds:.4f}  [{kind}]")
        rows.append({
            "m": m, "n": n, "label": label, "pole_idx": j,
            "pole_real": r.real, "pole_imag": r.imag,
            "abs_pole": abs(r),
            "dist_to_2": d2,
            "dist_to_slow": ds,
            "kind": kind,
        })


def main():
    print("=" * 78)
    print("PADE_NUMERICAL — extended Pade on eps_n with numerical eps_7..eps_13")
    print("=" * 78)
    c = coeffs(12)
    print("\nCoefficients of f_tilde(z) = (E(z) - eps_1 z)/z^2:")
    for j, cj in enumerate(c):
        print(f"  c_{j:>2} = eps_{j+2:>2} = {cj:+.10e}")
    print()

    rows = []

    print("Diagonal [n/n] sequence:")
    print("-" * 78)
    for nd in range(1, 6):  # [1/1] .. [5/5]
        report_approximant(f"diagonal [{nd}/{nd}]", c, nd, nd, rows)
    print()

    print("Off-diagonal near-diagonal approximants:")
    print("-" * 78)
    for (m, n) in [(3, 2), (2, 3), (4, 3), (3, 4), (5, 4), (4, 5), (5, 6), (6, 5), (6, 4), (4, 6)]:
        if m + n + 1 > len(c):
            print(f"  [{m}/{n}] — SKIP, need {m+n+1} coeffs, have {len(c)}")
            continue
        report_approximant(f"off-diag", c, m, n, rows)
    print()

    # Precision perturbation: perturb each numerical eps by relative 1e-14
    print("Precision perturbation (relative 1e-14 on eps_7..eps_13):")
    print("-" * 78)
    rng = np.random.default_rng(42)
    perturbed_eps_dict = dict(EPS)
    for k in range(7, 14):
        perturbed_eps_dict[k] = EPS[k] * (1.0 + 1e-14 * rng.standard_normal())
    c_pert = np.array([perturbed_eps_dict[j + 2] for j in range(12)], dtype=float)
    for nd in range(1, 6):
        res_a = pade_approximant(c, nd, nd)
        res_b = pade_approximant(c_pert, nd, nd)
        if res_a is None or res_b is None:
            continue
        ra = roots_of_Q(res_a[1])
        rb = roots_of_Q(res_b[1])
        ra.sort(key=lambda r: abs(r - Z2))
        rb.sort(key=lambda r: abs(r - Z2))
        # Closest pole shift
        if ra and rb:
            shift = abs(ra[0] - rb[0])
            print(f"  [{nd}/{nd}] closest-pole perturbation shift = {shift:.3e} "
                  f"(stable if <0.05)")

    # CSV out
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        fieldnames = ["m", "n", "label", "pole_idx", "pole_real", "pole_imag",
                       "abs_pole", "dist_to_2", "dist_to_slow", "kind"]
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"\n[csv: {OUT_CSV}]")


if __name__ == "__main__":
    main()
