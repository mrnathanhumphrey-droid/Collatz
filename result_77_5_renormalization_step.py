"""
result_77_5_renormalization_step.py — Stage 2 of R77.5

For each k = 2..5:
  - lift T_{k-1 -> k}(R_{k-1})
  - regression coeff c_k = <R_k, T(R_{k-1})> / <T(R_{k-1}), T(R_{k-1})>
  - orthogonal residual R_k^perp = R_k - c_k * T(R_{k-1})
  - report norm fractions

If c_k stabilizes across k, the inter-level renormalization Phi acts
diagonally on the leading mode with eigenvalue c. If c_∞^2 ~ 1/4 then
‖R_k‖^2 ratio ~ 1/4 ⇒ rate-1/2 envelope.

Reuses Stage 1 module (result_77_5_compute_R_k.py).

Outputs:
  C:\\Collatz\\result_77_5_phi_correlations.csv
"""
import sys
import os
import csv
from fractions import Fraction

sys.path.insert(0, r"C:\Collatz")
from result_77_5_compute_R_k import (
    pi_dict, lift_pi, sum_entries, squared_l2_norm
)

sys.stdout.reconfigure(encoding="utf-8")
OUT_CSV = r"C:\Collatz\result_77_5_phi_correlations.csv"


def lift_residual(R_km1, k_prev):
    """Lift R_{k-1} (vector on coprime Z/3^{k_prev+1}) to coprime Z/3^{k_prev+2}.

    R_{k_prev} = π_{k_prev+1} − T(π_{k_prev}), so R_{k_prev} lives on Z/3^{k_prev+1}.
    Lift target matches level of R_{k_prev+1} which lives on Z/3^{k_prev+2}.
    """
    Nkm1 = 3 ** (k_prev + 1)
    Nk = 3 ** (k_prev + 2)
    third = Fraction(1, 3)
    out = {}
    for rp in range(Nk):
        if rp % 3 == 0:
            continue
        r = rp % Nkm1
        out[rp] = R_km1[r] * third
    return out


def inner_product(u, v):
    """Σ_r u[r] * v[r] over the common index set."""
    s = Fraction(0)
    for r in u:
        s += u[r] * v[r]
    return s


def main():
    print("# R77.5 Stage 2: regress R_k onto T(R_{k-1})")
    print()

    # Compute pi_k and R_k for k = 1..5 (need 1..5 for R_k, 6 for k=5 case)
    pi_at = {}
    coprime_at = {}
    for k in [1, 2, 3, 4, 5, 6]:
        pi, cop = pi_dict(k)
        pi_at[k] = pi
        coprime_at[k] = cop
        print(f"  pi_{k}: {len(cop)} states")

    R = {}
    norms_sq = {}
    for k in [1, 2, 3, 4, 5]:
        T_pi_k = lift_pi(pi_at[k], k)
        R_k = {}
        for rp in pi_at[k + 1]:
            R_k[rp] = pi_at[k + 1][rp] - T_pi_k[rp]
        R[k] = R_k
        norms_sq[k] = squared_l2_norm(R_k)
        print(f"  R_{k}: ‖R_{k}‖² ≈ {float(norms_sq[k]):.6e}")

    print()
    print("# Compute c_k = <R_k, T(R_{k-1})> / <T(R_{k-1}), T(R_{k-1})>")
    rows = []
    for k in [2, 3, 4, 5]:
        T_R_prev = lift_residual(R[k - 1], k - 1)
        # T_R_prev lives on coprime Z/3^k same as R_k
        num = inner_product(R[k], T_R_prev)
        den = squared_l2_norm(T_R_prev)
        if den == 0:
            print(f"  k={k}: T(R_{k-1}) is zero — skipping.")
            continue
        c_k = num / den

        # Orthogonal residual
        R_perp = {r: R[k][r] - c_k * T_R_prev[r] for r in R[k]}
        norm_perp_sq = squared_l2_norm(R_perp)
        frac_explained = 1 - float(norm_perp_sq / norms_sq[k])

        # Note: T(R_{k-1}) has the same scaled-by-1/3 norm property:
        # ‖T(R_{k-1})‖^2 = (1/9) * 3 * ‖R_{k-1}‖^2 = ‖R_{k-1}‖^2 / 3
        # because each entry scales by 1/3 and there are 3x more entries
        # each contributing the same R_{k-1}[r]/3 squared.
        check_T = norms_sq[k - 1] / 3
        assert den == check_T, f"T-norm mismatch at k={k}: {den} vs {check_T}"

        print(f"  k={k}:")
        print(f"     <R_k, T(R_{k-1})>  = {num.numerator}/{num.denominator}  ≈ {float(num):+.6e}")
        print(f"     ‖T(R_{k-1})‖²      = {den.numerator}/{den.denominator}  ≈ {float(den):+.6e}")
        print(f"     c_k                = {c_k.numerator}/{c_k.denominator}  ≈ {float(c_k):+.6f}")
        print(f"     ‖R_k‖²             ≈ {float(norms_sq[k]):+.6e}")
        print(f"     ‖R_k^perp‖²        ≈ {float(norm_perp_sq):+.6e}")
        print(f"     frac_explained     ≈ {frac_explained:.6f}")
        print(f"     c_k^2              ≈ {float(c_k * c_k):.6f}  (1/4 target = 0.25)")

        rows.append({
            'k': k,
            'c_num': str(c_k.numerator),
            'c_den': str(c_k.denominator),
            'c_decimal': float(c_k),
            'R_k_norm_sq': float(norms_sq[k]),
            'R_k_perp_norm_sq': float(norm_perp_sq),
            'frac_explained': frac_explained,
            'c_k_squared': float(c_k * c_k),
        })

    print()
    with open(OUT_CSV, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["k", "c_k_num", "c_k_den", "c_k_decimal",
                    "R_k_norm_sq_decimal", "R_k_perp_norm_sq_decimal",
                    "frac_explained", "c_k_squared"])
        for row in rows:
            w.writerow([
                row['k'], row['c_num'], row['c_den'],
                f"{row['c_decimal']:.15e}",
                f"{row['R_k_norm_sq']:.15e}",
                f"{row['R_k_perp_norm_sq']:.15e}",
                f"{row['frac_explained']:.10f}",
                f"{row['c_k_squared']:.10f}",
            ])
    print(f"[save] {OUT_CSV}")


if __name__ == "__main__":
    main()
