"""result_77_5_d_R_identity_check.py — verify ‖R_k‖² = ‖d_{k+1}‖² exactly over Q, k=1..5.

Two identities to verify (over Q via fractions.Fraction):
  (A) ‖R_k‖² == ‖d_{k+1}‖²  where
        ‖R_k‖²       := Σ_{r'} (π_{k+1}(r') − T(π_k)(r'))²    (R77.5 definition)
        ‖d_{k+1}‖²   := Σ_{r'} π_{k+1}(r')² − (1/3) Σ_r π_k(r)²    (R74 definition)
  (B) ‖R_k‖² · 3^k == S_{k+1}/3  where  S_{k+1} := X_{k+1} − X_k, X_j := 3^j · Σ π_j²
        (R74's identity S_{k+1} = 3^{k+1} · ‖d_{k+1}‖² in different rearrangement)

Algebraic proof of (A):
  ‖R_k‖² = Σ π_{k+1}² − 2 Σ π_{k+1}·T(π_k) + Σ T(π_k)²
  Σ T(π_k)² = (1/9) · 3 · Σ π_k² = Σ π_k² / 3
  Σ π_{k+1}·T(π_k) = (1/3) Σ_r π_k(r) · [Σ_{r' lifts of r} π_{k+1}(r')]
                   = (1/3) Σ_r π_k(r) · π_k(r)    (marginal consistency)
                   = (1/3) Σ π_k²
  So ‖R_k‖² = Σ π_{k+1}² − (2/3) Σ π_k² + (1/3) Σ π_k² = Σ π_{k+1}² − (1/3) Σ π_k² = ‖d_{k+1}‖² ✓

This script verifies (A) and (B) as exact rational equalities at k = 1, 2, 3, 4, 5.
"""
import sys, os, csv
from fractions import Fraction

sys.path.insert(0, r"C:\Collatz")
from result_77_5_compute_R_k import pi_dict, lift_pi, squared_l2_norm

sys.stdout.reconfigure(encoding="utf-8")
OUT_CSV = r"C:\Collatz\result_77_5_d_R_norms.csv"


def main():
    print("# R77.5 follow-up: identity check ‖R_k‖² == ‖d_{k+1}‖² over Q\n")

    pis = {}
    for k in [1, 2, 3, 4, 5, 6]:
        pi, _ = pi_dict(k)
        pis[k] = pi
        print(f"  computed π_{k}: {len(pi)} states")
    print()

    # Compute X_k = 3^k · Σ π_k² for each k (R74's pre-S_k quantity)
    X = {0: Fraction(1)}  # X_0 := 1 (one trivial state)
    for k in [1, 2, 3, 4, 5, 6]:
        X[k] = Fraction(3 ** k) * sum(p * p for p in pis[k].values())

    # S_{k+1} = X_{k+1} − X_k
    S = {}
    for k in [1, 2, 3, 4, 5, 6]:
        S[k] = X[k] - X[k - 1]

    rows = []
    print(f"{'k':>2}  {'‖R_k‖² (decimal)':>20}  {'‖d_{k+1}‖² (decimal)':>22}  "
          f"{'(A) Test':>9}  {'‖R_k‖²·3^k':>14}  {'S_{k+1}/3':>14}  {'(B) Test':>9}")
    for k in [1, 2, 3, 4, 5]:
        # R_k as vector
        T_pi_k = lift_pi(pis[k], k)
        R_k = {rp: pis[k + 1][rp] - T_pi_k[rp] for rp in pis[k + 1]}
        norm_R = squared_l2_norm(R_k)

        # d_{k+1} squared via R74 definition
        sum_pi_k_sq = sum(p * p for p in pis[k].values())
        sum_pi_kp1_sq = sum(p * p for p in pis[k + 1].values())
        norm_d = sum_pi_kp1_sq - sum_pi_k_sq / 3

        test_A = (norm_R == norm_d)

        norm_R_scaled = norm_R * Fraction(3 ** k)
        S_kp1_over_3 = S[k + 1] / 3
        test_B = (norm_R_scaled == S_kp1_over_3)

        print(f"{k:>2}  {float(norm_R):>20.10e}  {float(norm_d):>22.10e}  "
              f"{'PASS' if test_A else 'FAIL':>9}  {float(norm_R_scaled):>14.10f}  "
              f"{float(S_kp1_over_3):>14.10f}  {'PASS' if test_B else 'FAIL':>9}")

        rows.append({
            "k": k,
            "norm_R_num": norm_R.numerator, "norm_R_den": norm_R.denominator,
            "norm_d_num": norm_d.numerator, "norm_d_den": norm_d.denominator,
            "test_A": test_A,
            "norm_R_3k_num": norm_R_scaled.numerator, "norm_R_3k_den": norm_R_scaled.denominator,
            "S_kp1_3_num": S_kp1_over_3.numerator, "S_kp1_3_den": S_kp1_over_3.denominator,
            "test_B": test_B,
            "decimal_R_3k": float(norm_R_scaled),
        })

    # Limit: 7/45
    target = Fraction(7, 45)
    print(f"\n  target as k → ∞:  7/45 = {float(target):.10f}")
    print(f"  ‖R_5‖² · 3^5    = {float(rows[-1]['decimal_R_3k']):.10f}")
    print(f"  diff             = {float(Fraction(rows[-1]['norm_R_3k_num'], rows[-1]['norm_R_3k_den']) - target):+.10e}")
    print()

    # CSV output
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "k",
            "norm_R_sq_num", "norm_R_sq_den",
            "norm_d_sq_num", "norm_d_sq_den",
            "test_A_norm_R_eq_norm_d",
            "norm_R_sq_times_3k_num", "norm_R_sq_times_3k_den",
            "S_kp1_over_3_num", "S_kp1_over_3_den",
            "test_B_norm_R_3k_eq_S_kp1_3",
            "decimal_R_sq_3k",
        ])
        for r in rows:
            w.writerow([
                r["k"],
                r["norm_R_num"], r["norm_R_den"],
                r["norm_d_num"], r["norm_d_den"],
                r["test_A"],
                r["norm_R_3k_num"], r["norm_R_3k_den"],
                r["S_kp1_3_num"], r["S_kp1_3_den"],
                r["test_B"],
                f"{r['decimal_R_3k']:.16e}",
            ])
    print(f"[save] {OUT_CSV}")

    n_pass_A = sum(1 for r in rows if r["test_A"])
    n_pass_B = sum(1 for r in rows if r["test_B"])
    print(f"\n  (A) ‖R_k‖² == ‖d_{{k+1}}‖²  : {n_pass_A}/5 passed")
    print(f"  (B) ‖R_k‖² · 3^k == S_{{k+1}}/3 : {n_pass_B}/5 passed")
    if n_pass_A == 5 and n_pass_B == 5:
        print("\n  Outcome: (IDENTITY). R_k = d_{k+1} algebraically; R77.5 = R74 in geometric basis.")


if __name__ == "__main__":
    main()
