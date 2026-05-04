"""R79b s-class structure of (G(a)/√q − ψ_lead(a)) deviation.

For each r ∈ {4, 6, 8, 10}: partition supp by s_star(a) ∈ {0, 1, 2} and compute
mean + variance of complex deviation per class. If class means differ → partition
is structurally meaningful with class-specific corrections; if not → leading-order
partition is an artifact at finite r.

Output: per-class mean(deviation), |mean|, std of |deviation|, n.
"""
import sys, math, time, csv
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

# Reuse machinery from r79b_scenario_comparison.py
sys.path.insert(0, r"C:\Collatz")
from r79b_scenario_comparison import (
    J_for_p3, truncated_3adic_log_mod_q, precompute_P_a_coeff_table,
    compute_dirichlet_arr, compute_G_arr, compute_psi_lead_arr,
)

PI = math.pi


def main():
    out_csv = r"C:\Collatz\r79b_s_class_deviation.csv"
    fieldnames = ["r", "j", "n_j", "mean_dev_re", "mean_dev_im",
                  "abs_mean_dev", "mean_abs_dev", "std_abs_dev", "max_abs_dev"]
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        # Warmup numba
        threes_, Ls_ = precompute_P_a_coeff_table(2, 27, J_for_p3(3))
        small = np.array([1, 4, 7], dtype=np.int64)
        _ = compute_G_arr(2, small, threes_, Ls_, 1, 9)
        _ = compute_psi_lead_arr(2, small, 1, 9, threes_, Ls_)

        print("="*72)
        print("# s-class structure of (G(a)/√q − ψ_lead(a)) deviation")
        print("# Per-class j ∈ {0, 1, 2}: mean(D_j), |mean(D_j)|, mean(|D_j|), std(|D_j|), n")
        print("="*72)
        print()

        for r in [4, 6, 8, 10]:
            t0 = time.time()
            q = 3 ** (r + 1)
            p_mm1 = 3 ** r
            J = J_for_p3(r + 1)
            sqrt_q = math.sqrt(q)
            print(f"## r = {r}: q = {q}, |supp| = {p_mm1//3 if p_mm1 % 3 == 0 else 0}, J = {J}")

            # Build supp + s_star mapping
            supp = np.array([a for a in range(p_mm1) if a % 3 == 1], dtype=np.int64)
            L4 = truncated_3adic_log_mod_q(1, J, q)
            L_tilde = L4 // 3
            L_tilde_inv = pow(L_tilde, -1, p_mm1)

            threes, Ls = precompute_P_a_coeff_table(r, q, J)

            # Compute G, ψ_lead
            g_re, g_im = compute_G_arr(r, supp, threes, Ls, L_tilde_inv, p_mm1)
            psi_lead_re, psi_lead_im = compute_psi_lead_arr(r, supp, L_tilde_inv, p_mm1, threes, Ls)

            # ψ_true = G / √q
            psi_true_re = g_re / sqrt_q
            psi_true_im = g_im / sqrt_q

            # Deviation D(a) = ψ_true − ψ_lead  (complex)
            D_re = psi_true_re - psi_lead_re
            D_im = psi_true_im - psi_lead_im

            # Determine s_star(a) for each a (mod 9 of C_a)
            s_star = np.array(
                [(((int(a) * L_tilde_inv) % p_mm1 - 1) // 3) % 3 for a in supp],
                dtype=np.int64,
            )

            print(f"   {'j':>3} {'n_j':>5} {'mean(D)_re':>12} {'mean(D)_im':>12} "
                  f"{'|mean D|':>10} {'mean|D|':>10} {'std|D|':>10} {'max|D|':>10}")
            for j in [0, 1, 2]:
                mask = (s_star == j)
                n_j = int(mask.sum())
                if n_j == 0:
                    continue
                D_re_j = D_re[mask]
                D_im_j = D_im[mask]
                mean_D_re = float(D_re_j.mean())
                mean_D_im = float(D_im_j.mean())
                abs_mean_D = math.sqrt(mean_D_re ** 2 + mean_D_im ** 2)
                abs_D = np.sqrt(D_re_j ** 2 + D_im_j ** 2)
                mean_abs_D = float(abs_D.mean())
                std_abs_D = float(abs_D.std())
                max_abs_D = float(abs_D.max())

                print(f"   {j:>3} {n_j:>5} {mean_D_re:>+12.6f} {mean_D_im:>+12.6f} "
                      f"{abs_mean_D:>10.6f} {mean_abs_D:>10.6f} {std_abs_D:>10.6f} {max_abs_D:>10.6f}")
                writer.writerow({
                    "r": r, "j": j, "n_j": n_j,
                    "mean_dev_re": f"{mean_D_re:.8f}",
                    "mean_dev_im": f"{mean_D_im:.8f}",
                    "abs_mean_dev": f"{abs_mean_D:.8f}",
                    "mean_abs_dev": f"{mean_abs_D:.8f}",
                    "std_abs_dev": f"{std_abs_D:.8f}",
                    "max_abs_dev": f"{max_abs_D:.8f}",
                })
            print(f"   elapsed {time.time()-t0:.2f}s")
            print()
            f.flush()

    print(f"\n[done] csv = {out_csv}")


if __name__ == "__main__":
    main()
