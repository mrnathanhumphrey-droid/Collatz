"""saddle_class_subsum_analysis.py — measure per-saddle-class subsums S_j(r).

For each r ∈ {6, 8, 10} and each saddle class j ∈ {0, 1, 2}:
  S_j(r) = |Σ_{a ∈ supp_j(r)} D(-3a) · ψ_true(a)|
where:
  D(-3a)  := Σ_{u=0}^{N-1} e_q(-3·a·u)        (Dirichlet kernel, length N=3^{r-1})
  ψ_true(a) := G(a)/√q                         (true phase, |ψ| = 1)
  G(a)    := Σ_{s=0}^{period-1} e_q(P_a(s))    (full inner Gauss sum, length 3^r)
  supp_j(r) := {a ∈ supp(r) : s*(C_a) = j},  s*(C_a) = (C_a − 1)/3 mod 3

Q1: power-law fit log(S_j) vs log(|supp_j|). Slope ≈ 0.5 ⟹ saturation; slope ≈ 0
or growth as log ⟹ controlled.

Q2 (only if Outcome A: j=1, j=2 controlled): characterize ψ on supp_0 — linearity,
periodicity in a-index, Fourier-mass concentration, |Σ ψ|/|supp_0| ratio.
"""
import sys, os, math, time, csv
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"C:\Collatz")
from r79b_scenario_comparison import (
    J_for_p3, truncated_3adic_log_mod_q, precompute_P_a_coeff_table,
    compute_dirichlet_arr, compute_G_arr, compute_psi_lead_arr,
)

PI = math.pi


def main():
    out_csv = r"C:\Collatz\saddle_class_subsum_data.csv"
    fieldnames = ["r", "j", "n_j", "S_j", "log_S_j", "log_n_j",
                  "S_j_over_sqrt_n_j", "psi_only_sum_abs",
                  "psi_only_over_n_j", "elapsed_s"]
    rows_for_csv = []

    print("="*72)
    print("# Q1: Per-saddle-class subsum S_j(r) at r ∈ {6, 8, 10}")
    print("="*72)
    print()

    # Warmup numba
    threes_, Ls_ = precompute_P_a_coeff_table(2, 27, J_for_p3(3))
    small = np.array([1, 4, 7], dtype=np.int64)
    _ = compute_G_arr(2, small, threes_, Ls_, 1, 9)
    _ = compute_dirichlet_arr(2, small)

    # Storage for per-r data (will be needed for Q2)
    full_data = {}

    for r in [6, 8, 10]:
        t0 = time.time()
        q = 3 ** (r + 1)
        p_mm1 = 3 ** r
        N = 3 ** (r - 1)
        J = J_for_p3(r + 1)
        sqrt_q = math.sqrt(q)

        print(f"## r = {r}: q = {q}, N = {N}, period = {3**r}, J = {J}")

        supp = np.array([a for a in range(p_mm1) if a % 3 == 1], dtype=np.int64)
        L4 = truncated_3adic_log_mod_q(1, J, q)
        L_tilde = L4 // 3
        L_tilde_inv = pow(L_tilde, -1, p_mm1)

        threes, Ls = precompute_P_a_coeff_table(r, q, J)

        # Compute G(a), then ψ_true = G/√q
        t_g = time.time()
        g_re, g_im = compute_G_arr(r, supp, threes, Ls, L_tilde_inv, p_mm1)
        print(f"   G(a) computed in {time.time()-t_g:.2f}s")
        psi_re = g_re / sqrt_q
        psi_im = g_im / sqrt_q

        # Verify |ψ_true| = 1
        psi_mag = np.sqrt(psi_re**2 + psi_im**2)
        max_dev = np.abs(psi_mag - 1.0).max()
        print(f"   max |ψ_true| − 1 = {max_dev:.2e}  (T78.3 saturation check)")

        # Compute D(-3a) for each a
        t_d = time.time()
        d_re, d_im = compute_dirichlet_arr(r, supp)
        print(f"   D(-3a) computed in {time.time()-t_d:.2f}s")

        # Determine s*(a) for each a
        s_star = np.array(
            [(((int(a) * L_tilde_inv) % p_mm1 - 1) // 3) % 3 for a in supp],
            dtype=np.int64,
        )

        # Per-class subsums
        print(f"   {'j':>3} {'n_j':>6} {'S_j':>14} {'log(S_j)':>10} {'log(n_j)':>10} "
              f"{'S_j/√n_j':>10} {'|Σ ψ_j|':>12} {'|Σ ψ_j|/n_j':>12}")
        for j in [0, 1, 2]:
            mask = (s_star == j)
            n_j = int(mask.sum())
            # S_j = |Σ D(-3a) · ψ(a)| over class j
            # Complex multiply: (d_re + i d_im) · (psi_re + i psi_im) = (d_re·psi_re − d_im·psi_im) + i(d_re·psi_im + d_im·psi_re)
            term_re = d_re[mask] * psi_re[mask] - d_im[mask] * psi_im[mask]
            term_im = d_re[mask] * psi_im[mask] + d_im[mask] * psi_re[mask]
            sum_re = float(term_re.sum())
            sum_im = float(term_im.sum())
            S_j = math.sqrt(sum_re ** 2 + sum_im ** 2)

            # Q2 prep: |Σ ψ_j| (without Dirichlet weighting)
            psi_only_re = float(psi_re[mask].sum())
            psi_only_im = float(psi_im[mask].sum())
            psi_only_abs = math.sqrt(psi_only_re ** 2 + psi_only_im ** 2)

            log_S = math.log(S_j) if S_j > 0 else float("-inf")
            log_nj = math.log(n_j)

            print(f"   {j:>3} {n_j:>6} {S_j:>14.4f} {log_S:>10.4f} {log_nj:>10.4f} "
                  f"{S_j/math.sqrt(n_j):>10.4f} {psi_only_abs:>12.4f} {psi_only_abs/n_j:>12.6f}")

            rows_for_csv.append({
                "r": r, "j": j, "n_j": n_j,
                "S_j": f"{S_j:.6f}",
                "log_S_j": f"{log_S:.6f}",
                "log_n_j": f"{log_nj:.6f}",
                "S_j_over_sqrt_n_j": f"{S_j/math.sqrt(n_j):.4f}",
                "psi_only_sum_abs": f"{psi_only_abs:.6f}",
                "psi_only_over_n_j": f"{psi_only_abs/n_j:.6f}",
                "elapsed_s": "",
            })

        # Save full per-a data for later structural tests on supp_0
        full_data[r] = {
            "supp": supp,
            "psi_re": psi_re,
            "psi_im": psi_im,
            "d_re": d_re,
            "d_im": d_im,
            "s_star": s_star,
            "L_tilde_inv": L_tilde_inv,
            "p_mm1": p_mm1,
            "q": q,
            "N": N,
        }
        print(f"   total elapsed: {time.time()-t0:.2f}s")
        print()

    # Power-law fit per j
    print("="*72)
    print("# Power-law fit: log(S_j) = a_j + β_j · log(n_j) per class")
    print("="*72)
    print()
    rs_test = [6, 8, 10]
    by_j = {0: [], 1: [], 2: []}
    for row in rows_for_csv:
        j = row["j"]
        log_S = float(row["log_S_j"])
        log_n = float(row["log_n_j"])
        by_j[j].append((log_n, log_S))

    print(f"  {'j':>3} {'slope β_j':>14} {'intercept a_j':>16} {'R²':>8} "
          f"{'interpretation':>50}")
    classification = {}
    for j in [0, 1, 2]:
        pts = by_j[j]
        xs = np.array([p[0] for p in pts])
        ys = np.array([p[1] for p in pts])
        n = len(pts)
        x_m = xs.mean()
        y_m = ys.mean()
        ssxy = ((xs - x_m) * (ys - y_m)).sum()
        ssxx = ((xs - x_m) ** 2).sum()
        beta = ssxy / ssxx if ssxx > 0 else 0
        alpha = y_m - beta * x_m
        y_pred = alpha + beta * xs
        ss_res = ((ys - y_pred) ** 2).sum()
        ss_tot = ((ys - y_m) ** 2).sum()
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

        if 0.45 <= beta <= 0.55:
            interp = "saturating (rate-1/2)"
        elif beta < 0.2:
            interp = "controlled (sub-poly, possibly O(log) or O(1))"
        elif beta < 0.45:
            interp = "intermediate (slow growth)"
        else:
            interp = f"fast growth (rate > 0.55)"
        classification[j] = (beta, interp)

        print(f"  {j:>3} {beta:>14.4f} {alpha:>16.4f} {r2:>8.4f} {interp:>50}")
    print()

    # Decision: A or B
    j1_beta = classification[1][0]
    j2_beta = classification[2][0]
    j0_beta = classification[0][0]

    print("="*72)
    print("# Decision: Outcome A vs B")
    print("="*72)
    if j1_beta < 0.45 and j2_beta < 0.45:
        outcome = "A"
        print("  Outcome A: j=1 and j=2 are controlled (sub-saturating).")
        print(f"    β_j=1 = {j1_beta:.4f}, β_j=2 = {j2_beta:.4f}")
        print(f"    Saddle-class partition is a real closure path.")
        print(f"    j=0 carries the load: β_j=0 = {j0_beta:.4f}")
    elif 0.45 <= j1_beta <= 0.55 and 0.45 <= j2_beta <= 0.55:
        outcome = "B"
        print("  Outcome B: all three classes saturate (rate-1/2 each).")
        print(f"    β_j=0 = {j0_beta:.4f}, β_j=1 = {j1_beta:.4f}, β_j=2 = {j2_beta:.4f}")
        print(f"    Saddle-class partition is structural but NOT a closure path.")
    else:
        outcome = "ambiguous"
        print(f"  Ambiguous: β_j=0={j0_beta:.4f}, β_j=1={j1_beta:.4f}, β_j=2={j2_beta:.4f}")
        print(f"  Need higher r to resolve.")
    print()

    # Save Q1 CSV
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows_for_csv:
            writer.writerow(row)
    print(f"[saved Q1 csv: {out_csv}]")
    print()

    # Q2: only if Outcome A
    if outcome == "A":
        print("="*72)
        print("# Q2: Structural tests on ψ over supp_0 (Outcome A — j=0 carries the load)")
        print("="*72)
        print()
        for r in [6, 8, 10]:
            d = full_data[r]
            mask0 = (d["s_star"] == 0)
            psi0_re = d["psi_re"][mask0]
            psi0_im = d["psi_im"][mask0]
            supp0 = d["supp"][mask0]
            n0 = len(supp0)
            print(f"## r = {r}, |supp_0| = {n0}")

            # Test 1: linear in a regression (real and imag separately).
            # Project to a-coordinate; fit linear ψ(a) ≈ α + β·a (mod 1 in phase).
            # Phase as angle in radians:
            angles = np.arctan2(psi0_im, psi0_re)  # in (-π, π]
            # Sort by a
            idx = np.argsort(supp0)
            a_sorted = supp0[idx]
            angles_sorted = angles[idx]

            # Fit: angle ≈ α + β·a (linear regression, ignoring 2π wrapping)
            # But angles wrap around — better to look at differences (angle_{i+1} - angle_i) mod 2π
            diffs = np.diff(angles_sorted)
            # wrap to (-π, π]
            diffs = (diffs + math.pi) % (2 * math.pi) - math.pi
            mean_step = float(np.mean(diffs))
            std_step = float(np.std(diffs))
            print(f"   linearity test: mean phase-step per a = {mean_step:+.4f} rad,  std = {std_step:.4f} rad")
            print(f"      (std ≪ |mean| ⟹ linear in a; std ≈ |mean| or larger ⟹ non-linear)")

            # Test 2: ψ Fourier transform over a-index (DFT over a as a sequence)
            # ψ as complex array (in supp_0 sorted-by-a order)
            psi0_sorted = (psi0_re[idx] + 1j * psi0_im[idx])
            psi_fft = np.fft.fft(psi0_sorted)
            mag = np.abs(psi_fft)
            top_idx = np.argsort(mag)[::-1][:5]
            print(f"   top 5 Fourier modes of ψ over a-index (mag, freq_idx, freq_normalized):")
            for k in top_idx:
                freq_norm = k / n0 if k <= n0 // 2 else (k - n0) / n0
                print(f"      mag = {mag[k]:>10.4f},  k = {k:>6} ,  k/n0 = {freq_norm:+.5f}")
            mass_top1 = float(mag[top_idx[0]])
            mass_total = float(np.sqrt((mag**2).sum()))  # = √(n_0 · n_0) = n_0 since |ψ|=1
            print(f"   top-1 mass / Parseval total ({mass_total:.2f}) = {mass_top1/mass_total:.4f}")

            # Test 3: |Σ ψ_0| / |supp_0|
            psi_sum = (psi0_re.sum() + 1j * psi0_im.sum())
            psi_sum_abs = abs(psi_sum)
            print(f"   |Σ_a ψ(a)| over supp_0 = {psi_sum_abs:.4f}, ratio to n_0 = {psi_sum_abs/n0:.6f}")
            print(f"      (small ratio ⟹ ψ cancels among itself; ratio → 1 ⟹ aligned)")

            # Test 4: cumulative sum behavior — random-walk-like vs structured
            cum = np.cumsum(psi0_sorted)
            cum_mag = np.abs(cum)
            cum_max = float(cum_mag.max())
            cum_end = float(cum_mag[-1])
            sqrt_n0 = math.sqrt(n0)
            print(f"   max|cumsum| = {cum_max:.2f}, end |cumsum| = {cum_end:.2f}, √n_0 = {sqrt_n0:.2f}")
            print(f"      max/√n_0 = {cum_max/sqrt_n0:.3f} (~1 ⟹ random-walk; >>1 ⟹ structured drift)")
            print()

    print("[done]")
    return outcome, classification


if __name__ == "__main__":
    main()
