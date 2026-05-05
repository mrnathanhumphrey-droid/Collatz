"""
d_scout_k5_a4_b3.py
===================
Scout probe for joint 2-3-adic deviation matrix at depth k=5 of the
Syracuse 3x+1 trajectory ensemble. Computes
    H(r_2, r_3)  = empirical joint histogram of (T^5(n) mod 16, T^5(n) mod 27)
                   over N = 10^6 odd starts uniformly in [1, 10^12]
    M2(r_2)      = marginal mod 16 (8 odd residues)
    M3(r_3)      = marginal mod 27 (18 coprime-to-3 residues)
    D(r_2, r_3)  = H/N_eff - M2*M3                (signed deviation from CRT
                                                   product structure)
    chi2         = sum_{r_2,r_3} (H/N_eff - M2*M3)^2 / (M2*M3 + eps)  / 144
    noise floor  = #cells / N_eff = 144 / N_eff

Outputs:
    C:\\Collatz\\d_scout_k5_a4_b3.csv     (8x18 deviation matrix)
    C:\\Collatz\\d_scout_k5_a4_b3_log.txt (full report)
"""
from __future__ import annotations
import csv
import os
import sys
import time

import numpy as np
from numba import njit

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz"
SEED = 42
N = 1_000_000
DEPTH_K = 5
A = 4              # mod 2^a = 16
B = 3              # mod 3^b = 27
LOG2_MAX = 60      # safety: |T^5(n)| should stay well under 2^60 from n in [1, 1e12]


# --------------- Syracuse map T(n) = (3n+1)/2^v iterated k times ---------------

@njit(cache=True)
def walk_k_steps(starts, k, n_eff_bool):
    """For each odd start in `starts`, walk k Syracuse 3x+1 compressed steps.
    Write result into `starts` in place. Set n_eff_bool[i] = False if trajectory
    collapsed to 1 (or we overflowed) before depth k.
    """
    n = len(starts)
    for i in range(n):
        m = starts[i]
        collapsed = False
        for _ in range(k):
            if m == 1:
                collapsed = True
                break
            x = 3 * m + 1
            # strip 2's
            while (x & 1) == 0:
                x >>= 1
            m = x
        starts[i] = m
        if collapsed or m == 1:
            n_eff_bool[i] = False
        else:
            n_eff_bool[i] = True


# --------------- Joint histogram + deviation matrix ---------------

def odd_residues_mod_2a(a):
    return [r for r in range(2 ** a) if r % 2 == 1]


def coprime_3_residues_mod_3b(b):
    return [r for r in range(3 ** b) if r % 3 != 0]


def main():
    print("=" * 78)
    print(f"Scout probe: D(r_2, r_3) at k={DEPTH_K}, (a, b)=({A}, {B}), "
          f"N={N:,}")
    print("=" * 78)

    # -- Sample odd starts
    rng = np.random.default_rng(SEED)
    starts = (2 * rng.integers(1, 5 * 10**11, size=N, dtype=np.int64) + 1)
    # values in [3, ~10^12+1]
    print(f"  sampled {N:,} odd starts in [3, ~10^12]")

    # -- Walk k steps
    t0 = time.time()
    n_eff_bool = np.zeros(N, dtype=np.bool_)
    walk_k_steps(starts, DEPTH_K, n_eff_bool)
    elapsed = time.time() - t0
    n_eff = int(n_eff_bool.sum())
    n_collapsed = N - n_eff
    print(f"  walked depth {DEPTH_K} in {elapsed:.2f}s")
    print(f"  N_eff = {n_eff:,}  ({100 * n_eff / N:.4f}% survived)")
    print(f"  collapsed/excluded: {n_collapsed:,}")

    # -- Sanity: all surviving values odd, all coprime to 3
    surv = starts[n_eff_bool]
    n_odd = int(((surv & 1) == 1).sum())
    n_coprime3 = int(((surv % 3) != 0).sum())
    print(f"  sanity: {n_odd:,} of {n_eff:,} surviving values are odd "
          f"({'OK' if n_odd == n_eff else 'FAIL'})")
    print(f"  sanity: {n_coprime3:,} of {n_eff:,} surviving values are "
          f"coprime to 3 ({'OK' if n_coprime3 == n_eff else 'FAIL'})")

    # -- Joint histogram
    M2 = 2 ** A
    M3 = 3 ** B
    r2_arr = surv % M2
    r3_arr = surv % M3
    flat = r2_arr * M3 + r3_arr
    H_full = np.bincount(flat, minlength=M2 * M3).reshape(M2, M3)

    # restrict to odd r_2 (rows where r_2 % 2 == 1) and coprime r_3 (r_3 % 3 != 0)
    odd_r2 = odd_residues_mod_2a(A)        # 8 values
    cop_r3 = coprime_3_residues_mod_3b(B)  # 18 values
    H = H_full[np.ix_(odd_r2, cop_r3)].astype(np.float64)
    H_sum = float(H.sum())
    assert int(H_sum) == n_eff, (
        f"H sum {H_sum} != N_eff {n_eff} (sanity)"
    )
    print(f"  joint histogram: {H.shape} grid, total {int(H_sum):,}")

    # Empirical joint probability and marginals
    P = H / n_eff
    M2_marg = P.sum(axis=1)   # length 8
    M3_marg = P.sum(axis=0)   # length 18
    print(f"  marginals: sum M2 = {M2_marg.sum():.6f}  "
          f"(expect 1), sum M3 = {M3_marg.sum():.6f}  (expect 1)")

    # CRT-product prediction
    CRT = np.outer(M2_marg, M3_marg)
    D = P - CRT
    print(f"  D matrix sum = {D.sum():.4e}  (expect ~0)")

    # chi-squared (per cell normalization), noise floor
    eps = 1e-300
    chi2_cells = (D ** 2) / (CRT + eps)
    chi2 = float(chi2_cells.sum()) / 144
    noise_floor = 144 / n_eff
    chi2_ratio = chi2 / noise_floor
    print()
    print(f"  CHI-SQ (per cell) = {chi2:.6e}")
    print(f"  noise floor       = {noise_floor:.6e}  (= 144 / N_eff)")
    print(f"  ratio (chi2 / noise floor) = {chi2_ratio:.2f}")

    # Dynamic range
    abs_D = np.abs(D)
    max_D = float(abs_D.max())
    p90_D = float(np.percentile(abs_D, 90))
    # per-cell sigma under independence (binomial-style)
    sigma_cell = np.sqrt(CRT * (1 - CRT) / n_eff + eps)
    sigma_mult = abs_D / sigma_cell
    n_above_3sigma = int((sigma_mult > 3).sum())
    print()
    print(f"  max |D|     = {max_D:.6e}")
    print(f"  90% |D|     = {p90_D:.6e}")
    print(f"  cells with |D| > 3 sigma_noise: {n_above_3sigma} of 144  "
          f"({100*n_above_3sigma/144:.1f}%)")

    # Top 5 deviating cells
    flat_idx = np.argsort(sigma_mult.ravel())[::-1][:5]
    print()
    print(f"  Top 5 deviating cells by sigma multiple:")
    print(f"    {'r_2':>4}  {'r_3':>4}  {'D':>14}  {'sigma_mult':>10}")
    for fi in flat_idx:
        i, j = int(fi // 18), int(fi % 18)
        print(f"    {odd_r2[i]:>4}  {cop_r3[j]:>4}  "
              f"{D[i,j]:>+14.6e}  {sigma_mult[i,j]:>10.2f}")

    # Verdict
    print()
    print("=" * 78)
    print("Verdict")
    print("=" * 78)
    if chi2_ratio < 1.5 and max_D < 3 * float(sigma_cell.max()):
        verdict = "CRT-INDEPENDENT (no structure detected at N=10^6)"
        recommendation = "skip full probe, or push N to 10^7-10^8 for tighter test"
    elif chi2_ratio >= 10:
        verdict = "STRUCTURED (joint Bohr deviation detected)"
        recommendation = "fire full multi-panel probe at N=10^7"
    else:
        verdict = "BORDERLINE (1.5 <= chi2/floor < 10)"
        recommendation = "single-panel N=10^7 retest before committing to full"
    print(f"  {verdict}")
    print(f"  recommendation: {recommendation}")
    print()

    # CSV
    out_csv = os.path.join(OUT_DIR, "d_scout_k5_a4_b3.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        # header: r_3 values
        w.writerow(["r_2 \\ r_3"] + cop_r3)
        for i, r2 in enumerate(odd_r2):
            row = [r2] + [f"{D[i, j]:.10e}" for j in range(len(cop_r3))]
            w.writerow(row)
    print(f"  saved D matrix: {out_csv}")

    # also save marginals + chi2 stats
    out_meta = os.path.join(OUT_DIR, "d_scout_k5_a4_b3_meta.csv")
    with open(out_meta, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["N", N])
        w.writerow(["N_eff", n_eff])
        w.writerow(["depth_k", DEPTH_K])
        w.writerow(["a", A])
        w.writerow(["b", B])
        w.writerow(["chi2_per_cell", chi2])
        w.writerow(["noise_floor_chi2", noise_floor])
        w.writerow(["chi2_ratio", chi2_ratio])
        w.writerow(["max_abs_D", max_D])
        w.writerow(["p90_abs_D", p90_D])
        w.writerow(["n_cells_above_3sigma", n_above_3sigma])
        w.writerow(["verdict", verdict])
        w.writerow(["recommendation", recommendation])
        w.writerow([])
        w.writerow(["odd_r2"] + list(odd_r2))
        w.writerow(["M2_marginal"] + [f"{x:.10e}" for x in M2_marg])
        w.writerow(["coprime_r3"] + list(cop_r3))
        w.writerow(["M3_marginal"] + [f"{x:.10e}" for x in M3_marg])
    print(f"  saved metadata: {out_meta}")


if __name__ == "__main__":
    main()
