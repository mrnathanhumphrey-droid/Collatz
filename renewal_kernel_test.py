"""
renewal_kernel_test.py — Result 53 candidate: build renewal kernel on
visit-event space at r=21 mod 32, marginalize via exit residue, test
whether stationary reproduces D_avg.

Framework (Family B, the renewal kernel candidate):
  - Each visit to r=21 mod 32 has state (m mod 64, V, band, G)
  - Per Result 49: V values are i.i.d. given band (kernel rows ~ identical,
    spectral gap 0.96)
  - So leading eigenvector of K_B is P(V | band) (marginal)
  - At each visit, exit residue m_next mod 32 is partially deterministic from
    (m mod 64, V): m=53 always exits to 5; m=21 exits depend on higher bits
  - D_predicted(r mod 32) = visit-rate-weighted exit residue distribution

Compare D_predicted to empirical D_avg from Result 50.

Two interpretations of D_predicted:
  (a) Exit residue distribution (just the immediate next residue after visit)
  (b) Renewal-reward: time fraction at each residue, computed by
      (visits per cycle) / (cycle length). For the simple framework,
      D_predicted(r) is the per-visit contribution to time-average.
"""
import sys
import time
from pathlib import Path
import numpy as np
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(r"C:\Collatz")
log_lines = []
def log(s):
    print(s, flush=True)
    log_lines.append(s)

LOG_2 = np.log(2.0)
LOG_3 = np.log(3.0)
K_H = 3.0 / np.log(4.0/3.0)
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_capture(starts, max_T, max_visits):
    """
    For each orbit, record at each r=21 mod 32 visit:
      [T_step, V, exit_mod32, m_mod_64, gap, full_residue_path_length]
    Plus: also record the residue path between consecutive visits
    (we'll record total time and residue distribution per band separately).

    Returns:
      visits[i, k, 0..4] = T, V, exit_mod32, m_mod_64, gap
      visit_count[i] = number of visits
      sigma_arr[i] = total sigma (full Collatz steps)
      T_arr[i] = total T (odd-skeleton steps)
      residue_counts[i, r] = #times orbit visits residue r mod 32 (all r)
    """
    n = len(starts)
    visits = np.full((n, max_visits, 5), -1, dtype=np.int32)
    visit_count = np.zeros(n, dtype=np.int32)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    T_arr = np.zeros(n, dtype=np.int32)
    residue_counts = np.zeros((n, 32), dtype=np.int64)

    for i in prange(n):
        m = np.int64(starts[i])
        sigma = 0
        T = 0
        last_visit_T = -1
        vc = 0
        failed = False
        while m != 1 and T < max_T:
            if (m & 1) == 0:
                m = m >> 1
                sigma += 1
                continue
            if m > MAX_VAL // 3:
                failed = True
                break
            r32 = m & 31
            residue_counts[i, r32] += 1
            x = 3 * m + 1
            v = 0
            while (x & 1) == 0:
                x >>= 1
                v += 1

            if r32 == 21 and vc < max_visits:
                gap = T - last_visit_T if last_visit_T >= 0 else -1
                visits[i, vc, 0] = T
                visits[i, vc, 1] = v
                visits[i, vc, 2] = x & 31    # m_next mod 32 (exit)
                visits[i, vc, 3] = m & 63    # m mod 64 at visit
                visits[i, vc, 4] = gap
                last_visit_T = T
                vc += 1

            sigma += 1 + v
            T += 1
            m = x

        if not failed and m == 1:
            sigma_arr[i] = sigma
            T_arr[i] = T
            visit_count[i] = vc
    return visits, visit_count, sigma_arr, T_arr, residue_counts


def chang_pi_32():
    from fractions import Fraction
    odd_residues = list(range(1, 64, 2))
    idx = {r: i for i, r in enumerate(odd_residues)}
    P = [[Fraction(0)] * 32 for _ in range(32)]
    for i, r in enumerate(odd_residues):
        counts = [0] * 32
        for k in range(128):
            n = r + 64 * k
            mm = 3 * n + 1
            while mm & 1 == 0:
                mm >>= 1
            counts[idx[mm % 64]] += 1
        for j in range(32):
            P[i][j] = Fraction(counts[j], 128)
    n = 32
    A = [[P[j][i] - (Fraction(1) if i == j else Fraction(0))
          for j in range(n)] for i in range(n)]
    A[n - 1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n - 1] = Fraction(1)
    for col in range(n):
        pivot = next(row for row in range(col, n) if A[row][col] != 0)
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            b[col], b[pivot] = b[pivot], b[col]
        piv = A[col][col]
        for j in range(col, n):
            A[col][j] /= piv
        b[col] /= piv
        for row in range(n):
            if row != col and A[row][col] != 0:
                factor = A[row][col]
                for j in range(col, n):
                    A[row][j] -= factor * A[col][j]
                b[row] -= factor * b[col]
    pi_32 = {r32: float(b[idx[r32]] + b[idx[r32 + 32]])
             for r32 in range(1, 32, 2)}
    return pi_32


def load_empirical_D_avg():
    D_avg = {}
    with open(OUT / "qsd_late_t_avg.csv") as f:
        f.readline()
        for line in f:
            parts = line.strip().split(",")
            D_avg[int(parts[0])] = float(parts[2])
    return D_avg


def main():
    log("=" * 78)
    log("Renewal kernel on visit-event space — D_predicted vs D_avg")
    log("=" * 78)

    pi_32 = chang_pi_32()
    D_avg = load_empirical_D_avg()
    odd_r32 = list(range(1, 32, 2))

    log("\nReference: Chang's pi_32 (essentially uniform 0.0625)")
    log("Reference: D_avg from Result 50 (late-t empirical)")
    for r in odd_r32:
        log(f"  r={r:>2}: pi={pi_32[r]:.6f}, D_avg={D_avg[r]:.4f}")

    # ------------------------------------------------------------------
    # Walk orbits with visit capture
    # ------------------------------------------------------------------
    log2N = 32
    N = 1 << log2N
    n_orbits_per_seed = 500_000
    seeds = [42, 137, 271]  # 1.5M orbits
    log(f"\nWalking {len(seeds) * n_orbits_per_seed:,} orbits at N=2^{log2N}")

    all_visits = []
    all_vc = []
    all_sigma = []
    all_logn = []
    all_residue_counts = []

    t0 = time.time()
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1) // 2,
                                  size=n_orbits_per_seed,
                                  dtype=np.int64) + 1
        visits, vc, sigma, T, rc = walk_capture(starts, 500, 60)
        ok = sigma > 0
        all_visits.append(visits[ok])
        all_vc.append(vc[ok])
        all_sigma.append(sigma[ok].astype(np.float64))
        all_logn.append(np.log(starts[ok].astype(np.float64)))
        all_residue_counts.append(rc[ok])

    visits = np.concatenate(all_visits, axis=0)
    vc = np.concatenate(all_vc)
    sigma = np.concatenate(all_sigma)
    log_n = np.concatenate(all_logn)
    rc_all = np.concatenate(all_residue_counts, axis=0)
    n_orbits = len(sigma)
    log(f"  walk: {time.time() - t0:.1f}s, n_orbits ok = {n_orbits:,}")

    # σ-bands
    sigma_resid = sigma - K_H * log_n
    edges = np.percentile(sigma_resid, [25, 50, 75, 95])
    band = np.digitize(sigma_resid, edges)
    band_names = ['0-25', '25-50', '50-75', '75-95', '95-100']

    # ------------------------------------------------------------------
    # Build flat visit table
    # ------------------------------------------------------------------
    flat_orbit = []
    flat_V = []
    flat_exit = []
    flat_mod64 = []
    flat_gap = []
    flat_band = []
    for i in range(n_orbits):
        nv = vc[i]
        for k in range(nv):
            flat_orbit.append(i)
            flat_V.append(visits[i, k, 1])
            flat_exit.append(visits[i, k, 2])
            flat_mod64.append(visits[i, k, 3])
            flat_gap.append(visits[i, k, 4])
            flat_band.append(band[i])
    flat_V = np.array(flat_V, dtype=np.int8)
    flat_exit = np.array(flat_exit, dtype=np.int8)
    flat_mod64 = np.array(flat_mod64, dtype=np.int8)
    flat_gap = np.array(flat_gap, dtype=np.int32)
    flat_band = np.array(flat_band, dtype=np.int8)

    n_visits = len(flat_V)
    log(f"  total r=21 visits: {n_visits:,}")

    # ------------------------------------------------------------------
    # Marginal P(V | band) and P(m mod 64 | band) — these are the eigenvectors
    # under Result 49's i.i.d. result
    # ------------------------------------------------------------------
    log("\n" + "=" * 78)
    log("Step 1: P(V | band) — leading eigenvector of K_B (Result 49)")
    log("=" * 78)
    log(f"  {'band':>6}  {'n_vis':>8}  " +
        "  ".join(f"V={v:>2}" for v in range(5, 13)))
    for b in range(5):
        mask = flat_band == b
        n_b = int(mask.sum())
        if n_b == 0:
            continue
        V_b = flat_V[mask]
        line = f"  {band_names[b]:>6}  {n_b:>8}"
        for v in range(5, 13):
            p = (V_b == v).mean()
            line += f"  {p:.3f}"
        log(line)

    log("\nP(m mod 64 = 21 | band):")
    for b in range(5):
        mask = flat_band == b
        if mask.sum() == 0:
            continue
        p21 = (flat_mod64[mask] == 21).mean()
        log(f"  {band_names[b]:>6}: p_21 = {p21:.4f}, p_53 = {1-p21:.4f}")

    # ------------------------------------------------------------------
    # Empirical exit-residue distribution conditional on (V, band)
    # ------------------------------------------------------------------
    log("\n" + "=" * 78)
    log("Step 2: P(exit_mod32 | V, band) — exit-residue distribution")
    log("=" * 78)

    # exit_distribution[band, V, residue]
    exit_dist = np.zeros((5, 13, 32))  # 5 bands, V in 5..12+, 32 residues
    for v_idx in range(n_visits):
        b = flat_band[v_idx]
        V = min(int(flat_V[v_idx]), 12)
        e = int(flat_exit[v_idx])
        exit_dist[b, V, e] += 1

    # Normalize
    for b in range(5):
        for V in range(13):
            tot = exit_dist[b, V, :].sum()
            if tot > 0:
                exit_dist[b, V, :] /= tot

    # ------------------------------------------------------------------
    # D_predicted(r) via renewal-kernel framework
    # Two computations:
    #   (a) Exit-residue marginal: D_pred_exit(r) = Σ_b P(b) Σ_V P(V|b) exit_dist[b,V,r]
    #   (b) Time-average residue: D_pred_time(r) = empirical residue counts / total time
    # ------------------------------------------------------------------
    log("\n" + "=" * 78)
    log("Step 3: D_predicted via renewal-kernel framework (exit marginal)")
    log("=" * 78)

    # Per-band weight = fraction of visits in band (= renewal-rate weighting)
    band_visit_count = np.array([
        (flat_band == b).sum() for b in range(5)
    ], dtype=np.float64)
    band_weight = band_visit_count / band_visit_count.sum()
    log(f"\n  Band weights (visit-frequency): " +
        ", ".join(f"{band_names[b]}={band_weight[b]:.3f}"
                  for b in range(5)))

    # P(V | band)
    P_V_given_b = np.zeros((5, 13))
    for b in range(5):
        mask = flat_band == b
        if mask.sum() == 0:
            continue
        V_b = flat_V[mask]
        for V in range(5, 13):
            P_V_given_b[b, V] = (V_b == min(V, 12)).mean()
        # Roll up V >= 12 into V=12
        P_V_given_b[b, 12] = (V_b >= 12).mean()
        # Recompute below-12 to actual count
        for V in range(5, 12):
            P_V_given_b[b, V] = (V_b == V).mean()

    # Exit-marginal D_predicted
    D_pred_exit = np.zeros(32)
    for b in range(5):
        for V in range(5, 13):
            D_pred_exit += band_weight[b] * P_V_given_b[b, V] * exit_dist[b, V, :]

    log("\n  D_predicted_exit(r mod 32) = visit-weighted exit residue distribution")
    log(f"\n  {'r':>3}  {'D_avg':>7}  {'D_exit':>7}  {'diff':>8}")
    total_dev_exit = 0
    for r in odd_r32:
        rho_exit = D_pred_exit[r]
        D_pred_r = rho_exit / pi_32[r]
        diff = D_pred_r - D_avg[r]
        total_dev_exit += abs(diff)
        log(f"  {r:>3}  {D_avg[r]:>7.4f}  {D_pred_r:>7.4f}  {diff:>+8.4f}")
    log(f"\n  Total |D_pred_exit - D_avg| = {total_dev_exit:.4f}")

    # ------------------------------------------------------------------
    # Time-average residue (proper trajectory measure from this orbit set)
    # ------------------------------------------------------------------
    log("\n" + "=" * 78)
    log("Step 4: D_traj_full = ALL residues from THIS orbit set (sanity check)")
    log("=" * 78)

    rho_full = rc_all.sum(axis=0)
    rho_full_odd = np.zeros(32)
    for r in odd_r32:
        rho_full_odd[r] = rho_full[r] / rho_full[1::2].sum()

    log(f"\n  Time-average residue (full orbit, including all m's, this dataset)")
    log(f"\n  {'r':>3}  {'D_avg(R50)':>11}  {'D_full':>7}  {'diff':>8}")
    total_dev_full = 0
    for r in odd_r32:
        D_full_r = rho_full_odd[r] / pi_32[r]
        diff = D_full_r - D_avg[r]
        total_dev_full += abs(diff)
        log(f"  {r:>3}  {D_avg[r]:>11.4f}  {D_full_r:>7.4f}  {diff:>+8.4f}")
    log(f"\n  Total |D_full - D_avg(R50)| = {total_dev_full:.4f}")
    log(f"  (D_full averaged over ENTIRE orbit, including transient — should differ"
        f"\n   from D_avg which averages over LATE-t survivors)")

    # ------------------------------------------------------------------
    # Late-t residue from THIS dataset (for an apples-to-apples reference)
    # ------------------------------------------------------------------
    # Need to walk and record m mod 32 at each step to compute late-t avg
    # Skipping for time; the comparison D_pred_exit vs D_avg is the key test
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Per-band D_pred_exit for diagnostic
    # ------------------------------------------------------------------
    log("\n" + "=" * 78)
    log("Step 5: D_pred per band (diagnostic)")
    log("=" * 78)

    line = f"\n  {'r':>3}  {'D_avg':>6}"
    for b in range(5):
        line += f"  {'b='+band_names[b]:>10}"
    log(line)

    D_pred_per_band = np.zeros((5, 32))
    for b in range(5):
        for V in range(5, 13):
            D_pred_per_band[b, :] += P_V_given_b[b, V] * exit_dist[b, V, :]

    for r in odd_r32:
        line = f"  {r:>3}  {D_avg[r]:>6.3f}"
        for b in range(5):
            d_pred = D_pred_per_band[b, r] / pi_32[r]
            line += f"  {d_pred:>10.4f}"
        log(line)

    # ------------------------------------------------------------------
    # Pearson correlation
    # ------------------------------------------------------------------
    D_avg_arr = np.array([D_avg[r] for r in odd_r32])
    D_pred_arr = np.array([D_pred_exit[r] / pi_32[r] for r in odd_r32])
    rho = np.corrcoef(D_avg_arr, D_pred_arr)[0, 1]

    log("\n" + "=" * 78)
    log("Step 6: Pearson correlation D_predicted vs D_avg")
    log("=" * 78)
    log(f"\n  Pearson ρ = {rho:.4f}")
    log(f"  Total |D_pred - D_avg| (16 residues) = {total_dev_exit:.4f}")
    log(f"  Mean |diff| per residue = {total_dev_exit / 16:.4f}")

    # ------------------------------------------------------------------
    # Save CSVs
    # ------------------------------------------------------------------
    with open(OUT / "renewal_kernel_predicted_D.csv", "w") as f:
        f.write("r,pi,D_avg,D_pred_exit,D_full,diff_exit\n")
        for r in odd_r32:
            d_pred = D_pred_exit[r] / pi_32[r]
            d_full = rho_full_odd[r] / pi_32[r]
            f.write(f"{r},{pi_32[r]:.6f},{D_avg[r]:.6f},{d_pred:.6f},"
                    f"{d_full:.6f},{d_pred - D_avg[r]:+.6f}\n")

    with open(OUT / "renewal_kernel_eigvecs.csv", "w") as f:
        f.write("band,V,P_V_given_band,p_21\n")
        for b in range(5):
            mask = flat_band == b
            if mask.sum() == 0:
                continue
            p21 = (flat_mod64[mask] == 21).mean()
            for V in range(5, 13):
                f.write(f"{band_names[b]},{V},{P_V_given_b[b, V]:.6f},{p21:.6f}\n")

    (OUT / "renewal_kernel_test_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )
    log("\n[wrote] renewal_kernel_predicted_D.csv, renewal_kernel_eigvecs.csv,"
        " renewal_kernel_test_log.txt")


if __name__ == "__main__":
    main()
