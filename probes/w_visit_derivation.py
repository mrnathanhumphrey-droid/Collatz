"""
w_visit_derivation.py — Result 78 candidate (R77 follow-up).

Test whether the empirical visit-frequency profile W_visit within each
(r mod 32, b) joint state matches R66's 3-adic Bohr stationary π_4 at
higher modulus (mod 81).

Hypothesis: W_visit's non-uniformity reflects 3-adic Bohr concentration
on m mod 81 (R65/R66 framework). If true: build K_full = K_dynamics
weighted by π_4-derived W_visit_predicted; check Perron(K_full)
recovers D_avg.

Steps:
1. Walk forward orbits, record (r mod 32, b, m mod 81) per Syracuse step
2. Compute empirical W_visit[(r, b)][m mod 81]
3. Compute R66's π_4 stationary on (Z/81Z)*
4. Compare W_visit conditional distribution mod 81 to π_4 across (r, b) cells
5. Build K_full with π_4-derived weights, check Perron → D_pred
"""
import math
import sys
import time
from pathlib import Path
from collections import defaultdict
import numpy as np
from numba import njit, prange
import scipy.sparse as sp

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(r"C:\Collatz")
log_lines = []
def log(s):
    print(s, flush=True)
    log_lines.append(s)

R = 16
B = 64
N_STATES = R * B
LOG_2 = math.log(2.0)
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_capture_extended(starts, max_T):
    """For each Syracuse step, record (r mod 32, b, m mod 81).
    Returns flat arrays per orbit, with sentinel -1 for unused slots."""
    n = len(starts)
    r32_seq = np.full((n, max_T), -1, dtype=np.int8)
    b_seq = np.full((n, max_T), -1, dtype=np.int8)
    m81_seq = np.full((n, max_T), -1, dtype=np.int8)
    count = np.zeros(n, dtype=np.int32)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    for i in prange(n):
        m = np.int64(starts[i])
        sigma = 0
        T = 0
        failed = False
        while (m & 1) == 0 and m > 1:
            m >>= 1
            sigma += 1
        if m == 1:
            sigma_arr[i] = sigma
            count[i] = 0
            continue
        while m != 1 and T < max_T:
            if m > MAX_VAL // 3:
                failed = True
                break
            r32 = m & 31
            # Compute b = floor(log_2 m)
            b = 0
            mm = m
            while mm > 1:
                mm >>= 1
                b += 1
            if b > 63:
                b = 63
            m81 = m % 81
            r32_seq[i, T] = r32
            b_seq[i, T] = b
            m81_seq[i, T] = m81
            x = 3 * m + 1
            v = 0
            while (x & 1) == 0:
                x >>= 1
                v += 1
            sigma += 1 + v
            T += 1
            m = x
        if not failed and m == 1:
            sigma_arr[i] = sigma
            count[i] = T
    return r32_seq, b_seq, m81_seq, count, sigma_arr


def build_markov_chain_3adic(k):
    """R66 chain on (Z/3^k Z)*."""
    N = 3 ** k
    M = 2 * 3 ** (k - 1)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime = [r for r in range(N) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(coprime)}
    n = len(coprime)
    K = np.zeros((n, n))
    Z_v = 1.0 - 2.0 ** (-M)
    for r in coprime:
        for v in range(1, M + 1):
            p = 2.0 ** (-v) / Z_v
            tgt = ((3 * r + 1) * powers_inv2[v - 1]) % N
            K[idx[r], idx[tgt]] += p
    return K, coprime, idx


def stationary_via_iteration(K, n_iter=2000):
    n = K.shape[0]
    pi = np.ones(n) / n
    for _ in range(n_iter):
        pi = pi @ K
        pi = pi / pi.sum()
    return pi


def main():
    log("=" * 78)
    log("W_visit derivation from R66 3-adic Bohr concentration (R77 follow-up)")
    log("=" * 78)

    # ---- Step 1: walk orbits ----
    log2N = 32
    N = 1 << log2N
    n_orbits = 1_500_000
    max_T = 250
    log(f"\nStep 1: walk {n_orbits:,} orbits at N=2^{log2N}, capture (r32, b, m81)")
    rng = np.random.default_rng(20260503)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_orbits, dtype=np.int64) + 1
    t0 = time.time()
    r32_seq, b_seq, m81_seq, count, sigma_arr = walk_capture_extended(starts, max_T)
    log(f"  walk: {time.time()-t0:.1f}s")
    ok = (sigma_arr > 0) & (count > 0)
    n_ok = int(ok.sum())
    log(f"  ok orbits: {n_ok:,}")

    # ---- Step 2: aggregate W_visit_empirical ----
    log("\nStep 2: aggregate W_visit_empirical[(r, b)][m mod 81]")
    # visit_counts[(r, b)] = array of length 81 with visit counts per mod-81 residue
    visit_counts = np.zeros((R, B, 81), dtype=np.int64)
    r32_seq_ok = r32_seq[ok]
    b_seq_ok = b_seq[ok]
    m81_seq_ok = m81_seq[ok]
    count_ok = count[ok]

    # Numba-accelerated aggregation
    @njit
    def aggregate(r32_seq, b_seq, m81_seq, count, visit_counts):
        for i in range(len(count)):
            c = count[i]
            for t in range(c):
                r = r32_seq[i, t]
                b = b_seq[i, t]
                m81 = m81_seq[i, t]
                if r > 0 and b >= 0 and m81 >= 0:
                    r_idx = (r - 1) // 2
                    if 0 <= r_idx < R and 0 <= b < B and 0 <= m81 < 81:
                        visit_counts[r_idx, b, m81] += 1

    aggregate(r32_seq_ok, b_seq_ok, m81_seq_ok, count_ok, visit_counts)
    log(f"  total visits aggregated: {visit_counts.sum():,}")

    # Marginal m mod 81 distribution (across all (r, b))
    marginal_m81 = visit_counts.sum(axis=(0, 1))
    marginal_m81_norm = marginal_m81 / marginal_m81.sum()
    log(f"\n  Marginal m mod 81 visit distribution:")
    log(f"  {'m81':>4}  {'count':>10}  {'p_emp':>9}  {'coprime?':>8}")
    for m81 in range(81):
        if marginal_m81[m81] > 0:
            cop = "yes" if m81 % 3 != 0 else "NO"
            log(f"  {m81:>4}  {marginal_m81[m81]:>10,}  "
                f"{marginal_m81_norm[m81]:>9.6f}  {cop:>8}")

    # ---- Step 3: compute R66 π_4 ----
    log("\n" + "=" * 78)
    log("Step 3: R66 chain stationary π_4 on (Z/81Z)*")
    log("=" * 78)
    K_4, coprime_4, idx_4 = build_markov_chain_3adic(4)
    pi_4 = stationary_via_iteration(K_4)
    log(f"  Chain on (Z/81Z)*: {len(coprime_4)} states, pi_4 sum = {pi_4.sum():.6f}")
    log(f"  pi_4 stats: min={pi_4.min():.6f}, max={pi_4.max():.6f}, "
        f"std={pi_4.std():.6f}")

    # Predicted marginal: pi_4 over coprime, ZERO at multiples of 3
    pi_4_full = np.zeros(81)
    for i, r in enumerate(coprime_4):
        pi_4_full[r] = pi_4[i]

    # ---- Step 4: compare empirical marginal to pi_4 ----
    log("\n" + "=" * 78)
    log("Step 4: compare empirical mod-81 marginal to R66 π_4")
    log("=" * 78)

    # The empirical marginal has BOTH coprime-to-3 and mod-3-zero mass
    # R66 π_4 is supported only on coprime residues
    coprime_emp = marginal_m81_norm[[r for r in range(81) if r % 3 != 0]]
    coprime_emp_renorm = coprime_emp / coprime_emp.sum()
    pi_4_renorm = pi_4 / pi_4.sum()

    pearson_marginal = float(np.corrcoef(coprime_emp_renorm, pi_4_renorm)[0, 1])
    mae_marginal = float(np.abs(coprime_emp_renorm - pi_4_renorm).mean())
    mass_zero_mod3 = sum(marginal_m81_norm[r] for r in range(81) if r % 3 == 0)
    log(f"\n  Mass at m ≡ 0 mod 3 residues: {mass_zero_mod3:.4f} "
        f"({100*mass_zero_mod3:.1f}%)")
    log(f"  (R66 chain assumes 0 here; reality has nonzero from forward orbits)")
    log(f"\n  Comparison on coprime-to-3 residues (54 values):")
    log(f"  Pearson(empirical mod-81 marginal, π_4) = {pearson_marginal:.4f}")
    log(f"  MAE = {mae_marginal:.6f}")

    # Show top-mass residues
    log(f"\n  Top 10 coprime residues by EMPIRICAL mod-81 mass:")
    log(f"  {'r81':>4}  {'p_emp':>9}  {'pi_4':>9}  {'ratio':>8}")
    coprime_idx_to_r = {i: r for i, r in enumerate(coprime_4)}
    sorted_idx = np.argsort(-coprime_emp_renorm)
    for i in sorted_idx[:10]:
        r81 = coprime_idx_to_r[i]
        pe = coprime_emp_renorm[i]
        pp = pi_4_renorm[i]
        ratio = pp / pe if pe > 0 else 0
        log(f"  {r81:>4}  {pe:>9.6f}  {pp:>9.6f}  {ratio:>8.4f}")

    # ---- Step 5: per-(r, b) conditional comparison ----
    log("\n" + "=" * 78)
    log("Step 5: per-(r mod 32, b) conditional distributions mod 81 vs π_4")
    log("=" * 78)
    log("\n  Test: within each (r mod 32, b) cell, is the m mod 81 distribution")
    log("        proportional to π_4 (independent of mod-32 conditioning)?")
    log("        Under CRT independence, conditional ≈ marginal.")

    # For each (r, b) with sufficient data, compute Pearson with π_4
    rhos = []
    cell_data = []
    for r_idx in range(R):
        r = 2 * r_idx + 1
        for b in range(B):
            cell_count = visit_counts[r_idx, b, :].sum()
            if cell_count < 200:
                continue
            cell_dist = visit_counts[r_idx, b, :] / cell_count
            cell_dist_coprime = cell_dist[[m81 for m81 in range(81) if m81 % 3 != 0]]
            s_cd = cell_dist_coprime.sum()
            if s_cd <= 0:
                continue
            cell_dist_coprime_renorm = cell_dist_coprime / s_cd
            rho = np.corrcoef(cell_dist_coprime_renorm, pi_4_renorm)[0, 1]
            if np.isfinite(rho):
                rhos.append(rho)
                cell_data.append((r, b, cell_count, rho, s_cd))

    log(f"\n  {len(rhos)} (r, b) cells with ≥ 200 visits")
    log(f"  Per-cell Pearson(empirical conditional, π_4): "
        f"mean={np.mean(rhos):.4f}, median={np.median(rhos):.4f}, "
        f"min={min(rhos):.4f}")

    # Sample cell printouts
    log(f"\n  Sample cells:")
    log(f"  {'r':>3}  {'b':>3}  {'visits':>9}  {'coprime_frac':>13}  {'ρ vs π_4':>10}")
    for r_target, b_target in [(3, 10), (5, 8), (1, 12), (15, 20),
                               (21, 8), (13, 10)]:
        for r, b, cnt, rho, s_cd in cell_data:
            if r == r_target and b == b_target:
                log(f"  {r:>3}  {b:>3}  {cnt:>9,}  {s_cd:>13.4f}  {rho:>10.4f}")
                break

    # ---- Step 6: build K_full with π_4-derived weights, compute Perron ----
    log("\n" + "=" * 78)
    log("Step 6: K_full = K_dynamics weighted by π_4-derived W_visit_predicted")
    log("=" * 78)

    # For each (r, b) joint state, treat the within-state m mod 81 distribution
    # as π_4 (renormalized to coprime support). Then for each m mod 81 within
    # the cell, compute (v, r', Δb) deterministically (where possible) and
    # aggregate weighted by π_4.

    # Strategy: for each "fine state" (r mod 32, m mod 81) with m mod 32 = r and
    # m mod 81 = q (q coprime to 3), there's a unique m mod 2592 by CRT. From
    # that, compute v, r' (deterministic for r ≠ 21; for r=21 use Geom shifted).
    # Aggregate over m mod 81 with weights π_4(q).

    # But for many m mod 81, the v depends on m mod 64 (or higher), not just
    # m mod 32. So we need to MC sample m with m ≡ r mod 32, m ≡ q mod 81,
    # and random higher bits, to get the v + r' distribution.

    log("\n  MC sampling for joint (m mod 32, m mod 81) → (v, r', Δb)")
    rng2 = np.random.default_rng(20260504)
    n_mc = 5000  # per (r, q) pair

    # Pre-compute MC transition table[(r, q)] = list of (v, r_next, prob)
    transition_table = {}
    for r in range(1, 32, 2):
        for q in range(81):
            if q % 3 == 0:
                continue
            # Find m mod 2592 such that m ≡ r mod 32, m ≡ q mod 81
            # CRT: 32 * inv(32 mod 81) * q + 81 * inv(81 mod 32) * r mod 2592
            inv32_mod81 = pow(32, -1, 81)
            inv81_mod32 = pow(81 % 32, -1, 32)  # 81 mod 32 = 17, inv 17 mod 32
            m_base = (32 * inv32_mod81 * q + 81 * inv81_mod32 * r) % 2592
            # Sample higher bits
            high = rng2.integers(0, 1 << 18, size=n_mc)
            ms = m_base + 2592 * high
            counts = defaultdict(int)
            for m in ms:
                x = 3 * int(m) + 1
                v = 0
                while x & 1 == 0:
                    x >>= 1
                    v += 1
                r_next = x & 31
                counts[(v, r_next)] += 1
            triples = [(v, r_next, c / n_mc)
                       for (v, r_next), c in counts.items()]
            transition_table[(r, q)] = triples

    log(f"  Built transition table for {len(transition_table)} (r, q) pairs")

    # Build K_full[s, s'] = Σ_q π_4(q) · Σ_{v,r'} P(v, r' | r, q) · P(Δb | v)
    # weighted within each (r, b) cell
    K_full = np.zeros((N_STATES, N_STATES), dtype=np.float64)
    for r in range(1, 32, 2):
        r_idx = (r - 1) // 2
        for b in range(B):
            s = r_idx * B + b
            # Sum over q coprime to 3
            for q_idx, q in enumerate(coprime_4):
                pi_q = pi_4[q_idx]
                if pi_q <= 0:
                    continue
                triples = transition_table.get((r, q))
                if triples is None:
                    continue
                for (v, r_next, p_vrn) in triples:
                    for db, p_db in [(-v + 1, 1/3), (-v + 2, 2/3)]:
                        b_next = b + db
                        if not (0 <= b_next < B):
                            continue
                        r_next_idx = (r_next - 1) // 2
                        if not (0 <= r_next_idx < R):
                            continue
                        s_next = r_next_idx * B + b_next
                        K_full[s, s_next] += pi_q * p_vrn * p_db
    # Renormalize rows (so sum to 1)
    rsum = K_full.sum(axis=1)
    nz = rsum > 1e-12
    K_full[nz, :] = K_full[nz, :] / rsum[nz, None]

    log(f"  K_full sparsity: {(K_full > 0).sum()} nonzero / {N_STATES**2}")

    # ---- Compute Perron and compare to D_avg ----
    log("\n  Computing Perron eigvec of K_full")
    pi_full = stationary_via_iteration(K_full, n_iter=2000)
    log(f"  pi_full sum = {pi_full.sum():.6f}")
    rho_pred = np.zeros(R)
    for r_idx in range(R):
        rho_pred[r_idx] = pi_full[r_idx * B:(r_idx + 1) * B].sum()
    if rho_pred.sum() > 0:
        rho_pred = rho_pred / rho_pred.sum()

    # Load D_avg
    D_avg = {}
    pi_32 = {}
    with open(OUT / "qsd_late_t_avg.csv") as f:
        f.readline()
        for line in f:
            parts = line.strip().split(",")
            r = int(parts[0])
            pi_32[r] = float(parts[1])
            D_avg[r] = float(parts[2])
    odd_r32 = list(range(1, 32, 2))

    log(f"\n  D_pred(r) from K_full Perron eigvec:")
    log(f"  {'r':>3}  {'rho_pred':>10}  {'D_pred':>8}  {'D_avg':>8}  {'diff':>9}")
    total_dev = 0.0
    D_pred_arr = []
    D_avg_arr = []
    for i, r in enumerate(odd_r32):
        d_pred = rho_pred[i] / pi_32[r]
        d_emp = D_avg[r]
        diff = d_pred - d_emp
        total_dev += abs(diff)
        D_pred_arr.append(d_pred)
        D_avg_arr.append(d_emp)
        log(f"  {r:>3}  {rho_pred[i]:>10.6f}  {d_pred:>8.4f}  {d_emp:>8.4f}  "
            f"{diff:>+9.4f}")
    pearson = float(np.corrcoef(D_pred_arr, D_avg_arr)[0, 1])
    log(f"\n  Total |D_pred - D_avg| = {total_dev:.4f}")
    log(f"  Pearson rho = {pearson:.4f}")
    log(f"  R60 empirical reference: total_dev = 3.40, ρ = +0.80")
    log(f"  R77 K_dynamics-only reference: total_dev = 5.40, ρ = -0.57")

    # ---- Save outputs ----
    log("\n" + "=" * 78)
    log("Save outputs")
    log("=" * 78)
    sp.save_npz(OUT / "K_full.npz", sp.csr_matrix(K_full))
    log("  [wrote] K_full.npz")

    with open(OUT / "K_full_perron_results.csv", "w") as f:
        f.write("r,pi_32,D_pred,D_avg,diff\n")
        for i, r in enumerate(odd_r32):
            f.write(f"{r},{pi_32[r]:.6f},{D_pred_arr[i]:.6f},"
                    f"{D_avg_arr[i]:.6f},{D_pred_arr[i]-D_avg_arr[i]:+.6f}\n")
    log("  [wrote] K_full_perron_results.csv")

    with open(OUT / "w_visit_predicted_vs_empirical.csv", "w") as f:
        f.write("r,b,n_visits,coprime_frac,pearson_pi4\n")
        for r, b, cnt, rho, s_cd in cell_data:
            f.write(f"{r},{b},{cnt},{s_cd:.6f},{rho:.6f}\n")
    log("  [wrote] w_visit_predicted_vs_empirical.csv")

    # ---- Verdict ----
    log("\n" + "=" * 78)
    log("VERDICT")
    log("=" * 78)
    log(f"\n  Per-cell Pearson(W_visit_emp conditional, π_4): "
        f"mean = {np.mean(rhos):.4f}")
    log(f"  Marginal Pearson(W_visit_emp marginal mod 81, π_4): "
        f"{pearson_marginal:.4f}")
    log(f"  K_full Perron → D_pred: total_dev = {total_dev:.4f}, "
        f"ρ = {pearson:.4f}")
    log(f"")
    if total_dev < 3.5 and pearson > 0.75 and pearson_marginal > 0.85:
        log(f"  Outcome (α): W_visit ≈ π_4-derived; K_full Perron recovers D_avg.")
        log(f"               Trajectory measure ↔ Perron(K_full) is FIRST-PRINCIPLES.")
    elif total_dev < 5.0 and pearson_marginal > 0.6:
        log(f"  Outcome (β): π_4 captures partial structure; residuals identify")
        log(f"               additional non-3-adic mechanism.")
    else:
        log(f"  Outcome (γ): π_4 alone doesn't explain W_visit; need additional")
        log(f"               structure (2-adic or joint 2-3-adic).")

    (OUT / "w_visit_derivation_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )
    log("\n  [wrote] w_visit_derivation_log.txt")


if __name__ == "__main__":
    main()
