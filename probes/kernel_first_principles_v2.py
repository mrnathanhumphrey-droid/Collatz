"""
kernel_first_principles_v2.py — corrected derivation.

Bug in v1: assumed r' depends only on (r mod 32, v). Actually for v=1 (r mod 4
= 3), the next residue r' alternates between two values by m mod 64. Similarly
for v=2 (r mod 8 = 1), r' takes 4 values by m mod 128. Fix: for each r mod 32,
MC sample m with high bits and build the proper conditional joint
distribution of (v, r') given r mod 32.

For the b transition: under within-state-uniform-m assumption, log_2(m) is
approx uniform in [b, b+1), so the b' = floor(log_2(m')) split is governed by
the fractional part of log_2(3m+1) - v. This is the "first-principles" piece.
"""
import math
import sys
from pathlib import Path
from collections import defaultdict
import numpy as np
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


def state_idx(r, b):
    return ((r - 1) // 2) * B + b


def derive_v_residue_table_mc(n_samples=500_000, max_high_bits=25):
    """For each r mod 32, MC sample m ≡ r mod 32 with random higher bits,
    walk one Syracuse step, build joint conditional P(v, r' | r mod 32)."""
    rng = np.random.default_rng(20260503)
    table = {}
    for r in range(1, 32, 2):
        high_bits = rng.integers(0, 1 << max_high_bits, size=n_samples)
        ms = r + 32 * high_bits.astype(np.int64)
        counts = defaultdict(int)
        for m in ms:
            x = 3 * int(m) + 1
            v = 0
            while x & 1 == 0:
                x >>= 1
                v += 1
            r_next = x & 31
            counts[(v, r_next)] += 1
        triples = []
        for (v, r_next), c in sorted(counts.items()):
            prob = c / n_samples
            triples.append((v, r_next, prob))
        table[r] = triples
    return table


def derive_kernel(table):
    """Build derived K from the per-residue (v, r') joint table.
    For each (r, b) state and (v, r') outcome, the Δb distribution comes from
    the within-state uniform-m assumption: Δb = -v + 1 with prob 1/3, -v + 2
    with prob 2/3 (since u = m/2^b is uniform in [1, 2), and the breakpoint is
    u = 4/3)."""
    K = np.zeros((N_STATES, N_STATES), dtype=np.float64)
    for r in range(1, 32, 2):
        for b in range(B):
            s = state_idx(r, b)
            for (v, r_next, p_v_rn) in table[r]:
                for db, p_db in [(-v + 1, 1/3), (-v + 2, 2/3)]:
                    b_next = b + db
                    if not (0 <= b_next < B):
                        continue
                    s_next = state_idx(r_next, b_next)
                    K[s, s_next] += p_v_rn * p_db
    return K


def stationary_via_iteration(K, n_iter=2000):
    n = K.shape[0]
    pi = np.ones(n) / n
    for _ in range(n_iter):
        pi = pi @ K
        s = pi.sum()
        if s > 0:
            pi = pi / s
    return pi


def main():
    log("=" * 78)
    log("Derive K from first principles — v2 (proper joint (v, r') conditional)")
    log("=" * 78)

    log("\nStep 1: build conditional P(v, r' | r mod 32) via MC")
    table = derive_v_residue_table_mc()
    log(f"\n  For each r, top (v, r', prob) entries:")
    for r in sorted(table.keys()):
        top = sorted(table[r], key=lambda t: -t[2])[:5]
        line = f"  r={r:>2}:  "
        for v, r_next, p in top:
            line += f"v={v}→r'={r_next}({p:.3f})  "
        log(line)

    log("\n  Number of distinct (v, r') outcomes per r:")
    log(f"  {'r':>3}  {'#outcomes':>10}  {'sum probs':>10}")
    for r in sorted(table.keys()):
        triples = table[r]
        sum_p = sum(t[2] for t in triples)
        log(f"  {r:>3}  {len(triples):>10}  {sum_p:>10.6f}")

    log("\nStep 2: build derived K_derived (1024×1024)")
    K_derived = derive_kernel(table)
    log(f"  Sparsity: {(K_derived > 0).sum()} nonzero / {N_STATES**2} = "
        f"{100*(K_derived>0).sum()/N_STATES**2:.2f}%")
    rsum = K_derived.sum(axis=1)
    log(f"  Row sums: min={rsum.min():.4f}, max={rsum.max():.4f}, "
        f"mean={rsum.mean():.4f}")

    # ---- Compare to R60 empirical ----
    log("\n" + "=" * 78)
    log("Step 3: compare K_derived (v2) to R60 empirical kernel")
    log("=" * 78)
    K_emp_sparse = sp.load_npz(OUT / "size_stratified_kernel.npz")
    keep_idx = np.load(OUT / "size_stratified_keep_idx.npy")
    K_emp_full = np.zeros((N_STATES, N_STATES), dtype=np.float64)
    K_emp_dense = K_emp_sparse.toarray()
    for i, s in enumerate(keep_idx):
        for j, s_next in enumerate(keep_idx):
            K_emp_full[s, s_next] = K_emp_dense[i, j]

    log(f"\n  Sample state-transition checks (top 5 transitions):")
    for r, b in [(3, 10), (5, 8), (1, 12), (15, 20), (21, 8)]:
        s = state_idx(r, b)
        d_row = K_derived[s, :]
        e_row = K_emp_full[s, :]
        if d_row.sum() == 0 and e_row.sum() == 0:
            continue
        top_d = np.argsort(-d_row)[:5]
        top_e = np.argsort(-e_row)[:5]
        log(f"\n    State (r={r}, b={b}):")
        for tag, top, row in [("DERIVED", top_d, d_row), ("EMPIRIC", top_e, e_row)]:
            top_str = ", ".join(
                f"({(idx//B)*2+1}|{idx%B}):{row[idx]:.4f}"
                for idx in top if row[idx] > 1e-4
            )
            log(f"      {tag}: {top_str}")

    # Per-state row Pearson
    rhos = []
    common_states = []
    for s in keep_idx:
        d_row = K_derived[s, :]
        e_row = K_emp_full[s, :]
        if d_row.sum() > 0.5 and e_row.sum() > 0.5:
            rho = np.corrcoef(d_row, e_row)[0, 1]
            if np.isfinite(rho):
                rhos.append(rho)
                common_states.append(s)
    log(f"\n  Per-state row Pearson rho (K_derived vs K_emp):")
    log(f"    n_common = {len(rhos)}, mean = {np.mean(rhos):.4f}, "
        f"median = {np.median(rhos):.4f}, min = {min(rhos):.4f}")

    # MAE on transitions
    diff = K_derived - K_emp_full
    nonzero_mask = (K_derived > 1e-9) | (K_emp_full > 1e-9)
    n_nonzero = int(nonzero_mask.sum())
    mae = np.abs(diff[nonzero_mask]).mean() if n_nonzero > 0 else 0
    log(f"  MAE over transitions with mass: {mae:.6f}")
    log(f"  Total |K_derived - K_emp|: {np.abs(diff).sum():.4f}")

    # ---- Perron eigvec ----
    log("\n" + "=" * 78)
    log("Step 4: Perron eigvec of K_derived → D_pred(r)")
    log("=" * 78)
    pi_derived = stationary_via_iteration(K_derived)
    log(f"  pi_derived sum = {pi_derived.sum():.6f}")
    rho_pred = np.zeros(R)
    for r_idx in range(R):
        rho_pred[r_idx] = pi_derived[r_idx * B:(r_idx + 1) * B].sum()

    # Load empirical D_avg + Chang pi_32
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

    log(f"\n  {'r':>3}  {'rho_pred':>10}  {'D_pred':>8}  {'D_avg':>8}  {'diff':>9}")
    total_dev = 0.0
    D_pred_arr = []
    D_avg_arr = []
    for i, r in enumerate(odd_r32):
        rho = rho_pred[i] / rho_pred.sum() if rho_pred.sum() > 0 else 0
        d_pred = rho / pi_32[r]
        d_emp = D_avg[r]
        diff = d_pred - d_emp
        total_dev += abs(diff)
        D_pred_arr.append(d_pred)
        D_avg_arr.append(d_emp)
        log(f"  {r:>3}  {rho:>10.6f}  {d_pred:>8.4f}  {d_emp:>8.4f}  "
            f"{diff:>+9.4f}")
    pearson = float(np.corrcoef(D_pred_arr, D_avg_arr)[0, 1])
    log(f"\n  Total |D_pred - D_avg| = {total_dev:.4f}")
    log(f"  Pearson rho = {pearson:.4f}")
    log(f"")
    log(f"  Reference points:")
    log(f"    R60 empirical (size-stratified, 1.5M orbits):  total_dev=3.40, ρ=0.80")
    log(f"    Trivial null (entire-orbit time average):       total_dev=4.72")

    # ---- Save outputs ----
    log("\n" + "=" * 78)
    log("Save outputs")
    log("=" * 78)
    sp.save_npz(OUT / "K_derived_v2.npz", sp.csr_matrix(K_derived))
    log("  [wrote] K_derived_v2.npz")

    with open(OUT / "D_predicted_derived_v2.csv", "w") as f:
        f.write("r,pi_32,D_pred_derived,D_avg,diff\n")
        for i, r in enumerate(odd_r32):
            f.write(f"{r},{pi_32[r]:.6f},{D_pred_arr[i]:.6f},"
                    f"{D_avg_arr[i]:.6f},{D_pred_arr[i]-D_avg_arr[i]:+.6f}\n")
    log("  [wrote] D_predicted_derived_v2.csv")

    # K_derived_vs_empirical comparison CSV
    with open(OUT / "K_derived_vs_empirical_v2.csv", "w") as f:
        f.write("from_state,to_state,K_derived,K_empirical,diff\n")
        for s in keep_idx:
            for s_next in keep_idx:
                d = K_derived[s, s_next]
                e = K_emp_full[s, s_next]
                if d > 1e-6 or e > 1e-6:
                    r_from = (s // B) * 2 + 1
                    b_from = s % B
                    r_to = (s_next // B) * 2 + 1
                    b_to = s_next % B
                    f.write(f"({r_from}|{b_from}),({r_to}|{b_to}),"
                            f"{d:.6f},{e:.6f},{d-e:+.6f}\n")
    log("  [wrote] K_derived_vs_empirical_v2.csv")

    # ---- Verdict ----
    log("\n" + "=" * 78)
    log("VERDICT")
    log("=" * 78)
    log(f"\n  K_derived v2 vs D_avg: total_dev = {total_dev:.4f}, ρ = {pearson:.4f}")
    log(f"  Per-state row Pearson (K_derived vs K_emp): mean = {np.mean(rhos):.4f}")

    if total_dev < 3.5 and pearson > 0.75 and np.mean(rhos) > 0.85:
        log(f"  Outcome (α): K_derived ≈ K_empirical structurally")
    elif total_dev < 5.0 and pearson > 0.6 and np.mean(rhos) > 0.7:
        log(f"  Outcome (β): K_derived captures core structure; some residuals")
    else:
        log(f"  Outcome (γ): K_derived differs systematically from K_empirical")
        log(f"               → joint state Markov assumption is missing structure")

    (OUT / "kernel_first_principles_v2_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )
    log("\n  [wrote] kernel_first_principles_v2_log.txt")


if __name__ == "__main__":
    main()
