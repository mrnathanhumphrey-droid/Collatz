"""
kernel_first_principles.py — Result 77 candidate.

Derive R60's 1024-state size-stratified Markov kernel K from Syracuse first
principles. Compare to R60's empirical kernel.

State: s = (r mod 32, b) with r ∈ {1, 3, ..., 31} (16 odd residues) and
b = floor(log_2 m) (64 size bins). 16 × 64 = 1024 states.

Construction (per the brief):
  Step 1: residue-dim transitions from arithmetic-deterministic v rule (R68)
  Step 2: size-dim transitions from log_2 random walk with v controlling drift
  Step 3: combine into joint kernel
  Step 4: compare to R60 empirical (size_stratified_kernel.npz)
  Step 5: Perron eigenvec → D_pred

Conditional v distribution given r mod 32 (from arithmetic determinism):
  r ∈ {3,7,11,15,19,23,27,31}: v=1 deterministic         (8 residues)
  r ∈ {1,9,17,25}:             v=2 deterministic          (4 residues)
  r ∈ {13,29}:                 v=3 deterministic          (2 residues)
  r = 5:                       v=4 deterministic          (1 residue)
  r = 21:                      v ≥ 5, geometric tail      (1 residue, the
                                                          singular boundary)

Conditional Δb given (v, b) under uniform-m-within-state:
  Δb = -v + 1 with prob 1/3 (when m fractional part u < 4/3 in u ∈ [1, 2))
  Δb = -v + 2 with prob 2/3 (when u ≥ 4/3)
"""
import math
import sys
import time
from pathlib import Path
from collections import defaultdict
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigs

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(r"C:\Collatz")
log_lines = []
def log(s):
    print(s, flush=True)
    log_lines.append(s)

R = 16     # 16 odd residues mod 32
B = 64     # log_2 size bins
N_STATES = R * B
LOG_2 = math.log(2.0)


def state_idx(r, b):
    """r is the actual odd residue mod 32 (1, 3, 5, ..., 31); b in 0..63."""
    return ((r - 1) // 2) * B + b


def derive_v_residue_table():
    """
    For each odd r ∈ {1, 3, ..., 31}, return list of (v, r_next, prob)
    triples for the residue-dim transitions, conditional on r alone.
    For r = 21 (singular), the (v, r_next) pairs come from Monte-Carlo
    sampling m ≡ 21 mod 32 with high-precision higher bits.
    """
    table = {}
    for r in range(1, 32, 2):
        if r == 21:
            continue
        # Deterministic v from r:
        if r % 4 == 3:
            v = 1
        elif r % 8 == 1:
            v = 2
        elif r % 16 == 5:  # r=5, 13, 21, 29 mod 16 — but only 5 and 13 here
            # Actually r mod 16 = 5: m mod 16 = 5 gives v ≥ 4; need more bits
            # For r mod 32 = 5 alone:
            if r == 5:
                v = 4
            else:  # r mod 16 = 5 but r ≠ 5 in mod 32 → only r = 21 already handled
                v = None
        elif r % 16 == 13:
            v = 3
        else:
            v = None
        if v is None:
            raise ValueError(f"residue r={r} not classified")
        # Compute r_next deterministically: r_next = (3r + 1) >> v mod 32
        r_next_int = (3 * r + 1) >> v
        r_next = r_next_int & 31  # mod 32
        # Verify r_next is odd
        if r_next % 2 == 0:
            raise ValueError(f"r={r}, v={v} → r_next={r_next} even")
        table[r] = [(v, r_next, 1.0)]

    # r = 21: Monte-Carlo sample m ≡ 21 mod 32 with various higher bits
    rng = np.random.default_rng(20260503)
    n_samples = 200_000
    high_bits = rng.integers(0, 1 << 25, size=n_samples)
    ms = 21 + 32 * high_bits  # m ≡ 21 mod 32
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
    table[21] = triples
    return table


def derive_kernel():
    """
    Build the derived kernel K of shape (N_STATES, N_STATES).
    K[s, s'] = P(s_{t+1} = s' | s_t = s)

    For each (r, b):
      - get conditional v distribution P(v | r) from arithmetic determinism
      - r_next determined by (r, v); for r = 21 from MC table
      - b_next ∈ {b - v + 1, b - v + 2} with probs 1/3, 2/3
    """
    table = derive_v_residue_table()
    K = np.zeros((N_STATES, N_STATES), dtype=np.float64)
    for r in range(1, 32, 2):
        for b in range(B):
            s = state_idx(r, b)
            for (v, r_next, p_v) in table[r]:
                # Δb distribution
                for db, p_db in [(-v + 1, 1/3), (-v + 2, 2/3)]:
                    b_next = b + db
                    if not (0 <= b_next < B):
                        continue
                    s_next = state_idx(r_next, b_next)
                    K[s, s_next] += p_v * p_db
    return K


def stationary_via_iteration(K, n_iter=1000):
    """Power iteration on K^T to get left eigenvector of K with eigenvalue 1."""
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
    log("Derive R60's size-stratified kernel K from Syracuse first principles")
    log("=" * 78)

    # ------ Step 1+2+3: derive K ------
    log("\nStep 1+2+3: derive K from arithmetic-deterministic v rule + log_2 random walk")
    table = derive_v_residue_table()
    log(f"\n  Conditional v distribution per r (from R45/R68 arithmetic determinism):")
    log(f"  {'r':>3}  v distribution → r_next")
    for r in sorted(table.keys()):
        triples = table[r]
        if len(triples) == 1:
            v, r_next, _ = triples[0]
            log(f"  {r:>3}  v={v} (det) → r'={r_next}")
        else:
            # r=21 case: show top entries
            top = sorted(triples, key=lambda t: -t[2])[:6]
            line = f"  {r:>3}  "
            for v, r_next, p in top:
                line += f"v={v}→r'={r_next}({p:.4f})  "
            log(line)
    log(f"\n  r=21 detail (singular boundary):")
    triples_21 = table[21]
    log(f"  {'v':>3}  {'r_next':>6}  {'prob':>10}  {'2^-(v-4)':>10}")
    for v, r_next, p in sorted(triples_21):
        target = 2.0 ** (-(v - 4)) if v >= 5 else 0
        log(f"  {v:>3}  {r_next:>6}  {p:>10.6f}  {target:>10.6f}")

    log(f"\n  Building 1024×1024 derived kernel K_derived...")
    K_derived = derive_kernel()
    log(f"  Sparsity: {(K_derived > 0).sum()} nonzero / {N_STATES**2} = "
        f"{100*(K_derived>0).sum()/N_STATES**2:.2f}%")

    # ------ Step 4: compare to R60 empirical ------
    log("\n" + "=" * 78)
    log("Step 4: compare K_derived to R60 empirical kernel")
    log("=" * 78)
    K_emp_path = OUT / "size_stratified_kernel.npz"
    keep_idx_path = OUT / "size_stratified_keep_idx.npy"
    if K_emp_path.exists() and keep_idx_path.exists():
        K_emp_sparse = sp.load_npz(K_emp_path)
        keep_idx = np.load(keep_idx_path)
        log(f"  Loaded R60 kernel: {K_emp_sparse.shape}, "
            f"{K_emp_sparse.nnz} nonzeros; keep_idx size {len(keep_idx)}")
        # K_emp is on the kept-states subspace only; expand back
        K_emp_full = np.zeros((N_STATES, N_STATES), dtype=np.float64)
        K_emp_dense = K_emp_sparse.toarray()
        for i, s in enumerate(keep_idx):
            for j, s_next in enumerate(keep_idx):
                K_emp_full[s, s_next] = K_emp_dense[i, j]

        # Compare row by row on kept states
        # MAE on transitions where either has mass
        diff = K_derived - K_emp_full
        nonzero_mask = (K_derived > 1e-9) | (K_emp_full > 1e-9)
        n_nonzero = int(nonzero_mask.sum())
        mae = np.abs(diff[nonzero_mask]).mean() if n_nonzero > 0 else 0
        log(f"\n  Transitions with mass in either: {n_nonzero}")
        log(f"  MAE over those entries: {mae:.6f}")
        log(f"  Fraction of derived mass also in empirical: "
            f"{((K_derived > 1e-9) & (K_emp_full > 1e-9)).sum() / max(1, (K_derived > 1e-9).sum()):.4f}")
        log(f"  Fraction of empirical mass also in derived: "
            f"{((K_derived > 1e-9) & (K_emp_full > 1e-9)).sum() / max(1, (K_emp_full > 1e-9).sum()):.4f}")

        # Per-state correlation
        rows_with_data = []
        for s in keep_idx:
            row_d = K_derived[s, :]
            row_e = K_emp_full[s, :]
            if row_d.sum() > 0 and row_e.sum() > 0:
                rho = np.corrcoef(row_d, row_e)[0, 1]
                rows_with_data.append(rho)
        if rows_with_data:
            log(f"\n  Per-state Pearson rho of K_derived[s] vs K_emp[s]: "
                f"mean={np.mean(rows_with_data):.4f}, "
                f"median={np.median(rows_with_data):.4f}")

        # Sample comparison: r=3 b=10 (r=3 → r'=5 deterministic v=1)
        log(f"\n  Sample state-transition checks:")
        for r, b in [(3, 10), (5, 8), (21, 8), (1, 12), (15, 20)]:
            s = state_idx(r, b)
            d_row = K_derived[s, :]
            e_row = K_emp_full[s, :]
            top_d = np.argsort(-d_row)[:5]
            top_e = np.argsort(-e_row)[:5]
            log(f"\n    State (r={r}, b={b}):")
            for tag, top, row in [("DERIVED", top_d, d_row), ("EMPIRIC", top_e, e_row)]:
                top_str = ", ".join(
                    f"({(idx//B)*2+1},{idx%B}):{row[idx]:.4f}"
                    for idx in top if row[idx] > 1e-4
                )
                log(f"      {tag}: {top_str}")
    else:
        log("  R60 empirical kernel not found; skipping comparison")
        K_emp_full = None

    # ------ Step 5: Perron eigenvec of K_derived → D_pred ------
    log("\n" + "=" * 78)
    log("Step 5: Perron eigenvec of K_derived → D_pred(r)")
    log("=" * 78)
    pi_derived = stationary_via_iteration(K_derived, n_iter=2000)
    log(f"  pi_derived sum = {pi_derived.sum():.6f}")
    # Marginalize to residue
    rho_pred = np.zeros(R)
    for r_idx in range(R):
        rho_pred[r_idx] = pi_derived[r_idx * B:(r_idx + 1) * B].sum()

    # Load empirical D_avg
    D_avg = {}
    with open(OUT / "qsd_late_t_avg.csv") as f:
        f.readline()
        for line in f:
            parts = line.strip().split(",")
            D_avg[int(parts[0])] = float(parts[2])
    odd_r32 = list(range(1, 32, 2))

    # Load Chang pi_32 for normalization (precomputed in qsd_late_t_avg.csv)
    pi_32 = {}
    with open(OUT / "qsd_late_t_avg.csv") as f:
        f.readline()
        for line in f:
            parts = line.strip().split(",")
            pi_32[int(parts[0])] = float(parts[1])

    log(f"\n  D_pred(r) from K_derived Perron eigvec / pi_32(r):")
    log(f"  {'r':>3}  {'rho_pred':>10}  {'pi_32':>8}  {'D_pred':>8}  "
        f"{'D_avg':>8}  {'diff':>+8}")
    total_dev = 0.0
    devs = []
    D_pred_arr = []
    D_avg_arr = []
    for i, r in enumerate(odd_r32):
        rho = rho_pred[i] / rho_pred.sum() if rho_pred.sum() > 0 else 0
        d_pred = rho / pi_32[r]
        d_emp = D_avg[r]
        diff = d_pred - d_emp
        total_dev += abs(diff)
        devs.append(diff)
        D_pred_arr.append(d_pred)
        D_avg_arr.append(d_emp)
        log(f"  {r:>3}  {rho:>10.6f}  {pi_32[r]:>8.5f}  "
            f"{d_pred:>8.4f}  {d_emp:>8.4f}  {diff:>+8.4f}")
    pearson = float(np.corrcoef(D_pred_arr, D_avg_arr)[0, 1])
    log(f"\n  Total |D_pred - D_avg| = {total_dev:.4f}")
    log(f"  MAD per residue = {total_dev/16:.4f}")
    log(f"  Pearson ρ = {pearson:.4f}")
    log(f"  R60 empirical reference: total_dev = 3.40, Pearson = 0.80")

    # ------ Save outputs ------
    log("\n" + "=" * 78)
    log("Save outputs")
    log("=" * 78)
    sp.save_npz(OUT / "K_derived.npz", sp.csr_matrix(K_derived))
    log("  [wrote] K_derived.npz")

    with open(OUT / "D_predicted_derived.csv", "w") as f:
        f.write("r,pi_32,D_pred_derived,D_avg,diff\n")
        for i, r in enumerate(odd_r32):
            f.write(f"{r},{pi_32[r]:.6f},{D_pred_arr[i]:.6f},"
                    f"{D_avg_arr[i]:.6f},{D_pred_arr[i]-D_avg_arr[i]:+.6f}\n")
    log("  [wrote] D_predicted_derived.csv")

    if K_emp_full is not None:
        with open(OUT / "K_derived_vs_empirical.csv", "w") as f:
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
        log("  [wrote] K_derived_vs_empirical.csv")

    # ------ Verdict ------
    log("\n" + "=" * 78)
    log("VERDICT")
    log("=" * 78)
    log(f"\n  K_derived total deviation from D_avg: {total_dev:.4f}")
    log(f"  R60 empirical total deviation: 3.40")
    log(f"  Match quality: {'BETTER' if total_dev < 3.40 else 'WORSE' if total_dev > 4.0 else 'COMPARABLE'}")
    log(f"")
    log(f"  Per-state row Pearson (when comparable): see above")

    if total_dev < 3.5 and pearson > 0.75:
        log(f"  Outcome (α): K_derived ≈ K_empirical structurally; trajectory")
        log(f"               measure ↔ Perron(K_derived) is rigorous identification.")
    elif total_dev < 5.0 or pearson > 0.6:
        log(f"  Outcome (β): K_derived captures some structure but residuals")
        log(f"               suggest missing components.")
    else:
        log(f"  Outcome (γ): K_derived differs systematically from K_empirical;")
        log(f"               joint state assumption likely missing structure.")

    (OUT / "kernel_first_principles_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )
    log("\n  [wrote] kernel_first_principles_log.txt")


if __name__ == "__main__":
    main()
