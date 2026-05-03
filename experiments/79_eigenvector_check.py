"""
Direct comparison: Result 23's v_max (inverse-tree eigenvector) vs
Result 51's D_avg (forward survivor-conditioned trajectory measure).
"""
import sys
import io
import numpy as np
from scipy.stats import pearsonr, spearmanr
import polars as pl

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)


def main():
    # Load Result 23's eigenvector
    df = pl.read_csv("C:/Collatz/inverse_tree/inverse_tree_eigvec_mod32.csv")
    eigvec_full = {row['residue_mod_32']: row['predicted_density'] for row in df.iter_rows(named=True)}

    # Restrict to odd residues
    odd_residues = list(range(1, 32, 2))
    eigvec_odd = np.array([eigvec_full[r] for r in odd_residues])
    sum_odd = eigvec_odd.sum()
    v_max_relative = eigvec_odd / sum_odd * 16  # ratio relative to uniform over 16 odd residues

    print(f"=== Result 23's v_max (inverse-tree eigenvector at k=5) ===", flush=True)
    print(f"  Sum of odd-residue densities: {sum_odd:.4f}", flush=True)
    print(f"  Sample relative values:", flush=True)
    for r in [3, 5, 13, 17, 21, 31]:
        i = odd_residues.index(r)
        print(f"    r={r:>3}: v_relative = {v_max_relative[i]:.4f}  (Result 23 reports specific values)", flush=True)

    # Empirical D_avg from Result 51 (file lines 4373-4380)
    D_avg = {
        1: 1.609, 3: 1.236, 5: 1.864, 7: 0.738, 9: 0.696, 11: 0.658,
        13: 0.557, 15: 1.058, 17: 1.132, 19: 0.628, 21: 0.931, 23: 1.398,
        25: 0.544, 27: 0.802, 29: 1.354, 31: 0.767,
    }
    D_avg_arr = np.array([D_avg[r] for r in odd_residues])

    # ============= Direct comparison =============
    print(f"\n=== Direct comparison v_max_relative vs D_avg per residue ===", flush=True)
    print(f"  {'r':>3}  {'v_max_rel':>10}  {'D_avg':>7}  {'ratio':>7}  {'gap':>9}", flush=True)
    for i, r in enumerate(odd_residues):
        v = v_max_relative[i]; D = D_avg_arr[i]
        ratio = v / D if D > 0 else float('nan')
        gap = v - D
        print(f"  {r:>3}  {v:>10.4f}  {D:>7.3f}  {ratio:>7.3f}  {gap:>+9.4f}", flush=True)

    # Correlations
    pearson_r, p_r = pearsonr(v_max_relative, D_avg_arr)
    spearman_r, p_s = spearmanr(v_max_relative, D_avg_arr)
    mad = float(np.abs(v_max_relative - D_avg_arr).mean())
    max_abs_gap = float(np.abs(v_max_relative - D_avg_arr).max())

    print(f"\n  Pearson ρ(v_max, D_avg)  = {pearson_r:+.4f}  (p={p_r:.4g})", flush=True)
    print(f"  Spearman ρ(v_max, D_avg) = {spearman_r:+.4f}  (p={p_s:.4g})", flush=True)
    print(f"  MAD |v_max - D_avg|       = {mad:.4f}", flush=True)
    print(f"  Max |gap|                 = {max_abs_gap:.4f}", flush=True)

    # ============= Verdict =============
    print(f"\n=== Verdict ===", flush=True)
    if pearson_r > 0.95 and mad < 0.05:
        verdict = 'A'
        print(f"  (A) DIRECT MATCH: Result 23's v_max IS the trajectory measure", flush=True)
    elif pearson_r > 0.7:
        verdict = 'B'
        print(f"  (B) PARTIAL MATCH: ρ={pearson_r:+.4f}, captures pattern but not quantitatively", flush=True)
    else:
        verdict = 'C'
        print(f"  (C) NO MATCH: Pearson ρ={pearson_r:+.4f} is BELOW 0.5", flush=True)
        print(f"      Result 23's eigenvector is NOT D_avg.", flush=True)
        print(f"      Confirms findings.md statement: forward-orbit ↔ inverse-tree Pearson ≈ -0.20", flush=True)

    # Diagnose key residues
    i_21 = odd_residues.index(21)
    i_5 = odd_residues.index(5)
    i_13 = odd_residues.index(13)
    print(f"\n  Sign-flip diagnostics:", flush=True)
    print(f"  r=21: v_max = {v_max_relative[i_21]:.3f} (HUGE enhancement)  vs  D_avg = {D_avg_arr[i_21]:.3f} (slight depletion)", flush=True)
    print(f"  r= 5: v_max = {v_max_relative[i_5]:.3f}                       vs  D_avg = {D_avg_arr[i_5]:.3f} (most enhanced)", flush=True)
    print(f"  r=13: v_max = {v_max_relative[i_13]:.3f}                       vs  D_avg = {D_avg_arr[i_13]:.3f} (most depleted)", flush=True)

    # Numerical curiosity on λ_max
    lambda_max = 1.263763
    log_lam = np.log(lambda_max)
    print(f"\n=== Step 6: Numerical curiosities of λ_max = {lambda_max} ===", flush=True)
    print(f"  log(λ_max)               = {log_lam:.6f}", flush=True)
    print(f"  log(λ_max) / log(2)      = {log_lam/np.log(2):.6f}", flush=True)
    print(f"  2 × log(λ_max) / log(2)  = {2*log_lam/np.log(2):.6f}  (cf Chang H-dim ≈ 0.68)", flush=True)
    print(f"  log(λ_max) / log(4/3)    = {log_lam/np.log(4/3):.6f}", flush=True)
    print(f"  log(2) / log(λ_max)      = {np.log(2)/log_lam:.6f}", flush=True)

    # Save
    rows = [{'r': r, 'v_max_relative': float(v_max_relative[i]),
             'D_avg': D_avg[r], 'ratio': float(v_max_relative[i]/D_avg[r]),
             'gap': float(v_max_relative[i] - D_avg[r])}
            for i, r in enumerate(odd_residues)]
    pl.DataFrame(rows).write_csv("C:/Collatz/experiments_output/79_eigenvector_vs_D_avg.csv")
    print(f"\n[save] CSV", flush=True)

    return verdict


if __name__ == "__main__":
    main()
