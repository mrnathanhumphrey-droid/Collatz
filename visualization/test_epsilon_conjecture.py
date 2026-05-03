"""
test_epsilon_conjecture.py - Test whether ε(σ) decomposes via per-j gap.

Conjecture: ε(σ) = per_j_gap × ⟨ℓ⟩  (with per_j_gap = -0.26 invariant in k).

Discriminator: k-invariance.
  Empirically ε(σ) ≈ -2.45 invariant in k.
  ⟨ℓ⟩ at k = 6,8,10,12,14 = 9.5, 12.5, 15.5, 18.5, 21.5 (grows linearly).
  If ε(σ) = -0.26 × ⟨ℓ⟩, ε(σ) would grow with k. It doesn't. Rules out simple form.

But the per-j ε(j) structure could still be linear in j with slope -0.26 *if* the
intercept shifts with k to keep ⟨ε⟩ invariant. This is the actual test:
  1. Fit ε(j) = α_k + β_k · j at each k.
  2. Is β_k ≈ -0.26 invariant in k? (predicts per-j gap is per-class linear)
  3. Does α_k shift with k by exactly 0.13 per unit k? (predicts ⟨ε⟩ stays -2.45)
"""
import importlib.util
import math
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[1]
EXP29 = ROOT / "experiments" / "29_qx1_cycle_classification.py"
spec = importlib.util.spec_from_file_location("exp29", EXP29)
exp29 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exp29)
qx1_prefix = exp29.qx1_prefix

K_H = 3.0 / (math.log(4.0) - math.log(3.0))


def precompute_j_table(k):
    j_arr = np.zeros(1 << k, dtype=np.int8)
    for r in range(1, 1 << k, 2):
        _, a_star, _ = qx1_prefix(r, k, 3)
        j = 0
        a = a_star
        while a > 1 and a % 3 == 0:
            a //= 3
            j += 1
        j_arr[r] = j if a == 1 else -1
    return j_arr


def main():
    df = pl.read_csv(Path(__file__).parent / "viz_outputs" / "descent_b_enlarged.csv")
    n_arr = df["n"].to_numpy()
    log_n = df["log_n"].to_numpy()
    sigma = df["sigma"].to_numpy().astype(np.float64)
    print(f"Loaded {len(df):,} orbits.  N range: [3, {n_arr.max():,}]  (log N_max = {math.log(n_arr.max()):.2f})")
    print(f"K_H = 3/log(4/3) = {K_H:.5f}")
    print(f"Heuristic per-j shift (1 + K_H·log(3)) = {1 + K_H*math.log(3):.5f}")
    print()

    print(f"{'='*100}")
    print(f"{'k':>3}  {'<ℓ>':>6}  {'<eps>_global':>13}  {'OLS slope on j':>15}  "
          f"{'eps(j) intercept':>17}  {'eps(j) slope β_k':>17}  {'predicted const_k':>18}")
    print(f"{'='*100}")

    rows = []
    for k in [6, 8, 10, 12, 14, 16]:
        j_table = precompute_j_table(k)
        j = j_table[n_arr % (1 << k)].astype(np.int64)
        ell = k + j  # prefix length per orbit
        ell_mean = ell.mean()

        # alpha_det per orbit (deterministic Tao prediction)
        alpha_det = ell + K_H * (j * math.log(3) - k * math.log(2))
        # epsilon per orbit
        eps = sigma - alpha_det - K_H * log_n
        eps_global = eps.mean()

        # OLS for per-j gap (matches sweep_k.py output)
        X = np.column_stack([np.ones(len(eps)), log_n, j.astype(np.float64)])
        beta, *_ = np.linalg.lstsq(X, sigma, rcond=None)
        ols_j_slope = beta[2]

        # Per-class eps(j): mean ε per j
        unique_j = sorted(set(j.tolist()))
        eps_per_j = []
        n_per_j = []
        for jv in unique_j:
            mask = j == jv
            if mask.sum() > 50:
                eps_per_j.append(eps[mask].mean())
                n_per_j.append(int(mask.sum()))
            else:
                eps_per_j.append(None)
                n_per_j.append(int(mask.sum()))

        # Linear fit eps(j) = α_k + β_k · j (only on classes with enough orbits)
        js_use = [jv for jv, v in zip(unique_j, eps_per_j) if v is not None]
        eps_use = [v for v in eps_per_j if v is not None]
        if len(js_use) >= 3:
            slope, intercept = np.polyfit(js_use, eps_use, 1)
        else:
            slope, intercept = float("nan"), float("nan")

        # Predicted const_k if ε(j) = -0.26·j + const_k and ⟨ε⟩=-2.45
        # const_k = -2.45 + 0.26·⟨j⟩ where ⟨j⟩ = (k+1)/2 for binomial
        predicted_const_k = -2.45 + 0.26 * (k + 1) / 2

        print(f"{k:>3}  {ell_mean:>6.2f}  {eps_global:>13.4f}  {ols_j_slope:>+15.4f}  "
              f"{intercept:>17.4f}  {slope:>+17.4f}  {predicted_const_k:>+18.4f}")
        rows.append((k, ell_mean, eps_global, ols_j_slope, intercept, slope, predicted_const_k, list(zip(js_use, eps_use))))

    print(f"{'='*100}")
    print(f"\nPer-class ε(j) tables:")
    print(f"{'='*60}")
    for k, ell_mean, eps_g, ols_s, ic, sl, pc, eps_table in rows:
        print(f"\nk={k}: ε(j) per class  (ε_global = {eps_g:.4f}, slope = {sl:+.4f})")
        print(f"  {'j':>3}  {'ε(j)':>9}  {'predicted (-0.26·(j-<j>) + ε_global)':>40}")
        for jv, ev in eps_table:
            pred = -0.26 * (jv - (k+1)/2) + eps_g
            print(f"  {jv:>3}  {ev:>+9.4f}  {pred:>+40.4f}")

    # Direct test of conjecture: ε(σ) vs predicted -0.26·⟨ℓ⟩
    print(f"\n{'='*100}")
    print(f"Direct conjecture test: ε(σ) vs -0.26 · ⟨ℓ⟩")
    print(f"{'='*100}")
    print(f"  {'k':>3}  {'<ℓ>':>6}  {'-0.26·<ℓ>':>11}  {'observed ε':>12}  {'gap':>10}")
    for k, ell_mean, eps_g, *_ in rows:
        pred = -0.26 * ell_mean
        gap = eps_g - pred
        print(f"  {k:>3}  {ell_mean:>6.2f}  {pred:>+11.4f}  {eps_g:>+12.4f}  {gap:>+10.4f}")
    print(f"\nIf conjecture holds: gap should be small AND k-invariant.")
    print(f"If conjecture fails (k-invariance): gap will grow linearly with k.")


if __name__ == "__main__":
    main()
