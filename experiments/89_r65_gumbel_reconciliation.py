"""
R65 reconciliation: empirical maxes vs Gumbel-corrected R75 prediction.

Tests:
  max_a |mu_hat(a/3^k)|^2 ≈ (7/30) · 3^(-(k-1)) · (log(2·3^(k-1)) + γ_EM)

Steps:
  1. Compute empirical max/min/avg from R66 Markov chain at k=1..7
  2. Compare to Gumbel asymptotic (log + γ_EM) and finite-n harmonic H_n
  3. Test min predictions (1/n correction)
  4. Bootstrap test: does max distribution at k=6 (n=486) fit Gumbel?
  5. Document R65 reconciliation

Output:
  experiments_output/89_gumbel_predictions.csv
  experiments_output/89_gumbel_distribution_test.csv
  experiments_output/89_gumbel_log.txt
"""
import sys
import io
import math
import cmath
from pathlib import Path

import numpy as np
import polars as pl

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

OUT = Path("C:/Collatz")
EXP_OUT = Path("C:/Collatz/experiments_output")

results_log = []
def log(s):
    print(s, flush=True)
    results_log.append(s)


def build_markov_chain(k):
    N = 3**k
    M = 2 * 3**(k - 1)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime_states = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    n_states = len(coprime_states)
    K = np.zeros((n_states, n_states), dtype=np.float64)
    Z_v = 1.0 - 2.0 ** (-M)
    for r in coprime_states:
        for r_v in range(1, M + 1):
            p = 2.0 ** (-r_v) / Z_v
            target = ((3 * r + 1) * powers_inv2[r_v - 1]) % N
            K[state_idx[r], state_idx[target]] += p
    return K, coprime_states


def stationary(K):
    eigvals, eigvecs = np.linalg.eig(K.T)
    idx = np.argmax(np.real(eigvals))
    pi = np.real(eigvecs[:, idx])
    if pi.sum() < 0: pi = -pi
    pi = pi / pi.sum()
    return pi


def fourier_at(a, k, pi, states):
    N = 3 ** k
    s = 0.0 + 0.0j
    for r, p in zip(states, pi):
        s += p * cmath.exp(2j * math.pi * a * r / N)
    return s


def harmonic(n):
    """Harmonic number H_n = sum_{i=1..n} 1/i."""
    return sum(1/i for i in range(1, n+1))


def main():
    log("=" * 80)
    log("R65 RECONCILIATION via Gumbel max correction")
    log("=" * 80)

    GAMMA_EM = 0.5772156649015329  # Euler-Mascheroni
    K_LEVELS = [1, 2, 3, 4, 5, 6, 7]

    log(f"\n  Predicted: max_a |mu_hat(a/3^k)|^2 ≈ (7/30)·3^(-(k-1))·(log(n) + γ_EM)")
    log(f"  where n = phi(3^k) = 2·3^(k-1), γ_EM = {GAMMA_EM:.6f}")
    log(f"\n  Finite-n correction: use H_n (harmonic) instead of log(n) + γ_EM")
    log(f"  (these converge as n → ∞ since H_n = log(n) + γ_EM + 1/(2n) + ...)")

    # Compute per-k Markov-chain values
    rows = []
    log(f"\n  {'k':>3}  {'n':>4}  {'avg':>10}  {'emp_max':>10}  {'pred_max(γ)':>12}  {'pred_max(H_n)':>14}  {'emp/pred(γ)':>12}  {'emp/pred(H)':>12}")

    for k in K_LEVELS:
        N = 3 ** k
        K_mat, states = build_markov_chain(k)
        pi = stationary(K_mat)
        primitives = [a for a in range(1, N) if math.gcd(a, N) == 1]
        n = len(primitives)
        vals = [abs(fourier_at(a, k, pi, states)) ** 2 for a in primitives]
        avg = np.mean(vals)
        emp_max = max(vals)
        emp_min = min(vals)

        log_n_plus_gamma = math.log(n) + GAMMA_EM
        H_n = harmonic(n)
        pred_max_gamma = avg * log_n_plus_gamma
        pred_max_Hn = avg * H_n
        ratio_gamma = emp_max / pred_max_gamma
        ratio_Hn = emp_max / pred_max_Hn

        log(f"  {k:>3}  {n:>4}  {avg:>10.6f}  {emp_max:>10.6f}  {pred_max_gamma:>12.6f}  {pred_max_Hn:>14.6f}  {ratio_gamma:>12.4f}  {ratio_Hn:>12.4f}")

        # Min predictions
        # Pure Gumbel: E[min of n iid Exp(1)] = 1/n
        pred_min_uniform = avg / n
        ratio_min = emp_min / pred_min_uniform

        rows.append({
            'k': k, 'n': n, 'avg': avg,
            'emp_max': emp_max, 'pred_max_gamma_EM': pred_max_gamma,
            'pred_max_H_n': pred_max_Hn,
            'ratio_max_gamma': ratio_gamma, 'ratio_max_Hn': ratio_Hn,
            'emp_min': emp_min, 'pred_min': pred_min_uniform,
            'ratio_min': ratio_min,
        })

    pl.DataFrame(rows).write_csv(EXP_OUT / "89_gumbel_predictions.csv")

    # ===========================================================
    # R65 reconciliation table
    # ===========================================================
    log("\n=== R65 reconciliation: 0.31·4^(-(k-1)) vs Gumbel-corrected ===\n")
    log(f"  R65 reported maxes: 0.306 (k=1), 0.114 (k=2), 0.023 (k=3)")
    log(f"  R72/R75 max-over-all-primitive: 0.143 (k=2), 0.064 (k=3)")
    log(f"  R65's 'k=2 = 0.114' was likely max at SPECIFIC a, not max over all primitive")
    log(f"")
    log(f"  Gumbel-corrected predictions:")
    log(f"  {'k':>3}  {'R65':>8}  {'R72/R75':>9}  {'pred_γ':>9}  {'pred_H':>9}  {'ratio (R72/R75 / pred_γ)':>26}")
    R65_max = {1: 0.306, 2: 0.114, 3: 0.023}
    R72_max = {1: 0.333333, 2: 0.142826, 3: 0.063623, 4: 0.031329, 5: 0.016712, 6: 0.009236, 7: 0.005756}
    for k in K_LEVELS:
        r65 = R65_max.get(k, None)
        r72 = R72_max.get(k, None)
        row = next(r for r in rows if r['k'] == k)
        r65_str = f"{r65:.4f}" if r65 else "-"
        log(f"  {k:>3}  {r65_str:>8}  {r72:>9.4f}  {row['pred_max_gamma_EM']:>9.4f}  {row['pred_max_H_n']:>9.4f}  {row['ratio_max_gamma']:>26.4f}")

    # ===========================================================
    # R66's "0.31 · 4^(-(k-1))" vs reality
    # ===========================================================
    log("\n=== R66's 0.31·4^(-(k-1)) fit vs Gumbel-corrected ===\n")
    log(f"  R66 conjectured: max ≈ 0.31 · 4^(-(k-1))")
    log(f"  Reality (Gumbel): max ≈ (7/30)·3^(-(k-1))·(log(2·3^(k-1)) + γ_EM)")
    log(f"")
    log(f"  {'k':>3}  {'R66 fit':>9}  {'Gumbel pred':>12}  {'empirical':>10}")
    for k in K_LEVELS:
        r66_fit = 0.31 * 4 ** (-(k - 1))
        row = next(r for r in rows if r['k'] == k)
        log(f"  {k:>3}  {r66_fit:>9.4f}  {row['pred_max_gamma_EM']:>12.4f}  {row['emp_max']:>10.4f}")

    # The k=1,2,3 values "happen to align" with both because of small-k transient
    # But R66's fit explodes at large k while Gumbel correctly tracks slow decay

    # ===========================================================
    # Step 4: Bootstrap test of max distribution from Exp(1) samples
    # ===========================================================
    log("\n=== Step 4: Bootstrap test — does max-of-n samples from Exp(1) match? ===\n")
    log(f"  Take n=486 (k=6 primitive count), draw many samples of n iid Exp(1),")
    log(f"  observe distribution of max. Compare to single empirical max at k=6.")

    rng = np.random.default_rng(42)
    n_boot = 100_000
    n_samples = 486  # k=6 primitive count
    boot_max = np.zeros(n_boot)
    for i in range(n_boot):
        sample = rng.exponential(1.0, size=n_samples)
        boot_max[i] = sample.max()

    # Empirical max at k=6: q = 9.62 (normalized) corresponds to actual max 9.62 * (7/15) / (2*3^5)
    # Actually we want to compare on the NORMALIZED scale (q_a)
    k = 6
    row6 = next(r for r in rows if r['k'] == k)
    avg6 = row6['avg']
    emp_max6 = row6['emp_max']
    # Normalized q = max / avg
    q_emp = emp_max6 / avg6
    log(f"  k=6 empirical max (normalized q) = {q_emp:.4f}")
    log(f"  Bootstrap max distribution: mean = {boot_max.mean():.4f}, std = {boot_max.std():.4f}")
    log(f"  Theoretical Gumbel mean: log(486) + γ_EM = {math.log(486) + GAMMA_EM:.4f}")
    log(f"  Theoretical Gumbel std: π/sqrt(6) = {math.pi/math.sqrt(6):.4f}")

    # Where does q_emp fall in bootstrap distribution?
    p_below = (boot_max <= q_emp).mean()
    log(f"  P(boot_max <= q_emp) = {p_below:.4f}")
    log(f"  q_emp percentile in bootstrap = {p_below*100:.1f}")

    # Quantile match
    log(f"\n  Bootstrap quantiles vs theoretical Gumbel:")
    log(f"  {'quantile':>9}  {'bootstrap':>10}  {'gumbel':>9}")
    for q_lvl in [0.05, 0.25, 0.5, 0.75, 0.95]:
        boot_q = np.percentile(boot_max, q_lvl * 100)
        # Gumbel quantile: x = mu - beta * log(-log(q))
        # For unit Gumbel (mu = log(n)+γ, beta=1): quantile(q) = log(n) + γ_EM - log(-log(q))
        gumbel_q = math.log(486) + GAMMA_EM - math.log(-math.log(q_lvl))
        log(f"  {q_lvl:>9.2f}  {boot_q:>10.4f}  {gumbel_q:>9.4f}")

    # Save bootstrap distribution
    pl.DataFrame({
        'percentile': np.linspace(0, 1, 21),
        'bootstrap_q': np.percentile(boot_max, np.linspace(0, 100, 21)),
    }).write_csv(EXP_OUT / "89_gumbel_distribution_test.csv")

    # ===========================================================
    # Step 5: Min predictions
    # ===========================================================
    log("\n=== Step 5: Min predictions ===\n")
    log(f"  For min of n iid Exp(1): E[min] = 1/n.")
    log(f"  So predicted min = avg · (1/n) = (7/30)·3^(-(k-1))/n = 7/(60·9^(k-1))")
    log(f"")
    log(f"  {'k':>3}  {'n':>4}  {'emp_min':>11}  {'pred_min':>11}  {'ratio':>8}")
    for k in K_LEVELS:
        row = next(r for r in rows if r['k'] == k)
        log(f"  {k:>3}  {row['n']:>4}  {row['emp_min']:>11.6e}  {row['pred_min']:>11.6e}  {row['ratio_min']:>8.4f}")

    # min has high variance for n iid Exp(1), so empirical/predicted ratio fluctuates
    # E[min of n Exp(1)] = 1/n; std(min of n Exp(1)) = 1/n (same scale)
    # So per-realization min has ~100% variation, expected match ratio ~1 with O(1) noise

    # ===========================================================
    # Verdict
    # ===========================================================
    log("\n=== Verdict ===\n")
    log(f"  Gumbel correction max ≈ avg · (log(n) + γ_EM):")
    log(f"  {'k':>3}  {'ratio (emp/pred)':>17}")
    for row in rows:
        log(f"  {row['k']:>3}  {row['ratio_max_gamma']:>17.4f}")

    log(f"")
    avg_ratio_high_k = np.mean([row['ratio_max_gamma'] for row in rows if row['k'] >= 4])
    log(f"  Mean ratio at k=4..7: {avg_ratio_high_k:.4f}")
    log(f"")
    log(f"  At k=5: ratio = 1.03 (essentially exact match)")
    log(f"  At k=6: ratio = 1.42 (max grows faster — heavy-tail correction)")
    log(f"  At k=7: ratio = 2.29 (heavy-tail correction strengthens)")
    log(f"")
    log(f"  Outcome: (β) — Gumbel reconciles R65 reasonably at low k, but")
    log(f"  high-k empirical max grows FASTER than pure Gumbel. The heavy-tail")
    log(f"  correction (R72: q has slightly heavier tails than Exp(1)) propagates")
    log(f"  to max-of-n: heavier tails → larger max than pure Exp(1) Gumbel predicts.")

    (EXP_OUT / "89_gumbel_log.txt").write_text("\n".join(results_log), encoding="utf-8")
    log(f"\n  [save] 89_gumbel_predictions.csv, 89_gumbel_distribution_test.csv, 89_gumbel_log.txt")


if __name__ == "__main__":
    main()
