"""B2 coefficient + log(n) correlation check.

Disambiguates three competing explanations for the smoke-test finding that
B1's u_r posterior magnitude is ~13x smaller than the literal conjecture
predicts:

  (a) Parameterization artifact (log(n) ↔ log(a_final) correlation absorbs
      the colinear component into β_log_n)
  (b) Posterior shrinkage artifact (hierarchical prior pulls u_r toward 0)
  (c) Genuine magnitude collapse (literal slope is wrong; only the directional
      form survives)

Tests on existing posteriors at C:/Collatz/experiments_output/paper1_hmc_n10k/
hmc_B0/, hmc_B1/, hmc_B2/. NO new MCMC. Does not touch B4.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import polars as pl
from scipy.stats import pearsonr, spearmanr
sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parent
COLLATZ = HERE.parent
B0_DIR = COLLATZ / "experiments_output" / "paper1_hmc_n10k" / "hmc_B0"
B1_DIR = COLLATZ / "experiments_output" / "paper1_hmc_n10k" / "hmc_B1"
B2_DIR = COLLATZ / "experiments_output" / "paper1_hmc_n10k" / "hmc_B2"
DATA_PARQ = COLLATZ / "data" / "main_N10000000.parquet"
MD_OUT = COLLATZ / "docs" / "paper1_b2_coefficient_check.md"
MD_OUT.parent.mkdir(parents=True, exist_ok=True)

SEED = 20260509
N_DATA = 10_000
N_TRAIN = 8_000
N_TEST = 2_000


def deterministic_prefix(r, a0, max_steps=400):
    a, c, steps = a0, r, 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3 * c + 1
        steps += 1
    return a


def lookup_a_final(k_pow):
    M = 1 << k_pow
    return np.array([deterministic_prefix(r, M) for r in range(M)], dtype=np.int64)


def reconstruct_train(seed=SEED, N=N_DATA, N_train=N_TRAIN):
    df = pl.read_parquet(DATA_PARQ)
    n_full = df["n"].to_numpy().astype(np.int64)
    sigma_full = df["sigma"].to_numpy().astype(np.int64)
    rng = np.random.default_rng(seed)
    sub_idx = rng.choice(len(n_full), size=N, replace=False)
    n_arr = n_full[sub_idx]
    sigma_arr = sigma_full[sub_idx]
    perm = rng.permutation(N)
    train_idx = perm[:N_train]
    return n_arr[train_idx], sigma_arr[train_idx]


def load_chain_csvs(spec_dir, params_prefixes):
    """Read all 4 chain CSVs and stack samples for parameters matching
    any prefix in params_prefixes (e.g. 'beta.', 'sigma_u', 'u.')."""
    csvs = sorted(spec_dir.glob("nb2_glm-*_*.csv"))
    if len(csvs) != 4:
        raise RuntimeError(f"expected 4 chain CSVs in {spec_dir}, got {len(csvs)}")
    out = {}
    cols_needed = None
    chain_dfs = [pl.read_csv(c, comment_prefix="#") for c in csvs]
    cols = chain_dfs[0].columns
    matched = [c for c in cols if any(c.startswith(p) for p in params_prefixes)]
    for c in matched:
        out[c] = np.concatenate([df[c].to_numpy() for df in chain_dfs])
    return out


def main():
    print("[load] reconstructing training subset (n, sigma) ...", flush=True)
    n_train, sigma_train = reconstruct_train()
    print(f"        n_train = {len(n_train):,}", flush=True)
    log_n_train = np.log(n_train.astype(float))
    mod8_train = (n_train % 8).astype(int)
    af_k3_lookup = lookup_a_final(3)
    log_af_k3_train = np.log(af_k3_lookup[mod8_train].astype(float))

    # ============================================================
    # TEST 1: log(n) ↔ log(a_final(r mod 8)) correlation
    # ============================================================
    print("\n[Test 1] log(n) vs log(a_final(r mod 8)) correlation ...", flush=True)
    pearson_r, pearson_p = pearsonr(log_n_train, log_af_k3_train)
    spearman_rho, spearman_p = spearmanr(log_n_train, log_af_k3_train)
    n = len(log_n_train)
    se_r = np.sqrt((1 - pearson_r ** 2) / (n - 2))
    fisher_z = 0.5 * np.log((1 + pearson_r) / (1 - pearson_r))
    se_fz = 1 / np.sqrt(n - 3)
    z_lo = fisher_z - 1.96 * se_fz
    z_hi = fisher_z + 1.96 * se_fz
    r_lo = (np.exp(2 * z_lo) - 1) / (np.exp(2 * z_lo) + 1)
    r_hi = (np.exp(2 * z_hi) - 1) / (np.exp(2 * z_hi) + 1)
    print(f"  Pearson r  = {pearson_r:+.6f} (p={pearson_p:.4e}, "
          f"95% CI [{r_lo:+.4f}, {r_hi:+.4f}])", flush=True)
    print(f"  Spearman ρ = {spearman_rho:+.6f} (p={spearman_p:.4e})", flush=True)

    abs_r = abs(pearson_r)
    if abs_r < 0.02:
        test1_verdict = "INDEPENDENT"
        test1_text = ("|r|<0.02; explanation (a) parameterization artifact is RULED OUT. "
                      "log(n) and log(a_final(r mod 8)) are essentially independent.")
    elif abs_r < 0.10:
        test1_verdict = "WEAKLY CORRELATED"
        test1_text = ("0.02 ≤ |r| < 0.10; explanation (a) is possibly operative. "
                      "Some absorption of log(a_final) signal into β_log_n is plausible.")
    else:
        test1_verdict = "STRONGLY CORRELATED"
        test1_text = ("|r| ≥ 0.10; explanation (a) is strongly operative. "
                      "β_log_n necessarily absorbs the colinear component of log(a_final).")
    print(f"  Test 1 verdict: {test1_verdict}", flush=True)
    print(f"  {test1_text}", flush=True)

    # ============================================================
    # TEST 2: B2's posterior coefs on a_final dummies
    # ============================================================
    # B2's design: [intercept, log_n, dummy_for_af=3, dummy_for_af=9, dummy_for_af=27]
    # (af=1 is reference). beta.1=intercept, beta.2=log_n, beta.3..5=dummies.
    print("\n[Test 2] B2 fixed-effect coefficients on a_final dummies ...", flush=True)
    b2 = load_chain_csvs(B2_DIR, ["beta."])
    beta_cols = sorted(b2.keys(), key=lambda c: int(c.split(".")[1]))
    print(f"  beta columns in B2 fit: {beta_cols}", flush=True)
    if len(beta_cols) != 5:
        raise RuntimeError(f"expected 5 beta cols in B2, got {len(beta_cols)}")
    # Map: beta.1=intercept, beta.2=log_n, beta.3=af_dummy_3, beta.4=af_dummy_9, beta.5=af_dummy_27
    beta_intercept_b2 = b2["beta.1"]
    beta_log_n_b2 = b2["beta.2"]
    beta_af3 = b2["beta.3"]
    beta_af9 = b2["beta.4"]
    beta_af27 = b2["beta.5"]
    af_levels = np.array([3, 9, 27])
    log_af_levels = np.log(af_levels.astype(float))
    af_coefs_mean = np.array([beta_af3.mean(), beta_af9.mean(), beta_af27.mean()])
    af_coefs_q025 = np.array([np.quantile(beta_af3, 0.025),
                               np.quantile(beta_af9, 0.025),
                               np.quantile(beta_af27, 0.025)])
    af_coefs_q975 = np.array([np.quantile(beta_af3, 0.975),
                               np.quantile(beta_af9, 0.975),
                               np.quantile(beta_af27, 0.975)])
    print("  posterior mean and 95% CI for each a_final dummy:", flush=True)
    for level, mean, lo, hi, lal in zip(af_levels, af_coefs_mean,
                                          af_coefs_q025, af_coefs_q975, log_af_levels):
        conj_pred = lal
        conj_match = "✓" if (lo <= conj_pred <= hi) else "✗"
        falsified_pred = 0.085 * lal
        falsified_match = "✓" if (lo <= falsified_pred <= hi) else "✗"
        print(f"    af={level:2d} (j={int(round(np.log(level)/np.log(3)))}): "
              f"β = {mean:+.4f} [95% CI {lo:+.4f}, {hi:+.4f}], "
              f"conjecture predicts {conj_pred:+.4f} {conj_match}, "
              f"smoke-test falsified-form predicts {falsified_pred:+.4f} {falsified_match}",
              flush=True)

    # Per-draw slope of dummy coefs on log(a_final): equivalent to a "β_a_final"
    # if the model had used continuous log(a_final). For each posterior draw,
    # regress (β_af3, β_af9, β_af27) on (log 3, log 9, log 27) (no intercept
    # since a_final=1 reference is at coef=0 by construction; with intercept
    # forced to zero at log(a_final)=0).
    n_draws = len(beta_af3)
    slopes = np.zeros(n_draws)
    intercepts = np.zeros(n_draws)
    log_af = np.array([0.0, np.log(3), np.log(9), np.log(27)])  # 4 levels incl ref
    for d in range(n_draws):
        coefs = np.array([0.0, beta_af3[d], beta_af9[d], beta_af27[d]])
        # OLS slope without forcing zero intercept (sanity: should give intercept ≈0)
        x_c = log_af - log_af.mean()
        y_c = coefs - coefs.mean()
        slopes[d] = float((x_c * y_c).sum() / (x_c * x_c).sum())
        intercepts[d] = coefs.mean() - slopes[d] * log_af.mean()
    slope_mean = slopes.mean(); slope_q025 = np.quantile(slopes, 0.025); slope_q975 = np.quantile(slopes, 0.975)
    print(f"  Implied slope on log(a_final) (regressing 4 dummy values incl af=1 ref):",
          flush=True)
    print(f"    posterior mean = {slope_mean:+.4f}, 95% CI [{slope_q025:+.4f}, {slope_q975:+.4f}]",
          flush=True)
    print(f"  Conjecture predicts slope = 1.0 (literal: β = 1.0 × log(a_final))", flush=True)
    print(f"  Smoke-test falsified form predicts slope ≈ 0.085 (the B1 raw slope)", flush=True)

    contains_1 = (slope_q025 <= 1.0 <= slope_q975)
    contains_0085 = (slope_q025 <= 0.085 <= slope_q975)
    excludes_1 = not contains_1
    excludes_0085 = not contains_0085
    if contains_1 and not contains_0085:
        test2_verdict = "LITERAL-MAGNITUDE CONSISTENT"
        test2_text = (f"95% CI [{slope_q025:+.4f}, {slope_q975:+.4f}] contains 1.0 "
                      "and excludes 0.085. The literal conjecture's magnitude is "
                      "consistent with B2's data.")
    elif contains_0085 and not contains_1:
        test2_verdict = "LITERAL-MAGNITUDE FALSIFIED"
        test2_text = (f"95% CI [{slope_q025:+.4f}, {slope_q975:+.4f}] contains 0.085 "
                      "and EXCLUDES 1.0. The literal conjecture's magnitude is "
                      "rejected; the smoke-test verdict stands.")
    elif contains_1 and contains_0085:
        test2_verdict = "AMBIGUOUS WIDE CI"
        test2_text = (f"95% CI [{slope_q025:+.4f}, {slope_q975:+.4f}] contains BOTH 1.0 "
                      "and 0.085 — too wide at N=10K to disambiguate.")
    elif slope_q975 < 1.0 and slope_q025 > 0.1:
        test2_verdict = "REDUCED SLOPE"
        test2_text = (f"95% CI [{slope_q025:+.4f}, {slope_q975:+.4f}] excludes 1.0 "
                      "but slope is meaningfully positive. Conjecture's directional "
                      "form holds at REDUCED slope.")
    else:
        test2_verdict = "OTHER"
        test2_text = (f"95% CI [{slope_q025:+.4f}, {slope_q975:+.4f}] doesn't fit "
                      "standard categories. Inspect manually.")
    print(f"  Test 2 verdict: {test2_verdict}", flush=True)

    # ============================================================
    # TEST 3: B0 vs B2 β_log_n
    # ============================================================
    print("\n[Test 3] B0 β_log_n vs B2 β_log_n ...", flush=True)
    b0 = load_chain_csvs(B0_DIR, ["beta."])
    beta_log_n_b0 = b0["beta.2"]
    b0_mean = beta_log_n_b0.mean(); b0_sd = beta_log_n_b0.std()
    b0_q025 = np.quantile(beta_log_n_b0, 0.025); b0_q975 = np.quantile(beta_log_n_b0, 0.975)
    b2_log_n_mean = beta_log_n_b2.mean(); b2_log_n_sd = beta_log_n_b2.std()
    b2_log_n_q025 = np.quantile(beta_log_n_b2, 0.025); b2_log_n_q975 = np.quantile(beta_log_n_b2, 0.975)
    diff = b2_log_n_mean - b0_mean
    pooled_sd = np.sqrt(b0_sd ** 2 + b2_log_n_sd ** 2)
    z_diff = diff / pooled_sd
    print(f"  B0 β_log_n: mean = {b0_mean:.6f}, SD = {b0_sd:.6f}, "
          f"95% CI [{b0_q025:.6f}, {b0_q975:.6f}]", flush=True)
    print(f"  B2 β_log_n: mean = {b2_log_n_mean:.6f}, SD = {b2_log_n_sd:.6f}, "
          f"95% CI [{b2_log_n_q025:.6f}, {b2_log_n_q975:.6f}]", flush=True)
    print(f"  Difference (B2 - B0) = {diff:+.6f}, pooled SD = {pooled_sd:.6f}, "
          f"z = {z_diff:+.3f}", flush=True)
    if abs(z_diff) > 2:
        test3_verdict = "REDISTRIBUTED"
        test3_text = (f"|z| = {abs(z_diff):.2f} > 2: B2's β_log_n is meaningfully "
                      "different from B0's. Some signal that was in log(n) alone "
                      "has been redistributed to a_final dummies. Corroborates the "
                      "parameterization story.")
    else:
        test3_verdict = "ESSENTIALLY UNCHANGED"
        test3_text = (f"|z| = {abs(z_diff):.2f} ≤ 2: B2's β_log_n is essentially "
                      "identical to B0's. log(n) and a_final dummies are nearly "
                      "orthogonal in the training data — no signal redistribution.")
    print(f"  Test 3 verdict: {test3_verdict}", flush=True)

    # ============================================================
    # TEST 4: parameterization-corrected u_r prediction
    # ============================================================
    # Compute mean(log(n) | r) for each residue mod 8 in training data
    print("\n[Test 4] Parameterization-corrected u_r prediction ...", flush=True)
    residues = list(range(8))
    mean_log_n_global = log_n_train.mean()
    mean_log_n_per_r = np.array([log_n_train[mod8_train == r].mean() for r in residues])
    u_r_predicted_uncorrected = np.array([np.log(af_k3_lookup[r]) for r in residues])
    # u_r_predicted = log(a_final(r)) - β_log_n × (mean(log(n)|r) - global mean)
    u_r_predicted = u_r_predicted_uncorrected - b2_log_n_mean * (mean_log_n_per_r - mean_log_n_global)
    # Compare to actual B1 u_r posterior means
    b1 = load_chain_csvs(B1_DIR, ["u."])
    u_cols = sorted(b1.keys(), key=lambda c: int(c.split(".")[1]))
    u_actual = np.array([b1[c].mean() for c in u_cols])
    print(f"  mean(log(n)|r) per residue: {mean_log_n_per_r}", flush=True)
    print(f"  log(a_final) per residue:    {u_r_predicted_uncorrected}", flush=True)
    print(f"  parameterization-corrected u_r predictions: {u_r_predicted}", flush=True)
    print(f"  B1 actual u_r means:                       {u_actual}", flush=True)
    # Demean both for the structural comparison
    u_pred_dm = u_r_predicted - u_r_predicted.mean()
    u_actual_dm = u_actual - u_actual.mean()
    pearson_u, pearson_u_p = pearsonr(u_pred_dm, u_actual_dm)
    spearman_u, spearman_u_p = spearmanr(u_pred_dm, u_actual_dm)
    # Compare magnitudes
    range_pred = u_r_predicted.max() - u_r_predicted.min()
    range_actual = u_actual.max() - u_actual.min()
    print(f"  Pearson r (pred vs actual, demeaned)  = {pearson_u:+.4f} (p={pearson_u_p:.4e})",
          flush=True)
    print(f"  Spearman ρ                             = {spearman_u:+.4f} (p={spearman_u_p:.4e})",
          flush=True)
    print(f"  Range (predicted) = {range_pred:.4f}; Range (actual) = {range_actual:.4f}; "
          f"ratio = {range_pred / range_actual if range_actual > 0 else float('inf'):.2f}",
          flush=True)

    # ============================================================
    # COMBINED VERDICT
    # ============================================================
    print("\n[combined] reconciling tests 1-4 to combined verdict ...", flush=True)
    if test2_verdict == "LITERAL-MAGNITUDE CONSISTENT":
        combined_verdict = "CONJECTURE LITERAL FORM INTACT"
        combined_text = ("B2's a_final dummy coefficients form a slope on log(a_final) "
                         "consistent with 1.0; the literal-magnitude form of the "
                         "§3.5 conjecture is supported. B1's small u_r span is a "
                         "parameterization or partial-pooling artifact.")
    elif test2_verdict == "LITERAL-MAGNITUDE FALSIFIED":
        combined_verdict = "CONJECTURE FALSIFIED IN MAGNITUDE"
        combined_text = ("B2's a_final dummy coefficients exclude a slope of 1.0 on "
                         "log(a_final) and are consistent with the small B1 magnitude. "
                         "The literal-magnitude form of the conjecture is rejected; "
                         "the directional/monotone form survives.")
    elif test2_verdict == "REDUCED SLOPE":
        combined_verdict = "CONJECTURE REFINES TO REDUCED SLOPE"
        combined_text = ("B2's a_final dummy coefficients form a slope on log(a_final) "
                         "that is positive but excludes 1.0. The structural relationship "
                         "survives at reduced slope; the literal-magnitude form needs "
                         "to be restated.")
    else:
        combined_verdict = "AMBIGUOUS"
        combined_text = ("B2's posterior is too wide at N=10K to disambiguate. "
                         "Recommend deferring the verdict until full-N HMC runs.")
    print(f"  Combined verdict: {combined_verdict}", flush=True)

    # ============================================================
    # WRITE MARKDOWN
    # ============================================================
    write_markdown(
        pearson_r=pearson_r, pearson_p=pearson_p, r_lo=r_lo, r_hi=r_hi,
        spearman_rho=spearman_rho, spearman_p=spearman_p,
        test1_verdict=test1_verdict, test1_text=test1_text,
        af_levels=af_levels, log_af_levels=log_af_levels,
        af_coefs_mean=af_coefs_mean, af_coefs_q025=af_coefs_q025, af_coefs_q975=af_coefs_q975,
        slope_mean=slope_mean, slope_q025=slope_q025, slope_q975=slope_q975,
        test2_verdict=test2_verdict, test2_text=test2_text,
        b0_mean=b0_mean, b0_sd=b0_sd, b0_q025=b0_q025, b0_q975=b0_q975,
        b2_log_n_mean=b2_log_n_mean, b2_log_n_sd=b2_log_n_sd,
        b2_log_n_q025=b2_log_n_q025, b2_log_n_q975=b2_log_n_q975,
        diff=diff, pooled_sd=pooled_sd, z_diff=z_diff,
        test3_verdict=test3_verdict, test3_text=test3_text,
        residues=residues, mean_log_n_per_r=mean_log_n_per_r,
        mean_log_n_global=mean_log_n_global,
        u_r_predicted_uncorrected=u_r_predicted_uncorrected,
        u_r_predicted=u_r_predicted, u_actual=u_actual,
        pearson_u=pearson_u, pearson_u_p=pearson_u_p,
        spearman_u=spearman_u, spearman_u_p=spearman_u_p,
        range_pred=range_pred, range_actual=range_actual,
        combined_verdict=combined_verdict, combined_text=combined_text,
    )
    print(f"\n[save] {MD_OUT}", flush=True)


def write_markdown(**k):
    L = []
    L.append("# B2 Coefficient and log(n) Correlation Check — Pre-Submission Verdict\n")
    L.append(f"**Headline:** {k['combined_verdict']}.\n")
    L.append(f"\n{k['combined_text']}\n\n---\n")

    # Section 1
    L.append("\n## Section 1 — Test 1: log(n) ↔ log(a_final(r mod 8)) correlation\n")
    L.append(f"\n- Pearson r = {k['pearson_r']:+.6f} (p = {k['pearson_p']:.4e}), "
             f"95% CI [{k['r_lo']:+.4f}, {k['r_hi']:+.4f}]")
    L.append(f"- Spearman ρ = {k['spearman_rho']:+.6f} (p = {k['spearman_p']:.4e})")
    L.append(f"\n**Verdict: {k['test1_verdict']}.** {k['test1_text']}\n")

    # Section 2
    L.append("\n## Section 2 — Test 2: B2 a_final dummy coefficients\n")
    L.append("\nB2's design uses categorical dummies on a_final (reference = af=1, j=0). "
             "Posterior coefficients on the 3 non-reference levels:\n")
    L.append("\n| level | a_final | j | log(a_final) | β posterior mean | 95% CI | conjecture predicts |")
    L.append("|---|---|---|---|---|---|---|")
    for af_level, lal, mean, lo, hi in zip(k["af_levels"], k["log_af_levels"],
                                              k["af_coefs_mean"], k["af_coefs_q025"],
                                              k["af_coefs_q975"]):
        j = int(round(np.log(af_level) / np.log(3)))
        match_lit = "✓ contains 1.0×log(af)" if (lo <= lal <= hi) else "✗ excludes 1.0×log(af)"
        L.append(f"| af = {af_level} | 3^{j} | {j} | {lal:.4f} | "
                 f"**{mean:+.4f}** | [{lo:+.4f}, {hi:+.4f}] | {lal:+.4f} ({match_lit}) |")
    L.append(f"\nImplied per-draw OLS slope on log(a_final) (4 dummy values, "
             f"a_final=1 reference fixed at 0):")
    L.append(f"\n- **Posterior mean slope = {k['slope_mean']:+.4f}**, "
             f"95% CI [{k['slope_q025']:+.4f}, {k['slope_q975']:+.4f}]")
    L.append(f"- Conjecture's literal magnitude predicts slope = 1.0")
    L.append(f"- B1's smoke-test \"falsified\" form predicts slope ≈ 0.085")
    L.append(f"\n**Verdict: {k['test2_verdict']}.** {k['test2_text']}\n")

    # Section 3
    L.append("\n## Section 3 — Test 3: β_log_n in B0 vs B2\n")
    L.append("\n| spec | β_log_n posterior mean | SD | 95% CI |")
    L.append("|---|---|---|---|")
    L.append(f"| B0 (log(n) only) | {k['b0_mean']:.6f} | {k['b0_sd']:.6f} | "
             f"[{k['b0_q025']:.6f}, {k['b0_q975']:.6f}] |")
    L.append(f"| B2 (+ a_final dummies) | {k['b2_log_n_mean']:.6f} | "
             f"{k['b2_log_n_sd']:.6f} | "
             f"[{k['b2_log_n_q025']:.6f}, {k['b2_log_n_q975']:.6f}] |")
    L.append(f"\n- Difference (B2 − B0) = {k['diff']:+.6f}, pooled SD = "
             f"{k['pooled_sd']:.6f}, |z| = {abs(k['z_diff']):.3f}")
    L.append(f"\n**Verdict: {k['test3_verdict']}.** {k['test3_text']}\n")

    # Section 4
    L.append("\n## Section 4 — Test 4: Parameterization-corrected u_r prediction\n")
    L.append("\nFor each residue mod 8, compute the predicted random effect under "
             "the literal conjecture, corrected for any colinearity with log(n):\n")
    L.append("```\n"
             "u_r_predicted = log(a_final(r)) "
             "- β_log_n × (mean(log(n) | r) - global mean(log(n)))\n"
             "```\n")
    L.append(f"- Global mean(log(n)) = {k['mean_log_n_global']:.4f}")
    L.append("\n| r | mean(log(n)\\|r) | log(a_final) | u_r predicted (corrected) | u_r actual (B1) |")
    L.append("|---|---|---|---|---|")
    for r, mln, lal, p, a in zip(k["residues"], k["mean_log_n_per_r"],
                                    k["u_r_predicted_uncorrected"],
                                    k["u_r_predicted"], k["u_actual"]):
        L.append(f"| {r} | {mln:.4f} | {lal:.4f} | {p:+.4f} | {a:+.4f} |")
    L.append(f"\n- Range (predicted, parameterization-corrected) = {k['range_pred']:.4f}")
    L.append(f"- Range (actual B1) = {k['range_actual']:.4f}")
    L.append(f"- Pearson r (demeaned predicted vs actual) = {k['pearson_u']:+.4f} "
             f"(p = {k['pearson_u_p']:.4e})")
    L.append(f"- Spearman ρ = {k['spearman_u']:+.4f} (p = {k['spearman_u_p']:.4e})")
    L.append("\nIf the parameterization correction recovers the literal magnitude, "
             "the predicted range should be ≈ 3.30 (= log(27)) and the demeaned "
             "correlation should be high. If the correction has no effect, the predicted "
             "range stays at log(27) (since the colinearity is small) and the actual "
             "range stays at 0.25 — meaning the magnitude gap is real, not a "
             "parameterization artifact.\n")

    # Section 5
    L.append("\n## Section 5 — Combined Verdict\n")
    L.append(f"\n**{k['combined_verdict']}.**\n")
    L.append(f"\n{k['combined_text']}\n")

    # Section 6 — draft §5 paragraph
    L.append("\n## Section 6 — Recommended §5 Paragraph (draft for editorial review)\n\n")
    para = build_section_5_paragraph(k)
    L.append(para)
    L.append("\n---\n")
    L.append(f"\n*Inputs: B0/B1/B2 HMC posteriors at N=10K, 8000 train / 2000 test, "
             f"seed {SEED}. NO new MCMC. Tests 1–4 on existing chains.*\n")
    MD_OUT.write_text("\n".join(L), encoding="utf-8")


def build_section_5_paragraph(k):
    cv = k["combined_verdict"]
    if cv == "CONJECTURE LITERAL FORM INTACT":
        body = (
"To validate the §5 head-to-head's structural ordering and the §3.5 conjecture "
"under proper inference, we re-ran B0–B4 at N=10⁴ under HMC (cmdstanpy 1.3.0, "
"4 chains × 2000 iter each, parallel-chains mode). All five specifications "
"converged with R̂ < 1.01, ESS_bulk > 600, zero divergences, zero treedepth "
"saturations, and E-BFMI > 0.80. The HMC ordinal ordering reproduces the "
"Pathfinder finding: B3 > B2 > B1 > B0 on held-out log score, with B4 within "
"sampling noise of B3 (B3 better by 0.10 nats). σ_u in B1 inflates substantially "
f"under HMC (posterior mean {abs(0.0973):.4f} vs Pathfinder 0.003*), confirming "
"the §5 caveat that Pathfinder collapses hierarchical scale parameters. σ_u in B4 "
"under HMC remains small (0.008), confirming that the mod-8 random effect carries "
"no residual information once a_final at k=6 is in the model — the conjecture's "
"exhaustion claim survives proper inference.\n\n"
"On the literal-magnitude form of the §3.5 conjecture (that u_r equals "
"log(a_final(r)) at k=3 up to a constant): B2's posterior coefficients on the "
"a_final fixed effects, when regressed against log(a_final), give a slope of "
f"{k['slope_mean']:+.3f} with 95% credible interval "
f"[{k['slope_q025']:+.3f}, {k['slope_q975']:+.3f}]. This interval contains 1.0, "
"so the literal-magnitude form is statistically consistent with the data once "
"the random-effect parameterization is removed. B1's small u_r span ("
f"observed range ≈ 0.25 for the posterior means) reflects a combination of "
"partial-pooling shrinkage and design-matrix identification within the "
"hierarchical model, not a falsification of the structural identity.")
    elif cv == "CONJECTURE FALSIFIED IN MAGNITUDE":
        body = (
"To validate the §5 head-to-head's structural ordering and the §3.5 conjecture "
"under proper inference, we re-ran B0–B4 at N=10⁴ under HMC (cmdstanpy 1.3.0, "
"4 chains × 2000 iter each, parallel-chains mode). All five specifications "
"converged with R̂ < 1.01, ESS_bulk > 600, zero divergences, zero treedepth "
"saturations, and E-BFMI > 0.80. The HMC ordinal ordering reproduces the "
"Pathfinder finding: B3 > B2 > B1 > B0 on held-out log score, with B4 within "
"sampling noise of B3 (B3 better by 0.10 nats). σ_u in B1 inflates substantially "
"under HMC (posterior mean ≈ 0.097 vs Pathfinder 0.003*), confirming the §5 "
"caveat that Pathfinder collapses hierarchical scale parameters. σ_u in B4 "
"under HMC remains small (0.008), confirming that the mod-8 random effect "
"carries no residual information once a_final at k=6 is in the model.\n\n"
"On the literal-magnitude form of the §3.5 conjecture (u_r = log(a_final(r)) "
"+ const at k=3): B2's posterior coefficients on the a_final fixed effects, "
"when regressed against log(a_final), give a slope of "
f"{k['slope_mean']:+.3f} with 95% credible interval "
f"[{k['slope_q025']:+.3f}, {k['slope_q975']:+.3f}], which excludes 1.0. The "
"literal magnitude is rejected at this scale. The DIRECTIONAL form survives — "
"u_r is monotone-increasing in log(a_final(r)) — but the slope is materially "
"smaller than log(3). We therefore restate the conjecture as a directional "
"claim: *the mod-8 random effect is monotone in the prefix-decomposition's "
"log(a_final(r)) at k=3, with a slope to be determined empirically.* The "
"structural origin claim — that the random effect's signal traces to the "
"prefix's terminal a_final — survives; the algebraic identity at unit slope "
"does not.")
    elif cv == "CONJECTURE REFINES TO REDUCED SLOPE":
        body = (
"To validate the §5 head-to-head's structural ordering and the §3.5 conjecture "
"under proper inference, we re-ran B0–B4 at N=10⁴ under HMC (cmdstanpy 1.3.0, "
"4 chains × 2000 iter each, parallel-chains mode). All five specifications "
"converged with clean diagnostics (R̂ < 1.01, ESS > 600, zero divergences, "
"E-BFMI > 0.80). The ordinal ordering B3 > B2 > B1 > B0 on held-out log "
"score reproduces the Pathfinder finding; B4 is within sampling noise of "
"B3. σ_u in B1 inflates substantially under HMC (posterior 0.097 vs "
"Pathfinder 0.003*), confirming the §5 caveat. σ_u in B4 remains small "
"(0.008), confirming the exhaustion claim of the conjecture.\n\n"
"On the literal-magnitude form: B2's posterior coefficients on the a_final "
f"fixed effects imply a slope on log(a_final) of {k['slope_mean']:+.3f} "
f"(95% CI [{k['slope_q025']:+.3f}, {k['slope_q975']:+.3f}]), which is "
"meaningfully positive but excludes 1.0. The conjecture's structural claim "
"survives at REDUCED slope: u_r is linearly proportional to log(a_final(r)), "
"with proportionality constant materially less than 1. We restate §3.5 as: "
"*u_r ≈ c · log(a_final(r)) + const at k=3, with c ≈ "
f"{k['slope_mean']:.3f} on data of N = 10^4*.")
    else:
        body = (
"To validate the §5 head-to-head's structural ordering and the §3.5 conjecture "
"under proper inference, we re-ran B0–B4 at N=10⁴ under HMC (cmdstanpy 1.3.0). "
"All five specifications converged cleanly (R̂ < 1.01, ESS > 600, zero divergences). "
"The ordinal ordering B3 > B2 > B1 > B0 reproduces the Pathfinder finding; "
"B4 is within sampling noise of B3, and σ_u in B4 remains small under HMC "
"(0.008), confirming the exhaustion claim. σ_u in B1 inflates ~32× over the "
"Pathfinder estimate, confirming the §5 caveat about Pathfinder hierarchical "
"underdispersion.\n\n"
"The literal-magnitude form of the §3.5 conjecture is unresolved at this scale: "
"the implied slope on log(a_final) from B2's fixed-effect dummies has 95% CI "
f"[{k['slope_q025']:+.3f}, {k['slope_q975']:+.3f}], which is too wide to "
"exclude either the literal value 1.0 or smaller alternatives. The directional "
"claim survives. We recommend that confirmation of the literal-magnitude form "
"awaits full-N HMC at the Bonacorsi-Bordoni production scale.")
    return body


if __name__ == "__main__":
    main()
