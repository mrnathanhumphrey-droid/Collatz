"""prefix_steps Confound Check — Residual Gap Analysis for B2 Slope.

Per the locked brief: tests whether the 14% residual gap between empirical
B2 slope (0.0774) and Pred A (K_h / E[σ] = 0.0675) is explained by the
prefix_steps(r) component of α_det = prefix_steps(r) + K_h · log(a_final(r)/2^k).

NO new MCMC. Restricted to 4 odd residues r ∈ {1, 3, 5, 7} per the brief.
"""
from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

K_h = 3 / np.log(4 / 3)
E_SIGMA = 154.4339
k = 3
MD_OUT = Path("C:/Collatz/docs/paper1_prefix_steps_confound_check.md")
MD_OUT.parent.mkdir(parents=True, exist_ok=True)


def prefix(r, a0, max_steps=400):
    a, c, s = a0, r, 0
    while a % 2 == 0 and s < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3 * c + 1
        s += 1
    return a, s


def main():
    odd_residues = [1, 3, 5, 7]
    rows = []
    for r in odd_residues:
        af, steps = prefix(r, 1 << k)
        log_af = float(np.log(af))
        alpha_det = steps + K_h * np.log(af / (1 << k))
        j = int(round(np.log(af) / np.log(3)))
        rows.append({"r": r, "a_final": af, "j": j, "prefix_steps": steps,
                     "log_a_final": log_af, "alpha_det": alpha_det})

    print("[Section 1] Per-residue triples (k=3, K_h = 3/log(4/3) ≈ 10.4282):")
    print(f"  {'r':>2} {'a_final':>8} {'j':>2} {'prefix_steps':>13} {'log(a_final)':>13} {'α_det':>10}")
    for row in rows:
        print(f"  {row['r']:>2} {row['a_final']:>8} {row['j']:>2} {row['prefix_steps']:>13} "
              f"{row['log_a_final']:>13.4f} {row['alpha_det']:>10.4f}")

    prefix_steps_arr = np.array([row["prefix_steps"] for row in rows], dtype=float)
    log_af_arr = np.array([row["log_a_final"] for row in rows], dtype=float)
    alpha_det_arr = np.array([row["alpha_det"] for row in rows], dtype=float)

    # Section 2: correlation
    pearson_r = float(np.corrcoef(prefix_steps_arr, log_af_arr)[0, 1])
    print(f"\n[Section 2] Pearson correlation prefix_steps vs log(a_final): "
          f"r = {pearson_r:+.6f}")

    # Section 3: Predictions
    # B2 dummies are per-a_final-level (a_final ∈ {3, 9, 27}, reference = 1).
    # For odd-residue predictions: each af class maps to specific odd residues.
    # The brief asks for Pred A and Pred B at the dummy level.
    af_levels = [3, 9, 27]
    log_af_levels = np.log(af_levels)
    observed = np.array([0.1037, 0.1806, 0.2578])

    # Pred A: β · log(a_final) with β = K_h / E[σ]
    beta_A = K_h / E_SIGMA
    pred_A = beta_A * log_af_levels

    # Pred B: α_det(r) / E[σ], aggregated within a_final class.
    # For odd residues only:
    #   af=3: r=5, α_det = 4 + K_h·log(3/8)
    #   af=9: r=1,3 (both have prefix_steps=5, α_det = 5 + K_h·log(9/8))
    #   af=27: r=7, α_det = 6 + K_h·log(27/8)
    # The dummy is a CLASS effect relative to af=1 reference. For odd-only
    # restriction the af=1 class has no odd residues. Per the brief's Pred B,
    # use α_det of the class (mean over odd residues mapping to that class).
    # Since dummies are RELATIVE to af=1, also subtract α_det(af=1) = α_det(r=0)
    # = 3 + K_h·log(1/8) = 3 - K_h·k·log(2).
    af_to_alpha_det = {}
    for row in rows:
        af_to_alpha_det.setdefault(row["a_final"], []).append(row["alpha_det"])
    af_to_alpha_det_mean = {af: float(np.mean(v)) for af, v in af_to_alpha_det.items()}
    # Reference α_det for a_final=1 (only r=0, even, but explicit per-formula):
    af1_steps = 3  # (8,0)→(4,0)→(2,0)→(1,0)
    alpha_det_af1 = af1_steps + K_h * np.log(1 / (1 << k))  # = 3 - K_h·k·log(2)

    pred_B = np.array([(af_to_alpha_det_mean[L] - alpha_det_af1) / E_SIGMA
                       for L in af_levels])

    print("\n[Section 3] Predictions and observations:")
    print(f"  β_A = K_h / E[σ]               = {beta_A:.6f}")
    print(f"  α_det(af=1) reference          = {alpha_det_af1:.6f}")
    print(f"  α_det per af class (mean over odd residues):")
    for L in af_levels:
        print(f"    af={L:>2d}: α_det = {af_to_alpha_det_mean[L]:.4f}")
    print()
    print(f"  {'af':>3} {'log(af)':>9} {'Pred A':>9} {'Pred B':>9} "
          f"{'Observed':>10} {'A ratio':>8} {'B ratio':>8}")
    print(f"  {'-'*3} {'-'*9} {'-'*9} {'-'*9} {'-'*10} {'-'*8} {'-'*8}")
    A_ratios = []; B_ratios = []
    for L, lal, pA, pB, obs in zip(af_levels, log_af_levels, pred_A, pred_B, observed):
        rA = obs / pA; rB = obs / pB
        A_ratios.append(rA); B_ratios.append(rB)
        print(f"  {L:>3d} {lal:>9.4f} {pA:>9.4f} {pB:>9.4f} "
              f"{obs:>10.4f} {rA:>8.4f} {rB:>8.4f}")

    mae_A = float(np.mean(np.abs(observed - pred_A)))
    mae_B = float(np.mean(np.abs(observed - pred_B)))
    print(f"\n  MAE Pred A = {mae_A:.6f}")
    print(f"  MAE Pred B = {mae_B:.6f}")

    # Slope-level comparison
    # Empirical: 0.0774 [0.0676, 0.0869]
    # Pred A slope: K_h / E[σ] = 0.0675
    # Pred B slope: regress α_det on log(a_final) across the residues
    # Across all residues prefix_steps = j + k (linear in log(af) with slope 1/log(3))
    # so slope of α_det on log(af) = 1/log(3) + K_h
    slope_A = K_h / E_SIGMA
    slope_B = (1 / np.log(3) + K_h) / E_SIGMA
    print(f"\n  Slope Pred A on log(a_final): {slope_A:.6f}")
    print(f"  Slope Pred B on log(a_final): {slope_B:.6f}")
    print(f"  Empirical slope:              0.0774 (CI [0.0676, 0.0869])")

    # Section 4: Verdict
    A_dist = [abs(r - 1.0) for r in A_ratios]
    B_dist = [abs(r - 1.0) for r in B_ratios]
    pp_diff = [(a - b) * 100 for a, b in zip(A_dist, B_dist)]
    print(f"\n[Section 4] Verdict criterion:")
    print(f"  Per-residue (|A_ratio - 1| - |B_ratio - 1|) × 100, in pp:")
    for L, d in zip(af_levels, pp_diff):
        print(f"    af={L:>2d}: {d:+.2f} pp ({'Pred B closer' if d > 0 else 'Pred A closer'})")
    all_B_closer_5pp = all(d >= 5 for d in pp_diff)
    all_within_5pp = all(abs(d) < 5 for d in pp_diff)
    Pred_B_worse = all(d < 0 for d in pp_diff)
    if Pred_B_worse:
        verdict = "PRED B WORSE — HALT"
        verdict_text = ("Pred B is meaningfully worse than Pred A. Something is "
                        "wrong with the α_det formula or the reasoning. Halt.")
    elif all_B_closer_5pp:
        verdict = "PREFIX_STEPS CONFOUND DOMINATES"
        verdict_text = ("Pred B's per-class ratios are closer to 1.0 than Pred A's "
                        "by ≥ 5 percentage points on all three classes. The prefix_steps "
                        "term is the dominant residual source. The conjecture should "
                        "reference α_det as a whole, not log(a_final) alone.")
    elif all_within_5pp:
        verdict = "WITHIN LINEARIZATION NOISE"
        verdict_text = ("Pred A and Pred B agree within 5 percentage points on all "
                        "three classes. prefix_steps doesn't materially shift the "
                        "prediction; the 14% residual is linearization-approximation "
                        "+ finite-sample noise. The conjecture as currently stated "
                        "(log(a_final) only) is the right scope.")
    else:
        verdict = "MIXED"
        verdict_text = ("Per-class differences mixed (some ≥5pp, some <5pp). "
                        "Mixed evidence; document and inspect manually.")
    print(f"  -> {verdict}")
    print(f"  {verdict_text}")

    write_markdown(rows, pearson_r, af_levels, log_af_levels, pred_A, pred_B,
                   observed, A_ratios, B_ratios, mae_A, mae_B,
                   af_to_alpha_det_mean, alpha_det_af1, slope_A, slope_B,
                   pp_diff, verdict, verdict_text)
    print(f"\n[save] {MD_OUT}")


def write_markdown(rows, pearson_r, af_levels, log_af_levels, pred_A, pred_B,
                   observed, A_ratios, B_ratios, mae_A, mae_B,
                   af_to_alpha_det_mean, alpha_det_af1, slope_A, slope_B,
                   pp_diff, verdict, verdict_text):
    L = []
    L.append("# prefix_steps Confound Check — Residual Gap Analysis for B2 Slope\n")
    L.append(f"\n**Headline:** {verdict}.\n")
    L.append(f"\n{verdict_text}\n")
    L.append("\n---\n")

    L.append("\n## Section 1 — Four (r, prefix_steps, log(a_final), α_det) triples\n")
    L.append(f"\nUsing K_h = 3/log(4/3) ≈ {K_h:.6f}, k=3.")
    L.append(f"\nα_det(r) = prefix_steps(r) + K_h · log(a_final(r) / 2^k).\n")
    L.append("\n| r (odd) | a_final | j | prefix_steps | log(a_final) | α_det |")
    L.append("|---|---|---|---|---|---|")
    for row in rows:
        L.append(f"| {row['r']} | {row['a_final']} | {row['j']} | "
                 f"{row['prefix_steps']} | {row['log_a_final']:.4f} | "
                 f"{row['alpha_det']:+.4f} |")

    L.append("\n## Section 2 — Pearson correlation between prefix_steps and log(a_final)\n")
    L.append(f"\nAcross the 4 odd residues:\n")
    L.append(f"\n- **Pearson r = {pearson_r:+.6f}**\n")
    L.append("\nReason: prefix_steps(r) = j(r) + k where a_final(r) = 3^j(r). At fixed k, "
             "prefix_steps is linear in j(r), and log(a_final(r)) = j(r) · log(3). So "
             "prefix_steps and log(a_final) are *exactly* linearly related (slope = 1/log(3) "
             "≈ 0.910), giving Pearson r = ±1.0 across any residue set whose distinct "
             "a_final values span more than one j level.\n")

    L.append("\n## Section 3 — Predictions vs observations\n")
    L.append(f"\nE[σ] (training set) = {E_SIGMA:.4f}.\n")
    L.append(f"\n- **Pred A** (literal conjecture, properly link-translated): "
             f"dummy = β · log(a_final) with β = K_h / E[σ] = {K_h / E_SIGMA:.6f}.")
    L.append(f"\n- **Pred B** (full α_det rescaled): "
             f"dummy = [α_det(class) − α_det(af=1)] / E[σ].\n")
    L.append("\n### α_det per class (mean over odd residues mapping to that class)\n")
    L.append("\n| a_final | α_det (mean over odd r) |")
    L.append("|---|---|")
    for af in af_levels:
        L.append(f"| {af} | {af_to_alpha_det_mean[af]:+.4f} |")
    L.append(f"| 1 (reference, even r=0) | {alpha_det_af1:+.4f} |")

    L.append("\n### Predictions table\n")
    L.append("\n| a_final | log(a_final) | Pred A | Pred B | Observed | A ratio | B ratio |")
    L.append("|---|---|---|---|---|---|---|")
    for af, lal, pA, pB, obs, rA, rB in zip(af_levels, log_af_levels, pred_A, pred_B,
                                              observed, A_ratios, B_ratios):
        L.append(f"| {af} | {lal:.4f} | {pA:.4f} | {pB:.4f} | "
                 f"{obs:.4f} | {rA:.4f} | {rB:.4f} |")
    L.append(f"\n- MAE Pred A = {mae_A:.6f}")
    L.append(f"- MAE Pred B = {mae_B:.6f}\n")

    L.append("\n### Slope-level comparison (across all 4 dummy points incl af=1 ref at 0)\n")
    L.append(f"\n- Pred A slope on log(a_final) = K_h / E[σ] = {slope_A:.6f}")
    L.append(f"- Pred B slope on log(a_final) = (1/log(3) + K_h) / E[σ] = {slope_B:.6f}")
    L.append(f"- Empirical OLS slope across 4 dummy points = 0.0774 (95% CI [0.0676, 0.0869])\n")

    L.append("\n## Section 4 — Verdict\n")
    L.append(f"\n**{verdict}.**\n")
    L.append(f"\n{verdict_text}\n")
    L.append("\nPer-class (|A_ratio − 1| − |B_ratio − 1|) × 100, in percentage points:\n")
    L.append("\n| a_final | pp difference (positive = Pred B closer) |")
    L.append("|---|---|")
    for af, d in zip(af_levels, pp_diff):
        L.append(f"| {af} | {d:+.2f} |")

    L.append("\n## Section 5 — Implication for §3.5 conjecture statement\n")
    if verdict == "PREFIX_STEPS CONFOUND DOMINATES":
        L.append("\nThe conjecture in §3.5 should reference **α_det(r)** as a whole, not "
                 "log(a_final(r)) alone. The prefix_steps(r) term contributes a "
                 "non-trivial component to the per-class effect on σ that the empirical "
                 "B2 dummies absorb. Stating the conjecture purely on log(a_final) "
                 "leaves a systematic residual that matches the prefix_steps slope "
                 "1 / (E[σ] · log(3)).\n")
        L.append("\nRecommended restatement:\n")
        L.append("\n*The Bonacorsi-Bordoni mod-8 random effect is, up to a global "
                 "additive constant, α_det(r) / E[σ] at k=3, where α_det(r) = "
                 "prefix_steps(r) + K_h · log(a_final(r)/2^k) is the prefix "
                 "decomposition's per-class structural offset on σ-scale.*\n")
    elif verdict == "WITHIN LINEARIZATION NOISE":
        L.append("\nPred A and Pred B agree within 5 percentage points across all three "
                 "classes. The prefix_steps term does not materially shift the "
                 "prediction at the precision the B2 dummies achieve at N=10⁴. The "
                 "§3.5 conjecture as currently stated (log(a_final) only) is the right "
                 "scope; the residual 14% gap in the slope comparison is "
                 "linearization-approximation noise plus finite-sample posterior "
                 "uncertainty, not a missing model term.\n")
    elif verdict == "PRED B WORSE — HALT":
        L.append("\n**Halt and report**: Pred B is worse than Pred A. Something is "
                 "wrong with the α_det formula or the derivation. Do not proceed with "
                 "§5 changes until the discrepancy is understood.\n")
    else:  # MIXED
        L.append("\nMixed evidence across the three classes. The verdict is not "
                 "decisive at this precision; document the per-class differences and "
                 "consider whether full-N HMC at the Bonacorsi-Bordoni production "
                 "scale would resolve the ambiguity.\n")

    L.append("\n## Section 6 — Recommended §5 paragraph language\n")
    if verdict == "PREFIX_STEPS CONFOUND DOMINATES":
        L.append(
"\nTo validate the §5 head-to-head's structural ordering and the §3.5 "
"conjecture under proper inference, we re-ran B0–B4 at N=10⁴ under HMC "
"(cmdstanpy 1.3.0, 4 chains × 2000 iter each, parallel-chains mode). All "
"five specifications converged with R̂ < 1.01, ESS_bulk > 600, zero "
"divergences, zero treedepth saturations, and E-BFMI > 0.80. The HMC ordinal "
"ordering reproduces the Pathfinder finding: B3 > B2 > B1 > B0 on held-out "
"log score, with B4 within sampling noise of B3 (B3 better by 0.10 nats). "
"σ_u in B1 inflates substantially under HMC (posterior mean 0.097 vs "
"Pathfinder 0.003*), confirming the §5 caveat that Pathfinder collapses "
"hierarchical scale parameters in non-Gaussian regions. σ_u in B4 under HMC "
"remains small (0.008), confirming the conjecture's exhaustion claim — once "
"a_final at k=6 is in the model, the mod-8 random effect carries no residual "
"information.\n\n"
"On the literal-magnitude form of the §3.5 conjecture: the empirical OLS "
"slope of B2's a_final fixed-effect dummies on log(a_final(r)) is 0.0774 "
"with 95% credible interval [0.0676, 0.0869]. Translating the §3.5 "
"conjecture through the NB2 log link, the predicted GLM-scale slope on the "
"full structural offset α_det(r) = prefix_steps(r) + K_h · log(a_final(r)/2^k) "
f"is (1/log(3) + K_h)/E[σ] = {slope_B:.4f}, against (K_h/E[σ]) = "
f"{slope_A:.4f} on log(a_final) alone. The empirical slope falls between "
"these two predictions, with the α_det form 5.4% below empirical against "
"the log(a_final)-only form 14.7% below. We therefore restate the "
"conjecture in §3.5 to reference α_det(r) as a whole rather than "
"log(a_final(r)) alone, since prefix_steps(r) contributes a non-negligible "
"component to the per-class effect that the empirical B2 dummies absorb at "
"this precision."
        )
    elif verdict == "WITHIN LINEARIZATION NOISE":
        L.append(
"\nTo validate the §5 head-to-head's structural ordering and the §3.5 "
"conjecture under proper inference, we re-ran B0–B4 at N=10⁴ under HMC. "
"All five specifications converged cleanly. The ordinal ordering B3 > B2 > "
"B1 > B0 reproduces the Pathfinder finding; σ_u in B1 inflates under HMC "
"(0.097 vs Pathfinder 0.003*) and σ_u in B4 remains small (0.008), "
"confirming the §5 caveat about Pathfinder hierarchical underdispersion "
"and the conjecture's exhaustion claim. The empirical OLS slope of B2's "
"a_final dummies on log(a_final(r)) is 0.0774 (95% CI [0.0676, 0.0869]); "
f"the §3.5 conjecture's literal form, properly link-translated, predicts "
f"K_h/E[σ] = {slope_A:.4f}, which sits at the lower edge of the empirical "
"95% CI. The remaining 14% gap is at the level of the linearization "
"approximation (the global-E[σ] denominator is an approximation to the "
"per-class-conditional E[σ|class]), so the conjecture's literal form is "
"consistent with the data at this scale."
        )
    elif verdict == "PRED B WORSE — HALT":
        L.append("\n*Section 5 paragraph deferred: verdict halt.*\n")
    else:
        L.append(
"\nThe empirical OLS slope of B2's a_final dummies on log(a_final(r)) is "
f"0.0774 (95% CI [0.0676, 0.0869]). The §3.5 conjecture's literal form, "
f"properly link-translated, predicts K_h/E[σ] ≈ {slope_A:.4f}; the full "
f"α_det form predicts {slope_B:.4f}. Both lie at or near the lower end of "
"the empirical CI; the data at N=10⁴ does not decisively favor one over "
"the other. We defer the magnitude-form question to full-N HMC at the "
"Bonacorsi-Bordoni production scale."
        )
    L.append("\n\n---\n")
    L.append("\n*Inputs: B2 HMC posterior at N=10K, 8000 train / 2000 test, seed 20260509. "
             "NO new MCMC. Restricted to 4 odd residues r ∈ {1, 3, 5, 7} per the brief.*\n")
    MD_OUT.write_text("\n".join(L), encoding="utf-8")


if __name__ == "__main__":
    main()
