"""Generate one-sheet PDF for Prof. Lin (Bayesian stats, Maryland)."""
from fpdf import FPDF
from pathlib import Path

OUT = Path(r"C:\Users\Nate\Desktop\collatz_one_sheet_lin_2026-05-01.pdf")

pdf = FPDF(orientation="portrait", unit="mm", format="letter")
pdf.set_auto_page_break(auto=True, margin=12)
pdf.add_page()
pdf.set_margins(left=15, top=12, right=15)

# Title
pdf.set_font("Helvetica", "B", 13)
pdf.cell(0, 6, "Prefix Decomposition for Collatz Trajectory Functionals",
         new_x="LMARGIN", new_y="NEXT")
pdf.set_font("Helvetica", "", 9)
pdf.cell(0, 4, "Empirical note. Nathan Humphrey. May 2026.",
         new_x="LMARGIN", new_y="NEXT")
pdf.ln(2)

def section(title):
    pdf.set_font("Helvetica", "B", 9.5)
    pdf.cell(0, 4.5, title, new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 8.5)

def body(text):
    pdf.multi_cell(0, 3.6, text, new_x="LMARGIN", new_y="NEXT")
    pdf.ln(0.5)

# --- Setup ---
section("1. Setup")
body(
    "sigma(n) denotes the Collatz total stopping time. The standard random-walk "
    "heuristic predicts E[sigma | n] ~ K log n, with K = 2/(ln 4 - ln 3) ~ 6.9521 "
    "for all n and K = 3/(ln 4 - ln 3) ~ 10.4282 for odd n. Bonacorsi & Bordoni "
    "(arXiv:2603.04479, Mar 2026) fit this distribution with a Bayesian NB2 GLM "
    "(log link), adding hierarchical random effects on the residue class n mod 2^k."
)

# --- Prior art ---
section("2. Prior art (Terras, 1976)")
body(
    "Terras Lemma 3: the parity vector P_k(n) of the first k Collatz iterates "
    "depends only on n mod 2^k. Lemma 4 (asymptotic): after k steps, "
    "S_k ~ S_0 * 3^{d(k)} * 2^{-k}, where d(k) is the count of odd steps among "
    "the first k iterations. This is the asymptotic version of the prefix "
    "decomposition used here, with the additive correction dropped."
)

# --- Refinement ---
section("3. Refinement (this work)")
body(
    "Track the symbolic state (a, c) such that state = a*m + c, where "
    "m = (n - r)/2^k is the integer 'tail' of n. Apply Collatz forward, with "
    "branch parity forced by c, until a becomes odd (variable prefix length, "
    "generally < k, deterministic per residue class). Result: post-prefix state "
    "= a_final * m + c_final, where a_final in {3^j : 1 <= j <= k}, and c_final "
    "is an integer offset retained explicitly. Differs from Terras Lemma 4 in "
    "(i) the stopping rule (first odd a, not fixed k), and (ii) retention of c_final."
)

# --- Empirical check ---
section("4. Empirical check (N = 10^7, k = 6, 32 odd classes mod 64)")
body(
    "Per-class OLS of sigma vs log(n): mean slope mu_beta = 10.4014 "
    "(odd-only heuristic 10.4282; 0.26% relative deviation). Between-class slope "
    "variance Var(beta) = 0.0118, mean within-class SE^2 = 0.0232; "
    "moment-corrected tau_beta^2 = max(0, Var - mean SE^2) = 0. Slope universality "
    "holds across classes."
)
body(
    "Predicted intercept alpha_det(r) = prefix_steps(r) + (3/(ln 4 - ln 3)) * "
    "log(a_final(r) / 2^k). Linear fit of actual alpha vs alpha_det: R^2 = 0.9854."
)
body(
    "Anderson-Darling 2-sample test on per-class residuals (496 pairs, "
    "n = 20K per class): same-a_final pairs (110/496) have median AD statistic "
    "0.79, max 13.0; different-a_final pairs (386/496) have median AD 15.7, "
    "max 162.6. Ratio of medians 19.9x: a_final is the primary clustering "
    "axis distributionally, not just for the mean. Note: at this n, AD is a "
    "sensitive test, and ~33% of same-a_final pairs still reject the same-"
    "distribution null at alpha = 0.05 -- a_final captures most but not all of "
    "the per-class distribution; a smaller secondary modulation (consistent "
    "with c_final structure) remains."
)

# --- Generalization across statistics ---
section("5. Does the decomposition hold for an independent functional?")
body(
    "sigma, syracuse, odd_steps, even_steps are linearly dependent by definition "
    "(sigma = odd_steps + even_steps, syracuse = odd_steps), so a single fit "
    "transfers to all four mechanically. The genuinely independent test is "
    "log(max_excursion), the log of the trajectory peak, which has no algebraic "
    "tie to the step-count statistics:"
)
pdf.set_font("Courier", "", 7.5)
pdf.cell(0, 3.4, "  outcome             mu_beta  heuristic  tau^2_corr  alpha_R^2  AD_ratio",
         new_x="LMARGIN", new_y="NEXT")
pdf.cell(0, 3.4, "  log(max_excursion)    1.002      --       0.000      0.8243      1.5x",
         new_x="LMARGIN", new_y="NEXT")
pdf.cell(0, 3.4, "  sigma (for reference) 10.401   10.428      0.000      0.9854     19.9x",
         new_x="LMARGIN", new_y="NEXT")
pdf.set_font("Helvetica", "", 8.5)
body(
    "The prefix decomposition partially explains log(max_excursion): R^2 = 0.82 "
    "between predicted-alpha and per-class peak intercept, AD ratio 1.5x. "
    "Substantial but smaller than for the step-count functionals -- peak amplitude "
    "depends on post-prefix trajectory chaos beyond the prefix algebra's reach, "
    "so the closed-form covariate captures a large but incomplete fraction of "
    "the per-class variation."
)

# --- Bayesian apples-to-apples ---
section("6. Bayesian replication in B&B framework (Pathfinder VI, full N_train = 500K)")
body(
    "Stan NB2 GLM with log link, priors beta ~ N(0,5), phi ~ LogN(log 5, 1) "
    "(phi >= 0.1), sigma_u ~ HN(0,2). Inference via cmdstanpy 1.3.0 Pathfinder "
    "(Zhang, Carpenter, Gelman, Vehtari, JMLR 2022), 4 paths, 1000 importance-"
    "resampled draws. Production HMC at this scale hit a Stan 2.36 unified-mode "
    "multi-chain lockup (sync issue with num_chains=N in shared-thread-pool "
    "process, distinct from the well-known reduce_sum grainsize pathology); "
    "Pathfinder converged cleanly. Held-out test set N_test = 50K."
)
pdf.set_font("Courier", "", 7.5)
pdf.cell(0, 3.4, "  spec  n_params   log_score    delta_vs_B0   sigma_u_mean",
         new_x="LMARGIN", new_y="NEXT")
pdf.cell(0, 3.4, "  B0    3          -274,150.3   --            --",
         new_x="LMARGIN", new_y="NEXT")
pdf.cell(0, 3.4, "  B1    12         -274,138.4   +12           0.0025  (RE on n mod 8)",
         new_x="LMARGIN", new_y="NEXT")
pdf.cell(0, 3.4, "  B2    6          -273,288.4   +862          --      (FE on a_final, k=3)",
         new_x="LMARGIN", new_y="NEXT")
pdf.cell(0, 3.4, "  B3    9          -272,435.3   +1,715        --      (FE on a_final, k=6)",
         new_x="LMARGIN", new_y="NEXT")
pdf.cell(0, 3.4, "  B4    18         -272,438.4   +1,712        0.0094  (B3 + RE on mod 8)",
         new_x="LMARGIN", new_y="NEXT")
pdf.set_font("Helvetica", "", 8.5)
body(
    "Pathfinder caveat: it is a quasi-Newton VI approximation that systematically "
    "underestimates posterior dispersion in non-Gaussian hierarchical regions, "
    "and so collapses sigma_u in B1/B4. Absolute log scores for hierarchical "
    "specs are below what HMC would produce; fixed-effect specs (B0, B2, B3) "
    "are not subject to this bias. HMC validation at production scale is the "
    "natural next step on the Bonacorsi side. Within Pathfinder, B2 and B3 "
    "(a_final fixed effects at k=3 and k=6) replace mod-8 RE with a closed-form "
    "covariate at the same modular resolution, with comparable predictive "
    "content and fewer parameters; B3 with 9 parameters scores 477 nats above "
    "B&B's reported HMC NB2-GLM, suggestive that finer modular resolution "
    "captures more structure than mod 8 (testable directly under HMC)."
)

# --- Status ---
section("7. Status")
body(
    "Code, model, scripts, summaries, and figures: private GitHub repo "
    "(github.com/mrnathanhumphrey-droid/Collatz). In correspondence with "
    "N. Bonacorsi (Columbia) regarding production-scale HMC validation and "
    "potential joint follow-up. Prior art audit (literature_check.md) identifies "
    "Terras 1976 Lemma 4 as the asymptotic predecessor; the contribution is "
    "(i) the exact (a_final, c_final) symbolic state with variable prefix length, "
    "and (ii) supplying the closed-form covariate that grounds Bonacorsi-Bordoni's "
    "mod-8 random effect in the symbolic Collatz prefix. A companion analysis "
    "extends the prefix decomposition to the qx+1 generalization (q in {5, 7, 9, 11}); "
    "reported separately."
)

# --- Refs ---
section("8. References")
body(
    "Terras, R. (1976). A stopping time problem on the positive integers. "
    "Acta Arithmetica 30, 241-252.\n"
    "Bonacorsi, N. & Bordoni, F. (2026). Bayesian modeling of Collatz stopping "
    "times. arXiv:2603.04479.\n"
    "Zhang, L., Carpenter, B., Gelman, A., Vehtari, A. (2022). Pathfinder: "
    "Parallel quasi-Newton variational inference. JMLR 23(306), 1-49. "
    "arXiv:2108.03782."
)

pdf.output(str(OUT))
print(f"[save] {OUT}")
