"""Generate one-sheet PDF for Yosef — qx+1 convergence prediction via Cramer."""
from fpdf import FPDF
from pathlib import Path

OUT = Path(r"C:\Users\Nate\Desktop\collatz_qx1_one_sheet_yosef_2026-05-01.pdf")

pdf = FPDF(orientation="portrait", unit="mm", format="letter")
pdf.set_auto_page_break(auto=True, margin=12)
pdf.add_page()
pdf.set_margins(left=15, top=12, right=15)

pdf.set_font("Helvetica", "B", 12.5)
pdf.cell(0, 5.5, "qx+1 Convergence Prediction via Cramer's Theorem on Prefix Complexity",
         new_x="LMARGIN", new_y="NEXT")
pdf.set_font("Helvetica", "", 8.5)
pdf.cell(0, 4, "Empirical note. Nathan Humphrey. May 2026.",
         new_x="LMARGIN", new_y="NEXT")
pdf.ln(1.5)

def section(title):
    pdf.set_font("Helvetica", "B", 9)
    pdf.cell(0, 4.2, title, new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 8)

def body(text):
    pdf.multi_cell(0, 3.4, text, new_x="LMARGIN", new_y="NEXT")
    pdf.ln(0.3)

def mono(text):
    pdf.set_font("Courier", "", 7.3)
    pdf.cell(0, 3.2, text, new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 8)

# 1 Setup
section("1. Setup")
body(
    "qx+1 dynamics: T_q(n) = n/2 if n even, q*n + 1 if n odd. For q=3 the orbits are "
    "conjectured to all reach 1 (Collatz). For q in {5, 7, 9, 11, ...}, heuristic factor "
    "q/4 > 1 implies divergence on average; only a small fraction of orbits return to 1."
)

# 2 Prefix decomposition
section("2. Prefix decomposition")
body(
    "Symbolic prefix algorithm starting at state (a=2^k, c=r) for residue r mod 2^k: "
    "while a even, if c even then a/=2 and c/=2 (c-even branch), else a*=q and c = q*c+1 "
    "(c-odd branch). Terminates when a becomes odd (after exactly k c-even branches and j "
    "= odd_steps_in_prefix c-odd branches; j varies by residue class). Post-prefix state = "
    "a_final * m + c_final, where a_final = q^j and m = (n-r)/2^k. Equivalent labels at "
    "fixed k: odd_steps = log_q(a_final) = prefix_total_steps - k. They are 100% collinear "
    "by construction."
)

# 3 Empirical law
section("3. Empirical law (per-class convergence rate)")
body(
    "Group residue classes by j and compute conv_rate(j) = fraction reaching 1 within max_steps."
)
mono("  q=5 N=10^8  k=8  32,785 conv  slope=-0.5619  R^2=0.999")
mono("  q=7 N=10^8  k=6     258 conv  slope=-1.3685  R^2=0.999")
mono("  q=9 N=10^8  k=6     104 conv  slope=-2.0529  R^2=0.994")
mono("  q=11 N=10^9 k=6      36 conv  slope=-1.6458  R^2=0.993")
body(
    "Each: log(conv_rate(j)) is essentially exactly linear in j with R^2 >= 0.99."
)

# 4 Trajectory measure
section("4. Trajectory measure (unconditional v_2 distribution)")
body(
    "Sampled 5,000 random odd n in [1, 10^7], ran q=5 trajectories for 200 Syracuse steps each, "
    "recorded v_2 = nu_2(qx+1) per Syracuse step. Pooled 999,695 step records."
)
mono("  E[v_2]    = 2.0028   (Geom(1/2) prediction: 2.0,   match to 0.14%)")
mono("  Var[v_2]  = 2.0056   (Geom(1/2) prediction: 2.0,   match to 0.28%)")
mono("  empirical drift = log(q) - E[v]*log(2) = 0.2212  (predicted log(5/4) = 0.2231)")
body(
    "Trajectory measure on v_2 IS i.i.d. Geom(1/2) at q=5 to ~0.3% precision."
)

# 5 Cramer derivation
section("5. Cramer's theorem (exact, not Gaussian small-step)")
body(
    "Random walk in log-coords with steps X = log(q) - v*log(2), v ~ Geom(1/2). MGF:"
)
mono("  E[exp(-theta*X)] = q^(-theta) * (2^theta/2) / (1 - 2^theta/2)")
body(
    "Cramer's theorem: P(reach 1 from log_value V) = exp(-theta * V), where theta is the unique "
    "positive root of E[exp(-theta*X)] = 1. Setting MGF = 1 and simplifying:"
)
mono("            q^(-theta) = 2^(1-theta) - 1                  (*)")
body(
    "Starting log-value at j-class is V ~ j*log(q) + const(N, k). Slope of log(conv_rate) "
    "with respect to j: -theta(q) * log(q). Predicted at q=5: theta(5) ~= 0.34908, slope_pred = "
    "-0.5618. Empirical slope = -0.5619. Match to four significant figures."
)

# 6 Comparison table
section("6. Cramer prediction vs empirical")
mono("  q  theta_Cramer  slope_pred  slope_emp  ratio_emp/pred  n_conv")
mono("  5     0.34908      -0.5618    -0.5619       1.0001      32,785")
mono("  7     0.62650      -1.2191    -1.3685       1.1226         258")
mono("  9     0.74189      -1.6301    -2.0529       1.2594         104")
mono(" 11     0.80415      -1.9283    -1.6458       0.8535          36")
body(
    "q=5 matches to 0.01%. Deviations at q=7,9,11 are likely finite-sample (n_conv ranges "
    "32,785 down to 36; slope fits at q=9,11 are dominated by single-orbit j-bins). The exact "
    "Cramer formula (*) is q-independent in form; the empirical 'C ~= 2.5' clustering at q in {5,7,9} "
    "is the value of theta(q)*log(q)/log(q/4) in that q range, not a universal constant."
)

# 7 Substantive claim
section("7. Substantive claim")
body(
    "P(orbit starting in residue class with prefix odd-step count j reaches 1 in qx+1) = "
    "A(q) * exp(-theta(q) * log(q) * j), where theta(q) solves equation (*). The trajectory "
    "v_2 distribution is i.i.d. Geom(1/2), verified empirically. The constant theta(q) is "
    "derivable from the discrete-step large-deviation Cramer theorem applied to the qx+1 "
    "single-step log-change, with no free parameters."
)

# 8 Open questions
section("8. Open questions")
body(
    "(i) Verify Cramer prediction at q=7,9,11 with N high enough for 10^4 convergent orbits "
    "per q (q=11 at N=10^11 would suffice). (ii) The q=9 deviation (composite q) may indicate "
    "structure beyond i.i.d. Geom(1/2); worth checking with sampled v_2 trajectory measure "
    "at q=9 specifically. (iii) Connection to Sinai (2003) statistical 3x+1 result, which "
    "treats q=3 (mean-reverting) regime; current result is for q >= 5 (transient) regime."
)

# 9 References
section("9. References")
body(
    "Cramer, H. (1938). Sur un nouveau theoreme-limite de la theorie des probabilites. "
    "Actualites Scientifiques 736.\n"
    "Terras, R. (1976). A stopping time problem on the positive integers. "
    "Acta Arithmetica 30, 241-252.\n"
    "Sinai, Ya. G. (2003). Statistical (3x+1) problem. Comm. Pure Appl. Math 56(7), 1016-1028.\n"
    "Bonacorsi, N. & Bordoni, F. (2026). Bayesian modeling of Collatz stopping times. "
    "arXiv:2603.04479."
)

pdf.output(str(OUT))
print(f"[save] {OUT}")
