"""
PROBE 41 -- does our r_q=1 boundary coincide with a named condition in Siegel's dissertation?
(Lead 1 from the Siegel read: identify the d=2 / r_3=1 / exceptional-point boundary in his
(p,q)-adic framework.)

SIEGEL Ch4 (Def 4.1, 4.2, eq 4.4/4.5, and the alpha_H(0)=1 degeneracy at eq 4.144-4.145).
For a p-Hydra map H (p = halving base), with branch derivatives H'_j(0) = a_j/d_j:
    alpha_H(t) = (1/p) sum_{j=0}^{p-1} H'_j(0) e^{-2 pi i j t}     (the Fourier SYMBOL)
    beta_H(t)  = (1/p) sum_{j=0}^{p-1} H_j(0)  e^{-2 pi i j t}
    gamma_H(t) = beta_H(t)/alpha_H(t)
    sigma_H    = log_p( sum_j H'_j(0) ) = 1 + log_p(alpha_H(0)),   alpha_H(0) = p^{sigma_H - 1}
  NON-SINGULARITY (Def 4.2): alpha_H(j/p) != 0 for all j in Z/pZ.
  DEGENERATE CASE (Siegel, eq 4.144-4.145, p277): at p=2, when alpha_H(0)=1, chi_hat_H picks up
    a factor gamma_H(1/2)*A_hat_H(t) that Siegel calls a "DEGENERATE MEASURE."

OUR MAP: qx+1 Syracuse as a 2-Hydra map (p=2): branch 0 = x/2 (H'_0=1/2, H_0(0)=0),
  branch 1 = (qx+1)/2 (H'_1=q/2, H_1(0)=1/2). So:
    alpha_H(0)   = (1+q)/4
    alpha_H(1/2) = (1-q)/4
    sigma_H      = log_2((q+1)/2)
    gamma_H(1/2) = 1/(q-1)

PRE-REGISTRATION (falsifier-first; my initial 'non-singularity = boundary' guess stated to LOSE).
------------------------------------------------------------------
H_NONSING (my ORIGINAL guess -- stated to lose): the boundary is Siegel non-singularity
    (alpha_H(j/2)=0 at q=3). PRED (honest): FALSE -- alpha_H(0)=(1+q)/4 and alpha_H(1/2)=(1-q)/4
    are NEVER 0 for odd q, so non-singularity holds for ALL q and does NOT single out q=3.
    Reporting this as a LOSS for my initial framing.
H_DEGEN (*** the corrected identification ***): our r_q=1 boundary = Siegel's alpha_H(0)=1
    (equivalently sigma_H=1) DEGENERATE case, which he flags explicitly at p=2. PRED: alpha_H(0)=1
    and sigma_H=1 EXACTLY at q=3, and nowhere else (odd q); this matches r_3=1 (gap closed) and
    r_q<1 <=> alpha_H(0)>1 (sigma_H>1) for q>=5.
    FALSIFIER: if alpha_H(0)=1 at some q!=3, or !=1 at q=3, the identification fails.
H_QUANT (measurement, NO verdict): is there a closed relation r_q vs alpha_H(0)/sigma_H? Report
    r_q vs 1/alpha_H(0)=2^{1-sigma_H}; NO fit committed (R28: r_q has no elementary closed form).

DECISION: H_DEGEN CONFIRMED iff {q : alpha_H(0)=1} = {3} = {q : r_q=1}. Exact.

NOT AT STAKE: R1-R40. This NAMES the boundary in Siegel's framework; it does not change r_q.
"""
from fractions import Fraction
from math import log2

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def siegel_quantities(q):
    """qx+1 as a 2-Hydra map. Exact rationals where possible."""
    a0 = Fraction(1, 2)          # H'_0(0)
    a1 = Fraction(q, 2)          # H'_1(0)
    sumder = a0 + a1             # = (q+1)/2
    alpha0 = sumder / 2          # alpha_H(0) = (1+q)/4
    alpha_half = (a0 - a1) / 2   # alpha_H(1/2) = (1-q)/4
    sigma = log2(float(sumder))  # sigma_H = log_2((q+1)/2)
    gamma_half = Fraction(1, q - 1)
    return alpha0, alpha_half, sigma, gamma_half


def main():
    log("# PROBE 41 -- our r_q=1 boundary in Siegel's (p,q)-adic framework")
    log("# H_NONSING (my guess, stated to lose): non-singularity = boundary?  H_DEGEN: alpha_H(0)=1 = boundary?")
    log("")
    R = {3: 1.0, 5: 0.62, 7: 0.39}   # r_q from R27/R32 (direct + operator)
    log(f"   {'q':>4} {'alpha_H(0)=(1+q)/4':>18} {'alpha_H(1/2)':>13} {'sigma_H':>9} "
        f"{'gamma_H(1/2)':>13} {'r_q':>6} {'2^(1-sigma)':>11}")
    degen_qs = []
    for q in [3, 5, 7, 11, 13, 17, 19]:
        a0, ah, sig, gh = siegel_quantities(q)
        two_1msig = 2.0 ** (1 - sig)
        rq = f"{R[q]:.2f}" if q in R else "-"
        if a0 == 1:
            degen_qs.append(q)
        log(f"   {q:>4} {str(a0):>18} {str(ah):>13} {sig:>9.5f} "
            f"{str(gh):>13} {rq:>6} {two_1msig:>11.5f}")
    log("")

    # ---- H_NONSING (my original guess) ----
    log("## H_NONSING (my ORIGINAL guess) -- is the boundary Siegel non-singularity alpha_H(j/2)=0?")
    nonsing_fail = [q for q in [3, 5, 7, 11, 13] if siegel_quantities(q)[0] == 0 or siegel_quantities(q)[1] == 0]
    log(f"   q with alpha_H(0)=0 or alpha_H(1/2)=0 (non-singularity FAILS): {nonsing_fail if nonsing_fail else 'NONE'}")
    log("   => alpha_H(0)=(1+q)/4 and alpha_H(1/2)=(1-q)/4 are NEVER 0 for odd q.")
    log("   => H_NONSING REFUTED (as I pre-committed it would be): non-singularity holds for ALL q,")
    log("      it does NOT single out q=3. My initial 'non-singularity = boundary' framing was WRONG.")
    log("")

    # ---- H_DEGEN (the corrected identification) ----
    log("## H_DEGEN (*** corrected ***) -- is the boundary Siegel's alpha_H(0)=1 (sigma_H=1) DEGENERACY?")
    log(f"   {{q : alpha_H(0)=1}} = {degen_qs}   (Siegel's degenerate-measure case, eq 4.144-4.145)")
    log(f"   {{q : r_q=1 (gap closed)}} = [3]   (r_3=1, R25/R32)")
    match = (degen_qs == [3])
    if match:
        log("   => H_DEGEN CONFIRMED: our r_q=1 boundary IS Siegel's alpha_H(0)=1 (sigma_H=1) case.")
        log("      alpha_H(0) = (1+q)/4 = 1  <=>  q=3  <=>  sigma_H = log_2((q+1)/2) = 1  <=>  r_3=1.")
        log("      Siegel EXPLICITLY flags alpha_H(0)=1 at p=2 as a 'degenerate measure' -- that is")
        log("      exactly our critical point / exceptional point (R39). For q>=5, alpha_H(0)>1")
        log("      (sigma_H>1, 'marginally expanding') and the gap opens (r_q<1).")
    else:
        log(f"   => H_DEGEN FAILED: degeneracy locus {degen_qs} != {{3}}. Identification wrong.")
    log("")

    # ---- H_QUANT (no verdict) ----
    log("## H_QUANT (measurement, NO verdict) -- r_q vs 2^{1-sigma_H} = 1/alpha_H(0):")
    for q in [5, 7]:
        a0, _, sig, _ = siegel_quantities(q)
        log(f"   q={q}: r_q={R[q]:.3f}  vs  1/alpha_H(0)={float(1/a0):.3f}  (2^(1-sigma)={2.0**(1-sig):.3f})  "
            f"-> close but NOT equal (R28: no elementary closed form)")
    log("")
    log("## READ -- the boundary now has a NAME in Siegel:")
    log("   r_q=1 (gap closed, EP) <=> sigma_H=1 (MARGINAL growth exponent) <=> alpha_H(0)=1")
    log("   (Siegel's degenerate measure). q>=5: sigma_H>1 (expanding), gap opens. This is a")
    log("   DYNAMICAL reading of the boundary (marginal vs expanding) + a citable home. It does")
    log("   NOT close L3 (the DECAY RATE / gap SIZE is still Siegel's open q-adic problem), but it")
    log("   identifies WHERE the transition sits in his framework. (Non-singularity was the wrong")
    log("   object; the marginal-exponent degeneracy is the right one.)")
    with open("result_41_siegel_alpha_boundary_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
