"""
PROBE R19 -- TYPICAL OR EXCEPTIONAL. Reuses R7/R9/R10. Per-r shape statistics only (NO cross-r rates/periods).

Two spikes to disambiguate:
  ADDITIVE:  mu_hat_r(xi), xi in Z/3^r  (R18-D found max/typical GROWS).
  A-SIDE:    A_r(m) = gamma_r(tau_m) - gamma_{r-1}(tau_m) = C_{r+1}(m)/3  (REAL exact rational; the moments of the
             dlog angular profile |theta_hat|^2, R12-B). A_r(0)=S_r (DC/uniform, excluded from delta).
R19-A additive argmax + coherence.  R19-B fixed-m typical-or-exceptional (THE decider).  R19-C within-stratum
spiking.  R19-D does the additive spike land ON the fixed small-m coefficients (same obstruction) or away.
R19-E max|mu_hat| six points + ratios column (labelled ratios, NOT a rate).
"""
import os, sys, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import probe_engine_R7 as R7
import probe_gamma_R9 as R9
import probe_charledger_R10 as R10

v3 = R9.v3


def mu_hat(mu_r, r, xi):
    N = 3 ** r
    return sum(complex(p) * cmath.exp(2j * math.pi * (xi * a % N) / N) for a, p in mu_r.items())


def A_spectrum(mu, r):
    """A_r(m) = gamma_r(tau_m) - gamma_{r-1}(tau_m), m=1..3^r-1 (REAL exact). Returns dict m->A_r(m)."""
    A = {}
    for m in range(1, 3 ** r):
        A[m] = R9.gamma(mu[r], r, R9.tau(m, r)) - R9.gamma(mu[r - 1], r - 1, R9.tau(m, r - 1))
    return A


def main():
    print("# PROBE R19 -- TYPICAL OR EXCEPTIONAL. Additive spike vs A-side fixed-m coefficients.\n")
    mu = {1: R7.mu1()}
    for k in range(2, 8):
        mu[k] = R7.build_mu(mu[k - 1], k)

    # additive coefficient magnitudes over primitive xi
    addmag = {}
    for r in range(2, 8):
        N = 3 ** r
        addmag[r] = {xi: abs(mu_hat(mu[r], r, xi)) for xi in range(1, N) if xi % 3 != 0}

    # ================= R19-A =================
    print("## R19-A  THE ARGMAX (measurement): additive argmax_{3|/xi} |mu_hat_r(xi)|, coherence, top-5")
    argmax = {}
    for r in range(2, 8):
        N = 3 ** r
        xs = addmag[r]
        xstar = max(xs, key=xs.get); argmax[r] = xstar
        coh = "n/a"
        if r > 2:
            coh = f"argmax_r mod 3^{r-1} = {xstar % 3**(r-1)}  (argmax_{{r-1}}={argmax[r-1]}; congruent? {xstar % 3**(r-1) == argmax[r-1] % 3**(r-1)})"
        top5 = sorted(xs.items(), key=lambda kv: -kv[1])[:5]
        print(f"   r={r} (N={N}): argmax xi*={xstar}  xi*/N={xstar/N:.5f}  v3(xi*-1)={v3(xstar-1)} v3(xi*+1)={v3(xstar+1)}")
        print(f"        coherence: {coh}")
        print(f"        top5 (xi:|mu_hat|): " + "  ".join(f"{xi}:{val:.5f}" for xi, val in top5))
    print("   [Q: fixed frequency / tower-coherent / wandering?]\n")

    # A-spectra (exact rationals) for r=3..7
    print("   (building A-spectra A_r(m), m=1..3^r-1, exact rationals, r=3..7 ...)")
    Aspec = {r: A_spectrum(mu, r) for r in range(3, 8)}
    normA = {r: sum(a * a for a in Aspec[r].values()) for r in Aspec}   # ||delta||^2_A = Sum_{m!=0} A_r(m)^2

    # ================= R19-B =================
    print("\n## R19-B  TYPICAL-OR-EXCEPTIONAL (measurement; THE decider): |A_r(m)|^2 vs stratum-typical")
    print("   typical = (||delta||^2_A / r) / N_{v3(m)},  N_j = 2*3^{r-1-j} (count of m in stratum v3=j)")
    MS = [1, 2, 3, 4, 9, 27]
    print(f"   {'r':>2} " + " ".join(f"m={m}(v3={v3(m)})".rjust(15) for m in MS))
    for r in range(3, 8):
        cells = []
        for m in MS:
            j = v3(m)
            if m >= 3 ** r:
                cells.append("--".rjust(15)); continue
            A2 = float(Aspec[r][m]) ** 2
            Nj = 2 * 3 ** (r - 1 - j)
            typ = (float(normA[r]) / r) / Nj
            cells.append(f"{A2/typ:>15.3f}")
        print(f"   {r:>2} " + " ".join(cells) + "   (ratios |A|^2/typical)")
    print(f"   raw |A_r(m)|^2 (float):")
    for r in range(3, 8):
        vals = " ".join(f"m{m}:{float(Aspec[r][m])**2:.3e}" if m < 3**r else f"m{m}:--" for m in MS)
        print(f"   r={r}: ||d||^2_A={float(normA[r]):.4f}  {vals}")
    print("   [Q: ratio O(1) stable in r  => equipartition survives;  growing => fixed-m coeffs EXCEPTIONAL, route dies.]\n")

    # ================= R19-C =================
    print("## R19-C  WITHIN-STRATUM DISTRIBUTION (measurement, NO fit): does A-side spike like R18-D additive?")
    print(f"   {'r':>2} {'j':>2} {'#members':>9} {'max/typ (within-stratum |A|^2)':>32} {'top3 frac of stratum mass':>26}")
    for r in range(5, 8):
        for j in (0, 1, 2):
            members = [m for m in range(1, 3 ** r) if v3(m) == j]
            a2 = [float(Aspec[r][m]) ** 2 for m in members]
            tot = sum(a2); mx = max(a2); typ = tot / len(a2)
            top3 = sum(sorted(a2, reverse=True)[:3])
            print(f"   {r:>2} {j:>2} {len(members):>9} {mx/typ:>32.3f} {top3/tot:>26.4f}")
    print("   [Q: does additive-side spiking (R18-D) reproduce on the A-side, and in which strata?]\n")

    # ================= R19-D =================
    print("## R19-D  SPIKE LOCATION vs FIXED m: does the additive spike land ON the small-m coefficients?")
    print("   U support law (R11-A/R12-A): U(k,xi)=0 unless v3(k)=v3(xi) AND k=xi mod 3 => transported stratum = v3(xi*).")
    for r in range(3, 8):
        xstar = argmax[r]
        # A-side argmax over ALL m, and over the small set
        mstar = max(Aspec[r], key=lambda m: abs(Aspec[r][m]))
        small_present = [m for m in [1, 2, 3, 4, 9, 27] if m < 3 ** r]
        msmall = max(small_present, key=lambda m: abs(Aspec[r][m]))
        jxi = v3(xstar)
        print(f"   r={r}: additive xi*={xstar} v3(xi*)={jxi} (=transported A-stratum)  |  "
              f"A-side argmax m*={mstar} v3(m*)={v3(mstar)} |A|={float(abs(Aspec[r][mstar])):.4f}  |  "
              f"in small set {small_present}? {mstar in small_present}")
    print("   [Q: transported additive-stratum v3(xi*) == v3 of the dominant/fixed small-m? on=same obstruction, away=harmless.]\n")

    # ================= R19-E =================
    print("## R19-E  max_xi |mu_hat_r(xi)| -- six points + RATIOS column (NOT a rate; do not extrapolate)")
    print(f"   {'r':>2} {'max|mu_hat_r| (float; algebraic, not rational)':>46} {'ratios':>10}")
    prev = None
    for r in range(2, 8):
        mx = max(addmag[r].values())
        rat = f"{mx/prev:.5f}" if prev else "--"
        print(f"   {r:>2} {mx:>46.8f} {rat:>10}")
        prev = mx
    print("   [six points reported; no exponent named, no extrapolation.]")


if __name__ == "__main__":
    main()
