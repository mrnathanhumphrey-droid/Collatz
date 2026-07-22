"""
PROBE R17 -- THE SLOW MODE. Reuses R7/R10. Transport symbol D and the deviation field in additive Fourier.

Renewal (R16, additive-Y): Y_r = 2^{-v}(1+3 Y_{r-1}) mod 3^r  =>
   mu_hat_r(xi) = E_v[ e(xi*2^{-v}/3^r) * mu_hat_{r-1}(xi*2^{-v} mod 3^{r-1}) ].
So the v-averaged transport phase is   D(xi) := E_v[ e(xi*2^{-v}/3^r) ],  v~Geom(1/2), p_v=2^{-v}.
Deviation field in Fourier: delta_r(xi) = mu_hat_r(xi) for xi != 0 (uniform has only DC).

R17-A: compute D directly, check (i) |D(xi)|=|D(-xi)|, (ii) mean|D|^2 vs 1/3, (iii) |D|^2 vs 1/(5-4cos(pi xi/3^r)).
   NO closed form assumed -- report match/mismatch and the true form.
"""
import os, sys, math, cmath, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import probe_engine_R7 as R7
import probe_charledger_R10 as R10

S = dict(R10.S)


def D_add(xi, r, Vmax=240):
    """DIRECT transport phase (additive Fourier): D(xi)=E_v[e(xi*(2^{-v} mod 3^r)/3^r)]."""
    N = 3 ** r; inv2 = pow(2, -1, N); acc = 0j; p = 0.5; u = inv2
    for v in range(1, Vmax + 1):
        acc += p * cmath.exp(2j * math.pi * (xi * u % N) / N)
        p *= 0.5; u = (u * inv2) % N
    return acc


def D_ideal2(xi, r):
    """IDEALIZED slow-mode symbol |D|^2 = 1/(5-4cos(pi xi/3^r)) (linear-in-v dlog shift -v/2, exact geom series)."""
    return 1.0 / (5 - 4 * math.cos(math.pi * xi / 3 ** r))


def mu_hat(mu_r, r, xi):
    N = 3 ** r
    return sum(complex(p) * cmath.exp(2j * math.pi * (xi * a % N) / N) for a, p in mu_r.items())


def main():
    print("# PROBE R17 -- THE SLOW MODE. Transport symbol D + deviation field (additive Fourier).\n")
    mu = {1: R7.mu1()}
    for k in range(2, 8):
        mu[k] = R7.build_mu(mu[k - 1], k)

    # ---- R17-A ----
    print("## R17-A  D CONVENTION (forced): D(xi)=E_v[e(xi 2^{-v}/3^r)] direct; check (i)(ii)(iii)")
    for r in range(2, 6):
        N = 3 ** r
        prim = [xi for xi in range(1, N) if xi % 3 != 0]
        Dvals = {xi: D_add(xi, r) for xi in range(N)}
        # (i) conjugation symmetry
        sym = max(abs(abs(Dvals[xi]) - abs(Dvals[(-xi) % N])) for xi in prim)
        # (ii) mean over ALL xi of |D|^2
        m2_all = sum(abs(Dvals[xi]) ** 2 for xi in range(N)) / N
        # (iii) closed forms
        dev_pi = max(abs(abs(Dvals[xi]) ** 2 - 1.0 / (5 - 4 * math.cos(math.pi * xi / N))) for xi in prim)
        dev_2pi = max(abs(abs(Dvals[xi]) ** 2 - 1.0 / (5 - 4 * math.cos(2 * math.pi * xi / N))) for xi in prim)
        print(f"   r={r} (N={N}): (i)|D(xi)|=|D(-xi)| max dev={sym:.2e}  "
              f"(ii)<|D|^2>_all={m2_all:.6f} (1/3={1/3:.6f}, dev={m2_all-1/3:+.2e})")
        print(f"         (iii) max| |D|^2 - 1/(5-4cos(pi xi/N)) |={dev_pi:.3e};  "
              f"vs 1/(5-4cos(2pi xi/N))={dev_2pi:.3e}")
    print("   [report: which closed form (if any) matches; is <|D|^2> exactly 1/3 or 1/3+O(2^-ord)?]\n")

    # ---- R17-D ----  (cheap forced, do early)
    print("## R17-D  FIXED-POINT FORM (forced): T self-map + invariance at achieved resolution")
    okD = True
    for r in range(2, 6):
        N = 3 ** r
        # self-map: 1+3*2^{-v}*X in 1+3Z always (trivially X=1 mod3 => 1+3*unit =1 mod3). verify by construction:
        selfmap = all((1 + 3 * pow(2, -1, 3 ** (r + 1)) * (1 + 3 * a)) % 3 == 1 for a in range(3))
        # invariance: mu_{r+1} folded mod 3^r == mu_r
        fold = {}
        for a, p in mu[r + 1].items():
            fold[a % N] = fold.get(a % N, F(0)) + p
        inv = (fold == mu[r])
        okD = okD and selfmap and inv
        print(f"   r={r}: 1+3*2^-v X in 1+3Z? {selfmap}   mu_{{r+1}} mod 3^r == mu_r (invariance)? {inv}")
    print(f"   => R17-D {'PASS -- self-map + achieved-resolution invariance certified' if okD else 'FAIL'}\n")

    # ---- R17-B ----
    print("## R17-B  QUASI-STATIONARY SELF-CONSISTENCY (measurement, NO fit): <|D|^2>_delta_r")
    print("   (D_ideal = closed-form slow-mode symbol |D|^2=1/(5-4cos(pi xi/3^r)); D_add = direct transport phase)")
    print(f"   {'r':>2} {'<|D_ideal|^2>_dr':>16} {'<|D_add|^2>_dr':>16} {'plain avg |D_ideal|^2':>22}")
    for r in range(3, 8):
        N = 3 ** r
        prim = [xi for xi in range(1, N) if xi % 3 != 0]
        num_i = num_a = den = 0.0; avg_i = 0.0
        for xi in prim:
            d2 = abs(mu_hat(mu[r], r, xi)) ** 2
            Di = D_ideal2(xi, r); Da = abs(D_add(xi, r)) ** 2
            num_i += d2 * Di; num_a += d2 * Da; den += d2; avg_i += Di
        print(f"   {r:>2} {num_i/den:>16.6f} {num_a/den:>16.6f} {avg_i/len(prim):>22.6f}")
    print("   [shape: flat in r (quasi-stationary) / drift up / down? per-r stat, no cross-r rate.]\n")

    # ---- R17-C ----
    print("## R17-C  LOCALIZATION SCALE (measurement, NO fit): |delta_r|^2-weighted x=xi/3^r moments")
    print(f"   {'r':>2} {'<x>':>10} {'<x^2>':>12} {'sqrt Var':>10} {'frac ||delta||^2 in x<1/12':>26}")
    for r in range(4, 8):
        N = 3 ** r
        prim = [xi for xi in range(1, N) if xi % 3 != 0]
        w = [abs(mu_hat(mu[r], r, xi)) ** 2 for xi in prim]
        xs = [min(xi / N, 1 - xi / N) for xi in prim]      # distance to 0 on the circle
        tot = sum(w)
        mx = sum(wi * xi for wi, xi in zip(w, xs)) / tot
        mx2 = sum(wi * xi ** 2 for wi, xi in zip(w, xs)) / tot
        frac = sum(wi for wi, xi in zip(w, xs) if xi < 1 / 12) / tot
        print(f"   {r:>2} {mx:>10.5f} {mx2:>12.6f} {math.sqrt(max(0,mx2-mx**2)):>10.5f} {frac:>26.4f}")
    print("   [width shrinks (squeezed into slow region) / holds (quasi-stationary) / spreads (leaking)?]\n")

    # ---- R17-E ----
    print("## R17-E  OUTSTANDING (one line each)")
    # symbol identity chi(4)w(chi) = 1/(4 - chi(2)^{-2}); verify numerically at r=3, few k
    r = 3; N = 3 ** r; ok_sym = True
    for k in (1, 2, 4):
        z4 = cmath.exp(2j * math.pi * k / N)            # chi_k(4)=zeta^k
        w = 1 / (4 * z4 - 1); lhs = z4 * w
        z2 = cmath.exp(2j * math.pi * k / (2 * N))       # chi_k(2)=zeta^{k/2}=e(k/2N); chi(2)^{-2}=e(-k/N)=z4^{-1}
        rhs = 1 / (4 - z2 ** (-2))
        if abs(lhs - rhs) > 1e-9: ok_sym = False
    print(f"   symbol identity chi(4)w(chi) = 1/(4 - chi(2)^-2): {ok_sym} (=1/(4-chi(4)^-1), verified numerically)")
    print("   R85 rung-1 n=8: one Bluestein/support-pruned probe, DEFER (unchanged R13-E/R16-E).")
    print("   branch weights / lambda-sweep: not specced here; carried as owed, no run.")


if __name__ == "__main__":
    main()
