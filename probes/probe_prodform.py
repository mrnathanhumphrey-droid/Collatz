"""
PROBE PRODFORM -- test Wilson's load-bearing conjecture: is the accumulated modulation the lacunary product
   Chat(m) =? const * Prod_{j=0}^{r-1} w(2^{-j} m / N),   w(u)=1/(5-4 cos 2pi u),
over 3 nmid m?  If yes, A(K) = sum_{sum n_j 2^{-j} == K} Prod w_hat(n_j) has the RESONANCE structure Wilson uses.
If no, the resonance interpretation is void (though A(K)>=0 still holds -- see below).

Also:
 (2) A(K) >= 0 WITHOUT the product form: A(K) = invFFT(Chat)(K) = autocorrelation C(K) = sum_s rho(s)rho(s+K),
     and rho >= 0 (a density) => C(K) >= 0 trivially. Confirm min_K C(K) >= 0. This is the positivity we HAVE,
     and it does NOT need the product form -- only the resonance labeling of A(K) does.
 (3) Aggregate 3 nmid m mass decay vs r (Wilson's contribution): is it exponential (beating Tao's pointwise
     superpolynomial |nu_hat(xi)| << n^{-A}), or a normalization artifact?
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"


def prof_of(r):
    rho = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy"))
    p = np.abs(np.fft.fft(rho)) ** 2
    return p, len(rho)


def logcorr(a, b):
    la = np.log(a); lb = np.log(b)
    la -= la.mean(); lb -= lb.mean()
    return float((la * lb).sum() / np.sqrt((la * la).sum() * (lb * lb).sum()))


def main():
    t0 = time.time()
    print("# PROBE PRODFORM -- is Chat(m) = const * Prod_j w(2^{-j}m/N)?  (Wilson's load-bearing conjecture)\n")

    # ---------- (1) PRODUCT FORM TEST ----------
    print("## (1) product form: corr(log Chat, log Prod_{j<J} w(2^{-j}m/N)) over 3 nmid m, vs #factors J")
    for r in (12, 16):
        N = 3 ** r
        prof, _ = prof_of(r)
        m = np.arange(N)
        prim = (m % 3 != 0)
        y = prof[prim]                                    # measured Chat on 3 nmid m
        inv2 = pow(2, -1, N)
        # build the running product, sampling correlation at a few J
        logP = np.zeros(N)
        checks = sorted(set([1, 2, 4, r // 2, r - 1, r, r + 2, 2 * r]))
        row = []
        for J in range(1, max(checks) + 1):
            dil = pow(inv2, J - 1, N)                      # 2^{-(J-1)}
            mp = (dil * m) % N
            logP += np.log(1.0 / (5.0 - 4.0 * np.cos(2 * np.pi * mp / N)))
            if J in checks:
                c = logcorr(y, np.exp(logP[prim]))
                row.append((J, c))
        print(f"   r={r}:  " + "  ".join(f"J={J}:{c:+.3f}" for J, c in row) + f"   ({time.time()-t0:.1f}s)")
        del prof, y, m, prim, logP
    print("   [J=1 = MODREG single-w (~0). If corr climbs to ~1 near J=r, product form holds; if it stalls low, void.]\n")

    # ---------- (2) A(K) = C(K) >= 0  (positivity WITHOUT the product form) ----------
    print("## (2) A(K) = autocorr C(K) = invFFT(Chat)(K); rho>=0 => C(K)>=0 trivially. Confirm min_K C(K).")
    for r in (12, 16):
        prof, N = prof_of(r)
        C = np.fft.ifft(prof).real                        # = autocorrelation of rho (times N)
        C0 = C[0]
        print(f"   r={r}: min_K C(K)/C(0) = {C.min()/C0:+.3e}   (>=0 expected)   #neg = {int((C < -1e-9*C0).sum())}")
        del prof, C
    print("   [A(K)>=0 confirmed => Wilson's Num = (2N/3)[A(1)+A(-1)] - (N/3)*partials is a nonneg-vs-nonneg race,")
    print("    and this needs ONLY rho>=0, NOT the product form. The product form is only for the resonance labels.]\n")

    # ---------- (3) aggregate 3 nmid m mass decay ----------
    print("## (3) aggregate fluctuation mass vs r: S=sum_{3nmid m}Chat,  C0=sum_m Chat=N||rho||^2,  frac=S/C0")
    print(f"   {'r':>2} {'S':>13} {'C0':>13} {'S/C0':>12} {'S/C0 rate':>10}")
    prev = None
    for r in (12, 13, 14, 15, 16):
        prof, N = prof_of(r)
        m = np.arange(N)
        S = float(prof[m % 3 != 0].sum())
        C0tot = float(prof.sum())
        frac = S / C0tot
        rate = frac / prev if prev else float('nan')
        prev = frac
        print(f"   {r:>2} {S:>13.6e} {C0tot:>13.6e} {frac:>12.6e} {rate:>10.4f}")
        del prof, m
    print("   [if S/C0 (or S) decays geometrically (~const rate<1) => exponential aggregate, beats Tao superpoly n^-A;")
    print("    check if it's S itself or only the RATIO S/C0 that decays (normalization artifact test).]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
