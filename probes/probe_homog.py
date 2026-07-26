"""
PROBE HOMOG -- is the ln|pi-hat| spectral shape EXACTLY stationary (forced) or merely settling (emergent)? (2026-07-26)

Measurement only, pi-hat (ADDITIVE) throughout -- NO rho / channel / U^2 quantity enters (that crossing caused the
U^2/ell4 mix-up). L = ln|pi-hat_k(a)| over 3-nmid a (= units mod 3^k, count 2*3^{k-1}).

A GATE: reproduce banked skew -0.65 / exkurt +1.4 (k=8..16). Report a-set, count, mean, var per k.
B DRIFT (decider): standardized cumulants k3,k4,k5,k6 vs k; fit {constant, linear (slope/SE), geometric c_inf+B rho^k}.
   Pre-reg: |slope/SE|<2 all + geom not preferred => FORCED. |slope/SE|>3 any w/ resolved rho<1 => EMERGENT (rho=observable).
   drift resolved but rho not <1 => NEITHER. Deparity too (period-2 killer): c'_k=(c_k+c_{k+1})/2.
C HOMOGENEITY: partition units by v2(a), a mod 9, a mod 27; skew/kurt per class. Forced => class-independent.
   (<2> orbit: 2 is a PRIMITIVE root mod 3^k => <2>=all units => that partition is DEGENERATE; noted, use dlog2 bands.)
D LEFT TAIL: bottom 1% |pi-hat| at k=12(,16); v2(a), a mod 9, a/N (near 0 / 1/2 / 1?). Structured or diffuse?

Reuses fwd_hat. HEAVY (k up to 15/16). Run in background.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_singlerec import fwd_hat


def v2arr(a):
    a = a.astype(np.int64).copy()
    v = np.zeros(len(a), dtype=np.int64)
    for _ in range(30):
        even = (a % 2 == 0) & (a > 0)
        if not even.any():
            break
        v[even] += 1
        a[even] //= 2
    return v


def cumulants(L):
    mu = L.mean(); d = L - mu
    m2 = (d ** 2).mean(); m3 = (d ** 3).mean(); m4 = (d ** 4).mean()
    m5 = (d ** 5).mean(); m6 = (d ** 6).mean()
    s = np.sqrt(m2)
    k3 = m3 / s ** 3
    k4 = m4 / s ** 4 - 3
    k5 = (m5 - 10 * m3 * m2) / s ** 5
    k6 = (m6 - 15 * m4 * m2 - 10 * m3 ** 2 + 30 * m2 ** 3) / s ** 6
    return mu, m2, k3, k4, k5, k6


def linfit_se(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); n = len(x)
    xm = x.mean(); Sxx = ((x - xm) ** 2).sum()
    slope = (((x - xm) * (y - y.mean())).sum()) / Sxx
    inter = y.mean() - slope * xm
    resid = y - (inter + slope * x)
    s2 = (resid ** 2).sum() / max(n - 2, 1)
    se = np.sqrt(s2 / Sxx)
    return slope, se


def geomfit(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    best = None
    for rho in np.arange(0.20, 0.985, 0.005):
        A = np.column_stack([np.ones_like(x), rho ** x])
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        ss = ((y - A @ coef) ** 2).sum()
        if best is None or ss < best[0]:
            best = (ss, rho, coef[0], coef[1])
    ss, rho, cinf, B = best
    sst = ((y - y.mean()) ** 2).sum()
    return cinf, rho, B, (1 - ss / sst if sst > 0 else 0)


def fit_block(name, ks, seq):
    ks = np.asarray(ks, float); seq = np.asarray(seq, float)
    cmean, cscat = seq.mean(), seq.std()
    sl, se = linfit_se(ks, seq)
    cinf, rho, B, r2 = geomfit(ks, seq)
    # deparitied
    dep = 0.5 * (seq[:-1] + seq[1:]); kd = 0.5 * (ks[:-1] + ks[1:])
    sld, sed = linfit_se(kd, dep)
    print(f"   {name:>4}: const {cmean:+.4f}(sc {cscat:.4f}) | lin slope {sl:+.5f} SE {se:.5f} "
          f"|s/SE| {abs(sl/se) if se>0 else 0:.2f} | dep |s/SE| {abs(sld/sed) if sed>0 else 0:.2f} "
          f"| geom c_inf {cinf:+.4f} rho {rho:.3f} r2 {r2:.3f}")
    return abs(sl / se) if se > 0 else 0, abs(sld / sed) if sed > 0 else 0, rho, r2


def main():
    t0 = time.time()
    print("# PROBE HOMOG -- is ln|pi-hat| shape forced (exactly stationary) or emergent (settling)?\n")
    KS = list(range(8, 17))
    store = {}
    cum = {n: [] for n in ('k3', 'k4', 'k5', 'k6')}
    kdone = []
    print("## HOMOG-A GATE: L=ln|pi-hat(a)|, a-set = {1<=a<3^k : 3-nmid a} (units, count 2*3^{k-1})")
    print(f"   {'k':>2} {'count':>10} {'mu':>9} {'var':>7} {'skew':>7} {'exkurt':>7}")
    for k in KS:
        try:
            ph, N = fwd_hat(k)
        except MemoryError:
            print(f"   k={k}: MemoryError -- stop A/B at k={k-1}")
            break
        a = np.arange(N)
        U = a[a % 3 != 0]
        mag = np.abs(ph[U]).astype(np.float64)
        L = np.log(mag)
        mu, m2, k3, k4, k5, k6 = cumulants(L)
        print(f"   {k:>2} {len(U):>10} {mu:>9.4f} {m2:>7.4f} {k3:>7.3f} {k4:>7.3f}")
        cum['k3'].append(k3); cum['k4'].append(k4); cum['k5'].append(k5); cum['k6'].append(k6)
        kdone.append(k)
        if k in (12, 16):
            store[k] = (U.copy(), mag.copy())
        del ph, a, U, mag, L
    print("   [reproduce banked skew ~ -0.65 / exkurt ~ +1.4 at high k -- if not, STOP.]")

    # --- B DRIFT ---
    print("\n## HOMOG-B DRIFT (decider): standardized cumulants vs k; |slope/SE| is the number")
    print("   pre-reg: |s/SE|<2 all & geom not preferred => FORCED ; |s/SE|>3 any w/ rho<1 resolved => EMERGENT")
    flags = {}
    for n in ('k3', 'k4', 'k5', 'k6'):
        sSE, sSEd, rho, r2 = fit_block(n, kdone, cum[n])
        flags[n] = (sSE, sSEd, rho, r2)
    maxsSE = max(f[1] for f in flags.values())      # use DEPARITIED slope/SE
    verdict = ("FORCED (stationary shape -- theorem to prove)" if maxsSE < 2 else
               "EMERGENT (settling -- rho is the observable)" if maxsSE > 3 else
               "AMBIGUOUS (2<=|s/SE|<=3) -- extend k")
    print(f"   => max deparitied |slope/SE| over cumulants = {maxsSE:.2f}  ::  {verdict}")

    # --- C HOMOGENEITY ---
    print("\n## HOMOG-C HOMOGENEITY: shape per class (forced => class-independent). k=12")
    if 12 in store:
        U, mag = store[12]; L = np.log(mag)
        v2 = v2arr(U)
        print("   by v2(a):")
        for vv in range(0, 7):
            m = v2 == vv
            if m.sum() >= 50:
                _, _, k3, k4, _, _ = cumulants(L[m])
                print(f"      v2={vv}: n={m.sum():>8}  skew {k3:+.3f}  exkurt {k4:+.3f}")
        print("   by a mod 9:")
        for r in (1, 2, 4, 5, 7, 8):
            m = (U % 9) == r
            _, _, k3, k4, _, _ = cumulants(L[m])
            print(f"      a=={r} mod9: n={m.sum():>8}  skew {k3:+.3f}  exkurt {k4:+.3f}")
        print("   by a mod 27 (agg |skew| spread): ", end="")
        sk = []
        for r in range(27):
            if r % 3 == 0:
                continue
            m = (U % 27) == r
            if m.sum() >= 50:
                _, _, k3, _, _, _ = cumulants(L[m])
                sk.append(k3)
        print(f"mean skew {np.mean(sk):+.3f}, spread(std) {np.std(sk):.3f} over {len(sk)} classes")

    # --- D LEFT TAIL ---
    print("\n## HOMOG-D LEFT TAIL: bottom 1% of |pi-hat| -- structured or diffuse?")
    for kk in (12, 16):
        if kk not in store:
            continue
        U, mag = store[kk]; N = 3 ** kk
        thr = np.quantile(mag, 0.01)
        tail = mag <= thr
        Ut = U[tail]; v2t = v2arr(Ut)
        frac = Ut / N
        print(f"   k={kk}: bottom 1% = {tail.sum()} freqs")
        print(f"      v2(a) dist: " + ", ".join(f"v2={v}:{(v2t==v).mean()*100:.1f}%" for v in range(0, 6))
              + f", v2>=6:{(v2t>=6).mean()*100:.1f}%   (baseline: v2=0 is 50%, v2=1 25%, ...)")
        print(f"      a mod 9 dist: " + ", ".join(f"{r}:{((Ut%9)==r).mean()*100:.0f}%" for r in (1, 2, 4, 5, 7, 8)))
        print(f"      a/N: near0(<.1) {np.mean(frac<0.1)*100:.0f}%  mid(.4-.6) {np.mean((frac>0.4)&(frac<0.6))*100:.0f}%  "
              f"near1(>.9) {np.mean(frac>0.9)*100:.0f}%  median {np.median(frac):.3f}")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
