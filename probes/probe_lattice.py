"""
PROBE LATTICE -- is the Syracuse tower a Lapidus self-similar p-adic string? (2026-07-27)

Tests whether the log-oscillation in Lambda_i is EXACTLY log-3-periodic (lattice, a theorem: bounded,
no crossing) or drifts (non-lattice / irrational-rotation, quasi-periodic). Zeta from SCALING RATIOS,
NOT operator diagonalization (that was the R29 mush). Certified T_i machinery only.

L-A rational-zeta test:  finite Hankel rank of Lambda_i  <=>  rational zeta in 3^-s  <=>  finite mode set.
L-B log-period, decisively: fit omega on full i=8..20 + sub-windows; lattice => stable at 2pi/log3.
L-C zeta poles vs spectrum: D from rho1 = 3^-D; recognizable?
L-D boundedness corollary if lattice.
"""
import numpy as np
from numba import njit, prange

@njit(cache=False)
def build_dl(q, Nn):
    DL = np.empty(q, np.uint32)
    for i in range(q): DL[i] = np.uint32(0xFFFFFFFF)
    g = 1
    for s in range(Nn): DL[g] = np.uint32(s); g = (g * 4) % q
    return DL

@njit(parallel=True, cache=False)
def scat(base, p, Nn):
    n = base.shape[0]; d = np.empty(n, np.int64); pu = np.uint64(p); Nu = np.uint64(Nn)
    for j in prange(n):
        idx = (np.uint64(base[j]) * pu) % Nu; ii = np.int64(idx); d[j] = ii - ii // 3 - 1
    return d

def T_level(n, q=3, vm=64):
    qN = q ** (n + 1); Nn = q ** n; inv2 = pow(2, -1, Nn)
    r = np.arange(Nn); cpr = r[r % q != 0]; ncp = len(cpr); base = ((q * cpr + 1) % Nn).astype(np.uint32)
    ps = []; p = 1
    for v in range(1, vm + 1): p = (p * inv2) % Nn; ps.append((p, 0.5 ** v))
    pi = np.full(ncp, 1.0 / ncp)
    for it in range(500):
        nxt = np.zeros(ncp)
        for p, w in ps: nxt += np.bincount(scat(base, p, Nn), weights=w * pi, minlength=ncp)
        nxt /= nxt.sum()
        if np.abs(nxt - pi).max() < 1e-14: pi = nxt; break
        pi = nxt
    nu = pi / pi.sum(); DL = build_dl(qN, Nn); Y = (q * cpr + 1) % qN
    rho = np.bincount(DL[Y].astype(np.int64), weights=nu, minlength=Nn)
    return q ** n * sum((4.0 ** -k) * float(np.dot(rho, np.roll(rho, k))) for k in range(1, 64))


def main():
    T = {0: 1.0 / 3}
    for n in range(1, 15): T[n] = T_level(n)
    T.update({15: 0.23567582169638104, 16: 0.23591007771310188, 17: 0.2361167285902292,
              18: 0.23629629695785445, 19: 0.23645299645774767, 20: 0.2365885345278137})
    I = np.arange(1, 21); Lam = np.array([T[i] - T[i - 1] for i in I])
    L3 = np.log(3); L2 = np.log(2)
    print("# PROBE LATTICE\n")
    print("Lambda_1..20 (x1e-3):", " ".join(f"{x*1e3:.4f}" for x in Lam), "\n")

    # ---------- L-A: rational-zeta test via Hankel rank ----------
    print("## L-A  rational zeta in 3^-s  <=>  finite Hankel rank of Lambda_i")
    seq = Lam[7:]                      # i=8..20, the clean-monotone-decay regime
    h = len(seq) // 2
    H = np.array([[seq[i + j] for j in range(h)] for i in range(len(seq) - h + 1)])
    sv = np.linalg.svd(H, compute_uv=False); sv = sv / sv[0]
    print("  Hankel singular values (norm):", " ".join(f"{s:.1e}" for s in sv))
    for tol in (1e-4, 1e-6, 1e-8):
        print(f"   numerical rank @ {tol:.0e}: {int(np.sum(sv > tol))}  (= # complex dimensions if finite)")
    ratios = Lam[8:] / Lam[7:-1]
    print(f"  Lambda ratios i=9..20: {[round(r,4) for r in ratios]}")
    print(f"  => per-level ratio ~{np.mean(ratios[-6:]):.4f}; is log(ratio)/log3 rational/nice?  "
          f"log(ratio)/log3 = {np.log(np.mean(ratios[-6:]))/L3:.4f}")
    print(f"  1/3={1/3:.4f} (r=1/3 => ratio should be 0.333). Measured ~0.87 != 1/3 => NOT a bare r=1/3 string.\n")

    # ---------- L-B: log-period, decisively ----------
    print("## L-B  measured log-period + units reconciliation")
    print(f"  Lapidus lattice period in log-modulus: 2pi/log3 = {2*np.pi/L3:.3f}")
    print(f"  base-2 candidate:                      2pi/log2 = {2*np.pi/L2:.3f}")
    print(f"  *** UNITS: a lattice-3 string sampled at INTEGER levels i (modulus x3/level) has")
    print(f"      oscillation exp(2pi i n * i * log3/log3) = exp(2pi i n i) = 1 -- period EXACTLY 1 level,")
    print(f"      ALIASED TO CONSTANT / INVISIBLE. So an observed ~9-level period is NOT a base-3 lattice mode.")
    print(f"      It matches 2pi/log2={2*np.pi/L2:.2f} (base-2, the /2^v multiplier), and log3/log2={L3/L2:.4f} is IRRATIONAL.\n")
    def fit_omega(idx):
        y = np.log(Lam[idx]); x = I[idx]
        b, a = np.polyfit(x, y, 1); res = y - (a + b * x)          # detrend at data's own rho1
        # scan omega, fit A cos(w i)+B sin(w i)
        best = (1e9, None)
        for w in np.linspace(0.2, 1.5, 4000):
            Dm = np.c_[np.cos(w * x), np.sin(w * x)]
            c, *_ = np.linalg.lstsq(Dm, res, rcond=None); r = res - Dm @ c
            ss = np.sum(r ** 2)
            if ss < best[0]: best = (ss, w)
        return np.exp(b), best[1], 2 * np.pi / best[1]
    for name, idx in (("full i=8..20", I >= 8), ("early i=8..14", (I >= 8) & (I <= 14)),
                       ("late  i=14..20", (I >= 14) & (I <= 20))):
        rho1, w, P = fit_omega(idx)
        print(f"  {name:15s}: rho1={rho1:.4f}  omega={w:.4f}  period_in_i={P:.2f}")
    print("  LATTICE => omega stable across windows at the predicted value; DRIFT => irrational-rotation live.")
    print("  (caveat: 13 pts fit omega poorly; STABLE omega = strong evidence, DRIFT = weak/noise.)\n")

    # ---------- L-C: D from the real rate ----------
    print("## L-C  dimension from the real decay rate")
    rho1 = np.exp(np.polyfit(I[I >= 12], np.log(Lam[I >= 12]), 1)[0])
    D = -np.log(rho1) / L3
    print(f"  rho1 (per level, i>=12) = {rho1:.4f}  =>  rho1 = 3^-D  =>  D = {D:.4f}")
    print(f"  recognizable? 1/8={1/8:.4f}  2/15={2/15:.4f}  1/7={1/7:.4f}  (2-log2/log3-1)?={round(2-L2/L3-1,4)}")
    print()

    # ---------- L-D ----------
    print("## L-D  boundedness corollary")
    print("  IF lattice: oscillatory part = sum_n c_n 3^-(D+2pi i n/log3) i = BOUNDED periodic modulation")
    print("     of the real 3^-Di decay => provably no sign change once real part dominates")
    print("     => upgrades i=20 'no crossing' from OBSERVATION to THEOREM.")
    print("     Does NOT cover a hidden REAL mode rho3 in (0.999,1) -- a real pole near 1, not on the lattice.")


if __name__ == "__main__":
    main()
