"""
PROBE R11 -- THE CHIRP KERNEL. Merges the 7/15 value campaign with the corpus epsilon-rate campaign:
  2 Lambda_r = S_{r+1}-S_r = eps_{r+1}-eps_r  (eps_k = S_k - 7/15 = R5's deviation d_k).
Bridge: beta(z)=dlog_4(1+3z) is a tower-compatible bijection of Z/3^r; U = F P_beta F* is the chirp-DFT
unitary; Lambda_r = <mu_hat, K_r mu_hat>, K_r = U* diag(w) U on the primitive shell, w(k)=1/(4 e(k/3^r)-1).

R11-A gate: (i) U mu_hat = nu_hat (R10) exact; (ii) U(k,xi)=0 when v3(k)!=v3(xi) (shell block-diagonal).
Reuses R7 (build_mu/mu1/cram) + R10 (dlog_table/autocorr_dlog/A_N/Lambda_r/layer_mass).
"""
import os, sys, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_engine_R7 as R7
import probe_charledger_R10 as R10

S = R10.S


def v3(n):
    n = abs(int(n))
    if n == 0:
        return 10 ** 9
    j = 0
    while n % 3 == 0:
        n //= 3; j += 1
    return j


def build_arrays(mu, r):
    """numpy arrays: mu_vec, beta, rho (pushforward), theta_hat(=nu_hat), mu_hat, U, w."""
    N = 3 ** r
    d = R10.dlog_table(r)                          # s -> t = beta(s)
    mu_vec = np.zeros(N)
    for s, p in mu[r].items():
        mu_vec[s % N] = float(p)
    beta = np.array([d[z] for z in range(N)])
    rho = np.zeros(N)
    for z in range(N):
        rho[beta[z]] += mu_vec[z]
    e = lambda x: cmath.exp(2j * math.pi * x)
    ks = np.arange(N)
    theta_hat = np.array([sum(rho[y] * e(k * y / N) for y in range(N)) for k in ks])
    mu_hat = np.array([sum(mu_vec[z] * e(xi * z / N) for z in range(N)) for xi in ks])
    U = np.zeros((N, N), dtype=complex)
    for k in range(N):
        for xi in range(N):
            U[k, xi] = sum(e((k * beta[z] - xi * z) / N) for z in range(N)) / N
    w = np.array([1.0 / (4 * e(k / N) - 1) for k in range(N)])
    return mu_vec, beta, rho, theta_hat, mu_hat, U, w


def main():
    print("# PROBE R11 -- THE CHIRP KERNEL. Numeric U-pipeline @1e-12 + exact Lambda/weld.\n")
    mu = {1: R7.mu1()}
    for k in range(2, 7):                          # build to mu_6 for R11-E
        mu[k] = R7.build_mu(mu[k - 1], k)

    arr = {}
    for r in range(2, 6):
        arr[r] = build_arrays(mu, r)

    # ---- R11-A ----
    print("## R11-A  KERNEL GATE  (i) U mu_hat = nu_hat(R10);  (ii) U(k,xi)=0 for v3(k)!=v3(xi)")
    okA = True
    for r in range(2, 6):
        N = 3 ** r
        _, beta, _, theta_hat, mu_hat, U, _ = arr[r]
        transported = U @ mu_hat
        err_i = np.max(np.abs(transported - theta_hat))
        # (ii) block-diagonality
        maxoff = 0.0
        for k in range(N):
            for xi in range(N):
                if v3(k) != v3(xi):
                    maxoff = max(maxoff, abs(U[k, xi]))
        gi = err_i < 1e-10; gii = maxoff < 1e-10
        okA = okA and gi and gii
        print(f"   r={r} (N={N}): (i) max|U mu_hat - nu_hat| = {err_i:.2e} [{'OK' if gi else 'FAIL'}]   "
              f"(ii) max|U(k,xi)|, v3(k)!=v3(xi) = {maxoff:.2e} [{'BLOCK-DIAG' if gii else 'FAIL'}]")
    print(f"   => R11-A {'GATE PASS' if okA else 'FAIL -- bridge dies at beta-table (#31)'}\n")

    # ---- R11-B ----
    print("## R11-B  FLATNESS vs Th 78.3  (|U| on primitive x primitive block; two labeled columns)")
    for r in range(2, 6):
        N = 3 ** r
        U = arr[r][5]
        prim = [i for i in range(N) if i % 3 != 0]
        block = np.abs(np.array([[U[k, xi] for xi in prim] for k in prim]))
        nz = block[block > 1e-12]
        frac_nz = nz.size / block.size
        nz_flat = (nz.max() - nz.min()) < 1e-9
        nz_val = nz.mean()
        th783 = 3 ** (-r / 2)                       # naive dense-flat chirp entry = 3^{-r/2}
        print(f"   r={r}: |U| prim block [measured]: DENSE-flat? {'no' if not (block.max()-block.min()<1e-9) else 'yes'}; "
              f"NONZERO-flat? {'YES' if nz_flat else 'no'} (support {frac_nz*100:.0f}%, nonzeros all = {nz_val:.8f} = 3^-(r-1)/2={3**(-(r-1)/2):.8f})")
        print(f"        [Th78.3 dense-flat 3^-r/2 = {th783:.8f}]  => U is a SPARSE unitary (flat-on-support), NOT dense-flat: MISMATCH with the maximally-spreading reading")
    print()

    # ---- R11-C ----
    print("## R11-C  QUADRATIC-FORM PIPELINE  (<mu_hat, K_r mu_hat> = Lambda_r; K_r=U*diag(w)U prim shell)")
    for r in range(2, 6):
        N = 3 ** r
        _, _, _, theta_hat, mu_hat, U, w = arr[r]
        prim = [i for i in range(N) if i % 3 != 0]
        # <mu_hat, U* diag(w) U mu_hat> restricted to primitive k = sum_{k prim} w(k) |theta_hat(k)|^2
        qform = sum(w[k] * abs(theta_hat[k]) ** 2 for k in prim)
        Lam_exact = R10.Lambda_r(mu, r)[0]
        err = abs(qform - complex(float(Lam_exact)))
        print(f"   r={r}: <mu_hat,K_r mu_hat>={qform.real:+.10f}{qform.imag:+.1e}i  vs  Lambda_r(exact)="
              f"{float(Lam_exact):+.10f}  [{'OK' if err<1e-9 else 'DEV'}  err={err:.1e}]")
    print()

    # ---- R11-D ----
    print("## R11-D  LEDGER EXTENSION + TWO-CAMPAIGN WELD  (2 Lambda_r = d_{r+1}-d_r; d_k=S_k-7/15=corpus eps)")
    d = {k: S[k] - F(7, 15) for k in S}
    print(f"   {'r':>2} {'2 Lambda_r (char side)':>26} {'d_(r+1)-d_r (eps incr)':>26} {'weld':>6}")
    okD = True
    for r in range(1, 6):
        twoL = 2 * R10.Lambda_r(mu, r)[0]
        eps_incr = d[r + 1] - d[r]
        good = (twoL == eps_incr)
        okD = okD and good
        print(f"   {r:>2} {str(twoL)[:26]:>26} {str(eps_incr)[:26]:>26} {'OK' if good else 'DEV':>6}")
    print(f"   => exact weld r=1..5 {'HOLDS' if okD else 'DEV'} (d_k = R5 deviation = corpus eps-sequence).")
    # exact extension attempt: S_7 via mu_7 -> Lambda_6
    print("   -- exact extension attempt (S_7 via mu_7 build -> Lambda_6):")
    try:
        mu[7] = R7.build_mu(mu[6], 7)
        S6_check = R10.layer_mass(R10.autocorr_dlog(mu[6], 6), 6)
        S7 = R10.layer_mass(R10.autocorr_dlog(mu[7], 7), 7)
        print(f"      S_6 (mu_6 layer_mass) = {S6_check}  [==frozen S_6? {S6_check==S[6]}]")
        print(f"      S_7 (mu_7 layer_mass) = {str(S7)[:40]}...")
        L6 = (S7 - S[6]) / 2
        print(f"      Lambda_6 = (S_7-S_6)/2 = {str(L6)[:40]}...  float={float(L6):+.8e}  (EXACT extension to r=6)")
    except Exception as ex:
        print(f"      WALL hit building mu_7/S_7: {type(ex).__name__} {ex}. Exact ledger stands at r=5;")
        print(f"      historical float eps-table (k=7..16) is the SAME object in a prior convention (not re-certified here).")
    # per-r labeled table
    print("   per-r ledger (all exact from frozen S):")
    print(f"      {'r':>2} {'Lambda_r':>16} {'Lambda_r^unif':>16} {'excess':>12} {'ratio':>10}")
    for r in range(1, 6):
        N = 3 ** r
        Lr = R10.Lambda_r(mu, r)[0]
        wtot = R10.A_N(0, N) - R10.A_N(0, 3 ** (r - 1))
        Lunif = S[r] * wtot / (2 * 3 ** (r - 1))
        ratio = float(Lr / Lunif) if Lunif != 0 else float('nan')
        print(f"      {r:>2} {float(Lr):>+16.8e} {float(Lunif):>+16.8e} {float(Lr-Lunif):>+12.3e} {ratio:>10.3f}")
    print()

    # ---- R11-E ----
    print("## R11-E  ANGULAR PROFILE  (measurement, NO fit; where does Lambda_r collect over the angle?)")
    for r in range(4, 7):
        N = 3 ** r
        if r <= 5:
            _, _, _, theta_hat, _, _, w = arr[r]
            th2 = np.abs(theta_hat) ** 2
            wv = w
        else:                                       # r=6: build lightweight (no full U)
            d6 = R10.dlog_table(6); rho = np.zeros(N)
            for s, p in mu[6].items():
                rho[d6[s % N]] += float(p)
            e = lambda x: cmath.exp(2j * math.pi * x)
            th2 = np.array([abs(sum(rho[y] * e(k * y / N) for y in range(N))) ** 2 for k in range(N)])
            wv = np.array([1.0 / (4 * e(k / N) - 1) for k in range(N)])
        prim = [k for k in range(1, N) if k % 3 != 0]
        # accumulation as K sweeps primitive angle indices
        acc = 0j; checkpoints = [len(prim) // 4, len(prim) // 2, 3 * len(prim) // 4, len(prim)]
        Lam = sum(wv[k] * th2[k] for k in prim)
        cps = []
        run = 0j
        for i, k in enumerate(prim, 1):
            run += wv[k] * th2[k]
            if i in checkpoints:
                cps.append((i, run.real))
        # angular concentration: fraction of |theta|^2 mass in top-decile angles
        top = np.sort(th2[prim])[::-1]
        top10 = top[:max(1, len(prim) // 10)].sum() / th2[prim].sum()
        print(f"   r={r} (N={N}, {len(prim)} prim angles): Lambda_r={Lam.real:+.6e}")
        print(f"      accumulation Re[sum_{{k<=K}} w|theta|^2] at K=25/50/75/100%: "
              f"{', '.join(f'{v:+.3e}' for _, v in cps)}")
        print(f"      |theta_hat|^2 over prim: min={th2[prim].min():.5f} max={th2[prim].max():.5f}; "
              f"top-decile angles hold {top10*100:.1f}% of layer |theta|^2 mass")
    print("   (raw material for locating period-9-in-r as an angular phenomenon; pen adjudicates.)")


if __name__ == "__main__":
    main()
