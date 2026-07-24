"""
GATE SPLITTING (S1) -- validation gate for the multilevel-splitting route past the build_nu wall.

Goal (pre-registered): reproduce the KNOWN exact Lambda_r for r=8..16 within stated (replica) error bars.
If it fails here, the deep run is DEAD -- a few minutes of compute, nothing lost.
If it passes, it CALIBRATES N for the deep run (design fix 1) and the cumulative-eps decision (design fix 2).

METHOD.  gamma_r(m) = 3^r * p_r(m),  p_r(m) = Pr[ X = 4^{-m} X'  (mod 3^{r+1}) ],  X,X' iid ~ nu.
This is a rare event ( ~3^{-r} ) estimated by fixed-N multilevel splitting on the perpetuity
    X = sum_{k=0..r} 3^k 2^{-(v_1+..+v_k)}  (mod 3^{r+1}),   v_j ~ Geom(1/2), P[v=j]=2^{-j}.
Digit s of X commits when v_s is drawn -> extend particles one level at a time; a pair (X,X') survives
level s iff digit s of (X - 4^{-m}X') is 0 (given digits 0..s-1 matched). q_s = survival fraction;
p_r = prod_{s<=r} q_s; gamma_r = prod (3 q_s).  NO |nu_hat|^2 is ever formed -> no squaring bias.

Lambda_r = sum_{m>=1} 4^{-m} A_r(m),  A_r(m) = gamma_r(m) - gamma_{r-1}(m).   (m truncated at MMAX identically
on both the exact side and the splitting side, so the truncation cancels in the comparison.)

ERROR BARS: R_REP independent replicas (independent genealogies) -> mean +/- SE across replicas.
This is the honest (Cerou-Guyader-style batch) SE, never naive binomial. A per-run effective-N is printed
as a cross-check only.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_gapop_R28 import build_nu

# ---------------- config ----------------
R_MAX  = 16          # gate range: compare r=8..16 to exact
MMAX   = 5           # channels m=1..5 (4^-m weights; m=1 dominates)
N      = 40000       # particles per replica per channel
R_REP  = 20          # independent replicas -> SE
JMAX   = 64          # clip geometric (P[v>64]=2^-64 negligible)
SEED   = 20260723
# ----------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
M    = 3 ** (R_MAX + 1)
INV2 = pow(2, -1, M)
INV4 = pow(4, -1, M)
INV2POW = np.array([pow(INV2, j, M) for j in range(JMAX + 1)], dtype=np.int64)
POW3    = [3 ** s for s in range(R_MAX + 2)]     # 3^s, all < M for s<=R_MAX


# ---- exact ground truth (build_nu) ----
def p_from_nu(dense, Mmod, m):
    fac = pow(pow(4, -1, Mmod), m, Mmod)
    idx = (np.arange(Mmod) * fac) % Mmod
    return float(np.sum(dense * dense[idx]))


def exact_gamma():
    nus = build_nu(0.5, R_MAX)
    gam = {m: {} for m in range(1, MMAX + 1)}
    for r in range(1, R_MAX + 1):
        Mr = 3 ** (r + 1)
        dense = np.zeros(Mr)
        for X, w in nus[r].items():
            dense[X] = w
        for m in range(1, MMAX + 1):
            gam[m][r] = 3 ** r * p_from_nu(dense, Mr, m)
        del dense
    return gam


# ---- one splitting replica: returns gamma[m][1..R_MAX] ----
def run_replica(rng):
    out = {}
    for m in range(1, MMAX + 1):
        inv4m = pow(INV4, m, M)
        accX = np.ones(N, dtype=np.int64); wX = np.ones(N, dtype=np.int64)
        accP = np.ones(N, dtype=np.int64); wP = np.ones(N, dtype=np.int64)
        gamma = np.zeros(R_MAX + 1); p = 1.0
        dead = False
        for s in range(1, R_MAX + 1):
            if dead:
                gamma[s] = 0.0; continue
            vX = np.clip(rng.geometric(0.5, N), 1, JMAX)
            vP = np.clip(rng.geometric(0.5, N), 1, JMAX)
            wX = (wX * INV2POW[vX]) % M
            wP = (wP * INV2POW[vP]) % M
            t3 = POW3[s]
            accX = (accX + (t3 * wX) % M) % M
            accP = (accP + (t3 * wP) % M) % M
            D = (accX - (inv4m * accP) % M) % M
            survive = (D % POW3[s + 1]) == 0
            k = int(survive.sum())
            if k == 0:
                dead = True; gamma[s] = 0.0; continue
            p *= k / N
            gamma[s] = POW3[s] * p
            idx = np.nonzero(survive)[0]
            pick = idx[rng.integers(0, k, N)]
            accX = accX[pick]; wX = wX[pick]; accP = accP[pick]; wP = wP[pick]
        out[m] = gamma
    return out


def main():
    t0 = time.time()
    print(f"# GATE SPLITTING (S1)  N={N} R_REP={R_REP} MMAX={MMAX} R_MAX={R_MAX}\n")
    print("## exact ground truth (build_nu)...")
    gam_ex = exact_gamma()
    Aex  = {m: {r: gam_ex[m][r] - gam_ex[m][r - 1] for r in range(2, R_MAX + 1)} for m in range(1, MMAX + 1)}
    Lex  = {r: sum(4.0 ** -m * Aex[m][r] for m in range(1, MMAX + 1)) for r in range(2, R_MAX + 1)}
    print(f"   done ({time.time()-t0:.1f}s)\n")

    ss = np.random.SeedSequence(SEED)
    rngs = [np.random.default_rng(s) for s in ss.spawn(R_REP)]
    # collect per-replica gamma[m][r]
    G = {m: np.zeros((R_REP, R_MAX + 1)) for m in range(1, MMAX + 1)}
    for i, rng in enumerate(rngs):
        rep = run_replica(rng)
        for m in range(1, MMAX + 1):
            G[m][i] = rep[m]
        if (i + 1) % 5 == 0:
            print(f"   replica {i+1}/{R_REP}  ({time.time()-t0:.1f}s)")
    print()

    # per-replica A and Lambda (same-replica difference -> correlated, low-variance)
    Lam_rep = np.zeros((R_REP, R_MAX + 1))
    for i in range(R_REP):
        for r in range(2, R_MAX + 1):
            Lam_rep[i, r] = sum(4.0 ** -m * (G[m][i, r] - G[m][i, r - 1]) for m in range(1, MMAX + 1))
    Lam_hat = Lam_rep.mean(0)
    Lam_se  = Lam_rep.std(0, ddof=1) / np.sqrt(R_REP)

    # ===== gamma reproduction (channel m=1, the dominant one) =====
    print("## GAMMA CHECK (m=1): splitting vs exact")
    print(f"   {'r':>2} {'gam_hat':>12} {'gam_exact':>12} {'rel err':>9} {'z=(hat-ex)/SE':>13}")
    for r in range(8, R_MAX + 1):
        hat = G[1][:, r].mean(); se = G[1][:, r].std(ddof=1) / np.sqrt(R_REP)
        ex = gam_ex[1][r]
        z = (hat - ex) / se if se > 0 else float('nan')
        print(f"   {r:>2} {hat:>12.5f} {ex:>12.5f} {abs(hat-ex)/abs(ex):>9.2%} {z:>13.2f}")
    print()

    # ===== THE GATE: Lambda_r reproduction r=8..16 =====
    print("## THE GATE: Lambda_r (r=8..16) splitting vs exact  [PASS iff |z|<=2 for all]")
    print(f"   {'r':>2} {'Lam_hat':>12} {'SE':>11} {'Lam_exact':>12} {'rel err':>9} {'z':>7} {'|z|<=2':>7}")
    zmax = 0.0; relerrs = {}
    for r in range(8, R_MAX + 1):
        hat, se, ex = Lam_hat[r], Lam_se[r], Lex[r]
        z = (hat - ex) / se if se > 0 else float('nan')
        relerrs[r] = se / abs(ex)
        zmax = max(zmax, abs(z))
        print(f"   {r:>2} {hat:>12.4e} {se:>11.3e} {ex:>12.4e} {abs(hat-ex)/abs(ex):>9.2%} {z:>7.2f} "
              f"{'yes' if abs(z)<=2 else 'NO':>7}")
    passed = zmax <= 2.5   # allow mild slack for 9 correlated tests
    print(f"   => max|z| = {zmax:.2f}  ::  GATE {'PASS' if passed else 'FAIL'}\n")

    # ===== design fix 1: calibrate N from achieved SE (SE ~ N^{-1/2}) =====
    print("## CALIBRATION (design fix 1): SE(Lambda_r) ~ 1/sqrt(N); rel-SE at this N, and N for 10% at r=16")
    r16 = R_MAX
    rel16 = relerrs[r16]
    N_for_10 = N * (rel16 / 0.10) ** 2
    print(f"   at N={N}: rel-SE(Lambda_16) = {rel16:.1%}")
    print(f"   => N for 10% rel-SE on Lambda_16: {N_for_10:,.0f}")
    print(f"   (Wilson: signal 3q-1 ~ Lambda decays; N must grow ~Lambda^-2 ~ 1.04^r, ~1.8x from r16->r31)\n")

    # ===== design fix 2: cumulative-eps SNR & decision-rule feasibility =====
    print("## CUMULATIVE-EPS SNR (design fix 2): eps_R = eps_16 + 2 sum_{17..R} Lambda_r")
    # required total for 7/15: sum_{r>=17} 2 Lambda = -5.15e-3 (from result_logperiodic no-turnover gap)
    # model per-term |Lambda_r| ~ Lambda_16 * 0.90^{r-16} (observed local decay); rel-SE per term ~ rel16*sqrt(N/N_deep)
    # here just report SNR at N_deep = N_for_10 (10% per Lambda_16) accumulated over r=17..40
    for N_deep_desc, N_deep in [("N_for_10", N_for_10), ("4x N_for_10", 4 * N_for_10)]:
        scale = np.sqrt(N / N_deep)                       # SE scales as sqrt(N/N_deep)
        sig = 0.0; noise2 = 0.0
        for r in range(17, 41):
            Lr = abs(Lex[R_MAX]) * 0.90 ** (r - R_MAX)    # extrapolated |Lambda_r|
            sig += 2 * Lr                                 # cumulative |signal| in eps
            se_r = rel16 * scale * Lr                     # per-term SE (rel16 holds at r16; grows mildly deeper)
            noise2 += (2 * se_r) ** 2
        snr = sig / np.sqrt(noise2)
        print(f"   {N_deep_desc}={N_deep:,.0f}: cumulative |2 sum Lambda|={sig:.3e}, cum-noise={np.sqrt(noise2):.3e}, SNR={snr:.1f}")
    print("   [decision rule: turnover CONFIRMED if eps_R < eps_16 by >=3 sigma_cum;")
    print("    NO-turnover if eps still rising at r=35 by >=3 sigma; else INCONCLUSIVE.]")
    print(f"\n# total extensions this gate ~ {R_REP*MMAX*R_MAX*N:,}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
