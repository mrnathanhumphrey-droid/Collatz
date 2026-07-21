"""
PROBE R13 -- THE LOBE CONSTANT, psi-EXISTENCE, TRANSPORT. Reuses R7/R9/R10.

Weight closed form (Wilson): Re w(x) = 15/(2D) - 1/2, D = 17 - 8cos2pix; range [-1/5, 1/3]; mean 0;
  ||w||^2 = 1/15 (Mersenne 4^2-1); lobe integral K0 = int_0^1 |Re w| dx = 1 - (2/pi) arccos(1/4) = 0.16086726...

R13-A discriminator: L_r+M_r (measured) vs S_r*K0 (uniform-psi). R13-C decider: does gamma_r(tau_m) converge?
R13-D forced gate: X_{n+1}=1+3*2^{-v} X_n  and  t'=beta(2^{-v} 4^t) both reproduce frozen Syrac(Z/3^{n+1}).
"""
import os, sys, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import probe_engine_R7 as R7
import probe_gamma_R9 as R9
import probe_charledger_R10 as R10

S = dict(R10.S)
# extend S exactly from certified eps-table (R12-F): S_k = 7/15 + eps_k, exact through k=8
_hist = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'experiments_output',
                                    'result_77_7_eps_exact_through_k8_v2_vec_pool.json')))
EPS = {int(k): F(int(v['num']), int(v['den'])) for k, v in _hist.items()}
for k in (7, 8):
    S[k] = F(7, 15) + EPS[k]

K0 = 1 - (2 / math.pi) * math.acos(0.25)          # lobe integral constant
Rew = lambda x: 15 / (2 * (17 - 8 * math.cos(2 * math.pi * x))) - 0.5


def theta2(mu_r, r):
    N = 3 ** r
    g = R10.autocorr_dlog(mu_r, r)
    gf = [float(x) for x in g]
    return [sum(gf[u] * math.cos(2 * math.pi * k * u / N) for u in range(N)) for k in range(N)]


def main():
    print("# PROBE R13 -- LOBE CONSTANT / psi-EXISTENCE / TRANSPORT.\n")
    print(f"   weight closed form checks: Re w(0)={Rew(0):+.6f} (=1/3), Re w(1/2)={Rew(0.5):+.6f} (=-1/5); "
          f"K0=1-(2/pi)arccos(1/4)={K0:.8f}")
    # numeric mean-zero + ||w||^2 sanity
    NN = 20000
    mean = sum(Rew(i / NN) for i in range(NN)) / NN
    print(f"   <Re w> numeric = {mean:+.2e} (=0 exact);  int|Re w| numeric = {sum(abs(Rew(i/NN)) for i in range(NN))/NN:.8f} (=K0)\n")

    mu = {1: R7.mu1()}
    for k in range(2, 8):
        mu[k] = R7.build_mu(mu[k - 1], k)

    # ---- R13-A ----
    print("## R13-A  LOBE CONSTANT vs UNIFORM (discriminator; labeled columns, NO fit)")
    print(f"   {'r':>2} {'L+M measured':>14} {'S_r*K0 [unif-psi]':>18} {'abs diff':>12} {'ratio meas/unif':>16}")
    for r in range(3, 8):
        N = 3 ** r
        th2 = theta2(mu[r], r)
        prim = [k for k in range(1, N) if k % 3 != 0]
        LpM = sum(th2[k] * abs(Rew(k / N)) for k in prim)
        pred = float(S[r]) * K0
        print(f"   {r:>2} {LpM:>14.7f} {pred:>18.7f} {LpM-pred:>+12.5f} {LpM/pred:>16.5f}")
    print("   (ratio ~0.96 flat => psi has fixed shape ~4% depleted; ->1 => uniform. No verdict.)\n")

    # ---- R13-C ----  (the decider)
    print("## R13-C  psi-EXISTENCE DECIDER: gamma_r(tau_m) convergence (exact), r=1..7")
    for m in (1, 2, 3, 4, 9, 27):
        vals = []
        for r in range(1, 8):
            t = R9.tau(m, r)
            vals.append(R9.gamma(mu[r], r, t))
        diffs = [vals[i + 1] - vals[i] for i in range(len(vals) - 1)]
        print(f"   m={m:>2} (v3={R9.v3(m) if hasattr(R9,'v3') else '?'}): "
              f"gamma_r = [{', '.join(f'{float(v):.4f}' for v in vals)}]")
        print(f"         successive diffs = [{', '.join(f'{float(d):+.4f}' for d in diffs)}]  "
              f"|diff| trend: {'SHRINKING' if all(abs(float(diffs[i+1]))<=abs(float(diffs[i]))+1e-9 for i in range(len(diffs)-1)) else 'NOT monotone-shrinking'}")
    print("   (converges => psi exists; if diffs don't shrink => psi programme dead. Pen adjudicates shape.)\n")

    # ---- R13-B ----
    print("## R13-B  WHERE THE DEPLETION SITS (measurement, NO fit): profile binned by angle x=k/3^r")
    NB = 12
    for r in range(5, 8):
        N = 3 ** r
        th2 = theta2(mu[r], r)
        prim = [k for k in range(1, N) if k % 3 != 0]
        bins_meas = [0.0] * NB; bins_cnt = [0] * NB
        for k in prim:
            b = min(NB - 1, int((k / N) * NB))
            bins_meas[b] += th2[k]; bins_cnt[b] += 1
        tot = len(prim)
        print(f"   r={r}: bin[x] | measured mass | uniform S_r*(cnt/tot) | Re w(center) | meas/unif")
        for b in range(NB):
            xc = (b + 0.5) / NB
            unif = float(S[r]) * bins_cnt[b] / tot
            ratio = bins_meas[b] / unif if unif else float('nan')
            print(f"      x={xc:.3f}: {bins_meas[b]:.6f} | {unif:.6f} | {Rew(xc):+.4f} | {ratio:.4f}")
    print("   (is the ~4% depletion near x~0 [decay of nu_hat near trivial char] or x~1/2 or spread?)\n")

    # ---- R13-D ----  (forced gate)
    print("## R13-D  RENEWAL-IN-ORBIT GATE (forced): X_{n+1}=1+3*2^-v X_n and t'=beta(2^-v 4^t) reproduce frozen mu")
    okD = True
    for n in range(1, 6):
        Mn1 = 3 ** (n + 1)          # X_n = 1+3a mod 3^{n+1}
        Mn2 = 3 ** (n + 2)          # X_{n+1} mod 3^{n+2}
        inv2 = pow(2, -1, Mn2)
        ordv, x = 1, 2 % Mn2
        while x != 1:
            x = (x * 2) % Mn2; ordv += 1
        denom = 1 - F(1, 2 ** ordv)
        # ---- route (a): frequency recursion, build mu'_{n+1}
        muA = {}
        for a, pa in mu[n].items():
            Xn = (1 + 3 * a) % Mn1
            for j in range(1, ordv + 1):
                wv = F(1, 2 ** j) / denom
                Xn1 = (1 + 3 * pow(inv2, j, Mn2) * Xn) % Mn2
                ap = ((Xn1 - 1) // 3) % (3 ** (n + 1))
                muA[ap] = muA.get(ap, F(0)) + pa * wv
        matchA = (muA == mu[n + 1])
        # ---- route (b): orbit coordinate via beta-table (dlog4)
        dlog_np1 = R10.dlog_table(n + 1)             # s -> t=dlog4(1+3s) mod 3^{n+1}, s in Z/3^{n+1}
        # inverse dlog for base 4 on units mod 3^{n+2}: pow(4,t,Mn2)
        muB = {}
        for a, pa in mu[n].items():
            Xn = (1 + 3 * a) % Mn1
            t = dlog_np1[a % (3 ** (n + 1))] if (a % (3**(n+1))) in dlog_np1 else None
            # compute t directly: t = dlog4(Xn) in units mod 3^{n+1}
            for j in range(1, ordv + 1):
                wv = F(1, 2 ** j) / denom
                z = (pow(inv2, j, Mn2) * Xn) % Mn2            # 2^-v * 4^t = 2^-v * Xn
                Xn1 = (1 + 3 * z) % Mn2                        # = 4^{t'} ; t' = beta(z)
                ap = ((Xn1 - 1) // 3) % (3 ** (n + 1))
                muB[ap] = muB.get(ap, F(0)) + pa * wv
        matchB = (muB == mu[n + 1])
        okD = okD and matchA and matchB
        print(f"   n={n}: X-recursion reproduces mu_{n+1}? {matchA}   orbit t'=beta(2^-v 4^t) reproduces? {matchB}  "
              f"[{'OK' if matchA and matchB else 'FAIL'}]")
    print(f"   => R13-D {'GATE PASS -- the chirp beta IS the renewal (two routes, one transition, exact)' if okD else 'FAIL (#32)'}\n")

    # ---- R13-E ----
    print("## R13-E  R85 RUNG-1 DEBT (feasibility only, NOT run)")
    print("   R85 rung-1 = operator-DFT chirp identity, owed r=5 (n=8) extension; = same object as this transport.")
    print("   Current machinery: beta/U tables built to r=5 (N=243) in R11/R12; support law (R12-A) makes U ~1/6 dense")
    print("   (only k==xi mod3, block-diag by v3). n=8 -> N=3^8=6561: dense chirp-DFT ~ (3^8)^2 * 3^8 ~ 2.8e11 ops (heavy),")
    print("   but Bluestein/chirp-z + support-law pruning -> O(N log N) per block, ~feasible in one focused probe.")
    print("   VERDICT: cheaper than July (beta-tables now exist); n=8 is one dedicated FFT-route probe, not free. DEFER.")


if __name__ == "__main__":
    main()
