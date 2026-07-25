"""
PROBE MODES -- the k=1 vs k=2 (vs k=3 fiber-mean) mode decomposition of the coupling.

g_r = <delta_r, Re w> = sum_{n>=1} 4^{-n} Re dhat_r(n),  dhat_r(n) = sum_k delta_r(k) e^{-2pi i n k / 3^r}
   = lag-n autocorrelation of nu in the orbit coordinate (delta_r = |nu_hat|^2/S - uniform).
Re w weights fall x1/4 per mode; 3|n modes are the fiber-mean, 3-nmid n the fiber-fluctuation.
 - fluctuation (3-nmid n): dominated by n=1 (weight 1/4) and n=2 (1/16); n=4 and up < 3%.
 - fiber-mean (3|n): n=3 (weight 1/64), ...

Three checks (Wilson):
 1. CLOSURE: 4^{-1}Re dhat(1) + 4^{-2}Re dhat(2) reproduces the FLUCTUATION coupling to ~3%.
 2. THE CROSSINGS (r=2, r=6): is Re dhat(1) itself changing sign, or is n=1 positive throughout with
    n=2 (or the n=3 fiber-mean) momentarily overpowering it? -- THE decisive question.
 3. RATE: does Re dhat(1) decay at ~0.9, and does Re dhat(2)/Re dhat(1) drift?
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu

Rew = lambda x: 15.0 / (2 * (17 - 8 * np.cos(2 * np.pi * x))) - 0.5
JMAX = 12


def delta_from_nu(nuj, j):
    N = 3 ** j
    mu = np.zeros(N)
    for X, w in nuj.items():
        mu[(X - 1) // 3 % N] += float(w)
    d = R10.dlog_table(j)
    rho = np.zeros(N)
    np.add.at(rho, np.array([d[a] for a in range(N)]), mu)
    prof = np.abs(np.fft.fft(rho)) ** 2
    prim = np.array([k for k in range(1, N) if k % 3 != 0]); M = len(prim)
    S = float(prof[prim].sum())
    dd = np.zeros(N); dd[prim] = prof[prim] / S - 1.0 / M
    return dd, S


def main():
    t0 = time.time()
    print(f"# PROBE MODES -- k=1/k=2/k=3 decomposition of the coupling, r=2..{JMAX}.\n")
    nus = build_nu(0.5, JMAX)
    dlt = {}
    for j in range(2, JMAX + 1):
        dlt[j], _ = delta_from_nu(nus[j], j)
    print(f"  delta built ({time.time()-t0:.1f}s)\n")

    dh = {}
    for r in range(2, JMAX + 1):
        N = 3 ** r
        H = np.fft.fft(dlt[r])                         # dhat(n) = sum_k delta(k) e^{-2pi i n k/N}
        dh[r] = {n: H[n].real for n in (1, 2, 3, 4, 6)}

    # ---- CHECK 1: closure ----
    print("## CHECK 1  CLOSURE: 4^-1 Re dhat(1)+4^-2 Re dhat(2) vs fluctuation coupling (3-nmid n) & g_r")
    print(f"   {'r':>2} {'g_r':>12} {'fluct(3!|n)':>12} {'1/4 d1 +1/16 d2':>16} {'fibermean(3|n)':>14}")
    for r in range(2, JMAX + 1):
        N = 3 ** r
        H = np.fft.fft(dlt[r]).real
        n = np.arange(N)
        w = np.where(n >= 1, 4.0 ** (-np.minimum(n, N - n).astype(float)), 0.0)   # 4^-|n| for n>=1 (fold)
        # careful: use n from 1..N/2 with 4^-n; contributions from n and N-n both ~ same (real)
        full = float(np.sum([4.0 ** (-m) * H[m] for m in range(1, min(40, N))]))  # g_r ~ sum 4^-n Re dhat(n)
        fluct = float(np.sum([4.0 ** (-m) * H[m] for m in range(1, min(40, N)) if m % 3 != 0]))
        fibm = float(np.sum([4.0 ** (-m) * H[m] for m in range(1, min(40, N)) if m % 3 == 0]))
        two = 0.25 * H[1] + (1.0 / 16) * H[2]
        gr = float(np.sum(dlt[r] * Rew(np.arange(N) / N)))
        print(f"   {r:>2} {gr:>+12.4e} {fluct:>+12.4e} {two:>+16.4e} {fibm:>+14.4e}")
    print("   [full sum_{n>=1}4^-n Re dhat(n) should = g_r; fluct=3-nmid part; two-mode approx vs fluct to ~3%.]\n")

    # ---- CHECK 2: the crossings ----
    print("## CHECK 2  THE CROSSINGS: sign of Re dhat(1), (2), (3) per r -- does n=1 flip, or is it overpowered?")
    print(f"   {'r':>2} {'Re dhat(1)':>13} {'Re dhat(2)':>13} {'Re dhat(3)':>13} {'sgn d1':>7} {'sgn d2':>7} {'sgn d3':>7} {'sgn g_r':>8}")
    for r in range(2, JMAX + 1):
        d1, d2, d3 = dh[r][1], dh[r][2], dh[r][3]
        N = 3 ** r
        gr = float(np.sum(dlt[r] * Rew(np.arange(N) / N)))
        print(f"   {r:>2} {d1:>+13.5e} {d2:>+13.5e} {d3:>+13.5e} {'+' if d1>0 else '-':>7} {'+' if d2>0 else '-':>7} "
              f"{'+' if d3>0 else '-':>7} {'+' if gr>0 else '-':>8}")
    print("   [if Re dhat(1) flips at r=2 & r=6 => the eventuality is unproved. if n=1 stays + and crossings are")
    print("    n=2/n=3 excursions => single-signed dominant mode + bounded perturbation; question = |d2|<4|d1|?]\n")

    # ---- CHECK 3: rates ----
    print("## CHECK 3  RATES: Re dhat(1) decay, Re dhat(2)/Re dhat(1) drift")
    print("   Re dhat(1) ratio r/(r-1):")
    print("   " + " ".join(f"{dh[r][1]/dh[r-1][1]:+.3f}" for r in range(4, JMAX + 1) if abs(dh[r-1][1]) > 1e-15))
    print("   Re dhat(2)/Re dhat(1) per r:")
    print("   " + " ".join(f"r{r}:{dh[r][2]/dh[r][1]:+.3f}" for r in range(2, JMAX + 1) if abs(dh[r][1]) > 1e-15))
    print("   |Re dhat(2)|/|Re dhat(1)| vs the 4x threshold (fluct sign safe iff <4):")
    print("   " + " ".join(f"r{r}:{abs(dh[r][2])/abs(dh[r][1]):.2f}" for r in range(2, JMAX + 1) if abs(dh[r][1]) > 1e-15))
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
