"""
PROBE MAXMODE -- gate Hank's load-bearing claim: the binding channel relaxation rate = nu's MAX Fourier mode, not
the AVERAGE (which is 1/3, what every external theorem quotes). (2026-07-26)

Hank: q_j(k)-1/3 for the binding channels k=3,4 relaxes at the rate of nu's SLOWEST (max) primitive Fourier
coefficient (~0.707 amplitude, R66's ~2^-k in |rho-hat|^2), NOT the average (rate 1/3). If true, the relaxation
bound is a spectral-gap lemma on the max mode, and the decay literature (avg rate) is the wrong mode.

TEST: power spectrum ps_j(m)=|FFT(rho_j)|^2 over primitive m (3-nmid m). Track:
  P_max(j)=max, P_avg(j)=mean; per-level POWER ratio and AMPLITUDE rate (sqrt).
Compare to the channel increment rate |A_j(k)/A_{j-1}(k)| for k=1,3,4 (A_j(k)=gamma_j(k)-gamma_{j-1}(k)).
Prediction (Hank): P_max amplitude rate ~ 0.70-0.73 ~ channel k=3,4 rate; P_avg power rate ~ 1/3.
Reuses probe_channelfam rho-build. Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"


def main():
    t0 = time.time()
    print("# PROBE MAXMODE -- is the binding channel relaxation = nu's MAX Fourier mode (not the average 1/3)?\n")

    nus = build_nu(0.5, 11)
    rho = {}
    for r in range(1, 12):
        N = 3 ** r
        mu = np.zeros(N)
        for X, w in nus[r].items():
            mu[(X - 1) // 3 % N] += float(w)
        d = R10.dlog_table(r)
        rr = np.zeros(N); rr[np.fromiter((d[a] for a in range(N)), np.int64, N)] = mu
        rho[r] = rr / rr.sum()
    del nus
    for r in range(12, 17):
        rr = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy")); rho[r] = rr / rr.sum()

    # power spectrum max/avg over primitive m
    Pmax, Pavg = {}, {}
    for j in range(2, 17):
        N = 3 ** j
        ps = np.abs(np.fft.fft(rho[j])) ** 2
        m = np.arange(N)
        prim = (m % 3 != 0)
        Pmax[j] = float(ps[prim].max())
        Pavg[j] = float(ps[prim].mean())

    print("## nu power spectrum |rho-hat(m)|^2 over primitive m: max vs avg decay")
    print(f"   {'j':>2} {'P_max':>11} {'max pow-ratio':>13} {'max ampl':>9} | {'P_avg':>11} {'avg pow-ratio':>13} {'avg ampl':>9}")
    for j in range(3, 17):
        mr = Pmax[j] / Pmax[j - 1]; ar = mr ** 0.5
        vr = Pavg[j] / Pavg[j - 1]; va = vr ** 0.5
        print(f"   {j:>2} {Pmax[j]:>11.4e} {mr:>13.4f} {ar:>9.4f} | {Pavg[j]:>11.4e} {vr:>13.4f} {va:>9.4f}")
    # late-window medians
    mr_late = np.median([Pmax[j] / Pmax[j - 1] for j in range(12, 17)])
    vr_late = np.median([Pavg[j] / Pavg[j - 1] for j in range(12, 17)])
    print(f"   late (j=12..16) median: max pow-ratio {mr_late:.4f} (ampl {mr_late**0.5:.4f}) ; "
          f"avg pow-ratio {vr_late:.4f} (ampl {vr_late**0.5:.4f})")
    print(f"   [avg pow-ratio should be ~1/3=0.333 = the external-theorem rate; max is the binding one]\n")

    # channel increment rates
    def gam(r, k):
        return 3.0 ** r * float(np.dot(rho[r], np.roll(rho[r], -k)))
    print("## channel increment |A_j(k)/A_{j-1}(k)| for k=1,3,4  (A_j(k)=gamma_j-gamma_{j-1}); vs P_max amplitude")
    for k in (1, 3, 4):
        A = {j: gam(j, k) - gam(j - 1, k) for j in range(2, 17)}
        rates = [abs(A[j] / A[j - 1]) for j in range(13, 17) if abs(A[j - 1]) > 1e-18]
        print(f"   k={k}: |A_j/A_{{j-1}}| late = {['%.3f'%x for x in rates]}  median {np.median(rates):.3f}")
    print()

    print("## VERDICT")
    amp_max = mr_late ** 0.5
    print(f"   P_max amplitude rate = {amp_max:.3f} ; P_avg power rate = {vr_late:.3f} (external theorems quote ~1/3)")
    print(f"   channel k=3,4 relaxation ~0.66-0.73. Match to MAX mode: {'YES' if 0.60 <= amp_max <= 0.80 else 'NO'}; "
          f"AVG mode (0.33) would be {'too fast (wrong mode confirmed)' if vr_late < 0.45 else 'comparable'}.")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
