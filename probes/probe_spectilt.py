"""
PROBE SPECTILT -- where the non-3-divisible power spectrum sits (Wilson's spectral identity).

Wilson's exact identity (from the replica-peak variables):
  Re dhat_r(1) = [ sum_{3 nmid m} Chat(m) cos(2pi m/N) ] / [ sum_{3 nmid m} Chat(m) ],  Chat(m)=|rho_hat(m)|^2 = prof[m] >= 0.
The dominant mode is the cos-weighted MEAN of the nonnegative fluctuation spectrum (restricted to 3 nmid m),
normalized to a probability measure. Automatically in [-1,1]; positive iff the spectrum is low-frequency weighted.
  denominator = C(0)-C(N/3) (x2) = total fluctuation spectral mass; its vanishing = forced-uniformity (R14), not cancellation.

No build_nu. Item 1 + consistency from the GATED AC-B normalized lags (r=2..16). Items 2,3 from the saved spectra
rho_12..16 (scratchpad). Wilson's flag: the SIGN obstruction (Chat>=0 constrains MASS not PLACEMENT) is unchanged;
what changed is placement is now ONE measurable distribution.

Item 1  per-level denominator C(N/3)/C(0) and its rate; re-read the numerator (margin) rate against it (identity check).
Item 2  cumulative distribution of sum_{3nmid m} Chat(m) vs |m|=min(m,N-m), deciles, r=12..16: mass vs N/4, migrating?
Item 3  numerator by band: |m|<N/4 (cos>0, +), |m|>N/4 (cos<0, -); is the +margin a residue of large opposing bands
        or a genuine small tilt of a spread measure?  (two completely different proof strategies)
"""
import os, sys, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"

# gated AC-B normalized lags (C./C0), r -> (C1, Cn3, C(n3-1), C(n3+1))
ACB = {
    2: (0.323810, 0.666667, 0.228571, 0.380952), 3: (0.269936, 0.734160, 0.310548, 0.218782),
    4: (0.230408, 0.773085, 0.201215, 0.255679), 5: (0.201118, 0.802422, 0.218242, 0.180818),
    6: (0.178584, 0.825201, 0.165052, 0.189408), 7: (0.160572, 0.843649, 0.167882, 0.151454),
    8: (0.145843, 0.858289, 0.139625, 0.150739), 9: (0.133573, 0.870334, 0.136892, 0.129244),
    10: (0.123216, 0.880480, 0.120165, 0.125407), 11: (0.114346, 0.889129, 0.115806, 0.112174),
    12: (0.106668, 0.896601, 0.105103, 0.107622), 13: (0.099955, 0.903139, 0.100578, 0.098809),
    14: (0.094032, 0.908898, 0.093188, 0.094431), 15: (0.088767, 0.914016, 0.089015, 0.088143),
    16: (0.084055, 0.918592, 0.083589, 0.084205),
}


def main():
    t0 = time.time()
    print("# PROBE SPECTILT -- placement of the 3-nmid power spectrum (Wilson's identity)\n")

    # ================= ITEM 1  denominator rate + numerator(margin) rate + identity check =================
    print("## ITEM 1  Den = C(0)-C(N/3), Num = 2C(1)-C(N/3-1)-C(N/3+1), d1 = Num/(2 Den)  [all /C0]")
    print(f"   {'r':>2} {'Den/C0':>10} {'DenRate':>8} {'Num/C0':>12} {'NumRate':>8} {'d1':>12} {'d1Rate':>8} {'d1*2Den vs Num':>15}")
    Den = {}; Num = {}; d1 = {}
    for r in range(2, 17):
        C1, Cn, Cm, Cp = ACB[r]
        Den[r] = 1.0 - Cn
        Num[r] = 2 * C1 - Cm - Cp
        d1[r] = Num[r] / (2 * Den[r])
    for r in range(2, 17):
        dR = Den[r] / Den[r - 1] if r > 2 else float('nan')
        nR = Num[r] / Num[r - 1] if r > 2 else float('nan')
        eR = d1[r] / d1[r - 1] if r > 2 else float('nan')
        chk = 2 * d1[r] * Den[r]                        # must equal Num[r]
        print(f"   {r:>2} {Den[r]:>10.6f} {dR:>8.4f} {Num[r]:>12.4e} {nR:>8.4f} {d1[r]:>12.4e} {eR:>8.4f} "
              f"{chk:>10.4e}={'OK' if abs(chk-Num[r])<1e-9 else 'X'}")
    def gm(seq, lo, hi):
        return (seq[hi] / seq[lo]) ** (1.0 / (hi - lo))
    print(f"\n   full-span (r2->16) geo rates: Den={gm(Den,2,16):.4f}  Num={gm(Num,2,16):.4f}  d1={gm(d1,2,16):.4f}")
    print(f"   late (r12->16)     geo rates: Den={gm(Den,12,16):.4f}  Num={gm(Num,12,16):.4f}  d1={gm(d1,12,16):.4f}")
    print(f"   identity says NumRate = d1Rate * DenRate; late: {gm(d1,12,16):.4f}*{gm(Den,12,16):.4f} = "
          f"{gm(d1,12,16)*gm(Den,12,16):.4f}  vs measured Num {gm(Num,12,16):.4f}")
    print("   [Den rate SLOWS (forced-uniformity easing); the margin/Num rate is what it must be by the identity.]\n")

    # ================= ITEMS 2,3  spectral placement from saved rho =================
    print("## ITEMS 2,3  the nonnegative fluctuation spectrum Chat(m)=|rho_hat(m)|^2, 3 nmid m\n")
    for r in (12, 13, 14, 15, 16):
        N = 3 ** r
        rho = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy"))
        prof = np.abs(np.fft.fft(rho)) ** 2                       # Chat(m)
        del rho
        m = np.arange(N)
        prim = (m % 3 != 0)                                       # 3 nmid m (excludes m=0)
        pm = prof[prim]; mm = m[prim]
        fm = np.minimum(mm, N - mm)                               # |m| folded into [1, N/2]
        cosw = np.cos(2 * np.pi * mm / N)
        S = pm.sum()
        Num_spec = float((pm * cosw).sum())
        d1_spec = Num_spec / S
        # gate: spectral d1 vs banked
        rel = abs(d1_spec - d1[r]) / abs(d1[r])
        # item 3: bands by sign of cos  (|m|<N/4 -> cos>0 ; |m|>N/4 -> cos<0)
        pos = cosw > 0
        pcon = float((pm[pos] * cosw[pos]).sum())                # positive-band numerator contribution
        ncon = float((pm[~pos] * cosw[~pos]).sum())              # negative-band (<=0)
        mpos = float(pm[pos].sum()) / S                          # mass fraction |m|<N/4
        mneg = float(pm[~pos].sum()) / S
        # item 2: deciles of mass vs |m|  (weighted histogram over [0,N/2])
        NB = 4000
        hist, edges = np.histogram(fm, bins=NB, range=(0, N / 2), weights=pm)
        cum = np.cumsum(hist) / S
        dec_frac = []
        for q in (0.1, 0.25, 0.5, 0.75, 0.9):
            idx = int(np.searchsorted(cum, q))
            mval = edges[min(idx, NB)]
            dec_frac.append(mval / N)                            # |m|/N  (N/4 -> 0.25)
        print(f"   r={r} (N={N}):  d1_spec={d1_spec:+.6e}  (vs banked {d1[r]:+.6e}, rel {rel:.1e})")
        print(f"      MASS: frac(|m|<N/4)={mpos:.4f}  frac(|m|>N/4)={mneg:.4f}   <cos>=Num/S={d1_spec:+.5e}")
        print(f"      BANDS(Num): pos(|m|<N/4)={pcon/S:+.5e}  neg(|m|>N/4)={ncon/S:+.5e}  "
              f"sum={ (pcon+ncon)/S:+.5e}  |residue|/|pos| = {abs(pcon+ncon)/abs(pcon):.4f}")
        print(f"      DECILES |m|/N (N/4=0.25): 10%={dec_frac[0]:.3f} 25%={dec_frac[1]:.3f} 50%={dec_frac[2]:.3f} "
              f"75%={dec_frac[3]:.3f} 90%={dec_frac[4]:.3f}   ({time.time()-t0:.1f}s)")
        del prof, pm, mm, fm, cosw
    print("\n   [item2 read: are deciles (|m|/N) MIGRATING outward r12->16 or stable-with-mass-shrinking?]")
    print("   [item3 read: |residue|/|pos| ~1 => genuine small tilt of spread measure; <<1 => small residue of large opposing bands.]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
