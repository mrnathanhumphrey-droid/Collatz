"""
PROBE P4 (Wilson) -- which frequencies carry gamma_inf(1)? (2026-07-26). FREE, data in hand.

gamma_r(k)-1 = Sum_{a!=0} |rho-hat(a)|^2 e(a k/3^r)  is LINEAR + EXPLICIT (Wiener-Khinchin on the CERTIFIED
gamma_r(k)=3^r <rho, shift_k rho>). No perturbation: the per-frequency influence on channel k IS the summand.
rho = cached dlog profile (CERTIFIED, scratchpad/rho_r{r}.npy, normalized sum 1). rho-hat = fft(rho).
Contribution of the pair {a, N-a} to gamma_r(1)-1 is  2|rho-hat(a)|^2 cos(2pi a/3^r)  (real; N-a is the conjugate).

GATE: 1 + Sum_{a!=0}|rho-hat|^2 e(+2pi i a/N) must reproduce 3^r <rho,roll(rho,-1)> to machine precision (sign conv).
Then sort the pair-contributions and read off WHICH frequencies carry gamma_inf(1)~0.733 (i.e. gamma-1~-0.267).
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"


def v3(n):
    n = int(n); v = 0
    while n and n % 3 == 0:
        n //= 3; v += 1
    return v


def main():
    t0 = time.time()
    print("# PROBE P4 -- which frequencies carry gamma_inf(1)? (linear summand read-off)\n")
    for r in (12, 14, 16):
        N = 3 ** r
        rho = np.load(os.path.join(SCRATCH, f"rho_r{r}.npy")).astype(np.float64)
        rho = rho / rho.sum()
        # certified gamma_r(1)
        g1_cert = 3.0 ** r * float(np.dot(rho, np.roll(rho, -1)))
        rhat = np.fft.fft(rho)
        P = (rhat.real ** 2 + rhat.imag ** 2)                 # |rho-hat(a)|^2
        a = np.arange(N)
        phase1 = np.cos(2 * np.pi * a / N)                    # real part of e(+a/N); imag cancels by symmetry
        contrib = P * phase1                                  # per-frequency real contribution to gamma_r(1)-1
        g1_fft = 1.0 + contrib[1:].sum()                      # a!=0
        # gate
        ok = abs(g1_fft - g1_cert) < 1e-6

        # pair up a and N-a: signed pair contribution = 2 P[a] cos(2pi a/N), a=1..(N-1)/2
        half = (N - 1) // 2
        aa = np.arange(1, half + 1)
        paircon = 2.0 * P[aa] * np.cos(2 * np.pi * aa / N)
        order = np.argsort(-np.abs(paircon))                  # by |contribution|
        total = g1_cert - 1.0                                 # = Sum pair contributions

        print(f"## r={r}: gamma_{r}(1) cert={g1_cert:.9f} fft={g1_fft:.9f} [{'GATE OK' if ok else 'MISMATCH'}]  "
              f"gamma-1 = {total:+.6f}")
        # cumulative: how many top pairs to reach 50/90/99% of |total|
        cum = np.cumsum(paircon[order])
        fr = cum / total
        for target in (0.5, 0.9, 0.99):
            idx = np.searchsorted(fr, target)
            print(f"     top {idx+1} pairs ({100*(idx+1)/half:.4f}% of freqs) reach {target*100:.0f}% of gamma-1")
        # top carriers structure
        print(f"     {'rank':>4} {'a':>10} {'a/N':>8} {'v3(a)':>5} {'cos':>7} {'|rhat|^2':>10} {'paircontrib':>12} {'cum/total':>9}")
        for i in range(10):
            j = order[i]; av = int(aa[j])
            print(f"     {i+1:>4} {av:>10} {av/N:>8.4f} {v3(av):>5} {np.cos(2*np.pi*av/N):>7.3f} "
                  f"{P[aa[j]]:>10.3e} {paircon[j]:>12.3e} {fr[i]:>9.4f}")
        # max-|rho-hat| frequency (the 'sup') -- does it carry?
        amax = 1 + int(np.argmax(P[1:]))
        print(f"     SUP freq a*={amax} (a/N={amax/N:.4f}, v3={v3(amax)}): |rhat|^2={P[amax]:.3e}, "
              f"contribution 2P cos = {2*P[amax]*np.cos(2*np.pi*amax/N):+.3e}  "
              f"({100*2*P[amax]*np.cos(2*np.pi*amax/N)/total:.4f}% of gamma-1)")
        # structure summary of the top 1% carriers
        top1 = order[:max(1, half // 100)]
        av1 = aa[top1]
        v3s = np.array([v3(int(x)) for x in av1[:2000]])      # sample v3 on up to 2000
        print(f"     top-1% carriers: median a/N={np.median(av1/N):.4f}, "
              f"v3 profile (sample) v3=0:{(v3s==0).mean()*100:.0f}% v3=1:{(v3s==1).mean()*100:.0f}% "
              f"v3>=2:{(v3s>=2).mean()*100:.0f}%\n")
    print(f"# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
