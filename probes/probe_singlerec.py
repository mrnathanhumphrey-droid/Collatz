"""
PROBE SINGLEREC -- Wilson's single-level recursion + the measure-free P_k(xi) (2026-07-26).

INDEPENDENCE FORK (resolved by reading stationary_trunc): it transitions r -> (3r+1)*2^-v with weight (0.5^v)/Z
INDEPENDENT of r (NOT the true valuation v_2(3r+1)) = the Syracuse-random-variable model = Wilson's assumption. Match.

SINGLE-LEVEL RECURSION (no level change):
    pi-hat_k(xi) = Sum_{a>=1} 2^-a e(xi 2^-a / 3^k) pi-hat_k(3 xi 2^-a mod 3^k).
TRAP: weight 2^-a is REAL (0.5^a); the 2^-a in the phase AND the frequency map is inv(2^a) mod 3^k. Same symbol,
different objects. Freq map xi -> 3 xi 2^-a raises v_3 by 1 => terminates in <=k steps (pi-hat(0)=1) => FINITE expansion.

MEASURE-FREE CARRIER:  P_k(xi) = |Sum_{a>=1} 2^-a e(xi 2^-a/3^k)| < 1 for xi != 0 (deterministic, no nu).
Escape/contraction: |pi-hat_k(xi)| <= |c| P_k(xi) + V, V=Sum_a 2^-a|c_a-c|, c_a=pi-hat(3 xi 2^-a). Contracts if V < S(1-P).
Tabulate P_k over 3-nmid xi -- the room 1-P at the frequencies that matter.

Reuses stationary_trunc. Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import gcd
import numpy as np
from probe_27_high_k_rho_q5 import stationary_trunc


def fwd_hat(k):
    N = 3 ** k
    pi, _ = stationary_trunc(3, k)
    cp = np.array([r for r in range(N) if gcd(r, 3) == 1], dtype=np.int64)
    dense = np.zeros(N); dense[cp] = pi; dense /= dense.sum()
    return np.conj(np.fft.fft(dense)), N          # pi-hat(xi)=E[e(+2pi i xi X/N)] = conj(numpy fft, which is e^-)


def main():
    t0 = time.time()
    print("# PROBE SINGLEREC -- single-level recursion + measure-free P_k(xi)\n")
    print("## independence fork: stationary_trunc weight = (0.5^v)/Z INDEPENDENT of source r (not v_2(3r+1))")
    print("   => Syracuse-random-variable model = Wilson's assumption. MATCH.\n")

    print("## SINGLE-LEVEL RECURSION gate: pi_k(xi) =?= Sum_a 2^-a e(xi 2^-a/3^k) pi_k(3 xi 2^-a mod 3^k)")
    print("   (weight 0.5^a REAL; phase/freq 2^-a = inv(2^a) mod 3^k)")
    for k in range(3, 7):
        ph, N = fwd_hat(k)
        inv2 = pow(2, -1, N)
        prim = [xi for xi in range(1, N) if xi % 3 != 0]
        errs = []
        for xi in prim:
            rhs = 0j; p = 1
            for a in range(1, 64):                        # real-weight tail 0.5^a -> machine precision by a~54
                p = (p * inv2) % N                        # inv(2^a) mod N
                w = xi * p % N                            # xi*2^-a mod 3^k
                rhs += (0.5 ** a) * np.exp(2j * np.pi * w / N) * ph[(3 * w) % N]
            # real-weight tail beyond 3k+1 is 2^-(3k+1) ~ negligible
            errs.append(abs(rhs - ph[xi]))
        rel = np.mean(errs) / np.mean([abs(ph[xi]) for xi in prim])
        print(f"   k={k}: rel = {rel:.3e}  [{'REPRODUCES' if rel<1e-9 else ('close' if rel<1e-4 else 'NO')}]")
    print()

    print("## MEASURE-FREE CARRIER P_k(xi) = |Sum_a 2^-a e(xi 2^-a/3^k)|  over 3-nmid xi  (deterministic, no nu)")
    print(f"   {'k':>2} {'max P':>8} {'min 1-P':>9} {'median P':>9} {'P at xi=1':>10} {'P at xi=2':>10}")
    for k in range(3, 9):
        N = 3 ** k
        inv2 = pow(2, -1, N)
        # P_k(xi) for all 3-nmid xi
        xis = np.array([xi for xi in range(1, N) if xi % 3 != 0], dtype=np.int64)
        # accumulate sum_a 0.5^a e(xi 2^-a/N)
        P = np.zeros(len(xis), dtype=complex)
        p = 1
        for a in range(1, 3 * k + 40):
            p = (p * inv2) % N
            P += (0.5 ** a) * np.exp(2j * np.pi * (xis * p % N) / N)
        Pabs = np.abs(P)
        # P at xi=1,2
        def Pxi(xi):
            s = 0j; q = 1
            for a in range(1, 3 * k + 40):
                q = (q * inv2) % N
                s += (0.5 ** a) * np.exp(2j * np.pi * (xi * q % N) / N)
            return abs(s)
        print(f"   {k:>2} {Pabs.max():>8.5f} {1-Pabs.max():>9.5f} {np.median(Pabs):>9.5f} "
              f"{Pxi(1):>10.5f} {Pxi(2):>10.5f}")
    print("   [1-P = the deterministic room; contraction needs V < S(1-P). max P closest to 1 = tightest frequencies.]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
