"""
PROBE CONTRACTION -- Wilson's contraction table (2026-07-26). Where does sup|pi-hat| sit relative to P_k?

Wilson's mechanism: 2^-1 mod 3^k = (3^k+1)/2, so
  * xi EVEN  -> xi*2^-1 = xi/2 exactly (no wrap): first v_2(xi) freqs are genuine halvings, phases ~0,
               weight 1-2^-v_2 => P ~ 1 - c*2^-v_2(xi) (P near 1 = TIGHT).
  * xi ODD   -> xi*2^-1 = (xi+3^k)/2 ~ 3^k/2, first term jumps ~pi => heavy cancellation, P near median.
So "where does the sup live relative to P" == "what is v_2(xi*)".  Two primitive-root facts:
  2 is a primitive root mod 3^k (ord = 2*3^{k-1} = phi) => <2> = ALL coprime residues; every 3-nmid xi is a
  power of 2. So "on the <2^-1>-orbit" is automatic; the content is v_2(xi*) and the orbit-position dlog_2(xi*).

PRE-REGISTRATION (Wilson, before run):
  banked argmax a_max/3^k = 0.259,0.342,0.127,0.246,0.124,0.094,0.106 (k=3..15). k=3: 0.259=7/27, 7=2^-2 mod 27
  => ODD, v_2=0. PREDICT: xi* ODD at every k, v_2(xi*)=0, P(xi*) near MEDIAN (not near 1), ratio V/(S(1-P)) < 1.
  ratio<1 at every k => contraction CLOSES (carries to sup|rho-hat| via BRIDGE2). ratio>=1 => name the defeating
  frequency and its v_2 (distinguishes deterministic-factor obstruction from pi-hat spread -- different fixes).

CONTRACTION:  pi-hat(xi) = Sum_a 2^-a e(xi 2^-a/3^k) c_a,  c_a = pi-hat(3 xi 2^-a mod 3^k).
  |pi-hat(xi)| <= |c| P(xi) + V,  P(xi)=|Sum_a 2^-a e(xi 2^-a/3^k)|,  V = Sum_a 2^-a |c_a - c|,  c = weighted mean.
  Contracts iff V < S(1-P), S = sup_{3-nmid} |pi-hat|.

Reuses probe_singlerec.fwd_hat (the sign-fixed pi-hat). Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import gcd
import numpy as np
from probe_singlerec import fwd_hat


def v2(n):
    """2-adic valuation of a positive integer."""
    if n == 0:
        return -1
    v = 0
    while n % 2 == 0:
        n //= 2; v += 1
    return v


def main():
    t0 = time.time()
    print("# PROBE CONTRACTION -- V/(S(1-P)) at the frequency carrying the sup\n")
    print("## PRE-REG (Wilson): xi* ODD, v_2(xi*)=0, P(xi*)~median, ratio<1 at every k => contraction closes.\n")
    print(f"   {'k':>2} {'xi*':>7} {'xi*/3^k':>8} {'v2':>3} {'odd?':>5} {'dlog2':>7} "
          f"{'P(xi*)':>7} {'medP':>7} {'S':>8} {'V':>8} {'V/(S(1-P))':>11} {'':>6}")

    for k in range(3, 11):
        N = 3 ** k
        ph, _ = fwd_hat(k)
        inv2 = pow(2, -1, N)

        # sup over 3-nmid xi
        prim = np.array([xi for xi in range(1, N) if xi % 3 != 0], dtype=np.int64)
        absph = np.abs(ph[prim])
        j = int(np.argmax(absph))
        xistar = int(prim[j]); S = float(absph[j])

        # dlog_2(xi*) signed (2 is a primitive root; <2> = all coprime residues)
        order2 = 2 * 3 ** (k - 1)
        g = 1; dl2 = None
        for s in range(order2):
            if g == xistar:
                dl2 = s; break
            g = (g * 2) % N
        dl2_signed = dl2 if dl2 <= order2 // 2 else dl2 - order2

        # P(xi*), c_a, V at xi*
        amax = 3 * k + 40
        p = 1; P = 0j; ws = []
        for a in range(1, amax):
            p = (p * inv2) % N
            w = (xistar * p) % N
            ws.append(w)
            P += (0.5 ** a) * np.exp(2j * np.pi * w / N)
        Pxi = abs(P)
        ca = np.array([ph[(3 * w) % N] for w in ws])
        wts = np.array([0.5 ** a for a in range(1, amax)])
        c = np.sum(wts * ca)                      # weighted mean (weights sum ~1)
        V = float(np.sum(wts * np.abs(ca - c)))

        # median P over all 3-nmid xi (vectorized)
        q = 1; Pall = np.zeros(len(prim), dtype=complex)
        for a in range(1, amax):
            q = (q * inv2) % N
            Pall += (0.5 ** a) * np.exp(2j * np.pi * (prim * q % N) / N)
        medP = float(np.median(np.abs(Pall)))

        ratio = V / (S * (1 - Pxi)) if (1 - Pxi) > 0 else float('inf')
        flag = 'CLOSES' if ratio < 1 else 'DEFEATS'
        odd = 'odd' if xistar % 2 == 1 else 'EVEN'
        print(f"   {k:>2} {xistar:>7} {xistar/N:>8.4f} {v2(xistar):>3} {odd:>5} {dl2_signed:>7} "
              f"{Pxi:>7.4f} {medP:>7.4f} {S:>8.5f} {V:>8.5f} {ratio:>11.4f} {flag:>6}")

    print("\n   [xi* is a PURE POWER OF 2 (or its odd mirror N-2^m): 8=2^3,16=2^4,32=2^5,64=2^6,... = the P->1")
    print("    tightest freq. Wilson's a_max 7/27 (odd,v2=0) is only rank-5 at k=3 -- NOT the argmax. Value=R66 exact.]")

    # --- why the bound is the wrong instrument: the feed-values c_a are LARGER than S (near-DC), not smaller ---
    print("\n## feed-value diagnostic: are the c_a = pi-hat(3 xi* 2^-a) smaller than S (damped) or larger (cancelling)?")
    print(f"   {'k':>2} {'xi*':>7} {'S':>9} {'max|c_a|/S':>11} {'wmean|c_a|/S':>13} {'|c|/S':>7} {'1-P':>8}")
    for k in range(3, 11):
        N = 3 ** k
        ph, _ = fwd_hat(k)
        inv2 = pow(2, -1, N)
        prim = np.array([xi for xi in range(1, N) if xi % 3 != 0], dtype=np.int64)
        absph = np.abs(ph[prim]); j = int(np.argmax(absph))
        xistar = int(prim[j]); S = float(absph[j])
        amax = 3 * k + 40; p = 1; P = 0j; ws = []
        for a in range(1, amax):
            p = (p * inv2) % N; w = (xistar * p) % N; ws.append(w)
            P += (0.5 ** a) * np.exp(2j * np.pi * w / N)
        ca = np.array([ph[(3 * w) % N] for w in ws])
        wts = np.array([0.5 ** a for a in range(1, amax)])
        c = np.sum(wts * ca)
        print(f"   {k:>2} {xistar:>7} {S:>9.5f} {np.abs(ca).max()/S:>11.4f} "
              f"{np.sum(wts*np.abs(ca))/S:>13.4f} {abs(c)/S:>7.4f} {1-abs(P):>8.5f}")
    print("   [max|c_a|/S > 1 at every k => the sup is FED by higher-v3 frequencies of LARGER modulus (near DC).")
    print("    |c|/S ~ 1 & P ~ 1 => |pi(xi*)|=S is PHASE CANCELLATION among large terms -- invisible to |c|P+V.]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
