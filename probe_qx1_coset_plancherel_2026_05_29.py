"""
probe_qx1_coset_plancherel_2026_05_29.py

NEXT PROBE (teed up in VALIDATION_CHECKLIST_2026_05_29.md):
Does the qx+1 char-fn Plancherel mass reflect the ord_2(q) mixing degree?

Setup. qx+1 Syracuse offset over Z/q^n (generalize the q=3 builder, 3->q):
    X_n^(q) = sum_{j=1}^{n} q^{j-1} * 2^{-(a_1+...+a_j)}  mod q^n,  a iid Geom(2).
char-fn  mu_hat(xi) = E[ e(-xi X_n / q^n) ],  xi in (Z/q^n)* .

Mechanism. mu_hat resonates (large |mu_hat|) only where xi stays "aligned" with
powers of 2: xi * 2^{-S} must repeatedly land near 0 mod q^{...}, which happens iff
xi lives in the cyclic subgroup H = <2> of (Z/q^n)*. So:
  - q=3, q=5: 2 is a primitive root mod q^n  => H = all units (index 1) => NO dark coset.
  - q=7: ord_2(7^n)=phi/2 (index 2) => half the frequency space is a "dark" coset that
    powers of 2 never reach => predicted SUPPRESSED Plancherel mass there.

HYPOTHESIS: Plancherel mass fraction in H is ~1 for q=3,5 (trivially, index 1) and,
for q=7, the in-H coset dominates while the complementary coset is suppressed.
A measurable in-H >> out-of-H asymmetry = the ord_2(q) half-mixing signature.
A ~50/50 split would be the NULL.

Control: q=3,5 are index 1 (cannot show asymmetry); q=7 is the test.
"""
from __future__ import annotations
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

A_MAX = 100

def offset_distribution(q, n):
    """Exact distribution of X_n^(q) on Z/q^n (length q^n probability vector). 3->q generalization."""
    N = q ** n
    inv2 = pow(2, -1, N)
    arange = np.arange(N)
    P_U = np.zeros(N, dtype=np.float64)
    p = inv2
    for a in range(1, A_MAX + 1):
        P_U[p] += 2.0 ** (-a)
        p = (p * inv2) % N
    P_U /= P_U.sum()
    u_support = np.nonzero(P_U)[0]
    u_weight = P_U[u_support]
    P_S = P_U.copy()
    for j in range(n - 1, 0, -1):
        v_idx = (1 + q * arange) % N          # Horner step 1 + q*S
        P_V = np.zeros(N, dtype=np.float64)
        np.add.at(P_V, v_idx, P_S)
        P_new = np.zeros(N, dtype=np.float64)
        for u, w in zip(u_support, u_weight):
            perm = (u * arange) % N
            P_new[perm] += w * P_V
        P_S = P_new
    return P_S

def subgroup_pow2(q, n):
    """H = <2> in (Z/q^n)* as a boolean mask over [0,N)."""
    N = q ** n
    H = np.zeros(N, dtype=bool)
    x = 1 % N
    seen = set()
    while x not in seen:
        seen.add(x); H[x] = True
        x = (x * 2) % N
    return H, len(seen)

def main():
    print("=" * 96)
    print("qx+1 coset-resolved Plancherel: in-<2> vs complementary coset mass")
    print("=" * 96)
    print(f"{'q':>2} {'n':>2} {'N=q^n':>8} {'ord2':>6} {'phi':>7} {'idx':>4} "
          f"{'massfrac_inH':>13} {'massfrac_out':>13} {'maxH/maxOut':>12} {'argmax_inH':>10}")
    for q in (3, 5, 7):
        for n in range(2, 7):
            N = q ** n
            P = offset_distribution(q, n)
            assert abs(P.sum() - 1.0) < 1e-9
            assert P[::q].sum() < 1e-11   # supported on units
            mu = np.fft.fft(P)
            xi = np.arange(N)
            units = (xi % q != 0)
            H, ord2 = subgroup_pow2(q, n)
            phi = q ** (n - 1) * (q - 1)
            idx = phi // ord2
            inH = units & H
            outH = units & (~H)
            a2 = np.abs(mu) ** 2
            mass_tot = a2[units].sum()
            mass_inH = a2[inH].sum()
            mass_out = a2[outH].sum()
            fH = mass_inH / mass_tot
            fO = mass_out / mass_tot
            maxH = np.abs(mu[inH]).max() if inH.any() else 0.0
            maxO = np.abs(mu[outH]).max() if outH.any() else float('nan')
            ratio = (maxH / maxO) if (outH.any() and maxO > 0) else float('inf')
            # is the global argmax over units in H?
            amu = np.abs(mu).copy(); amu[~units] = -1
            argmax_xi = int(np.argmax(amu))
            argmax_inH = bool(H[argmax_xi])
            print(f"{q:>2} {n:>2} {N:>8} {ord2:>6} {phi:>7} {idx:>4} "
                  f"{fH:>13.6f} {fO:>13.6f} {ratio:>12.4g} {str(argmax_inH):>10}")
        print()

    print("Reading:")
    print("  q=3,5 -> idx=1, out-coset empty (fH=1.0 by construction): full mixing, no dark coset.")
    print("  q=7   -> idx=2: compare fH vs fO and maxH/maxO. fH>>fO + argmax_inH=True = half-mixing signature.")
    print("          fH~=fO would be the NULL (coset structure not visible in Plancherel mass).")

if __name__ == "__main__":
    main()
