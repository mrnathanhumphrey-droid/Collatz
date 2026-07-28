"""
PROBE CF-CANCEL -- integer-level corroboration of the log3/log2 rotation number (2026-07-27)

Wilson thought 1, item 2 (the ONE honest new test; low priority). The base-2 oscillation sampled at
integer base-3 levels is an irrational rotation with number  alpha = log3/log2 - 1 = 0.58496...
The DMT continued-fraction convergents of log3/log2 -- numerators 1,2,3,8,19,65 -- are the levels
where the phase best realigns (best cancellation). This is a statement about INTEGER levels (which
exist), unlike the half-integer strobe (which does not -- the tower is 3-to-1 with no intermediate
group; M^{1/2} is a different operator whose (-1)^n would be a branch-cut artifact). So we test the
convergents on the banked Lambda_i, i<=20.

HONEST caveat up front: 13 points, only i in {8,19} (h_k numerators) and {12} (q_k denominators)
land in range. This CORROBORATES an already-confirmed rotation number (it is the 2pi/log2=9.06 period
restated); it decides nothing. A clean alignment is weak positive; a miss is near-meaningless.
"""
import numpy as np
from fractions import Fraction

L3 = np.log(3); L2 = np.log(2)

# banked T_i (0..20): T_0=1/3, T_1..14 recomputed elsewhere = exact-ish; use the certified banked tail.
T = {0: 1.0 / 3,
     15: 0.23567582169638104, 16: 0.23591007771310188, 17: 0.2361167285902292,
     18: 0.23629629695785445, 19: 0.23645299645774767, 20: 0.2365885345278137}
# T_1..14 from the same machinery (probe_lattice recomputes them); import to avoid drift.
from probe_lattice import T_level
for n in range(1, 15):
    T[n] = T_level(n)
I = np.arange(1, 21)
Lam = np.array([T[i] - T[i - 1] for i in I])


def convergents(x, n):
    """first n continued-fraction convergents (h_k, k_k) of x."""
    a = []; y = x
    for _ in range(n):
        ai = int(np.floor(y)); a.append(ai)
        fr = y - ai
        if fr < 1e-12: break
        y = 1.0 / fr
    h0, h1, k0, k1 = 1, a[0], 0, 1
    out = [(a[0], 1)]
    for ai in a[1:]:
        h0, h1 = h1, ai * h1 + h0
        k0, k1 = k1, ai * k1 + k0
        out.append((h1, k1))
    return a, out


def main():
    print("# PROBE CF-CANCEL -- integer-level rotation-number corroboration\n")
    alpha = L3 / L2 - 1.0
    print(f"rotation number alpha = log3/log2 - 1 = {alpha:.6f}   (=0.585, mescaline; = 2pi/log2 period restated)")
    a, cv = convergents(L3 / L2, 8)
    print(f"CF(log3/log2) = {a}")
    print(f"convergents h_k/k_k = {cv}")
    hnum = [h for h, k in cv]; kden = [k for h, k in cv]
    print(f"  DMT numerators h_k = {hnum}  (predicted best-cancellation levels)")
    print(f"  denominators  k_k = {kden}  (irrational-rotation return times)\n")

    # detrend Lambda by its own log-linear rho1 (i=8..20 clean regime), get oscillatory residual
    idx = I >= 8
    b, aa = np.polyfit(I[idx], np.log(Lam[idx]), 1)
    res = np.log(Lam) - (aa + b * I)                     # residual defined for all i, fit on i>=8
    print(f"detrended log-residual res_i (rho1 = e^{b:.4f} = {np.exp(b):.4f}/level):")
    print("  i :  res_i")
    for i, r in zip(I, res):
        mark = ""
        if i in hnum: mark += "  <- h_k (DMT best-cancel)"
        if i in kden: mark += "  <- k_k (return time)"
        print(f"  {i:2d}: {r:+.4f}{mark}")

    # phase model: predicted phase phi_i = 2pi * frac(i * (log3/log2)); best cancel when phi near 0 mod 2pi
    print("\nphase check: phi_i = 2*pi*frac(i*log3/log2); |sin(phi_i/2)| small => near a cancellation node")
    print("  i :  frac(i*a)   |dist to integer|   (small at return times)")
    for i in I:
        f = (i * (L3 / L2)) % 1.0
        d = min(f, 1 - f)
        mark = "  <-- h_k" if i in hnum else ("  <-- k_k" if i in kden else "")
        print(f"  {i:2d}:  {f:.4f}      {d:.4f}{mark}")

    # is |res| locally minimal at the in-range convergent levels?
    inrange_h = [i for i in hnum if 1 <= i <= 20]
    inrange_k = [i for i in kden if 1 <= i <= 20]
    absres = np.abs(res)
    print(f"\nin-range h_k levels: {inrange_h};  k_k levels: {inrange_k}")
    print(f"mean |res| all i>=8 = {absres[I >= 8].mean():.4f}")
    for i in sorted(set(inrange_h + inrange_k)):
        if i >= 8:
            print(f"  |res_{i}| = {absres[i-1]:.4f}  (vs mean {absres[I>=8].mean():.4f}) "
                  f"{'below-mean (aligns)' if absres[i-1] < absres[I>=8].mean() else 'above-mean (no align)'}")
    print("\nVERDICT NOTE: corroboration only. Rotation number already confirmed by 2pi/log2=9.06.")


if __name__ == "__main__":
    main()
