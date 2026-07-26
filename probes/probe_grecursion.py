"""
PROBE GRECURSION -- Wilson's high-v_2 recursion in closed form (2026-07-26).

On the lattice xi = 3^i 2^j the genuine-halving branch (a<=j, no wraparound) has UNIT-MODULUS coefficients:
  F(i,j) := pi-hat_k(3^i 2^j),  G(i,j) := 2^j F(i,j)
  G(i,j) = Sum_{b<j} e(2^b / 3^{k-i}) G(i+1,b)   [+ wrap terms a>j, bounded by 2^-j S vs main 2^-j G]
  |pi-hat(3^i 2^j)| = 2^-j |G(i,j)|.

Three checks (Wilson's run order -- #3 FIRST):
 3. PIN m(k): argmax|pi-hat| = 2^m (or mirror N-2^m). Tabulate m(k), dm/dk, and TEST sup rate vs 2^-dm.
    If dm/dk accounts for the observed 0.80/step, the decay is a COUNTING fact (where the argmax sits),
    not a cancellation fact. Also |G(0,m)| = 2^m S -- is it O(1) (Wilson) and does it converge?
 1. GATE the (truncated, genuine-halving) G-recursion vs pi-hat on the lattice -> machine precision?
    If not, the wrap terms are NOT negligible.
 2. |G(0,j)| profile across j -- bounded, or growing?

Reuses probe_singlerec.fwd_hat (sign-fixed pi-hat). Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_singlerec import fwd_hat


def v2(n):
    n = int(n); v = 0
    while n and n % 2 == 0:
        n //= 2; v += 1
    return v


def main():
    t0 = time.time()
    print("# PROBE GRECURSION -- high-v_2 recursion in closed form (G(i,j)=2^j pi-hat(3^i 2^j))\n")

    KS = list(range(3, 13))
    phs = {k: fwd_hat(k)[0] for k in KS}

    # ---------- #3 PIN m(k) ----------
    print("## #3  PIN m(k): argmax|pi-hat| sits at 2^m (or mirror N-2^m). m = exponent of the power of 2.")
    print(f"   {'k':>2} {'xi*':>8} {'m':>3} {'dm':>4} {'S=|pi(xi*)|':>12} {'Srate':>7} {'2^-dm':>7} "
          f"{'|G(0,m)|=2^m S':>14} {'Grate':>7}")
    rows = []
    prev_m = None; prev_S = None; prev_G = None
    for k in KS:
        N = 3 ** k
        ph = phs[k]
        prim = np.array([xi for xi in range(1, N) if xi % 3 != 0], dtype=np.int64)
        absp = np.abs(ph[prim]); j = int(np.argmax(absp))
        xistar = int(prim[j]); S = float(absp[j])
        # identify m: the power of 2 in the degenerate pair {2^m, N-2^m}
        even = xistar if xistar % 2 == 0 else N - xistar
        m = v2(even)
        assert even == 2 ** m, f"k={k}: argmax {xistar} not a power-of-2 pair (even twin {even} != 2^{m})"
        G = (2 ** m) * S
        dm = '' if prev_m is None else f"{m-prev_m:+d}"
        srate = '' if prev_S is None else f"{S/prev_S:.4f}"
        twodm = '' if prev_m is None else f"{2.0**(-(m-prev_m)):.4f}"
        grate = '' if prev_G is None else f"{G/prev_G:.4f}"
        print(f"   {k:>2} {xistar:>8} {m:>3} {dm:>4} {S:>12.6f} {srate:>7} {twodm:>7} {G:>14.4f} {grate:>7}")
        rows.append((k, m, S, G)); prev_m, prev_S, prev_G = m, S, G
    print("   [Wilson: if Srate ~ 2^-dm at every step, decay is COUNTING (m grows), not cancellation.")
    print("    if |G(0,m)| stays O(1) and converges, that constant is the coeff in sup|pi| ~ c 2^-m.]")

    # ---------- #1 GATE the truncated G-recursion ----------
    print("\n## #1  GATE truncated G-recursion  G(i,j) =?= Sum_{b<j} e(2^b/3^{k-i}) G(1?,b)  at (i=0, j=m)")
    print(f"   {'k':>2} {'m':>3} {'|G(0,m)| direct':>15} {'|RHS trunc|':>12} {'rel(trunc)':>11} {'rel(+wrap)':>11}")
    for (k, m, S, Gval) in rows:
        N = 3 ** k
        ph = phs[k]
        i = 0
        # direct: G(0,m) = 2^m pi-hat(2^m mod N)
        Gdirect = (2 ** m) * ph[(2 ** m) % N]
        # truncated RHS: Sum_{b=0}^{m-1} e(2^b/3^{k-i}) * 2^b pi-hat(3^{i+1} 2^b mod N)
        rhs_trunc = 0j
        for b in range(0, m):
            ph_e = np.exp(2j * np.pi * (2 ** b) / (3 ** (k - i)))
            Gib = (2 ** b) * ph[(3 ** (i + 1) * 2 ** b) % N]
            rhs_trunc += ph_e * Gib
        # full RHS: include wrap terms a>m  (b=j-a<0 -> 2^-a with a>m; xi 2^-a wraps mod N)
        inv2 = pow(2, -1, N)
        rhs_full = rhs_trunc
        # a from m+1 upward: contribution 2^{m-a} e(xi* 2^-a/N) pi-hat(3 xi* 2^-a mod N), xi*=2^m
        p = pow(inv2, m, N)   # 2^-m
        xistar_e = (2 ** m) % N
        for a in range(m + 1, m + 60):
            p = (p * inv2) % N          # 2^-a
            w = (xistar_e * p) % N      # xi* 2^-a
            rhs_full += (2.0 ** (m - a)) * np.exp(2j * np.pi * w / N) * ph[(3 * w) % N]
        rel_t = abs(rhs_trunc - Gdirect) / (abs(Gdirect) + 1e-30)
        rel_f = abs(rhs_full - Gdirect) / (abs(Gdirect) + 1e-30)
        print(f"   {k:>2} {m:>3} {abs(Gdirect):>15.6f} {abs(rhs_trunc):>12.6f} {rel_t:>11.2e} {rel_f:>11.2e}")
    print("   [rel(+wrap) ~ machine 0 => recursion exact; rel(trunc) ~ 2^-k => wrap negligible as Wilson claims.]")

    # ---------- #2 |G(0,j)| profile ----------
    print("\n## #2  |G(0,j)| = 2^j |pi-hat(2^j)| profile across j  (bounded, or growing?)")
    for k in (7, 10, 12):
        N = 3 ** k; ph = phs[k]
        vals = []
        for j in range(1, 2 * k + 6):
            xi = (2 ** j) % N
            vals.append((j, (2.0 ** j) * abs(ph[xi])))
        jstar = max(vals, key=lambda t: t[1])
        prof = " ".join(f"{g:.2f}" for _, g in vals[:2 * k + 2])
        print(f"   k={k:>2}: j*={jstar[0]} |G|max={jstar[1]:.3f}  profile |G(0,j)| j=1..: {prof}")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
