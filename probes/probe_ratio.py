"""
PROBE RATIO -- the successive-ratio sequence of the dominant fluctuation mode Re dhat_r(1).

Wilson's spec (2026-07-25). The MODES probe left four ratio rows; this computes all eleven
levels (r=2..12) and reads the successive ratio rho_r := Re dhat_{r+1}(1)/Re dhat_r(1), r=2..11.

The conditional (Wilson): if Re dhat_r(1)=sum_i c_i rho_i^r with finitely many REAL rho_i>0,
and the (positive) sequence has successive ratios RISING to a limit, the largest-rho component
carries c>0 -- a negative dominant coefficient would force the sequence through zero, at which
point successive ratios FALL to zero and go negative rather than rising. So:
  MONOTONE RISE + positive sequence  ==>  no approaching sign change of the dominant mode.
  OSCILLATION in rho_r               ==>  a COMPLEX PAIR ==> conditional void, 7/15 branch live.

Three reads:
 1. Monotone rise, or oscillation? (sign pattern of the successive differences of rho_r)
 2. Does rho_r turn DOWN anywhere? (the standing monitor: while rho rises, Re dhat(1) is not
    about to flip; a downturn is the leading indicator of a dominant-mode crossing)
 3. What LIMIT does rho_r curve toward -- 0.91 or the banked ~0.984? (Aitken dt^2; tail sum ~1/(1-rho))
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import probe_charledger_R10 as R10
from probe_gapop_R28 import build_nu

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
    return dd


def aitken(seq):
    """Aitken delta^2 extrapolation: x_hat = x_n - (dx)^2/(ddx). Returns array of accelerated estimates."""
    s = np.asarray(seq, float)
    out = []
    for i in range(len(s) - 2):
        d1 = s[i + 1] - s[i]; d2 = s[i + 2] - 2 * s[i + 1] + s[i]
        out.append(s[i] - d1 * d1 / d2 if abs(d2) > 1e-18 else float('nan'))
    return np.array(out)


def main():
    t0 = time.time()
    print(f"# PROBE RATIO -- successive ratio of Re dhat_r(1), r=2..{JMAX}\n")
    nus = build_nu(0.5, JMAX)
    d1 = {}
    for j in range(2, JMAX + 1):
        H = np.fft.fft(delta_from_nu(nus[j], j))
        d1[j] = H[1].real
    print(f"  Re dhat_r(1) built ({time.time()-t0:.1f}s)\n")

    print("## Re dhat_r(1) (the dominant fluctuation mode)")
    for r in range(2, JMAX + 1):
        print(f"   r={r:>2}  {d1[r]:>+14.7e}  {'(+)' if d1[r] > 0 else '(-)  <-- SIGN FLIP'}")
    print()

    rho = {r: d1[r + 1] / d1[r] for r in range(2, JMAX)}
    rs = list(range(2, JMAX))
    print("## READ 1&2 -- rho_r = Re dhat_{r+1}(1)/Re dhat_r(1), and its successive change")
    print(f"   {'r':>2} {'rho_r':>10} {'d(rho)':>11}  monotone?")
    prev = None; downturns = []; signs = []
    for r in rs:
        drho = None if prev is None else rho[r] - prev
        if drho is not None:
            signs.append(1 if drho > 0 else -1)
            if drho < 0:
                downturns.append(r)
        tag = "" if drho is None else ("RISE" if drho > 0 else "**DOWNTURN**")
        ds = "" if drho is None else f"{drho:>+11.4f}"
        print(f"   {r:>2} {rho[r]:>10.4f} {ds:>11}  {tag}")
        prev = rho[r]
    # oscillation = any sign change in the successive-difference sign pattern
    flips = sum(1 for i in range(1, len(signs)) if signs[i] != signs[i - 1])
    print()
    print(f"   sign pattern of d(rho): {' '.join('+' if s > 0 else '-' for s in signs)}")
    print(f"   downturns at r = {downturns if downturns else 'NONE'}")
    print(f"   direction flips in d(rho): {flips}   "
          f"=> {'MONOTONE (real dominant mode consistent)' if flips == 0 else 'OSCILLATION (complex pair -> conditional VOID, 7/15 live)'}")
    print()

    # ---- READ 3: limit ----
    print("## READ 3 -- what limit does rho_r curve toward?")
    rvals = [rho[r] for r in rs]
    ait = aitken(rvals)
    print("   raw rho tail:   " + " ".join(f"{v:.4f}" for v in rvals))
    print("   Aitken dt^2:    " + " ".join(f"{v:.4f}" for v in ait))
    # tail rate on the clean monotone tail (last 4 points)
    tail = rvals[-4:]
    print(f"   last raw rho = {rvals[-1]:.4f};  last Aitken = {ait[-1]:.4f} (if finite)")
    for cand in (0.907, 0.984):
        print(f"   tail-sum factor 1/(1-rho) at rho={cand}: {1.0/(1.0-cand):.2f}")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
