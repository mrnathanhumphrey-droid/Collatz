"""
PROBE FOURCELL -- Wilson's direct four-cell decomposition + rate question (2026-07-25).

Direct (no w=-u-v substitution, no cell-mixing): q_r(1)-1/3 = U0 + U3 + V1 + W2, where each cell uses its OWN
rotation c mod 3 (c=0,3 -> rot 0 = u; c=1 -> rot 1 = v; c=2 -> rot 2 = w):
    U_c = p_c * E[own-rotation coherence | c],   p_c = mu-measure of cell c,   mu(x) prop nu_low(x) nu_low(4x mod 3^r).
    u(x)=<dpi_x,dpi_L>, v=<dpi_x,T1 dpi_L>, w=<dpi_x,T2 dpi_L>;  U0=p0 E[u|0], U3=p3 E[u|3], V1=p1 E[v|1], W2=p2 E[w|2].

THE QUESTION (Wilson): which of U0,U3,V1,W2 decays at the EXCESS rate 0.89 (vs faster)? Two-step rates only
(the parity alias -- 3^r mod 4 = 1 (r even) / 3 (r odd) -- fooled the single-step rate). Three outcomes:
  all four ~0.89 same-signed => bulk, per-cell positivity, tractable;
  individually alternating, only sum steady => carry conspiracy, harder;
  mixed (one slow) => sign lives in that cell.
Plus: parity split (even-r vs odd-r subsequence per cell) -- does the alternation subtract out to smooth bulk?
"""
import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_gapop_R28 import build_nu

RTOP = 16
SCRATCH = r"C:\Users\Nate\AppData\Local\Temp\claude\c--As-Above-So-Below-Master\7870203b-c92a-4380-9d18-22b6077a7418\scratchpad"


def four_cells(nu_r, r):
    Mp = 3 ** (r + 1); M = 3 ** r
    hi = np.zeros(Mp)
    for X, w in nu_r.items():
        hi[X % Mp] += float(w)
    hi /= hi.sum()
    low = hi[:M] + hi[M:2 * M] + hi[2 * M:3 * M]
    xs = np.nonzero(low > 0)[0]
    c = (4 * xs) // M
    L = (4 * xs) % M
    wgt = low[xs] * low[L]
    keep = wgt > 0
    xs, c, L, wgt = xs[keep], c[keep], L[keep], wgt[keep]
    Z = wgt.sum()
    dpx = np.stack([hi[xs + d * M] / low[xs] for d in range(3)], 1) - 1.0 / 3
    dpL = np.stack([hi[L + d * M] / low[L] for d in range(3)], 1) - 1.0 / 3
    coh = {0: (dpx * dpL).sum(1),
           1: (dpx * np.roll(dpL, -1, axis=1)).sum(1),
           2: (dpx * np.roll(dpL, -2, axis=1)).sum(1)}
    out = {}
    for cell in range(4):
        sel = (c == cell)
        pc = float(wgt[sel].sum() / Z)
        rot = cell % 3
        Ucoh = float((wgt[sel] * coh[rot][sel]).sum() / wgt[sel].sum()) if sel.any() else 0.0
        out[cell] = (pc, Ucoh, pc * Ucoh)      # (p_c, E[own-coh|c], U_c)
    ex = out[0][2] + out[3][2] + out[1][2] + out[2][2]
    return out, ex, Z


def two_step(seq, r):
    """per-step rate from (|x_r|/|x_{r-2}|)^(1/2), sign-tracking."""
    if r - 2 not in seq or seq[r][2] == 0 or seq[r - 2][2] == 0:
        return float('nan')
    return (abs(seq[r][2]) / abs(seq[r - 2][2])) ** 0.5


def main():
    t0 = time.time()
    print("# PROBE FOURCELL -- direct four-cell decomposition + rate question\n")
    print(f"building build_nu to {RTOP} ... (~7 min)")
    nus = build_nu(0.5, RTOP)
    print(f"  built ({time.time()-t0:.1f}s)\n")

    data = {}
    for r in range(4, RTOP + 1):
        data[r], ex, _ = four_cells(nus[r], r)

    # ---- the four terms per r ----
    print("## the four cell terms  U0 U3 V1 W2  (sum = q-1/3);  q-1/3 for reference")
    print(f"   {'r':>2} {'U0':>11} {'U3':>11} {'V1':>11} {'W2':>11} | {'q-1/3':>11} {'U0+U3':>11} {'V1+W2':>11}")
    for r in range(4, RTOP + 1):
        U0, U3, V1, W2 = data[r][0][2], data[r][3][2], data[r][1][2], data[r][2][2]
        ex = U0 + U3 + V1 + W2
        print(f"   {r:>2} {U0:>+11.3e} {U3:>+11.3e} {V1:>+11.3e} {W2:>+11.3e} | {ex:>+11.3e} {U0+U3:>+11.3e} {V1+W2:>+11.3e}")
    print()

    # ---- signs across r (does each term alternate, or hold sign?) ----
    print("## sign pattern across r=4..16 (does each cell term ALTERNATE or HOLD sign?)")
    for name, cell in (('U0', 0), ('U3', 3), ('V1', 1), ('W2', 2)):
        signs = ''.join('+' if data[r][cell][2] > 0 else '-' for r in range(4, RTOP + 1))
        print(f"   {name}: {signs}")
    print("   excess: " + ''.join('+' if (data[r][0][2] + data[r][3][2] + data[r][1][2] + data[r][2][2]) > 0
                                   else '-' for r in range(4, RTOP + 1)))
    print()

    # ---- two-step rates (the decisive measurement) ----
    print("## TWO-STEP per-step rates (|x_r/x_{r-2}|)^.5 vs EXCESS 0.89  [which term decays at the slow rate?]")
    seqs = {'U0': {r: data[r][0] for r in data}, 'U3': {r: data[r][3] for r in data},
            'V1': {r: data[r][1] for r in data}, 'W2': {r: data[r][2] for r in data}}
    exseq = {r: (0, 0, data[r][0][2] + data[r][3][2] + data[r][1][2] + data[r][2][2]) for r in data}
    print(f"   {'r':>2} {'U0':>7} {'U3':>7} {'V1':>7} {'W2':>7} {'excess':>7}")
    for r in range(12, RTOP + 1):
        row = [two_step(seqs[n], r) for n in ('U0', 'U3', 'V1', 'W2')]
        er = two_step(exseq, r)
        print(f"   {r:>2} " + " ".join(f"{x:>7.3f}" for x in row) + f" {er:>7.3f}")
    print("   [~0.89 = decays with the excess (asymptotic carrier); << 0.89 (e.g. ~0.76) = transient, vanishes vs excess.]")
    print()

    # ---- parity / bulk-boundary split ----
    print("## PARITY SPLIT (bulk vs boundary): even-r and odd-r subsequences per term (two-step within parity)")
    for name, cell in (('U0', 0), ('U3', 3), ('V1', 1), ('W2', 2)):
        ev = [data[r][cell][2] for r in range(4, RTOP + 1) if r % 2 == 0]
        od = [data[r][cell][2] for r in range(5, RTOP + 1) if r % 2 == 1]
        rev = (abs(ev[-1]) / abs(ev[-2])) ** 0.5 if len(ev) >= 2 and ev[-2] != 0 else float('nan')
        rod = (abs(od[-1]) / abs(od[-2])) ** 0.5 if len(od) >= 2 and od[-2] != 0 else float('nan')
        se = '+' if ev[-1] > 0 else '-'; so = '+' if od[-1] > 0 else '-'
        print(f"   {name}: even-r rate {rev:.3f} (sign {se}) | odd-r rate {rod:.3f} (sign {so})"
              f"  {'PARITY-SMOOTH (bulk, boundary subtracts)' if abs(rev-rod) < 0.05 and se == so else 'parity structure remains'}")
    print()

    # save
    tbl = {r: {n: data[r][c] for n, c in (('U0', 0), ('U3', 3), ('V1', 1), ('W2', 2))} for r in data}
    json.dump({str(r): {k: list(v) for k, v in tbl[r].items()} for r in tbl},
              open(os.path.join(SCRATCH, 'fourcell.json'), 'w'))
    print(f"  [four-cell table saved -> scratchpad/fourcell.json]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
