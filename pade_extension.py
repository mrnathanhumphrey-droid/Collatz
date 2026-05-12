"""
pade_extension.py
=================
PADE_EXTENSION probe (2026-05-12, Wilson).

Re-examines R77.6's Padé approximants for FULL pole structure (not just
closest-to-z=2), runs the ratio diagnostic |eps_n|/|eps_{n-1}|, and
characterizes secondary-singularity structure.

INPUT:  experiments_output/result_77_7_eps_exact_through_k7.json (cached
        eps_n exact rationals, n=1..6).
OUTPUT: stdout tables matching the PADE_EXTENSION_*.md deliverables.

Exact rationals via fractions.Fraction throughout Padé construction.
numpy.roots used for pole-finding (floating point; documented precision).

Run from C:/Collatz/ as:  python pade_extension.py
"""
from __future__ import annotations

import json
import math
from fractions import Fraction

import numpy as np


# --------------------------------------------------------------------------- #
# Load cached eps_n (exact rationals)                                         #
# --------------------------------------------------------------------------- #

def load_eps(path="experiments_output/result_77_7_eps_exact_through_k7.json"):
    with open(path) as f:
        data = json.load(f)
    return {int(k): Fraction(int(v["num"]), int(v["den"])) for k, v in data.items()}


# --------------------------------------------------------------------------- #
# Padé construction over Q (copied from result_77_6_pade_construction.py)     #
# --------------------------------------------------------------------------- #

def gauss_solve_Q(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        pivot = -1
        for row in range(col, n):
            if M[row][col] != 0:
                pivot = row
                break
        if pivot == -1:
            return None
        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]
        piv = M[col][col]
        for j in range(col, n + 1):
            M[col][j] /= piv
        for row in range(n):
            if row != col and M[row][col] != 0:
                factor = M[row][col]
                for j in range(col, n + 1):
                    M[row][j] -= factor * M[col][j]
    return [M[i][n] for i in range(n)]


def pade_approximant(coeffs, m, n):
    if len(coeffs) < m + n + 1:
        raise ValueError(f"Need at least {m + n + 1} coefficients; got {len(coeffs)}")
    c = coeffs
    if n == 0:
        return [c[k] for k in range(m + 1)], [Fraction(1)]
    A, b = [], []
    for i in range(1, n + 1):
        row = []
        for j in range(1, n + 1):
            idx = m + i - j
            row.append(c[idx] if idx >= 0 else Fraction(0))
        A.append(row)
        b.append(-c[m + i])
    sol = gauss_solve_Q(A, b)
    if sol is None:
        return None
    Q = [Fraction(1)] + sol
    P = []
    for k in range(m + 1):
        pk = Fraction(0)
        for j in range(0, min(k, n) + 1):
            pk += Q[j] * c[k - j]
        P.append(pk)
    return P, Q


def find_roots(coeffs):
    cf = [float(c) for c in coeffs[::-1]]
    while cf and cf[0] == 0:
        cf.pop(0)
    if len(cf) <= 1:
        return []
    return list(np.roots(cf))


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #

def main():
    print("=" * 78)
    print("PADE_EXTENSION: extended Padé re-examination + ratio diagnostic")
    print("Wilson, 2026-05-12. R77.6 follow-up.")
    print("=" * 78)
    print()

    eps = load_eps()
    print(f"Loaded eps_n for n in {sorted(eps.keys())}")
    print()

    # ---- Phase 3: ratio diagnostic ---- #
    print("-" * 78)
    print("PHASE 3 — RATIO DIAGNOSTIC")
    print("-" * 78)
    print()
    print("Empirical |eps_n| values (exact):")
    print(f"{'n':>3} {'|eps_n| (float)':>18}")
    for n in sorted(eps.keys()):
        print(f"{n:3d} {abs(float(eps[n])):18.10e}")
    print()

    print("Raw ratios r_n := |eps_n|/|eps_{n-1}| (target: 0.5 for clean rate-1/2):")
    print(f"{'n':>3} {'r_n':>14} {'r_n - 0.5':>14} {'cumulative |r_n - 0.5|':>22}")
    ks = sorted(eps.keys())
    cum_abs = 0.0
    ratios = {}
    for i in range(1, len(ks)):
        n_prev, n_cur = ks[i - 1], ks[i]
        r = abs(eps[n_cur]) / abs(eps[n_prev])
        r_f = float(r)
        ratios[n_cur] = r_f
        dev = r_f - 0.5
        cum_abs += abs(dev)
        print(f"{n_cur:3d} {r_f:14.6f} {dev:+14.6f} {cum_abs:22.6f}")
    print()
    print("Differences for n=3..6 (excluding ratio_2 = 1/21 which is from transient eps_1):")
    print(f"{'n':>3} {'r_n':>14} {'r_n - 0.5':>14}")
    for n in [3, 4, 5, 6]:
        print(f"{n:3d} {ratios[n]:14.6f} {ratios[n] - 0.5:+14.6f}")
    print()
    print("Sign pattern of (r_n - 0.5): + - - - (peak above 0.5 at n=3, then below)")
    print()
    print("Geometric/oscillatory check: differences (r_n - 0.5):")
    diffs = [ratios[n] - 0.5 for n in [3, 4, 5, 6]]
    print(f"  n=3 -> n=4: ratio of consecutive diffs = {diffs[1]/diffs[0]:.4f}")
    print(f"  n=4 -> n=5: ratio of consecutive diffs = {diffs[2]/diffs[1]:.4f}")
    print(f"  n=5 -> n=6: ratio of consecutive diffs = {diffs[3]/diffs[2]:.4f}")
    print()
    print("Test: r_n = 0.5 + c/n fit (monotone correction):")
    print(f"{'n':>3} {'r_n':>14} {'(r_n - 0.5)*n':>16}")
    for n in [3, 4, 5, 6]:
        print(f"{n:3d} {ratios[n]:14.6f} {(ratios[n] - 0.5) * n:+16.6f}")
    print()
    print("Test: r_n = 0.5 + c/n^2 fit:")
    for n in [3, 4, 5, 6]:
        print(f"  n={n}: (r_n - 0.5)*n^2 = {(ratios[n] - 0.5) * n * n:+10.6f}")
    print()

    # ---- Phase 1: extended Padé table with FULL pole structure ---- #
    print("-" * 78)
    print("PHASE 1 — EXTENDED PADÉ POLE TABLE (full pole list per approximant)")
    print("-" * 78)
    print()
    print("Working series: f_tilde(z) = (E(z) - eps_1 z)/z^2")
    print("Coefficients: c_j = eps_{j+2} for j=0..4 (5 input coefficients)")
    print("Padé constraint: m + n <= 4 (no NEW approximants computable from N=5 data)")
    print()
    print("R77.6 already enumerated all (m,n) with m+n<=4. Re-examining for")
    print("FULL pole structure (secondary poles, complex poles, cluster patterns).")
    print()

    coeffs = [eps[2], eps[3], eps[4], eps[5], eps[6]]
    targets = [
        (1, 1, "diagonal"),
        (2, 1, "sub-diag m+n=3"),
        (1, 2, "super-diag m+n=3"),
        (3, 1, "sub-diag m+n=4"),
        (2, 2, "diagonal m+n=4"),
        (1, 3, "super-diag m+n=4"),
        (4, 0, "Taylor"),
        (0, 4, "all-pole"),
    ]

    all_poles = []  # (m, n, role, idx, complex pole)
    for m, n, role in targets:
        if m + n + 1 > len(coeffs):
            continue
        out = pade_approximant(coeffs, m, n)
        if out is None:
            print(f"  [{m}/{n}]: SINGULAR")
            continue
        P, Q = out
        roots = find_roots(Q)
        print(f"  [{m}/{n}] ({role}):")
        if not roots:
            print(f"    no poles (n=0 approximant)")
        else:
            # sort by distance to z=2
            roots_sorted = sorted(roots, key=lambda r: abs(r - 2.0))
            for j, r in enumerate(roots_sorted):
                d2 = abs(r - 2.0)
                tag = "PRIMARY" if j == 0 else f"sec_{j}"
                marker = " (complex)" if abs(r.imag) > 1e-9 else ""
                print(f"    {tag}: z = {r.real:+10.6f} {r.imag:+10.6f}j   |z| = {abs(r):8.4f}"
                      f"   |z-2| = {d2:8.4f}{marker}")
                all_poles.append((m, n, role, j, complex(r)))
        print()

    # ---- Phase 2: pattern across approximants ---- #
    print("-" * 78)
    print("PHASE 2 — PATTERN ACROSS APPROXIMANTS")
    print("-" * 78)
    print()
    print("(a) Closest-pole-to-z=2 across all approximants:")
    print(f"{'(m,n)':>8} {'role':>20} {'pole z':>30} {'|z-2|':>10}")
    by_mn = {}
    for m, n, role, idx, z in all_poles:
        by_mn.setdefault((m, n, role), []).append((idx, z))
    for (m, n, role), pairs in sorted(by_mn.items()):
        pairs.sort(key=lambda x: abs(x[1] - 2.0))
        idx, z = pairs[0]
        print(f"  ({m},{n}) {role:>20} {z.real:+10.4f}{z.imag:+10.4f}j {abs(z - 2.0):10.4f}")
    print()

    print("(b) All poles with |z - 2| < 1 (potential primary cluster):")
    real_primary = []
    imag_primary = []
    for m, n, role, idx, z in all_poles:
        if abs(z - 2.0) < 1.0:
            print(f"  ({m},{n}) idx={idx}: z = {z.real:+10.6f} {z.imag:+10.6f}j  |z-2|={abs(z - 2.0):.4f}")
            real_primary.append(z.real)
            imag_primary.append(z.imag)
    print()
    print(f"  Primary cluster: real parts in [{min(real_primary):.4f}, {max(real_primary):.4f}]")
    print(f"                   imag parts in [{min(imag_primary):.4e}, {max(imag_primary):.4e}]")
    print(f"  Range of |z-2|: 0.05 .. 0.35  (NOT within 0.05 radius — pre-reg consistency check FAILS)")
    print()

    print("(c) Diagonal [n/n] subsequence (Stahl-natural probe):")
    for m, n, role, idx, z in all_poles:
        if m == n:
            d2 = abs(z - 2.0)
            print(f"  [{m}/{n}] pole idx={idx}: z = {z.real:+10.6f} {z.imag:+10.6f}j  |z-2| = {d2:.4f}")
    print()
    print("  [1/1] primary pole at 2.0764, [2/2] primary at 2.0513 — monotone descent.")
    print("  Ratio 0.0513/0.0764 = 0.67. Per R77.6: inconsistent with simple pole.")
    print()

    print("(d) Off-axis poles (|imag| > 0.01):")
    off_axis = [(m, n, role, idx, z) for (m, n, role, idx, z) in all_poles if abs(z.imag) > 0.01]
    if off_axis:
        for m, n, role, idx, z in off_axis:
            print(f"  ({m},{n}) idx={idx} {role}: z = {z.real:+10.4f} {z.imag:+10.4f}j"
                  f"  |z| = {abs(z):.4f}  arg = {math.degrees(math.atan2(z.imag, z.real)):+.1f}deg")
    else:
        print("  (none)")
    print()

    print("(e) Secondary poles per approximant:")
    for (m, n, role), pairs in sorted(by_mn.items()):
        if len(pairs) <= 1:
            continue
        pairs.sort(key=lambda x: abs(x[1] - 2.0))
        sec_poles = pairs[1:]
        sec_strs = [f"z = {z.real:+.4f}{z.imag:+.4f}j (|z-2|={abs(z - 2.0):.4f})"
                    for _, z in sec_poles]
        print(f"  ({m},{n}) secondary poles: {sec_strs}")
    print()
    print("  Notes:")
    print("  - [2/2] secondary at z=+0.6878 (REAL, inside |z|<1). Likely artifact:")
    print("    not present in [1/1] (no room) nor [3/1] (one pole only); does NOT")
    print("    correspond to expected next singularity at z=4 (R76 §10).")
    print("  - [1/2] spurious at z=+155 — classical Padé artifact (pole at large |z|).")
    print("  - [1/3] spurious at z=-12, z=+8 — Padé instability with 1 num + 3 denom.")
    print("  - [0/4] all-pole has TWO cc-pairs:")
    print("       z = -0.95 ± 0.76j on |z|≈1.22")
    print("       z = +0.45 ± 1.01j on |z|≈1.10  <-- INSIDE expected disk |z|=2")
    print("    These are [0/4] artifacts (numerator forced to constant); off-axis")
    print("    poles on circles of |z| approaching ρ are the classical 'Froissart")
    print("    doublet'-like behavior for [0/n] approximants of functions with")
    print("    branch structure, NOT a real secondary singularity.")
    print()

    print("(f) Pre-registered consistency check:")
    print("  Definition: closest poles across approximants within radius 0.05 of a")
    print("  common point => 'consistent'.")
    primary_zs = []
    for (m, n, role), pairs in sorted(by_mn.items()):
        pairs.sort(key=lambda x: abs(x[1] - 2.0))
        primary_zs.append(pairs[0][1])
    # All primary closest-to-z=2 are real (off-axis only in [0/4]); compute spread
    pure_primary = [z for z in primary_zs if abs(z.imag) < 1e-9 and z.real < 5]
    if pure_primary:
        reals = [z.real for z in pure_primary]
        spread = max(reals) - min(reals)
        center = sum(reals) / len(reals)
        max_dev = max(abs(r - center) for r in reals)
        print(f"  Primary real-axis poles: {[f'{r:.4f}' for r in reals]}")
        print(f"  Spread: {spread:.4f}, center: {center:.4f}, max deviation: {max_dev:.4f}")
        print(f"  Consistency (within radius 0.05 of common point): "
              f"{'PASS' if max_dev < 0.05 else 'FAIL'}")
        print()
        # Restrict to diagonals only
        diag_reals = []
        for (m, n, role), pairs in by_mn.items():
            if m == n:
                pairs.sort(key=lambda x: abs(x[1] - 2.0))
                if abs(pairs[0][1].imag) < 1e-9:
                    diag_reals.append(pairs[0][1].real)
        if diag_reals:
            print(f"  Diagonal-only [{1},{1}], [{2},{2}] primary poles: {[f'{r:.4f}' for r in diag_reals]}")
            print(f"  Diagonal spread: {max(diag_reals) - min(diag_reals):.4f}")
    print()

    # ---- Phase 4 synthesis ---- #
    print("-" * 78)
    print("PHASE 4 — SYNTHESIS")
    print("-" * 78)
    print()
    print("Padé pattern (Phase 1+2):")
    print("  - Closest poles span [2.05, 2.35] across approximants (NOT within 0.05).")
    print("  - Diagonal subsequence is tight & monotone (2.076 -> 2.051), ratio 0.67.")
    print("  - Off-diagonals push closest pole AWAY from z=2 (2.13 .. 2.35).")
    print("  - Only off-axis poles appear in [0/4] (cc-pairs); other approximants")
    print("    have all closest+real-secondary poles on the positive real axis.")
    print("  - NO complex-conjugate primary poles. Secondary structure shows NO")
    print("    evidence of complex secondary singularities on the principal sheet.")
    print()
    print("Ratio pattern (Phase 3):")
    print("  - r_3 = +0.5347 (ABOVE 0.5 by 0.035)")
    print("  - r_4 = +0.4816 (below by 0.018)")
    print("  - r_5 = +0.4697 (below by 0.030)")
    print("  - r_6 = +0.4323 (below by 0.068)")
    print()
    print("  Sign pattern: + - - -.  NOT monotone.")
    print("  NOT geometric (ratio of successive differences: -0.53, +1.63, +2.24).")
    print("  NOT 1/n fit ((r_n - 0.5)*n drifts down from +0.10 to -0.41).")
    print("  NOT 1/n^2 fit either ((r_n-0.5)*n^2: +0.31, -0.30, -0.76, -2.44).")
    print()
    print("  KEY FINDING: ratios are LEAVING 0.5 going downward through n=4,5,6,")
    print("  with departure ACCELERATING (0.018 -> 0.030 -> 0.068 absolute).")
    print()
    print("  This is a strong signal: the SECONDARY structure is causing |eps_n|")
    print("  to fall faster than (1/2)^n at n=4..6. Two consistent readings:")
    print("    (i) Secondary pole at |z|=ρ<2 with NEGATIVE coefficient subtracting")
    print("        from |eps_n|.")
    print("    (ii) Branch-cut at z=2 with subleading correction of form")
    print("         (1/2)^n * (negative-coefficient * n^{-α} or log(n) term).")
    print()
    print("Disposition (see PADE_EXTENSION_DISPOSITION.md for full reasoning):")
    print("  -- Phase 1 Padé: no NEW approximants computable from N=5; existing")
    print("     ones show diagonal monotone descent + off-diagonal scatter,")
    print("     ALL primary poles real-axis above z=2, no complex secondary.")
    print("  -- Phase 3 Ratio: deviations from 1/2 are non-monotone (sign + then -)")
    print("     and ACCELERATING departure downward for n=4..6.")
    print("  -- Both diagnostics agree: branch-cut/mixed structure at z=2,")
    print("     subleading correction with negative drift, no complex secondary.")
    print("  -- DISPOSITION: H_AMBIGUOUS / INCONCLUSIVE (consistent with R77.6")
    print("     baseline). The probe RULES OUT H_COMPLEX_SECONDARY (no oscillation")
    print("     in ratios, no complex primary Padé poles) but does NOT close the")
    print("     H_BRANCH_CUT_ORDER question. Needs eps_7, eps_8 (Route A).")
    print()


if __name__ == "__main__":
    main()
