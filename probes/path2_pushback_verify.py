"""
path2_pushback_verify.py — Adversarial extension of path2_family_verify.py.

Implements Check 3 (Hensel tightness at r=4,5,6) and Check 5 (extended primes
p ∈ {13,17,19,23,29,31} at r ∈ {2,3}). All checks operate on the same
F̂_p / G_p machinery as the original Phase 2 script.

What is computed per cell (p,r):
  - C1 bijection a ↔ C_a (exhaustive count)
  - C2 |G_p(a)| = p^{(r+1)/2} (max deviation)
  - C3 r=2 uniform Gaussian factor (only when r=2)
  - C4 r=3 saddle exactness (only when r=3); recorded as max_saddle_diff
  - C5 |S_partial| empirical = |Σ_{a∈supp} 1̂(p·a) · G_p(a)| (DIRECT bilinear sum)
  - C6 ratio |S_partial| / √N (the empirical "constant")
  - bound2sqrtN = 2·√N (claimed at r ≤ 3)
  - boundHenselLog = 2·√N · (1 + log N)  (claimed at r ≥ 4)

Output: PATH2_PUSHBACK_EXTENDED.csv
"""

import sys
import os
import math
import cmath
from fractions import Fraction

# Force UTF-8 stdout so the √ and ↔ characters in print statements don't crash
# the default Windows cp1252 codepage.
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

OUTDIR = r"C:\Collatz"
CSV_PATH = os.path.join(OUTDIR, "PATH2_PUSHBACK_EXTENDED.csv")


def J_for_p(p, m):
    """Max j s.t. for all k≤j, k − v_p(k) < m. Same as original script."""
    j = 1
    while True:
        x = j + 1
        v = 0
        while x % p == 0:
            x //= p
            v += 1
        if (j + 1) - v >= m:
            return j
        j += 1


def truncated_padic_log(p, s, J):
    L = Fraction(0)
    for j in range(1, J + 1):
        L += Fraction((-1) ** (j - 1), j) * Fraction(p * s) ** j
    return L


def L_mod_q(p, s, J, q):
    if s == 0:
        return 0
    L_frac = truncated_padic_log(p, s, J)
    num = L_frac.numerator
    den = L_frac.denominator
    g = math.gcd(den, q)
    if g == 1:
        return (num * pow(den, -1, q)) % q
    v = 0
    d = den
    while d % p == 0:
        d //= p
        v += 1
    u = d
    if num % (p ** v) != 0:
        raise ValueError(f"L_p(1+{p}·{s}) not p-adic int: num={num}, den={den}, v={v}")
    num_reduced = num // (p ** v)
    return (num_reduced * pow(u, -1, q)) % q


def F_hat_full_at_pa(p, r, a, c=1):
    """F̂_p^full(p·a) via length-period DFT directly.
       G[a] = Σ_{s=0}^{period-1} f(s) · e_{period}(-a·s), where f(s) = e_q(c·(1+p)^s).
       Then F̂_full(p·a) = p · G[a]."""
    q = p ** (r + 1)
    period = p ** r
    base = (1 + p) % q
    pw = 1
    total = complex(0, 0)
    for s in range(period):
        phase = (c * pw) % q
        # multiply by e_{period}(-a·s) = exp(-2πi a s / period)
        total += cmath.exp(2j * cmath.pi * (phase / q - a * s / period))
        pw = (pw * base) % q
    return p * total  # F̂_full(p·a) = p · G[a]


def G_p_a(p, r, a, c=1):
    """Direct G[a] (length-period DFT of f at frequency a)."""
    q = p ** (r + 1)
    period = p ** r
    base = (1 + p) % q
    pw = 1
    total = complex(0, 0)
    for s in range(period):
        phase = (c * pw) % q
        total += cmath.exp(2j * cmath.pi * (phase / q - a * s / period))
        pw = (pw * base) % q
    return total


def one_hat_at_pa(p, r, a, N):
    """1̂(p·a) = Σ_{u=0}^{N-1} e_q(p·a·u), q = p^{r+1}.
       = Σ_u exp(2πi·p·a·u/q) = Σ_u exp(2πi·a·u/p^r)."""
    q = p ** (r + 1)
    period = p ** r  # = N for r-1 normalization... actually N = p^{r-1} per doc.
    # Wait: 1̂(ξ) := Σ_{u=0}^{N-1} e_q(ξ·u) with N = p^{r-1}.
    # So e_q(p·a·u) = e_{p^r}(a·u). Sum from u=0 to N-1 = p^{r-1}-1.
    total = complex(0, 0)
    for u in range(N):
        total += cmath.exp(2j * cmath.pi * a * u / (p ** r))
    return total


def L_tilde_setup(p, r):
    q = p ** (r + 1)
    J = J_for_p(p, r + 1)
    L1 = L_mod_q(p, 1, J, q)
    if L1 % p != 0:
        raise ValueError(f"L_p(1+p) at p={p} r={r} not divisible by p")
    L_tilde = (L1 // p) % (p ** r)
    if L_tilde % p != 1:
        raise ValueError(f"L̃_p ≢ 1 mod p at p={p} r={r}: L_tilde={L_tilde}")
    L_tilde_inv = pow(L_tilde, -1, p ** r)
    return J, L_tilde, L_tilde_inv


def verify_cell(p, r):
    """All checks at cell (p, r)."""
    q = p ** (r + 1)
    period = p ** r
    N = p ** (r - 1)
    sqrt_q = math.sqrt(q)
    pred_mag = p ** ((r + 1) / 2.0)
    J, L_tilde, L_tilde_inv = L_tilde_setup(p, r)

    # Support: a ∈ Z/p^r with a ≡ 1 mod p
    support = [a for a in range(period) if a % p == 1]
    assert len(support) == N

    # C1 + C2 + bilinear S_partial
    C_set = set()
    saddle_diffs = []
    saddle_ratios = []
    mag_devs = []
    S_partial = complex(0, 0)
    one_hats = []
    G_vals = []
    for a in support:
        C_a = (a * L_tilde_inv) % period
        C_set.add(C_a)
        G = G_p_a(p, r, a, c=1)
        G_vals.append(G)
        mag_devs.append(abs(abs(G) - pred_mag))
        oh = one_hat_at_pa(p, r, a, N)
        one_hats.append(oh)
        S_partial += oh * G

        # Saddle: s* = (C_a - 1)/p mod p
        s_star = ((C_a - 1) // p) % p
        L_at_sstar = L_mod_q(p, s_star, J, q)
        P_at_sstar = (p * s_star - C_a * L_at_sstar) % q
        psi_pred = cmath.exp(2j * cmath.pi * P_at_sstar / q)
        psi_emp = G / pred_mag
        saddle_diffs.append(abs(psi_emp - psi_pred))
        saddle_ratios.append(psi_emp / psi_pred)

    bijection_ok = (len(C_set) == len(support) and all((c_val % p == 1) for c_val in C_set))
    max_mag_dev = max(mag_devs)
    mag_ok = max_mag_dev < 1e-9 * pred_mag

    max_saddle_diff = max(saddle_diffs)
    saddle_ok_r3 = (r == 3 and max_saddle_diff < 1e-9)

    # r=2 uniform-Gaussian-factor check
    if r == 2:
        avg = sum(saddle_ratios) / len(saddle_ratios)
        ratio_std = max(abs(rv - avg) for rv in saddle_ratios)
        r2_uniform_factor = ratio_std < 1e-9
        avg_re, avg_im = avg.real, avg.imag
    else:
        r2_uniform_factor = None
        avg_re = avg_im = None

    abs_S = abs(S_partial)
    sqrtN = math.sqrt(N) if N > 0 else 1.0
    ratio_S_over_sqrtN = abs_S / sqrtN
    bound_2sqrtN = 2 * sqrtN
    bound_HenselLog = 2 * sqrtN * (1 + math.log(N)) if N > 1 else 2 * sqrtN

    return {
        'p': p, 'r': r, 'q': q, 'period': period, 'N': N, 'J': J,
        'supp_size': len(support),
        'bijection_ok': bijection_ok,
        'max_mag_dev': max_mag_dev, 'mag_ok': mag_ok,
        'max_saddle_diff': max_saddle_diff,
        'saddle_ok_r3': saddle_ok_r3 if r == 3 else None,
        'r2_uniform_factor': r2_uniform_factor,
        'r2_avg_re': avg_re, 'r2_avg_im': avg_im,
        'S_partial_abs': abs_S,
        'ratio_S_over_sqrtN': ratio_S_over_sqrtN,
        'bound_2sqrtN': bound_2sqrtN,
        'bound_HenselLog': bound_HenselLog,
        'sqrtN': sqrtN,
    }


def main():
    print("# Path 2 pushback: extended numerical verification")
    print()
    # Check 5 cells: p ∈ {13..31} × r ∈ {2,3}
    check5_cells = [(p, r) for p in [13, 17, 19, 23, 29, 31] for r in (2, 3)]
    # Check 3 cells: p ∈ {3,5,7} × r ∈ {4,5,6}
    check3_cells = [(p, r) for p in [3, 5, 7] for r in (4, 5, 6)]
    # Also redo p∈{3,5,7,11} r∈{2,3} for ratio_S_over_sqrtN baseline.
    baseline_cells = [(p, r) for p in [3, 5, 7, 11] for r in (2, 3)]

    all_cells = baseline_cells + check5_cells + check3_cells
    csv_rows = [['phase', 'p', 'r', 'q', 'period', 'N', 'J', 'supp_size',
                 'bijection_ok', 'max_mag_dev', 'mag_ok',
                 'max_saddle_diff', 'saddle_ok_r3', 'r2_uniform_factor',
                 'r2_avg_re', 'r2_avg_im',
                 'S_partial_abs', 'ratio_S_over_sqrtN',
                 'bound_2sqrtN', 'bound_HenselLog', 'sqrtN']]

    for (p, r) in all_cells:
        if (p, r) in baseline_cells:
            phase = 'baseline'
        elif (p, r) in check5_cells:
            phase = 'check5'
        else:
            phase = 'check3'
        print(f"## phase={phase}  p={p}  r={r}")
        try:
            res = verify_cell(p, r)
        except Exception as e:
            print(f"  ERROR: {e}")
            csv_rows.append([phase, p, r] + ['ERR'] * (len(csv_rows[0]) - 3))
            continue
        print(f"  q={res['q']} period={res['period']} N={res['N']} J={res['J']} |supp|={res['supp_size']}")
        print(f"  C1 bijection_ok={res['bijection_ok']}")
        print(f"  C2 max_mag_dev={res['max_mag_dev']:.2e} mag_ok={res['mag_ok']}")
        if r == 2:
            print(f"  r=2 uniform_factor={res['r2_uniform_factor']} avg=({res['r2_avg_re']:.6f}, {res['r2_avg_im']:.6f}i)")
        if r == 3:
            print(f"  C4 max_saddle_diff={res['max_saddle_diff']:.2e} saddle_ok_r3={res['saddle_ok_r3']}")
        if r >= 4:
            print(f"  (r≥4) max_saddle_diff={res['max_saddle_diff']:.4f} (Hensel deviation)")
        print(f"  |S_partial|={res['S_partial_abs']:.6f}  ratio/√N={res['ratio_S_over_sqrtN']:.4f}")
        print(f"  bound 2√N={res['bound_2sqrtN']:.4f}  bound 2√N·(1+log N)={res['bound_HenselLog']:.4f}")
        print()
        csv_rows.append([phase, p, r, res['q'], res['period'], res['N'], res['J'], res['supp_size'],
                         res['bijection_ok'], f"{res['max_mag_dev']:.6e}", res['mag_ok'],
                         f"{res['max_saddle_diff']:.6e}", res['saddle_ok_r3'], res['r2_uniform_factor'],
                         f"{res['r2_avg_re']:.6f}" if res['r2_avg_re'] is not None else '',
                         f"{res['r2_avg_im']:.6f}" if res['r2_avg_im'] is not None else '',
                         f"{res['S_partial_abs']:.6f}", f"{res['ratio_S_over_sqrtN']:.6f}",
                         f"{res['bound_2sqrtN']:.6f}", f"{res['bound_HenselLog']:.6f}",
                         f"{res['sqrtN']:.6f}"])

    import csv
    with open(CSV_PATH, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerows(csv_rows)
    print(f"# CSV written to {CSV_PATH}")


if __name__ == "__main__":
    main()
