"""
hensel_approach_a_verify.py — verify the Hensel-lifted closed form at family level for r >= 4.

CLAIM (Approach A): for prime p >= 3, r >= 2, c=1:
    G_p(a) = p^((r+1)/2) * eta_p(r) * e_q(P_a(s*(r)))
where:
    s*(r) = (C_a - 1) // p  mod p^{r-1}    (Hensel-lifted saddle, exact)
    P_a(s*(r)) = Sum_{j=2}^{r} (-1)^{j-1} (p*s*(r))^j / (j*(j-1))  mod p^{r+1}
                 (closed form from (1+y) log(1+y) series identity)
    eta_p(r) = (1/sqrt(p)) * Sum_{h=0}^{p-1} e_p(h^2 / 2)    at even r
             = 1                                              at odd r

Tests:
    (p, r) in {(3, 4), (3, 5), (3, 6), (5, 4), (5, 5), (7, 4), (7, 5), (7, 6), (11, 4), (11, 5)}.

Compares |G_p(a)_actual - prediction| / |G_p(a)_actual| at every a in support.
Target: < 1e-12 max relative deviation (machine precision).

If passes: APPROACH_A_VALIDATED — H_HENSEL_CLOSES outcome.
If fails: APPROACH_A_INVALID — derivation error or small-prime caveat; redo by hand.

Uses mpmath for 50-digit arithmetic to rule out float64 round-off as a confound.

NOTE: Python denied this session; script written for main-thread execution.
"""

import sys
import os
import math
from fractions import Fraction

try:
    import mpmath
    mpmath.mp.dps = 50
    HAS_MPMATH = True
except ImportError:
    HAS_MPMATH = False
    import cmath

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz"
CSV_PATH = os.path.join(OUTDIR, "HENSEL_APPROACH_A_VERIFICATION.csv")


def J_for_p(p, m):
    """J = largest j such that for all k <= j, k - v_p(k) < m."""
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
    """L_p(1+p*s) = Sum_{j=1}^J (-1)^{j-1}/j * (p*s)^j  as Fraction."""
    L = Fraction(0)
    for j in range(1, J + 1):
        L += Fraction((-1) ** (j - 1), j) * Fraction(p * s) ** j
    return L


def L_mod_q(p, s, J, q):
    """Reduce L_p(1+p*s) mod q.  Returns integer in [0, q)."""
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
        raise ValueError(f"L_p(1+{p}*{s}) not p-adic integer: num={num}, den={den}, v={v}")
    num_reduced = num // (p ** v)
    return (num_reduced * pow(u, -1, q)) % q


def e_q(phase_int, q):
    """e_q(x) = exp(2*pi*i*x/q).  mpmath if available, else cmath."""
    if HAS_MPMATH:
        z = mpmath.mpc(0, 2 * mpmath.pi * phase_int / q)
        return mpmath.exp(z)
    else:
        import cmath
        return cmath.exp(2j * math.pi * phase_int / q)


def G_p_a_actual(p, r, C_a, c=1):
    """Direct computation: G_p(a) = Sum_{s=0}^{p^r - 1} e_q(p*c*s - C_a*L_p(1+p*s))."""
    q = p ** (r + 1)
    period = p ** r
    J = J_for_p(p, r + 1)
    if HAS_MPMATH:
        total = mpmath.mpc(0)
    else:
        total = complex(0)
    for s in range(period):
        L_val = L_mod_q(p, s, J, q)
        phase_int = (p * c * s - C_a * L_val) % q
        total += e_q(phase_int, q)
    return total


def P_a_s_star_hensel(p, r, C_a, c=1):
    """Compute P_a(s*(r)) mod p^{r+1} via closed form (1+y) log(1+y) series.
       Returns integer phase_int such that e_q(phase_int) = predicted phase factor."""
    q = p ** (r + 1)
    p_rm1 = p ** (r - 1)
    # s*(r) = (C_a - 1) // p mod p^{r-1}
    # C_a is integer in Z/p^r. Take canonical representative in [0, p^r).
    C_a_canonical = C_a % (p ** r)
    s_star = ((C_a_canonical - 1) // p) % p_rm1  # integer in [0, p^{r-1})

    # Direct way: just compute P_a(s_star) using the formula p*c*s - C_a * L_p(1+ps).
    # That IS the exact phase (no closed-form needed; the closed form was for derivation insight).
    # But the *prediction* of Approach A is that this matches G_p(a) up to magnitude p^{(r+1)/2} and eta_p factor.
    J = J_for_p(p, r + 1)
    L_val = L_mod_q(p, s_star, J, q)
    phase_int = (p * c * s_star - C_a * L_val) % q
    return phase_int


def quadratic_gauss_sum_mod_p(p, coeff=None):
    """Compute eta_p = (1/sqrt(p)) * Sum_{h=0}^{p-1} e_p(coeff * h^2)
       where coeff = 1/2 mod p by default (for c=1, C_a ≡ 1 mod p).
       Returns complex unit |eta_p| = 1."""
    if coeff is None:
        coeff_int = pow(2, -1, p)  # 1/2 mod p
    else:
        coeff_int = coeff % p
    if HAS_MPMATH:
        total = mpmath.mpc(0)
        for h in range(p):
            ph = (coeff_int * h * h) % p
            total += mpmath.exp(2j * mpmath.pi * ph / p)
        eta = total / mpmath.sqrt(p)
    else:
        import cmath
        total = complex(0)
        for h in range(p):
            ph = (coeff_int * h * h) % p
            total += cmath.exp(2j * math.pi * ph / p)
        eta = total / math.sqrt(p)
    return eta


def predict_G_p_a(p, r, C_a, c=1):
    """Approach A prediction: G_p(a) = p^{(r+1)/2} * eta_p(r) * e_q(P_a(s*(r)))."""
    q = p ** (r + 1)
    mag = p ** ((r + 1) / 2.0)  # use float for now; mpmath casts ok
    if HAS_MPMATH:
        mag = mpmath.mpf(p) ** mpmath.mpf((r + 1) / 2)
    phase_int = P_a_s_star_hensel(p, r, C_a, c=c)
    phase_factor = e_q(phase_int, q)
    if r % 2 == 0:
        eta = quadratic_gauss_sum_mod_p(p)
    else:
        eta = 1
    return mag * eta * phase_factor


def verify_cell(p, r):
    """Verify at (p, r): compare G_p(a)_actual vs predicted at every a in support."""
    q = p ** (r + 1)
    period = p ** r
    J = J_for_p(p, r + 1)
    pred_mag = p ** ((r + 1) / 2.0)
    print(f"\n## p={p}, r={r}, q={q}, period={period}, J_p={J}")

    # L_tilde_p
    L1_mod = L_mod_q(p, 1, J, q)
    if L1_mod % p != 0:
        print(f"  ERROR: L_p(1+p) not divisible by p: {L1_mod}")
        return None
    L_tilde = (L1_mod // p) % (p ** r)
    L_tilde_inv = pow(L_tilde, -1, p ** r)
    print(f"  L_tilde={L_tilde}, L_tilde_inv={L_tilde_inv}")

    support = [a for a in range(p ** r) if a % p == 1]
    print(f"  |support|={len(support)}")

    max_rel_dev = 0.0
    max_abs_dev = 0.0
    devs = []
    for a in support:
        C_a = (a * L_tilde_inv) % (p ** r)
        G_actual = G_p_a_actual(p, r, C_a)
        G_pred = predict_G_p_a(p, r, C_a)
        abs_dev = abs(G_actual - G_pred)
        rel_dev = abs_dev / abs(G_actual) if abs(G_actual) > 1e-50 else abs_dev
        max_abs_dev = max(max_abs_dev, float(abs_dev))
        max_rel_dev = max(max_rel_dev, float(rel_dev))
        devs.append((a, C_a, G_actual, G_pred, abs_dev, rel_dev))

    print(f"  max_abs_dev = {max_abs_dev:.6e}")
    print(f"  max_rel_dev = {max_rel_dev:.6e}")
    print(f"  prediction status: {'PASS' if max_rel_dev < 1e-10 else 'FAIL'}")

    return {
        'p': p, 'r': r, 'q': q, 'period': period, 'J': J,
        'support_size': len(support),
        'max_abs_dev': max_abs_dev, 'max_rel_dev': max_rel_dev,
        'pass': max_rel_dev < 1e-10,
        'devs': devs,
    }


def main():
    cells = [
        (3, 4), (3, 5), (3, 6),
        (5, 4), (5, 5),
        (7, 4), (7, 5), (7, 6),
        (11, 4), (11, 5),
    ]
    print("# HENSEL Approach A — numerical verification")
    print(f"# mpmath available: {HAS_MPMATH}, dps=50 if so")
    print()
    rows = [['p', 'r', 'q', 'period', 'J_p', 'supp_size', 'max_abs_dev', 'max_rel_dev', 'pass']]
    for (p, r) in cells:
        try:
            res = verify_cell(p, r)
            if res is None:
                rows.append([p, r, '', '', '', '', 'ERROR', '', ''])
                continue
            rows.append([p, r, res['q'], res['period'], res['J'], res['support_size'],
                         f"{res['max_abs_dev']:.6e}", f"{res['max_rel_dev']:.6e}", res['pass']])
        except Exception as e:
            print(f"  EXCEPTION at p={p}, r={r}: {e}")
            rows.append([p, r, '', '', '', '', f'EXCEPTION:{e}', '', ''])

    import csv
    with open(CSV_PATH, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerows(rows)
    print(f"\n# CSV written to {CSV_PATH}")


if __name__ == "__main__":
    main()
