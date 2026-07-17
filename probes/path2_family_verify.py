"""
path2_family_verify.py — verify family-level extension of R78.4-78.6 at p ∈ {3,5,7,11}, r ∈ {2,3}.

Tests:
  C1: bijection a ↔ C_a (exhaustive)
  C2: |G_p(a)| = p^{(r+1)/2} for all a in support
  C3: at r=2, magnitude only; document phase
  C4: at r=3, saddle prediction G_p(a) = p^{(r+1)/2} · e_q(P_a(s*(C_a)))

Writes CSV per-cell + console table.
"""
import sys
import os
import math
import cmath
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz"
CSV_PATH = os.path.join(OUTDIR, "PATH2_FAMILY_EXTENSION_VERIFICATION.csv")


def J_for_p(p, m):
    """J = max j such that for all k ≤ j, k - v_p(k) < m. Returns largest such j."""
    j = 1
    while True:
        # check if (j+1) still satisfies (j+1) - v_p(j+1) < m
        x = j + 1
        v = 0
        while x % p == 0:
            x //= p
            v += 1
        if (j + 1) - v >= m:
            return j
        j += 1


def truncated_padic_log(p, s, J):
    """L_p(1+p·s) = Σ_{j=1}^J (-1)^{j-1}/j · (p·s)^j as Fraction."""
    L = Fraction(0)
    for j in range(1, J + 1):
        L += Fraction((-1) ** (j - 1), j) * Fraction(p * s) ** j
    return L


def L_mod_q(p, s, J, q):
    """Reduce L_p(1+p·s) mod q. Returns integer in [0, q)."""
    if s == 0:
        return 0
    L_frac = truncated_padic_log(p, s, J)
    num = L_frac.numerator
    den = L_frac.denominator
    # Need den invertible mod q: gcd(den, p) must be 1 (which holds since denoms are j ≤ J < p^? and we picked J such that no j divisible by p^large).
    # Actually denominators include j up to J. For p=3 J=3 r=2, denoms include 2, 3 — and 3 isn't invertible mod 3^k!
    # The Cochrane formula handles this by: the j with p|j contributes only via the v_p(j) reduction. Need careful handling.
    #
    # Approach: compute num · den^{-1} mod q, but den may share factor with p.
    # The truncated log's denominators are LCM of {1,2,...,J}. For p odd and J small, denom is coprime to p as long as no j ≤ J has p|j.
    # For p=3, J=3: j ∈ {1,2,3}, j=3 has v_3(j)=1, denom includes 3 → NOT coprime to 3.
    #
    # Resolution: the (p·s)^j factor has v_p = j·v_p(s) + j ≥ j. So at j=3 in p=3 case, (3s)^3 has v_3 ≥ 3.
    # The j=3 term is ± (3s)^3 / 3 = ± 9·s^3, which has v_3 ≥ 2. So it's a well-defined p-adic integer.
    #
    # Concretely: each term t_j = (-1)^{j-1}/j · (p·s)^j has v_p = j - v_p(j) ≥ 1 (by construction).
    # So the entire L is a p-adic integer; the rational form just has cancellations.
    #
    # Safe reduction mod q = p^{r+1}: compute each t_j as fraction, multiply out the denominator
    # if gcd(den, p) > 1, factor out the p-power.
    #
    # Cleanest: convert L to integer mod q via the following algorithm:
    # L_frac = num/den. Let v = v_p(den). Multiply both num and den by p^v won't work directly.
    # Better: compute L · q' where q' = q·den_unit. Use Fraction integer reduction.
    #
    # Practical: use pow(den, -1, q) but only if gcd(den, q) = 1.
    g = math.gcd(den, q)
    if g == 1:
        return (num * pow(den, -1, q)) % q
    # gcd > 1, meaning p | den. Need to factor it out.
    # Decompose den = p^v · u with gcd(u, p) = 1.
    v = 0
    d = den
    while d % p == 0:
        d //= p
        v += 1
    u = d
    # L = num / (p^v · u). For L to be a p-adic integer, p^v must divide num.
    # Check:
    if num % (p ** v) != 0:
        raise ValueError(f"L_p(1+{p}·{s}) is not a p-adic integer: num={num}, den={den}, v={v}")
    num_reduced = num // (p ** v)
    # Now L = num_reduced / u, u coprime to p, so invertible mod q.
    return (num_reduced * pow(u, -1, q)) % q


def G_p_a(p, r, C_a, c=1):
    """G_p(a) = Σ_{s=0}^{p^r - 1} e_q(P_a(s)) with P_a(s) = p·c·s − C_a · L_p(1+ps)."""
    q = p ** (r + 1)
    period = p ** r
    J = J_for_p(p, r + 1)
    total = complex(0, 0)
    for s in range(period):
        L_val = L_mod_q(p, s, J, q)
        phase_int = (p * c * s - C_a * L_val) % q
        total += cmath.exp(2j * cmath.pi * phase_int / q)
    return total


def F_hat_direct(p, r, ksi, c=1):
    """F̂_p(ξ) = Σ_{u=0}^{q-1} e_q(c·(1+p)^u − ξ·u) directly."""
    q = p ** (r + 1)
    total = complex(0, 0)
    pw = 1  # (1+p)^0
    base = 1 + p
    for u in range(q):
        phase_int = (c * pw - ksi * u) % q
        total += cmath.exp(2j * cmath.pi * phase_int / q)
        pw = (pw * base) % q
    return total


def verify_cell(p, r):
    """Verify family-level extension at (p, r). Returns dict of results."""
    c = 1
    q = p ** (r + 1)
    period = p ** r
    J = J_for_p(p, r + 1)
    p_r = period
    sqrt_q = math.sqrt(q)
    pred_mag = p ** ((r + 1) / 2)

    # Compute L_p(1+p) mod q, extract L̃_p
    L1_mod = L_mod_q(p, 1, J, q)
    # L_p(1+p) has v_p = 1 (since leading term is p), so L1_mod = p · L̃_p with L̃_p unit mod p^r.
    assert L1_mod % p == 0, f"L_p(1+p) at p={p}, r={r} not divisible by p: {L1_mod}"
    L_tilde = (L1_mod // p) % (p ** r)
    # Sanity: L̃_p ≡ 1 mod p
    assert L_tilde % p == 1, f"L̃_p at p={p}, r={r} not ≡ 1 mod p: {L_tilde}"
    L_tilde_inv = pow(L_tilde, -1, p ** r)

    # Support: a ∈ {a ∈ Z/p^r : a ≡ 1 mod p}
    support = [a for a in range(p_r) if a % p == 1]
    assert len(support) == p ** (r - 1)

    # C1: bijection a ↔ C_a
    C_set = set()
    rows = []
    for a in support:
        C_a = (a * L_tilde_inv) % (p ** r)
        C_set.add(C_a)
        rows.append((a, C_a))
    bijection_ok = (len(C_set) == len(support)) and all((c_val % p == 1) for c_val in C_set)

    # C2: |G_p(a)| = p^{(r+1)/2}
    mag_devs = []
    G_vals = {}
    F_hat_vals = {}
    for (a, C_a) in rows:
        G = G_p_a(p, r, C_a, c=c)
        G_vals[a] = G
        mag_devs.append(abs(abs(G) - pred_mag))
        # Cross-check vs direct F̂_p(p·a) = p · e_q(c) · G_p(a)
        F_hat = F_hat_direct(p, r, (p * a) % q, c=c)
        F_hat_vals[a] = F_hat
        F_pred = p * cmath.exp(2j * cmath.pi * c / q) * G
        # Difference
        diff = abs(F_hat - F_pred)
        # store
        rows_with_G = rows  # placeholder
    max_mag_dev = max(mag_devs)
    mag_ok = max_mag_dev < 1e-9 * pred_mag

    # C3 / C4: saddle prediction at r=2 (mag only) and r=3 (mag + phase)
    saddle_phase_matches = []
    saddle_phase_full_check = []
    for (a, C_a) in rows:
        G = G_vals[a]
        # s* = (C_a - 1)/p mod p
        s_star = ((C_a - 1) // p) % p
        # Predicted P_a(s*) mod q
        L_at_sstar = L_mod_q(p, s_star, J, q)
        P_at_sstar = (p * c * s_star - C_a * L_at_sstar) % q
        # ψ predicted = e_q(P_a(s*))
        psi_pred = cmath.exp(2j * cmath.pi * P_at_sstar / q)
        # ψ empirical = G / |G| (or G / pred_mag since |G|=pred_mag)
        psi_emp = G / pred_mag
        # Phase difference (in units of 2π)
        # For r=3 we expect exact match; for r=2 there's an extra Gaussian factor (per R78.6 at q=3)
        phase_diff_full = (cmath.phase(psi_emp) - cmath.phase(psi_pred)) % (2 * cmath.pi)
        if phase_diff_full > cmath.pi:
            phase_diff_full -= 2 * cmath.pi
        saddle_phase_matches.append(abs(psi_emp - psi_pred))
        saddle_phase_full_check.append({
            'a': a, 'C_a': C_a, 's_star': s_star,
            'P_at_sstar': P_at_sstar,
            'psi_pred': psi_pred,
            'psi_emp': psi_emp,
            'phase_diff_rad': phase_diff_full,
            'abs_diff': abs(psi_emp - psi_pred),
        })

    max_saddle_diff = max(saddle_phase_matches)
    # At r=3 expect tiny; at r=2 expect potentially non-zero
    saddle_ok_r3 = (r == 3 and max_saddle_diff < 1e-9)
    saddle_ok_r2_mag = (r == 2)  # we don't check phase at r=2

    # At r=2, check that the ratio ψ_emp / ψ_pred is the SAME constant for all a
    # (the "Gaussian factor" should be a-independent)
    if r == 2:
        ratios = [d['psi_emp'] / d['psi_pred'] for d in saddle_phase_full_check]
        # Compare std deviation of ratios
        avg_ratio = sum(ratios) / len(ratios)
        ratio_std = max(abs(r_val - avg_ratio) for r_val in ratios)
        r2_uniform_factor = ratio_std < 1e-9
    else:
        avg_ratio = None
        r2_uniform_factor = None

    return {
        'p': p, 'r': r, 'q': q, 'period': period, 'J': J,
        'support_size': len(support),
        'L_tilde': L_tilde, 'L_tilde_inv': L_tilde_inv,
        'bijection_ok': bijection_ok,
        'max_mag_dev': max_mag_dev,
        'mag_ok': mag_ok,
        'max_saddle_diff': max_saddle_diff,
        'saddle_ok_r3': saddle_ok_r3 if r == 3 else None,
        'r2_uniform_factor': r2_uniform_factor,
        'avg_ratio_r2': avg_ratio,
        'saddle_phase_data': saddle_phase_full_check,
    }


def main():
    print("# Phase 2: family-level extension empirical verification")
    print()
    cells = [(3, 2), (3, 3), (5, 2), (5, 3), (7, 2), (7, 3), (11, 2), (11, 3)]
    csv_rows = []
    csv_rows.append(['p', 'r', 'q', 'period', 'J', 'supp_size', 'bijection_ok', 'max_mag_dev', 'mag_ok',
                     'max_saddle_diff', 'saddle_ok_r3', 'r2_uniform_factor', 'r2_avg_ratio_re', 'r2_avg_ratio_im'])
    for (p, r) in cells:
        print(f"## p = {p}, r = {r}")
        try:
            res = verify_cell(p, r)
        except Exception as e:
            print(f"  ERROR: {e}")
            csv_rows.append([p, r, '', '', '', '', 'ERROR', '', '', '', '', '', '', ''])
            continue
        print(f"  q={res['q']}, period={res['period']}, J={res['J']}, |supp|={res['support_size']}")
        print(f"  L̃_p = {res['L_tilde']}, L̃_p^{{-1}} = {res['L_tilde_inv']}")
        print(f"  C1 bijection: {res['bijection_ok']}")
        print(f"  C2 |G|=p^{{(r+1)/2}}: max_mag_dev = {res['max_mag_dev']:.2e}, OK = {res['mag_ok']}")
        if r == 2:
            print(f"  C3 saddle (r=2): max diff = {res['max_saddle_diff']:.4f} (expected non-zero from Gaussian factor)")
            if res['r2_uniform_factor'] is not None:
                avg = res['avg_ratio_r2']
                print(f"    Uniform factor across all a: {res['r2_uniform_factor']}")
                if avg is not None:
                    print(f"    Avg ratio ψ_emp/ψ_pred = {avg.real:.6f} + {avg.imag:.6f}i (|ratio|={abs(avg):.6f})")
        elif r == 3:
            print(f"  C4 saddle (r=3): max diff = {res['max_saddle_diff']:.2e}, exact OK = {res['saddle_ok_r3']}")
            # Show per-a table
            print(f"    {'a':>4}  {'C_a':>6}  {'s*':>4}  {'P_a(s*)':>10}  {'psi_pred':>30}  {'psi_emp':>30}  {'|diff|':>10}")
            for d in res['saddle_phase_data']:
                pp = d['psi_pred']; pe = d['psi_emp']
                print(f"    {d['a']:>4}  {d['C_a']:>6}  {d['s_star']:>4}  {d['P_at_sstar']:>10}  "
                      f"{pp.real:>+12.4f}{pp.imag:>+12.4f}i  {pe.real:>+12.4f}{pe.imag:>+12.4f}i  {d['abs_diff']:>10.2e}")
        print()
        avg = res.get('avg_ratio_r2')
        csv_rows.append([p, r, res['q'], res['period'], res['J'], res['support_size'],
                         res['bijection_ok'], f"{res['max_mag_dev']:.6e}", res['mag_ok'],
                         f"{res['max_saddle_diff']:.6e}", res['saddle_ok_r3'],
                         res['r2_uniform_factor'],
                         f"{avg.real:.6f}" if avg is not None else '',
                         f"{avg.imag:.6f}" if avg is not None else ''])

    # Write CSV
    import csv
    with open(CSV_PATH, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerows(csv_rows)
    print(f"# CSV written to {CSV_PATH}")


if __name__ == "__main__":
    main()
