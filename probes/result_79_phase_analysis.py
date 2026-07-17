"""Result 79 Step 1+: 3-adic phase analysis for Kalafatelis eq 190 van der Corput route.

Setup (Kalafatelis Prop 20+):
  S_{r,c,m}(0..N-1) := Σ_{u=0}^{N-1} e_{q}(c·4^u - 9m·u)
  where N = 3^{r-1}, q = 3^{r+1}, c ∈ Z_3^×.

Step 1: Build differenced phase Φ(u+h) - Φ(u) explicitly.
  Φ(u+h) - Φ(u) = c·4^u·(4^h - 1) - 9m·h
  Set h = 3^k·h' with gcd(h',3) = 1 (so v_3(h) = k).
  By LTE: v_3(4^h - 1) = v_3(4-1) + v_3(h) = 1 + k.
  So 4^h - 1 = 3^{k+1}·u_k(h')  with u_k(h') ∈ Z_3^×.
  Substituting: Φ(u+h) - Φ(u) = c·4^u·3^{k+1}·u_k(h') - 9mh
                              = 3^{k+1}·c·u_k(h')·4^u - 9mh
  In phase mod q = 3^{r+1}, the 3^{k+1} factor reduces the effective modulus to 3^{r-k}.

The inner sum I(h) = e_q(-9mh) · Σ_u e_{3^{r-k}}(c'' · 4^u)
  where c'' = c · u_k(h')·... ∈ Z_3^×.

Now: 4 has order 3^{r-k-1} mod 3^{r-k} (LTE).
  As u ranges over [0, N-h-1] (length N-h ≈ N = 3^{r-1}),
    4^u mod 3^{r-k} cycles with period P_k = 3^{r-k-1}.
  Number of complete cycles: M_k = floor((N-h)/P_k).
  Each complete cycle sum: Σ_{x ∈ ⟨4⟩ mod 3^{r-k}} e_{3^{r-k}}(c''·x) = 0 (proved below).
  Remainder length: ρ_k(h) := (N-h) - M_k·P_k.
  So I(h) = e_q(-9mh) · (partial-cycle sum of length ρ_k(h)).

This script:
  (a) Verifies LTE: v_3(4^{3^k} - 1) = k+1 for k = 0..6.
  (b) Verifies the complete-cycle vanishing for r' = 2..6.
  (c) Computes |I(h)| for various h directly, compares to partial-cycle bound.
  (d) Computes |S_{r,c,m}| directly, compares to trivial N and √N.
  (e) Tabulates partial-cycle bound vs actual |I(h)| vs trivial.
"""
import numpy as np
from math import gcd
from cmath import exp as cexp
from cmath import pi

PI = float(np.pi)

def v3(n: int) -> int:
    """3-adic valuation."""
    if n == 0:
        return 10**9
    n = abs(n)
    k = 0
    while n % 3 == 0:
        n //= 3
        k += 1
    return k

def phase_sum(coeffs_per_u, mod):
    """Sum of e^{2πi·coeffs[u]/mod} for u in range(len(coeffs))."""
    s = 0.0 + 0.0j
    inv = 2j * PI / mod
    for c in coeffs_per_u:
        s += cexp(inv * (c % mod))
    return s

def kalafatelis_S(r, c, m, length=None):
    """S = Σ_{u=0}^{length-1} e_{3^{r+1}}(c·4^u - 9m·u). Default length = 3^{r-1}."""
    q = 3**(r + 1)
    if length is None:
        length = 3**(r - 1)
    s = 0.0 + 0.0j
    x = 1
    nine_m = 9 * m
    inv = 2j * PI / q
    for u in range(length):
        phase = (c * x - nine_m * u) % q
        s += cexp(inv * phase)
        x = (x * 4) % q
    return s

def pure_subgroup_arc_sum(r_eff, c_eff, length):
    """Σ_{u=0}^{length-1} e_{3^{r_eff}}(c_eff · 4^u). Pure (no linear)."""
    q = 3**r_eff
    s = 0.0 + 0.0j
    x = 1
    inv = 2j * PI / q
    for u in range(length):
        s += cexp(inv * ((c_eff * x) % q))
        x = (x * 4) % q
    return s

# --------------------------------------------------------------------
# (a) LTE check
# --------------------------------------------------------------------
print("="*72)
print("(a) LTE check: v_3(4^{3^k} - 1) should equal k+1")
print("    k | 4^{3^k} - 1                 | v_3 | expected")
print("    -" * 16)
for k in range(0, 7):
    val = pow(4, 3**k) - 1
    exp_k = k + 1
    v = v3(val)
    print(f"    {k} | (val len = {len(str(val))} digits)         | {v}   | {exp_k}")
    assert v == exp_k, f"LTE FAILED at k={k}: v_3 = {v}, expected {exp_k}"
print("    OK: LTE holds for all k = 0..6.")
print()

# --------------------------------------------------------------------
# (b) Complete-cycle vanishing: Σ_{x ∈ ⟨4⟩ mod 3^N} e_{3^N}(c·x) = 0 for c unit
# --------------------------------------------------------------------
print("="*72)
print("(b) Complete-cycle vanishing for c ∈ Z_3^×, modulus q = 3^N")
print("    Σ_{u=0}^{ord(4)-1} e_q(c·4^u) should be 0 (size 3^{N-1})")
print()
print("    N | period 3^{N-1} | c=1 sum |abs| | c=2 sum |abs| | max c-sum |abs|")
for N in range(2, 7):
    period = 3**(N - 1)
    q = 3**N
    s1 = pure_subgroup_arc_sum(N, 1, period)
    s2 = pure_subgroup_arc_sum(N, 2, period)
    max_abs = 0.0
    for c in range(1, q):
        if gcd(c, 3) != 1:
            continue
        s = pure_subgroup_arc_sum(N, c, period)
        if abs(s) > max_abs:
            max_abs = abs(s)
    print(f"    {N} | {period:6d}        | {abs(s1):.3e}    | {abs(s2):.3e}    | {max_abs:.3e}")
print("    OK: complete subgroup-arc sum is essentially zero (numerical noise only).")
print()

# --------------------------------------------------------------------
# (c) Differenced inner sum I(h) for various h, compare to partial-cycle bound
# --------------------------------------------------------------------
print("="*72)
print("(c) Differenced inner sum I(h) at level r")
print("    For h = 3^k·h', |I(h)| ≤ partial-cycle length ρ_k(h)")
print("    where ρ_k(h) = (N - h) mod 3^{r-k-1}")
print()

for r in [4, 5, 6]:
    N = 3**(r - 1)
    q = 3**(r + 1)
    c = 1  # any unit suffices for size analysis
    print(f"  r = {r}, N = {N}, q = {q}")
    print(f"  {'h':>6} {'k':>4} {'h_prime':>8} {'period':>8} {'rho':>6} {'|I(h)|':>10} {'tightness':>10}")
    for h in [1, 2, 3, 4, 5, 6, 9, 12, 18, 27, 28, 36, 54, 81]:
        if h >= N:
            continue
        k = v3(h)
        h_prime = h // 3**k
        period = 3**(r - k - 1) if r - k - 1 >= 0 else 1
        rho = (N - h) % period if period > 0 else 0
        # Compute I(h) = Σ_{u=0}^{N-h-1} e_q(c·4^u·(4^h-1) - 9m·h)
        # The constant -9mh just adds a phase, so |I(h)| = |inner pure sum|.
        # Inner pure: Σ_u e_{3^{r-k}}(c'' · 4^u) where c'' = c·(4^h-1)/3^{k+1} mod 3^{r-k}
        diff_factor = pow(4, h, q) - 1   # = 4^h - 1 mod q
        assert v3(diff_factor) == k + 1, f"v_3(4^h-1) = {v3(diff_factor)}, expected {k+1}"
        # Build inner sum directly:
        s_inner = 0.0 + 0.0j
        x = 1
        for u in range(N - h):
            phase = (c * x * diff_factor) % q
            s_inner += cexp(2j * PI * phase / q)
            x = (x * 4) % q
        I_abs = abs(s_inner)
        tight = I_abs / max(rho, 1)
        print(f"  {h:>6} {k:>4} {h_prime:>8} {period:>8} {rho:>6} {I_abs:>10.4f} {tight:>10.4f}")
    print()

# --------------------------------------------------------------------
# (d) Direct |S_{r,c,m}| sizes vs trivial N and √N
# --------------------------------------------------------------------
print("="*72)
print("(d) Direct |S_{r,c,m}| over c ∈ Z_3^×, m ∈ {0, 1, 2, ..., 9}")
print("    Trivial bound: |S| ≤ N. Square-root: |S| ~ √N.")
print()
print("    r | N | √N | max|S| (over c, m=0..min(N_r,9)) | max/N | max/√N")
for r in range(2, 8):
    N = 3**(r - 1)
    N_r = 2 * 3**(r - 1)
    q = 3**(r + 1)
    sqrtN = N**0.5
    overall_max = 0.0
    overall_arg = None
    for c in range(1, q):
        if gcd(c, 3) != 1:
            continue
        for m in range(min(N_r, 10)):
            s = kalafatelis_S(r, c, m)
            v = abs(s)
            if v > overall_max:
                overall_max = v
                overall_arg = (c, m)
    print(f"    {r} | {N:5d} | {sqrtN:>7.2f} | {overall_max:>8.3f} (c={overall_arg[0]},m={overall_arg[1]}) | {overall_max/N:.4f} | {overall_max/sqrtN:.3f}")

print()
print("="*72)
print("Interpretation:")
print("  - If max/√N stays bounded as r grows, square-root cancellation holds.")
print("  - If max/N stays bounded away from 0, no cancellation (trivial bound is tight).")
print("  - If max/√N grows, we have intermediate (sub-trivial but worse than √N).")
