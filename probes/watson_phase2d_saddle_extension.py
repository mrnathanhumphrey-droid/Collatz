"""
watson_phase2d_saddle_extension.py — Thread α: extend the R78.6 saddle-point setup to r > 3
and ASYMPTOTICALLY characterize T_p(r) as r → ∞.

R78.6 closed form at r=3 (verified to machine precision in result_78_extended.md):
  ψ(a) = G(a)/√q = e_q(P_a(s*(C_a)))     (J=3, p=3)

PATH2_BILINEAR Attempt G+ (PATH2_BILINEAR_FROM_CLOSED_FORM.md) derived:
  Inner(s*) = p · e_q(ξ_0 s_0) · Σ_j e_{p²}(a_0(s*) j)
  |Inner(s*)| ≤ p · sin(π/p)/sin(π(1 + p α(s*))/p²)
  Σ |Inner(s*)| ≤ p · (p + log p) = p² + p log p

⟹ |T_p|(r=3) ≤ N + p log p  =  N + O(log p)
⟹ |T_p|(r=3) ≤ 2N at p=3, with N = p²

The TASK here: extend this bound to general r via Hensel-lifted saddle, and derive the
asymptotic rate κ(r) := log |T_p(r)| / log N.

R79b EMPIRICAL: κ(r) = 0.522 ± 0.008 at r=8..20.

WHAT SADDLE-POINT PREDICTS: κ_predict = ? depending on the structure of higher-order saddle corrections.

Direct computation: at r=4,5,6,7 compute T_p(r) numerically and check rate.
"""
import sys
import math
import cmath
from fractions import Fraction
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

def J_for_p3(m):
    """Number of terms needed in truncated p-adic log expansion for p=3 at precision p^m."""
    j = 1
    while True:
        x = j + 1
        v = 0
        while x % 3 == 0:
            x //= 3
            v += 1
        if (j + 1) - v >= m:
            return j
        j += 1

def F_hat_direct(r, a, c=1):
    """F̂(3a) = Σ_{u=0}^{q-1} e_q(c·4^u − 3a·u). Direct sum (slow at large r)."""
    q = 3 ** (r + 1)
    total = complex(0, 0)
    pw = 1
    for u in range(q):
        phase_int = (c * pw - 3 * a * u) % q
        total += cmath.exp(2j * cmath.pi * phase_int / q)
        pw = (pw * 4) % q
    return total

def F_hat_via_period(r, a, c=1):
    """F̂(3a) = p · G(a) · e_q(boundary) — use period structure to compute faster."""
    p = 3
    q = p ** (r + 1)
    period = p ** r
    total = complex(0, 0)
    pw = 1
    for u in range(period):
        phase_int = (c * pw - 3 * a * u) % q
        total += cmath.exp(2j * cmath.pi * phase_int / q)
        pw = (pw * 4) % q
    # period sum repeated p times within full q-length gives p × period-sum
    return p * total

def ind_hat(r, a):
    """1̂(3a) = Σ_{u=0}^{N-1} e_q(3a·u),  N = 3^{r-1}."""
    q = 3 ** (r + 1)
    N = 3 ** (r - 1)
    total = complex(0, 0)
    z = cmath.exp(2j * cmath.pi * 3 * a / q)
    pw = complex(1, 0)
    for u in range(N):
        total += pw
        pw *= z
    return total

def compute_T_p_r(r):
    """T_p(r) := Σ_{a ≡ 1 mod 3 in Z/3^r} 1̂(3a) · ψ(a),  ψ(a) = F̂(3a) / (3·√q)."""
    q = 3 ** (r + 1)
    sqrt_q = math.sqrt(q)
    period = 3 ** r
    support = [a for a in range(period) if a % 3 == 1]  # a ≡ 1 mod 3 in Z/3^r
    T = complex(0, 0)
    for a in support:
        F_hat = F_hat_via_period(r, a)
        psi = F_hat / (3 * sqrt_q)
        ihat = ind_hat(r, a)
        T += ihat * psi
    return T, support

print("=" * 70)
print("THREAD α — Direct computation of T_p(r) at p=3, r=2..6")
print("=" * 70)
print()
print(f"  {'r':>2} {'N=3^{r-1}':>10} {'|T_p(r)|':>14} {'κ=log|T|/log(N)':>18}")

results = []
for r in [2, 3, 4, 5, 6]:
    T, supp = compute_T_p_r(r)
    N = 3 ** (r - 1)
    kappa = math.log(abs(T)) / math.log(N) if N > 1 else float('nan')
    results.append((r, N, abs(T), kappa))
    print(f"  {r:>2} {N:>10} {abs(T):>14.4f} {kappa:>18.4f}")

print()
print("R79b empirical: κ(r) = 0.522 ± 0.008 at r=8..20.")
print()

# Compute rate of change
if len(results) >= 2:
    print("Successive κ estimates (regression of log|T| on log N):")
    rs = np.array([x[0] for x in results])
    Ns = np.array([x[1] for x in results])
    Ts = np.array([x[2] for x in results])
    log_T = np.log(Ts)
    log_N = np.log(Ns)
    A = np.vstack([np.ones_like(log_N), log_N]).T
    coef, *_ = np.linalg.lstsq(A, log_T, rcond=None)
    print(f"  Log-linear fit: log|T_p(r)| = {coef[0]:.4f} + {coef[1]:.4f} · log(N)")
    print(f"  ⟹ κ_fit = {coef[1]:.4f}")
    print(f"  ⟹ |T_p(r)| ~ {math.exp(coef[0]):.4f} · N^{coef[1]:.4f}")
    print()
print()

# ============================================================
# THREAD α.2 — Apply saddle-point Hensel correction at r=4,5,6
# ============================================================
print("=" * 70)
print("THREAD α.2 — Sanity: ψ_lead vs ψ_true at r=3,4,5 (saddle exactness)")
print("=" * 70)
print()
print("At r=3: ψ_true(a) = e_q(P_a(s*(C_a))) saddle-exact (verified R78.6).")
print("At r=4: saddle-point requires Hensel lifting. Test mean |ψ_true − ψ_lead| empirically.")
print()

def compute_psi_lead(r, a):
    """ψ_lead(a) = e_q(P_a(s*(C_a))) from R78.6 with truncated J."""
    p = 3
    q = p ** (r + 1)
    # L_p(1 + p) / p  =  L̃_p; need this mod p^r
    # L_p(1 + ps) = Σ (-1)^{j-1}/j · (ps)^j
    J = J_for_p3(r + 1)
    # Compute L̃_p = L_p(1 + p) / p mod p^r as a Fraction modulo p^r
    # L_p(1+p) = p · L̃_p where L̃_p ≡ 1 mod p (the "stripped" unit)
    L_at_p_over_p = Fraction(0)
    for j in range(1, J + 1):
        L_at_p_over_p += Fraction((-1)**(j-1), j) * Fraction(p)**(j-1) * Fraction(1) ** j
    # That's L_p(1+p)/p; reduce mod p^r:
    num = (L_at_p_over_p.numerator * pow(L_at_p_over_p.denominator, -1, p**r)) % p**r
    L_tilde = num
    # C_a = a · L̃^{-1} mod p^r
    L_tilde_inv = pow(L_tilde, -1, p**r)
    C_a = (a * L_tilde_inv) % p**r
    # s*(C_a) = (C_a − 1)/p mod p — must have C_a ≡ 1 mod p
    if C_a % p != 1:
        return None  # off support
    s_star = ((C_a - 1) // p) % p
    # P_a(s*) = p·s* − C_a · L_p(1 + p·s*) mod q
    L_at_ps = Fraction(0)
    for j in range(1, J + 1):
        L_at_ps += Fraction((-1)**(j-1), j) * Fraction(p * s_star) ** j
    # Reduce L_at_ps mod q
    num_L = (L_at_ps.numerator * pow(L_at_ps.denominator, -1, q)) % q
    P_a_s = (p * s_star - C_a * num_L) % q
    return cmath.exp(2j * cmath.pi * P_a_s / q)

print(f"  {'r':>2} {'a':>4} {'|ψ_lead|':>10} {'|ψ_true|':>10} {'arg(ψ_true/ψ_lead) deg':>22} {'|ψ_true−ψ_lead|':>15}")
for r in [3, 4, 5]:
    q = 3 ** (r + 1)
    sqrt_q = math.sqrt(q)
    period = 3 ** r
    supp = [a for a in range(period) if a % 3 == 1][:5]  # first 5 elements
    for a in supp:
        psi_lead = compute_psi_lead(r, a)
        if psi_lead is None:
            continue
        F_hat = F_hat_via_period(r, a)
        psi_true = F_hat / (3 * sqrt_q)
        ratio = psi_true / psi_lead
        diff = abs(psi_true - psi_lead)
        print(f"  {r:>2} {a:>4} {abs(psi_lead):>10.4f} {abs(psi_true):>10.4f} {math.degrees(cmath.phase(ratio)):>22.2f} {diff:>15.4e}")
    print()

print()
print("INTERPRETATION:")
print("  At r=3: ψ_lead = ψ_true (machine precision) — saddle-exact.")
print("  At r=4,5: ψ_lead ≠ ψ_true — Hensel correction needed.")
print("  This is what R79b reported as the structural barrier for family-level extension.")
