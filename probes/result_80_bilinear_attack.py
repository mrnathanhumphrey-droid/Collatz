"""
result_80_bilinear_attack.py — empirical bilinear-sum attack on eq 190 closure.

Goal: |Σ_{a ∈ supp} 1̂(3a) · ψ(a)| ≪ q^{1−δ}, equivalently |S_partial| ≪ q^{1/2−δ}.

After the Milićević / Banks-Shparlinski non-applicability finding
(milicevic_banks_verification.md), the path forward identified is:

(A) Partition Σ_a by saddle class j = s*(C_a) ∈ {0, 1, 2}
    — within each class at r = 3, ψ(a) is a linear additive character.
(B) Empirically measure scaling at r = 3, 4, 5, 6 to extract δ.
(C) Test whether the saddle-class partition gives within-class cancellation,
    cross-class cancellation, or neither.

Setup recap (from result_78_extended.md):
  q = 3^{r+1}, supp = {a ∈ Z/3^r : a ≡ 1 mod 3}, |supp| = 3^{r-1}
  C_a = a · L̃⁻¹ mod 3^r,  L̃ = L(4)/3
  s*(C_a) = (C_a − 1)/3 mod 3
  ψ(a) = F̂(3a) / (3√q),  empirically |ψ(a)| = 1
  1̂(3a) = Σ_{u=0}^{N−1} e_q(3au), N = 3^{r-1}

At r = 3, P_a(s*=j) is linear in C_a (and hence in a):
  j = 0:  P_a(0) = 0
  j = 1:  P_a(1) = 3 − (15/2) C_a   (mod 81)
  j = 2:  P_a(2) = 6 − 60 C_a       (mod 81)

So ψ(a)|_{class j} = e_q(P_a(j)) = e_q(3j − α_j · a)  with  α_j = L(1+3j)/L̃.
This is the structure to exploit at r = 3 (and to break/preserve at r ≥ 4).
"""
import sys
import os
import math
import cmath
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def J_for_p3(m):
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


def truncated_3adic_log(s, J):
    L = Fraction(0)
    for j in range(1, J + 1):
        L += Fraction((-1) ** (j - 1), j) * (Fraction(3 * s)) ** j
    return L


def F_hat_direct(r, a, c=1):
    """F̂(3a) = Σ_{u=0}^{q−1} e_q(c·4^u − 3a·u). Direct computation."""
    q = 3 ** (r + 1)
    total = complex(0, 0)
    pw = 1
    for u in range(q):
        phase_int = (c * pw - 3 * a * u) % q
        total += cmath.exp(2j * cmath.pi * phase_int / q)
        pw = (pw * 4) % q
    return total


def ind_hat(r, a):
    """1̂(3a) = Σ_{u=0}^{N−1} e_q(3a·u),  N = 3^{r-1}."""
    q = 3 ** (r + 1)
    N = 3 ** (r - 1)
    total = complex(0, 0)
    for u in range(N):
        total += cmath.exp(2j * cmath.pi * 3 * a * u / q)
    return total


def C_a_of(a, r):
    """C_a = a · L̃⁻¹ mod 3^r,  L̃ = L(4)/3."""
    q = 3 ** (r + 1)
    p_mm1 = 3 ** r
    J = J_for_p3(r + 1)
    L4_frac = truncated_3adic_log(1, J)
    L4_mod = (L4_frac.numerator * pow(L4_frac.denominator, -1, q)) % q
    L_tilde = L4_mod // 3
    L_tilde_inv = pow(L_tilde, -1, p_mm1)
    return (a * L_tilde_inv) % p_mm1


def saddle_class(a, r):
    """j = s*(C_a) = (C_a − 1)/3 mod 3."""
    C_a = C_a_of(a, r)
    return ((C_a - 1) // 3) % 3


def supp_iter(r):
    period = 3 ** r
    return [a for a in range(period) if a % 3 == 1]


def bilinear_total(r):
    """Σ_a 1̂(3a) · ψ(a) where ψ(a) = F̂(3a)/(3√q)."""
    q = 3 ** (r + 1)
    sqrt_q = math.sqrt(q)
    total = complex(0, 0)
    for a in supp_iter(r):
        ind = ind_hat(r, a)
        F = F_hat_direct(r, a)
        psi = F / (3 * sqrt_q)
        total += ind * psi
    return total


def bilinear_by_class(r):
    """Same sum but grouped by saddle class j = s*(C_a)."""
    q = 3 ** (r + 1)
    sqrt_q = math.sqrt(q)
    by_class = {0: complex(0, 0), 1: complex(0, 0), 2: complex(0, 0)}
    for a in supp_iter(r):
        j = saddle_class(a, r)
        ind = ind_hat(r, a)
        F = F_hat_direct(r, a)
        psi = F / (3 * sqrt_q)
        by_class[j] += ind * psi
    return by_class


def main():
    print("# Empirical bilinear-sum attack on eq 190")
    print("# |Σ_{a ∈ supp} 1̂(3a) · ψ(a)|, target ≪ q^{1−δ}")
    print()

    rs = [3, 4, 5, 6]
    rows = []

    print(f"  {'r':>2}  {'q':>8}  {'|supp|':>7}  {'|Σ_total|':>14}  "
          f"{'log_q(|Σ|)':>12}  "
          f"{'|Σ_j=0|':>12}  {'|Σ_j=1|':>12}  {'|Σ_j=2|':>12}  "
          f"{'sum |Σ_j|':>12}  {'cancel ratio':>14}")

    for r in rs:
        q = 3 ** (r + 1)
        N_supp = 3 ** (r - 1)
        total = bilinear_total(r)
        by_cls = bilinear_by_class(r)
        mag_total = abs(total)
        mag_by_cls = {j: abs(v) for j, v in by_cls.items()}
        sum_mags = sum(mag_by_cls.values())
        log_q_total = math.log(mag_total) / math.log(q) if mag_total > 0 else -math.inf
        cancel_ratio = mag_total / sum_mags if sum_mags > 0 else 0

        rows.append({
            "r": r, "q": q, "N_supp": N_supp,
            "mag_total": mag_total, "log_q_total": log_q_total,
            "by_cls": mag_by_cls, "sum_mags": sum_mags,
            "cancel_ratio": cancel_ratio,
            "by_cls_complex": by_cls,
        })

        print(f"  {r:>2}  {q:>8}  {N_supp:>7}  {mag_total:>14.4f}  "
              f"{log_q_total:>12.4f}  "
              f"{mag_by_cls[0]:>12.4f}  {mag_by_cls[1]:>12.4f}  {mag_by_cls[2]:>12.4f}  "
              f"{sum_mags:>12.4f}  {cancel_ratio:>14.4f}")

    print()
    print("# Empirical δ from log_q(|Σ|): saving δ is 1 − log_q(|Σ|)")
    print(f"  {'r':>2}  {'log_q(|Σ|)':>12}  {'δ_emp = 1 − that':>18}")
    for row in rows:
        delta_emp = 1.0 - row["log_q_total"]
        print(f"  {row['r']:>2}  {row['log_q_total']:>12.4f}  {delta_emp:>18.4f}")

    print()
    print("# Class-by-class scaling: log_q(|Σ_j|)")
    print(f"  {'r':>2}  {'log_q|Σ_0|':>12}  {'log_q|Σ_1|':>12}  {'log_q|Σ_2|':>12}")
    for row in rows:
        cs = []
        for j in (0, 1, 2):
            v = row["by_cls"][j]
            cs.append(math.log(v) / math.log(row["q"]) if v > 0 else -math.inf)
        print(f"  {row['r']:>2}  {cs[0]:>12.4f}  {cs[1]:>12.4f}  {cs[2]:>12.4f}")

    print()
    print("# Phases of class sums (look for cross-class cancellation pattern)")
    for row in rows:
        print(f"  r = {row['r']}:")
        for j in (0, 1, 2):
            z = row["by_cls_complex"][j]
            phase_norm = (cmath.phase(z) / (2 * math.pi)) % 1.0 if abs(z) > 1e-10 else 0.0
            print(f"    Σ_{j} = {z.real:>+10.4f}{z.imag:>+10.4f}i,  "
                  f"|·| = {abs(z):>10.4f},  phase/2π = {phase_norm:.6f}")
        z_tot = sum(row["by_cls_complex"].values())
        print(f"    sum  = {z_tot.real:>+10.4f}{z_tot.imag:>+10.4f}i,  |·| = {abs(z_tot):>10.4f}")

    print()
    print("# Trivial bound comparison: Σ_a |1̂(3a)| (no ψ-cancellation)")
    print(f"  {'r':>2}  {'q':>8}  {'Σ|1̂|':>14}  {'log_q(Σ|1̂|)':>14}  {'q^1':>10}  {'|Σ|/Σ|1̂|':>12}")
    for r, row in zip(rs, rows):
        q = row["q"]
        triv = sum(abs(ind_hat(r, a)) for a in supp_iter(r))
        log_q_triv = math.log(triv) / math.log(q) if triv > 0 else -math.inf
        ratio = row["mag_total"] / triv if triv > 0 else 0
        print(f"  {r:>2}  {q:>8}  {triv:>14.4f}  {log_q_triv:>14.4f}  {q:>10}  {ratio:>12.6f}")


if __name__ == "__main__":
    main()
