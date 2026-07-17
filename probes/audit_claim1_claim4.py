"""
Audit script for CLAIM 1 and CLAIM 4 in milicevic_banks_verification.md.

CLAIM 1: F̂(3a) vs 1̂(3a) at r=3.
  - F̂(3a) = Σ_{u=0}^{q-1} e_q(c·4^u - 3a·u)   (full sum, c=1)
  - 1̂(3a) = Σ_{u=0}^{N-1} e_q(3a·u)            (indicator FT, no character)
  with q = 81, N = 3^{r-1} = 9 (per result_78_extended line 56:
    "1̂(3a) = Σ_{u=0}^{N−1} e_q(3au) (short-window character sum, N = 3^{r-1})")

CLAIM 4: P_a(s*) at r=3 for s* in {0,1,2}.
  P_a(s) = 3s - C_a · L(1+3s) mod 81.
  Doc claims:
    P_a(s*=0) = 0
    P_a(s*=1) = 3 - (15/2) C_a (mod 81)
    P_a(s*=2) = 6 - 60 C_a (mod 81)

  Also verify L̃ = L(4)/3 claim.
"""
import cmath, math
from fractions import Fraction

def trunc_log(s, J):
    L = Fraction(0)
    for j in range(1, J+1):
        L += Fraction((-1)**(j-1), j) * (Fraction(3*s)) ** j
    return L

def J_for_p3(m):
    j = 1
    while True:
        x = j+1
        v = 0
        while x % 3 == 0:
            x //= 3
            v += 1
        if (j+1) - v >= m:
            return j
        j += 1

r = 3
q = 3**(r+1)            # 81
N = 3**(r-1)            # 9
period = 3**r           # 27
m = r+1                 # 4
J = J_for_p3(m)         # should be 3 (per doc)

print(f"r={r}, q={q}, N={N}, period={period}, m={m}, J={J}")
supp = [a for a in range(period) if a % 3 == 1]
print(f"supp (a mod 3^r = {period}): {supp}")
print()

# CLAIM 1: F̂(3a) vs 1̂(3a)
print("=== CLAIM 1: F̂(3a) vs 1̂(3a) ===")
print()
print(f"{'a':>3} {'|F̂(3a)|':>12} {'|1̂(3a)|':>12} {'F̂':>30} {'1̂':>30}")
for a in supp:
    # F̂(3a) full
    Fhat = 0+0j
    pw = 1
    for u in range(q):
        ph = (1*pw - 3*a*u) % q
        Fhat += cmath.exp(2j*math.pi*ph/q)
        pw = (pw * 4) % q
    # 1̂(3a) over u=0..N-1
    Ihat = 0+0j
    for u in range(N):
        ph = (3*a*u) % q
        Ihat += cmath.exp(2j*math.pi*ph/q)
    print(f"{a:>3} {abs(Fhat):>12.6f} {abs(Ihat):>12.6f} "
          f"{Fhat.real:>+10.4f}{Fhat.imag:>+10.4f}i {Ihat.real:>+10.4f}{Ihat.imag:>+10.4f}i")

print()
print(f"3*sqrt(q) = {3*math.sqrt(q):.6f}")
print()

# CLAIM 4: saddle-class linear phases
print("=== CLAIM 4: P_a(s*) values mod 81 ===")
print()

# L(1) = L(1+3*0) = 0 (sum starts at j=1, with 3*0=0)
L1 = trunc_log(0, J)
L4 = trunc_log(1, J)
L7 = trunc_log(2, J)
print(f"L(1) = L(1+3*0) = {L1}  (should be 0)")
print(f"L(4) = L(1+3*1) = {L4} = {float(L4):.6f}")
print(f"L(7) = L(1+3*2) = {L7} = {float(L7):.6f}")
print()

# Reduce mod q
def red(F, mod):
    if F == 0:
        return 0
    return (F.numerator * pow(F.denominator, -1, mod)) % mod

L1_mod = red(L1, q)
L4_mod = red(L4, q)
L7_mod = red(L7, q)
print(f"L(1) mod {q} = {L1_mod}")
print(f"L(4) mod {q} = {L4_mod}")
print(f"L(7) mod {q} = {L7_mod}")
print()

# Doc says L̃ = L(4)/3 = 5/2.  Let's check:
# L(4) as Fraction is = ?
print(f"L(4) numerator/denominator: {L4.numerator}/{L4.denominator}")
print(f"L(4)/3 = {L4 / 3} (does this equal 5/2?)")
print(f"L̃ from script: L4_mod // 3 = {L4_mod // 3}")
print()

# Now compute P_a(s) = 3s - C_a · L(1+3s) symbolically in C_a
# P_a(0) = 0 - C_a · 0 = 0
# P_a(1) = 3 - C_a · L(4)
# P_a(2) = 6 - C_a · L(7)

# Check claim: P_a(1) = 3 - (15/2) C_a (mod 81)?
# That would mean L(4) ≡ 15/2 mod 81
# 15/2 mod 81: 2^{-1} mod 81 = 41 (since 2*41 = 82 ≡ 1)
# 15 * 41 mod 81 = 615 mod 81 = 615 - 7*81 = 615 - 567 = 48
print(f"15/2 mod 81 = {(15 * pow(2,-1,q)) % q}")
print(f"  ⟹ doc claim 'P_a(1) = 3 - (15/2)C_a mod 81' means L(4) ≡ {(15*pow(2,-1,q))%q} mod 81")
print(f"  Actual L(4) mod 81 = {L4_mod}")
print()

# Check claim: P_a(2) = 6 - 60 C_a (mod 81)?
# That means L(7) ≡ 60 mod 81
print(f"  Doc claim 'P_a(2) = 6 - 60·C_a mod 81' means L(7) ≡ 60 mod 81")
print(f"  Actual L(7) mod 81 = {L7_mod}")
print()

# Empirical verification: for each a, compute P_a(s) for s=s*(C_a) and check vs result_78_extended table
print("=== Verification: P_a(s*) per support a ===")
print(f"{'a':>3} {'C_a':>5} {'s*':>3} {'P_a(s*) computed':>18} {'doc 78 table':>15}")
# C_a from path B: L̃ * C_a ≡ a mod 3^{m-1} = 27
p_mm1 = 3**(m-1)  # 27
L_tilde = L4_mod // 3
L_tilde_inv = pow(L_tilde, -1, p_mm1)

# Doc Table at r=3 (from result_78_extended.md):
doc_table = {1:(22,1,0), 4:(7,2,72), 7:(19,0,0), 10:(4,1,54),
             13:(16,2,18), 16:(1,0,0), 19:(13,1,27),
             22:(25,2,45), 25:(10,0,0)}

for a in supp:
    C_a = (a * L_tilde_inv) % p_mm1
    # s*(C_a) = (C_a - 1)/3 mod 3 (per Theorem 78.6)
    s_star = ((C_a - 1) // 3) % 3
    # Compute P_a(s*) mod 81
    L_at_s = [L1_mod, L4_mod, L7_mod][s_star]
    P_val = (3*s_star - C_a * L_at_s) % q
    Cd, sd, Pd = doc_table[a]
    print(f"{a:>3} {C_a:>5} {s_star:>3} {P_val:>18} {Pd:>15} {'✓' if (C_a==Cd and s_star==sd and P_val==Pd) else '✗'}")

# Now express P_a as linear-in-C_a closed forms mod 81:
print()
print("=== Linear-in-C_a check ===")
print(f"P_a(0) = 0 (independent of C_a)  ✓")
print(f"P_a(1) = 3 - C_a · L(4) mod 81 = 3 - {L4_mod}·C_a mod 81")
print(f"P_a(2) = 6 - C_a · L(7) mod 81 = 6 - {L7_mod}·C_a mod 81")
print()
print(f"Doc claims:")
print(f"  P_a(1) = 3 - (15/2)·C_a mod 81 = 3 - 48·C_a mod 81")
print(f"  P_a(2) = 6 - 60·C_a mod 81")
