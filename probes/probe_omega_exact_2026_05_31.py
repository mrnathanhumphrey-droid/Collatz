"""
Compute |<omega, p_m>|^2 as EXACT RATIONALS for m=0,1,2.

omega = chi_4 = primitive order-4 character of (Z/17)*.
omega(sigma) values pre-computed via dlog_3.

|<omega, p_m>|^2 = Re^2 + Im^2 of the complex sum sum_sigma omega(sigma) p_m(sigma).

Both Re and Im are in Q (since p_m is rational and omega values are in {±1, ±i}).
Hence |<omega, p_m>|^2 ∈ Q.

Compare to c(m) (also exact rational) and to c_inf (50-digit reference).

Three possible outcomes:
  A. |<omega, p_m>|^2 = c(m) at every m → depth-wise exact identity
  B. |<omega, p_m>|^2 ≠ c(m) but converges to c_inf rapidly → asymptotic identity
  C. |<omega, p_m>|^2 far from c(m) and from c_inf → coincidence
"""
from __future__ import annotations
import sys, json
from fractions import Fraction
from mpmath import mp, mpf
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 60

# omega = chi_4 in our k-indexing (k=4 → ω(σ) = i^(dlog_3(σ) mod 4))
omega_re = {1: 1, 2: -1, 4: 1, 8: -1, 9: -1, 13: 1, 15: -1, 16: 1}
omega_im = {3: 1, 5: 1, 6: -1, 7: -1, 10: -1, 11: -1, 12: 1, 14: 1}

q = 17

# Load exact p_0, p_1, p_2
with open("C:/Collatz/pm_distributions_2026_05_31.json") as f:
    data = json.load(f)

p_0 = {int(s): Fraction(*v) for s, v in data["p_0"].items()}
p_1 = {int(s): Fraction(*v) for s, v in data["p_1_rational"].items()}
p_2 = {int(s): Fraction(*v) for s, v in data["p_2_rational"].items()}

# c(m) exact
c0_exact = Fraction(19, 127)
c1_num = 265011804960406635465672455997699
c1_den = 1730087916969634762193659498034425
c1_exact = Fraction(c1_num, c1_den)
# c(2) we don't have exact rational on hand, but can compute via Σ chi_2(σ) p_2(σ)

# Legendre chi_2 (= χ_L):
def chi_L(s):
    s %= q
    if s == 0: return 0
    return 1 if pow(s, (q-1)//2, q) == 1 else -1

# c(m) computed from p_m (will match exact reference)
def compute_c(p_m):
    return sum(Fraction(chi_L(s)) * p_m[s] for s in range(1, q))

c0_check = compute_c(p_0)
c1_check = compute_c(p_1)
c2_check = compute_c(p_2)

print("=== Verify c(m) from p_m ===")
print(f"  c(0) = {c0_check}  expected 19/127 = {c0_exact}")
print(f"  match: {c0_check == c0_exact}")
print(f"  c(1) = {c1_check}")
print(f"  match: {c1_check == c1_exact}")
print(f"  c(2) = {c2_check}")
print(f"  c(2) float = {float(c2_check)}")

# Compute <omega, p_m>:
def omega_inner_product(p_m):
    """Return (Re_part, Im_part) as exact Fractions."""
    Re = Fraction(0)
    Im = Fraction(0)
    for s in range(1, q):
        if s in omega_re:
            Re += omega_re[s] * p_m[s]
        if s in omega_im:
            Im += omega_im[s] * p_m[s]
    return (Re, Im)

print("\n=== Computing <omega, p_m> exactly ===")
for m, p_m, c_m in [(0, p_0, c0_exact), (1, p_1, c1_exact), (2, p_2, c2_check)]:
    Re, Im = omega_inner_product(p_m)
    mag_sq = Re * Re + Im * Im
    print(f"\n  m={m}:")
    print(f"    Re<omega, p_{m}> = {Re}  (= {float(Re):.15f})")
    print(f"    Im<omega, p_{m}> = {Im}  (= {float(Im):.15f})")
    print(f"    |<omega, p_{m}>|^2 = {mag_sq}")
    print(f"    |<omega, p_{m}>|^2 (float) = {float(mag_sq):.30f}")
    print(f"    c({m}) = {c_m}")
    print(f"    c({m}) (float)            = {float(c_m):.30f}")
    diff = mag_sq - c_m
    print(f"    DIFF |<omega>|^2 - c({m}) = {diff}")
    print(f"    DIFF (float)               = {float(diff):.2e}")
    print(f"    Equal as rationals: {mag_sq == c_m}")

# Also compare to c_inf
c_inf_ref = mpf("0.15298912060588517527891674877413229926086222622334")
print(f"\n=== Compare to c_inf = {c_inf_ref} ===")
for m, p_m in [(0, p_0), (1, p_1), (2, p_2)]:
    Re, Im = omega_inner_product(p_m)
    mag_sq = Re * Re + Im * Im
    mag_sq_mp = mpf(mag_sq.numerator) / mpf(mag_sq.denominator)
    diff = mag_sq_mp - c_inf_ref
    print(f"  m={m}: |<omega, p_{m}>|^2 - c_inf = {float(diff):+.6e}")

# Decay ratios
print("\n=== Decay ratios |<omega, p_m>|^2 - c_inf ===")
diffs = []
for m, p_m in [(0, p_0), (1, p_1), (2, p_2)]:
    Re, Im = omega_inner_product(p_m)
    mag_sq = Re * Re + Im * Im
    mag_sq_mp = mpf(mag_sq.numerator) / mpf(mag_sq.denominator)
    diffs.append(mag_sq_mp - c_inf_ref)

for i in range(1, len(diffs)):
    if diffs[i-1] != 0:
        print(f"  diff_{i}/diff_{i-1} = {float(diffs[i]/diffs[i-1]):+.6e}")

print("\n=== Decay ratios c(m) - c_inf ===")
cm_diffs = []
for c_m in [c0_exact, c1_exact, c2_check]:
    cm_mp = mpf(c_m.numerator) / mpf(c_m.denominator)
    cm_diffs.append(cm_mp - c_inf_ref)
for i in range(1, len(cm_diffs)):
    if cm_diffs[i-1] != 0:
        print(f"  diff_{i}/diff_{i-1} = {float(cm_diffs[i]/cm_diffs[i-1]):+.6e}")

print("\n=== Done ===")
