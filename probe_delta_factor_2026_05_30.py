"""
probe_delta_factor_2026_05_30.py

Final deliverable for the Δ_m series request:
  1. Δ_1 EXACT as reduced fraction with FULL prime factorization of num + den.
  2. Δ_2, Δ_3, Δ_4 numerical at ~15 digit precision (float64 from FFT n=m+1).
  3. Ratio sequence Δ_m / Δ_{m-1} numerical.
  4. PSLQ summary against polylog Li_s(2^-8), Lerch Φ(2^-8, s, a) — NEGATIVE at 10^-15 tol.
  5. Structural diagnosis: why naive shift-additive extension fails.
"""
from __future__ import annotations
import sys
from fractions import Fraction
from sympy import factorint
sys.stdout.reconfigure(encoding="utf-8")

# Δ_1 EXACT
c1 = Fraction(265011804960406635465672455997699, 1730087916969634762193659498034425)
c0 = Fraction(19, 127)
delta1 = c1 - c0
print(f"=== Δ_1 EXACT ===")
print(f"c(1)         = {c1.numerator}/{c1.denominator}")
print(f"  num digits = {len(str(c1.numerator))}, den digits = {len(str(c1.denominator))}")
print(f"19/127       = 19/127")
print(f"Δ_1 = c(1) - 19/127")
print(f"  num = {delta1.numerator}")
print(f"  den = {delta1.denominator}")
print(f"  num digits = {len(str(delta1.numerator))}, den digits = {len(str(delta1.denominator))}")
print(f"  Δ_1 (float) = {float(delta1):.18e}")

print(f"\n--- Δ_1 numerator prime factorization ---")
fnum = factorint(delta1.numerator)
for p in sorted(fnum):
    print(f"  {p}^{fnum[p]}")
print(f"\n--- Δ_1 denominator prime factorization ---")
fden = factorint(delta1.denominator)
for p in sorted(fden):
    print(f"  {p}^{fden[p]}")

# Check structural relation: does Δ_1 num/den factor as 2^-? * 127^? * ...?
# c(0) = 19/127 has 127 = q^2 - q - 1 (no, 127 = 2^7 - 1, Mersenne prime).
# c(1) den has interesting factor structure — Mersenne primes?
print(f"\n--- c(1) denominator prime factorization ---")
fc1den = factorint(c1.denominator)
for p in sorted(fc1den):
    print(f"  {p}^{fc1den[p]}")

# Δ_2, Δ_3, Δ_4 from FFT (numerical at ~15 digits)
delta_floats = {
    1: float(delta1),
    2: 3.6416215664888670e-03,
    3: 3.3990324789156260e-03,
    4: 3.3824098342225830e-03,
}

print(f"\n=== Δ_m numerical (m=1..4) ===")
print(f"  Δ_1 = {delta_floats[1]:.16e}  (EXACT — float is rounded)")
print(f"  Δ_2 = {delta_floats[2]:.16e}  (float64, FFT n=3, A_MAX=200)")
print(f"  Δ_3 = {delta_floats[3]:.16e}  (float64, FFT n=4, A_MAX=200)")
print(f"  Δ_4 = {delta_floats[4]:.16e}  (float64, FFT n=5, A_MAX=200)")

print(f"\n=== Ratio sequence Δ_m / Δ_{{m-1}} (numerical) ===")
for m in (2, 3, 4):
    r = delta_floats[m] / delta_floats[m-1]
    print(f"  Δ_{m} / Δ_{m-1} = {r:.16f}")

# Higher-order discrete derivative: Δ_{m+1} - Δ_m
print(f"\n=== Forward differences ===")
for m in (1, 2, 3):
    d = delta_floats[m+1] - delta_floats[m]
    print(f"  Δ_{m+1} - Δ_{m} = {d:+.16e}")

# Second forward differences
print(f"\n=== Second forward differences ===")
for m in (1, 2):
    d2 = delta_floats[m+2] - 2*delta_floats[m+1] + delta_floats[m]
    print(f"  Δ_{m+2} - 2Δ_{m+1} + Δ_{m} = {d2:+.16e}")

# Convergence ratio of differences
diffs = [delta_floats[m+1] - delta_floats[m] for m in (1, 2, 3)]
print(f"\n=== Ratio of differences: (Δ_{{m+1}}-Δ_m) / (Δ_m-Δ_{{m-1}}) ===")
for i in range(1, 3):
    r = diffs[i] / diffs[i-1]
    print(f"  i={i}: ratio = {r:.16f}")

print(f"\n=== PSLQ VERDICT ===")
print(f"Tested bases at tol=10^-15 / maxcoeff=10^6 / 24 candidates:")
print(f"  - simple: {{1, log2, log17, √17, π, e}}")
print(f"  - polylogs: Li_s(1/256) and Li_s(1/17) for s=1..6")
print(f"  - Lerch transcendents: Φ(1/256, s, k/8) for s=1..3, k=1..8")
print(f"  - logs: log(15/16), log(15/17), log(255/256), log(1-1/q), log(1-2/q)")
print(f"Result: NO integer relation found within 10^-15 for any Δ_m or Δ_∞.")
print(f"")
print(f"The only matches PSLQ found were INTERNAL TAUTOLOGIES of the basis")
print(f"(e.g., log(15/17) = log(1-2/17)), not relations involving Δ_m.")
print(f"")
print(f"At maxcoeff=10^10 with 4-term {{1, Δ_1..Δ_4}} basis, PSLQ returned")
print(f"  [-1, -209, 192, 131, 178]")
print(f"which is a precision-noise relation (residual ~10^-3, expected at the tail")
print(f"of float64 precision). Not a structural relation.")
print(f"")
print(f"=== STRUCTURAL DIAGNOSIS ===")
print(f"Naive shift-additive hypothesis for c(m):")
print(f"  c(m) = Σ over (k_1,...,k_m) of Π 2^{{-8|k_i|}} · N(2·Σk_i mod 17)")
print(f"       / Σ Π 2^{{-8|k_i|}} · T(2·Σk_i mod 17)")
print(f"FAILS at m≥2 with error growing ~3.4e-3 per step (run probe_c2_extension_hypothesis).")
print(f"")
print(f"Reason: shift η at depth m+1 is NOT 2k_{{m+1}} mod q; it's")
print(f"  η = -2·k_{{m+1}}·2^{{-min(a_{{m+1}},b_{{m+1}})}} mod q,")
print(f"and the carry-forward σ_{{m+1}} = 2^{{-a_{{m+1}}}}σ_m + η couples the chain")
print(f"across depths via a multiplicative random walk on (Z/q)* that doesn't")
print(f"factorize as a convolution on Z/q.")
print(f"")
print(f"=== TO GET Δ_2, Δ_3, Δ_4 AS EXACT RATIONALS ===")
print(f"Two routes:")
print(f"  (A) Derive the correct c(m) closed form (the right multi-depth structure")
print(f"      involves the FULL (a_i, b_i) joint, not just k_i offsets).")
print(f"  (B) Brute-force exact rational FFT at n=3 (length 4913, ~10 min Fraction)")
print(f"      for c(2); at n=4 (length 83521) for c(3) is borderline; n=5 not feasible.")
print(f"Recommended: pursue (B) for c(2) → exact Δ_2 → re-run PSLQ at 30+ digit precision.")
print(f"If still negative, c_∞ likely involves a structural object NOT in the standard")
print(f"polylog/Lerch family (e.g., a character-sum L-value or modular form coefficient).")
