"""
nisoli_closure_corrected.py — verify Nisoli closure inequality at T_lead's
corrected rate lambda = 43/45 for c = 7/45.

Phase 1: verify M_3'' = ||(I - T_lead)^{-1}|| (operator l^2 norm).
Phase 2-3: tabulate closure inequality |K| * K^{-A} * M_3'' < 1
           over (r, A, K) grids.

Outputs:
- M_3'' exact via Fractions + SVD via float.
- A table of |K|*K^{-A}*M_3'' < 1 over candidate (r, A, K).
"""

from fractions import Fraction
import math

# ---- Phase 1: T_lead and M_3'' ----

# T_lead = (1/45) * [[7, 9], [28, 36]]
T_lead = [[Fraction(7, 45), Fraction(9, 45)],
          [Fraction(28, 45), Fraction(36, 45)]]

# I - T_lead = (1/45) * [[38, -9], [-28, 9]]
I_minus_T = [[Fraction(1, 1) - T_lead[0][0], -T_lead[0][1]],
             [-T_lead[1][0], Fraction(1, 1) - T_lead[1][1]]]

print("I - T_lead =")
for row in I_minus_T:
    print(" ", [str(x) for x in row])

# det(I - T_lead)
a, b = I_minus_T[0]
c, d = I_minus_T[1]
det = a * d - b * c
print(f"\ndet(I - T_lead) = {det} = {float(det)}")

# (I - T_lead)^{-1} = (1/det) * [[d, -b], [-c, a]]
inv = [[d / det, -b / det],
       [-c / det, a / det]]
print("\n(I - T_lead)^{-1} =")
for row in inv:
    print(" ", [str(x) for x in row])
    print("    ", [float(x) for x in row])

# Operator l^2 norm: sqrt of max eigenvalue of A^T A
# A = inv
A = inv
ATA = [[A[0][0]**2 + A[1][0]**2, A[0][0]*A[0][1] + A[1][0]*A[1][1]],
       [A[0][0]*A[0][1] + A[1][0]*A[1][1], A[0][1]**2 + A[1][1]**2]]

print("\nA^T A (exact rationals):")
for row in ATA:
    print(" ", [str(x) for x in row])

trace = ATA[0][0] + ATA[1][1]
det_ATA = ATA[0][0] * ATA[1][1] - ATA[0][1] * ATA[1][0]
print(f"\ntrace(A^T A) = {trace} = {float(trace)}")
print(f"det(A^T A) = {det_ATA} = {float(det_ATA)}")

# Eigenvalues of 2x2: (trace +- sqrt(trace^2 - 4*det)) / 2
disc = float(trace)**2 - 4*float(det_ATA)
sigma2_max = (float(trace) + math.sqrt(disc)) / 2
sigma2_min = (float(trace) - math.sqrt(disc)) / 2

sigma_max = math.sqrt(sigma2_max)
sigma_min = math.sqrt(sigma2_min)

print(f"\nsigma^2_max = {sigma2_max}")
print(f"sigma^2_min = {sigma2_min}")
print(f"\nsigma_max = ||(I - T_lead)^{{-1}}||_2 = {sigma_max}")
print(f"sigma_min = {sigma_min}")
print(f"\nCondition number kappa = sigma_max/sigma_min = {sigma_max/sigma_min}")

# Also check spectral radius
# Eigenvalues of (I - T_lead)^{-1}: 1/(1-0) = 1 and 1/(1 - 43/45) = 45/2
print(f"\nSpectral radius of (I - T_lead)^{{-1}} = 45/2 = {45/2}")
print(f"  (l^2 op norm exceeds spectral radius due to non-orthogonal eigenbasis)")

M3_pp = sigma_max
print(f"\n=== M_3'' = ||(I - T_lead)^{{-1}}||_op = {M3_pp:.6f} ===")

# ---- Phase 2-3: Closure inequality tabulation ----

print("\n" + "="*70)
print("Phase 2-3: Closure inequality |K_bil| * K^{-A} * M_3'' < 1")
print("="*70)

p = 3  # q = p = 3 throughout for c = 7/45
sqrt_p = math.sqrt(p)

def bilinear_bound(r):
    """|K| at level r. r <= 3: 2*sqrt(N) strict. r >= 4: 2*sqrt(p)*sqrt(N)."""
    N = p**(r - 1)
    sqrt_N = math.sqrt(N)
    if r <= 3:
        return 2 * sqrt_N, "strict 2*sqrt(N)"
    else:
        return 2 * sqrt_p * sqrt_N, "polylog-free 2*sqrt(p)*sqrt(N)"

# Tabulate
print(f"\nM_3'' = {M3_pp:.4f}")
print(f"\nLegend: |K_bil| = bilinear bound on |S_partial| at level r")
print(f"        K = truncation level (separate from r; Nisoli's K param)")
print(f"        Closure fires iff |K_bil| * K^{{-A}} * M_3'' < 1")
print()

header = f"{'r':>3} {'N':>10} {'sqrt(N)':>10} {'|K_bil|':>12} {'A':>3} {'K':>5} {'K^-A':>14} {'|Kbil|*K^-A*M3':>18} {'fires':>8}"
print(header)
print("-" * len(header))

firing_cells = []

r_values = [2, 3, 4, 5, 6, 8, 10]
A_values = [1, 2, 3, 5, 10, 20]
K_values = [6, 10, 20, 50, 100, 500, 1000]

for r in r_values:
    N = p**(r - 1)
    sqrt_N = math.sqrt(N)
    K_bil, note = bilinear_bound(r)
    for A in A_values:
        for K in K_values:
            KmA = K**(-A)
            product = K_bil * KmA * M3_pp
            fires = product < 1
            mark = "YES" if fires else "no"
            if fires:
                firing_cells.append((r, A, K, K_bil, KmA, product))
            print(f"{r:>3} {N:>10} {sqrt_N:>10.3f} {K_bil:>12.3f} {A:>3} {K:>5} {KmA:>14.4e} {product:>18.6e} {mark:>8}")

# Summary
print("\n" + "="*70)
print(f"Firing cells (closure inequality satisfied): {len(firing_cells)}")
print("="*70)

if firing_cells:
    print(f"\n{'r':>3} {'A':>3} {'K':>5} {'|K_bil|':>12} {'K^-A':>14} {'product':>14}")
    print("-" * 60)
    for (r, A, K, K_bil, KmA, product) in firing_cells[:30]:
        print(f"{r:>3} {A:>3} {K:>5} {K_bil:>12.3f} {KmA:>14.4e} {product:>14.6e}")
    print(f"\n[showing up to 30 of {len(firing_cells)} firing cells]")

# Smallest firing A at each (r, K)
print("\n" + "="*70)
print("Minimum A for closure at each (r, K)")
print("="*70)

print(f"\n{'r':>3} {'K':>5} {'|K_bil|':>12} {'min A':>8} {'K^-A':>14} {'product':>14}")
print("-" * 60)

for r in r_values:
    K_bil, _ = bilinear_bound(r)
    for K in [6, 10, 20, 50, 100, 1000]:
        # Solve K_bil * K^{-A} * M3 < 1  =>  K^A > K_bil * M3  =>  A > log(K_bil*M3)/log(K)
        thr = K_bil * M3_pp
        if K > 1:
            A_min = math.log(thr) / math.log(K)
            A_int = math.ceil(A_min) if A_min > 0 else 1
            prod = K_bil * K**(-A_int) * M3_pp
            print(f"{r:>3} {K:>5} {K_bil:>12.3f} {A_min:>8.3f} {K**(-A_int):>14.4e} {prod:>14.6e}")

# Specific scenario: r=3 (strict 2*sqrt(N)), K matched to truncation level
print("\n" + "="*70)
print("FOCUS: r=3 (strict 2*sqrt(N) bilinear bound, N=9, |K_bil|=6)")
print("="*70)
r = 3
K_bil, _ = bilinear_bound(r)
print(f"  |K_bil| = {K_bil}")
print(f"  |K_bil| * M_3'' = {K_bil * M3_pp:.4f}  (required K^A > this)")
print()
print(f"{'A':>3} {'K_min for closure':>20}")
print("-" * 30)
thr = K_bil * M3_pp
for A in [1, 2, 3, 5, 10]:
    # K^A > thr => K > thr^(1/A)
    K_min = thr**(1/A)
    print(f"{A:>3} {K_min:>20.4f}")

# Also: r=4 polylog-free
print("\n" + "="*70)
print("FOCUS: r=4 (polylog-free 2*sqrt(p)*sqrt(N) bilinear bound, N=27, |K_bil|=2*sqrt(3)*sqrt(27))")
print("="*70)
r = 4
K_bil, _ = bilinear_bound(r)
print(f"  |K_bil| = {K_bil:.4f}")
print(f"  |K_bil| * M_3'' = {K_bil * M3_pp:.4f}  (required K^A > this)")
print()
print(f"{'A':>3} {'K_min for closure':>20}")
print("-" * 30)
thr = K_bil * M3_pp
for A in [1, 2, 3, 5, 10]:
    K_min = thr**(1/A)
    print(f"{A:>3} {K_min:>20.4f}")

print("\n" + "="*70)
print("Summary done.")
print("="*70)
