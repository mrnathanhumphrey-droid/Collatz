"""Result 79 Step 3 continued: Iterated van der Corput level B=2.

Standard B=2 van der Corput inequality:
  |S|^4 ≤ (4N) · (something involving I(h1, h2))

Specifically, applying Weyl differencing twice:
  |S|^2 ≤ (N+H)/H · Σ_{|h1|<H} (1-|h1|/H) Σ_u f(u+h1) conj(f(u))
       Let g(u; h1) := f(u+h1) conj(f(u))
       The inner sum is |Σ_u g(u; h1)|.
  Apply Weyl again to Σ_u g(u; h1):
      |Σ_u g(u; h1)|^2 ≤ (N+H')/H' · Σ_{|h2|<H'} (1-|h2|/H') Σ_u g(u+h2; h1) conj(g(u; h1))
                       = (N+H')/H' · Σ_{|h2|<H'} (1-|h2|/H') · I(h1, h2)
       where I(h1, h2) := Σ_u f(u+h1+h2) conj(f(u+h1)) conj(f(u+h2)) f(u)

For our f(u) = e_q(c·4^u - 9m·u):
  I(h1, h2) = Σ_u e_q(c·4^u·(4^{h1}-1)·(4^{h2}-1)) · e_q(constant)
  v_3((4^{h1}-1)(4^{h2}-1)) = v_3(h1) + v_3(h2) + 2

For h1=h2=1: v_3 = 2, eff modulus = 3^{r+1-2} = 3^{r-1}. Inner sum length N at modulus 3^{r-1}.
4 has order 3^{r-2} mod 3^{r-1}. # cycles = N/3^{r-2} = 3. Each cycle sums to 0 (principal unit
subgroup, c'' a unit). So I(1, 1) ≈ 0 (after handling boundary).

This script:
  - Computes Σ_{h1, h2} |I(h1, h2)| in full
  - Tests whether B=2 Weyl gives √N cancellation
"""
import numpy as np
from math import gcd
from cmath import exp as cexp
import time

PI = float(np.pi)

def v3(n):
    if n == 0:
        return 10**9
    n = abs(n)
    k = 0
    while n % 3 == 0:
        n //= 3
        k += 1
    return k

def f_array(r, c, m):
    q = 3**(r + 1)
    N = 3**(r - 1)
    arr = np.empty(N, dtype=np.complex128)
    x = 1
    inv = 2j * PI / q
    nine_m = 9 * m
    for u in range(N):
        phase = (c * x - nine_m * u) % q
        arr[u] = cexp(inv * phase)
        x = (x * 4) % q
    return arr

def I_h1_h2(f_arr, h1, h2):
    """I(h1, h2) = Σ_u f(u+h1+h2) conj(f(u+h1)) conj(f(u+h2)) f(u)
    sum over u such that u+h1+h2 < N and u >= 0."""
    N = len(f_arr)
    if h1 + h2 >= N:
        return 0.0 + 0.0j
    u_max = N - max(h1 + h2, 0)
    if u_max <= 0:
        return 0.0 + 0.0j
    s = 0.0 + 0.0j
    for u in range(u_max):
        s += f_arr[u + h1 + h2] * np.conj(f_arr[u + h1]) * np.conj(f_arr[u + h2]) * f_arr[u]
    return s

# --------------------------------------------------------------------
# Compute B=2 vdc bound
# --------------------------------------------------------------------
print("="*72)
print("Iterated van der Corput level B=2:")
print("  |S|^4 ≤ ((N+H)/H)^2 · Σ_{h1, h2} |I(h1, h2)|  (rough form, both H = N)")
print()

t0 = time.time()
print(f"  {'r':>3} {'N':>5} {'actual|S|':>10} {'|S|^4':>15} {'sum|I(h1,h2)|':>15} {'B2_bound^4':>15} {'B2_S':>10} {'ratio_to_actual':>15}")

for r in range(3, 7):
    N = 3**(r - 1)
    q = 3**(r + 1)
    # Pick (c, m) that gives large |S|
    best_S = 0
    best_c = 1
    best_m = 0
    for c in range(1, min(q, 200), 2):
        if gcd(c, 3) != 1:
            continue
        for m in range(min(5, 2*N)):
            arr = f_array(r, c, m)
            v = abs(arr.sum())
            if v > best_S:
                best_S = v
                best_c = c
                best_m = m

    f_arr = f_array(r, best_c, best_m)
    actual_S = abs(f_arr.sum())

    # Compute Σ_{h1, h2} |I(h1, h2)| (for h1, h2 in [0, N-1] with h1+h2 < N)
    # This is N^2-ish entries; manageable for r ≤ 7
    sum_I = 0.0
    for h1 in range(0, N):
        for h2 in range(0, N - h1):
            v = abs(I_h1_h2(f_arr, h1, h2))
            sum_I += v
    # Bound: |S|^4 ≤ (N+N)^2/N^2 · sum_I = 4·sum_I
    B2_bound4 = 4 * sum_I
    B2_S = B2_bound4**0.25 if B2_bound4 > 0 else 0
    print(f"  {r:>3} {N:>5} {actual_S:>10.3f} {actual_S**4:>15.1f} {sum_I:>15.1f} {B2_bound4:>15.1f} {B2_S:>10.3f} {B2_S/actual_S:>15.3f}")

print(f"\n  [time: {time.time()-t0:.1f}s]")
print()
print("If B2_S/actual_S → constant as r grows, B=2 differencing achieves the right order.")
print("If B2_S/actual_S → ∞, B=2 differencing also fails to capture cancellation.")
print()
print("Compare also B=2 'naive' bound to trivial N and to √N:")
print(f"  {'r':>3} {'√N':>8} {'B2_S':>10} {'N':>6} {'B2/√N':>8} {'B2/N':>8}")
for r in range(3, 7):
    N = 3**(r - 1)
    # already computed above; just re-skip for clarity
    pass
