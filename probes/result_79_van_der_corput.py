"""Result 79 Step 2-3: van der Corput differencing — does it gain a factor?

Tests:
  (1) Empirical: |S_{r,c,m}|/√N at r = 2..8 (confirm √N cancellation extends)
  (2) For each (r, c, m=0), compute FULL Σ_{h=1}^{N-1} |I(h)| (van der Corput input)
  (3) Compute Weyl bound: |S|² ≤ N + 2·Σ_h |I(h)| (taking H = N for simplicity)
  (4) Compare to direct |S|², see if Weyl is sharp or loose
  (5) Try iterated differencing (van der Corput level B = 2): |S|^4 ≤ N²·H + N·Σ_{h1,h2} |I(h1,h2)|
"""
import numpy as np
from math import gcd
from cmath import exp as cexp
from cmath import pi
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

def kalafatelis_S(r, c, m, length=None):
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

def kalafatelis_S_array(r, c, m):
    """Return array of f(u) = e_q(c·4^u - 9m·u) for u=0..N-1."""
    q = 3**(r + 1)
    N = 3**(r - 1)
    arr = np.empty(N, dtype=np.complex128)
    x = 1
    nine_m = 9 * m
    inv = 2j * PI / q
    for u in range(N):
        phase = (c * x - nine_m * u) % q
        arr[u] = cexp(inv * phase)
        x = (x * 4) % q
    return arr

def autocorr(f_arr):
    """I(h) = Σ_u f(u+h)·conj(f(u)) for h = 0..N-1.
    Computed via FFT-style auto-correlation. Returns array of length N.
    Note: standard def, h-th lag with full overlap (length N - |h|)."""
    N = len(f_arr)
    # Compute correlation via FFT for efficiency
    F = np.fft.fft(f_arr, 2*N)
    # autocorrelation = ifft(|F|^2), but need only positive lags
    R = np.fft.ifft(F * F.conj())
    # R[h] = sum_{u} f(u+h) conj(f(u)) for h = 0..N-1, with circular wrap
    # We want the LINEAR autocorrelation: sum_{u: u+h<N, u>=0} f(u+h) conj(f(u))
    # That's: R[h] for h in 0..N-1 (taking the first N entries from the 2N FFT).
    return R[:N]

def vdc_weyl_bound(I_arr, H):
    """Standard van der Corput bound:
       |S|² ≤ (N+H-1)/H · Σ_{|h|<H} (1 - |h|/H) Re(I(h))
    Returns (bound², bound_S)."""
    N = len(I_arr)
    # Σ_{h=-(H-1)}^{H-1} (1-|h|/H) Re I(h)
    # = (1) Re I(0) + 2·Σ_{h=1}^{H-1} (1-h/H) Re I(h)
    if H > N:
        H = N
    s = float(I_arr[0].real)
    for h in range(1, H):
        s += 2 * (1 - h/H) * float(I_arr[h].real)
    s *= (N + H - 1) / H
    if s < 0:
        s = 0  # numerical
    return s, s**0.5

def vdc_naive_l1(I_arr, H):
    """Loose bound: Σ_{|h|<H} |I(h)|, then |S|² ≤ (N+H-1)/H · this."""
    N = len(I_arr)
    if H > N:
        H = N
    s = float(abs(I_arr[0]))
    for h in range(1, H):
        s += 2 * (1 - h/H) * float(abs(I_arr[h]))
    s *= (N + H - 1) / H
    return s, s**0.5

# --------------------------------------------------------------------
# (1) Direct sizes at r = 2..8
# --------------------------------------------------------------------
print("="*72)
print("(1) Direct |S_{r,c,m}| sizes — square-root scaling check")
print()
print("    r | N | √N | max|S| (over c-units, m=0..min(9,N_r-1)) | max/√N")

t0 = time.time()
for r in range(2, 9):
    N = 3**(r - 1)
    q = 3**(r + 1)
    sqrtN = N**0.5
    max_S = 0.0
    arg = None
    # Sample c uniformly: too many c's at high r, sample 200
    cs = list(range(1, q)) if q < 1000 else [c for c in np.linspace(1, q-1, 200, dtype=int) if gcd(int(c), 3) == 1]
    if q < 1000:
        cs = [c for c in cs if gcd(c, 3) == 1]
    for c in cs[:300]:  # cap at 300 c-values
        for m in range(min(10, 2*N)):
            s = kalafatelis_S(r, int(c), m)
            v = abs(s)
            if v > max_S:
                max_S = v
                arg = (int(c), m)
    print(f"    {r} | {N:5d} | {sqrtN:>7.2f} | {max_S:>8.3f} (c={arg[0]:>5d}, m={arg[1]}) | {max_S/sqrtN:.3f}")

print(f"  [time: {time.time()-t0:.1f}s]")
print()

# --------------------------------------------------------------------
# (2-3) Compute autocorrelation and Weyl bound
# --------------------------------------------------------------------
print("="*72)
print("(2-3) Autocorrelation Σ_h |I(h)| and van der Corput bound")
print("       Standard form: |S|² ≤ (1) Re I(0) + 2 Σ (1-h/H) Re I(h), with H = N")
print("       Loose form (l1):     replace Re I(h) with |I(h)|")
print()
print("    r | N | √N | actual |S| | sum |I(h)| | √N·sumI/N | Weyl-real | Weyl-l1")

for r in range(3, 9):
    N = 3**(r - 1)
    q = 3**(r + 1)
    sqrtN = N**0.5
    # Pick a (c, m) that gave large |S| from panel (1)
    # For simplicity, scan a few c's
    best_S = 0
    best_c = 1
    best_m = 0
    cs_to_try = list(range(1, min(q, 1000)))
    cs_to_try = [c for c in cs_to_try if gcd(c, 3) == 1][:50]
    for c in cs_to_try:
        for m in range(min(5, 2*N)):
            s = kalafatelis_S(r, c, m)
            if abs(s) > best_S:
                best_S = abs(s)
                best_c = c
                best_m = m

    f_arr = kalafatelis_S_array(r, best_c, best_m)
    actual_S = abs(f_arr.sum())

    I_arr = autocorr(f_arr)
    sum_abs_I = float(np.abs(I_arr).sum())  # raw sum of |I(h)|

    # Weyl with H = N
    weyl_real_bound2, weyl_real_S = vdc_weyl_bound(I_arr, N)
    weyl_l1_bound2, weyl_l1_S = vdc_naive_l1(I_arr, N)

    # heuristic: if sum |I(h)| is O(N), Weyl gives √N cancellation
    # Specifically, Weyl_l1_S ≤ √(2·N·B) where B = max |I(h)|. Let's see.
    print(f"    {r} | {N:5d} | {sqrtN:>6.2f} | {actual_S:>8.3f}    | {sum_abs_I:>8.1f}    | {sum_abs_I/N:.2f}      | {weyl_real_S:>8.2f}  | {weyl_l1_S:>8.2f}")

print()

# --------------------------------------------------------------------
# (4) Per-h breakdown: distribution of |I(h)| by v_3(h)
# --------------------------------------------------------------------
print("="*72)
print("(4) Per-h |I(h)| distribution by k = v_3(h)")
print("    Theory: I(h) = e^{phase} · partial_cycle of length ρ_k(h)")
print("    Expectation: |I(h)|² has mean ≈ ρ_k (random-phase model)")
print()
for r in [5, 6, 7]:
    N = 3**(r - 1)
    q = 3**(r + 1)
    f_arr = kalafatelis_S_array(r, 1, 0)
    I_arr = autocorr(f_arr)
    print(f"  r = {r}, N = {N}")
    print(f"    {'k':>3} {'#h':>5} {'<rho>':>8} {'<|I|>':>10} {'max|I|':>10} {'sum|I|':>10}")
    by_k = {}
    for h in range(1, N):
        k = v3(h)
        period = 3**(r - k - 1) if r - k - 1 >= 0 else 1
        rho = (N - h) % period if period > 0 else 0
        absI = abs(I_arr[h])
        by_k.setdefault(k, []).append((rho, absI))
    for k in sorted(by_k.keys()):
        items = by_k[k]
        n = len(items)
        avg_rho = np.mean([x[0] for x in items])
        avg_I = np.mean([x[1] for x in items])
        max_I = max(x[1] for x in items)
        sum_I = sum(x[1] for x in items)
        print(f"    {k:>3} {n:>5d} {avg_rho:>8.1f} {avg_I:>10.3f} {max_I:>10.3f} {sum_I:>10.2f}")
    print()

# --------------------------------------------------------------------
# (5) Iterated van der Corput — level B = 2 differencing
# --------------------------------------------------------------------
print("="*72)
print("(5) Iterated differencing (level B=2)")
print("    |S|^4 ≤ N²·H + N²/H · Σ_{h1,h2} |I(h1,h2)|  (rough form)")
print("    where I(h1,h2) = Σ_u f(u+h1+h2) conj(f(u+h1)) conj(f(u+h2)) f(u)")
print("    For our f(u) = e_q(c·4^u - 9m·u):")
print("      I(h1,h2) = e_q(... constant ...) · Σ_u e_q(c·(4^{h1}-1)·(4^{h2}-1)·4^u)")
print("    LTE: v_3((4^{h1}-1)(4^{h2}-1)) = (k1+1)+(k2+1) = k1+k2+2")
print("    So inner sum has effective modulus 3^{r+1-k1-k2-2} = 3^{r-k1-k2-1}.")
print()

def iter_vdc_inner(r, c, m, h1, h2):
    """Return Σ_{u=0}^{N-h1-h2-1} f(u+h1+h2)·conj(f(u+h1))·conj(f(u+h2))·f(u)
    where f(u) = e_q(c·4^u - 9m·u)."""
    N = 3**(r - 1)
    q = 3**(r + 1)
    # f(u+a) f(u+b) form: phase(u+a) - phase(u+b) - phase(u+c) + phase(u+d) etc.
    # Total: phase(u+h1+h2) - phase(u+h1) - phase(u+h2) + phase(u)
    # = c·(4^{u+h1+h2} - 4^{u+h1} - 4^{u+h2} + 4^u) - 9m(0)  [linear cancels!]
    # = c·4^u·(4^{h1+h2} - 4^{h1} - 4^{h2} + 1)
    # = c·4^u·(4^{h1} - 1)·(4^{h2} - 1)
    diff = (pow(4, h1, q) - 1) * (pow(4, h2, q) - 1) % q
    s = 0.0 + 0.0j
    x = 1
    inv = 2j * PI / q
    for u in range(N - h1 - h2):
        phase = (c * x * diff) % q
        s += cexp(inv * phase)
        x = (x * 4) % q
    return s, diff

print("  Sample (r, h1, h2, k1=v3(h1), k2=v3(h2), v3(diff), |I(h1,h2)|, partial_cycle bound)")
for r in [5, 6]:
    N = 3**(r - 1)
    print(f"  r = {r}, N = {N}")
    print(f"    {'h1':>4} {'h2':>4} {'k1+k2':>5} {'v3':>4} {'eff_mod':>10} {'period':>8} {'rho':>6} {'|I|':>10}")
    for h1, h2 in [(1, 1), (1, 3), (3, 3), (3, 9), (9, 9), (1, 9), (2, 5), (3, 6), (1, 2)]:
        if h1 + h2 >= N:
            continue
        s, diff = iter_vdc_inner(r, 1, 0, h1, h2)
        k_sum = v3(h1) + v3(h2)
        v_diff = v3(diff) if diff != 0 else -1
        eff_mod_exp = (r + 1) - v_diff if v_diff > 0 else r + 1
        eff_period = 3**(eff_mod_exp - 1) if eff_mod_exp >= 1 else 1
        rho = (N - h1 - h2) % eff_period if eff_period > 0 else 0
        print(f"    {h1:>4} {h2:>4} {k_sum:>5d} {v_diff:>4d} {3**eff_mod_exp:>10d} {eff_period:>8d} {rho:>6d} {abs(s):>10.4f}")
    print()
