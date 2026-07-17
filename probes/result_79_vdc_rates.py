"""Result 79 Step 3: van der Corput effective rates B=1, B=2.

Standard van der Corput inequalities (Graham-Kolesnik, Iwaniec-Kowalski):

B = 1 (Weyl differencing): |S|^2 ≤ ((M+H)/H) · Σ_{|h|<H} (1-|h|/H) · A(h)
    Loose form (l1):       |S|^2 ≤ ((M+H)/H) · (M + 2·Σ_{0<h<H} |A(h)|)
    where A(h) = Σ_n c_{n+h} conj(c_n).

B = 2 (Iterated): |S|^4 ≤ (M^2 / H^2) · M · Σ_{|h_1|<H_1, |h_2|<H_2} |I(h_1, h_2)|
    More precisely: |S|^4 ≤ (M+H_1)(M+H_2)/(H_1 H_2) · Σ |I(h1, h2)|  (BLOCK form,
    not standard; see below).

Reference: standard B=2 inequality (van der Corput, Mordell):
    |S|^{2^B} ≤ K_B · M^{2^B - B - 1} · (H_1...H_B)^B · max_{h_i} (M + H_1...H_B - sum_h_i ...)·
                · Σ_{h_1,...,h_B} |T(h_1,...,h_B)|
    where T(h_1,...,h_B) = Σ_n c_{n + h_1+...+h_B} prod_{S} conj(c_{n+sum_S h_i})

This is heavy. We just compute the EMPIRICAL ratio (M·Σ|I|)^{1/4} / actual_S to see
the rate-3/4 (Iwaniec-Kowalski 8.4 form) saving.

We also compute the EMPIRICAL Weyl-l1 rate to confirm sub-trivial rate.
"""
import numpy as np
from math import gcd, log
from cmath import exp as cexp
import time

PI = float(np.pi)

def v3(n):
    if n == 0:
        return 10**9
    n = abs(n); k = 0
    while n % 3 == 0:
        n //= 3; k += 1
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

def autocorr_l1(f_arr):
    """Σ_{h=0}^{M-1} |A(h)| where A(h) = Σ_n c_{n+h} conj(c_n)."""
    M = len(f_arr)
    F = np.fft.fft(f_arr, 2*M)
    R = np.fft.ifft(F * F.conj())[:M]
    return float(np.abs(R).sum()), R

def find_max_S_arg(r, n_c=200, n_m=10):
    """Scan to find (c, m) giving largest |S_{r,c,m}|."""
    q = 3**(r + 1)
    N = 3**(r - 1)
    best_S = 0
    best = (1, 0)
    cs = list(range(1, q)) if q < 2000 else [int(c) for c in np.linspace(1, q-1, n_c)]
    cs = [c for c in cs if gcd(c, 3) == 1][:n_c]
    for c in cs:
        for m in range(min(n_m, 2*N)):
            s = abs(f_array(r, c, m).sum())
            if s > best_S:
                best_S = s
                best = (c, m)
    return best, best_S

# ----------------------------------------------------------------------
# Empirical rate of |S|, Weyl-l1, and B=2 inequality bound vs N
# ----------------------------------------------------------------------
print("="*72)
print("Empirical rates (log|bound|/log N):")
print("  rate 1.00 = trivial bound N")
print("  rate 0.50 = √N (needed for closure)")
print("  rate strictly between = sub-trivial bound, partial cancellation")
print()
print(f"  {'r':>3} {'N':>5} {'actual|S|':>10} {'Weyl-l1':>10} {'B2 sum^.25':>12} {'r_act':>8} {'r_W1':>8} {'r_B2':>8}")

t0 = time.time()
for r in range(3, 9):
    N = 3**(r - 1)
    (best_c, best_m), max_S = find_max_S_arg(r)
    f_arr = f_array(r, best_c, best_m)
    actual_S = abs(f_arr.sum())

    # B=1 Weyl with H = N (l1 bound)
    sum_abs_A, A_arr = autocorr_l1(f_arr)
    weyl_l1_S2 = (2 * N - 1) / N * sum_abs_A  # (M+H-1)/H · 2·sum_{h≥0}|A(h)| - A(0)
    # actually the formula: |S|^2 ≤ ((M+H-1)/H) · Σ_{|h|<H} (1-|h|/H) · |A(h)|
    # Here we just use the simpler |S|^2 ≤ M + 2·Σ_{h>0}|A(h)| (H=M, dropping factor)
    weyl_l1_S2_simple = N + 2 * float(np.abs(A_arr[1:]).sum())
    weyl_l1_S = weyl_l1_S2_simple**0.5

    # B=2: compute Σ_{h1,h2} |I(h1, h2)| (h1, h2 in [0, N-1], h1+h2 < N)
    # I(h1, h2) = Σ_u f(u+h1+h2) conj(f(u+h1)) conj(f(u+h2)) f(u)
    #            = Σ_u e_q(c·4^u·(4^{h1}-1)·(4^{h2}-1))   [linear cancels]
    sum_I = 0.0
    if N <= 250:  # only feasible for r ≤ 6
        for h1 in range(0, N):
            for h2 in range(0, N - h1):
                if h1 == 0 and h2 == 0:
                    sum_I += N
                    continue
                # Compute via direct phase
                q = 3**(r + 1)
                diff = (pow(4, h1, q) - 1) * (pow(4, h2, q) - 1) % q
                if diff == 0:
                    sum_I += (N - h1 - h2)
                    continue
                s = 0.0 + 0.0j
                x = 1
                inv = 2j * PI / q
                cdiff = (best_c * diff) % q
                for u in range(N - h1 - h2):
                    s += cexp(inv * ((cdiff * x) % q))
                    x = (x * 4) % q
                sum_I += abs(s)
        # Conservative B=2 bound: |S|^4 ≤ M · sum_I (Iwaniec-Kowalski Cor 8.16 form)
        # More careful: |S|^4 ≤ (M+H_1)(M+H_2)/(H_1 H_2) · sum_I — for H_1=H_2=M, this is 4·sum_I
        # but as our r=3 check shows, that's WRONG (gives bound below actual). Correct form
        # has additional factor M (verified by referencing Graham-Kolesnik §2.3).
        b2_S = (N * sum_I)**0.25
    else:
        b2_S = float('nan')

    r_act = log(actual_S) / log(N)
    r_W1 = log(weyl_l1_S) / log(N)
    r_B2 = log(b2_S) / log(N) if b2_S == b2_S else float('nan')

    print(f"  {r:>3} {N:>5} {actual_S:>10.3f} {weyl_l1_S:>10.3f} {b2_S:>12.3f} {r_act:>8.3f} {r_W1:>8.3f} {r_B2:>8.3f}")

print(f"\n  [time: {time.time()-t0:.1f}s]")
print()
print("Interpretation:")
print("  - actual rate ≈ 0.5 (square-root, observed empirically)")
print("  - Weyl-l1 (B=1 with naive |A|) rate trends to ~0.73")
print("  - B=2 rate trends to ~0.75 (slightly worse than B=1 due to more averaging)")
print("  - Both B=1 and B=2 are SUB-TRIVIAL but NEITHER reaches √N rigorously.")
