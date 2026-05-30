"""
probe_wirsching_check_2026_05_30.py

Test whether c_inf ~= 0.152988994 (q=17 deep-collision shift moment) equals the chi-bar
moment of the Wirsching density phi_17 or its self-convolution phi_17 * phi_17^rev
projected to first 17-digit.

Wirsching W_q operator: (W_q f)(x) = (q/(q-1)) * integral of f over [qx-(q-1), qx] cap [0,1].
phi_q is the fixed point with int phi_q = 1, supp [0,1].

First 17-digit projection: P(floor(17 x) = j) = int_{j/17}^{(j+1)/17} phi_17(x) dx.

Compare moments:
  c(0) = 19/127 = 0.149606 (Tao Syracuse self-difference at depth 0)
  c_inf ~= 0.152989 (deep-collision shift)
"""
from __future__ import annotations
import numpy as np
import sys
sys.stdout.reconfigure(encoding="utf-8")

def legendre(x, q):
    x = int(x) % q
    if x == 0: return 0
    return 1 if pow(x, (q - 1) // 2, q) == 1 else -1

def iterate_W(q, N=20001, n_iter=200, verbose=False):
    """Iterate W_q to find phi_q. Returns (phi, x_grid)."""
    x = np.linspace(0, 1, N)
    dx = 1.0 / (N - 1)
    phi = np.ones(N)  # uniform initial
    # cumulative trapezoidal integral
    for it in range(n_iter):
        # F[i] = integral_0^{x[i]} phi
        F = np.zeros(N)
        F[1:] = np.cumsum((phi[:-1] + phi[1:]) / 2 * dx)
        hi = np.clip(q * x, 0, 1)
        lo = np.clip(q * x - (q - 1), 0, 1)
        F_hi = np.interp(hi, x, F)
        F_lo = np.interp(lo, x, F)
        phi_new = (q / (q - 1)) * (F_hi - F_lo)
        norm = np.trapezoid(phi_new, x)
        if norm > 0:
            phi_new /= norm
        # convergence
        diff = np.max(np.abs(phi_new - phi))
        phi = phi_new
        if verbose and (it < 5 or it % 50 == 0):
            print(f"  iter {it}: max diff = {diff:.2e}, norm = {norm:.6f}")
        if diff < 1e-12:
            if verbose: print(f"  converged at iter {it}, diff = {diff:.2e}")
            break
    return phi, x

def digit_dist(phi, x, q):
    """First q-digit distribution: P(floor(qx) = j) for j = 0..q-1."""
    F = np.zeros(len(x))
    dx = x[1] - x[0]
    F[1:] = np.cumsum((phi[:-1] + phi[1:]) / 2 * dx)
    p = np.zeros(q)
    for j in range(q):
        lo, hi = j / q, (j + 1) / q
        p[j] = np.interp(hi, x, F) - np.interp(lo, x, F)
    return p

def chi_bar_moment(p, q):
    """sum_j legendre(j, q) * p[j] / sum_j p[j] over j != 0."""
    chi = np.array([legendre(j, q) for j in range(q)])
    units = chi != 0
    return np.sum(chi[units] * p[units]) / np.sum(p[units])

def self_conv_digit_dist(phi, x, q):
    """Distribution of (U - V) mod q where U, V iid ~ phi on [0,1].
    Returns P((U - V) mod 1 in [j/q, (j+1)/q)) projected to Z/q (digit of mod-1 diff)."""
    # Use FFT to convolve phi with phi_reverse
    # phi_rev(x) = phi(1 - x); convolve phi * phi_rev with circular mod-1
    # Actually we want P(U - V mod 1) which is autocorrelation: phi conv phi_rev
    # Take phi on [0,1], extend by zero, but for mod-1 diff just use circular conv
    # P(D = d) = int phi(d + v) phi(v) dv for d in [0,1) (mod 1)
    # = (phi * phi_reflected)(d) circularly on torus
    F = np.fft.fft(phi)
    PD = np.real(np.fft.ifft(F * np.conj(F))) * (x[1] - x[0])  # autocorrelation
    # PD lives on x grid, normalize
    norm = np.trapezoid(PD, x)
    if norm > 0: PD /= norm
    # First q-digit of (U-V mod 1)
    F_PD = np.zeros(len(x))
    dx = x[1] - x[0]
    F_PD[1:] = np.cumsum((PD[:-1] + PD[1:]) / 2 * dx)
    p = np.zeros(q)
    for j in range(q):
        lo, hi = j / q, (j + 1) / q
        p[j] = np.interp(hi, x, F_PD) - np.interp(lo, x, F_PD)
    return p

for q in (3, 17):
    print(f"\n{'='*70}")
    print(f"q = {q}")
    print(f"{'='*70}")
    phi, x = iterate_W(q, N=20001, n_iter=300, verbose=(q == 17))
    print(f"  phi_{q} support / shape: min={phi.min():.4f}, max={phi.max():.4f}, "
          f"phi(1/2)={np.interp(0.5, x, phi):.4f}, phi(1/q)={np.interp(1/q, x, phi):.4f}")
    p_digit = digit_dist(phi, x, q)
    print(f"  first-{q}-digit distribution of phi_{q}:")
    for j in range(q):
        print(f"    P(digit={j}) = {p_digit[j]:.8f}  (Legendre={legendre(j,q):+d})")
    if q > 2:
        chi_phi = chi_bar_moment(p_digit, q)
        print(f"  Legendre chi-bar moment of phi_{q} first digit = {chi_phi:.10f}")
    # Self-conv test
    p_conv = self_conv_digit_dist(phi, x, q)
    if q > 2:
        chi_conv = chi_bar_moment(p_conv, q)
        print(f"  Legendre chi-bar moment of (phi_{q} * phi_{q}^rev) first digit = {chi_conv:.10f}")

print(f"\n{'='*70}")
print("Comparison to known q=17 values:")
print(f"{'='*70}")
print(f"  c(0)  = 19/127           = {19/127:.10f}  (Tao Syracuse self-diff depth 0)")
print(f"  c_inf ~= 0.152988994400  (deep-collision shift moment)")
print()
print("If phi_17 first-digit chi-bar moment matches c_inf or c(0), there is a structural")
print("identification. If neither matches, Wirsching's framework does not crack pi.")
