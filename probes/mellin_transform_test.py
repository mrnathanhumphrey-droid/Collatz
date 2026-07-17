"""
mellin_transform_test.py — Mellin transform of σ-distribution and σ-band-conditional
T-distribution.

Mellin transform of a discrete distribution on positive integers:
  M[P](s) = Σ_x P(x) · x^(s-1) = E[x^(s-1)]

For σ-distribution P(σ|r) per residue r mod 32:
  M(s, r) = Σ_σ P(σ|r) · σ^(s-1)

Probe at:
  - Real s ∈ [-2, 4] (moments)
  - Imaginary s = it for t ∈ [0, 30] (oscillatory)
  - Complex s = a + it grid (look for zeros)

Test:
  (1) Pole structure (none expected for bounded exp-tail distribution)
  (2) Zero structure on imaginary line / critical line s = 1/2 + it
  (3) Functional equations M(s) = f · M(c-s)
  (4) Asymptotic expansion at large |s|
  (5) Universality across residues / σ-bands

5 seeds × 1M orbits at N=2^32. Compute: ~30s.
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from scipy import stats, special

sys.stdout.reconfigure(encoding="utf-8")

LOG2 = np.log(2.0)
LOG3 = np.log(3.0)
LOG43 = 2 * LOG2 - LOG3
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk(starts, max_steps=400000):
    n = len(starts)
    sigma_arr = np.zeros(n, dtype=np.int32)
    j_arr = np.full(n, -1, dtype=np.int8)
    for i in prange(n):
        m = starts[i]
        steps = 0
        j_attr = -1
        while m != 1 and steps < max_steps:
            if m > MAX_VAL // 3:
                break
            three_m = 3 * m + 1
            v = 0
            while (three_m & 1) == 0:
                three_m >>= 1
                v += 1
            steps += 1
            if three_m == 1:
                if v % 2 == 0:
                    j_attr = v // 2
            m = three_m
        sigma_arr[i] = steps
        j_arr[i] = j_attr
    return sigma_arr, j_arr


def mellin_discrete(P, x_values, s):
    """M(s) = Σ_x P(x) · x^(s-1) for complex s.
       x_values is array of σ values (positive integers)."""
    x = x_values.astype(np.complex128)
    return np.sum(P.astype(np.complex128) * x ** (s - 1))


def mellin_grid(P, x_values, s_grid):
    """Vectorized Mellin transform over a grid of complex s values."""
    out = np.zeros(len(s_grid), dtype=np.complex128)
    for i, s in enumerate(s_grid):
        out[i] = mellin_discrete(P, x_values, s)
    return out


def find_zeros_on_line(P, x_values, t_grid, line_re=0.5):
    """Search for zeros of M(line_re + it) along t-grid."""
    s_grid = line_re + 1j * t_grid
    M_vals = np.array([mellin_discrete(P, x_values, s) for s in s_grid])
    abs_M = np.abs(M_vals)
    log_abs = np.log(abs_M + 1e-30)
    # Find local minima
    minima = []
    for i in range(1, len(log_abs) - 1):
        if log_abs[i] < log_abs[i - 1] and log_abs[i] < log_abs[i + 1]:
            if log_abs[i] < log_abs.mean() - 2:  # significant dip
                minima.append((t_grid[i], abs_M[i], log_abs[i]))
    return minima, M_vals, abs_M


def main():
    log2N = 32
    N = 1 << log2N
    seeds = [42, 137, 271, 314, 1729]
    n_starts = 1_000_000

    print(f"# Mellin transform test for σ-distribution at N=2^{log2N}")
    print(f"# {len(seeds)} seeds × {n_starts:,} orbits\n")

    # Walk and aggregate per residue mod 32
    n_residues = 16  # residues 1, 3, ..., 31
    sigma_per_r = {r: [] for r in range(1, 32, 2)}
    sigma_per_band = {q: [] for q in [0.125, 0.375, 0.625, 0.875, 0.975]}
    all_sigma = []
    t0_total = time.perf_counter()
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1) // 2, size=n_starts, dtype=np.int64) + 1
        t0 = time.perf_counter()
        sigma, j = walk(starts)
        elapsed = time.perf_counter() - t0
        valid = sigma > 0
        sigma = sigma[valid]
        starts = starts[valid]
        all_sigma.append(sigma)
        # Per residue mod 32
        for r in range(1, 32, 2):
            mask = (starts & 31) == r
            sigma_per_r[r].append(sigma[mask])
        # Per σ-band
        edges = np.quantile(sigma, [0.25, 0.50, 0.75, 0.95])
        band_idx = np.digitize(sigma, edges)
        for b_idx, q in enumerate([0.125, 0.375, 0.625, 0.875, 0.975]):
            mask = band_idx == b_idx
            sigma_per_band[q].append(sigma[mask])
        print(f"  seed={seed:>5} ({elapsed:.1f}s)  valid={valid.sum():,}")
    print(f"\n# Total: {time.perf_counter()-t0_total:.1f}s\n")

    # Concatenate
    all_sigma = np.concatenate(all_sigma)
    for r in sigma_per_r:
        sigma_per_r[r] = np.concatenate(sigma_per_r[r])
    for q in sigma_per_band:
        sigma_per_band[q] = np.concatenate(sigma_per_band[q])

    # Build empirical P(σ) global
    sigma_max = max(int(all_sigma.max()), 200)
    sigma_values = np.arange(1, sigma_max + 1, dtype=np.float64)
    counts_global = np.bincount(all_sigma, minlength=sigma_max + 1)[1:]  # skip 0
    P_global = counts_global / counts_global.sum()

    # === Step 1: Mellin of global σ-distribution at real s ===
    print(f"# Mellin M(s) at real s (moments of σ):")
    print(f"#   {'s':>6}  {'M(s)_real':>12}  {'interpretation':>30}")
    for s_real in [-1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]:
        M = mellin_discrete(P_global, sigma_values, s_real)
        if s_real == 1:
            interp = "Σ P(σ) = 1 (sanity)"
        elif s_real == 2:
            interp = f"⟨σ⟩ = {M.real:.2f}"
        elif s_real == 3:
            interp = f"⟨σ²⟩ = {M.real:.2f}"
        else:
            interp = f"E[σ^({s_real-1})]"
        print(f"#   {s_real:>6.2f}  {M.real:>12.6f}  {interp:>30}")

    # === Step 2: Imaginary axis - look for oscillation / zeros ===
    print(f"\n# Mellin M(it) on imaginary axis (oscillatory structure):")
    t_grid = np.linspace(0, 30, 301)
    s_grid_imag = 1j * t_grid
    M_imag = np.array([mellin_discrete(P_global, sigma_values, s) for s in s_grid_imag])
    abs_M_imag = np.abs(M_imag)
    print(f"#   |M(it)| at sample points:")
    for t in [0, 1, 2, 5, 10, 20, 30]:
        idx = int(t / 30 * 300)
        print(f"#     t={t:>4}: |M(it)| = {abs_M_imag[idx]:.4e}, arg = {np.angle(M_imag[idx]):.4f}")
    # Local minima search
    print(f"#   Local minima of log|M(it)| (potential zeros):")
    log_abs = np.log(abs_M_imag + 1e-30)
    found_min = 0
    for i in range(2, len(log_abs) - 2):
        if (log_abs[i] < log_abs[i - 1] and log_abs[i] < log_abs[i + 1] and
            log_abs[i] < log_abs.mean() - 1):
            print(f"#     t={t_grid[i]:.4f}: |M(it)| = {abs_M_imag[i]:.4e}, ratio to mean = "
                  f"{abs_M_imag[i] / abs_M_imag.mean():.4f}")
            found_min += 1
            if found_min >= 5:
                break
    if found_min == 0:
        print(f"#     (no significant minima found; |M(it)| is smooth)")
    # Decay rate
    decay_slope, decay_intercept = np.polyfit(t_grid[10:], log_abs[10:], 1)
    print(f"#   Asymptotic decay |M(it)| ~ exp({decay_slope:+.4f}·t) for t > 1")

    # === Step 3: Critical-line search s = 1/2 + it ===
    print(f"\n# Mellin M(1/2 + it) on critical line (Riemann-style):")
    s_grid_crit = 0.5 + 1j * t_grid
    M_crit = np.array([mellin_discrete(P_global, sigma_values, s) for s in s_grid_crit])
    abs_M_crit = np.abs(M_crit)
    print(f"#   |M(1/2 + it)| at sample points:")
    for t in [0, 1, 2, 5, 10, 20, 30]:
        idx = int(t / 30 * 300)
        print(f"#     t={t:>4}: |M(1/2+it)| = {abs_M_crit[idx]:.4e}")
    # Look for zeros on critical line
    log_abs_c = np.log(abs_M_crit + 1e-30)
    found_min_c = 0
    print(f"#   Local minima on critical line (potential zeros):")
    for i in range(2, len(log_abs_c) - 2):
        if (log_abs_c[i] < log_abs_c[i - 1] and log_abs_c[i] < log_abs_c[i + 1] and
            log_abs_c[i] < log_abs_c.mean() - 1):
            print(f"#     t={t_grid[i]:.4f}: |M| = {abs_M_crit[i]:.4e}")
            found_min_c += 1
            if found_min_c >= 5:
                break
    if found_min_c == 0:
        print(f"#     (no significant minima on critical line)")

    # === Step 4: Per-residue Mellin uniformity test ===
    print(f"\n# Per-residue Mellin uniformity at fixed s values:")
    print(f"#   Compare M(s, r) across r ∈ {{1, 3, ..., 31}}")
    test_s_values = [0.0, 1.0, 2.0, 1j, 5j, 0.5+10j]
    print(f"#   {'r':>3}  ", "  ".join(f"|M({s:.1f})|" for s in test_s_values[:4]) + "  ⟨σ|r⟩")
    Ms_per_r = {s: [] for s in test_s_values}
    for r in range(1, 32, 2):
        sg = sigma_per_r[r]
        if len(sg) < 100:
            continue
        smax = int(sg.max())
        sv = np.arange(1, smax + 1, dtype=np.float64)
        cnt = np.bincount(sg, minlength=smax + 1)[1:]
        P_r = cnt / cnt.sum()
        line = f"#   {r:>3}  "
        for s in test_s_values[:4]:
            M_rs = mellin_discrete(P_r, sv, s)
            Ms_per_r[s].append(abs(M_rs))
            line += f"{abs(M_rs):>7.4f}  "
        line += f"{sg.mean():>6.2f}"
        print(line)
    # Spread of |M| across residues for each s
    print(f"\n#   Spread of |M(s, r)| across residues (CV = sd/mean):")
    for s in test_s_values:
        Ms = np.array(Ms_per_r.get(s, []))
        if len(Ms) > 0:
            cv = float(Ms.std() / Ms.mean()) if Ms.mean() > 0 else 0
            print(f"#     s={s}: CV = {cv:.4f}, range = [{Ms.min():.4f}, {Ms.max():.4f}]")

    # === Step 5: Per-σ-band Mellin ===
    print(f"\n# Per-σ-band Mellin signature:")
    print(f"#   {'q':>6}  {'⟨σ|q⟩':>8}  {'M(2)/⟨σ⟩':>10}  {'|M(it)| decay rate':>20}")
    for q in [0.125, 0.375, 0.625, 0.875, 0.975]:
        sg = sigma_per_band[q]
        if len(sg) < 100:
            continue
        smax = int(sg.max())
        sv = np.arange(1, smax + 1, dtype=np.float64)
        cnt = np.bincount(sg, minlength=smax + 1)[1:]
        P_q = cnt / cnt.sum()
        # Decay rate
        s_g = 1j * t_grid
        M_q = np.array([mellin_discrete(P_q, sv, s) for s in s_g])
        log_a = np.log(np.abs(M_q) + 1e-30)
        decay, _ = np.polyfit(t_grid[10:], log_a[10:], 1)
        print(f"#   {q:>6.3f}  {sg.mean():>8.2f}  {mellin_discrete(P_q, sv, 2.0).real / sg.mean():>10.4f}  "
              f"{decay:>+12.4f}")

    # === Step 6: Compare to known Mellin structures ===
    # Gaussian (with mean μ, sd σ): M(s) ≈ μ^(s-1) · exp((s-1)²σ²/(2μ²)) for ⟨σ⟩ ≫ sd
    # Try ratio M(s) / [⟨σ⟩^(s-1) · "Gaussian-like factor"]
    print(f"\n# Compare empirical Mellin to Gaussian baseline:")
    mu = float(all_sigma.mean())
    sigma_sd = float(all_sigma.std(ddof=1))
    print(f"#   ⟨σ⟩ = {mu:.4f}, sd = {sigma_sd:.4f}, CV = {sigma_sd/mu:.4f}")
    print(f"#   For Gaussian with these moments: M(s) ≈ μ^(s-1) · MGF[(s-1)/μ * Gaussian] ...")
    print(f"#   {'s':>4}  {'M_emp(s)':>14}  {'M_Gauss(s)':>14}  {'ratio':>10}")
    for s_real in [0.5, 1.5, 2.0, 2.5, 3.0]:
        M_emp = mellin_discrete(P_global, sigma_values, s_real).real
        # Gaussian Mellin (continuous approx via lognormal MGF or direct):
        # M_gauss(s) = E[X^(s-1)] for X~N(μ, σ): no closed form for general s
        # Use numerical: integrate Normal PDF * x^(s-1)
        from scipy.integrate import quad
        from scipy.stats import norm
        def integrand(x, s=s_real, m=mu, sd=sigma_sd):
            return norm.pdf(x, m, sd) * x ** (s - 1)
        try:
            M_g, _ = quad(integrand, max(1, mu - 5*sigma_sd), mu + 5*sigma_sd)
        except Exception:
            M_g = float('nan')
        ratio = M_emp / M_g if M_g != 0 and not np.isnan(M_g) else float('nan')
        print(f"#   {s_real:>4.2f}  {M_emp:>14.4e}  {M_g:>14.4e}  {ratio:>10.4f}")

    # Save summary
    rows = []
    for r in range(1, 32, 2):
        sg = sigma_per_r[r]
        if len(sg) < 100:
            continue
        smax = int(sg.max())
        sv = np.arange(1, smax + 1, dtype=np.float64)
        cnt = np.bincount(sg, minlength=smax + 1)[1:]
        P_r = cnt / cnt.sum()
        rows.append({
            'r_mod_32': r,
            'n_orbits': len(sg),
            'mean_sigma': float(sg.mean()),
            'M_at_s_0': float(mellin_discrete(P_r, sv, 0.0).real),
            'M_at_s_1': float(mellin_discrete(P_r, sv, 1.0).real),
            'M_at_s_2': float(mellin_discrete(P_r, sv, 2.0).real),
            'absM_at_5i': float(abs(mellin_discrete(P_r, sv, 5j))),
            'absM_at_10i': float(abs(mellin_discrete(P_r, sv, 10j))),
            'absM_at_half_plus_10i': float(abs(mellin_discrete(P_r, sv, 0.5 + 10j))),
        })
    pl.DataFrame(rows).write_csv("experiments_output/mellin_per_residue.csv")
    print(f"\n[save] experiments_output/mellin_per_residue.csv")


if __name__ == "__main__":
    main()
