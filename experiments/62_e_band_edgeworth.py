"""
Test third-moment Edgeworth correction for E_band(q).

Linear Gaussian: E[V | band] = μ_V + ρ·σ_V · E[Z_σ | band]
With Edgeworth correction: + σ_V · C · E[Z_σ² − 1 | band]

If the correction term captures the asymmetric Geom-skew residual, E_band(q)
closes in two terms.

Method:
  1. Walk orbits at N=2^36 (and 2^32, 2^34, 2^38 for stability)
  2. Compute per-orbit (V, σ) → joint moments (means, vars, cov, cross-cumulants)
  3. For each band q, compute:
       E[V | band]_emp
       Linear pred: μ_V + ρ·σ_V·E[Z_σ|band]
       Edgeworth correction term: σ_V·C·E[Z_σ²-1|band]
       Empirical "correction needed" = emp - linear
       Test: is empirical correction ∝ E[Z²-1|band]? (slope = σ_V·C)
  4. Fit C from empirical and test stability/structure
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from scipy.stats import norm

import io
sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

LOG_2 = np.log(2.0); LOG_3 = np.log(3.0)
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_v_only(starts, max_value, max_syr):
    n = len(starts)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    n_odd_arr = np.zeros(n, dtype=np.int32)
    ok_arr = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        m = np.int64(starts[i])
        sigma_total = 0; syr = 0
        failed = False
        while m != 1 and syr < max_syr:
            if (m & 1) == 0:
                m = m >> 1; sigma_total += 1; continue
            if m > max_value // 3:
                failed = True; break
            x = 3*m + 1; v = 0
            while (x & 1) == 0:
                x >>= 1; v += 1
            m = x; sigma_total += 1 + v; syr += 1
        if not failed and m == 1:
            sigma_arr[i] = sigma_total
            n_odd_arr[i] = syr
            ok_arr[i] = True
    return sigma_arr, n_odd_arr, ok_arr


def trunc_z_moments(a, b):
    """E[Z|Z in (a,b)] and E[Z^2-1|Z in (a,b)] for standard Gaussian."""
    Phi_a = norm.cdf(a); Phi_b = norm.cdf(b)
    phi_a = norm.pdf(a); phi_b = norm.pdf(b)
    P = Phi_b - Phi_a
    if P < 1e-15: return 0.0, 0.0
    # E[Z | band] = (phi_a - phi_b)/P
    EZ = (phi_a - phi_b) / P
    # E[Z^2 | band] = 1 + (a*phi_a - b*phi_b)/P
    a_eff = a if np.isfinite(a) else 0.0
    b_eff = b if np.isfinite(b) else 0.0
    a_phi_a = a_eff * phi_a if np.isfinite(a) else 0.0
    b_phi_b = b_eff * phi_b if np.isfinite(b) else 0.0
    EZ2 = 1.0 + (a_phi_a - b_phi_b) / P
    return EZ, EZ2 - 1.0  # (E[Z|band], E[Z^2 - 1|band])


def analyze_N(N, n_per_seed=100_000, seeds=(42,137,271,314,1729)):
    """Walk and compute (V, σ) per-orbit moments + per-band E[V] for 10 bands."""
    log2N = int(np.log2(N))
    print(f"\n# N = 2^{log2N} = {N:,}, walking orbits", flush=True)

    V_all = []; sigma_all = []; logn_all = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        starts = 2 * rng.integers(1, (N - 1)//2, size=n_per_seed, dtype=np.int64) + 1
        sigma, n_odd, ok = walk_v_only(starts, MAX_VAL, 1_000_000)
        starts_ok = starts[ok]
        sigma_ok = sigma[ok].astype(np.float64)
        n_odd_ok = n_odd[ok].astype(np.float64)
        V = (sigma_ok - n_odd_ok) / np.maximum(n_odd_ok, 1)
        V_all.append(V); sigma_all.append(sigma_ok)
        logn_all.append(np.log(starts_ok.astype(np.float64)))

    V = np.concatenate(V_all); sigma = np.concatenate(sigma_all); log_n = np.concatenate(logn_all)
    print(f"  total ok: {len(V):,}", flush=True)

    # Detrend σ by log_n → σ_resid
    log_n_c = log_n - log_n.mean()
    sigma_c = sigma - sigma.mean()
    beta = float((log_n_c * sigma_c).sum() / (log_n_c * log_n_c).sum())
    alpha = float(sigma.mean() - beta * log_n.mean())
    sigma_resid = sigma - (alpha + beta * log_n)

    # Joint moments of (V, σ_resid)
    mu_V = float(V.mean()); mu_S = float(sigma_resid.mean())
    Vc = V - mu_V; Sc = sigma_resid - mu_S
    var_V = float((Vc**2).mean()); var_S = float((Sc**2).mean())
    sd_V = np.sqrt(var_V); sd_S = np.sqrt(var_S)
    cov_VS = float((Vc*Sc).mean())
    rho = cov_VS / (sd_V * sd_S)

    # Third cumulants (in standardized coords for clarity)
    kappa_300_V = float((Vc**3).mean()) / sd_V**3   # skew[V]
    kappa_030_S = float((Sc**3).mean()) / sd_S**3   # skew[σ]
    kappa_111 = float((Vc * Sc**2).mean()) / (sd_V * sd_S**2)  # E[Z_V · Z_σ^2]
    kappa_210 = float((Vc**2 * Sc).mean()) / (sd_V**2 * sd_S)  # E[Z_V^2 · Z_σ]

    print(f"  μ_V = {mu_V:.4f}  σ_V = {sd_V:.4f}  σ_σ = {sd_S:.4f}", flush=True)
    print(f"  ρ(V, σ_resid) = {rho:+.4f}", flush=True)
    print(f"  Std third cumulants:", flush=True)
    print(f"    skew[V] = κ_300 = {kappa_300_V:+.4f}", flush=True)
    print(f"    skew[σ] = κ_030 = {kappa_030_S:+.4f}", flush=True)
    print(f"    κ_111 (E[Z_V·Z_σ²]) = {kappa_111:+.4f}", flush=True)
    print(f"    κ_210 (E[Z_V²·Z_σ]) = {kappa_210:+.4f}", flush=True)

    # Per-band E[V | band] empirical
    band_defs = [
        (0.05,  0.00, 0.10),
        (0.125, 0.10, 0.15),
        (0.20,  0.15, 0.25),
        (0.375, 0.25, 0.50),
        (0.50,  0.45, 0.55),
        (0.625, 0.50, 0.75),
        (0.80,  0.75, 0.85),
        (0.875, 0.85, 0.90),
        (0.95,  0.90, 1.00),
        (0.975, 0.95, 1.00),
    ]
    band_results = []
    for q, lo, hi in band_defs:
        if lo == 0.0:
            lo_val = -np.inf; a_z = -np.inf
        else:
            lo_val = float(np.percentile(sigma_resid, lo*100))
            a_z = (lo_val - mu_S)/sd_S
        if hi == 1.0:
            hi_val = np.inf; b_z = np.inf
        else:
            hi_val = float(np.percentile(sigma_resid, hi*100))
            b_z = (hi_val - mu_S)/sd_S
        mask = (sigma_resid > lo_val) & (sigma_resid <= hi_val)
        n_band = int(mask.sum())
        if n_band == 0: continue

        E_V_band = float(V[mask].mean())

        # Truncated-Gaussian moments
        EZ_band, EZ2m1_band = trunc_z_moments(a_z, b_z)

        # Linear Gaussian prediction
        E_V_lin = mu_V + rho * sd_V * EZ_band

        # "Empirical correction needed" = E[V | band] - linear
        corr_emp = E_V_band - E_V_lin

        # Edgeworth correction term: σ_V · C · E[Z²-1|band], where C involves κ_111
        # The conditional mean Edgeworth: E[Z_V | Z_σ = z] = ρ z + (κ_111 - ρ κ_021)/(2(1-ρ²)) · (z²-1) + ...
        # κ_021 is the standardized E[Z_σ²·Z_V]... wait that's what I called κ_111.
        # Conventions vary; let me just use κ_111 = E[Z_V · Z_σ²] and try simple form:
        # corr_pred = σ_V · (κ_111 / 2) · E[Z²-1|band]
        # (skipping the (1-ρ²) factor for simplicity; will fit slope empirically)
        edgeworth_factor = sd_V * (kappa_111 / 2.0)
        corr_pred_simple = edgeworth_factor * EZ2m1_band

        band_results.append({
            'q': q, 'lo': lo, 'hi': hi,
            'a_z': a_z, 'b_z': b_z,
            'EZ_band': EZ_band, 'EZ2m1_band': EZ2m1_band,
            'n_band': n_band,
            'E_V_band_emp': E_V_band,
            'E_V_lin_Gauss': E_V_lin,
            'corr_emp': corr_emp,
            'corr_pred_edgeworth': corr_pred_simple,
            'E_V_pred_with_edge': E_V_lin + corr_pred_simple,
        })

    return {
        'log2N': log2N, 'mu_V': mu_V, 'sd_V': sd_V, 'sd_S': sd_S, 'rho': rho,
        'skew_V': kappa_300_V, 'skew_S': kappa_030_S,
        'kappa_111': kappa_111, 'kappa_210': kappa_210,
        'bands': band_results,
    }


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"

    all_results = {}
    for log2N in [32, 34, 36, 38]:
        N = 1 << log2N
        r = analyze_N(N)
        all_results[log2N] = r

    # Per-N report
    print(f"\n\n=== Per-N joint cumulants ===", flush=True)
    print(f"  {'log2N':>6}  {'μ_V':>7}  {'σ_V':>7}  {'σ_σ':>8}  {'ρ':>8}  "
          f"{'skew_V':>9}  {'skew_σ':>9}  {'κ_111':>9}  {'κ_210':>9}", flush=True)
    for log2N, r in all_results.items():
        print(f"  {log2N:>6}  {r['mu_V']:>7.4f}  {r['sd_V']:>7.4f}  {r['sd_S']:>8.2f}  "
              f"{r['rho']:>+7.4f}  {r['skew_V']:>+8.4f}  {r['skew_S']:>+8.4f}  "
              f"{r['kappa_111']:>+8.4f}  {r['kappa_210']:>+8.4f}", flush=True)

    # Per-band Edgeworth verification at N=2^36 (most precise)
    r_main = all_results[36]
    print(f"\n\n=== Per-band Edgeworth correction at N=2^36 ===", flush=True)
    print(f"  {'q':>6}  {'a_z':>8}  {'b_z':>8}  {'EZ_band':>9}  {'EZ²-1':>9}  "
          f"{'E_V_emp':>8}  {'E_V_lin':>8}  {'corr_emp':>10}  {'corr_pred':>10}  {'gap':>9}", flush=True)
    for b in r_main['bands']:
        gap = b['corr_emp'] - b['corr_pred_edgeworth']
        print(f"  {b['q']:>6.3f}  {b['a_z']:>+7.3f}  {b['b_z']:>+7.3f}  "
              f"{b['EZ_band']:>+8.3f}  {b['EZ2m1_band']:>+8.3f}  "
              f"{b['E_V_band_emp']:>8.4f}  {b['E_V_lin_Gauss']:>8.4f}  "
              f"{b['corr_emp']:>+9.4f}  {b['corr_pred_edgeworth']:>+9.4f}  {gap:>+8.4f}", flush=True)

    # Empirical fit of correction slope C
    print(f"\n=== Empirical fit: corr_emp = C · E[Z²-1|band] ===", flush=True)
    EZ2m1 = np.array([b['EZ2m1_band'] for b in r_main['bands']])
    corr_emp = np.array([b['corr_emp'] for b in r_main['bands']])
    # OLS through origin
    C_origin = float((EZ2m1 * corr_emp).sum() / (EZ2m1**2).sum())
    pred = C_origin * EZ2m1
    R2_origin = 1 - float(np.sum((corr_emp - pred)**2) / np.sum((corr_emp - corr_emp.mean())**2))
    # Free intercept OLS
    xc = EZ2m1 - EZ2m1.mean(); yc = corr_emp - corr_emp.mean()
    slope = float((xc*yc).sum() / (xc*xc).sum())
    intercept = float(corr_emp.mean() - slope*EZ2m1.mean())
    pred_free = intercept + slope*EZ2m1
    R2_free = 1 - float(np.sum((corr_emp - pred_free)**2) / np.sum((corr_emp - corr_emp.mean())**2))
    print(f"  Through-origin: C = {C_origin:+.5f}  R² = {R2_origin:.4f}", flush=True)
    print(f"  Free intercept: slope = {slope:+.5f}  intercept = {intercept:+.5f}  R² = {R2_free:.4f}", flush=True)
    # Predicted slope from edgeworth: σ_V · κ_111 / 2
    C_pred = r_main['sd_V'] * r_main['kappa_111'] / 2.0
    print(f"  Edgeworth prediction: σ_V·κ_111/2 = {C_pred:+.5f}", flush=True)
    print(f"  Empirical / prediction ratio: {C_origin / C_pred if abs(C_pred) > 1e-6 else 'inf':>+.3f}", flush=True)

    # Stability: same C fit at different N
    print(f"\n=== Stability of empirical C across N ===", flush=True)
    print(f"  {'log2N':>6}  {'σ_V':>7}  {'κ_111':>9}  {'σ_V·κ_111/2':>13}  "
          f"{'fit C':>11}  {'R²':>7}", flush=True)
    for log2N, r in all_results.items():
        EZ2m1_N = np.array([b['EZ2m1_band'] for b in r['bands']])
        corr_N = np.array([b['corr_emp'] for b in r['bands']])
        C_N = float((EZ2m1_N * corr_N).sum() / (EZ2m1_N**2).sum())
        pred_N = C_N * EZ2m1_N
        R2_N = 1 - float(np.sum((corr_N - pred_N)**2) / np.sum((corr_N - corr_N.mean())**2))
        C_pred_N = r['sd_V'] * r['kappa_111'] / 2.0
        print(f"  {log2N:>6}  {r['sd_V']:>7.4f}  {r['kappa_111']:>+8.4f}  "
              f"{C_pred_N:>+13.5f}  {C_N:>+10.5f}  {R2_N:>7.4f}", flush=True)

    # Save
    rows = []
    for log2N, r in all_results.items():
        for b in r['bands']:
            rows.append({
                'log2N': log2N,
                'q': b['q'], 'a_z': b['a_z'], 'b_z': b['b_z'],
                'EZ_band': b['EZ_band'], 'EZ2m1_band': b['EZ2m1_band'],
                'mu_V': r['mu_V'], 'sd_V': r['sd_V'], 'rho': r['rho'],
                'kappa_111': r['kappa_111'],
                'E_V_emp': b['E_V_band_emp'],
                'E_V_lin_Gauss': b['E_V_lin_Gauss'],
                'corr_emp': b['corr_emp'],
                'corr_pred_edgeworth': b['corr_pred_edgeworth'],
                'E_V_pred_with_edge': b['E_V_pred_with_edge'],
                'gap_edge_minus_emp': b['corr_pred_edgeworth'] - b['corr_emp'],
            })
    pl.DataFrame(rows).write_csv(out_dir / "62_e_band_edgeworth.csv")
    print(f"\n[save] CSV", flush=True)


if __name__ == "__main__":
    main()
