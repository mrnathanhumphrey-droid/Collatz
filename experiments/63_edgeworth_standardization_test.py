"""
Test standardization correction for Edgeworth coefficient C in E_band(q).

Result 27/29 established empirical C ~ 0.04, with C/(σ_V·κ_111/2) ratio INCREASING
in N (8.30 → 10.81 → 13.31 → 15.92 across N = 2^32, 2^34, 2^36, 2^38).

Hypothesis: standardization is the issue. σ has trivial K_h·log(n) drift PLUS
structural fluctuation. Standardization A uses raw σ; Standardization B uses
σ_resid = σ − K_h·log(n) (theoretical K_h, not empirical β).

This test:
  - Recompute κ_111 under both standardizations at each N
  - Recompute predicted C and empirical C for each
  - If C_emp / C_pred_B is N-stable: standardization is the issue (outcome a or b)
  - If still N-varying: standardization isn't the issue (outcome c)
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
K_H = 3.0 / np.log(4.0/3.0)  # ~ 10.4185, theoretical Syracuse drift
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
    EZ = (phi_a - phi_b) / P
    a_phi_a = a * phi_a if np.isfinite(a) else 0.0
    b_phi_b = b * phi_b if np.isfinite(b) else 0.0
    EZ2 = 1.0 + (a_phi_a - b_phi_b) / P
    return EZ, EZ2 - 1.0


def fit_band_C(EZ_arr, EZ2m1_arr, corr_emp_arr):
    """Fit corr_emp = ρ_term + C·EZ2m1 (we already subtracted ρ_term)."""
    # Through origin
    denom = float((EZ2m1_arr**2).sum())
    if denom < 1e-15:
        return 0.0, 0.0
    C = float((EZ2m1_arr * corr_emp_arr).sum() / denom)
    pred = C * EZ2m1_arr
    ss_res = float(((corr_emp_arr - pred)**2).sum())
    ss_tot = float(((corr_emp_arr - corr_emp_arr.mean())**2).sum())
    R2 = 1 - ss_res/ss_tot if ss_tot > 1e-15 else 1.0
    return C, R2


def analyze_N_with_both_standardizations(N, n_per_seed=100_000, seeds=(42,137,271,314,1729)):
    """Walk and compute joint moments + per-band E[V] under BOTH standardizations."""
    log2N = int(np.log2(N))
    print(f"\n# N = 2^{log2N} = {N:,}, walking orbits", flush=True)

    V_all = []; sigma_all = []; logn_all = []
    t0 = time.time()
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

    V = np.concatenate(V_all); sigma_raw = np.concatenate(sigma_all); log_n = np.concatenate(logn_all)
    print(f"  total ok: {len(V):,}  walk-time: {time.time()-t0:.1f}s", flush=True)

    # ============== STANDARDIZATION A: raw σ ==============
    mu_sA = float(sigma_raw.mean()); sd_sA = float(sigma_raw.std())
    Z_sA = (sigma_raw - mu_sA) / sd_sA

    # ============== STANDARDIZATION B: σ - K_h·log_n ==============
    sigma_residB = sigma_raw - K_H * log_n
    # Allow intercept (theoretical K_h gives the slope; intercept is set by data)
    mu_sB = float(sigma_residB.mean()); sd_sB = float(sigma_residB.std())
    Z_sB = (sigma_residB - mu_sB) / sd_sB

    # Also: B' empirical β regression (= what exp 62 did) for cross-check
    log_n_c = log_n - log_n.mean()
    sigma_c = sigma_raw - sigma_raw.mean()
    beta_emp = float((log_n_c * sigma_c).sum() / (log_n_c * log_n_c).sum())
    intercept_emp = float(sigma_raw.mean() - beta_emp * log_n.mean())
    sigma_residBp = sigma_raw - (intercept_emp + beta_emp * log_n)
    mu_sBp = float(sigma_residBp.mean()); sd_sBp = float(sigma_residBp.std())
    Z_sBp = (sigma_residBp - mu_sBp) / sd_sBp

    # V moments (same regardless of standardization)
    mu_V = float(V.mean()); sd_V = float(V.std())
    Z_V = (V - mu_V) / sd_V

    print(f"  V: μ={mu_V:.4f}  σ={sd_V:.4f}", flush=True)
    print(f"  Standardization A (raw σ): mean={mu_sA:.2f}  SD={sd_sA:.2f}", flush=True)
    print(f"  Standardization B (σ-K_h·log_n, K_h={K_H:.4f}): mean={mu_sB:.2f}  SD={sd_sB:.2f}", flush=True)
    print(f"  Standardization B' (σ-β_emp·log_n, β_emp={beta_emp:.4f}): mean={mu_sBp:.2f}  SD={sd_sBp:.2f}", flush=True)

    # Joint third cumulants under each standardization
    kappa_111_A = float((Z_V * Z_sA**2).mean())
    kappa_111_B = float((Z_V * Z_sB**2).mean())
    kappa_111_Bp = float((Z_V * Z_sBp**2).mean())

    # Correlations
    rho_A = float((Z_V * Z_sA).mean())
    rho_B = float((Z_V * Z_sB).mean())
    rho_Bp = float((Z_V * Z_sBp).mean())

    print(f"  ρ(V, Z_σ): A={rho_A:+.4f}  B={rho_B:+.4f}  B'={rho_Bp:+.4f}", flush=True)
    print(f"  κ_111: A={kappa_111_A:+.5f}  B={kappa_111_B:+.5f}  B'={kappa_111_Bp:+.5f}", flush=True)

    # Edgeworth predicted C
    C_pred_A = sd_V * kappa_111_A / 2.0
    C_pred_B = sd_V * kappa_111_B / 2.0
    C_pred_Bp = sd_V * kappa_111_Bp / 2.0

    # Per-band analysis under each standardization
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

    def band_corr_arrays(sigma_proxy, sd_proxy, mu_proxy, rho_use):
        """For bands defined in sigma_proxy quantiles, return EZ_band, EZ2m1_band, corr_emp arrays."""
        EZ_list = []; EZ2m1_list = []; corr_list = []; emp_list = []; lin_list = []
        for q, lo, hi in band_defs:
            if lo == 0.0:
                lo_val = -np.inf; a_z = -np.inf
            else:
                lo_val = float(np.percentile(sigma_proxy, lo*100))
                a_z = (lo_val - mu_proxy) / sd_proxy
            if hi == 1.0:
                hi_val = np.inf; b_z = np.inf
            else:
                hi_val = float(np.percentile(sigma_proxy, hi*100))
                b_z = (hi_val - mu_proxy) / sd_proxy
            mask = (sigma_proxy > lo_val) & (sigma_proxy <= hi_val)
            n_band = int(mask.sum())
            if n_band < 100: continue
            EZ_band, EZ2m1_band = trunc_z_moments(a_z, b_z)
            E_V_band = float(V[mask].mean())
            E_V_lin = mu_V + rho_use * sd_V * EZ_band
            corr_emp = E_V_band - E_V_lin
            EZ_list.append(EZ_band); EZ2m1_list.append(EZ2m1_band)
            corr_list.append(corr_emp); emp_list.append(E_V_band); lin_list.append(E_V_lin)
        return (np.array(EZ_list), np.array(EZ2m1_list), np.array(corr_list),
                np.array(emp_list), np.array(lin_list))

    # Per-band under A (bands defined in raw σ quantiles)
    EZ_A, EZ2m1_A, corr_A, emp_A, lin_A = band_corr_arrays(sigma_raw, sd_sA, mu_sA, rho_A)
    # Per-band under B (bands defined in σ_residB quantiles)
    EZ_B, EZ2m1_B, corr_B, emp_B, lin_B = band_corr_arrays(sigma_residB, sd_sB, mu_sB, rho_B)
    # Per-band under B' (bands defined in σ_residBp quantiles)
    EZ_Bp, EZ2m1_Bp, corr_Bp, emp_Bp, lin_Bp = band_corr_arrays(sigma_residBp, sd_sBp, mu_sBp, rho_Bp)

    C_emp_A, R2_A = fit_band_C(EZ_A, EZ2m1_A, corr_A)
    C_emp_B, R2_B = fit_band_C(EZ_B, EZ2m1_B, corr_B)
    C_emp_Bp, R2_Bp = fit_band_C(EZ_Bp, EZ2m1_Bp, corr_Bp)

    return {
        'log2N': log2N, 'mu_V': mu_V, 'sd_V': sd_V,
        'A': {
            'sd_sigma': sd_sA, 'rho': rho_A, 'kappa_111': kappa_111_A,
            'C_pred': C_pred_A, 'C_emp': C_emp_A, 'R2': R2_A,
            'EZ': EZ_A, 'EZ2m1': EZ2m1_A, 'corr_emp': corr_A,
            'emp': emp_A, 'lin': lin_A,
        },
        'B': {
            'sd_sigma': sd_sB, 'rho': rho_B, 'kappa_111': kappa_111_B,
            'C_pred': C_pred_B, 'C_emp': C_emp_B, 'R2': R2_B,
            'EZ': EZ_B, 'EZ2m1': EZ2m1_B, 'corr_emp': corr_B,
            'emp': emp_B, 'lin': lin_B,
        },
        'Bp': {
            'sd_sigma': sd_sBp, 'rho': rho_Bp, 'kappa_111': kappa_111_Bp,
            'C_pred': C_pred_Bp, 'C_emp': C_emp_Bp, 'R2': R2_Bp,
            'EZ': EZ_Bp, 'EZ2m1': EZ2m1_Bp, 'corr_emp': corr_Bp,
            'emp': emp_Bp, 'lin': lin_Bp,
        },
        'beta_emp': beta_emp,
        'log_n_var': float(log_n.var()),
    }


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(exist_ok=True)

    print(f"K_h (theoretical) = 3/log(4/3) = {K_H:.6f}", flush=True)

    all_results = {}
    for log2N in [32, 34, 36, 38]:
        N = 1 << log2N
        r = analyze_N_with_both_standardizations(N)
        all_results[log2N] = r

    # ============== Master comparison table ==============
    print(f"\n\n=== Per-N comparison: A (raw σ) vs B (σ−K_h·log_n) vs B' (empirical β) ===", flush=True)
    hdr = (f"  {'log2N':>6}  {'σ_V':>7}  "
           f"{'β_emp':>8}  {'Var[ln]':>8}  "
           f"{'sd_A':>7}  {'sd_B':>7}  {'sd_Bp':>7}  "
           f"{'ρ_A':>8}  {'ρ_B':>8}  {'ρ_Bp':>8}")
    print(hdr, flush=True)
    for log2N, r in all_results.items():
        print(f"  {log2N:>6}  {r['sd_V']:>7.4f}  "
              f"{r['beta_emp']:>8.4f}  {r['log_n_var']:>8.4f}  "
              f"{r['A']['sd_sigma']:>7.2f}  {r['B']['sd_sigma']:>7.2f}  {r['Bp']['sd_sigma']:>7.2f}  "
              f"{r['A']['rho']:>+7.4f}  {r['B']['rho']:>+7.4f}  {r['Bp']['rho']:>+7.4f}", flush=True)

    print(f"\n=== κ_111 comparison ===", flush=True)
    print(f"  {'log2N':>6}  {'κ_111_A':>10}  {'κ_111_B':>10}  {'κ_111_Bp':>10}", flush=True)
    for log2N, r in all_results.items():
        print(f"  {log2N:>6}  {r['A']['kappa_111']:>+10.5f}  "
              f"{r['B']['kappa_111']:>+10.5f}  {r['Bp']['kappa_111']:>+10.5f}", flush=True)

    print(f"\n=== Predicted vs empirical C, with ratios ===", flush=True)
    print(f"  {'log2N':>6}  std  {'C_pred':>10}  {'C_emp':>10}  {'ratio':>8}  {'R²':>7}", flush=True)
    for log2N, r in all_results.items():
        for key in ['A', 'B', 'Bp']:
            d = r[key]
            ratio = d['C_emp']/d['C_pred'] if abs(d['C_pred']) > 1e-9 else float('nan')
            print(f"  {log2N:>6}   {key:<2}  {d['C_pred']:>+10.5f}  {d['C_emp']:>+10.5f}  "
                  f"{ratio:>8.3f}  {d['R2']:>7.4f}", flush=True)

    # ============== Decisive ratios across N ==============
    print(f"\n=== Ratio (C_emp / C_pred) across N — N-stable means standardization closes ===", flush=True)
    print(f"  {'log2N':>6}  {'ratio_A':>10}  {'ratio_B':>10}  {'ratio_Bp':>10}", flush=True)
    for log2N, r in all_results.items():
        rA = r['A']['C_emp']/r['A']['C_pred'] if abs(r['A']['C_pred']) > 1e-9 else float('nan')
        rB = r['B']['C_emp']/r['B']['C_pred'] if abs(r['B']['C_pred']) > 1e-9 else float('nan')
        rBp = r['Bp']['C_emp']/r['Bp']['C_pred'] if abs(r['Bp']['C_pred']) > 1e-9 else float('nan')
        print(f"  {log2N:>6}  {rA:>10.3f}  {rB:>10.3f}  {rBp:>10.3f}", flush=True)

    # ============== Verdict ==============
    print(f"\n=== VERDICT ===", flush=True)
    ratios_B = [r['B']['C_emp']/r['B']['C_pred'] if abs(r['B']['C_pred']) > 1e-9 else float('nan')
                for r in all_results.values()]
    ratios_A = [r['A']['C_emp']/r['A']['C_pred'] if abs(r['A']['C_pred']) > 1e-9 else float('nan')
                for r in all_results.values()]
    cv_A = np.std(ratios_A) / abs(np.mean(ratios_A))
    cv_B = np.std(ratios_B) / abs(np.mean(ratios_B))
    print(f"  ratio_A across N: {ratios_A}", flush=True)
    print(f"  ratio_B across N: {ratios_B}", flush=True)
    print(f"  CV(ratio_A) = {cv_A:.3f}", flush=True)
    print(f"  CV(ratio_B) = {cv_B:.3f}", flush=True)
    print(f"  mean ratio_B = {np.mean(ratios_B):.3f}", flush=True)

    if cv_B < 0.10 and abs(np.mean(ratios_B) - 1.0) < 0.10:
        print("  → outcome (a): standardization B closes C as σ_V·κ_111_B/2 (ratio ≈ 1, N-stable)", flush=True)
    elif cv_B < 0.10:
        print(f"  → outcome (b): standardization B closes C up to constant factor {np.mean(ratios_B):.3f}", flush=True)
        print(f"     Test if g(ρ) explains it.", flush=True)
    else:
        print("  → outcome (c): standardization is NOT the issue (ratio still N-varying)", flush=True)
        print("     Need full Stuart-Ord bivariate Edgeworth derivation with multiple cumulants", flush=True)

    # Test (1-ρ²)^k factor candidates if outcome (b) suggested
    if cv_B < 0.20:
        print(f"\n=== Test ρ-dependent factor candidates ===", flush=True)
        for log2N, r in all_results.items():
            rho = r['B']['rho']
            ratio = r['B']['C_emp']/r['B']['C_pred'] if abs(r['B']['C_pred']) > 1e-9 else float('nan')
            print(f"  log2N={log2N}  ρ={rho:+.4f}  ratio={ratio:.3f}  "
                  f"1/(1-ρ²)={1/(1-rho**2):.3f}  "
                  f"1/(1-ρ²)²={1/(1-rho**2)**2:.3f}  "
                  f"ρ²/(1-ρ²)={rho**2/(1-rho**2):.3f}", flush=True)

    # Save CSV
    rows = []
    for log2N, r in all_results.items():
        for key in ['A', 'B', 'Bp']:
            d = r[key]
            for i, q in enumerate([0.05, 0.125, 0.20, 0.375, 0.50, 0.625, 0.80, 0.875, 0.95, 0.975]):
                if i >= len(d['EZ']): continue
                rows.append({
                    'log2N': log2N, 'standardization': key, 'q': q,
                    'sd_V': r['sd_V'], 'mu_V': r['mu_V'],
                    'sd_sigma': d['sd_sigma'], 'rho': d['rho'], 'kappa_111': d['kappa_111'],
                    'C_pred': d['C_pred'], 'C_emp': d['C_emp'], 'R2': d['R2'],
                    'EZ_band': float(d['EZ'][i]), 'EZ2m1_band': float(d['EZ2m1'][i]),
                    'corr_emp': float(d['corr_emp'][i]),
                    'E_V_emp': float(d['emp'][i]), 'E_V_lin': float(d['lin'][i]),
                })
    df = pl.DataFrame(rows)
    out_csv = out_dir / "63_edgeworth_standardization_test.csv"
    df.write_csv(out_csv)
    print(f"\n[save] {out_csv}", flush=True)


if __name__ == "__main__":
    main()
