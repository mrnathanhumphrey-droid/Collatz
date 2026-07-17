"""
qx_plus_1_framework.py

Systematic bridge-equation pass for q in {3, 5, 7, 9}.
Loads existing data/q_main_q*_N*.parquet (already generated).

Per q, for CONVERGED-to-1 orbits only:
  Constant 1 analog: per-class alpha_det offset structure (R^2 of OLS sigma vs alpha_det+log)
  Constant 2 analog: empirical slope of <sigma|class> on log(n) -- the K_h(q) analog
  Closed-form candidate: K_h(q)? = (q+1)/(2*log(2)) -- inverse of average log-step magnitude

For q in {5, 7, 9}, only converged orbits exist (conv_rate ~ 10^-3 to 10^-6).
For q=3, all orbits converge.

Cap: <300 lines.
"""
import sys
import numpy as np
import polars as pl
from pathlib import Path

DATA = Path(r"C:\Collatz\data")
OUT = Path(r"C:\Collatz")

sys.stdout.reconfigure(encoding="utf-8")
log_lines = []


def log(s):
    print(s)
    log_lines.append(s)


# Largest available parquet per q
PARQUETS = {
    3: DATA / "q_main_q3_N1000000.parquet",
    5: DATA / "q_main_q5_N100000000.parquet",
    7: DATA / "q_main_q7_N100000000.parquet",
    9: DATA / "q_main_q9_N100000000.parquet",
}


def fit_per_class_slope(df_conv, k=6):
    """
    For converged orbits only, fit sigma ~ a + b*log(n) within each residue class mod 2^k.
    Return per-class (intercept, slope, n_class) and pooled slope on the full sample.
    """
    M = 1 << k
    df = df_conv.with_columns([
        (pl.col("n") % M).alias(f"r_{k}"),
        pl.col("n").log().alias("log_n"),
        pl.col("sigma_q").cast(pl.Float64).alias("sigma_f"),
    ])
    log_n = df["log_n"].to_numpy()
    sigma = df["sigma_f"].to_numpy()
    r_k = df[f"r_{k}"].to_numpy()

    # Pooled
    A_pool = np.vstack([np.ones_like(log_n), log_n]).T
    b_pool, slope_pool = np.linalg.lstsq(A_pool, sigma, rcond=None)[0]

    # Per-class
    classes = sorted(set(r_k.tolist()))
    rows = []
    for c in classes:
        mask = r_k == c
        n_c = mask.sum()
        if n_c < 50:
            continue
        ln_c = log_n[mask]
        sg_c = sigma[mask]
        if ln_c.std() < 1e-9:
            continue
        A = np.vstack([np.ones_like(ln_c), ln_c]).T
        try:
            sol = np.linalg.lstsq(A, sg_c, rcond=None)[0]
            rows.append({"class": int(c), "n": int(n_c), "intercept": float(sol[0]),
                         "slope": float(sol[1]), "mean_sigma": float(sg_c.mean()),
                         "mean_log_n": float(ln_c.mean())})
        except Exception:
            pass
    return rows, float(b_pool), float(slope_pool)


def main():
    log("=" * 78)
    log("qx+1 bridge equation framework — systematic pass for q in {3, 5, 7, 9}")
    log("=" * 78)

    LOG2 = np.log(2.0)
    results = {}

    for q in [3, 5, 7, 9]:
        path = PARQUETS[q]
        if not path.exists():
            log(f"\n[q={q}] parquet missing: {path}")
            continue
        log(f"\n{'='*60}")
        log(f"[q={q}] loading {path.name}")
        df = pl.read_parquet(path)
        n_total = df.height
        n_conv = df.filter(pl.col("converged")).height
        log(f"  total odd n: {n_total:,}; converged-to-1: {n_conv:,} "
            f"({100*n_conv/n_total:.4f}%)")

        if n_conv < 200:
            log(f"  WARNING: too few converged orbits for slope fit; skipping")
            continue
        df_c = df.filter(pl.col("converged"))

        # --- Constant 2 analog: K_h(q) ---
        log(f"\n  === Constant 2 analog: K_h(q) — slope of sigma vs log(n) on converged orbits ===")
        rows, intercept_pool, slope_pool = fit_per_class_slope(df_c, k=6)
        log(f"  Pooled OLS:  sigma = {intercept_pool:+.4f} + {slope_pool:.4f}*log(n)")

        # Closed-form candidate: random-walk slope = 1/|E[X]| with E[X]=log(q/4) and full step count.
        # For q=3: K_h = 3/log(4/3) = 10.43. Generalization?
        # For q=3: each odd step accompanied by E[v]=2 even steps; total slope = (1+E[v])/|E[X_odd]|
        #   E[X_odd] = log(q) - E[v]·log(2); for q=3: log(3) - 2log(2) = -log(4/3); so K = 3/log(4/3)
        # For q>=5 conditional-on-convergence: orbits diverge on AVERAGE, so converged orbits are
        #   those with anomalously LARGE v's (Esscher-tilted v-distribution per Result 18 / Cramer).
        # Empirical Cramer-tilted E*[v] for q=5: 2.92 (agent2 task 3). Then K(q=5; conv) approx =
        #   (1+2.92)/(2.92*log(2) - log(5)) = 3.92/(2.024 - 1.609) = 9.45.
        # For q=7 conv: tilted E*[v] would be larger. For q=9 even more.

        # Predicted K_h(q) under the "conditional Cramer-tilted E*[v]" model:
        # E*[v] approx (log m_start + J*log(q)) / (J*log(2)) for typical converged orbit
        # where J = mean odd_steps. Use empirical odd_steps/sigma to get E*[v].
        odd_s_mean = float(df_c["odd_steps"].mean())
        even_s_mean = float(df_c["even_steps"].mean())
        sigma_mean = float(df_c["sigma_q"].mean())
        E_v_emp = even_s_mean / odd_s_mean if odd_s_mean > 0 else float("nan")
        log(f"  Empirical mean odd_steps={odd_s_mean:.2f}, even_steps={even_s_mean:.2f}")
        log(f"  Implied E*[v]_conv = even/odd = {E_v_emp:.4f}")
        log(f"  Implied step drift mu* = E*[v]*log(2) - log(q) = "
            f"{E_v_emp*LOG2 - np.log(q):+.4f} (POSITIVE means descending)")
        K_pred = (1.0 + E_v_emp) / (E_v_emp * LOG2 - np.log(q)) if E_v_emp*LOG2 > np.log(q) else float("nan")
        log(f"  Closed-form prediction K_h(q;conv) = (1+E*[v]) / (E*[v]*log(2) - log(q)) = {K_pred:.4f}")
        log(f"  vs empirical pooled slope                                                = {slope_pool:.4f}")
        if np.isfinite(K_pred):
            log(f"  Match: gap = {slope_pool - K_pred:+.4f} ({100*(slope_pool-K_pred)/K_pred:+.2f}%)")

        # Per-class slope spread
        if rows:
            slopes = np.array([r["slope"] for r in rows])
            intercepts = np.array([r["intercept"] for r in rows])
            log(f"  Per-class slope (k=6, {len(rows)} classes ≥50 obs): "
                f"median={np.median(slopes):.4f}, std={slopes.std():.4f}, range=[{slopes.min():.3f}, {slopes.max():.3f}]")
            log(f"  Per-class intercept: median={np.median(intercepts):.4f}, std={intercepts.std():.4f}")

        # --- Constant 1 analog: alpha_det offset structure ---
        # Test universality: are per-class slopes ~ pooled, with all class-difference in intercepts?
        if rows and len(rows) >= 5:
            slope_consistency = slopes.std() / abs(slopes.mean()) if abs(slopes.mean()) > 1e-9 else float("nan")
            log(f"\n  === Constant 1 analog: per-class intercept structure ===")
            log(f"  Slope consistency CV = std/|mean| = {slope_consistency:.4f}")
            log(f"  Intercept range: {intercepts.max() - intercepts.min():.4f} (variation across classes)")
            if slope_consistency < 0.05:
                log(f"  --> Slopes are q-universal across classes; class structure lives in intercepts (CONST-1 analog HOLDS)")
            else:
                log(f"  --> Slopes vary ~{slope_consistency*100:.1f}% across classes; CONST-1 analog less clean")

        results[q] = {
            "n_conv": n_conv,
            "pooled_slope": slope_pool,
            "pooled_intercept": intercept_pool,
            "K_pred": K_pred,
            "E_v_emp": E_v_emp,
            "odd_steps_mean": odd_s_mean,
            "sigma_mean": sigma_mean,
        }

    # ---------- Cross-q comparison ----------
    log(f"\n{'='*78}")
    log(f"Cross-q summary")
    log(f"{'='*78}")
    log(f"\n{'q':>3} {'n_conv':>10} {'sigma_mean':>12} {'E*[v]_conv':>12} "
        f"{'K_h(q;conv)':>14} {'pooled_slope':>14} {'gap':>10}")
    for q, r in results.items():
        gap = r["pooled_slope"] - r["K_pred"] if np.isfinite(r["K_pred"]) else float("nan")
        log(f"{q:>3} {r['n_conv']:>10} {r['sigma_mean']:>12.2f} {r['E_v_emp']:>12.4f} "
            f"{r['K_pred']:>14.4f} {r['pooled_slope']:>14.4f} {gap:>+10.4f}")

    # Universal Cramer multiplier check (already shipped by user; cross-reference)
    log(f"\n[Cramer multiplier check (cross-reference, already shipped)]")
    log(f"  Universal C ~= 5/2 confirmed for q in {{5, 7, 9}}, R^2 >= 0.99 (findings.md 2026-05-01 late)")
    log(f"  This is the q-universal piece of the convergence-rate slope.")

    # Save
    rows_out = []
    for q, r in results.items():
        row = {"q": q}
        row.update(r)
        rows_out.append(row)
    pl.DataFrame(rows_out).write_csv(OUT / "qx_plus_1_framework.csv")

    log_path = OUT / "qx_plus_1_framework_log.txt"
    log_path.write_text("\n".join(log_lines), encoding="utf-8")
    log(f"\n[wrote] qx_plus_1_framework.csv, qx_plus_1_framework_log.txt")


if __name__ == "__main__":
    main()
