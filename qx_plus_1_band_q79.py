"""
qx_plus_1_band_q79.py — Constant 4 band test at q=7 (N=10^9) and q=9 (N=10^8).

q=7 N=10^9: expected ~2.5K conv orbits (8x more than N=10^8's 258)
q=9 N=10^8: 104 conv orbits — use 2-band stratification (median split)

Test if U-shape per-band slope pattern reproduces at q=7, q=9.
"""
import sys
import numpy as np
import polars as pl
from pathlib import Path

DATA = Path(r"C:\Collatz\data")
OUT = Path(r"C:\Collatz")
sys.stdout.reconfigure(encoding="utf-8")
log_lines = []
def log(s): print(s); log_lines.append(s)


def per_band_slope(df_conv, n_bands):
    log_n = df_conv["n"].log().to_numpy()
    sigma = df_conv["sigma_q"].cast(pl.Float64).to_numpy()
    qs = np.quantile(sigma, [i/n_bands for i in range(1, n_bands)])
    bands = np.digitize(sigma, qs)
    rows = []
    for b in range(n_bands):
        mask = bands == b
        n_b = mask.sum()
        if n_b < 20:
            continue
        ln = log_n[mask]; sg = sigma[mask]
        if ln.std() < 1e-9:
            continue
        A = np.vstack([np.ones_like(ln), ln]).T
        sol, _, _, _ = np.linalg.lstsq(A, sg, rcond=None)
        # SE on slope: sqrt(MSE / sum((x-x_mean)^2))
        resid = sg - A @ sol
        mse = (resid**2).sum() / max(1, len(sg) - 2)
        se_slope = float(np.sqrt(mse / ((ln - ln.mean())**2).sum()))
        rows.append({
            "band": b, "n_band": int(n_b),
            "intercept": float(sol[0]), "slope": float(sol[1]), "se_slope": se_slope,
            "sigma_mean": float(sg.mean()),
        })
    return rows


log("=" * 78)
log("Constant 4 analog at q=7 (N=10^9) and q=9 (N=10^8)")
log("=" * 78)

# q=7 at N=10^9
log(f"\n[q=7 N=10^9]")
df7 = pl.read_parquet(DATA / "q_main_q7_N1000000000.parquet").filter(pl.col("converged"))
log(f"  {df7.height:,} converged orbits")
for nb in (4, 3):
    if df7.height >= 80 * nb:
        rows7 = per_band_slope(df7, n_bands=nb)
        log(f"  {nb}-band stratification:")
        log(f"  band   n_band  σ_mean   slope ± SE")
        for r in rows7:
            log(f"  {r['band']:>4d}  {r['n_band']:>7d}  {r['sigma_mean']:>6.1f}  "
                f"{r['slope']:>+7.4f} ± {r['se_slope']:.4f}")
        slopes = np.array([r['slope'] for r in rows7])
        spread = slopes.max() - slopes.min()
        log(f"  spread = {spread:.3f}")
        # U-shape check: low+high higher than middle
        if nb == 4:
            avg_ends = (slopes[0] + slopes[3]) / 2
            avg_mid = (slopes[1] + slopes[2]) / 2
            log(f"  U-shape diagnostic: avg(ends) = {avg_ends:.3f}, avg(middle) = {avg_mid:.3f}, "
                f"diff = {avg_ends - avg_mid:+.3f}")
            log(f"  --> {'U-shape PRESENT' if avg_ends > avg_mid + 0.5 else 'U-shape ABSENT/WEAK'}")
        elif nb == 3:
            log(f"  Tri-band: ends {slopes[0]:.3f}, {slopes[2]:.3f}, middle {slopes[1]:.3f}")
            log(f"  --> {'U-shape PRESENT' if max(slopes[0], slopes[2]) > slopes[1] + 0.5 else 'U-shape ABSENT/WEAK'}")
        break

# q=9 at N=10^8
log(f"\n[q=9 N=10^8]")
df9 = pl.read_parquet(DATA / "q_main_q9_N100000000.parquet").filter(pl.col("converged"))
log(f"  {df9.height:,} converged orbits")
for nb in (4, 3, 2):
    if df9.height >= 30 * nb:
        rows9 = per_band_slope(df9, n_bands=nb)
        log(f"  {nb}-band stratification:")
        log(f"  band   n_band  σ_mean   slope ± SE")
        for r in rows9:
            log(f"  {r['band']:>4d}  {r['n_band']:>7d}  {r['sigma_mean']:>6.1f}  "
                f"{r['slope']:>+7.4f} ± {r['se_slope']:.4f}")
        slopes = np.array([r['slope'] for r in rows9])
        spread = slopes.max() - slopes.min()
        log(f"  spread = {spread:.3f}")
        if nb == 4:
            avg_ends = (slopes[0] + slopes[3]) / 2
            avg_mid = (slopes[1] + slopes[2]) / 2
            log(f"  U-shape diagnostic: avg(ends) = {avg_ends:.3f}, avg(middle) = {avg_mid:.3f}, "
                f"diff = {avg_ends - avg_mid:+.3f}")
            log(f"  --> {'U-shape PRESENT' if avg_ends > avg_mid + 0.5 else 'U-shape ABSENT/WEAK'}")
        elif nb == 3:
            log(f"  Tri-band: ends {slopes[0]:.3f}, {slopes[2]:.3f}, middle {slopes[1]:.3f}")
            log(f"  --> {'U-shape PRESENT' if max(slopes[0], slopes[2]) > slopes[1] + 0.5 else 'U-shape ABSENT/WEAK'}")
        elif nb == 2:
            log(f"  2-band: low {slopes[0]:.3f}, high {slopes[1]:.3f} (insufficient for U-shape test)")
        break

log(f"\n[Cross-q U-shape pattern (low+high vs middle)]")
log(f"  q=3: U-shape PRESENT (band 0,3 = +2.90, +5.48; band 1,2 = -0.10, +0.58; spread 5.58)")
log(f"  q=5: U-shape PRESENT (band 0,3 = +2.73, +3.80; band 1,2 = +0.01, +0.96; spread 3.80)")

(OUT / "qx_plus_1_band_q79_log.txt").write_text("\n".join(log_lines), encoding="utf-8")
log(f"\n[wrote] qx_plus_1_band_q79_log.txt")
