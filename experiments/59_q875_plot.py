"""Plot K_eff_band(q=0.875) and ξ_X trajectories with 5-seed bootstrap CIs."""
import sys
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")

K_H_USER = 3.0 / (np.log(4.0) - np.log(3.0))

here = Path(__file__).resolve().parent
out_dir = here.parent / "experiments_output"

df = pl.read_csv(out_dir / "59_q875_extension_summary.csv")
log2N = df["log2N"].to_numpy()
K_q875 = df["K_q875_mean"].to_numpy()
K_q875_se = df["K_q875_se"].to_numpy()
xi_X = df["xi_X_mean"].to_numpy()
xi_X_se = df["xi_X_se"].to_numpy()
K_full = df["K_full_mean"].to_numpy()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Top: K_eff_band(q=0.875) with CIs
ax1.errorbar(log2N, K_q875, yerr=1.96 * K_q875_se, fmt='o-',
             color='C0', capsize=4, markersize=7, label='K_eff_band(q=0.875), 95% CI')
ax1.axhline(14.5, ls='--', color='gray', alpha=0.5, label='v3.4 assumed asymptote (14.5)')
ax1.axhline(K_q875[-1], ls=':', color='C0', alpha=0.7,
            label=f'2^42 value = {K_q875[-1]:.3f}')
# Mark the two reversals
ax1.annotate(f'peak\n2^30: {K_q875[3]:.3f}', xy=(30, K_q875[3]), xytext=(28.5, 14.85),
             fontsize=8, ha='center',
             arrowprops=dict(arrowstyle='->', color='red', alpha=0.6))
ax1.annotate(f'low\n2^36: {K_q875[6]:.3f}', xy=(36, K_q875[6]), xytext=(34.5, 14.50),
             fontsize=8, ha='center',
             arrowprops=dict(arrowstyle='->', color='red', alpha=0.6))
ax1.annotate(f'2nd reversal\n2^38: {K_q875[7]:.3f}', xy=(38, K_q875[7]), xytext=(40, 14.50),
             fontsize=8, ha='center',
             arrowprops=dict(arrowstyle='->', color='green', alpha=0.6))
ax1.set_ylabel('K_eff_band(q=0.875)')
ax1.set_title('q=0.875 band reversal verification: rises, peaks at 2^30, dips to 2^36, '
              'reverses up and stabilizes ~14.63')
ax1.grid(True, alpha=0.3)
ax1.legend(loc='lower right', fontsize=9)

# Bottom: ξ_X trajectory (Weibull → Gumbel rotation)
ax2.errorbar(log2N, xi_X, yerr=1.96 * xi_X_se, fmt='s-',
             color='C2', capsize=4, markersize=7, label='ξ_X (USER frame), 95% CI')
ax2.axhline(0, ls='--', color='gray', alpha=0.5, label='ξ=0 (Gumbel limit)')
# Linear extrapolation 36→42
log2_pts = np.array([36, 38, 40, 42])
xi_pts = xi_X[-4:]
p = np.polyfit(log2_pts, xi_pts, 1)
xx = np.array([34, 50])
ax2.plot(xx, p[0] * xx + p[1], ':', color='C2', alpha=0.5,
         label=f'lin extrap: ξ=0 at log2N≈{-p[1]/p[0]:.1f}')
ax2.set_xlabel('log2 N')
ax2.set_ylabel('ξ_X (GPD shape)')
ax2.set_title('ξ_X rotation toward 0: Weibull-to-Gumbel transition; not yet complete at 2^42')
ax2.grid(True, alpha=0.3)
ax2.legend(loc='lower right', fontsize=9)

plt.tight_layout()
out_png = out_dir / "59_q875_extension.png"
plt.savefig(out_png, dpi=110, bbox_inches='tight')
print(f"[plot] {out_png}")
plt.close()
