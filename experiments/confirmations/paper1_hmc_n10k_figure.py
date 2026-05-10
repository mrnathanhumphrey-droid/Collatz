"""Figure: Pathfinder vs HMC posterior for sigma_u in B1 and B4 at N=10^4.

Reads sigma_u draws from experiments_output/paper1_hmc_n10k/sigma_u_draws.npz
and the metadata from data/paper1/hmc_n10k_results.parquet.

Two-panel histogram (B1 left, B4 right). Pathfinder vs HMC overlaid.
Posterior medians annotated.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import polars as pl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parent
COLLATZ = HERE.parent
NPZ = COLLATZ / "experiments_output" / "paper1_hmc_n10k" / "sigma_u_draws.npz"
PARQ = COLLATZ / "data" / "paper1" / "hmc_n10k_results.parquet"
FIG = COLLATZ / "figures" / "hmc_n10k_sigma_u_comparison.png"


def main():
    data = np.load(NPZ)
    hmc_B1 = data["hmc_B1"]; pf_B1 = data["pf_B1"]
    hmc_B4 = data["hmc_B4"]; pf_B4 = data["pf_B4"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, hmc, pf, title in [
        (axes[0], hmc_B1, pf_B1, r"$B_1$: $\log(n)$ + RE on $(n \bmod 8)$"),
        (axes[1], hmc_B4, pf_B4, r"$B_4$: $B_3$ + RE on $(n \bmod 8)$"),
    ]:
        all_vals = np.concatenate([hmc, pf]) if (len(hmc) and len(pf)) else hmc
        if len(all_vals) == 0:
            continue
        upper = max(np.percentile(all_vals, 99.5), 1e-3)
        bins = np.linspace(0, upper * 1.05, 60)
        ax.hist(pf, bins=bins, alpha=0.55, label="Pathfinder",
                color="tab:orange", density=True, edgecolor="white", linewidth=0.3)
        ax.hist(hmc, bins=bins, alpha=0.55, label="HMC",
                color="tab:blue", density=True, edgecolor="white", linewidth=0.3)
        med_pf = float(np.median(pf)) if len(pf) else float("nan")
        med_hmc = float(np.median(hmc)) if len(hmc) else float("nan")
        ax.axvline(med_pf, color="tab:orange", linestyle="--", linewidth=1.5,
                   label=f"Pathfinder median = {med_pf:.4f}")
        ax.axvline(med_hmc, color="tab:blue", linestyle="--", linewidth=1.5,
                   label=f"HMC median = {med_hmc:.4f}")
        ax.set_xlabel(r"$\sigma_u$")
        ax.set_ylabel("posterior density")
        ax.set_title(title)
        ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
        ax.grid(True, alpha=0.3)

    fig.suptitle(r"Posterior on $\sigma_u$: Pathfinder vs HMC at $N = 10^4$ "
                 "(8000 train / 2000 test)", y=1.00)
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=160, bbox_inches="tight")
    print(f"[save] {FIG}")
    print(f"  HMC B1 median = {np.median(hmc_B1):.6f}, n_draws = {len(hmc_B1)}")
    print(f"  PF  B1 median = {np.median(pf_B1):.6f}, n_draws = {len(pf_B1)}")
    print(f"  HMC B4 median = {np.median(hmc_B4):.6f}, n_draws = {len(hmc_B4)}")
    print(f"  PF  B4 median = {np.median(pf_B4):.6f}, n_draws = {len(pf_B4)}")


if __name__ == "__main__":
    main()
