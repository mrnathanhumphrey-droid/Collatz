"""
Experiment 03 — N-scaling of universality

For each N in a sweep, on full odd-only data:
  - Pooled OLS slope mu_beta (does it converge to heuristic 10.4282?)
  - Per-class slope SD (naive tau_beta)
  - Per-class sampling SE on slope (noise floor)
  - Moment-corrected tau_beta = sqrt(max(0, naive^2 - SE^2))
  - Mean GPD shape xi from per-class top-5% residual fits

Tracks asymptotic universality: slope (mu_beta -> 10.4282), tail shape
(xi -> 0), between-class slope variation (tau_beta at noise floor).

Usage:
    python 03_n_scaling.py --Ns 1048576 4194304 8388608 16777216 33554432 134217728
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from scipy.stats import genpareto

sys.stdout.reconfigure(encoding="utf-8")

LOG_FACTOR_ODD = 3.0 / (np.log(4.0) - np.log(3.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Ns", type=int, nargs="+",
                    default=[1<<20, 1<<22, 1<<23, 1<<24, 1<<25, 1<<27])
    ap.add_argument("--k", type=int, default=6, help="Modular resolution (default 6)")
    ap.add_argument("--data", type=str, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data) if args.data else here.parent / "data"
    out_dir = Path(args.out) if args.out else here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    M = 1 << args.k; K = M // 2
    print(f"[config] k={args.k}, K={K} odd classes")
    print()

    rows = []
    for N in args.Ns:
        path = data_dir / f"main_N{N}.parquet"
        if not path.exists():
            print(f"[skip] {path} (run generate.py --N {N} first)")
            continue
        df = pl.read_parquet(path).filter(
            (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
        n = df["n"].to_numpy().astype(np.int64)
        log_n = np.log(n.astype(np.float64))
        sigma = df["sigma"].to_numpy().astype(np.float64)
        res = (n % M).astype(np.int32)
        class_idx = ((res - 1) // 2).astype(np.int32)

        # Pooled OLS
        b_p, a_p = np.polyfit(log_n, sigma, 1)
        phi = (sigma - (a_p + b_p * log_n)).std()

        # Per-class
        betas = np.zeros(K); alphas = np.zeros(K); ses = np.zeros(K); xis = np.zeros(K)
        for kk in range(K):
            m = class_idx == kk
            if m.sum() < 50: continue
            bk, ak = np.polyfit(log_n[m], sigma[m], 1)
            rk = sigma[m] - (ak + bk * log_n[m])
            betas[kk] = bk; alphas[kk] = ak
            var_x = log_n[m].var()
            ses[kk] = rk.std() / np.sqrt(m.sum() * var_x)
            # GPD on top 5%
            thr = np.quantile(rk, 0.95)
            excess = rk[rk > thr] - thr
            if len(excess) >= 30:
                try:
                    p_gpd = genpareto.fit(excess, floc=0)
                    xis[kk] = p_gpd[0]
                except Exception:
                    xis[kk] = np.nan

        sd_betas = betas.std()
        mean_se = ses.mean()
        tau_beta_corr = np.sqrt(max(0, sd_betas**2 - mean_se**2))

        rows.append(dict(
            N=N, n_odd=int(df.height), per_class=int(df.height // K),
            mu_beta=b_p, phi=phi, sd_class_beta=sd_betas, mean_se_beta=mean_se,
            tau_beta_corr=tau_beta_corr,
            sd_class_alpha=alphas.std(),
            mean_xi=np.nanmean(xis), std_xi=np.nanstd(xis)))

    print(f"{'N':>13} {'odd':>13} {'per/cls':>10} {'mu_beta':>9} {'phi':>7} "
          f"{'SD(beta)':>9} {'SE(beta)':>9} {'tau_corr':>9} "
          f"{'SD(alpha)':>10} {'<xi>':>8}")
    for r in rows:
        print(f"{r['N']:>13,} {r['n_odd']:>13,} {r['per_class']:>10,} "
              f"{r['mu_beta']:>9.4f} {r['phi']:>7.2f} "
              f"{r['sd_class_beta']:>9.4f} {r['mean_se_beta']:>9.4f} {r['tau_beta_corr']:>9.4f} "
              f"{r['sd_class_alpha']:>10.3f} {r['mean_xi']:>+8.4f}")

    # Plot
    Ns = np.array([r["N"] for r in rows])
    log_Ns = np.log2(Ns)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    ax = axes[0, 0]
    ax.semilogy(log_Ns, [r["sd_class_beta"] for r in rows], "o-",
                label="SD(class slopes)")
    ax.semilogy(log_Ns, [r["mean_se_beta"] for r in rows], "^-",
                label="per-class SE", alpha=0.6)
    ax.set_xlabel("log2(N)"); ax.set_ylabel("tau_beta")
    ax.set_title(f"tau_beta vs N at k={args.k}\n(SD tracks SE = at noise floor)")
    ax.legend(); ax.grid(True, alpha=0.3, which="both")

    ax = axes[0, 1]
    ax.plot(log_Ns, [r["mu_beta"] for r in rows], "o-")
    ax.axhline(LOG_FACTOR_ODD, color="red", linestyle="--", alpha=0.6,
               label=f"heuristic = {LOG_FACTOR_ODD:.4f}")
    ax.set_xlabel("log2(N)"); ax.set_ylabel("mu_beta (pooled OLS slope)")
    ax.set_title("mu_beta vs N")
    ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(log_Ns, [r["phi"] for r in rows], "o-", color="C3")
    ax.set_xlabel("log2(N)"); ax.set_ylabel("phi")
    ax.set_title("Pooled residual SD vs N")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.errorbar(log_Ns, [r["mean_xi"] for r in rows],
                yerr=[r["std_xi"] for r in rows], fmt="o-", color="C2")
    ax.axhline(0, color="black", linestyle=":", alpha=0.5,
               label="exponential (xi=0)")
    ax.set_xlabel("log2(N)"); ax.set_ylabel("mean GPD shape xi")
    ax.set_title("xi vs N (across classes; mean +/- SD)")
    ax.legend(); ax.grid(True, alpha=0.3)

    fig.suptitle(f"N-scaling of asymptotic universality (k={args.k})", y=1.01)
    plt.tight_layout()
    out_png = out_dir / f"03_n_scaling_k{args.k}.png"
    plt.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close()
    print()
    print(f"[save] {out_png}")


if __name__ == "__main__":
    main()
