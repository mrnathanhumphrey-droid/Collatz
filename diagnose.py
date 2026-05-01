"""
Stage 4: diagnostics, posterior summaries, GPD tail layer, tail-probability
plots for a Collatz hierarchical Stan fit.

Inputs:
  --fit-dir   directory of cmdstan output CSVs (e.g. fits/k6_uniform_N1600000_full)
  --input     stage 3 input parquet (used to compute residuals against posterior)
  --out       output dir (default: stage4_results/<fit_dir.name>)

Outputs (PNG + CSV):
  diagnose.txt              raw cmdstan diagnose() output
  summary.csv               selected param summary (means, q05/50/95, R-hat)
  fig_hyperparams.png       hyperparameter posteriors
  fig_class_forest.png      per-class alpha & beta forest plots
  gpd_per_class.csv         per-class GPD fit (xi, scale, threshold)
  fig_gpd_per_class.png     histograms of GPD shape and scale across classes
  fig_tail_probabilities.png posterior intervals on P(sigma > c*ln(n) | class)
"""

import argparse
import sys
import warnings
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from cmdstanpy import from_csv
from scipy.stats import genpareto, norm


# Heuristic for ODD-n filter
LOG_FACTOR_ODD = 3.0 / (np.log(4.0) - np.log(3.0))


def diagnostics_summary(fit, out_dir: Path):
    diag = fit.diagnose()
    with open(out_dir / "diagnose.txt", "w", encoding="utf-8") as f:
        f.write(diag)
    print("[diag] saved diagnose.txt")

    # quick R-hat / ESS scan via stan_variables_xr-equivalent: use summary()
    try:
        s = fit.summary(percentiles=(5, 50, 95))
        s.to_csv(out_dir / "summary.csv")
        rhat = s["R_hat"].dropna()
        ess = s["ESS_bulk"].dropna()
        bad_rhat = (rhat > 1.01).sum()
        bad_ess = (ess < 400).sum()
        print(f"[diag] params: {len(rhat)}  | "
              f"R-hat > 1.01: {bad_rhat}  | ESS_bulk < 400: {bad_ess}")
        print(f"[diag] R-hat range: [{rhat.min():.4f}, {rhat.max():.4f}]")
    except Exception as e:
        print(f"[diag] summary() failed: {e}")


def hyperparam_summary(fit, out_dir: Path):
    params = ["mu_alpha", "mu_beta", "tau_alpha", "tau_beta", "phi"]
    print()
    print(f"  {'param':>10}  {'mean':>10}  {'std':>10}  {'q05':>10}  {'q50':>10}  {'q95':>10}")
    rows = []
    for p in params:
        d = fit.stan_variable(p).flatten()
        rows.append((p, d.mean(), d.std(),
                     np.quantile(d, 0.05), np.quantile(d, 0.50), np.quantile(d, 0.95)))
        print(f"  {p:>10}  {d.mean():>10.4f}  {d.std():>10.4f}  "
              f"{np.quantile(d, 0.05):>10.4f}  {np.quantile(d, 0.50):>10.4f}  "
              f"{np.quantile(d, 0.95):>10.4f}")

    fig, axes = plt.subplots(1, 5, figsize=(20, 4))
    for ax, p in zip(axes, params):
        d = fit.stan_variable(p).flatten()
        ax.hist(d, bins=50, color="C0", alpha=0.75)
        ax.axvline(d.mean(), color="red", linestyle="--", alpha=0.8,
                   label=f"mean={d.mean():.4g}")
        ax.set_xlabel(p)
        ax.set_title(p)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig.suptitle("Global hyperparameter posteriors", y=1.02)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_hyperparams.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    print("[hyper] saved fig_hyperparams.png")


def class_forest_plots(fit, out_dir: Path):
    alpha = fit.stan_variable("alpha")  # (draws, K)
    beta = fit.stan_variable("beta")
    K = alpha.shape[1]

    alpha_q = np.quantile(alpha, [0.025, 0.5, 0.975], axis=0)
    beta_q = np.quantile(beta, [0.025, 0.5, 0.975], axis=0)
    mu_alpha = fit.stan_variable("mu_alpha").mean()
    mu_beta = fit.stan_variable("mu_beta").mean()

    label_fontsize = max(4, min(10, 200 // K))
    height = max(6, K * 0.18)

    fig, axes = plt.subplots(1, 2, figsize=(14, height))

    for ax, q, mu, name in zip(axes, [alpha_q, beta_q], [mu_alpha, mu_beta],
                               ["alpha", "beta"]):
        order = np.argsort(q[1])
        ys = np.arange(K)
        ax.errorbar(q[1, order], ys,
                    xerr=[q[1, order] - q[0, order], q[2, order] - q[1, order]],
                    fmt="o", markersize=3, lw=0.8, color="C0")
        ax.set_yticks(ys)
        # show odd residue label: class_idx k -> 2k+1 mod 2^(k_power)
        ax.set_yticklabels([f"k={k} (n≡{2*k+1})" for k in order],
                           fontsize=label_fontsize)
        ax.axvline(mu, color="red", linestyle="--", alpha=0.6,
                   label=f"mu_{name}={mu:.3f}")
        ax.set_xlabel(f"{name}  (95% credible interval, sorted by median)")
        ax.set_title(f"per-class {name}")
        ax.legend(loc="best", fontsize=9)
        ax.grid(True, alpha=0.3, axis="x")

    fig.suptitle(f"Per-class hierarchical posteriors  (K={K} odd residue classes)",
                 y=1.005)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_class_forest.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"[forest] saved fig_class_forest.png  (K={K})")

    # also print top/bottom 5 for each
    print()
    print("[forest] beta extremes:")
    order_b = np.argsort(beta_q[1])
    for k in order_b[:5]:
        print(f"  shallowest k={k:>4} (n≡{2*k+1:>5} mod {2*K}):  "
              f"beta = [{beta_q[0,k]:.3f}, {beta_q[1,k]:.3f}, {beta_q[2,k]:.3f}]")
    for k in order_b[-5:][::-1]:
        print(f"  steepest   k={k:>4} (n≡{2*k+1:>5} mod {2*K}):  "
              f"beta = [{beta_q[0,k]:.3f}, {beta_q[1,k]:.3f}, {beta_q[2,k]:.3f}]")


def gpd_residual_fit(fit, data_path: Path, out_dir: Path, threshold_q: float = 0.95):
    """Returns (xi, scale, threshold, n_excess) per class. NaN where fit failed."""
    df = pl.read_parquet(data_path)
    log_n = df["log_n"].to_numpy()
    sigma = df["sigma"].to_numpy().astype(np.float64)
    class_idx = df["class_idx"].to_numpy()
    K = int(class_idx.max()) + 1

    alpha_mean = fit.stan_variable("alpha").mean(axis=0)
    beta_mean = fit.stan_variable("beta").mean(axis=0)

    pred = alpha_mean[class_idx] + beta_mean[class_idx] * log_n
    resid = sigma - pred

    xi = np.full(K, np.nan)
    scale = np.full(K, np.nan)
    threshold = np.full(K, np.nan)
    n_excess = np.zeros(K, dtype=int)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for k in range(K):
            m = class_idx == k
            rk = resid[m]
            thr = np.quantile(rk, threshold_q)
            excess = rk[rk > thr] - thr
            n_excess[k] = len(excess)
            if len(excess) < 30:
                continue
            try:
                # genpareto.fit with floc=0 fits (xi, scale)
                params = genpareto.fit(excess, floc=0)
                xi[k] = params[0]
                scale[k] = params[2]
                threshold[k] = thr
            except Exception:
                pass

    out_csv = out_dir / "gpd_per_class.csv"
    pl.DataFrame({
        "class_idx": np.arange(K),
        "n_excess": n_excess,
        "threshold": threshold,
        "xi": xi,
        "scale": scale,
    }).write_csv(out_csv)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    valid_xi = xi[~np.isnan(xi)]
    valid_scale = scale[~np.isnan(scale)]

    axes[0].hist(valid_xi, bins=20, color="C2", alpha=0.75, edgecolor="black")
    axes[0].axvline(0, color="k", linestyle="--", alpha=0.7,
                    label="xi=0 (exponential)")
    axes[0].axvline(valid_xi.mean(), color="red", linestyle=":", alpha=0.8,
                    label=f"mean={valid_xi.mean():.3f}")
    axes[0].set_xlabel("GPD shape parameter xi")
    axes[0].set_ylabel("# classes")
    axes[0].set_title(f"Per-class GPD shape (top {(1-threshold_q)*100:.0f}% of residuals)\n"
                      f"xi<0: bounded, xi=0: exponential, xi>0: power-law")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].hist(valid_scale, bins=20, color="C3", alpha=0.75, edgecolor="black")
    axes[1].axvline(valid_scale.mean(), color="red", linestyle=":", alpha=0.8,
                    label=f"mean={valid_scale.mean():.2f}")
    axes[1].set_xlabel("GPD scale parameter sigma")
    axes[1].set_ylabel("# classes")
    axes[1].set_title("Per-class GPD scale")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / "fig_gpd_per_class.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

    print()
    print(f"[gpd]   xi summary across {len(valid_xi)} classes: "
          f"mean={valid_xi.mean():.3f}, std={valid_xi.std():.3f}, "
          f"range=[{valid_xi.min():.3f}, {valid_xi.max():.3f}]")
    n_sub = (xi < -0.05).sum()
    n_exp = ((xi >= -0.05) & (xi <= 0.05)).sum()
    n_sup = (xi > 0.05).sum()
    print(f"[gpd]   sub-exp (xi<-.05): {n_sub}  |  exp (|xi|<=.05): {n_exp}  "
          f"|  super-exp (xi>.05): {n_sup}")
    print(f"[gpd]   saved gpd_per_class.csv + fig_gpd_per_class.png")
    return xi, scale, threshold, n_excess


def tail_probability_plots(fit, gpd_xi, gpd_scale, gpd_threshold, out_dir: Path,
                           threshold_q: float = 0.95):
    """
    Hybrid tail-probability sweep: Normal bulk below the per-class GPD threshold,
    GPD CCDF above it. GPD parameters fixed at posterior-mean (point estimate);
    Normal-bulk uncertainty propagates from posterior of alpha, beta, phi.

    For class k, residual threshold = c*ln_n - alpha_k - beta_k*ln_n. If this
    is below the per-class u_k threshold, use Normal CDF on the residual.
    Otherwise, P = (1 - threshold_q) * GPD_CCDF(residual - u_k; xi_k, scale_k).
    """
    alpha = fit.stan_variable("alpha")  # (draws, K)
    beta = fit.stan_variable("beta")
    phi = fit.stan_variable("phi").flatten()  # (draws,)
    K = alpha.shape[1]
    excess_prob = 1.0 - threshold_q  # nominal P(r > u) = 0.05

    # c grid chosen to span bulk -> tail regime. Residual = c*ln_n - alpha - beta*ln_n;
    # GPD threshold u ~ 1.65*phi ~ 104. Crossover happens around c - mu_beta ~ u/ln_n.
    # At ln_n=17 with phi=63: c >= 10.45 + 100/17 ~ 16.4 enters GPD regime.
    c_grid = np.array([11.0, 13.0, 15.0, 17.0, 20.0])
    ln_n_eval = 17.0

    # broadcast GPD params to (1, K)
    u_b = gpd_threshold[None, :]
    xi_b = gpd_xi[None, :]
    sc_b = gpd_scale[None, :]
    no_gpd = np.isnan(gpd_threshold)  # classes with too few exceedances

    fig, axes = plt.subplots(len(c_grid), 1, figsize=(14, 2.5 * len(c_grid)),
                             sharex=True)

    for ax, c in zip(axes, c_grid):
        # residual threshold per (draw, class)
        r_thr = c * ln_n_eval - alpha - beta * ln_n_eval  # (draws, K)

        # bulk: P(r > r_thr) = 1 - Phi(r_thr / phi)
        p_bulk = 1.0 - norm.cdf(r_thr / phi[:, None])

        # GPD tail: only valid where r_thr > u and GPD fitted
        with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
            excess = r_thr - u_b
            base = 1.0 + xi_b * excess / sc_b
            base = np.where(base > 0, base, 0.0)
            # ξ ≈ 0 falls back to exponential
            gpd_ccdf = np.where(
                np.abs(xi_b) < 1e-6,
                np.exp(-excess / sc_b),
                base ** (-1.0 / np.where(np.abs(xi_b) < 1e-6, 1.0, xi_b)),
            )
        p_tail = excess_prob * gpd_ccdf

        # choose hybrid: Normal where r_thr <= u OR no GPD; else GPD-scaled
        use_gpd = (~no_gpd[None, :]) & (r_thr > u_b)
        p_hybrid = np.where(use_gpd, p_tail, p_bulk)
        # numerical safety
        p_hybrid = np.clip(p_hybrid, 1e-30, 1.0)

        q = np.quantile(p_hybrid, [0.025, 0.5, 0.975], axis=0)
        order = np.argsort(q[1])
        xs = np.arange(K)

        ax.errorbar(xs, q[1, order],
                    yerr=[q[1, order] - q[0, order], q[2, order] - q[1, order]],
                    fmt="o", markersize=3, lw=0.7, color="C0",
                    label="Normal+GPD hybrid")

        # also plot Normal-only as a faded reference for comparison
        q_bulk = np.quantile(p_bulk, [0.5], axis=0)[0]
        ax.plot(xs, q_bulk[order], "x", color="gray", alpha=0.5, markersize=4,
                label="Normal-only (for reference)")

        ax.set_yscale("log" if q[1].min() > 0 else "linear")
        ax.set_ylabel(f"P(sigma > {c} * ln n)")
        ax.set_title(f"c = {c},  ln(n) = {ln_n_eval}  (n ~ {np.exp(ln_n_eval):.1e}),  "
                     f"classes sorted by hybrid median")
        ax.legend(loc="best", fontsize=8)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("class rank (left = lowest tail probability)")
    fig.suptitle("Posterior tail-probability sweep — Normal bulk + GPD tail hybrid",
                 y=1.005)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_tail_probabilities.png", dpi=120,
                bbox_inches="tight")
    plt.close(fig)
    print("[tail]  saved fig_tail_probabilities.png  (Normal+GPD hybrid)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit-dir", required=True,
                    help="cmdstan output dir (contains *.csv chain files)")
    ap.add_argument("--input", required=True,
                    help="Stage-3 input parquet (used for residual GPD)")
    ap.add_argument("--out", default=None,
                    help="Output dir. Default: stage4_results/<fit_dir.name>")
    ap.add_argument("--threshold-q", type=float, default=0.95,
                    help="Per-class threshold quantile for GPD fit (default 0.95)")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    fit_dir = Path(args.fit_dir)
    out_dir = (Path(args.out) if args.out
               else here / "stage4_results" / fit_dir.name)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load]  fit from {fit_dir}")
    fit = from_csv(str(fit_dir))
    K = fit.stan_variable("alpha").shape[1]
    n_draws = fit.stan_variable("mu_beta").shape[0]
    print(f"[load]  chains={fit.chains}, draws={n_draws}, K={K}")
    print(f"[out]   {out_dir}")
    print()

    print("=== Convergence diagnostics ===")
    diagnostics_summary(fit, out_dir)

    print()
    print("=== Hyperparameters ===")
    hyperparam_summary(fit, out_dir)

    print()
    print("=== Per-class forest plots ===")
    class_forest_plots(fit, out_dir)

    print()
    print("=== Post-hoc GPD on residuals ===")
    gpd_xi, gpd_scale, gpd_threshold, _ = gpd_residual_fit(
        fit, Path(args.input), out_dir, threshold_q=args.threshold_q)

    print()
    print("=== Posterior tail probabilities (Normal+GPD hybrid) ===")
    tail_probability_plots(fit, gpd_xi, gpd_scale, gpd_threshold, out_dir,
                           threshold_q=args.threshold_q)

    print()
    print(f"[done]  outputs -> {out_dir}")


if __name__ == "__main__":
    main()
