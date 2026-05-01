"""
Stage 2: Exploratory analysis. Reads Stage 1 parquets, produces 4 figures.

Figures (saved to figures/):
  1. fig1_sigma_vs_logn_by_mod16.png
       4x4 hexbin grid, one panel per residue mod 16, with heuristic
       reference line E[sigma|n] = (2/(log4-log3)) * log(n).
  2. fig2_v_distribution.png
       Empirical P(v=k) for k=1..max vs Geometric(1/2) prediction
       2^-k, semi-log y. Shows where Syracuse v deviates from heuristic.
  3. fig3_residual_tail.png
       CCDF of residual r(n) = sigma(n) - heuristic on (a) semi-log and
       (b) log-log. Linear in (a) => exponential tail; linear in (b) =>
       power-law tail. Reference lines fit by least squares to the
       upper-tail decade.
  4. fig4_residual_tail_by_mod64.png
       8x8 grid of per-class CCDFs (mod 64), log-log. Shows whether the
       heaviest tails live in particular residue classes.
"""

import argparse
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib as mpl


# Heuristic: E[sigma|n] = (2 / (log 4 - log 3)) * log(n)  ~ 6.9521 * log(n)
LOG_FACTOR = 2.0 / (np.log(4.0) - np.log(3.0))


def load(data_dir: Path, N: int):
    main = pl.read_parquet(data_dir / f"main_N{N}.parquet")
    vseq = pl.read_parquet(data_dir / f"v_seq_N{N}.parquet")
    return main, vseq


def fig1_sigma_vs_logn_by_mod16(df: pl.DataFrame, out: Path):
    n = df["n"].to_numpy()
    sigma = df["sigma"].to_numpy()
    res16 = df["res_mod_16"].to_numpy()
    log_n = np.log(n.astype(np.float64))

    fig, axes = plt.subplots(4, 4, figsize=(16, 12), sharex=True, sharey=True)
    log_max = log_n.max()
    sig_max = sigma.max()
    ref_x = np.linspace(0.0, log_max, 200)
    ref_y = LOG_FACTOR * ref_x

    for k in range(16):
        ax = axes[k // 4, k % 4]
        m = res16 == k
        if m.sum() > 0:
            ax.hexbin(
                log_n[m], sigma[m],
                gridsize=80, cmap="viridis", mincnt=1, bins="log",
                extent=(0, log_max, 0, sig_max),
            )
        ax.plot(ref_x, ref_y, color="red", lw=1.0, alpha=0.7,
                label="heuristic" if k == 0 else None)
        ax.set_title(f"n ≡ {k} (mod 16)   [count={int(m.sum()):,}]",
                     fontsize=9)
        ax.set_xlim(0, log_max)
        ax.set_ylim(0, sig_max)

    fig.supxlabel("log(n)")
    fig.supylabel("σ(n)")
    fig.suptitle(f"σ(n) vs log(n) stratified by residue class mod 16  "
                 f"(N = {len(df):,})", y=0.995)
    axes[0, 0].legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)


def fig2_v_distribution(vdf: pl.DataFrame, out: Path):
    counts = (
        vdf.group_by("v")
        .agg(pl.len().alias("count"))
        .sort("v")
    )
    v_vals = counts["v"].to_numpy().astype(np.int32)
    cnts = counts["count"].to_numpy().astype(np.float64)
    total = cnts.sum()
    emp = cnts / total

    th = np.where(v_vals >= 1, 0.5 ** v_vals, 0.0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    width = 0.4
    ax1.bar(v_vals - width / 2, emp, width=width, label="empirical",
            color="C0", alpha=0.85)
    ax1.bar(v_vals + width / 2, th, width=width, label="Geometric(1/2)",
            color="C3", alpha=0.85)
    ax1.set_yscale("log")
    ax1.set_xlabel("v = ν₂(3m+1)")
    ax1.set_ylabel("P(v)")
    ax1.set_title("Syracuse v-distribution: empirical vs Geometric(1/2)")
    ax1.legend()
    ax1.grid(True, alpha=0.3, which="both")

    ratio = np.where(th > 0, emp / th, np.nan)
    ax2.axhline(1.0, color="k", lw=0.8, alpha=0.5)
    ax2.bar(v_vals, ratio, color="C2", alpha=0.85)
    ax2.set_xlabel("v")
    ax2.set_ylabel("empirical / theoretical")
    ax2.set_title("Ratio (>1 = empirical heavier than geometric)")
    ax2.set_ylim(0, max(2.0, np.nanmax(ratio) * 1.1))
    ax2.grid(True, alpha=0.3)

    fig.suptitle(f"v-distribution diagnostics  (total v-steps = {int(total):,})")
    fig.tight_layout()
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)


def _ccdf(x: np.ndarray):
    """Return (sorted_values_desc, survival_prob) for the upper tail."""
    s = np.sort(x)
    n = len(s)
    # P(X > x_i) for sorted ascending = (n - i - 1) / n; we'll plot vs s
    p = np.arange(n - 1, -1, -1, dtype=np.float64) / n
    # drop the last point (p=0) for log scales
    return s, p


def _fit_powerlaw_tail(t, p, frac=0.05):
    """Least-squares fit of log p ~ a + b * log t on the top-frac of t > 0."""
    mask = (t > 0) & (p > 0)
    t = t[mask]; p = p[mask]
    if len(t) < 50:
        return None
    cutoff = np.quantile(t, 1.0 - frac)
    sel = t >= cutoff
    if sel.sum() < 20:
        return None
    lt = np.log(t[sel]); lp = np.log(p[sel])
    b, a = np.polyfit(lt, lp, 1)
    return a, b, cutoff


def _fit_exponential_tail(t, p, frac=0.05):
    """Least-squares fit of log p ~ a + b * t on the top-frac of t."""
    mask = (p > 0)
    t = t[mask]; p = p[mask]
    if len(t) < 50:
        return None
    cutoff = np.quantile(t, 1.0 - frac)
    sel = t >= cutoff
    if sel.sum() < 20:
        return None
    lp = np.log(p[sel])
    b, a = np.polyfit(t[sel], lp, 1)
    return a, b, cutoff


def fig3_residual_tail(df: pl.DataFrame, out: Path):
    n = df["n"].to_numpy()
    sigma = df["sigma"].to_numpy().astype(np.float64)
    log_n = np.log(n.astype(np.float64))
    r = sigma - LOG_FACTOR * log_n

    # Right tail: positive residuals
    r_pos = r[r > 0]

    t, p = _ccdf(r_pos)
    # drop zero-survival points for log scales
    nz = p > 0
    t, p = t[nz], p[nz]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # (a) semi-log: linear => exponential
    ax1.plot(t, p, color="C0", lw=1.2, label="empirical CCDF")
    fit_e = _fit_exponential_tail(t, p, frac=0.05)
    if fit_e is not None:
        a, b, cutoff = fit_e
        xs = np.linspace(cutoff, t.max(), 100)
        ax1.plot(xs, np.exp(a + b * xs), "r--", lw=1.0,
                 label=f"exp fit (top 5%): rate={-b:.3f}")
    ax1.set_yscale("log")
    ax1.set_xlabel("t = σ(n) - heuristic")
    ax1.set_ylabel("P(residual > t)")
    ax1.set_title("Right-tail CCDF, semi-log\n(linear ⇒ exponential tail)")
    ax1.legend()
    ax1.grid(True, alpha=0.3, which="both")

    # (b) log-log: linear => power-law
    ax2.plot(t, p, color="C0", lw=1.2, label="empirical CCDF")
    fit_p = _fit_powerlaw_tail(t, p, frac=0.05)
    if fit_p is not None:
        a, b, cutoff = fit_p
        xs = np.geomspace(cutoff, t.max(), 100)
        ax2.plot(xs, np.exp(a) * xs ** b, "r--", lw=1.0,
                 label=f"power fit (top 5%): exponent={b:.2f}")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("t = σ(n) - heuristic")
    ax2.set_ylabel("P(residual > t)")
    ax2.set_title("Right-tail CCDF, log-log\n(linear ⇒ power-law tail)")
    ax2.legend()
    ax2.grid(True, alpha=0.3, which="both")

    fig.suptitle(f"Residual r(n) = σ(n) − (2/(log 4−log 3))·log(n)  "
                 f"|  positive-residual sample size = {len(r_pos):,}")
    fig.tight_layout()
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)


def fig4_residual_tail_by_mod64(df: pl.DataFrame, out: Path):
    n = df["n"].to_numpy()
    sigma = df["sigma"].to_numpy().astype(np.float64)
    res64 = df["res_mod_64"].to_numpy()
    log_n = np.log(n.astype(np.float64))
    r = sigma - LOG_FACTOR * log_n

    fig, axes = plt.subplots(8, 8, figsize=(20, 18),
                             sharex=True, sharey=True)

    # global limits for shared axes
    r_pos_all = r[r > 0]
    if len(r_pos_all) == 0:
        plt.close(fig); return
    t_max = r_pos_all.max()

    # collect tail exponents to color panels
    exponents = np.full(64, np.nan)
    panel_data = []

    for k in range(64):
        m = (res64 == k) & (r > 0)
        rk = r[m]
        if len(rk) < 30:
            panel_data.append((k, None, None, None))
            continue
        t, p = _ccdf(rk)
        nz = p > 0
        t, p = t[nz], p[nz]
        fit = _fit_powerlaw_tail(t, p, frac=0.10)
        slope = fit[1] if fit is not None else np.nan
        exponents[k] = slope
        panel_data.append((k, t, p, slope))

    # color scale based on exponent (more negative = thinner tail)
    valid = exponents[~np.isnan(exponents)]
    if len(valid) > 0:
        vmin, vmax = np.percentile(valid, [5, 95])
    else:
        vmin, vmax = -3, -1
    cmap = mpl.colormaps["coolwarm"]
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    for (k, t, p, slope) in panel_data:
        ax = axes[k // 8, k % 8]
        if t is None:
            ax.text(0.5, 0.5, f"k={k}\nn<30", ha="center", va="center",
                    transform=ax.transAxes, fontsize=8, color="gray")
            ax.set_xticks([]); ax.set_yticks([])
            continue
        color = cmap(norm(slope)) if not np.isnan(slope) else "gray"
        ax.plot(t, p, color=color, lw=1.0)
        ax.set_xscale("log"); ax.set_yscale("log")
        title = f"k={k}"
        if not np.isnan(slope):
            title += f"  α={slope:.2f}"
        ax.set_title(title, fontsize=8)
        ax.tick_params(labelsize=6)

    fig.supxlabel("t = σ(n) - heuristic  (log)")
    fig.supylabel("P(residual > t)  (log)")
    fig.suptitle(
        f"Right-tail CCDF stratified by residue mod 64  (N = {len(df):,})\n"
        f"color = power-law tail exponent (blue thinner, red thicker)",
        y=0.995, fontsize=12,
    )

    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), shrink=0.6,
                        pad=0.02, location="right")
    cbar.set_label("power-law exponent (top-10% fit)")

    fig.savefig(out, dpi=110, bbox_inches="tight")
    plt.close(fig)

    # print top-5 thickest and thinnest tails
    valid_idx = np.where(~np.isnan(exponents))[0]
    if len(valid_idx) > 0:
        order = valid_idx[np.argsort(exponents[valid_idx])]
        thinnest = order[:5]
        thickest = order[-5:][::-1]
        print(f"[fig4]   thinnest tails (mod 64): "
              f"{[f'k={k}(a={exponents[k]:.2f})' for k in thinnest]}")
        print(f"[fig4]   thickest tails (mod 64): "
              f"{[f'k={k}(a={exponents[k]:.2f})' for k in thickest]}")


def main():
    ap = argparse.ArgumentParser(description="Stage 2 Collatz EDA plots.")
    ap.add_argument("--N", type=int, default=1 << 20)
    ap.add_argument("--data", type=str, default=None,
                    help="Data dir. Default: <script_dir>/data")
    ap.add_argument("--out", type=str, default=None,
                    help="Figures dir. Default: <script_dir>/figures")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data) if args.data else here / "data"
    fig_dir = Path(args.out) if args.out else here / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load]   reading parquets from {data_dir}", flush=True)
    df, vdf = load(data_dir, args.N)
    print(f"[load]   main: {df.shape}, v-seq: {vdf.shape}", flush=True)

    print("[fig1]   sigma(n) vs log(n) by mod 16 ...", flush=True)
    fig1_sigma_vs_logn_by_mod16(df, fig_dir / "fig1_sigma_vs_logn_by_mod16.png")

    print("[fig2]   v-distribution vs Geometric(1/2) ...", flush=True)
    fig2_v_distribution(vdf, fig_dir / "fig2_v_distribution.png")

    print("[fig3]   residual tail (semi-log + log-log) ...", flush=True)
    fig3_residual_tail(df, fig_dir / "fig3_residual_tail.png")

    print("[fig4]   residual tail by mod 64 (8x8 grid) ...", flush=True)
    fig4_residual_tail_by_mod64(df, fig_dir / "fig4_residual_tail_by_mod64.png")

    print(f"[done]   figures -> {fig_dir}", flush=True)


if __name__ == "__main__":
    main()
