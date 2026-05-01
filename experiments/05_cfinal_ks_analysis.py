"""
Experiment 05 — c_final structure: pairwise KS within and across a_final groups

For odd residue classes mod 64 grouped by their prefix terminal a_final value,
test whether classes sharing a_final but differing in c_final have
distinguishable per-class residual distributions (eps = sigma - alpha_class -
beta_class * ln n via per-class OLS).

Within-cluster pairwise KS detects c_final modulation.
Cross-cluster pairwise KS detects a_final clustering.

Also sweeps across N to test whether the within-cluster effect shrinks with
N (would imply asymptotic "exactly k distributions") or persists at fixed
ratio to the cross-cluster effect.

Usage:
    python 05_cfinal_ks_analysis.py --Ns 8388608 16777216 33554432 134217728
"""
import argparse
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

sys.stdout.reconfigure(encoding="utf-8")


def deterministic_prefix(r, a0=64):
    a, c = a0, r; steps = 0
    while a % 2 == 0:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3*c + 1
        steps += 1
    return steps, a, c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Ns", type=int, nargs="+",
                    default=[1<<23, 1<<24, 1<<25, 1<<27])
    ap.add_argument("--n_ks", type=int, default=50_000,
                    help="Subsample size per class for KS test (default 50000)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--data", type=str, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data) if args.data else here.parent / "data"
    out_dir = Path(args.out) if args.out else here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    # k=6 a_final grouping
    by_a = defaultdict(list)
    for k in range(32):
        r = 2*k + 1
        steps, a_f, c_f = deterministic_prefix(r)
        by_a[a_f].append((k, r, c_f))
    multi_a = [a for a in sorted(by_a) if len(by_a[a]) > 1]
    print("Class membership by a_final:")
    for a_f in sorted(by_a):
        members = by_a[a_f]
        print(f"  a_final={a_f:>3}: {len(members):>2} classes  "
              + ", ".join(f"r={r}(c={c})" for _, r, c in members))
    print()

    all_results = {}
    for N in args.Ns:
        path = data_dir / f"main_N{N}.parquet"
        if not path.exists():
            print(f"[skip] {path}")
            continue
        df = pl.read_parquet(path).filter(
            (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
        n = df["n"].to_numpy().astype(np.int64)
        log_n = np.log(n.astype(np.float64))
        sigma = df["sigma"].to_numpy().astype(np.float64)
        res64 = (n % 64).astype(np.int32)
        class_idx = ((res64 - 1) // 2).astype(np.int32)

        # Per-class residuals
        class_eps = {}
        for k in range(32):
            m = class_idx == k
            bk, ak = np.polyfit(log_n[m], sigma[m], 1)
            class_eps[k] = sigma[m] - (ak + bk * log_n[m])

        # Subsample
        rng = np.random.default_rng(args.seed)
        eps_sub = {}
        for k, eps in class_eps.items():
            if len(eps) >= args.n_ks:
                idx = rng.choice(len(eps), args.n_ks, replace=False)
                eps_sub[k] = eps[idx]
            else:
                eps_sub[k] = eps

        # Within
        within = []
        for a_f in multi_a:
            members = by_a[a_f]
            for i in range(len(members)):
                for j in range(i+1, len(members)):
                    k_i = members[i][0]; k_j = members[j][0]
                    stat, p = ks_2samp(eps_sub[k_i], eps_sub[k_j])
                    within.append((a_f, stat, p))
        # Cross (one rep per a_final group)
        cross = []
        a_list = list(multi_a)
        for i in range(len(a_list)):
            for j in range(i+1, len(a_list)):
                k_i = by_a[a_list[i]][0][0]
                k_j = by_a[a_list[j]][0][0]
                stat, p = ks_2samp(eps_sub[k_i], eps_sub[k_j])
                cross.append((a_list[i], a_list[j], stat, p))

        within_ks = np.array([w[1] for w in within])
        cross_ks = np.array([c[2] for c in cross])
        bonf_alpha = 0.05 / len(within)
        n_reject = sum(1 for w in within if w[2] < bonf_alpha)

        print(f"N={N:>13,}: within-KS median={np.median(within_ks):.4f} mean={within_ks.mean():.4f}, "
              f"cross-KS median={np.median(cross_ks):.4f} mean={cross_ks.mean():.4f}, "
              f"ratio={cross_ks.mean()/within_ks.mean():.2f}x, rejects={n_reject}/{len(within)}")
        all_results[N] = dict(within=within_ks, cross=cross_ks, n_reject=n_reject,
                              n_within=len(within))

    # Plot
    Ns_used = sorted(all_results.keys())
    log_Ns = np.log2(Ns_used)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.plot(log_Ns, [np.median(all_results[N]["within"]) for N in Ns_used], "o-",
            label="within-a_final (median)", color="C0")
    ax.plot(log_Ns, [all_results[N]["within"].mean() for N in Ns_used], "o--",
            label="within-a_final (mean)", color="C0", alpha=0.5)
    ax.plot(log_Ns, [np.median(all_results[N]["cross"]) for N in Ns_used], "s-",
            label="cross-a_final (median)", color="C3")
    ax.plot(log_Ns, [all_results[N]["cross"].mean() for N in Ns_used], "s--",
            label="cross-a_final (mean)", color="C3", alpha=0.5)
    ax.set_xlabel("log2(N)"); ax.set_ylabel("KS_stat"); ax.set_yscale("log")
    ax.set_title("KS effect size vs N")
    ax.legend(); ax.grid(True, alpha=0.3, which="both")

    ax = axes[1]
    n_rej = [all_results[N]["n_reject"] for N in Ns_used]
    n_within = all_results[Ns_used[0]]["n_within"]
    ax.plot(log_Ns, n_rej, "o-", color="C2", markersize=8)
    ax.set_xlabel("log2(N)"); ax.set_ylabel("within-group rejections at Bonferroni")
    ax.set_title(f"Rejection count out of {n_within} within-group pairs")
    ax.grid(True, alpha=0.3)

    fig.suptitle("Does c_final effect shrink with N?", y=1.02)
    plt.tight_layout()
    out_png = out_dir / "05_cfinal_ks_scaling.png"
    plt.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close()
    print()
    print(f"[save] {out_png}")

    if len(Ns_used) >= 2:
        log_w = np.log([np.median(all_results[N]["within"]) for N in Ns_used])
        b_w, _ = np.polyfit(np.log(Ns_used), log_w, 1)
        log_c = np.log([np.median(all_results[N]["cross"]) for N in Ns_used])
        b_c, _ = np.polyfit(np.log(Ns_used), log_c, 1)
        print(f"within-median KS scaling exponent: {b_w:.3f}")
        print(f"cross-median  KS scaling exponent: {b_c:.3f}")


if __name__ == "__main__":
    main()
