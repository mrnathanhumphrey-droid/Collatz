"""
Experiment 10 — qx+1 generalization of the prefix decomposition.

For q in {3, 5, 7} on the same N=10^6, k=6 setup:
  1. Per residue class r mod 2^k, compute the convergence rate (fraction of
     n in class that reach 1 under T_q) and statistics on convergent n only.
  2. Compute the q-parameterized prefix decomposition: while a even, branch
     parity forced by c, apply T_q. Stop when a becomes odd. a_final in {q^j}.
  3. Test whether a_final predicts:
        (a) per-class convergence rate (the open qx+1 question)
        (b) per-class sigma_q on convergent orbits (the prefix-decomp transfer)

Usage:
    python 10_q_decomposition.py --q 5 --N 1000000 --k 6
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import pearsonr

sys.stdout.reconfigure(encoding="utf-8")


def deterministic_prefix_q(r, q, a0, max_steps=400):
    """Symbolic prefix under T_q starting from state (a0, r). Stop when a odd."""
    a, c = a0, r
    steps = 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2
            c //= 2
        else:
            a *= q
            c = q * c + 1
        steps += 1
    return steps, a, c


def analyze_q(q, df, k):
    M = 1 << k
    K = M // 2  # number of odd residue classes

    n = df["n"].to_numpy().astype(np.int64)
    log_n = np.log(n.astype(np.float64))
    sigma_q = df["sigma_q"].to_numpy().astype(np.int64)
    converged = df["converged"].to_numpy().astype(bool)
    res = (n % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)

    # Per-class deterministic prefix
    a_final_per = np.zeros(K, dtype=np.int64)
    c_final_per = np.zeros(K, dtype=np.int64)
    prefix_per = np.zeros(K, dtype=np.int64)
    for kk in range(K):
        r_class = 2 * kk + 1
        s, a, c = deterministic_prefix_q(r_class, q, M)
        prefix_per[kk] = s
        a_final_per[kk] = a
        c_final_per[kk] = c

    # Per-class convergence rate
    n_total_per = np.zeros(K, dtype=np.int64)
    n_conv_per = np.zeros(K, dtype=np.int64)
    sigma_mean_conv = np.full(K, np.nan)
    sigma_std_conv = np.full(K, np.nan)
    beta_conv = np.full(K, np.nan)
    alpha_conv = np.full(K, np.nan)

    for kk in range(K):
        m_class = class_idx == kk
        n_total_per[kk] = m_class.sum()
        m_conv = m_class & converged
        n_conv_per[kk] = m_conv.sum()
        if n_conv_per[kk] >= 50:
            sigma_c = sigma_q[m_conv].astype(np.float64)
            log_n_c = log_n[m_conv]
            sigma_mean_conv[kk] = sigma_c.mean()
            sigma_std_conv[kk] = sigma_c.std(ddof=1)
            try:
                bk, ak = np.polyfit(log_n_c, sigma_c, 1)
                beta_conv[kk] = bk
                alpha_conv[kk] = ak
            except Exception:
                pass

    conv_rate = n_conv_per / np.maximum(n_total_per, 1)

    # Does a_final predict convergence rate?
    log_a_final = np.log(a_final_per.astype(np.float64))
    if (conv_rate > 0).sum() >= 3:
        try:
            r_a, p_a = pearsonr(log_a_final, conv_rate)
            r_pre, p_pre = pearsonr(prefix_per, conv_rate)
        except Exception:
            r_a = p_a = r_pre = p_pre = np.nan
    else:
        r_a = p_a = r_pre = p_pre = np.nan

    return dict(
        q=q,
        K=K,
        n_total=int(n_total_per.sum()),
        n_conv=int(n_conv_per.sum()),
        conv_rate_overall=float(n_conv_per.sum() / max(n_total_per.sum(), 1)),
        a_final_per=a_final_per,
        prefix_per=prefix_per,
        conv_rate_per=conv_rate,
        n_conv_per=n_conv_per,
        log_a_final=log_a_final,
        beta_conv=beta_conv,
        alpha_conv=alpha_conv,
        sigma_mean_conv=sigma_mean_conv,
        r_logafinal_convrate=r_a,
        p_logafinal_convrate=p_a,
        r_prefix_convrate=r_pre,
        p_prefix_convrate=p_pre,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=5)
    ap.add_argument("--N", type=int, default=1_000_000)
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--data", type=str, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data) if args.data else here.parent / "data"
    out_dir = Path(args.out) if args.out else here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] q={args.q}, N={args.N:,}, k={args.k}", flush=True)
    df = pl.read_parquet(data_dir / f"q_main_q{args.q}_N{args.N}.parquet")
    print(f"        rows={len(df):,}, converged={df['converged'].sum():,}", flush=True)

    r = analyze_q(args.q, df, args.k)

    print()
    print(f"=== qx+1 prefix decomposition: q={args.q}, k={args.k}, N={args.N:,} ===")
    print(f"  Total odd n: {r['n_total']:,}")
    print(f"  Convergent (reached 1): {r['n_conv']:,} ({100*r['conv_rate_overall']:.4f}%)")
    print()
    print(f"  Per-class table (top 10 by convergence rate):")
    print(f"  {'class':>5} {'n_tot':>8} {'n_conv':>7} {'rate':>8} {'a_final':>8} "
          f"{'prefix':>7} {'mean(sig)':>10} {'beta':>8}")
    order = np.argsort(r["conv_rate_per"])[::-1]
    for kk in order[:10]:
        sm = r["sigma_mean_conv"][kk]
        b = r["beta_conv"][kk]
        sm_s = f"{sm:.1f}" if not np.isnan(sm) else "  --"
        b_s = f"{b:.3f}" if not np.isnan(b) else "  --"
        print(f"  {kk:>5} {int(r['n_conv_per'][kk]) + (int(r['n_conv_per'][kk]) == 0):>8} "
              f"{int(r['n_conv_per'][kk]):>7} {r['conv_rate_per'][kk]:>8.5f} "
              f"{int(r['a_final_per'][kk]):>8} {int(r['prefix_per'][kk]):>7} "
              f"{sm_s:>10} {b_s:>8}")

    print()
    print(f"  Per-class summary stats:")
    print(f"    convergence rate range: [{r['conv_rate_per'].min():.5f}, {r['conv_rate_per'].max():.5f}]")
    n_classes_with_any_conv = (r['n_conv_per'] > 0).sum()
    print(f"    classes with >=1 convergent n: {n_classes_with_any_conv}/{r['K']}")
    n_classes_with_50_conv = (r['n_conv_per'] >= 50).sum()
    print(f"    classes with >=50 convergent n (usable for OLS): {n_classes_with_50_conv}/{r['K']}")

    print()
    print(f"  Predictive correlations (across {r['K']} odd residue classes):")
    print(f"    Pearson r(log(a_final), conv_rate)  = {r['r_logafinal_convrate']:>7.4f}  "
          f"(p={r['p_logafinal_convrate']:.4g})")
    print(f"    Pearson r(prefix_steps, conv_rate)  = {r['r_prefix_convrate']:>7.4f}  "
          f"(p={r['p_prefix_convrate']:.4g})")

    if n_classes_with_50_conv >= 5:
        # On convergent orbits only: does a_final predict slope / intercept?
        valid = (r['n_conv_per'] >= 50)
        beta_v = r['beta_conv'][valid]
        alpha_v = r['alpha_conv'][valid]
        log_a_v = r['log_a_final'][valid]
        prefix_v = r['prefix_per'][valid]
        try:
            r_b, _ = pearsonr(log_a_v, beta_v)
            r_a_int, _ = pearsonr(log_a_v, alpha_v)
            print(f"    Pearson r(log(a_final), beta_conv)  = {r_b:>7.4f}  (n={valid.sum()} usable classes)")
            print(f"    Pearson r(log(a_final), alpha_conv) = {r_a_int:>7.4f}")
            print(f"    mu_beta over usable classes: {beta_v.mean():.4f}")
        except Exception as e:
            print(f"    correlation skipped: {e}")
    else:
        print(f"    too few classes with >=50 convergent n for OLS comparison")

    # Save CSV summary
    out_csv = out_dir / f"10_q_decomposition_q{args.q}_k{args.k}_N{args.N}.csv"
    rows = []
    for kk in range(r["K"]):
        rows.append({
            "class_idx": kk,
            "residue": 2 * kk + 1,
            "a_final": int(r["a_final_per"][kk]),
            "prefix_steps": int(r["prefix_per"][kk]),
            "n_total": int(r["n_conv_per"][kk]) if r["n_conv_per"][kk] > 0 else 0,
            "n_conv": int(r["n_conv_per"][kk]),
            "conv_rate": float(r["conv_rate_per"][kk]),
            "sigma_mean_conv": float(r["sigma_mean_conv"][kk]) if not np.isnan(r["sigma_mean_conv"][kk]) else None,
            "beta_conv": float(r["beta_conv"][kk]) if not np.isnan(r["beta_conv"][kk]) else None,
            "alpha_conv": float(r["alpha_conv"][kk]) if not np.isnan(r["alpha_conv"][kk]) else None,
        })
    pl.DataFrame(rows).write_csv(out_csv)
    print()
    print(f"[save] {out_csv}")


if __name__ == "__main__":
    main()
