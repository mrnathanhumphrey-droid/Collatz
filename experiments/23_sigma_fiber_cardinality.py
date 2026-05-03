"""
Experiment 23 — sigma-fiber cardinality (Avenue A, cryptographic hardness candidate).

Forward map: f_{N,k}(n) = (sigma(n), n mod 2^k), for odd n in [1, N].
Inverse problem: given (sigma, r), recover n.

The candidate hard problem is: knowing (sigma(n), n mod 2^k), how many odd n in
[1, N] match? Tier-1 measurement is the size of these preimage classes
("fibers"). Tier-2 measurement is within-fiber structure: even if a fiber is
large, lattice / arithmetic-progression structure makes it easy to crack.

  Tier 1: For each k in {6, 8, 10, 12, 14}, group odd n by (sigma, r mod 2^k)
          and report the distribution of fiber sizes (median, max, P(F=1),
          benchmark vs uniform N/(2^(k-1) * sigma_range)).

  Tier 2: For the largest fibers at k=10, look at the m values (m = (n-r)/2^k).
          Pairwise gap distribution, GCD of gaps, sorted-m lag-1 autocorrelation.
          Hard fibers are unstructured; lattice-friendly fibers are detectable
          here.

Usage:
    python 23_sigma_fiber_cardinality.py --N 10000000
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.stdout.reconfigure(encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10_000_000)
    ap.add_argument("--ks", type=str, default="6,8,10,12,14")
    ap.add_argument("--top_fibers", type=int, default=20,
                    help="How many top-size fibers to analyze for tier-2 gap structure.")
    ap.add_argument("--gap_k", type=int, default=10,
                    help="Modular resolution to use for tier-2 gap analysis.")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"

    print(f"[load] N={args.N:,}", flush=True)
    df = pl.read_parquet(data_dir / f"main_N{args.N}.parquet")
    df = df.filter(pl.col("n") % 2 == 1)
    n = df["n"].to_numpy().astype(np.int64)
    sigma = df["sigma"].to_numpy().astype(np.int64)
    print(f"        odd rows: {len(n):,}  sigma range: [{sigma.min()}, {sigma.max()}]",
          flush=True)
    print()

    ks = [int(k) for k in args.ks.split(",")]

    # ============================================================
    # Tier 1: fiber-size distribution at each k
    # ============================================================
    print("=" * 78)
    print("Tier 1: fiber-size distribution at each k")
    print("=" * 78)
    print()
    print(f"{'k':>3} {'mod':>7} {'n_fibers':>11} {'median':>7} "
          f"{'mean':>8} {'max':>7} {'P(F=1)':>8} {'P(F>1)':>8} "
          f"{'uniform_E':>11}")

    for k in ks:
        M = 1 << k
        r = (n % M).astype(np.int64)
        # Encode (sigma, r) as a single key for fast grouping.
        # sigma fits in <2^11, r fits in <2^k. Use (sigma << k) | r.
        key = (sigma << k) | r
        # Count duplicates
        _, counts = np.unique(key, return_counts=True)
        n_fibers = len(counts)
        med = int(np.median(counts))
        mn = float(counts.mean())
        mx = int(counts.max())
        p_unit = float((counts == 1).sum()) / n_fibers
        p_multi = 1.0 - p_unit
        # Uniform expectation: total odd points / (2^(k-1) odd-residue classes * sigma-range)
        sigma_range = int(sigma.max() - sigma.min() + 1)
        uniform_E = len(n) / (M / 2 * sigma_range)
        print(f"{k:>3} {M:>7} {n_fibers:>11,} {med:>7} "
              f"{mn:>8.3f} {mx:>7} {p_unit:>8.4f} {p_multi:>8.4f} "
              f"{uniform_E:>11.3e}")

    print()
    # ============================================================
    # Tier 1b: max-fiber scaling vs k
    # ============================================================
    print("Tier 1b: how fast does max-fiber-size shrink with k?")
    print("If F_max ~ N / 2^k * g(sigma), the table above tells us g(sigma)'s")
    print("amplitude. Below: fitted exponent for max-fiber decay.")
    print()
    max_fibers = []
    for k in ks:
        M = 1 << k
        r = (n % M).astype(np.int64)
        key = (sigma << k) | r
        _, counts = np.unique(key, return_counts=True)
        max_fibers.append(int(counts.max()))
    log_max = np.log2(np.array(max_fibers, dtype=np.float64))
    log_k = np.array(ks, dtype=np.float64)
    if len(ks) >= 2:
        slope = np.polyfit(log_k, log_max, 1)[0]
        print(f"  log2(max_fiber) vs k slope: {slope:>+.3f}")
        print(f"  (uniform-grid prediction: -1.0; actual >  -1 means fibers")
        print(f"   shrink slower than 1/2^k; actual < -1 means structural")
        print(f"   thinning beyond uniform splitting.)")

    print()
    print("=" * 78)
    print(f"Tier 2: within-fiber dispersion at k={args.gap_k}")
    print("=" * 78)
    print()

    # ============================================================
    # Tier 2: structure within the top fibers at gap_k
    # ============================================================
    k = args.gap_k
    M = 1 << k
    r = (n % M).astype(np.int64)
    m = (n - r) // M  # m = integer tail
    key = (sigma << k) | r

    # Group by key, get largest fibers
    sort_idx = np.argsort(key)
    key_sorted = key[sort_idx]
    n_sorted = n[sort_idx]
    m_sorted = m[sort_idx]
    r_sorted = r[sort_idx]
    sigma_sorted = sigma[sort_idx]

    # Run-length encode
    unique_keys, run_starts, run_counts = np.unique(
        key_sorted, return_index=True, return_counts=True)
    # Top fibers
    top_idx = np.argsort(-run_counts)[:args.top_fibers]
    print(f"  top {args.top_fibers} fibers at k={k} (mod {M}):")
    print()
    print(f"  {'rank':>4} {'sigma':>5} {'r':>5} {'F':>6} "
          f"{'m_min':>10} {'m_max':>10} {'gap_mean':>11} {'gap_std':>11} "
          f"{'gap_GCD':>9} {'lag1_acf':>10}")

    pooled_gaps = []
    for rank, ti in enumerate(top_idx):
        start = run_starts[ti]
        cnt = run_counts[ti]
        sl = slice(start, start + cnt)
        ms_in_fiber = np.sort(m_sorted[sl])
        sig_val = sigma_sorted[start]
        r_val = r_sorted[start]
        gaps = np.diff(ms_in_fiber)
        if len(gaps) == 0:
            continue
        gap_mean = float(gaps.mean())
        gap_std = float(gaps.std())
        gap_gcd = int(np.gcd.reduce(gaps.astype(np.int64))) if len(gaps) > 0 else 0
        # Lag-1 autocorrelation of sorted m: detects arithmetic-progression-ness.
        # For a perfect AP, gaps are constant, so lag-1 acf of gaps is undefined
        # (zero variance). We use the gap-of-gaps test: consecutive gap diffs.
        if len(gaps) >= 2:
            gap_diffs = np.diff(gaps.astype(np.float64))
            if gaps.std() > 0:
                lag1 = float(np.corrcoef(gaps[:-1], gaps[1:])[0, 1])
            else:
                lag1 = float("nan")
        else:
            lag1 = float("nan")

        pooled_gaps.append(gaps)
        print(f"  {rank:>4} {sig_val:>5} {r_val:>5} {cnt:>6} "
              f"{int(ms_in_fiber.min()):>10,} {int(ms_in_fiber.max()):>10,} "
              f"{gap_mean:>11.2f} {gap_std:>11.2f} {gap_gcd:>9} "
              f"{lag1:>10.4f}")

    # Pool all gaps from top fibers
    if pooled_gaps:
        pooled = np.concatenate(pooled_gaps)
        print()
        print(f"  Pooled within-fiber gaps (top-{args.top_fibers} fibers, "
              f"n={len(pooled):,}):")
        print(f"    mean        = {pooled.mean():.2f}")
        print(f"    std         = {pooled.std():.2f}")
        print(f"    median      = {np.median(pooled):.0f}")
        print(f"    GCD         = {int(np.gcd.reduce(pooled.astype(np.int64)))}")
        # Histogram in log bins
        # If gaps look uniform on [1, N/2^k], expected mean = N / (2*2^k) ≈ N/(2^(k+1))
        N_eff = len(n) * 2  # we filtered to odd, so original space is 2x
        uniform_gap_mean = N_eff / (M * 2)  # rough; gaps within a fiber span a slice of 2^k
        print(f"    benchmark   = ~N/2^(k+1) ~ {uniform_gap_mean:.0f} "
              f"(if m's were uniform on [0, N/2^k))")

        # Periodicity test: how often does the same gap value appear?
        gap_unique, gap_freq = np.unique(pooled, return_counts=True)
        top_gaps = gap_unique[np.argsort(-gap_freq)[:5]]
        top_freqs = gap_freq[np.argsort(-gap_freq)[:5]]
        print(f"    top-5 most common gap values:")
        for gv, gf in zip(top_gaps, top_freqs):
            print(f"      gap = {int(gv):>10} occurs {int(gf):>6} times "
                  f"({100*gf/len(pooled):.2f}%)")

    print()
    print("=" * 78)
    print("Verdict heuristic")
    print("=" * 78)
    print()
    print("  Tier 1: if max(F) at k=14 is small (<10) for a 5M-row dataset, the")
    print("  problem is bounded by exhaustion at modest k. If max(F) stays large,")
    print("  the fiber structure has the cardinality to support hardness.")
    print()
    print("  Tier 2: if pooled-gap GCD > 1 OR top-5 gap values dominate the")
    print("  histogram, the m-values within fibers have arithmetic structure")
    print("  and lattice reduction is a plausible attack. Hard fibers should")
    print("  show GCD = 1 and a roughly continuous gap distribution.")


if __name__ == "__main__":
    main()
