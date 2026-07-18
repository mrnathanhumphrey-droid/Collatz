"""
test_tail_shape.py - Test 3: stopping-time tail shape vs K_eff deficit.

Hypothesis: K_h - K_eff = 0.486 reflects a tail-weight effect — the right tail
of σ at fixed log(n) shifts with log(n) faster than the median, so quantile
slopes differ from K_h.

Method:
  1. Walk 1M orbits at each of N=2^25 and N=2^27.
  2. Bin by octave of log_2(n).
  3. Per octave: compute σ quantiles {0.5, 0.75, 0.9, 0.95, 0.99}, fit GPD tail.
  4. Regress each quantile on log_2(n) octave midpoint -> per-quantile slope.
  5. Compare to K_h = 10.43 (in nat units; in log_2 units = 10.43 * ln(2) = 7.23).

Decisive:
  - If high-quantile slopes match K_h, median below: tail NOT mechanism.
  - If high-quantile slopes differ from K_h, integrate to K_eff: tail IS mechanism.
"""
import math
import time
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from scipy import stats

K_H_NAT = 3.0 / (math.log(4.0) - math.log(3.0))  # 10.43 per ln(n)
K_H_LOG2 = K_H_NAT * math.log(2.0)               # = K_h per log_2(n) ≈ 7.23

OUT_DIR = Path(__file__).parent / "viz_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)


@njit(cache=True, parallel=True)
def walk_sigma_batch(starts, max_steps):
    n_arr = len(starts)
    sigmas = np.zeros(n_arr, dtype=np.int64)
    for i in prange(n_arr):
        n = starts[i]
        s = -1
        for step in range(max_steps):
            if n == 1:
                s = step
                break
            if n % 2 == 0:
                n = n >> 1
            else:
                n = 3 * n + 1
        sigmas[i] = s
    return sigmas


def generate_sample(N_log2, n_samples, seed=42):
    """Sample n_samples odd values from [3, 2^N_log2], walk to get sigma."""
    out_path = OUT_DIR / f"tail_sample_2pow{N_log2}_n{n_samples}.csv"
    if out_path.exists():
        print(f"[gen] using cached {out_path.name}")
        return pl.read_csv(out_path)

    print(f"[gen] sampling {n_samples:,} odd orbits from [3, 2^{N_log2}] ...")
    rng = np.random.default_rng(seed)
    N = 1 << N_log2
    # Sample odd n uniformly from [3, N-1]
    odd_pool_size = (N - 2) // 2  # count of odd values in [3, N-1]
    indices = rng.choice(odd_pool_size, size=n_samples, replace=False)
    starts = (3 + 2 * indices).astype(np.int64)

    # JIT warmup
    _ = walk_sigma_batch(np.array([7], dtype=np.int64), 5000)

    t0 = time.perf_counter()
    sigmas = walk_sigma_batch(starts, max_steps=10_000)
    elapsed = time.perf_counter() - t0
    print(f"[gen] walk done in {elapsed:.1f}s")

    df = pl.DataFrame({
        "n": starts,
        "log2_n": np.log2(starts.astype(np.float64)),
        "sigma": sigmas,
    })
    df.write_csv(out_path)
    print(f"[gen] wrote {out_path}")
    return df


def analyze_octaves(df, N_log2, label):
    print(f"\n{'='*100}")
    print(f"  {label}  (N=2^{N_log2}, n_orbits={len(df):,})")
    print(f"{'='*100}")
    log2_n = df["log2_n"].to_numpy()
    sigma = df["sigma"].to_numpy().astype(np.float64)

    # Octave bins from 20 to N_log2 (lower octaves have few samples)
    octave_edges = np.arange(20, N_log2 + 1)
    octave_mid = (octave_edges[:-1] + octave_edges[1:]) / 2.0

    quantile_levels = [0.5, 0.75, 0.9, 0.95, 0.99]
    octave_quantiles = {q: [] for q in quantile_levels}
    octave_means = []
    octave_n = []
    octave_xi_gpd = []

    print(f"  {'octave':>8}  {'n_orb':>8}  {'mean_sigma':>10}  {'q50':>7}  {'q75':>7}  "
          f"{'q90':>7}  {'q95':>7}  {'q99':>7}  {'GPD ξ':>9}")
    for i in range(len(octave_edges) - 1):
        lo, hi = octave_edges[i], octave_edges[i+1]
        mask = (log2_n >= lo) & (log2_n < hi)
        n_oct = int(mask.sum())
        if n_oct < 200:
            print(f"  [{lo:>2d},{hi:>2d}]  {n_oct:>8d}  (insufficient, skip)")
            for q in quantile_levels:
                octave_quantiles[q].append(np.nan)
            octave_means.append(np.nan)
            octave_n.append(n_oct)
            octave_xi_gpd.append(np.nan)
            continue
        sub = sigma[mask]
        means = sub.mean()
        qs = {q: np.quantile(sub, q) for q in quantile_levels}
        # Fit GPD to top 10% tail
        threshold = qs[0.9]
        tail = sub[sub > threshold] - threshold
        if len(tail) > 50:
            try:
                xi, loc, scale = stats.genpareto.fit(tail, floc=0)
            except Exception:
                xi = np.nan
        else:
            xi = np.nan
        for q in quantile_levels:
            octave_quantiles[q].append(qs[q])
        octave_means.append(means)
        octave_n.append(n_oct)
        octave_xi_gpd.append(xi)

        print(f"  [{lo:>2d},{hi:>2d}]  {n_oct:>8d}  {means:>10.2f}  "
              f"{qs[0.5]:>7.1f}  {qs[0.75]:>7.1f}  {qs[0.9]:>7.1f}  "
              f"{qs[0.95]:>7.1f}  {qs[0.99]:>7.1f}  {xi:>+9.4f}")

    # Per-quantile slopes (regress quantile on octave midpoint, in log_2 units)
    print(f"\n  Per-quantile slopes (σ-quantile vs log_2(n)):")
    print(f"  {'quantile':>10}  {'slope (log_2)':>14}  {'slope (ln)':>11}  {'gap from K_h':>13}")
    print(f"  {'mean':>10}  ", end="")
    valid = ~np.isnan(octave_means)
    if valid.sum() >= 3:
        sl, _ = np.polyfit(octave_mid[valid], np.array(octave_means)[valid], 1)
        sl_nat = sl / math.log(2.0)
        print(f"{sl:>14.4f}  {sl_nat:>11.4f}  {sl_nat - K_H_NAT:>+13.4f}")
    else:
        print("(insufficient octaves)")

    quantile_slopes = {}
    for q in quantile_levels:
        vals = np.array(octave_quantiles[q])
        valid = ~np.isnan(vals)
        if valid.sum() >= 3:
            sl, _ = np.polyfit(octave_mid[valid], vals[valid], 1)
            sl_nat = sl / math.log(2.0)
            quantile_slopes[q] = (sl, sl_nat)
            print(f"  {q:>10}  {sl:>14.4f}  {sl_nat:>11.4f}  {sl_nat - K_H_NAT:>+13.4f}")

    # GPD ξ trend
    print(f"\n  GPD tail-shape ξ vs octave:")
    xi_arr = np.array(octave_xi_gpd)
    valid = ~np.isnan(xi_arr)
    if valid.sum() >= 3:
        xi_slope, _ = np.polyfit(octave_mid[valid], xi_arr[valid], 1)
        print(f"  ξ values: {xi_arr[valid]}")
        print(f"  ξ slope on log_2(n): {xi_slope:+.4f}")
    return {
        "N_log2": N_log2,
        "octave_mid": octave_mid,
        "octave_means": np.array(octave_means),
        "octave_quantiles": {q: np.array(v) for q, v in octave_quantiles.items()},
        "octave_xi": np.array(octave_xi_gpd),
        "quantile_slopes": quantile_slopes,
    }


def main():
    print(f"K_h = 3/log(4/3) = {K_H_NAT:.5f} per ln(n) = {K_H_LOG2:.5f} per log_2(n)")
    print(f"K_eff (from closed_form_findings.md) ≈ 9.94 per ln(n) = {9.94 * math.log(2):.5f} per log_2(n)")
    print(f"Slope K_h - K_eff = 0.486 per ln(n) = {0.486 * math.log(2):.5f} per log_2(n)")
    print()

    N_LIST = [25, 27]
    N_SAMPLES = 1_000_000

    results = []
    for N_log2 in N_LIST:
        df = generate_sample(N_log2, N_SAMPLES)
        result = analyze_octaves(df, N_log2, f"N=2^{N_log2}")
        results.append(result)

    # Cross-N comparison
    print(f"\n{'='*100}")
    print(f"  CROSS-N COMPARISON: per-quantile slopes (per ln(n))")
    print(f"{'='*100}")
    print(f"  {'quantile':>10}", end="")
    for r in results:
        print(f"  {'N=2^'+str(r['N_log2'])+' slope':>17}", end="")
    print(f"  {'K_h - q_slope':>15}")

    for q in [0.5, 0.75, 0.9, 0.95, 0.99]:
        print(f"  {q:>10}  ", end="")
        slopes_nat = []
        for r in results:
            if q in r["quantile_slopes"]:
                _, sl_nat = r["quantile_slopes"][q]
                slopes_nat.append(sl_nat)
                print(f"{sl_nat:>17.4f}  ", end="")
            else:
                slopes_nat.append(np.nan)
                print(f"{'--':>17}  ", end="")
        if slopes_nat:
            avg = np.nanmean(slopes_nat)
            print(f"{K_H_NAT - avg:>+15.4f}")


if __name__ == "__main__":
    main()
