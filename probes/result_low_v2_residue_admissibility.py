"""
result_low_v2_residue_admissibility.py — empirical probe of whether low-v_2
trajectories have shorter / longer 2-adic-bit-exhaustion horizons than
typical trajectories.

Definition (canonical 2-adic admissibility horizon):
  Starting from n mod 2^m known, after k Syracuse steps T^k(n) mod 2^(m - Σ v_i)
  is determined; once Σ v_i ≥ m, the residue trace at level 2^m is no longer
  determined by the initial residue alone — bits beyond 2^m are required.
  k_max(m) := first k with cumulative Σ_{i<=k} v_2(3n_i+1) ≥ m.

  TYPICAL trajectory: E[v_2] ≈ 2, so k_max(m) ≈ m/2.
  LOW-v_2 trajectory: E[v_2] < 2, so k_max(m) > m/2 (LARGER, not smaller).

The brief's hypothesis is that low-v_2 trajectories "run out of admissible
residues earlier" (smaller k_max). The bit-budget calculation predicts the
opposite. This probe computes both and reports.

Procedure:
  1. Iterate all odd n coprime to 3 in [1, 10^7] under Syracuse for up to
     1000 steps (most collapse far earlier).
  2. Per trajectory: track Σ v_2, count of steps with v_2 ≥ t (t=2..6),
     and on-the-fly compute k_max(m) for m ∈ {6,8,10,12,16,20}.
  3. Stratify: bottom q ∈ {0.001, 0.01, 0.05} by count of v_2 ≥ 2.
  4. Compare horizon distributions.
"""
from __future__ import annotations
import csv
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = Path(r"C:\Collatz")
OUT_MD = OUTDIR / "result_low_v2_residue_admissibility.md"
OUT_CSV = OUTDIR / "result_low_v2_residue_admissibility.csv"

N_MAX = 10_000_000
MAX_STEPS = 1000
THRESHOLDS = [2, 3, 4, 5, 6]
M_VALUES = [6, 8, 10, 12, 16, 20]
QUANTILES = [0.001, 0.01, 0.05]
INT64_GUARD = 1 << 60


def main():
    t0 = time.time()
    print(f"[1/4] generating starts; iterating up to {MAX_STEPS} Syracuse "
          f"steps", flush=True)
    arr = np.arange(1, N_MAX + 1, 2, dtype=np.int64)
    arr = arr[arr % 3 != 0]
    arr = arr[arr != 1]
    starts = arr
    N = len(starts)
    print(f"  {N:,} starting integers", flush=True)

    n_arr = starts.copy()
    alive = np.ones(N, dtype=bool)
    step = np.zeros(N, dtype=np.int32)
    length = np.zeros(N, dtype=np.int32)
    cum_v2 = np.zeros(N, dtype=np.int32)
    counts_thresh = np.zeros((N, len(THRESHOLDS)), dtype=np.int32)
    horizon = np.full((N, len(M_VALUES)), -1, dtype=np.int32)
    overflow_seen = False

    for s in range(MAX_STEPS):
        sub = np.where(alive)[0]
        if sub.size == 0:
            break
        ns = n_arr[sub]
        m = 3 * ns + 1
        if not overflow_seen and (m > INT64_GUARD).any():
            print(f"  WARN: int64 overflow risk at step {s}, "
                  f"max m = {m.max():,}", flush=True)
            overflow_seen = True
        v = np.zeros(len(sub), dtype=np.int32)
        m_work = m.copy()
        while True:
            even = (m_work & 1) == 0
            if not even.any():
                break
            v[even] += 1
            m_work[even] >>= 1
        m_final = m_work

        # Update cumulative v_2 and threshold counts
        cum_v2[sub] += v
        for ti, t in enumerate(THRESHOLDS):
            counts_thresh[sub, ti] += (v >= t).astype(np.int32)

        # On-the-fly horizon detection
        for mi, m_val in enumerate(M_VALUES):
            # for indices in sub: horizon not set AND cum_v2 just reached m_val
            sub_horizon_unset = horizon[sub, mi] == -1
            sub_reached = cum_v2[sub] >= m_val
            sub_newly = sub_horizon_unset & sub_reached
            if sub_newly.any():
                horizon[sub[sub_newly], mi] = step[sub[sub_newly]] + 1

        n_arr[sub] = m_final
        step[sub] += 1
        new_dead = m_final == 1
        if new_dead.any():
            dead_idx = sub[new_dead]
            length[dead_idx] = step[dead_idx]
            alive[dead_idx] = False

        if (s + 1) % 100 == 0:
            n_alive = int(alive.sum())
            print(f"  step {s+1}: alive={n_alive:,}", flush=True)
            if n_alive == 0:
                break

    still_alive = np.where(alive)[0]
    length[still_alive] = MAX_STEPS
    print(f"  iteration done in {time.time()-t0:.1f}s; "
          f"{still_alive.size:,} did not collapse within {MAX_STEPS} steps",
          flush=True)

    # Diagnostics on cum_v2 and length
    print(f"\n[2/4] descriptive stats on full ensemble", flush=True)
    print(f"  total v_2 (all trajectories): "
          f"mean={cum_v2.mean():.2f}, median={np.median(cum_v2):.0f}, "
          f"min={cum_v2.min()}, max={cum_v2.max()}", flush=True)
    print(f"  trajectory length: "
          f"mean={length.mean():.2f}, median={np.median(length):.0f}, "
          f"max={length.max()}", flush=True)
    avg_v2_per_step = cum_v2.astype(np.float64) / np.maximum(length, 1)
    print(f"  avg v_2 per step: "
          f"mean={avg_v2_per_step.mean():.4f}, "
          f"median={np.median(avg_v2_per_step):.4f}", flush=True)

    print(f"\n  fraction with horizon reached:", flush=True)
    for mi, m_val in enumerate(M_VALUES):
        reached = (horizon[:, mi] >= 0).sum()
        not_reached = N - reached
        if reached > 0:
            mean_h = horizon[:, mi][horizon[:, mi] >= 0].mean()
            median_h = np.median(horizon[:, mi][horizon[:, mi] >= 0])
        else:
            mean_h = float("nan")
            median_h = float("nan")
        print(f"    m={m_val:>2}: reached {100*reached/N:.2f}%  "
              f"(mean horizon {mean_h:.2f}, median {median_h:.0f}; "
              f"{not_reached:,} did NOT reach m bits)", flush=True)

    # Stratify by bottom-q on counts_thresh[:,0] (v_2 ≥ 2 count)
    # That's the most natural "low-valuation" criterion.
    print(f"\n[3/4] stratifying by bottom-q quantile on count(v_2 ≥ 2)",
          flush=True)
    primary_count = counts_thresh[:, 0]  # v_2 ≥ 2 count per trajectory
    print(f"  count(v_2 ≥ 2) summary: mean={primary_count.mean():.2f}, "
          f"median={np.median(primary_count):.0f}, "
          f"min={primary_count.min()}, max={primary_count.max()}",
          flush=True)

    # NOTE: count of "v_2 ≥ 2" steps is mostly correlated with trajectory
    # length. Bottom q ≈ "shortest trajectories with few high-v_2 steps".
    # That's actually starts that collapsed quickly. To get "long
    # trajectories with low-v_2", use AVG v_2 per step.
    primary_metric = avg_v2_per_step
    print(f"  using PRIMARY metric: avg v_2 per step (length-normalized)",
          flush=True)
    print(f"  avg_v2 summary: mean={primary_metric.mean():.4f}, "
          f"min={primary_metric.min():.4f}, "
          f"max={primary_metric.max():.4f}",
          flush=True)

    rows_out = []  # for CSV
    print(f"\n  bottom-q quantile cuts:", flush=True)
    for q in QUANTILES:
        cut = np.quantile(primary_metric, q)
        low_mask = primary_metric <= cut
        n_low = low_mask.sum()
        print(f"\n  q={q:.3f}: cut={cut:.4f}, {n_low:,} trajectories at "
              f"or below", flush=True)
        # Compare horizon distributions
        print(f"    {'m':>3} {'mean h (low)':>14} {'mean h (rest)':>14} "
              f"{'median h (low)':>16} {'median h (rest)':>16} "
              f"{'frac reached (low)':>20} {'frac reached (rest)':>20}",
              flush=True)
        for mi, m_val in enumerate(M_VALUES):
            low_h = horizon[low_mask, mi]
            rest_h = horizon[~low_mask, mi]
            low_reached = low_h[low_h >= 0]
            rest_reached = rest_h[rest_h >= 0]
            mean_low = low_reached.mean() if len(low_reached) > 0 else float("nan")
            mean_rest = rest_reached.mean() if len(rest_reached) > 0 else float("nan")
            med_low = np.median(low_reached) if len(low_reached) > 0 else float("nan")
            med_rest = np.median(rest_reached) if len(rest_reached) > 0 else float("nan")
            frac_low = len(low_reached) / max(n_low, 1)
            frac_rest = len(rest_reached) / max(N - n_low, 1)
            print(f"    {m_val:>3} {mean_low:>14.4f} {mean_rest:>14.4f} "
                  f"{med_low:>16.1f} {med_rest:>16.1f} "
                  f"{frac_low:>20.4f} {frac_rest:>20.4f}", flush=True)
            rows_out.append({
                "q": q, "m": m_val,
                "n_low": int(n_low),
                "mean_horizon_low": float(mean_low),
                "mean_horizon_rest": float(mean_rest),
                "median_horizon_low": float(med_low),
                "median_horizon_rest": float(med_rest),
                "frac_reached_low": float(frac_low),
                "frac_reached_rest": float(frac_rest),
                "ratio_means_low_over_rest": (
                    float(mean_low / mean_rest)
                    if mean_rest > 0 else float("nan")
                ),
            })

    # Verdict
    # Are mean horizons systematically larger or smaller in low-v group?
    # If low-v has CONSISTENTLY larger horizons across (q, m): bit-budget
    # calculation confirmed.
    # If low-v has CONSISTENTLY smaller horizons: brief's hypothesis confirmed.
    # If mixed or null: report.
    larger_count = sum(
        1 for r in rows_out
        if (not np.isnan(r["ratio_means_low_over_rest"])
            and r["ratio_means_low_over_rest"] > 1.0)
    )
    smaller_count = sum(
        1 for r in rows_out
        if (not np.isnan(r["ratio_means_low_over_rest"])
            and r["ratio_means_low_over_rest"] < 1.0)
    )
    total = len(rows_out)
    print(f"\n[4/4] verdict on horizon direction:", flush=True)
    print(f"  cells where low-v has LARGER mean horizon than rest: "
          f"{larger_count}/{total}", flush=True)
    print(f"  cells where low-v has SMALLER mean horizon than rest: "
          f"{smaller_count}/{total}", flush=True)

    if larger_count > smaller_count * 3:
        verdict = ("BIT-BUDGET CONFIRMED. Low-v_2 trajectories have "
                   "SYSTEMATICALLY LARGER admissibility horizons than "
                   "typical trajectories. The brief's hypothesis (smaller "
                   "k_max for low-v) is REFUTED — the math predicts and "
                   "the data confirms that lower per-step v_2 means more "
                   "steps to exhaust m bits.")
    elif smaller_count > larger_count * 3:
        verdict = ("BRIEF'S HYPOTHESIS SUPPORTED. Low-v_2 trajectories "
                   "have SYSTEMATICALLY SMALLER admissibility horizons. "
                   "Some non-trivial mechanism (beyond bit accumulation) "
                   "must explain this. Worth characterizing.")
    else:
        verdict = ("MIXED / NULL. No consistent direction across (q, m) "
                   "cells. The 2-adic admissibility horizon doesn't "
                   "stratify cleanly by v_2 deficit at this scale.")
    print(f"\n*** VERDICT: {verdict}", flush=True)

    # Write CSV
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["q", "m", "n_low", "mean_horizon_low",
                    "mean_horizon_rest", "median_horizon_low",
                    "median_horizon_rest", "frac_reached_low",
                    "frac_reached_rest", "ratio_means_low_over_rest"])
        for r in rows_out:
            w.writerow([r["q"], r["m"], r["n_low"],
                        f"{r['mean_horizon_low']:.4f}",
                        f"{r['mean_horizon_rest']:.4f}",
                        f"{r['median_horizon_low']:.1f}",
                        f"{r['median_horizon_rest']:.1f}",
                        f"{r['frac_reached_low']:.6f}",
                        f"{r['frac_reached_rest']:.6f}",
                        f"{r['ratio_means_low_over_rest']:.6f}"])
    print(f"  saved {OUT_CSV}", flush=True)

    # Markdown
    L = []
    L.append("# Low-v_2 Residue Admissibility — empirical horizon test")
    L.append("")
    L.append(f"**Date:** 2026-05-05.  N = {N:,} odd-coprime-to-3 starts in "
             f"[1, 10^7]; iterated up to {MAX_STEPS} Syracuse steps "
             f"(most collapse to 1 in <100 steps).")
    L.append("")
    L.append("**Admissibility horizon definition.**  k_max(m) = first step "
             "k where cumulative Σ_{i≤k} v_2(3n_i+1) ≥ m. Beyond k_max, the "
             "residue r_k mod 2^m is no longer determined by the initial "
             "residue mod 2^m alone — bits beyond 2^m are required to "
             "continue the trace deterministically.")
    L.append("")
    L.append("Bit-budget arithmetic predicts: low-v_2 trajectories take "
             "MORE steps to consume m bits → LARGER k_max. The brief's "
             "hypothesis is the opposite (smaller k_max for low-v).")
    L.append("")

    L.append("## Verdict")
    L.append("")
    L.append(verdict)
    L.append("")

    L.append("## Ensemble stats")
    L.append("")
    L.append("| metric | value |")
    L.append("|---|---:|")
    L.append(f"| total v_2 mean | {cum_v2.mean():.2f} |")
    L.append(f"| total v_2 median | {np.median(cum_v2):.0f} |")
    L.append(f"| trajectory length mean | {length.mean():.2f} |")
    L.append(f"| trajectory length median | {np.median(length):.0f} |")
    L.append(f"| avg v_2 / step (mean) | {avg_v2_per_step.mean():.4f} |")
    L.append(f"| avg v_2 / step (median) | "
             f"{np.median(avg_v2_per_step):.4f} |")
    L.append("")

    L.append("## Horizon reached fraction by m (full ensemble)")
    L.append("")
    L.append("| m | frac reached | mean horizon | median horizon |")
    L.append("|---:|---:|---:|---:|")
    for mi, m_val in enumerate(M_VALUES):
        h_col = horizon[:, mi]
        reached = h_col[h_col >= 0]
        mh = reached.mean() if len(reached) > 0 else float("nan")
        med_h = np.median(reached) if len(reached) > 0 else float("nan")
        L.append(f"| {m_val} | {len(reached)/N:.4f} | {mh:.2f} | "
                 f"{med_h:.1f} |")
    L.append("")

    L.append("## Horizons stratified by bottom-q on (avg v_2 per step)")
    L.append("")
    L.append("Definition: 'low-v_2' = trajectories with avg v_2 per step "
             "in the bottom-q quantile. Comparing mean k_max(m) for "
             "(low) vs (rest).")
    L.append("")
    L.append("Ratio = mean horizon (low) / mean horizon (rest).  "
             "Ratio > 1 ↔ low-v has LARGER horizons (bit-budget prediction). "
             " Ratio < 1 ↔ low-v has SMALLER horizons (brief's hypothesis).")
    L.append("")
    L.append("| q | m | n_low | mean h (low) | mean h (rest) | "
             "median h (low) | median h (rest) | ratio |")
    L.append("|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows_out:
        L.append(f"| {r['q']:.3f} | {r['m']} | {r['n_low']:,} | "
                 f"{r['mean_horizon_low']:.2f} | "
                 f"{r['mean_horizon_rest']:.2f} | "
                 f"{r['median_horizon_low']:.1f} | "
                 f"{r['median_horizon_rest']:.1f} | "
                 f"{r['ratio_means_low_over_rest']:.4f} |")
    L.append("")
    L.append(f"**Cells with ratio > 1 (low has larger horizon):** "
             f"{larger_count}/{total}")
    L.append(f"**Cells with ratio < 1 (low has smaller horizon):** "
             f"{smaller_count}/{total}")
    L.append("")

    L.append("## Files")
    L.append("")
    L.append("- `result_low_v2_residue_admissibility.py` — script")
    L.append("- `result_low_v2_residue_admissibility.csv` — per-(q, m) "
             "horizon stats")
    L.append("- `result_low_v2_residue_admissibility.md` — this writeup")
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"  saved {OUT_MD}", flush=True)
    print(f"\n[done] runtime = {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
