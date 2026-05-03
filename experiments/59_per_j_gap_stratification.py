"""
σ-quartile × j-class stratification of the −0.26 per-level a★ gap.

Heuristic Δα per a★-level = K_h · log(3) + 1 = 12.46.
Empirical (k-invariant): ≈ 12.20.
Gap: −0.26 per j.

Method:
  For k ∈ {6, 8, 10, 12, 14}, at parquet N=2^27:
    For each odd residue r mod 2^k, compute deterministic_prefix(r, M=2^k)
      → (prefix_steps, a_final). j_r = log_3(a_final) (the a★-level).
      α_pred(r) = prefix_steps + K_h · log(a_final / M)
    Stratify orbits by σ_resid quartile (q1-q4) globally.
    For each (r, q) pair:
      Compute α_emp(r, q) = OLS intercept of σ vs log(n) for orbits in class.
    Compute gap(r, q) = α_emp(r, q) − α_pred(r).

  Aggregate:
    gap_by_quartile(q) = mean gap(r, q) over r at each k → does it concentrate?
    gap_by_j(j) = mean gap(r, q) over (r with j_r = j) → does it concentrate at j?

Verdict criteria (from brief):
  (a) Gap localizes to upper σ-quartiles (q3, q4) → same as constants 3, 4
  (b) Gap localizes to lower σ-quartiles (q1) → same as E[v]_q125 mechanism
  (c) Gap concentrates at specific j ∈ {2, 4, 5} → absorbing-chain machinery
  (d) Gap uniform across both → decoupled mechanism
"""
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

import io
sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)

K_H = 3.0 / (np.log(4.0) - np.log(3.0))
LOG_3 = np.log(3.0)
DELTA_ALPHA_HEURISTIC = K_H * LOG_3 + 1.0  # = 12.46
print(f"# Heuristic per-a★-level Δα = K_h · log(3) + 1 = {DELTA_ALPHA_HEURISTIC:.4f}", flush=True)
print(f"# K_h = {K_H:.4f}", flush=True)


def deterministic_prefix(r, M, max_steps=400):
    """Walk (a, c) from (M, r) until a is odd. Return (steps, a_final)."""
    a, c = M, r
    steps = 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3 * c + 1
        steps += 1
    return steps, a


def trim_mean_quantile(arr, q_low=0.005, q_high=0.995):
    """Trimmed mean (1% trim each tail, default)."""
    if len(arr) < 10:
        return float(arr.mean()) if len(arr) > 0 else float('nan')
    lo = np.quantile(arr, q_low)
    hi = np.quantile(arr, q_high)
    keep = (arr >= lo) & (arr <= hi)
    if keep.sum() == 0:
        return float(arr.mean())
    return float(arr[keep].mean())


def compute_per_class_alpha(n_arr, log_n, sigma, class_idx, K_total, mask_q):
    """For each residue class kk and quartile mask, compute OLS intercept of σ vs log_n.
    Returns alpha_actual[K_total], se_alpha[K_total], n_per[K_total]."""
    alpha = np.full(K_total, np.nan)
    se = np.full(K_total, np.nan)
    npc = np.zeros(K_total, dtype=np.int64)

    # Apply mask
    if mask_q is not None:
        keep = mask_q
        ci = class_idx[keep]; ln = log_n[keep]; sg = sigma[keep]
    else:
        ci = class_idx; ln = log_n; sg = sigma

    if len(ci) == 0: return alpha, se, npc

    # Sort by class
    idx = np.argsort(ci, kind='stable')
    ci_s = ci[idx]; ln_s = ln[idx]; sg_s = sg[idx]
    # Find class boundaries
    boundaries = np.concatenate([[0], np.where(np.diff(ci_s) != 0)[0] + 1, [len(ci_s)]])

    for bi in range(len(boundaries) - 1):
        st, en = boundaries[bi], boundaries[bi+1]
        kk = int(ci_s[st])
        cnt = en - st
        if cnt < 30:  # reduced threshold for per-quartile sub-classes
            continue
        x = ln_s[st:en]; y = sg_s[st:en]
        var_x = x.var()
        if var_x == 0: continue
        b, a = np.polyfit(x, y, 1)
        r = y - (a + b * x)
        alpha[kk] = a
        se[kk] = r.std() * np.sqrt(1.0 / cnt + x.mean()**2 / (cnt * var_x))
        npc[kk] = cnt

    return alpha, se, npc


def analyze_k(k, df, log_n_global, sigma_resid_q25_q50_q75):
    """For modulus 2^k, compute per-(class, σ-quartile) gap."""
    M = 1 << k
    K_total = M // 2

    n = df["n"].to_numpy().astype(np.int64)
    log_n = log_n_global
    sigma = df["sigma"].to_numpy().astype(np.float64)

    res = (n % M).astype(np.int32)
    class_idx = ((res - 1) // 2).astype(np.int32)

    # σ_resid quartiles (computed globally for consistency)
    q25, q50, q75 = sigma_resid_q25_q50_q75
    sigma_resid = sigma - K_H * log_n  # detrend by K_h · log_n
    quartile = np.zeros(len(sigma), dtype=np.int8)  # 0 = q1, 1 = q2, 2 = q3, 3 = q4
    quartile[(sigma_resid > q25) & (sigma_resid <= q50)] = 1
    quartile[(sigma_resid > q50) & (sigma_resid <= q75)] = 2
    quartile[sigma_resid > q75] = 3

    # Compute prefix prediction per class
    print(f"  k={k}, M={M}, K_total={K_total}: computing prefixes...", flush=True)
    t0 = time.perf_counter()
    prefix_steps_arr = np.zeros(K_total, dtype=np.int64)
    a_final_arr = np.zeros(K_total, dtype=np.int64)
    for kk in range(K_total):
        r = 2 * kk + 1
        ps, af = deterministic_prefix(r, M)
        prefix_steps_arr[kk] = ps
        a_final_arr[kk] = af
    log_a_norm = np.log(a_final_arr / float(M))
    alpha_pred = prefix_steps_arr + K_H * log_a_norm
    print(f"    prefix done in {time.perf_counter()-t0:.1f}s", flush=True)

    # j-index per class: j_r = log_3(a_final) — the a★-level
    j_per_class = np.round(np.log(a_final_arr.astype(np.float64)) / LOG_3).astype(np.int32)

    # Per-quartile α_emp
    print(f"  computing α per (class, σ-quartile)...", flush=True)
    t0 = time.perf_counter()
    results_per_q = {}
    for q in range(4):
        mask_q = quartile == q
        alpha_emp, se_emp, npc_emp = compute_per_class_alpha(
            n, log_n, sigma, class_idx, K_total, mask_q)
        results_per_q[q] = (alpha_emp, se_emp, npc_emp)
        print(f"    q{q}: {(npc_emp > 0).sum()} populated classes, "
              f"min/median/max n = {npc_emp[npc_emp>0].min()} / "
              f"{int(np.median(npc_emp[npc_emp>0]))} / {npc_emp[npc_emp>0].max()}", flush=True)
    print(f"    α-fit done in {time.perf_counter()-t0:.1f}s", flush=True)

    # Also unstratified (all orbits)
    alpha_emp_all, se_emp_all, npc_all = compute_per_class_alpha(
        n, log_n, sigma, class_idx, K_total, None)

    # Compute gap = α_emp − α_pred per (q, class), aggregate
    summary = {'k': k, 'per_q': {}, 'per_j_per_q': {}}
    for q in range(4):
        ae, se_e, npc = results_per_q[q]
        gap = ae - alpha_pred
        valid = ~np.isnan(gap)
        summary['per_q'][q] = {
            'mean_gap': trim_mean_quantile(gap[valid]) if valid.sum() > 0 else float('nan'),
            'median_gap': float(np.median(gap[valid])) if valid.sum() > 0 else float('nan'),
            'sd_gap': float(np.std(gap[valid])) if valid.sum() > 0 else float('nan'),
            'n_classes': int(valid.sum()),
        }
        # Also per-j aggregation
        per_j = {}
        for j in np.unique(j_per_class):
            j_mask = (j_per_class == j) & valid
            if j_mask.sum() > 0:
                per_j[int(j)] = {
                    'mean_gap': trim_mean_quantile(gap[j_mask]),
                    'n_classes': int(j_mask.sum()),
                }
        summary['per_j_per_q'][q] = per_j

    # Unstratified
    gap_all = alpha_emp_all - alpha_pred
    valid_all = ~np.isnan(gap_all)
    summary['all'] = {
        'mean_gap': trim_mean_quantile(gap_all[valid_all]),
        'median_gap': float(np.median(gap_all[valid_all])),
        'n_classes': int(valid_all.sum()),
    }

    # Also expose per-j gap aggregated across all quartiles
    summary['per_j_all'] = {}
    for j in np.unique(j_per_class):
        j_mask = (j_per_class == j) & valid_all
        if j_mask.sum() > 0:
            summary['per_j_all'][int(j)] = {
                'mean_gap': trim_mean_quantile(gap_all[j_mask]),
                'n_classes': int(j_mask.sum()),
            }

    return summary


def main():
    here = Path(__file__).resolve().parent
    data_dir = here.parent / "data"
    out_dir = here.parent / "experiments_output"

    N = 1 << 27
    print(f"\n# Loading parquet at N=2^27 ({N:,})", flush=True)
    t0 = time.perf_counter()
    df = pl.read_parquet(data_dir / f"main_N{N}.parquet").filter(
        (pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    n_arr = df["n"].to_numpy().astype(np.int64)
    log_n_global = np.log(n_arr.astype(np.float64))
    sigma = df["sigma"].to_numpy().astype(np.float64)
    print(f"  Loaded {len(n_arr):,} odd orbits in {time.perf_counter()-t0:.1f}s", flush=True)

    # Global σ_resid quartiles
    sigma_resid = sigma - K_H * log_n_global
    q25 = float(np.percentile(sigma_resid, 25))
    q50 = float(np.percentile(sigma_resid, 50))
    q75 = float(np.percentile(sigma_resid, 75))
    print(f"  σ_resid quartiles: q25={q25:.2f}, q50={q50:.2f}, q75={q75:.2f}", flush=True)

    all_summaries = []
    for k in [6, 8, 10, 12, 14]:
        print(f"\n========== k = {k} ==========", flush=True)
        s = analyze_k(k, df, log_n_global, (q25, q50, q75))
        all_summaries.append(s)

        print(f"  {'k='+str(k):>4}: gap (mean over classes)", flush=True)
        for q in range(4):
            r = s['per_q'][q]
            print(f"    q{q+1}: mean gap = {r['mean_gap']:+8.4f}  median = {r['median_gap']:+8.4f}  "
                  f"sd = {r['sd_gap']:.4f}  n_classes = {r['n_classes']}", flush=True)
        print(f"    all: mean gap = {s['all']['mean_gap']:+8.4f}  median = {s['all']['median_gap']:+8.4f}  "
              f"n_classes = {s['all']['n_classes']}", flush=True)

        # Per-j summary (top-5 j-values by class count)
        per_j = s['per_j_all']
        if per_j:
            j_sorted = sorted(per_j.keys(), key=lambda j: -per_j[j]['n_classes'])
            print(f"    Top-5 j-values by # classes:", flush=True)
            print(f"      {'j':>4}  {'gap (all)':>10}  {'n_cl':>6}", end='', flush=True)
            for q in range(4):
                print(f"  {'gap q'+str(q+1):>9}", end='', flush=True)
            print('', flush=True)
            for j in j_sorted[:8]:
                d = per_j[j]
                print(f"      {j:>4}  {d['mean_gap']:>+10.4f}  {d['n_classes']:>6}", end='', flush=True)
                for q in range(4):
                    pq = s['per_j_per_q'][q].get(j, None)
                    if pq:
                        print(f"  {pq['mean_gap']:>+9.4f}", end='', flush=True)
                    else:
                        print(f"  {'—':>9}", end='', flush=True)
                print('', flush=True)

    # Cross-k summary
    print(f"\n\n========== Cross-k summary ==========", flush=True)
    print(f"\n# Mean gap by σ-quartile (across k):", flush=True)
    print(f"  {'k':>3}", end='', flush=True)
    for q in range(4):
        print(f"  {'q'+str(q+1):>10}", end='', flush=True)
    print(f"  {'all':>10}", flush=True)
    for s in all_summaries:
        print(f"  {s['k']:>3}", end='', flush=True)
        for q in range(4):
            print(f"  {s['per_q'][q]['mean_gap']:>+10.4f}", end='', flush=True)
        print(f"  {s['all']['mean_gap']:>+10.4f}", flush=True)

    # σ-quartile localization test
    print(f"\n# σ-quartile localization test (mean gap variability across q at each k):", flush=True)
    for s in all_summaries:
        gaps_q = np.array([s['per_q'][q]['mean_gap'] for q in range(4)])
        spread = float(gaps_q.max() - gaps_q.min())
        print(f"  k={s['k']}: q1={gaps_q[0]:+.3f}, q2={gaps_q[1]:+.3f}, "
              f"q3={gaps_q[2]:+.3f}, q4={gaps_q[3]:+.3f}; spread={spread:.3f}", flush=True)

    # Verdict
    print(f"\n# VERDICT (per brief):", flush=True)
    # Compute average q-pattern across k
    avg_per_q = np.zeros(4)
    for s in all_summaries:
        for q in range(4):
            avg_per_q[q] += s['per_q'][q]['mean_gap']
    avg_per_q /= len(all_summaries)
    print(f"  Avg gap by quartile across k: q1={avg_per_q[0]:+.3f}, q2={avg_per_q[1]:+.3f}, "
          f"q3={avg_per_q[2]:+.3f}, q4={avg_per_q[3]:+.3f}", flush=True)
    spread = avg_per_q.max() - avg_per_q.min()
    print(f"  Spread q-min to q-max: {spread:.3f}", flush=True)
    if spread < 0.05:
        print(f"  → (d) Gap is uniform across σ-quartiles (spread < 0.05)", flush=True)
    elif avg_per_q[2:4].mean() - avg_per_q[0:2].mean() < -0.1:
        print(f"  → (a) Gap localizes to LOWER σ-quartiles", flush=True)
    elif avg_per_q[2:4].mean() - avg_per_q[0:2].mean() > 0.1:
        print(f"  → Gap localizes to UPPER σ-quartiles (constants-3,4 phenomenon)", flush=True)
    else:
        print(f"  → Spread present but not strongly directional", flush=True)

    # Save CSVs
    rows = []
    for s in all_summaries:
        for q in range(4):
            r = s['per_q'][q]
            rows.append({
                'k': s['k'], 'sigma_quartile': q + 1,
                'mean_gap': r['mean_gap'], 'median_gap': r['median_gap'],
                'sd_gap': r['sd_gap'], 'n_classes': r['n_classes'],
            })
        rows.append({
            'k': s['k'], 'sigma_quartile': 0,
            'mean_gap': s['all']['mean_gap'],
            'median_gap': s['all']['median_gap'], 'sd_gap': float('nan'),
            'n_classes': s['all']['n_classes'],
        })
    pl.DataFrame(rows).write_csv(out_dir / "59_per_j_gap_stratification_quartile.csv")

    rows_j = []
    for s in all_summaries:
        for j, d in s['per_j_all'].items():
            row = {'k': s['k'], 'j': j, 'gap_all': d['mean_gap'], 'n_classes': d['n_classes']}
            for q in range(4):
                pq = s['per_j_per_q'][q].get(j, None)
                row[f'gap_q{q+1}'] = pq['mean_gap'] if pq else float('nan')
                row[f'n_q{q+1}'] = pq['n_classes'] if pq else 0
            rows_j.append(row)
    pl.DataFrame(rows_j).write_csv(out_dir / "59_per_j_gap_stratification_j.csv")
    print(f"\n[save] CSVs", flush=True)


if __name__ == "__main__":
    main()
