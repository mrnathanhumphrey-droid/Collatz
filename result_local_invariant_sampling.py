"""
result_local_invariant_sampling.py — exploratory soil sampling of marginal
local invariants on Syracuse trajectories.

Five tests on per-step data from all odd-coprime-to-3 starts in [1, 10^7]:
  T1: v_2(3n+1) distribution vs geometric null.
  T2: v_2 conditional on n mod 27 / mod 81 (residue-conditional 2-adic).
  T3: sibling-prime joint distributions (mod 5, 7, 11) × mod 27 vs CRT.
  T4: return-time τ(n) conditional on starting residue mod 27.
  T5: trajectory residue visit frequencies vs framework's π_3 stationary.

Histograms accumulated incrementally; per-step data is NOT stored on disk.
Bonferroni-corrected significance threshold p < 0.05/5 = 0.01 across tests.
"""
from __future__ import annotations
import csv
import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, r"C:\Collatz")
from result_77_5_compute_R_k import pi_dict

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = Path(r"C:\Collatz")
OUT_MD = OUTDIR / "result_local_invariant_sampling.md"
OUT_V2 = OUTDIR / "result_local_invariant_sampling_v2_dist.csv"
OUT_SIB = OUTDIR / "result_local_invariant_sampling_sibling.csv"
OUT_RT = OUTDIR / "result_local_invariant_sampling_returntime.csv"
OUT_VIS = OUTDIR / "result_local_invariant_sampling_visits.csv"

N_MAX = 10_000_000
MAX_STEPS = 10_000
INT64_GUARD = 1 << 60
ALPHA_BONFERRONI = 0.05 / 5  # five tests


def coprime_residues(M3):
    return [r for r in range(M3) if r % 3 != 0]


def chi2_geometric(observed_v2, k_min=1, k_max=29):
    """Chi² test of observed counts vs geometric P(v_2=k) = 2^-k for k>=1."""
    n_total = int(observed_v2[k_min:k_max+1].sum())
    if n_total < 100:
        return float("nan"), 0, float("nan")
    expected = np.array([n_total * 2.0**(-k) for k in range(k_min, k_max+1)])
    obs = observed_v2[k_min:k_max+1].astype(np.float64)
    # Drop bins with expected < 5 to keep chi² valid
    valid = expected >= 5
    chi2 = float(np.sum((obs[valid] - expected[valid])**2 / expected[valid]))
    df = int(valid.sum()) - 1
    p = float(stats.chi2.sf(chi2, df)) if df > 0 else float("nan")
    return chi2, df, p


def main():
    t0 = time.time()
    print("[1/5] generating odd-coprime-to-3 starts in [1, 10^7]", flush=True)
    arr = np.arange(1, N_MAX + 1, 2, dtype=np.int64)
    arr = arr[arr % 3 != 0]
    # Exclude 1 itself (it's the cycle vertex)
    arr = arr[arr != 1]
    starts = arr
    N_STARTS = len(starts)
    print(f"  {N_STARTS:,} starts", flush=True)

    n_arr = starts.copy()
    starts_mod27 = (starts % 27).astype(np.int64)
    alive = np.ones(N_STARTS, dtype=bool)
    step = np.zeros(N_STARTS, dtype=np.int32)
    return_times = np.zeros(N_STARTS, dtype=np.int32)

    H_v2 = np.zeros(40, dtype=np.int64)
    H_v2_mod27 = np.zeros((30, 27), dtype=np.int64)
    H_v2_mod81 = np.zeros((30, 81), dtype=np.int64)
    H_5_27 = np.zeros((5, 27), dtype=np.int64)
    H_7_27 = np.zeros((7, 27), dtype=np.int64)
    H_11_27 = np.zeros((11, 27), dtype=np.int64)
    H_visits = np.zeros((27, 27), dtype=np.int64)

    print("[2/5] iterating trajectories...", flush=True)
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
        v2 = np.zeros(len(sub), dtype=np.int64)
        m_work = m.copy()
        while True:
            even = (m_work & 1) == 0
            if not even.any():
                break
            v2[even] += 1
            m_work[even] >>= 1
        m_final = m_work

        ns_mod27 = (ns % 27).astype(np.int64)
        ns_mod81 = (ns % 81).astype(np.int64)
        ns_mod5 = (ns % 5).astype(np.int64)
        ns_mod7 = (ns % 7).astype(np.int64)
        ns_mod11 = (ns % 11).astype(np.int64)
        v2_clip40 = np.minimum(v2, 39).astype(np.int64)
        v2_clip30 = np.minimum(v2, 29).astype(np.int64)
        starts_alive_mod27 = starts_mod27[sub]

        H_v2 += np.bincount(v2_clip40, minlength=40)
        H_v2_mod27 += np.bincount(v2_clip30 * 27 + ns_mod27,
                                   minlength=30*27).reshape(30, 27)
        H_v2_mod81 += np.bincount(v2_clip30 * 81 + ns_mod81,
                                   minlength=30*81).reshape(30, 81)
        H_5_27 += np.bincount(ns_mod5 * 27 + ns_mod27,
                              minlength=5*27).reshape(5, 27)
        H_7_27 += np.bincount(ns_mod7 * 27 + ns_mod27,
                              minlength=7*27).reshape(7, 27)
        H_11_27 += np.bincount(ns_mod11 * 27 + ns_mod27,
                               minlength=11*27).reshape(11, 27)
        H_visits += np.bincount(starts_alive_mod27 * 27 + ns_mod27,
                                minlength=27*27).reshape(27, 27)

        n_arr[sub] = m_final
        step[sub] += 1
        new_dead = m_final == 1
        if new_dead.any():
            dead_idx = sub[new_dead]
            return_times[dead_idx] = step[dead_idx]
            alive[dead_idx] = False

        if (s + 1) % 50 == 0:
            n_alive = int(alive.sum())
            print(f"  step {s+1}: alive={n_alive:,}  "
                  f"({100*n_alive/N_STARTS:.2f}%)", flush=True)
            if n_alive == 0:
                break

    still_alive = np.where(alive)[0]
    return_times[still_alive] = step[still_alive]
    print(f"  iteration done in {time.time()-t0:.1f}s. "
          f"still alive after cap: {still_alive.size:,}", flush=True)

    coprime27 = coprime_residues(27)  # 18 classes
    coprime81 = coprime_residues(81)  # 54 classes

    print(f"\n[3/5] running tests...", flush=True)
    results = {}

    # === T1: v_2 distribution vs geometric ===
    print("\nT1: v_2(3n+1) distribution vs geometric null", flush=True)
    chi2_t1, df_t1, p_t1 = chi2_geometric(H_v2)
    total_steps = int(H_v2.sum())
    print(f"  total Syracuse steps recorded: {total_steps:,}", flush=True)
    print(f"  chi² = {chi2_t1:.2f}, df = {df_t1}, p = {p_t1:.4e}", flush=True)
    # Per-k deviations
    n_total = int(H_v2[1:30].sum())
    expected = np.array([n_total * 2.0**(-k) for k in range(1, 30)])
    obs = H_v2[1:30].astype(np.float64)
    rel_dev = (obs - expected) / np.maximum(expected, 1.0)
    print(f"  per-k relative deviation (k=1..10):", flush=True)
    for k in range(1, 11):
        print(f"    k={k:>2}: obs={int(obs[k-1]):>10,}  "
              f"exp={expected[k-1]:>14,.0f}  rel_dev={rel_dev[k-1]:+.4%}",
              flush=True)
    t1_verdict = "STRUCTURE" if p_t1 < ALPHA_BONFERRONI else "NULL"
    results["T1"] = {
        "name": "v_2 distribution vs geometric",
        "chi2": chi2_t1, "df": df_t1, "p": p_t1,
        "verdict": t1_verdict,
        "rel_dev": rel_dev,
    }
    print(f"  -> {t1_verdict}  (Bonferroni p<{ALPHA_BONFERRONI:.4f})",
          flush=True)

    # === T2: v_2 conditional on mod 27, mod 81 ===
    print("\nT2: v_2 conditional on n mod 27 (and mod 81)", flush=True)
    # Pooled homogeneity test (independence between v_2 levels and residues)
    sub27 = H_v2_mod27[1:30][:, coprime27]  # rows v_2=1..29, cols 18 residues
    sub81 = H_v2_mod81[1:30][:, coprime81]  # rows v_2=1..29, cols 54 residues
    # Drop empty rows (very rare v_2 levels)
    nonzero27 = sub27.sum(axis=1) >= 5
    nonzero81 = sub81.sum(axis=1) >= 5
    chi2_27, p_27, df_27, _ = stats.chi2_contingency(sub27[nonzero27])
    chi2_81, p_81, df_81, _ = stats.chi2_contingency(sub81[nonzero81])
    print(f"  mod 27: chi² = {chi2_27:.2f}, df = {df_27}, p = {p_27:.4e}",
          flush=True)
    print(f"  mod 81: chi² = {chi2_81:.2f}, df = {df_81}, p = {p_81:.4e}",
          flush=True)
    t2_verdict = ("STRUCTURE" if min(p_27, p_81) < ALPHA_BONFERRONI
                  else "NULL")
    results["T2"] = {
        "name": "v_2 conditional on mod 27/81 (independence test)",
        "chi2_mod27": chi2_27, "p_mod27": p_27,
        "chi2_mod81": chi2_81, "p_mod81": p_81,
        "verdict": t2_verdict,
    }
    print(f"  -> {t2_verdict}", flush=True)

    # === T3: sibling primes ===
    print("\nT3: sibling-prime joint distributions vs CRT product", flush=True)
    sib_results = {}
    for prime, H in [(5, H_5_27), (7, H_7_27), (11, H_11_27)]:
        table = H[:, coprime27]  # (prime) × 18
        chi2_p, p_p, df_p, _ = stats.chi2_contingency(table)
        sib_results[prime] = (chi2_p, df_p, p_p)
        print(f"  mod {prime} × mod 27: chi² = {chi2_p:.2f}, "
              f"df = {df_p}, p = {p_p:.4e}", flush=True)
    t3_verdict = ("STRUCTURE"
                  if min(r[2] for r in sib_results.values()) < ALPHA_BONFERRONI
                  else "NULL")
    results["T3"] = {
        "name": "sibling primes vs CRT independence",
        "sib": sib_results, "verdict": t3_verdict,
    }
    print(f"  -> {t3_verdict}", flush=True)

    # === T4: return time ===
    print("\nT4: return time τ(n) conditional on starting mod 27", flush=True)
    log_starts = np.log(starts.astype(np.float64))
    rt_float = return_times.astype(np.float64)
    slope_t4, intercept_t4, r_t4, p_lin, _ = stats.linregress(log_starts,
                                                                rt_float)
    print(f"  τ ≈ {slope_t4:.4f}·log(n) + {intercept_t4:.4f}  "
          f"(R = {r_t4:.4f})", flush=True)
    print(f"  Lagarias predict slope = 2/log(4/3) = {2/np.log(4/3):.4f}",
          flush=True)
    residuals_t4 = rt_float - (slope_t4 * log_starts + intercept_t4)
    # ANOVA: do conditional means by mod 27 differ?
    # Use one-way ANOVA F-test
    groups = [residuals_t4[starts_mod27 == r]
              for r in coprime27 if (starts_mod27 == r).any()]
    if len(groups) >= 2 and all(len(g) > 1 for g in groups):
        F_t4, p_t4 = stats.f_oneway(*groups)
    else:
        F_t4 = float("nan"); p_t4 = float("nan")
    print(f"  ANOVA F = {F_t4:.4f}, p = {p_t4:.4e}", flush=True)
    cond_means = {r: float(residuals_t4[starts_mod27 == r].mean())
                  for r in coprime27}
    cond_stds = {r: float(residuals_t4[starts_mod27 == r].std())
                 for r in coprime27}
    cond_counts = {r: int((starts_mod27 == r).sum()) for r in coprime27}
    t4_verdict = "STRUCTURE" if p_t4 < ALPHA_BONFERRONI else "NULL"
    results["T4"] = {
        "name": "return-time τ vs starting mod 27",
        "slope": slope_t4, "intercept": intercept_t4, "R": r_t4,
        "F": F_t4, "p": p_t4, "verdict": t4_verdict,
        "cond_means": cond_means, "cond_stds": cond_stds,
        "cond_counts": cond_counts,
    }
    print(f"  -> {t4_verdict}", flush=True)

    # === T5: visit frequencies vs π_3 ===
    print("\nT5: visit-frequency vs framework's π_3 stationary",
          flush=True)
    pi_3_dict, _ = pi_dict(3)
    pi_3_vec = np.zeros(27, dtype=np.float64)
    for r, p in pi_3_dict.items():
        pi_3_vec[r] = float(p)
    # H_visits[starting_r, current_r]: per starting residue, distribution
    # of visited residues. For each starting r in coprime27 with substantial
    # visit count, test whether empirical matches π_3.
    per_start_chi2 = {}
    for r in coprime27:
        row = H_visits[r, :].astype(np.float64)
        n_r = row.sum()
        if n_r < 100:
            continue
        # restrict to coprime classes
        emp = row[coprime27]
        expected_t5 = pi_3_vec[coprime27] * n_r
        valid = expected_t5 >= 5
        chi2_r = float(np.sum((emp[valid] - expected_t5[valid])**2
                                / expected_t5[valid]))
        df_r = int(valid.sum()) - 1
        p_r = float(stats.chi2.sf(chi2_r, df_r)) if df_r > 0 else float("nan")
        per_start_chi2[r] = (chi2_r, df_r, p_r, int(n_r))

    # Aggregate p (Fisher's combined): -2 Σ log(p) ~ chi²(2k)
    if per_start_chi2:
        all_p = np.array([v[2] for v in per_start_chi2.values()])
        all_p = np.maximum(all_p, 1e-300)  # avoid log(0)
        fisher_stat = -2 * np.sum(np.log(all_p))
        fisher_df = 2 * len(all_p)
        fisher_p = float(stats.chi2.sf(fisher_stat, fisher_df))
    else:
        fisher_p = float("nan")
        fisher_stat = float("nan")
        fisher_df = 0
    n_signif = sum(1 for v in per_start_chi2.values()
                   if v[2] < ALPHA_BONFERRONI / len(per_start_chi2))
    print(f"  per-starting-residue χ² tests (Bonferroni-corrected within "
          f"this test):", flush=True)
    print(f"    {n_signif} of {len(per_start_chi2)} starting residues "
          f"reject π_3 match at p<{ALPHA_BONFERRONI/max(len(per_start_chi2),1):.2e}",
          flush=True)
    print(f"  Fisher combined: stat = {fisher_stat:.2f}, df = {fisher_df}, "
          f"p = {fisher_p:.4e}", flush=True)
    t5_verdict = "STRUCTURE" if fisher_p < ALPHA_BONFERRONI else "NULL"
    results["T5"] = {
        "name": "visit frequency vs π_3",
        "per_start_chi2": per_start_chi2,
        "fisher_stat": fisher_stat, "fisher_df": fisher_df,
        "fisher_p": fisher_p,
        "verdict": t5_verdict,
    }
    print(f"  -> {t5_verdict}", flush=True)

    # ============ Write CSVs ============
    print(f"\n[4/5] writing CSVs", flush=True)
    with open(OUT_V2, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["k", "obs_pooled", "expected_geom", "rel_dev"]
                   + [f"obs_mod27_r={r}" for r in coprime27]
                   + [f"obs_mod81_r={r}" for r in coprime81[:10]])
        for k in range(1, 30):
            row = [k, int(H_v2[k]), float(expected[k-1]),
                   float(rel_dev[k-1])]
            for r in coprime27:
                row.append(int(H_v2_mod27[k, r]))
            for r in coprime81[:10]:
                row.append(int(H_v2_mod81[k, r]))
            w.writerow(row)

    with open(OUT_SIB, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["sibling_prime", "a_mod_p", "r_mod_27", "count"])
        for prime, H in [(5, H_5_27), (7, H_7_27), (11, H_11_27)]:
            for a in range(prime):
                for r in coprime27:
                    w.writerow([prime, a, r, int(H[a, r])])

    with open(OUT_RT, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["start_mod27", "n_starts", "mean_resid_tau",
                    "std_resid_tau"])
        for r in coprime27:
            cm = results["T4"]["cond_means"].get(r, float("nan"))
            cs = results["T4"]["cond_stds"].get(r, float("nan"))
            cn = results["T4"]["cond_counts"].get(r, 0)
            w.writerow([r, cn, f"{cm:.6f}", f"{cs:.6f}"])

    with open(OUT_VIS, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["start_mod27", "current_mod27", "visit_count",
                    "expected_pi3", "obs_normalized"])
        for r in coprime27:
            row_total = float(H_visits[r, :].sum())
            for r2 in coprime27:
                cnt = int(H_visits[r, r2])
                exp = float(pi_3_vec[r2] * row_total)
                obs_n = cnt / row_total if row_total > 0 else 0.0
                w.writerow([r, r2, cnt, f"{exp:.4f}", f"{obs_n:.6f}"])

    # ============ Markdown writeup ============
    print(f"\n[5/5] writing markdown", flush=True)
    L = []
    L.append("# Local-Invariant Soil Sampling on Syracuse Trajectories")
    L.append("")
    L.append(f"Date: 2026-05-05.  Ensemble: all odd integers coprime to 3 in "
             f"[1, 10^7] (excluding the cycle vertex n=1) → "
             f"**{N_STARTS:,} starting points**, iterated under accelerated "
             f"Syracuse map T(n) = (3n+1)/2^v_2(3n+1) until return to 1 or "
             f"{MAX_STEPS} step cap.")
    L.append("")
    L.append(f"Bonferroni-corrected significance threshold: "
             f"p < 0.05/5 = **{ALPHA_BONFERRONI:.4f}** (5 tests).")
    L.append("")
    L.append(f"Total Syracuse steps recorded: {total_steps:,}.")
    L.append(f"Trajectories that hit the step cap (didn't return to 1): "
             f"{still_alive.size:,}.")
    L.append(f"Iteration runtime: {time.time()-t0:.1f} s.")
    L.append("")

    L.append("## Verdict summary")
    L.append("")
    L.append("| test | invariant | null model | statistic | p-value | verdict |")
    L.append("|---|---|---|---:|---:|:---:|")
    L.append(f"| T1 | v_2(3n+1) marginal | geometric P(v_2=k) = 2^(-k) | "
             f"χ²={chi2_t1:.1f}, df={df_t1} | {p_t1:.2e} | "
             f"**{results['T1']['verdict']}** |")
    L.append(f"| T2 | v_2 \\| residue (mod 27) | independence | "
             f"χ²={chi2_27:.1f} | {p_27:.2e} | "
             f"**{results['T2']['verdict']}** |")
    L.append(f"| T2' | v_2 \\| residue (mod 81) | independence | "
             f"χ²={chi2_81:.1f} | {p_81:.2e} | (same verdict) |")
    for prime in (5, 7, 11):
        c, d, p = sib_results[prime]
        L.append(f"| T3.{prime} | (mod {prime}, mod 27) joint | "
                 f"CRT product | χ²={c:.1f}, df={d} | {p:.2e} | "
                 f"**{results['T3']['verdict']}** |")
    L.append(f"| T4 | return-time τ \\| start mod 27 | starting residue "
             f"is process variable | F={F_t4:.2f}, ANOVA | {p_t4:.2e} | "
             f"**{results['T4']['verdict']}** |")
    L.append(f"| T5 | visit-freq vs π_3 stationary | empirical → π_3 | "
             f"Fisher χ²={fisher_stat:.1f}, df={fisher_df} | "
             f"{fisher_p:.2e} | **{results['T5']['verdict']}** |")
    L.append("")

    L.append("## T1 — v_2 marginal vs geometric")
    L.append("")
    L.append("| k | observed | expected (geom) | rel dev |")
    L.append("|---:|---:|---:|---:|")
    for k in range(1, 11):
        L.append(f"| {k} | {int(obs[k-1]):,} | {expected[k-1]:,.0f} | "
                 f"{rel_dev[k-1]:+.3%} |")
    L.append("")
    L.append("(See `result_local_invariant_sampling_v2_dist.csv` for k=1..29.)")
    L.append("")

    L.append("## T2 — v_2 conditional on residue")
    L.append("")
    L.append(f"Test of independence between v_2 levels (1..29, dropping "
             f"empty bins) and residue class (coprime mod 27 / mod 81). "
             f"Under the geometric null with CRT independence, v_2 should "
             f"factor from any mod-3^k residue.")
    L.append("")
    L.append(f"- mod 27: χ² = {chi2_27:.2f}, df = {df_27}, "
             f"p = {p_27:.4e}")
    L.append(f"- mod 81: χ² = {chi2_81:.2f}, df = {df_81}, "
             f"p = {p_81:.4e}")
    L.append("")

    L.append("## T3 — sibling primes")
    L.append("")
    L.append("| prime p | χ² | df | p | verdict-this-pair |")
    L.append("|---:|---:|---:|---:|:---:|")
    for prime in (5, 7, 11):
        c, d, p = sib_results[prime]
        v = "STRUCTURE" if p < ALPHA_BONFERRONI else "NULL"
        L.append(f"| {prime} | {c:.2f} | {d} | {p:.4e} | {v} |")
    L.append("")

    L.append("## T4 — return-time conditional on starting residue")
    L.append("")
    L.append(f"τ(n) ≈ {slope_t4:.4f}·log(n) + {intercept_t4:.4f}  "
             f"(Pearson R = {r_t4:.4f})")
    L.append(f"Lagarias prediction: 2/log(4/3) = {2/np.log(4/3):.4f}")
    L.append("")
    L.append(f"ANOVA on residuals by starting mod 27 class: "
             f"F = {F_t4:.4f}, p = {p_t4:.4e}.")
    L.append("")
    L.append("| start mod 27 | N_starts | mean(τ - τ_pred) | std |")
    L.append("|---:|---:|---:|---:|")
    for r in coprime27:
        cm = cond_means.get(r, float("nan"))
        cs = cond_stds.get(r, float("nan"))
        cn = cond_counts.get(r, 0)
        L.append(f"| {r} | {cn:,} | {cm:+.4f} | {cs:.4f} |")
    L.append("")

    L.append("## T5 — visit frequencies vs π_3")
    L.append("")
    L.append(f"Fisher combined p over 18 starting residues: "
             f"{fisher_p:.4e}.")
    L.append("")
    L.append("| start mod 27 | visits total | χ² vs π_3 | p |")
    L.append("|---:|---:|---:|---:|")
    for r in coprime27:
        if r in per_start_chi2:
            c, d, p, n = per_start_chi2[r]
            L.append(f"| {r} | {n:,} | {c:.2f} | {p:.4e} |")
    L.append("")

    L.append("## Interpretation")
    L.append("")
    null_count = sum(1 for k in ("T1", "T2", "T3", "T4", "T5")
                     if results[k]["verdict"] == "NULL")
    struct_count = 5 - null_count
    L.append(f"Of 5 tests: **{null_count} NULL, {struct_count} STRUCTURE** "
             f"(Bonferroni α = {ALPHA_BONFERRONI:.4f}).")
    L.append("")

    if struct_count == 0:
        L.append("All marginal local invariants factor cleanly. The "
                 "framework's global tools (Plancherel mass, conservation "
                 "laws) are well-suited to the problem; no detectable "
                 "structure beyond first moments lives at the per-step or "
                 "per-residue local level. Negative result, but informative: "
                 "rules out hidden local arithmetic that the global "
                 "treatment would average over.")
    else:
        L.append("STRUCTURE was detected in some tests. See per-test sections "
                 "above for effect sizes and per-residue/per-prime patterns. "
                 "Note: with N ~ 10^7 starts and ~10^8 total steps, even "
                 "small relative deviations from the null can register at "
                 "very small p — interpret effect sizes (relative deviations "
                 "and per-residue magnitudes), not just p-values, before "
                 "claiming structural significance.")
    L.append("")
    L.append("**Caveat on power.** This test has very high power to detect "
             "small deviations (large N). A STRUCTURE flag at p < 0.01 may "
             "represent a relative effect of 0.1–1% that doesn't carry "
             "useful arithmetic information. The effect-size tables above "
             "should be inspected before drawing structural conclusions.")
    L.append("")

    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    print(f"\n[done] runtime = {time.time()-t0:.1f}s; outputs in {OUTDIR}",
          flush=True)


if __name__ == "__main__":
    main()
