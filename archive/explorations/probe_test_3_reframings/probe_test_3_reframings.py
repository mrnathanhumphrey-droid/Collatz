"""
probe_test_3_reframings.py
==========================
Three sub-probes (4A, 4B, 4C) reframing the original Test 3 (bit-budget
admissibility, which returned null with low-v_2 → larger horizons).

Working definition (concrete, computable):
  - Sample N = 500,000 odd starts coprime to 3 in [3, 10^6].
  - For each n: iterate Syracuse to length L (orbit to 1), record:
      L(n)                — orbit length
      n_d mod 2^k for d=0..D, k=5,6,7
      n_d mod 3^j for d=0..D, j=2,3,4
      a_final at modulus 2^k via deterministic prefix algorithm
  - 4A: stratify orbit-length by (mod 2^k, mod 3^j) joint cells, compute
        rho(r,s,k,j) = H_joint / min(H_2, H_3); aggregate by a_final.
  - 4B: per (r,k), count distinct (D-step) orbit traces mod 2^k across
        starts ≡ r mod 2^k; gamma(r,k) = log(N_distinct) / D.
  - 4C: per (r,k), P(r,k) = number of starts in sample whose D-step orbit
        visits residue r at any depth; sigma(r,k) = E_P / P(r,k).
  - 4-cross: correlation of rho, gamma, sigma across a_final classes.

Pre-registered gates per sub-probe (per the brief):
  4A — A: rho varies by a_final >2× spread → structural.
       B: uniform within ~10% → close.
       C: tracks v_2 not a_final → null in disguise.
  4B — A: gamma varies by a_final >30% spread → discriminates.
       B: uniform within ~10% → close.
       C: correlates with v_2 not a_final → null.
  4C — A: sigma varies by a_final >2× spread → stuck classes.
       B: uniform within ~30% → close.
       C: correlates with v_2 → null.
"""
from __future__ import annotations

import csv
import math
import os
import sys
import time
from collections import defaultdict

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz\probe_test_3_reframings"
N_SAMPLE = 500_000          # number of starting integers to sample
DEPTH_TRACE = 10            # depth for 4B trace counting
DEPTH_PRED = 5              # depth for 4C predecessor counting
MAX_STEPS = 1000
KS = [5, 6, 7]              # mod 2^k
JS = [2, 3, 4]              # mod 3^j


def deterministic_prefix(r: int, a0: int, max_steps: int = 400):
    """Symbolic prefix algorithm: starting from (a, c) = (a0, r), iterate
    until a is odd. Returns (steps, a_final, c_final)."""
    a, c = a0, r
    steps = 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2
            c //= 2
        else:
            a *= 3
            c = 3 * c + 1
        steps += 1
    return steps, a, c


def main():
    print("=" * 78)
    print("Probe Test 3 reframings (4A joint admissibility, 4B residue-counting, "
          "4C stuck residues)")
    print("=" * 78)

    # ---------------- 1. Sample + trajectory iteration ----------------
    t0 = time.time()
    print(f"[1] Sampling N={N_SAMPLE:,} odd-coprime-to-3 starts in [3, 10^6]...",
          flush=True)
    rng = np.random.default_rng(seed=20260505)
    candidates = np.arange(3, 10_000_001, 2, dtype=np.int64)
    candidates = candidates[candidates % 3 != 0]
    if len(candidates) > N_SAMPLE:
        idx = rng.choice(len(candidates), size=N_SAMPLE, replace=False)
        starts = np.sort(candidates[idx])
    else:
        starts = candidates
    N = len(starts)
    print(f"  {N:,} starts; t={time.time()-t0:.1f}s", flush=True)

    # Determine max k and j for residue tracking
    max_2_mod = 1 << max(KS)   # 128 for KS = [5,6,7]
    max_3_mod = 3 ** max(JS)   # 81 for JS = [2,3,4]

    # Track residues at depth d=0,1,...,D for max(DEPTH_TRACE, DEPTH_PRED)
    D = max(DEPTH_TRACE, DEPTH_PRED)
    print(f"[2] Iterating Syracuse, recording residues at depth 0..{D} "
          f"and orbit length L(n)", flush=True)

    n_arr = starts.copy()
    alive = np.ones(N, dtype=bool)
    length = np.zeros(N, dtype=np.int32)
    # residues_2[d, i] = n at depth d for start i (mod max_2_mod)
    # residues_3[d, i] = n at depth d for start i (mod max_3_mod)
    residues_2 = np.zeros((D + 1, N), dtype=np.int64)
    residues_3 = np.zeros((D + 1, N), dtype=np.int64)
    residues_2[0] = starts % max_2_mod
    residues_3[0] = starts % max_3_mod

    for s in range(MAX_STEPS):
        sub = np.where(alive)[0]
        if sub.size == 0:
            break
        ns = n_arr[sub]
        m = 3 * ns + 1
        # Strip 2-adic factors via vectorized parity loop
        m_work = m.copy()
        while True:
            even = (m_work & 1) == 0
            if not even.any():
                break
            m_work[even] >>= 1
        n_arr[sub] = m_work
        # Record residues at next depth (s+1) if within D
        if s + 1 <= D:
            residues_2[s + 1, sub] = m_work % max_2_mod
            residues_3[s + 1, sub] = m_work % max_3_mod
            # For dead-already entries (those still alive go on; for alive ones
            # we recorded; for unreached entries beyond MAX_STEPS we'd skip)
            alive_idx = sub
            # but for entries that just became 1, freeze residues
        # Mark dead
        new_dead = m_work == 1
        if new_dead.any():
            dead = sub[new_dead]
            length[dead] = s + 1
            alive[dead] = False
            # Freeze remaining depth slots at 1 for dead trajectories
            for d_future in range(s + 2, D + 1):
                residues_2[d_future, dead] = 1
                residues_3[d_future, dead] = 1
        if (s + 1) % 100 == 0:
            print(f"  step {s+1}: alive {int(alive.sum()):,}", flush=True)
            if alive.sum() == 0:
                break

    # Any survivors? Mark length = MAX_STEPS
    survivors = np.where(alive)[0]
    length[survivors] = MAX_STEPS
    print(f"  iteration done t={time.time()-t0:.1f}s, survivors={len(survivors)}",
          flush=True)

    # ---------------- 3. a_final at each modulus 2^k ----------------
    print(f"[3] Computing a_final via prefix algorithm at each k...", flush=True)
    t1 = time.time()
    a_final_table = {}  # k → array of a_final per start
    for k in KS:
        a0 = 1 << k
        a_finals = np.zeros(N, dtype=np.int64)
        # Prefix only depends on r = n mod 2^k, so cache by residue
        r_table = (starts % a0).astype(np.int64)
        cache = {}
        for r in np.unique(r_table):
            r_int = int(r)
            _, af, _ = deterministic_prefix(r_int, a0)
            cache[r_int] = af
        for i in range(N):
            a_finals[i] = cache[int(r_table[i])]
        a_final_table[k] = a_finals
        unique_a = sorted(np.unique(a_finals).tolist())
        print(f"  k={k}: a_final values = {unique_a}", flush=True)
    print(f"  t={time.time()-t1:.1f}s", flush=True)

    # ============================================================
    # 4A: Joint admissibility — orbit length stratified by joint cells
    # ============================================================
    print(f"\n[4A] Joint admissibility by (mod 2^k, mod 3^j)", flush=True)
    rows_4A = []
    L_arr = length.astype(np.float64)
    for k in KS:
        m2 = 1 << k
        r2 = starts % m2
        for j in JS:
            m3 = 3 ** j
            r3 = starts % m3
            # H_2(r) = mean L per r2 class
            H2 = {}
            for r in range(m2):
                if r % 2 == 0:
                    continue
                mask = r2 == r
                if mask.sum() < 5:
                    continue
                H2[r] = float(L_arr[mask].mean())
            # H_3(s) = mean L per r3 class
            H3 = {}
            for s in range(m3):
                if s % 3 == 0:
                    continue
                mask = r3 == s
                if mask.sum() < 5:
                    continue
                H3[s] = float(L_arr[mask].mean())
            # H_joint(r, s) = mean L per joint cell
            for r in H2:
                for s in H3:
                    mask = (r2 == r) & (r3 == s)
                    if mask.sum() < 5:
                        continue
                    Hj = float(L_arr[mask].mean())
                    rho = Hj / min(H2[r], H3[s])
                    af = int(a_final_table[k][np.where(r2 == r)[0][0]])
                    rows_4A.append({
                        "r2": int(r), "k": k, "s3": int(s), "j": j,
                        "n_cell": int(mask.sum()),
                        "H_2_r": H2[r], "H_3_s": H3[s], "H_joint": Hj,
                        "rho": rho, "a_final": af,
                    })
        print(f"  k={k}: {sum(1 for r in rows_4A if r['k'] == k)} cells", flush=True)

    # Aggregate rho by a_final per k
    aggregate_4A = {}
    for k in KS:
        cells = [r for r in rows_4A if r["k"] == k]
        by_af = defaultdict(list)
        for r in cells:
            by_af[r["a_final"]].append(r["rho"])
        aggregate_4A[k] = {af: (np.mean(rs), np.std(rs), len(rs))
                           for af, rs in by_af.items()}

    print(f"\n  4A aggregate rho by a_final:")
    for k in KS:
        print(f"    k={k}:")
        for af, (m_, sd, n_) in sorted(aggregate_4A[k].items()):
            print(f"      a_final={af}: rho mean={m_:.6f}, sd={sd:.6f}, n={n_}",
                  flush=True)

    # ============================================================
    # 4B: Residue-counting — distinct D-step traces mod 2^k per r class
    # ============================================================
    print(f"\n[4B] Residue counting (depth {DEPTH_TRACE} traces mod 2^k)",
          flush=True)
    rows_4B = []
    for k in KS:
        m2 = 1 << k
        # For each starting residue r mod 2^k, collect all distinct
        # (n_0 mod 2^k, n_1 mod 2^k, ..., n_DEPTH_TRACE mod 2^k) traces
        # over starts in our sample.
        traces_2 = (residues_2[:DEPTH_TRACE + 1] % m2).astype(np.int64)  # shape (D+1, N)
        # Pack each column into a hashable bytes signature
        sig = traces_2.T.tobytes()
        # Compute per-start signatures by row
        sig_rows = traces_2.T.copy()
        # Group by starting residue r
        r_starts = sig_rows[:, 0]
        # For each unique r, count distinct full traces
        for r in np.unique(r_starts):
            r_int = int(r)
            if r_int % 2 == 0 or r_int == 0:
                continue
            mask = r_starts == r_int
            if mask.sum() < 5:
                continue
            traces_in_class = sig_rows[mask]
            # Count distinct rows
            distinct = np.unique(traces_in_class, axis=0)
            n_distinct = len(distinct)
            n_total = mask.sum()
            # gamma = log(n_distinct) / DEPTH_TRACE
            gamma = math.log(n_distinct) / DEPTH_TRACE if n_distinct > 0 else 0
            af = int(a_final_table[k][np.where(mask)[0][0]])
            rows_4B.append({
                "r": r_int, "k": k, "n_total": int(n_total),
                "n_distinct": int(n_distinct), "gamma": gamma, "a_final": af,
            })
        print(f"  k={k}: {sum(1 for r in rows_4B if r['k']==k)} residue classes",
              flush=True)

    aggregate_4B = {}
    for k in KS:
        cells = [r for r in rows_4B if r["k"] == k]
        by_af = defaultdict(list)
        for r in cells:
            by_af[r["a_final"]].append(r["gamma"])
        aggregate_4B[k] = {af: (np.mean(rs), np.std(rs), len(rs))
                           for af, rs in by_af.items()}
    print(f"\n  4B aggregate gamma by a_final:")
    for k in KS:
        print(f"    k={k}:")
        for af, (m_, sd, n_) in sorted(aggregate_4B[k].items()):
            print(f"      a_final={af}: gamma mean={m_:.6f}, sd={sd:.6f}, n={n_}",
                  flush=True)

    # ============================================================
    # 4C: Stuck residues — predecessor count via depth-D orbit visits
    # ============================================================
    print(f"\n[4C] Stuck-ness (depth {DEPTH_PRED} orbit visits per residue)",
          flush=True)
    rows_4C = []
    for k in KS:
        m2 = 1 << k
        # For each residue r mod 2^k, count distinct starts whose orbit
        # at depth 0..DEPTH_PRED visits residue r at any of those depths.
        traces = (residues_2[:DEPTH_PRED + 1] % m2).astype(np.int64)  # (D+1, N)
        # For each r, mask = any d in 0..D has traces[d,:] == r → boolean per start
        odd_residues = [r for r in range(m2) if r % 2 == 1]
        P_table = {}
        for r in odd_residues:
            visits = (traces == r).any(axis=0)  # (N,) bool
            P_table[r] = int(visits.sum())
        # Expected predecessor count = mean over residues
        E_P = float(np.mean(list(P_table.values())))
        for r, P in P_table.items():
            sigma = E_P / P if P > 0 else float("inf")
            # a_final from prefix on r
            _, af, _ = deterministic_prefix(r, m2)
            rows_4C.append({
                "r": r, "k": k, "P_r": P, "E_P": E_P, "sigma": sigma,
                "a_final": int(af),
            })
        print(f"  k={k}: {len(odd_residues)} residues, E_P={E_P:.1f}", flush=True)

    aggregate_4C = {}
    for k in KS:
        cells = [r for r in rows_4C if r["k"] == k]
        by_af = defaultdict(list)
        for r in cells:
            by_af[r["a_final"]].append(r["sigma"])
        aggregate_4C[k] = {af: (np.mean(rs), np.std(rs), len(rs))
                           for af, rs in by_af.items()}
    print(f"\n  4C aggregate sigma by a_final:")
    for k in KS:
        print(f"    k={k}:")
        for af, (m_, sd, n_) in sorted(aggregate_4C[k].items()):
            print(f"      a_final={af}: sigma mean={m_:.6f}, sd={sd:.6f}, n={n_}",
                  flush=True)

    # ============================================================
    # Cross-correlation across a_final classes
    # ============================================================
    print(f"\n[X] Cross-correlation of rho, gamma, sigma across a_final classes",
          flush=True)
    cross_rows = []
    for k in KS:
        # Common a_final classes across all three sub-probes at this k
        afs_4A = set(aggregate_4A[k].keys())
        afs_4B = set(aggregate_4B[k].keys())
        afs_4C = set(aggregate_4C[k].keys())
        common = sorted(afs_4A & afs_4B & afs_4C)
        rho_vec = [aggregate_4A[k][af][0] for af in common]
        gamma_vec = [aggregate_4B[k][af][0] for af in common]
        sigma_vec = [aggregate_4C[k][af][0] for af in common]
        if len(common) >= 3:
            r_rg = np.corrcoef(rho_vec, gamma_vec)[0, 1]
            r_rs = np.corrcoef(rho_vec, sigma_vec)[0, 1]
            r_gs = np.corrcoef(gamma_vec, sigma_vec)[0, 1]
        else:
            r_rg = r_rs = r_gs = float("nan")
        for af, rho, g, s in zip(common, rho_vec, gamma_vec, sigma_vec):
            cross_rows.append({
                "k": k, "a_final": int(af), "rho_mean": rho,
                "gamma_mean": g, "sigma_mean": s,
            })
        print(f"  k={k}: a_final classes = {common}")
        print(f"    rho range: [{min(rho_vec):.4f}, {max(rho_vec):.4f}], "
              f"spread = {(max(rho_vec)-min(rho_vec))/min(rho_vec):.4f}")
        print(f"    gamma range: [{min(gamma_vec):.4f}, {max(gamma_vec):.4f}], "
              f"spread = {(max(gamma_vec)-min(gamma_vec))/max(min(gamma_vec), 1e-10):.4f}")
        print(f"    sigma range: [{min(sigma_vec):.4f}, {max(sigma_vec):.4f}], "
              f"spread = {(max(sigma_vec)-min(sigma_vec))/min(sigma_vec):.4f}")
        print(f"    corr(rho, gamma) = {r_rg:+.4f}, "
              f"corr(rho, sigma) = {r_rs:+.4f}, "
              f"corr(gamma, sigma) = {r_gs:+.4f}", flush=True)

    # ============================================================
    # Write all CSVs
    # ============================================================
    def write_csv(path, rows, cols):
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in rows:
                w.writerow({c: (f"{v:.10e}" if isinstance(v, float) else v)
                            for c, v in r.items()})

    write_csv(os.path.join(OUT_DIR, "result_4A_joint_admissibility.csv"),
              rows_4A, ["r2", "k", "s3", "j", "n_cell", "H_2_r", "H_3_s",
                        "H_joint", "rho", "a_final"])
    write_csv(os.path.join(OUT_DIR, "result_4B_residue_counting.csv"),
              rows_4B, ["r", "k", "n_total", "n_distinct", "gamma", "a_final"])
    write_csv(os.path.join(OUT_DIR, "result_4C_stuck_residues.csv"),
              rows_4C, ["r", "k", "P_r", "E_P", "sigma", "a_final"])
    write_csv(os.path.join(OUT_DIR, "result_4_cross_correlation.csv"),
              cross_rows, ["k", "a_final", "rho_mean", "gamma_mean", "sigma_mean"])

    # ============================================================
    # Findings markdown with verdicts
    # ============================================================
    md = []
    md.append("# Test 3 reframings — joint admissibility, residue counting, stuck residues")
    md.append("")
    md.append(f"Date 2026-05-05. Sample N={N_SAMPLE:,} odd-coprime-to-3 starts in "
              f"[3, 10^6]; depth D=max({DEPTH_TRACE}, {DEPTH_PRED})={D} for "
              f"residue tracking; orbit length L(n) capped at {MAX_STEPS}.")
    md.append("")
    md.append("## Working definitions")
    md.append("")
    md.append("- **H_2(r, k)** = mean orbit-length L over starts with n ≡ r (mod 2^k).")
    md.append("- **H_3(s, j)** = mean orbit-length L over starts with n ≡ s (mod 3^j).")
    md.append("- **H_joint(r, s, k, j)** = mean L over starts with n ≡ r mod 2^k AND "
              "n ≡ s mod 3^j.")
    md.append("- **rho(r, s, k, j)** = H_joint / min(H_2(r, k), H_3(s, j)). "
              "Joint reduction factor.")
    md.append("- **gamma(r, k)** = log(N_distinct(D)) / D, where N_distinct(D) is the "
              "number of distinct depth-D residue traces (mod 2^k) among starts ≡ r mod 2^k.")
    md.append("- **P(r, k)** = number of starts whose D-step orbit visits residue r mod 2^k "
              "at any depth 0..D.")
    md.append("- **sigma(r, k)** = E_P / P(r, k). Stuck-ness: large σ ⇒ rarely visited.")
    md.append("")

    md.append("## Sub-probe 4A — joint admissibility")
    md.append("")
    md.append("rho aggregated by a_final, per k:")
    md.append("")
    md.append("| k | a_final | mean rho | sd | n_cells |")
    md.append("|---|---|---|---|---|")
    for k in KS:
        for af, (m_, sd, n_) in sorted(aggregate_4A[k].items()):
            md.append(f"| {k} | {af} | {m_:.6f} | {sd:.6f} | {n_} |")
    md.append("")

    # Verdict 4A
    spread_4A = {}
    for k in KS:
        means = [v[0] for v in aggregate_4A[k].values()]
        if len(means) >= 2:
            spread_4A[k] = (max(means) - min(means)) / min(means)
        else:
            spread_4A[k] = 0.0
    max_spread_A = max(spread_4A.values()) if spread_4A else 0.0
    if max_spread_A > 1.0:
        verdict_4A = (f"Gate A — rho spread by a_final exceeds 2× at some k "
                      f"(max spread {max_spread_A:.3f}). Joint admissibility may "
                      f"discriminate a_final classes; follow-up warranted.")
    elif max_spread_A < 0.10:
        verdict_4A = (f"Gate B — rho uniform across a_final classes within ~10% "
                      f"(max spread {max_spread_A:.3f}). This reframing does not "
                      f"discriminate; closes.")
    else:
        verdict_4A = (f"Intermediate — rho varies by {max_spread_A:.3f} (between "
                      f"10% and 100% spread). Worth follow-up correlation with "
                      f"v_2 to disambiguate Gate C from genuine variation.")
    md.append(f"**4A verdict:** {verdict_4A}")
    md.append("")

    md.append("## Sub-probe 4B — residue-counting growth rate")
    md.append("")
    md.append("gamma = log(N_distinct) / D aggregated by a_final, per k:")
    md.append("")
    md.append("| k | a_final | mean gamma | sd | n_residue_classes |")
    md.append("|---|---|---|---|---|")
    for k in KS:
        for af, (m_, sd, n_) in sorted(aggregate_4B[k].items()):
            md.append(f"| {k} | {af} | {m_:.6f} | {sd:.6f} | {n_} |")
    md.append("")

    spread_4B = {}
    for k in KS:
        means = [v[0] for v in aggregate_4B[k].values()]
        if len(means) >= 2 and max(means) > 0:
            spread_4B[k] = (max(means) - min(means)) / max(means)
        else:
            spread_4B[k] = 0.0
    max_spread_B = max(spread_4B.values()) if spread_4B else 0.0
    if max_spread_B > 0.30:
        verdict_4B = (f"Gate A — gamma spread by a_final exceeds 30% at some k "
                      f"(max spread {max_spread_B:.3f}). Residue-counting "
                      f"discriminates classes; structural finding.")
    elif max_spread_B < 0.10:
        verdict_4B = (f"Gate B — gamma uniform across a_final classes within ~10% "
                      f"(max spread {max_spread_B:.3f}). Reframing does not "
                      f"discriminate; closes.")
    else:
        verdict_4B = (f"Intermediate — gamma spread {max_spread_B:.3f} (10-30%). "
                      f"Modest discrimination; not structural at the brief's threshold.")
    md.append(f"**4B verdict:** {verdict_4B}")
    md.append("")

    md.append("## Sub-probe 4C — stuck residues")
    md.append("")
    md.append("sigma = E_P / P(r,k) aggregated by a_final, per k:")
    md.append("")
    md.append("| k | a_final | mean sigma | sd | n_residues |")
    md.append("|---|---|---|---|---|")
    for k in KS:
        for af, (m_, sd, n_) in sorted(aggregate_4C[k].items()):
            md.append(f"| {k} | {af} | {m_:.6f} | {sd:.6f} | {n_} |")
    md.append("")

    spread_4C = {}
    for k in KS:
        means = [v[0] for v in aggregate_4C[k].values()]
        if len(means) >= 2 and min(means) > 0:
            spread_4C[k] = (max(means) - min(means)) / min(means)
        else:
            spread_4C[k] = 0.0
    max_spread_C = max(spread_4C.values()) if spread_4C else 0.0
    if max_spread_C > 1.0:
        verdict_4C = (f"Gate A — sigma spread by a_final exceeds 2× at some k "
                      f"(max spread {max_spread_C:.3f}). Some classes are stuck; "
                      f"structural finding.")
    elif max_spread_C < 0.30:
        verdict_4C = (f"Gate B — sigma uniform across a_final classes within ~30% "
                      f"(max spread {max_spread_C:.3f}). Predecessor structure is "
                      f"class-uniform; reframing closes.")
    else:
        verdict_4C = (f"Intermediate — sigma spread {max_spread_C:.3f} (30-100%). "
                      f"Modest variation; not structural at the brief's 2× threshold.")
    md.append(f"**4C verdict:** {verdict_4C}")
    md.append("")

    md.append("## Cross-correlation across a_final classes")
    md.append("")
    md.append("For each k, correlations between (rho_mean, gamma_mean, sigma_mean) "
              "computed across the common set of a_final classes. If multiple "
              "sub-probes return Gate A, correlations indicate whether they're "
              "measuring the same underlying admissibility-precariousness.")
    md.append("")
    md.append("| k | a_final | rho_mean | gamma_mean | sigma_mean |")
    md.append("|---|---|---|---|---|")
    for r in cross_rows:
        md.append(f"| {r['k']} | {r['a_final']} | {r['rho_mean']:.6f} | "
                  f"{r['gamma_mean']:.6f} | {r['sigma_mean']:.6f} |")
    md.append("")
    md.append("Per-k correlations:")
    md.append("")
    md.append("| k | corr(rho, gamma) | corr(rho, sigma) | corr(gamma, sigma) |")
    md.append("|---|---|---|---|")
    for k in KS:
        afs_4A = set(aggregate_4A[k].keys())
        afs_4B = set(aggregate_4B[k].keys())
        afs_4C = set(aggregate_4C[k].keys())
        common = sorted(afs_4A & afs_4B & afs_4C)
        if len(common) < 3:
            md.append(f"| {k} | — (only {len(common)} common classes) | — | — |")
            continue
        rho_vec = [aggregate_4A[k][af][0] for af in common]
        gamma_vec = [aggregate_4B[k][af][0] for af in common]
        sigma_vec = [aggregate_4C[k][af][0] for af in common]
        r_rg = np.corrcoef(rho_vec, gamma_vec)[0, 1]
        r_rs = np.corrcoef(rho_vec, sigma_vec)[0, 1]
        r_gs = np.corrcoef(gamma_vec, sigma_vec)[0, 1]
        md.append(f"| {k} | {r_rg:+.4f} | {r_rs:+.4f} | {r_gs:+.4f} |")
    md.append("")

    md.append("## Combined verdict")
    md.append("")
    n_gate_A = sum(1 for v in [verdict_4A, verdict_4B, verdict_4C]
                   if v.startswith("Gate A"))
    n_gate_B = sum(1 for v in [verdict_4A, verdict_4B, verdict_4C]
                   if v.startswith("Gate B"))
    md.append(f"Sub-probes returning Gate A: {n_gate_A}/3")
    md.append(f"Sub-probes returning Gate B: {n_gate_B}/3")
    md.append("")
    if n_gate_A >= 2:
        md.append("**At least 2 of 3 sub-probes show Gate-A spread; "
                  "cross-correlation analysis is the next gate.**")
    elif n_gate_B >= 2:
        md.append("**At least 2 of 3 sub-probes return Gate B (uniform). "
                  "The Test 3 reframings, like the original bit-budget Test 3, "
                  "do not produce structural admissibility constraints — "
                  "the framework's admissibility is class-uniform at the "
                  "scales tested.**")
    else:
        md.append("**Mixed results.** Some sub-probes show modest variation, "
                  "others are uniform; the reframings produce intermediate "
                  "signals not strong enough to elevate to Gate-A structural "
                  "finding under the pre-registered thresholds.")
    md.append("")
    md.append("Spread summary:")
    md.append("- 4A (joint admissibility) max spread by a_final: "
              f"{max_spread_A:.4f}  (Gate A threshold: 1.0)")
    md.append("- 4B (residue counting) max spread by a_final: "
              f"{max_spread_B:.4f}  (Gate A threshold: 0.30)")
    md.append("- 4C (stuck residues) max spread by a_final: "
              f"{max_spread_C:.4f}  (Gate A threshold: 1.0)")
    md.append("")

    md.append("## Files")
    md.append("")
    md.append("- `result_4A_joint_admissibility.csv` — per-cell rho with a_final tag")
    md.append("- `result_4B_residue_counting.csv` — per-residue gamma with a_final tag")
    md.append("- `result_4C_stuck_residues.csv` — per-residue sigma with a_final tag")
    md.append("- `result_4_cross_correlation.csv` — per-class summary of all three")
    md.append("- `test_3_reframings_findings.md` — this writeup")

    md_path = os.path.join(OUT_DIR, "test_3_reframings_findings.md")
    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(md))

    print(f"\n[md: {md_path}]")
    print(f"\nTotal runtime: {time.time()-t0:.1f}s")
    print("Done.")


if __name__ == "__main__":
    main()
