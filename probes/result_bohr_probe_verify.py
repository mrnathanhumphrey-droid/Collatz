"""
Verification follow-up to result_bohr_probe.py / result_bohr_probe_cliff.py.

Five tests at the (a=5, b=4) joint cell, N = 1,000,000:
  A. Survival count and fraction at each k.
  B. Controlled subset S_15: trajectories alive at k=15. Compute chi²/d.f. on
     this subset (members still alive at later k) at k = 15, 18, 20, 22, 25, 30.
  C. Geometric fit chi²(k) = chi²_0 · r^k on k ∈ {18, 19, 20, 22, 25, 30}.
  D. Saturation test: extrapolate; if predicted > 1000 at k=50, run k=40, 50.
  E. Top 10 deviating (r_2, r_3) cells at k=30 vs k=20.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

SEED = 20260504
N = 1_000_000
DEPTHS_SAMPLE = [5, 10, 15, 18, 19, 20, 22, 25, 30, 40, 50]
DEPTHS_FIT = [18, 19, 20, 22, 25, 30]
A, B = 5, 4
M2 = 1 << A         # 32
M3 = 3 ** B         # 81
ROWS = 1 << (A - 1) # 16
COLS = 2 * 3 ** (B - 1)  # 54
DF = (ROWS - 1) * (COLS - 1)
S15_DEPTH = 15

OUTDIR = Path(r"C:\Collatz")
CSV_PATH = OUTDIR / "result_bohr_probe_verify.csv"
TOP_PATH = OUTDIR / "result_bohr_probe_verify_top_cells.csv"
MD_PATH = OUTDIR / "result_bohr_probe_verify.md"


def coprime_residues_mod_3b(b):
    return np.array([r for r in range(3 ** b) if r % 3 != 0], dtype=np.int64)


def odd_residues_mod_2a(a):
    return np.array([r for r in range(1 << a) if r % 2 == 1], dtype=np.int64)


def histogram(vals, inv3):
    r2 = (vals & (M2 - 1)).astype(np.int64)
    row = (r2 >> 1)
    r3 = (vals % M3).astype(np.int64)
    col = inv3[r3]
    flat = row * COLS + col
    return np.bincount(flat, minlength=ROWS * COLS).reshape(ROWS, COLS)


def chi2_df_from_H(H, N_eff):
    if N_eff == 0:
        return 0.0, float("nan"), float("nan")
    P = H.astype(np.float64) / N_eff
    Mr = P.sum(axis=1, keepdims=True)
    Mc = P.sum(axis=0, keepdims=True)
    Q = Mr * Mc
    with np.errstate(divide="ignore", invalid="ignore"):
        contrib = np.where(Q > 0, N_eff * (P - Q) ** 2 / Q, 0.0)
    chi2 = float(contrib.sum())
    chi2_per_df = chi2 / DF
    z = (chi2 - DF) / np.sqrt(2.0 * DF)
    return chi2, chi2_per_df, z


def top_cells(H, N_eff, k_count=10):
    P = H.astype(np.float64) / N_eff
    Mr = P.sum(axis=1, keepdims=True)
    Mc = P.sum(axis=0, keepdims=True)
    Q = Mr * Mc
    D = P - Q
    r2_list = odd_residues_mod_2a(A)
    r3_list = coprime_residues_mod_3b(B)
    flat = np.argsort(np.abs(D).ravel())[::-1][:k_count]
    out = []
    for fi in flat:
        i = fi // COLS
        j = fi % COLS
        out.append({
            "r2": int(r2_list[i]),
            "r3": int(r3_list[j]),
            "P": float(P[i, j]),
            "Q": float(Q[i, j]),
            "D": float(D[i, j]),
            "z_cell": float(D[i, j] * np.sqrt(N_eff / max(Q[i, j], 1e-30))),
        })
    return out


def main():
    t0 = time.time()
    print(f"[t={time.time()-t0:6.1f}s] verify probe; N={N:,}; "
          f"depths={DEPTHS_SAMPLE}", flush=True)

    inv3 = -np.ones(M3, dtype=np.int64)
    cnt = 0
    for r in range(M3):
        if r % 3 != 0:
            inv3[r] = cnt
            cnt += 1

    rng = np.random.default_rng(seed=SEED)
    x = rng.integers(low=0, high=500_000_000_000, size=N, dtype=np.int64)
    n_arr = (2 * x + 1).astype(np.int64)
    alive = np.ones(N, dtype=bool)
    in_S15 = None  # set when depth=15 completes

    H_full = {k: np.zeros((ROWS, COLS), dtype=np.int64) for k in DEPTHS_SAMPLE}
    H_s15 = {k: np.zeros((ROWS, COLS), dtype=np.int64) for k in DEPTHS_SAMPLE
             if k >= S15_DEPTH}
    n_alive_full = {k: 0 for k in DEPTHS_SAMPLE}
    n_alive_s15 = {k: 0 for k in DEPTHS_SAMPLE if k >= S15_DEPTH}
    n_overflow = 0

    depth_set = set(DEPTHS_SAMPLE)
    depth_max = max(DEPTHS_SAMPLE)

    # Threshold to switch a trajectory to Python int
    INT64_GUARD = 1 << 60

    overflow_idx = set()  # indices using Python int
    py_vals = {}          # idx -> Python int

    for depth in range(1, depth_max + 1):
        # Vectorized step on int64-resident alive trajectories
        sub = np.where(alive & ~np.isin(np.arange(N), list(overflow_idx)
                                        if overflow_idx else np.array([],
                                                                      dtype=np.int64)))[0]
        if sub.size > 0:
            ns = n_arr[sub]
            m = 3 * ns + 1
            while True:
                even = (m & 1) == 0
                if not even.any():
                    break
                m[even] >>= 1
            # Detect entries that need bigint promotion (will overflow next step)
            big_mask = m > INT64_GUARD
            if big_mask.any():
                big_local = np.where(big_mask)[0]
                for li in big_local:
                    gi = sub[li]
                    overflow_idx.add(int(gi))
                    py_vals[int(gi)] = int(m[li])
                # zero out their int64 slot to avoid further use
                m[big_mask] = 0
            n_arr[sub] = m
            collapsed = (m == 1)
            if collapsed.any():
                alive[sub[collapsed]] = False

        # Python int step on overflow trajectories
        if overflow_idx:
            to_remove = []
            for gi in list(overflow_idx):
                if not alive[gi]:
                    to_remove.append(gi)
                    continue
                v = py_vals[gi]
                m = 3 * v + 1
                while not (m & 1):
                    m >>= 1
                if m <= INT64_GUARD:
                    # demote back to int64
                    n_arr[gi] = m
                    del py_vals[gi]
                    to_remove.append(gi)
                else:
                    py_vals[gi] = m
                if m == 1:
                    alive[gi] = False
            for gi in to_remove:
                overflow_idx.discard(gi)

        n_overflow_now = len(overflow_idx)
        if n_overflow_now > n_overflow:
            n_overflow = n_overflow_now

        if depth == S15_DEPTH:
            in_S15 = alive.copy()

        if depth in depth_set:
            # Full ensemble
            alive_idx = np.where(alive)[0]
            # values: from n_arr or py_vals
            if overflow_idx:
                # split: int64 alive + overflow alive
                int64_alive = np.array([i for i in alive_idx
                                        if i not in overflow_idx],
                                       dtype=np.int64)
                py_alive = [i for i in alive_idx if i in overflow_idx]
                vals = n_arr[int64_alive]
                # process int64 alive
                if vals.size > 0:
                    H_full[depth] += histogram(vals, inv3)
                # process py alive: residues
                for gi in py_alive:
                    v = py_vals[gi]
                    r2 = v & (M2 - 1)
                    r3 = v % M3
                    if r2 % 2 == 1 and r3 % 3 != 0:
                        i = r2 >> 1
                        j = inv3[r3]
                        H_full[depth][i, j] += 1
                n_alive_full[depth] = alive_idx.size
            else:
                vals = n_arr[alive_idx]
                if vals.size > 0:
                    H_full[depth] += histogram(vals, inv3)
                n_alive_full[depth] = alive_idx.size

            # Controlled subset S_15
            if depth >= S15_DEPTH and in_S15 is not None:
                # alive AND was alive at k=15
                ctrl_mask = alive & in_S15
                ctrl_idx = np.where(ctrl_mask)[0]
                if overflow_idx:
                    int64_ctrl = np.array([i for i in ctrl_idx
                                           if i not in overflow_idx],
                                          dtype=np.int64)
                    py_ctrl = [i for i in ctrl_idx if i in overflow_idx]
                    if int64_ctrl.size > 0:
                        H_s15[depth] += histogram(n_arr[int64_ctrl], inv3)
                    for gi in py_ctrl:
                        v = py_vals[gi]
                        r2 = v & (M2 - 1)
                        r3 = v % M3
                        if r2 % 2 == 1 and r3 % 3 != 0:
                            i = r2 >> 1
                            j = inv3[r3]
                            H_s15[depth][i, j] += 1
                else:
                    if ctrl_idx.size > 0:
                        H_s15[depth] += histogram(n_arr[ctrl_idx], inv3)
                n_alive_s15[depth] = ctrl_idx.size

            print(f"  [t={time.time()-t0:6.1f}s] depth={depth:>2}  "
                  f"alive_full={n_alive_full[depth]:>9,}  "
                  f"alive_s15={n_alive_s15.get(depth, '-'):>9}  "
                  f"overflow={n_overflow_now}", flush=True)

    print(f"[t={time.time()-t0:6.1f}s] computing chi2 and growth fit",
          flush=True)

    rows_data = []
    for k in DEPTHS_SAMPLE:
        chi2_f, cpdf_f, z_f = chi2_df_from_H(H_full[k], n_alive_full[k])
        if k >= S15_DEPTH:
            chi2_s, cpdf_s, z_s = chi2_df_from_H(H_s15[k], n_alive_s15[k])
            n_s = n_alive_s15[k]
        else:
            chi2_s = float("nan"); cpdf_s = float("nan"); z_s = float("nan")
            n_s = 0
        rows_data.append({
            "k": k,
            "n_alive_full": n_alive_full[k],
            "frac_alive_full": n_alive_full[k] / N,
            "chi2_full": chi2_f,
            "cpdf_full": cpdf_f,
            "z_full": z_f,
            "n_alive_s15": n_s,
            "chi2_s15": chi2_s,
            "cpdf_s15": cpdf_s,
            "z_s15": z_s,
        })
        print(f"  k={k:>2}  alive_full={n_alive_full[k]:>9,} "
              f"({100*n_alive_full[k]/N:5.2f}%)  "
              f"cpdf_full={cpdf_f:7.3f}  z_full={z_f:+8.2f}  "
              f"cpdf_s15={cpdf_s:7.3f}  z_s15={z_s:+8.2f}",
              flush=True)

    # Geometric fit on full ensemble
    fit_ks = np.array(DEPTHS_FIT, dtype=float)
    fit_y = np.array([r["cpdf_full"] for r in rows_data
                      if r["k"] in DEPTHS_FIT], dtype=float)
    fit_logy = np.log(fit_y)
    slope, intercept = np.polyfit(fit_ks, fit_logy, 1)
    r_growth = float(np.exp(slope))
    chi2_0 = float(np.exp(intercept))
    print(f"[fit] chi2/df ~ {chi2_0:.4e} * {r_growth:.4f}^k    "
          f"(log-linear LS on k={DEPTHS_FIT})", flush=True)

    candidates = {
        "4/3": 4 / 3,
        "log_2(3)": np.log2(3),
        "1/log_2(3)": 1 / np.log2(3),
        "3/2": 1.5,
        "exp(1/3)": np.exp(1 / 3),
    }
    closest = min(candidates.items(), key=lambda kv: abs(kv[1] - r_growth))
    print(f"[fit] closest candidate: {closest[0]} = {closest[1]:.4f} "
          f"(|diff|={abs(closest[1]-r_growth):.4f})", flush=True)

    # Predicted at k=40, 50
    pred_40 = chi2_0 * r_growth ** 40
    pred_50 = chi2_0 * r_growth ** 50
    print(f"[predict] k=40: chi2/df = {pred_40:.2f}  "
          f"actual = {rows_data[DEPTHS_SAMPLE.index(40)]['cpdf_full']:.2f}",
          flush=True)
    print(f"[predict] k=50: chi2/df = {pred_50:.2f}  "
          f"actual = {rows_data[DEPTHS_SAMPLE.index(50)]['cpdf_full']:.2f}",
          flush=True)

    # Top cells at k=20 and k=30
    top20 = top_cells(H_full[20], n_alive_full[20])
    top30 = top_cells(H_full[30], n_alive_full[30])

    # CSV
    with open(CSV_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["k", "n_alive_full", "frac_alive_full",
                    "chi2_full", "chi2_per_df_full", "z_full",
                    "n_alive_s15", "chi2_s15", "chi2_per_df_s15", "z_s15"])
        for r in rows_data:
            w.writerow([r["k"], r["n_alive_full"], f"{r['frac_alive_full']:.6f}",
                        f"{r['chi2_full']:.4f}", f"{r['cpdf_full']:.6f}",
                        f"{r['z_full']:.4f}",
                        r["n_alive_s15"],
                        f"{r['chi2_s15']:.4f}", f"{r['cpdf_s15']:.6f}",
                        f"{r['z_s15']:.4f}"])

    with open(TOP_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["k", "rank", "r_2", "r_3", "P", "Q", "P_minus_Q", "z_cell"])
        for rank, c in enumerate(top20, 1):
            w.writerow([20, rank, c["r2"], c["r3"],
                        f"{c['P']:.6e}", f"{c['Q']:.6e}",
                        f"{c['D']:+.6e}", f"{c['z_cell']:.4f}"])
        for rank, c in enumerate(top30, 1):
            w.writerow([30, rank, c["r2"], c["r3"],
                        f"{c['P']:.6e}", f"{c['Q']:.6e}",
                        f"{c['D']:+.6e}", f"{c['z_cell']:.4f}"])

    # Markdown
    lines = []
    lines.append("# Bohr Probe — Verification (survival, controlled subset, "
                 "growth fit, dominant cells)")
    lines.append("")
    lines.append(f"N = {N:,}; (a, b) = ({A}, {B}); d.f. = {DF}")
    lines.append("")

    lines.append("## A. Survival")
    lines.append("")
    lines.append("| k | N_alive | fraction |")
    lines.append("|---:|---:|---:|")
    for r in rows_data:
        lines.append(f"| {r['k']} | {r['n_alive_full']:,} | "
                     f"{r['frac_alive_full']:.4f} |")
    lines.append("")

    lines.append("## B. Full ensemble vs controlled subset S_15 "
                 "(alive-at-k=15)")
    lines.append("")
    lines.append("| k | N_full | chi²/df full | z full | N_S15 | "
                 "chi²/df S15 | z S15 |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows_data:
        ns = f"{r['n_alive_s15']:,}" if r["k"] >= S15_DEPTH else "—"
        cps = f"{r['cpdf_s15']:.3f}" if r["k"] >= S15_DEPTH else "—"
        zs = f"{r['z_s15']:+.2f}" if r["k"] >= S15_DEPTH else "—"
        lines.append(f"| {r['k']} | {r['n_alive_full']:,} | "
                     f"{r['cpdf_full']:.3f} | {r['z_full']:+.2f} | "
                     f"{ns} | {cps} | {zs} |")
    lines.append("")

    lines.append("## C. Geometric fit on k ∈ {18, 19, 20, 22, 25, 30}")
    lines.append("")
    lines.append(f"chi²/df ≈ {chi2_0:.4e} · r^k")
    lines.append("")
    lines.append(f"**r = {r_growth:.4f}** (log-linear LS)")
    lines.append("")
    lines.append("Candidates:")
    for name, val in candidates.items():
        lines.append(f"- {name} = {val:.4f}  (|diff| = "
                     f"{abs(val-r_growth):.4f})")
    lines.append(f"")
    lines.append(f"**Closest candidate:** {closest[0]} = {closest[1]:.4f}")
    lines.append("")

    lines.append("## D. Saturation test (k=40, 50)")
    lines.append("")
    lines.append("| k | predicted chi²/df | actual chi²/df | "
                 "actual / predicted |")
    lines.append("|---:|---:|---:|---:|")
    for k_ext in [40, 50]:
        pred = chi2_0 * r_growth ** k_ext
        actual = next(r for r in rows_data if r["k"] == k_ext)["cpdf_full"]
        ratio = actual / pred if pred > 0 else float("nan")
        lines.append(f"| {k_ext} | {pred:.2f} | {actual:.2f} | {ratio:.3f} |")
    lines.append("")

    lines.append("## E. Top 10 deviating cells: k=20 vs k=30")
    lines.append("")
    lines.append("### k = 20")
    lines.append("")
    lines.append("| rank | r_2 | r_3 | P | Q | P−Q | z_cell |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|")
    for rk, c in enumerate(top20, 1):
        lines.append(f"| {rk} | {c['r2']} | {c['r3']} | {c['P']:.4e} | "
                     f"{c['Q']:.4e} | {c['D']:+.3e} | {c['z_cell']:+.2f} |")
    lines.append("")
    lines.append("### k = 30")
    lines.append("")
    lines.append("| rank | r_2 | r_3 | P | Q | P−Q | z_cell |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|")
    for rk, c in enumerate(top30, 1):
        lines.append(f"| {rk} | {c['r2']} | {c['r3']} | {c['P']:.4e} | "
                     f"{c['Q']:.4e} | {c['D']:+.3e} | {c['z_cell']:+.2f} |")
    lines.append("")

    # Cell overlap
    set20 = {(c["r2"], c["r3"]) for c in top20}
    set30 = {(c["r2"], c["r3"]) for c in top30}
    overlap = set20 & set30
    lines.append(f"**Cell overlap (top 10):** {len(overlap)} / 10 cells "
                 f"shared: {sorted(overlap)}")
    lines.append("")

    MD_PATH.write_text("\n".join(lines), encoding="utf-8")
    print(f"[done] runtime = {time.time()-t0:.1f}s; outputs in {OUTDIR}",
          flush=True)


if __name__ == "__main__":
    main()
