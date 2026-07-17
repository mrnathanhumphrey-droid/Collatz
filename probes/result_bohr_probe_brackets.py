"""
result_bohr_probe_brackets.py — decompose cell mass at k=50 by integer-value
bracket to separate small-attractor funnel from residue-level joint structure.

Brackets at depth k=50:
  A: 1 ≤ v ≤ 100
  B: 100 < v ≤ 10^4
  C: 10^4 < v ≤ 10^6
  D: 10^6 < v ≤ 10^9
  E: v > 10^9

For each diagonal cell (r, r) with r ∈ {5, 7, 11, 13, 17, 19, 23, 25, 29, 31}
and 5 auto-selected non-diagonal control cells (|z| < 0.5, Q·N > 500),
report: per-bracket count and fraction.

Decision logic (in writeup): diagonal cells dominantly bracket A → small-
attractor funnel; diagonal distribution similar to control with elevated
totals → residue-level joint structure; mixed → both.
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
TARGET_K = 50
A_, B_ = 5, 4
M2, M3 = 32, 81
ROWS, COLS = 16, 54
DF = (ROWS - 1) * (COLS - 1)
LOW, HIGH = 500_000_000, 500_000_000_000
INT64_GUARD = 1 << 60
DIAGONAL_RS = [5, 7, 11, 13, 17, 19, 23, 25, 29, 31]
N_CONTROL = 5
BRACKET_NAMES = ["A:[1,100]", "B:(100,1e4]", "C:(1e4,1e6]",
                 "D:(1e6,1e9]", "E:>1e9"]
BRACKET_BOUNDS = [100, 10_000, 1_000_000, 1_000_000_000]

OUTDIR = Path(r"C:\Collatz")
CSV_PATH = OUTDIR / "result_bohr_probe_brackets.csv"
MD_PATH = OUTDIR / "result_bohr_probe_brackets.md"


def coprime_residues_mod_3b(b):
    return np.array([r for r in range(3 ** b) if r % 3 != 0], dtype=np.int64)


def odd_residues_mod_2a(a):
    return np.array([r for r in range(1 << a) if r % 2 == 1], dtype=np.int64)


def build_inv3():
    inv = -np.ones(M3, dtype=np.int64)
    cnt = 0
    for r in range(M3):
        if r % 3 != 0:
            inv[r] = cnt
            cnt += 1
    return inv


def histogram(vals, inv3):
    r2 = (vals & (M2 - 1)).astype(np.int64)
    row = (r2 >> 1)
    r3 = (vals % M3).astype(np.int64)
    col = inv3[r3]
    flat = row * COLS + col
    return np.bincount(flat, minlength=ROWS * COLS).reshape(ROWS, COLS)


def bracket_of(v):
    if v <= 100:
        return 0
    if v <= 10_000:
        return 1
    if v <= 1_000_000:
        return 2
    if v <= 1_000_000_000:
        return 3
    return 4


def main():
    t0 = time.time()
    print(f"[t={time.time()-t0:6.1f}s] bracket probe; N={N:,}; "
          f"k={TARGET_K}", flush=True)

    inv3 = build_inv3()
    odd_r2 = odd_residues_mod_2a(A_)
    cop_r3 = coprime_residues_mod_3b(B_)

    rng = np.random.default_rng(seed=SEED)
    x = rng.integers(low=LOW, high=HIGH, size=N, dtype=np.int64)
    n_arr = (2 * x + 1).astype(np.int64)
    alive = np.ones(N, dtype=bool)
    overflow = np.zeros(N, dtype=bool)
    py_vals = {}
    n_promoted = 0

    for depth in range(1, TARGET_K + 1):
        int64_mask = alive & ~overflow
        sub = np.where(int64_mask)[0]
        if sub.size > 0:
            ns = n_arr[sub]
            m = 3 * ns + 1
            while True:
                even = (m & 1) == 0
                if not even.any():
                    break
                m[even] >>= 1
            big_mask = m > INT64_GUARD
            if big_mask.any():
                for li in np.where(big_mask)[0]:
                    gi = sub[li]
                    overflow[gi] = True
                    py_vals[gi] = int(m[li])
                    n_promoted += 1
            n_arr[sub] = m
            collapsed = (m == 1)
            if collapsed.any():
                alive[sub[collapsed]] = False

        py_idx = [gi for gi in list(py_vals) if alive[gi]]
        for gi in py_idx:
            v = py_vals[gi]
            m = 3 * v + 1
            while not (m & 1):
                m >>= 1
            if m == 1:
                alive[gi] = False
                del py_vals[gi]
                overflow[gi] = False
                n_arr[gi] = 1
            elif m <= INT64_GUARD:
                n_arr[gi] = m
                del py_vals[gi]
                overflow[gi] = False
            else:
                py_vals[gi] = m

    n_alive_50 = int(alive.sum())
    print(f"[t={time.time()-t0:6.1f}s] depth=50 alive={n_alive_50:,} "
          f"({100*n_alive_50/N:.2f}%)  overflow={n_promoted}", flush=True)

    # Build H[50]
    H = np.zeros((ROWS, COLS), dtype=np.int64)
    int64_alive_idx = np.where(alive & ~overflow)[0]
    if int64_alive_idx.size > 0:
        H += histogram(n_arr[int64_alive_idx], inv3)
    for gi in [g for g in py_vals if alive[g]]:
        v = py_vals[gi]
        r2 = v & (M2 - 1)
        r3 = v % M3
        if r2 % 2 == 1 and r3 % 3 != 0:
            H[r2 >> 1, inv3[r3]] += 1

    # Chi² and z grid
    P = H.astype(np.float64) / n_alive_50
    Mr = P.sum(axis=1, keepdims=True)
    Mc = P.sum(axis=0, keepdims=True)
    Q = Mr * Mc
    D = P - Q
    with np.errstate(divide="ignore", invalid="ignore"):
        z_grid = np.where(Q > 0,
                          D * np.sqrt(n_alive_50 / np.maximum(Q, 1e-30)),
                          0.0)

    # Pick control cells: |z| < 0.5, Q*N > 500, non-diagonal
    Q_threshold = 500 / n_alive_50
    cands = []
    for i in range(ROWS):
        for j in range(COLS):
            r2 = int(odd_r2[i])
            r3 = int(cop_r3[j])
            if r2 == r3:
                continue
            if Q[i, j] < Q_threshold:
                continue
            if abs(z_grid[i, j]) >= 0.5:
                continue
            cands.append((abs(z_grid[i, j]), r2, r3, i, j))
    cands.sort()
    control_cells = cands[:N_CONTROL]
    print(f"[t={time.time()-t0:6.1f}s] picked {len(control_cells)} control "
          f"cells (|z|<0.5, Q·N>500)", flush=True)
    for _, r2, r3, _, _ in control_cells:
        print(f"  control: ({r2}, {r3})  z={z_grid[r2>>1, inv3[r3]]:+.3f}",
              flush=True)

    # Build target cell set
    target_cells = set()
    for r in DIAGONAL_RS:
        target_cells.add((r, r))
    for _, r2, r3, _, _ in control_cells:
        target_cells.add((r2, r3))

    # Bracket decomposition
    counts = {cell: [0] * 5 for cell in target_cells}
    print(f"[t={time.time()-t0:6.1f}s] bracket decomposition", flush=True)
    # int64 alive
    if int64_alive_idx.size > 0:
        vals = n_arr[int64_alive_idx]
        for v in vals:
            v_int = int(v)
            r2 = v_int & (M2 - 1)
            r3 = v_int % M3
            if (r2, r3) in target_cells:
                counts[(r2, r3)][bracket_of(v_int)] += 1
    # python int alive
    for gi in [g for g in py_vals if alive[g]]:
        v = py_vals[gi]
        r2 = v & (M2 - 1)
        r3 = v % M3
        if (r2, r3) in target_cells:
            counts[(r2, r3)][bracket_of(v)] += 1

    # Build summary rows
    rows = []
    for r in DIAGONAL_RS:
        cell = (r, r)
        cnts = counts[cell]
        total = sum(cnts)
        # cell-level chi² metrics
        i = r >> 1 if r < M2 else None
        j = inv3[r] if r < M3 else None
        z_c = float(z_grid[r >> 1, inv3[r]]) if (r % 2 == 1 and r % 3 != 0) \
            else float("nan")
        rows.append({
            "kind": "diag",
            "r2": r, "r3": r,
            "z_cell": z_c,
            "total": total,
            "counts": cnts,
            "fracs": [c / total if total > 0 else 0.0 for c in cnts],
        })

    for _, r2, r3, i, j in control_cells:
        cell = (r2, r3)
        cnts = counts[cell]
        total = sum(cnts)
        z_c = float(z_grid[i, j])
        rows.append({
            "kind": "control",
            "r2": r2, "r3": r3,
            "z_cell": z_c,
            "total": total,
            "counts": cnts,
            "fracs": [c / total if total > 0 else 0.0 for c in cnts],
        })

    # Write CSV
    with open(CSV_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["kind", "r2", "r3", "z_cell", "total",
                    "n_A", "n_B", "n_C", "n_D", "n_E",
                    "f_A", "f_B", "f_C", "f_D", "f_E"])
        for r in rows:
            w.writerow([r["kind"], r["r2"], r["r3"],
                        f"{r['z_cell']:.4f}", r["total"],
                        *r["counts"],
                        *[f"{f:.6f}" for f in r["fracs"]]])

    # Markdown writeup
    L = []
    L.append("# Bohr Probe — Bracket Decomposition at k=50")
    L.append("")
    L.append(f"N = {N:,}; (a, b) = ({A_}, {B_}); k = {TARGET_K}; "
             f"d.f. = {DF}; N_alive = {n_alive_50:,} "
             f"({100*n_alive_50/N:.2f}%)")
    L.append("")
    L.append("## Brackets")
    L.append("")
    L.append("| label | range |")
    L.append("|---|---|")
    L.append("| A | 1 ≤ v ≤ 100 |")
    L.append("| B | 100 < v ≤ 10⁴ |")
    L.append("| C | 10⁴ < v ≤ 10⁶ |")
    L.append("| D | 10⁶ < v ≤ 10⁹ |")
    L.append("| E | v > 10⁹ |")
    L.append("")

    L.append("## Diagonal cells (r, r) for r ∈ {5, 7, 11, 13, 17, 19, 23, "
             "25, 29, 31}")
    L.append("")
    L.append("### Counts")
    L.append("")
    L.append("| r | total | A | B | C | D | E | z_cell |")
    L.append("|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        if row["kind"] != "diag":
            continue
        L.append(f"| {row['r2']} | {row['total']} | "
                 f"{row['counts'][0]} | {row['counts'][1]} | "
                 f"{row['counts'][2]} | {row['counts'][3]} | "
                 f"{row['counts'][4]} | {row['z_cell']:+.2f} |")
    L.append("")

    L.append("### Fractions")
    L.append("")
    L.append("| r | total | f_A | f_B | f_C | f_D | f_E |")
    L.append("|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        if row["kind"] != "diag":
            continue
        L.append(f"| {row['r2']} | {row['total']} | "
                 f"{row['fracs'][0]:.3f} | {row['fracs'][1]:.3f} | "
                 f"{row['fracs'][2]:.3f} | {row['fracs'][3]:.3f} | "
                 f"{row['fracs'][4]:.3f} |")
    L.append("")

    L.append(f"## Control cells (non-diagonal, |z| < 0.5, Q·N > 500) — "
             f"{N_CONTROL} cells")
    L.append("")
    L.append("### Counts")
    L.append("")
    L.append("| (r2, r3) | total | A | B | C | D | E | z_cell |")
    L.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        if row["kind"] != "control":
            continue
        L.append(f"| ({row['r2']}, {row['r3']}) | {row['total']} | "
                 f"{row['counts'][0]} | {row['counts'][1]} | "
                 f"{row['counts'][2]} | {row['counts'][3]} | "
                 f"{row['counts'][4]} | {row['z_cell']:+.2f} |")
    L.append("")

    L.append("### Fractions")
    L.append("")
    L.append("| (r2, r3) | total | f_A | f_B | f_C | f_D | f_E |")
    L.append("|---|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        if row["kind"] != "control":
            continue
        L.append(f"| ({row['r2']}, {row['r3']}) | {row['total']} | "
                 f"{row['fracs'][0]:.3f} | {row['fracs'][1]:.3f} | "
                 f"{row['fracs'][2]:.3f} | {row['fracs'][3]:.3f} | "
                 f"{row['fracs'][4]:.3f} |")
    L.append("")

    # Summary aggregate
    diag_rows = [r for r in rows if r["kind"] == "diag"]
    ctrl_rows = [r for r in rows if r["kind"] == "control"]

    def aggregate_fracs(rs):
        total_count = sum(r["total"] for r in rs)
        if total_count == 0:
            return [0.0] * 5
        agg = [0] * 5
        for r in rs:
            for i in range(5):
                agg[i] += r["counts"][i]
        return [a / total_count for a in agg]

    diag_agg = aggregate_fracs(diag_rows)
    ctrl_agg = aggregate_fracs(ctrl_rows)

    L.append("## Aggregate comparison")
    L.append("")
    L.append("Pooled bracket fractions across all cells in each group:")
    L.append("")
    L.append("| group | total | f_A | f_B | f_C | f_D | f_E |")
    L.append("|---|---:|---:|---:|---:|---:|---:|")
    L.append(f"| diagonal (10 cells) | {sum(r['total'] for r in diag_rows):,} "
             f"| {diag_agg[0]:.3f} | {diag_agg[1]:.3f} | {diag_agg[2]:.3f} | "
             f"{diag_agg[3]:.3f} | {diag_agg[4]:.3f} |")
    L.append(f"| control (5 cells) | {sum(r['total'] for r in ctrl_rows):,} "
             f"| {ctrl_agg[0]:.3f} | {ctrl_agg[1]:.3f} | {ctrl_agg[2]:.3f} | "
             f"{ctrl_agg[3]:.3f} | {ctrl_agg[4]:.3f} |")
    L.append("")
    L.append("Excess factor (diagonal_frac / control_frac) per bracket:")
    L.append("")
    L.append("| bracket | diag/ctrl |")
    L.append("|---:|---:|")
    for i, name in enumerate(BRACKET_NAMES):
        ratio = (diag_agg[i] / ctrl_agg[i]) if ctrl_agg[i] > 0 \
            else float("inf")
        L.append(f"| {name} | {ratio:.3f} |")
    L.append("")

    L.append("## Decision read")
    L.append("")
    if diag_agg[0] > 0.5:
        verdict = ("**Small-attractor funnel dominates.** > 50% of diagonal-"
                   "cell mass is in bracket A (v ≤ 100). The 'joint Bohr' "
                   "finding is fundamentally about the Collatz descent funnel "
                   "structure; the residue diagonal is a consequence of "
                   "trajectories visiting small odd integers.")
    elif diag_agg[0] < 0.1 and abs(diag_agg[0] - ctrl_agg[0]) < 0.05:
        verdict = ("**Residue-level joint structure dominates.** Diagonal "
                   "cells distribute across brackets similarly to non-"
                   "diagonal cells but with elevated totals; the structure "
                   "is in residue configuration itself, independent of small-"
                   "attractor effects.")
    else:
        verdict = ("**Hybrid mechanism.** Diagonal cells show enriched "
                   "bracket A relative to controls but also significant mass "
                   "in larger brackets (D, E). Both small-attractor funnel "
                   "and residue-level joint structure contribute. The residue-"
                   "level structure can be characterized cleanly on the "
                   "bracket-D/E subset.")
    L.append(verdict)
    L.append("")
    L.append(f"Diagonal A-fraction: **{diag_agg[0]:.3f}**  |  "
             f"Control A-fraction: **{ctrl_agg[0]:.3f}**  |  "
             f"Diagonal/control ratio in A: **{diag_agg[0]/max(ctrl_agg[0], 1e-9):.1f}×**")
    L.append("")

    MD_PATH.write_text("\n".join(L), encoding="utf-8")
    print(f"\n[done] runtime = {time.time()-t0:.1f}s; outputs in {OUTDIR}",
          flush=True)


if __name__ == "__main__":
    main()
