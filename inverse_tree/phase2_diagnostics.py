"""Phase 2 — structural diagnostics on the d≤50 inverse Collatz tree."""
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl

OUT = Path(r"C:\Collatz\inverse_tree")
df = pl.read_parquet(OUT / "tree_d50.parquet")
DMAX = int(df["depth"].max())

# ===== 2.1 Branching by depth =====
br_by_depth = (
    df.group_by("depth")
    .agg(
        n_nodes=pl.len(),
        n_branching=pl.col("child_right").is_not_null().sum(),
    )
    .sort("depth")
    .with_columns(
        branching_ratio=(pl.col("n_branching") / pl.col("n_nodes")).cast(pl.Float64)
    )
)
br_by_depth.write_parquet(OUT / "phase2_branching.parquet")

# ===== 2.2 a★_6-stratified branching =====
astar_br = (
    df.group_by(["depth", "a_final_6"])
    .agg(
        n_nodes=pl.len(),
        n_branching=pl.col("child_right").is_not_null().sum(),
        mean_n=pl.col("n").mean().cast(pl.Float64),
    )
    .sort(["depth", "a_final_6"])
    .with_columns(
        branching_ratio=(pl.col("n_branching") / pl.col("n_nodes")).cast(pl.Float64)
    )
)
astar_br.write_parquet(OUT / "phase2_astar_branching.parquet")

# ===== 2.3 P(a★_6 | depth) =====
totals = df.group_by("depth").agg(total=pl.len())
astar_density = (
    df.group_by(["depth", "a_final_6"])
    .agg(c=pl.len())
    .sort(["depth", "a_final_6"])
    .join(totals, on="depth")
    .with_columns(p=(pl.col("c") / pl.col("total")).cast(pl.Float64))
)
astar_density.write_parquet(OUT / "phase2_astar_density.parquet")

# ===== 2.5 Residue-class density per k =====
for k in (6, 8, 10):
    rd = (
        df.group_by(["depth", f"r_{k}"])
        .agg(c=pl.len())
        .sort(["depth", f"r_{k}"])
        .join(totals, on="depth")
        .with_columns(p=(pl.col("c") / pl.col("total")).cast(pl.Float64))
    )
    rd.write_parquet(OUT / f"phase2_residue_density_k{k}.parquet")

# ===== 2.4 Self-similarity per a★_6 class =====
n_to_depth = dict(zip(df["n"].to_list(), df["depth"].to_list()))
n_set = set(n_to_depth)
children_map = defaultdict(list)
for n, cl, cr in zip(df["n"].to_list(), df["child_left"].to_list(), df["child_right"].to_list()):
    if cl is not None and cl in n_set:
        children_map[n].append(cl)
    if cr is not None and cr in n_set:
        children_map[n].append(cr)

smallest_per_class = (
    df.group_by("a_final_6").agg(min_n=pl.col("n").min()).sort("a_final_6")
)

def subtree_layers(root_n):
    layers = [[root_n]]
    seen = {root_n}
    while layers[-1]:
        nxt = []
        for nn in layers[-1]:
            for c in children_map.get(nn, []):
                if c not in seen:
                    seen.add(c)
                    nxt.append(c)
        layers.append(nxt)
    if not layers[-1]:
        layers.pop()
    return layers, seen

ss_rows = []
for row in smallest_per_class.iter_rows(named=True):
    av = row["a_final_6"]
    rn = row["min_n"]
    rd = n_to_depth[rn]
    sub_layers, sub_set = subtree_layers(rn)
    max_offset = len(sub_layers) - 1
    if max_offset >= 10:
        offsets = np.arange(max(0, max_offset - 20), max_offset + 1)
        ys = np.log(np.array([len(sub_layers[o]) for o in offsets], dtype=float))
        slope, intercept = np.polyfit(offsets, ys, 1)
    else:
        slope, intercept = float("nan"), float("nan")
    sub_df = df.filter(pl.col("n").is_in(list(sub_set)))
    n_br = sub_df.filter(pl.col("child_right").is_not_null()).height
    br_ratio = n_br / sub_df.height if sub_df.height else float("nan")
    ss_rows.append(
        {
            "a_final_6": av,
            "root_n": rn,
            "root_depth": rd,
            "max_offset_in_window": max_offset,
            "subtree_total": sum(len(L) for L in sub_layers),
            "slope_recent20": float(slope),
            "branching_ratio_subtree": float(br_ratio),
        }
    )
ss_df = pl.DataFrame(ss_rows)
ss_df.write_parquet(OUT / "phase2_selfsimilarity.parquet")

# ===== Aggregate stats for the report =====
br_asymptote = float(
    br_by_depth.filter(pl.col("depth") >= 30).select(pl.col("branching_ratio").mean()).item()
)

# a★ collapse test: spread of branching_ratio across classes at each large d
spread_rows = []
for d in range(30, DMAX + 1):
    sub = astar_br.filter((pl.col("depth") == d) & (pl.col("n_nodes") >= 10))
    if sub.height >= 2:
        ratios = sub["branching_ratio"].to_numpy()
        spread_rows.append(
            {
                "depth": d,
                "spread": float(ratios.max() - ratios.min()),
                "std": float(ratios.std()),
            }
        )
spread_df = pl.DataFrame(spread_rows)
mean_spread = float(spread_df["spread"].mean())
mean_std = float(spread_df["std"].mean())

# Equidistribution chi² at d=DMAX over observed residues
chi_results = {}
for k in (6, 8, 10):
    deep = df.filter(pl.col("depth") == DMAX)
    counts = deep.group_by(f"r_{k}").agg(c=pl.len())
    n_total = deep.height
    n_unique = counts.height
    obs = counts["c"].to_numpy().astype(float)
    expected = n_total / n_unique
    chi2 = float(((obs - expected) ** 2 / expected).sum())
    chi_results[k] = {
        "k": k,
        "M": 1 << k,
        "n_unique": n_unique,
        "n_total": n_total,
        "chi2": chi2,
        "chi2_per_dof": chi2 / max(1, n_unique - 1),
    }

# Self-similarity verdict: compare sub-tree slopes to full-tree slope.
# Restrict to well-sampled sub-trees (≥100 nodes) so small-sample classes
# don't dominate the verdict.
FULL_SLOPE = 0.2343
SS_MIN_NODES = 100
ss_well_sampled = ss_df.filter(pl.col("subtree_total") >= SS_MIN_NODES)
ss_slope_dev = ss_well_sampled["slope_recent20"].to_numpy() - FULL_SLOPE
ss_slope_dev_max = float(np.nanmax(np.abs(ss_slope_dev))) if ss_well_sampled.height else float("nan")
ss_n_well = ss_well_sampled.height
ss_n_total = ss_df.height
# Mass coverage of well-sampled classes
mass_total = ss_df["subtree_total"].sum()
mass_well = ss_well_sampled["subtree_total"].sum()
mass_pct = 100.0 * mass_well / mass_total if mass_total else 0.0

# ===== Report =====
lines = [
    "# Phase 2 Diagnostics — Inverse Collatz Tree (d ≤ 50)",
    "",
    "Source: `tree_d50.parquet` (379,600 nodes).",
    "",
    "## 2.1 Branching ratio by depth",
    "",
    "| depth | n_nodes | n_branching | branching_ratio |",
    "|---:|---:|---:|---:|",
]
for d in (5, 10, 15, 20, 25, 30, 35, 40, 45, 50):
    row = br_by_depth.filter(pl.col("depth") == d).row(0, named=True)
    lines.append(
        f"| {d} | {row['n_nodes']:,} | {row['n_branching']:,} | {row['branching_ratio']:.4f} |"
    )
lines += [
    "",
    f"Mean branching ratio over d=30..50: **{br_asymptote:.4f}**.",
    "",
    "(Heuristic: ~1/3 of nodes are ≡1 mod 3; of those, roughly 1/2 yield odd (n−1)/3, "
    "giving an upper-bound expectation ~1/6 ≈ 0.167. Observed equilibrates around the empirical value above.)",
    "",
    "## 2.2 a★_6 stratified branching — collapse test",
    "",
    f"Across d=30..50 (classes with ≥10 nodes): mean spread (max − min branching_ratio) = "
    f"**{mean_spread:.4f}**, mean std = **{mean_std:.4f}**.",
    "",
    "**Verdict:** ",
]
if mean_spread < 0.03:
    lines.append("Branching ratio is approximately a★-INDEPENDENT — lines collapse. "
                 "a★ adds no incremental signal for branching prediction beyond what's "
                 "implicit in the mod-3 structure of n.")
elif mean_spread < 0.10:
    lines.append("Branching ratio shows MODEST a★-stratification — partial separation but "
                 "much smaller than the asymptotic ratio itself.")
else:
    lines.append("Branching ratio is a★-STRATIFIED — clear separation across classes "
                 "indicates a★ encodes branching-relevant structure.")

lines += [
    "",
    f"Per-class branching ratios at d={DMAX} (≥10 nodes):",
    "",
    "| a★_6 | n_nodes | branching_ratio |",
    "|---:|---:|---:|",
]
for row in (
    astar_br.filter((pl.col("depth") == DMAX) & (pl.col("n_nodes") >= 10))
    .sort("a_final_6")
    .iter_rows(named=True)
):
    lines.append(f"| {row['a_final_6']} | {row['n_nodes']:,} | {row['branching_ratio']:.4f} |")

lines += [
    "",
    "## 2.3 a★_6 density at deepest layer",
    "",
    f"P(a★_6 | depth={DMAX}):",
    "",
    "| a★_6 | count | P |",
    "|---:|---:|---:|",
]
for row in astar_density.filter(pl.col("depth") == DMAX).sort("a_final_6").iter_rows(named=True):
    lines.append(f"| {row['a_final_6']} | {row['c']:,} | {row['p']:.4f} |")
lines += [
    "",
    "Full table (all d): `phase2_astar_density.parquet`.",
    "",
    "## 2.4 Self-similarity per a★_6 class",
    "",
    f"Sub-tree rooted at smallest n in each a★_6 class. Slope fit on last 20 offset levels. "
    f"Full-tree slope (d=10..50) = **0.2343** for comparison.",
    "",
    "| a★_6 | root_n | root_depth | sub_total | slope (last 20) | Δ vs full | branching_ratio_sub |",
    "|---:|---:|---:|---:|---:|---:|---:|",
]
for row in ss_df.iter_rows(named=True):
    delta = row["slope_recent20"] - FULL_SLOPE
    lines.append(
        f"| {row['a_final_6']} | {row['root_n']:,} | {row['root_depth']} | "
        f"{row['subtree_total']:,} | {row['slope_recent20']:.4f} | {delta:+.4f} | "
        f"{row['branching_ratio_subtree']:.4f} |"
    )
lines += [
    "",
    f"**Verdict** (restricted to well-sampled sub-trees, ≥{SS_MIN_NODES} nodes — "
    f"covers {ss_n_well}/{ss_n_total} classes, **{mass_pct:.2f}%** of total mass):",
    "",
    f"max |slope deviation vs full tree| = **{ss_slope_dev_max:.4f}**.",
]
if ss_slope_dev_max < 0.02:
    lines.append("Sub-trees match the full tree's growth rate within tight tolerance — "
                 "the tree is approximately self-similar with respect to a★_6 stratification "
                 "across well-sampled classes.")
elif ss_slope_dev_max < 0.05:
    lines.append("Sub-trees match the full tree's growth rate moderately well; small "
                 "class-dependent deviations exist.")
else:
    lines.append("Sub-tree growth rates deviate substantially from the full tree — "
                 "self-similarity claim is NOT supported even on well-sampled classes.")
lines.append("")
lines.append(f"Small-sample classes excluded from the verdict: "
             f"{ss_n_total - ss_n_well} ({100 - mass_pct:.2f}% of mass) — "
             "their sub-trees are too shallow inside d≤50 to fit a meaningful slope.")

lines += [
    "",
    f"## 2.5 Residue-class density at d={DMAX} — equidistribution",
    "",
    "Chi² distance from uniform across observed residues (residues that don't appear are excluded from dof).",
    "",
    "| k | M=2^k | unique residues observed | n_total | chi² | chi²/dof |",
    "|---:|---:|---:|---:|---:|---:|",
]
for k in (6, 8, 10):
    r = chi_results[k]
    lines.append(
        f"| {k} | {r['M']} | {r['n_unique']} / {r['M']} | {r['n_total']:,} | "
        f"{r['chi2']:.1f} | {r['chi2_per_dof']:.2f} |"
    )

lines += [
    "",
    "Per-k full density tables: `phase2_residue_density_k{6,8,10}.parquet`.",
    "",
    "Note: chi²/dof ≈ 1 is consistent with uniform; >> 1 indicates residue-dependent over- or under-representation.",
    "",
    "## Outputs",
    "",
    "- `phase2_branching.parquet` — depth, n_nodes, n_branching, branching_ratio",
    "- `phase2_astar_branching.parquet` — (depth × a_final_6) cross-tab",
    "- `phase2_astar_density.parquet` — long-format P(a★_6 | depth)",
    "- `phase2_residue_density_k{6,8,10}.parquet` — long-format P(r_k | depth)",
    "- `phase2_selfsimilarity.parquet` — sub-tree growth-rate per class",
]
(OUT / "phase2_findings.md").write_text("\n".join(lines), encoding="utf-8")

print(f"[ok] mean_spread={mean_spread:.4f} mean_std={mean_std:.4f}")
print(f"[ok] mean branching d>=30: {br_asymptote:.4f}")
print(f"[ok] self-sim max abs(slope_dev) (well-sampled, n>={SS_MIN_NODES}): {ss_slope_dev_max:.4f}")
print(f"[ok] well-sampled classes: {ss_n_well}/{ss_n_total} = {mass_pct:.2f}% of mass")
for k, r in chi_results.items():
    print(f"[ok] k={k} chi2={r['chi2']:.1f} chi2/dof={r['chi2_per_dof']:.2f}")
print(f"[wrote] {OUT/'phase2_findings.md'}")
