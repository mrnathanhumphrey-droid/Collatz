"""Phase 3 — visualizations: tree to d<=15, a-star density stacked-area, branching by stratum."""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = Path(r"C:\Collatz\inverse_tree")
df = pl.read_parquet(OUT / "tree_d50.parquet")

# Color palette for a-star_6 — fixed mapping across all 3 figures
ALL_ASTAR = sorted(df["a_final_6"].unique().to_list())
cmap = plt.colormaps["tab10"]
A_COLOR = {a: cmap(i) for i, a in enumerate(ALL_ASTAR)}

# ===== 3.1 Tree visualization, depth <= 15 =====
DTREE = 15
sub = df.filter(pl.col("depth") <= DTREE)
n_at = set(sub["n"].to_list())

children_map = defaultdict(list)
for n, cl, cr, d in zip(
    sub["n"].to_list(),
    sub["child_left"].to_list(),
    sub["child_right"].to_list(),
    sub["depth"].to_list(),
):
    if d == DTREE:
        continue
    if cl is not None and cl in n_at:
        children_map[n].append(cl)
    if cr is not None and cr in n_at:
        children_map[n].append(cr)
for parent in children_map:
    children_map[parent] = sorted(children_map[parent])

# Recursive layout: leaves get sequential y; internal nodes get midpoint of children
y_pos = {}
counter = [0]

def assign_y(n):
    kids = children_map.get(n, [])
    if not kids:
        y_pos[n] = counter[0]
        counter[0] += 1
        return y_pos[n]
    ys = [assign_y(c) for c in kids]
    y_pos[n] = (min(ys) + max(ys)) / 2.0
    return y_pos[n]

assign_y(1)

n_to_a = dict(zip(sub["n"].to_list(), sub["a_final_6"].to_list()))
n_to_d = dict(zip(sub["n"].to_list(), sub["depth"].to_list()))

fig, ax = plt.subplots(figsize=(14, 9))
# Edges
for parent, kids in children_map.items():
    px, py = n_to_d[parent], y_pos[parent]
    for c in kids:
        cx, cy = n_to_d[c], y_pos[c]
        ax.plot([px, cx], [py, cy], "-", color="0.6", lw=0.6, alpha=0.7, zorder=1)

# Nodes
ns = sub["n"].to_list()
xs = [n_to_d[n] for n in ns]
ys = [y_pos[n] for n in ns]
colors = [A_COLOR[n_to_a[n]] for n in ns]
sizes = [30 + 6 * np.log2(max(n, 1)) for n in ns]
ax.scatter(xs, ys, c=colors, s=sizes, edgecolors="black", linewidths=0.4, zorder=2)

# Label small-n nodes
for n in ns:
    if n <= 64:
        ax.annotate(
            str(n), (n_to_d[n], y_pos[n]),
            xytext=(4, 4), textcoords="offset points",
            fontsize=7, color="black",
        )

legend_handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=A_COLOR[a],
           markersize=10, markeredgecolor="black", label=f"a* = {a}")
    for a in sorted(set(n_to_a.values()))
]
ax.legend(handles=legend_handles, loc="upper left", fontsize=9,
          title="prefix terminal a* at k=6")
ax.set_xlabel("depth (= sigma(n))")
ax.set_ylabel("")
ax.set_title(f"Inverse Collatz tree, depth <= {DTREE} ({sub.height} nodes), colored by a*_6")
ax.set_xticks(range(0, DTREE + 1))
ax.grid(True, axis="x", alpha=0.3)
ax.set_yticks([])
plt.tight_layout()
plt.savefig(OUT / "phase3_tree_d15.svg", bbox_inches="tight")
plt.savefig(OUT / "phase3_tree_d15.png", bbox_inches="tight", dpi=150)
plt.close()
print(f"[wrote] phase3_tree_d15.{{svg,png}} ({sub.height} nodes)")

# ===== 3.2 a-star density stacked-area =====
density = pl.read_parquet(OUT / "phase2_astar_density.parquet")
depths = sorted(density["depth"].unique().to_list())
a_vals = sorted(density["a_final_6"].unique().to_list())
mat = np.zeros((len(a_vals), len(depths)))
for row in density.iter_rows(named=True):
    ai = a_vals.index(row["a_final_6"])
    di = depths.index(row["depth"])
    mat[ai, di] = row["p"]

fig, ax = plt.subplots(figsize=(12, 6))
colors_stack = [A_COLOR[a] for a in a_vals]
ax.stackplot(depths, mat, labels=[f"a* = {a}" for a in a_vals],
             colors=colors_stack, alpha=0.85, edgecolor="white", linewidth=0.4)
ax.set_xlabel("depth")
ax.set_ylabel("P(a*_6 | depth)")
ax.set_title("a*_6 density over depth (stacked-area)")
ax.set_xlim(0, 50)
ax.set_ylim(0, 1)
ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT / "phase3_astar_density.svg", bbox_inches="tight")
plt.savefig(OUT / "phase3_astar_density.png", bbox_inches="tight", dpi=150)
plt.close()
print(f"[wrote] phase3_astar_density.{{svg,png}}")

# ===== 3.3 Branching ratio per a-star stratum =====
br = pl.read_parquet(OUT / "phase2_astar_branching.parquet")
br_full = pl.read_parquet(OUT / "phase2_branching.parquet").filter(pl.col("depth") >= 5).sort("depth")

fig, ax = plt.subplots(figsize=(12, 6))
for a in a_vals:
    rows = (
        br.filter((pl.col("a_final_6") == a) & (pl.col("n_nodes") >= 10))
        .sort("depth")
    )
    if rows.height < 3:
        continue
    ax.plot(
        rows["depth"].to_numpy(),
        rows["branching_ratio"].to_numpy(),
        "-o", markersize=4, label=f"a* = {a}", color=A_COLOR[a],
    )
ax.plot(
    br_full["depth"].to_numpy(),
    br_full["branching_ratio"].to_numpy(),
    "--", color="black", lw=1.5, alpha=0.7, label="all nodes",
)
ax.set_xlabel("depth")
ax.set_ylabel("branching ratio")
ax.set_title("Branching ratio vs depth, per a*_6 stratum (cells with >=10 nodes)")
ax.set_xlim(0, 50)
ax.set_ylim(0, 0.5)
ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT / "phase3_branching_by_stratum.svg", bbox_inches="tight")
plt.savefig(OUT / "phase3_branching_by_stratum.png", bbox_inches="tight", dpi=150)
plt.close()
print(f"[wrote] phase3_branching_by_stratum.{{svg,png}}")

print("[done]")
