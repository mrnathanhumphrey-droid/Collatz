"""Phase 6 — supplementary visuals.
1. branching by parity per a*_6 class
2. Q3 transition matrix heatmap + predicted vs empirical equilibrium bars
3. P((n-1)/3 odd | n%3==1) over depth
4. even-fraction per a*_6 class
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from pathlib import Path

import numpy as np
import polars as pl
import matplotlib.pyplot as plt

OUT = Path(r"C:\Collatz\inverse_tree")
parity = pl.read_parquet(OUT / "phase5_q1_parity.parquet").sort("a_final_6")
q2 = pl.read_parquet(OUT / "phase5_q2_branching_decomp.parquet").sort("depth")
q3 = pl.read_parquet(OUT / "phase5_q3_eq_eigenvector.parquet").sort("a_final_6")
M = np.loadtxt(OUT / "phase5_q3_transition_matrix.csv", delimiter=",", skiprows=1)

classes = parity["a_final_6"].to_list()
class_labels = [str(c) for c in classes]
cmap = plt.colormaps["tab10"]
A_COLOR = {a: cmap(i) for i, a in enumerate(classes)}

# ===== 1. Branching by parity per a*_6 =====
fig, ax = plt.subplots(figsize=(11, 6))
xs = np.arange(len(classes))
w = 0.38
br_even = [r if r == r else 0.0 for r in parity["br_even"].to_list()]
br_odd = [r if r == r else 0.0 for r in parity["br_odd"].to_list()]
n_even = parity["n_even"].to_list()
n_odd = parity["n_odd"].to_list()
b1 = ax.bar(xs - w/2, br_even, w, label="even n", color="#2b8cbe")
b2 = ax.bar(xs + w/2, br_odd, w, label="odd n", color="#e34a33")
ax.axhline(1/3, color="black", linestyle="--", lw=1, alpha=0.6, label="1/3 reference")
for i, (be, ne) in enumerate(zip(br_even, n_even)):
    if ne > 0:
        ax.annotate(f"n={ne:,}", (xs[i] - w/2, be + 0.01), ha="center", fontsize=7)
for i, (bo, no) in enumerate(zip(br_odd, n_odd)):
    if no > 0:
        ax.annotate(f"n={no:,}", (xs[i] + w/2, max(bo, 0) + 0.01), ha="center", fontsize=7)
ax.set_xticks(xs)
ax.set_xticklabels(class_labels)
ax.set_xlabel("a*_6 class")
ax.set_ylabel("branching ratio")
ax.set_title("Branching ratio by parity, within each a*_6 class (d=50)")
ax.set_ylim(0, 0.6)
ax.legend(loc="upper right")
ax.grid(True, axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(OUT / "phase6_branching_by_parity.svg")
plt.savefig(OUT / "phase6_branching_by_parity.png", dpi=150)
plt.close()
print("[wrote] phase6_branching_by_parity.{svg,png}")

# ===== 2. Q3 transition matrix heatmap + predicted vs empirical bars =====
fig, axes = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={"width_ratios": [1, 1.3]})

ax = axes[0]
im = ax.imshow(M, cmap="viridis", aspect="auto")
ax.set_xticks(range(len(classes)))
ax.set_yticks(range(len(classes)))
ax.set_xticklabels(class_labels)
ax.set_yticklabels(class_labels)
ax.set_xlabel("child a*_6")
ax.set_ylabel("parent a*_6")
ax.set_title("Transition matrix M[i,j] = E[# class-j children per class-i parent]")
for i in range(len(classes)):
    for j in range(len(classes)):
        if M[i, j] > 0.001:
            txt_color = "white" if M[i, j] < M.max() * 0.5 else "black"
            ax.text(j, i, f"{M[i,j]:.3f}", ha="center", va="center",
                    color=txt_color, fontsize=8)
plt.colorbar(im, ax=ax, fraction=0.046)

ax = axes[1]
emp = q3["empirical_eq"].to_numpy()
pred = q3["predicted_eq"].to_numpy()
xs2 = np.arange(len(classes))
w2 = 0.38
ax.bar(xs2 - w2/2, emp, w2, label="empirical (d=50)", color="#2b8cbe")
ax.bar(xs2 + w2/2, pred, w2, label="predicted (left eigvec)", color="#e34a33")
ax.set_xticks(xs2)
ax.set_xticklabels(class_labels)
ax.set_xlabel("a*_6 class")
ax.set_ylabel("equilibrium probability")
ax.set_title("Predicted vs empirical equilibrium  (max |diff| = 0.0003)")
ax.set_yscale("log")
ax.legend()
ax.grid(True, axis="y", alpha=0.3, which="both")
plt.tight_layout()
plt.savefig(OUT / "phase6_transition_and_equilibrium.svg")
plt.savefig(OUT / "phase6_transition_and_equilibrium.png", dpi=150)
plt.close()
print("[wrote] phase6_transition_and_equilibrium.{svg,png}")

# ===== 3. P((n-1)/3 odd | n%3==1) vs depth =====
fig, ax = plt.subplots(figsize=(11, 6))
ax.plot(q2["depth"].to_numpy(), q2["p_cond_eligible"].to_numpy(),
        "-o", markersize=4, label="P((n-1)/3 odd | n%3==1)", color="#2b8cbe")
ax.axhline(0.5, color="black", linestyle="--", lw=1, alpha=0.6, label="0.5 baseline (independence)")
ax.axhline(0.792, color="red", linestyle=":", lw=1, alpha=0.7, label="0.792 (asymptotic)")
ax.plot(q2["depth"].to_numpy(), q2["actual_branching"].to_numpy(),
        "-s", markersize=3, label="actual branching ratio", color="#e34a33", alpha=0.7)
ax.plot(q2["depth"].to_numpy(), q2["product"].to_numpy(),
        "--", lw=1, label="product (P(mod3) * P(cond))", color="green", alpha=0.7)
ax.set_xlabel("depth")
ax.set_ylabel("probability / ratio")
ax.set_title("Q2 — Conditional eligibility drives the 0.264 branching asymptote")
ax.legend(loc="lower right")
ax.set_ylim(0, 1)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT / "phase6_q2_eligibility_over_depth.svg")
plt.savefig(OUT / "phase6_q2_eligibility_over_depth.png", dpi=150)
plt.close()
print("[wrote] phase6_q2_eligibility_over_depth.{svg,png}")

# ===== 4. Even-fraction per a*_6 class =====
fig, ax = plt.subplots(figsize=(11, 6))
n_even_arr = np.array(parity["n_even"].to_list(), dtype=float)
n_odd_arr = np.array(parity["n_odd"].to_list(), dtype=float)
totals = n_even_arr + n_odd_arr
even_frac = np.where(totals > 0, n_even_arr / totals, np.nan)
colors_bar = [A_COLOR[c] for c in classes]
ax.bar(xs, even_frac, color=colors_bar, edgecolor="black", linewidth=0.5)
ax.axhline(0.5, color="black", linestyle="--", lw=1, alpha=0.6, label="0.5 reference")
for i, (ef, t) in enumerate(zip(even_frac, totals)):
    if t > 0:
        ax.annotate(f"{ef:.3f}\n(n={int(t):,})",
                    (xs[i], ef), ha="center", va="bottom" if ef < 0.95 else "top",
                    fontsize=8)
ax.set_xticks(xs)
ax.set_xticklabels(class_labels)
ax.set_xlabel("a*_6 class")
ax.set_ylabel("fraction of class with even n")
ax.set_title("Even-n fraction per a*_6 class (d=50)")
ax.set_ylim(0, 1.1)
ax.legend()
ax.grid(True, axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(OUT / "phase6_even_fraction_per_class.svg")
plt.savefig(OUT / "phase6_even_fraction_per_class.png", dpi=150)
plt.close()
print("[wrote] phase6_even_fraction_per_class.{svg,png}")

print("[done]")
