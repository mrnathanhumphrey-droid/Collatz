"""Generate publication-quality plot for R79b empirical δ vs theoretical curves."""
import sys
import os
import math
import csv
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

# Use non-interactive backend
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


CSV_PATH = r"C:\Collatz\r79b_S_partial_data.csv"
OUT_PNG = r"C:\Collatz\r79b_empirical_delta.png"


def read_data(csv_path: str):
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "r": int(row["r"]),
                "q": int(row["q"]),
                "N": int(row["N"]),
                "K_max": float(row["K_max_abs"]),
                "K_c1m0": float(row["K_c1m0_abs"]),
                "rho": float(row["rho"]),
                "delta_emp": float(row["delta_emp"]),
                "baseline": float(row["random_baseline_abs"]),
            })
    rows.sort(key=lambda r: r["r"])
    return rows


def main():
    rows = read_data(CSV_PATH)
    if len(rows) < 3:
        print(f"Not enough data points: {len(rows)}")
        return

    rs = np.array([row["r"] for row in rows])
    qs = np.array([row["q"] for row in rows], dtype=np.float64)
    Ns = np.array([row["N"] for row in rows], dtype=np.float64)
    Ks_max = np.array([row["K_max"] for row in rows])
    Ks_c1 = np.array([row["K_c1m0"] for row in rows])
    bases = np.array([row["baseline"] for row in rows])

    # Two views: (1) δ_brief vs r, (2) log-log plot of |K| vs N
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # Panel A: empirical δ_brief vs r, with theoretical lines
    ax = axes[0, 0]
    delta_emp = -np.array([row["rho"] for row in rows]) / 2
    ax.plot(rs, delta_emp, "o-", label="Empirical δ (max over c, m sample)", color="#1f77b4")
    # Same for c=1, m=0 only:
    delta_c1 = (np.log(qs) - 2 * np.log(Ks_c1)) / (2 * np.log(qs))
    ax.plot(rs, delta_c1, "s--", label="Empirical δ (c=1, m=0 only)", color="#ff7f0e", alpha=0.7)
    # Theoretical predictions:
    ax.axhline(0, color="gray", linestyle=":", label="Trivial (Plancherel) δ = 0")
    ax.axhline(0.0855, color="green", linestyle="--", label="Sub-Weyl (Milićević ABA³B) δ ≈ 0.0855")
    ax.axhline(1/6, color="purple", linestyle="--", label="Weyl δ = 1/6 ≈ 0.167")
    ax.axhline(0.5, color="red", linestyle=":", label="Square-root cancellation δ = 1/2")
    ax.set_xlabel("r")
    ax.set_ylabel("δ (saving against q^{1/2})")
    ax.set_title("Empirical δ vs r (brief's parameterization)")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 0.55)

    # Panel B: log-log fit |K| vs N
    ax = axes[0, 1]
    ax.loglog(Ns, Ks_max, "o-", label="|K_max|", color="#1f77b4")
    ax.loglog(Ns, Ks_c1, "s--", label="|K_{c=1,m=0}|", color="#ff7f0e", alpha=0.7)
    # Square-root reference
    sqrt_N = np.sqrt(Ns)
    ax.loglog(Ns, sqrt_N, ":", color="red", label="√N (square-root)")
    ax.loglog(Ns, 2 * sqrt_N, "--", color="red", alpha=0.5, label="2√N")
    # Trivial
    ax.loglog(Ns, Ns, ":", color="gray", label="N (trivial)")
    # Random baseline
    ax.loglog(Ns, bases, "v", color="#2ca02c", alpha=0.5, label="Random baseline")
    ax.set_xlabel("N = 3^{r-1}")
    ax.set_ylabel("|K|")
    ax.set_title("|K_partial(r)| vs N (log-log)")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, which="both", alpha=0.3)

    # Panel C: |K|/√N ratio
    ax = axes[1, 0]
    ratios_max = Ks_max / sqrt_N
    ratios_c1 = Ks_c1 / sqrt_N
    ratios_base = bases / sqrt_N
    ax.plot(rs, ratios_max, "o-", label="|K_max| / √N", color="#1f77b4")
    ax.plot(rs, ratios_c1, "s--", label="|K_{c=1}| / √N", color="#ff7f0e", alpha=0.7)
    ax.plot(rs, ratios_base, "v:", label="random / √N", color="#2ca02c", alpha=0.5)
    ax.axhline(1.0, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("r")
    ax.set_ylabel("|K| / √N")
    ax.set_title("Pre-factor: |K|/√N (constant ⟺ exact rate-1/2)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel D: r mod 3 dependence
    ax = axes[1, 1]
    for k in [0, 1, 2]:
        mask = (rs % 3 == k)
        ax.plot(rs[mask], (Ks_max / sqrt_N)[mask], "o-", label=f"r ≡ {k} (mod 3)", alpha=0.8)
    ax.set_xlabel("r")
    ax.set_ylabel("|K|/√N")
    ax.set_title("r mod 3 oscillation check (Subtask C2)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.suptitle(f"R79b: Empirical |S_partial(r)| / Kalafatelis-eq-190 sum  (n = {len(rows)} points)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=140, bbox_inches="tight")
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()
