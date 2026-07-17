"""
inverse_tree_vs_davg.py — Result 52 candidate: test whether D_avg matches
inverse-tree visit-frequency distribution.

Hypothesis: surviving orbits at late t are "close to termination" (current m
has small remaining-time-to-1, equivalently: m is at small depth d from m=1
in the backward inverse Collatz tree). The residue distribution of such
m's is the inverse-tree visit-frequency.

Test multiple weighting schemes:
  (1) Uniform over odd nodes at depth d for various small d
  (2) Geometric weighting in d: w(d) ~ rho^d for various rho
  (3) Forward-orbit visit count: each starting m_0 contributes its full orbit
      residue counts -- this is the actual trajectory measure mechanism

Compare each candidate D_inv_tree(r) to empirical D_avg(r).
"""
import sys
from pathlib import Path
from fractions import Fraction
import numpy as np
import polars as pl

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(r"C:\Collatz")
log_lines = []
def log(s):
    print(s, flush=True)
    log_lines.append(s)


def chang_pi_32():
    """Chang's stationary mod 32 (essentially uniform 1/16)."""
    odd_residues = list(range(1, 64, 2))
    idx = {r: i for i, r in enumerate(odd_residues)}
    P = [[Fraction(0)] * 32 for _ in range(32)]
    for i, r in enumerate(odd_residues):
        counts = [0] * 32
        for k in range(128):
            n = r + 64 * k
            mm = 3 * n + 1
            while mm & 1 == 0:
                mm >>= 1
            counts[idx[mm % 64]] += 1
        for j in range(32):
            P[i][j] = Fraction(counts[j], 128)
    n = 32
    A = [[P[j][i] - (Fraction(1) if i == j else Fraction(0))
          for j in range(n)] for i in range(n)]
    A[n - 1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n - 1] = Fraction(1)
    for col in range(n):
        pivot = next(row for row in range(col, n) if A[row][col] != 0)
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            b[col], b[pivot] = b[pivot], b[col]
        piv = A[col][col]
        for j in range(col, n):
            A[col][j] /= piv
        b[col] /= piv
        for row in range(n):
            if row != col and A[row][col] != 0:
                factor = A[row][col]
                for j in range(col, n):
                    A[row][j] -= factor * A[col][j]
                b[row] -= factor * b[col]
    pi_32 = {r32: float(b[idx[r32]] + b[idx[r32 + 32]])
             for r32 in range(1, 32, 2)}
    return pi_32


def load_empirical_D_avg():
    D_avg = {}
    with open(OUT / "qsd_late_t_avg.csv") as f:
        f.readline()
        for line in f:
            parts = line.strip().split(",")
            D_avg[int(parts[0])] = float(parts[2])
    return D_avg


def main():
    log("=" * 78)
    log("Inverse-tree vs D_avg: candidate alternative kernel")
    log("=" * 78)

    pi_32 = chang_pi_32()
    D_avg = load_empirical_D_avg()

    # Load tree
    tree_path = OUT / "inverse_tree" / "tree_d50.parquet"
    log(f"\nLoading {tree_path}...")
    df = pl.read_parquet(tree_path)
    log(f"  {len(df):,} nodes, depth range {df['depth'].min()}-{df['depth'].max()}")

    # Add residue mod 32 column
    df = df.with_columns([
        (pl.col("n") % 32).alias("r32"),
        ((pl.col("n") % 2) == 1).alias("is_odd"),
    ])

    n_odd = df.filter(pl.col("is_odd")).height
    log(f"  Total odd nodes: {n_odd:,}")
    log(f"  Total even nodes: {len(df) - n_odd:,}")

    # ------------------------------------------------------------------
    # Diagnostic: odd-node count by depth
    # ------------------------------------------------------------------
    log("\nOdd-node count by depth d (in inverse tree from m=1):")
    by_d = (df.filter(pl.col("is_odd"))
              .group_by("depth")
              .len()
              .sort("depth"))
    for row in by_d.iter_rows(named=True):
        d, n = row["depth"], row["len"]
        if d <= 30 or d % 5 == 0:
            log(f"  d={d:>2}: {n:>7} odd nodes")

    # ------------------------------------------------------------------
    # Weighting scheme 1: residue distribution at depth d for various small d
    # ------------------------------------------------------------------
    log("\n" + "=" * 78)
    log("Scheme 1: residue distribution at fixed depth d (odd nodes only)")
    log("=" * 78)

    odd_r32 = list(range(1, 32, 2))

    def log_inline(s):
        # log() doesn't support end=; just print without newline
        sys.stdout.write(s)
        sys.stdout.flush()

    def D_at_depth(d_min, d_max):
        """D(r) = (count of odd nodes with r_32=r at depth in [d_min, d_max]) / pi(r)
        normalized so sum = 1 over odd residues, then divide by pi.
        Actually let's compute rho(r) = count(r) / total_count, then D = rho / pi.
        """
        sub = df.filter((pl.col("depth") >= d_min) &
                        (pl.col("depth") <= d_max) &
                        pl.col("is_odd"))
        n = sub.height
        if n == 0:
            return None, n
        counts = sub.group_by("r32").len().sort("r32")
        c_dict = {row["r32"]: row["len"] for row in counts.iter_rows(named=True)}
        rho = {r: c_dict.get(r, 0) / n for r in odd_r32}
        D = {r: rho[r] / pi_32[r] for r in odd_r32}
        return D, n

    test_ranges = [
        ("d=0..5", 0, 5),
        ("d=0..10", 0, 10),
        ("d=5..10", 5, 10),
        ("d=10..15", 10, 15),
        ("d=15..20", 15, 20),
        ("d=20..25", 20, 25),
        ("d=25..30", 25, 30),
        ("d=30..40", 30, 40),
        ("d=40..50", 40, 50),
    ]

    results = {}
    header = f"\n{'r':>3}  {'D_avg':>7}"
    for label, _, _ in test_ranges:
        header += f"  {label:>10}"
    log(header)

    # Compute all
    D_per_range = {}
    for label, d_min, d_max in test_ranges:
        D, n = D_at_depth(d_min, d_max)
        D_per_range[label] = (D, n)

    for r in odd_r32:
        line = f"  {r:>3}  {D_avg[r]:>7.4f}"
        for label, _, _ in test_ranges:
            D, n = D_per_range[label]
            if D is None:
                line += f"  {'(empty)':>10}"
            else:
                line += f"  {D[r]:>10.4f}"
        log(line)

    log(f"\nSample sizes (n_odd in range):")
    for label, _, _ in test_ranges:
        D, n = D_per_range[label]
        if D is not None:
            total_dev = sum(abs(D_avg[r] - D[r]) for r in odd_r32)
            log(f"  {label}: n={n:>6}, total |D - D_avg| = {total_dev:.4f}")

    # ------------------------------------------------------------------
    # Scheme 2: geometric weighting in d
    # ------------------------------------------------------------------
    log("\n" + "=" * 78)
    log("Scheme 2: geometric weighting w(d) = rho^d, sum over all odd nodes")
    log("=" * 78)

    odd_df = df.filter(pl.col("is_odd"))
    header = f"\n{'rho':>6}"
    for r in odd_r32:
        header += f"  r={r:>2}"
    header += "  total_dev"
    log(header)

    rho_values = [0.5, 0.7, 0.85, 0.94, 0.95, 0.97, 0.99]
    geom_results = {}
    for rho in rho_values:
        # weight each node by rho^d
        weights = rho ** odd_df["depth"].to_numpy()
        residues = odd_df["r32"].to_numpy()
        rho_dist = {}
        total_w = weights.sum()
        for r in odd_r32:
            mask = residues == r
            rho_dist[r] = weights[mask].sum() / total_w
        D = {r: rho_dist[r] / pi_32[r] for r in odd_r32}
        geom_results[rho] = D
        total_dev = sum(abs(D_avg[r] - D[r]) for r in odd_r32)
        line = f"  {rho:>6.2f}"
        for r in odd_r32:
            line += f"  {D[r]:>5.3f}"
        line += f"  {total_dev:>10.4f}"
        log(line)

    # ------------------------------------------------------------------
    # Scheme 3: forward-orbit visit count from inverse-tree nodes
    # For each odd node n at depth d, walk its forward orbit and count
    # the residues it visits. Weight each starting node uniformly.
    # ------------------------------------------------------------------
    log("\n" + "=" * 78)
    log("Scheme 3: forward-orbit visit count, weighted by depth window")
    log("=" * 78)
    log("  Each starting node n contributes its full forward-orbit residue counts.")
    log("  Sum over odd nodes weighted by 1 (uniform over inverse tree).")

    def forward_orbit_residues(n):
        """Return list of m mod 32 along forward odd-skeleton orbit until m=1."""
        residues = []
        m = n
        while m != 1:
            residues.append(m & 31)
            three_m = 3 * m + 1
            while (three_m & 1) == 0:
                three_m >>= 1
            m = three_m
        residues.append(1)  # final m=1
        return residues

    # Use odd nodes at d in [10, 30] as starting set (mid-depth, well-mixed)
    ranges_to_test = [(10, 30), (15, 25), (20, 30), (5, 20), (30, 50)]
    for d_min, d_max in ranges_to_test:
        starting = (odd_df.filter((pl.col("depth") >= d_min) &
                                  (pl.col("depth") <= d_max))
                          ["n"].to_list())
        log(f"\n  Starting set: odd nodes at d in [{d_min},{d_max}], "
            f"count={len(starting):,}")
        counts = {r: 0 for r in odd_r32}
        for n in starting:
            for r in forward_orbit_residues(n):
                if r in counts:
                    counts[r] += 1
        total = sum(counts.values())
        if total == 0:
            continue
        rho_visit = {r: counts[r] / total for r in odd_r32}
        D = {r: rho_visit[r] / pi_32[r] for r in odd_r32}
        total_dev = sum(abs(D_avg[r] - D[r]) for r in odd_r32)
        log(f"  {'r':>3}  {'D_avg':>7}  {'D_visit':>8}  {'diff':>8}")
        for r in odd_r32:
            log(f"  {r:>3}  {D_avg[r]:>7.4f}  {D[r]:>8.4f}  "
                f"{D_avg[r] - D[r]:>+8.4f}")
        log(f"  total |D_visit - D_avg| = {total_dev:.4f}")

    # ------------------------------------------------------------------
    # Scheme 4: forward-orbit visit count, but weighted to match late-t
    # The empirical trajectory measure conditions on survival to t large.
    # Each orbit of length L contributes 1/L visits per residue to residues
    # at "remaining time" L-t. For large t, only orbits with L > t contribute,
    # and they contribute weight at the tail of their orbit.
    # ------------------------------------------------------------------
    log("\n" + "=" * 78)
    log("Scheme 4: forward-orbit tail count (last K residues before m=1)")
    log("=" * 78)
    log("  For each starting node n, count residues in the LAST K steps of")
    log("  its forward orbit before reaching m=1. This isolates the descent")
    log("  endpoint structure.")

    starting = (odd_df.filter((pl.col("depth") >= 5) &
                              (pl.col("depth") <= 30))
                      ["n"].to_list())
    log(f"\n  Starting set: odd nodes at d in [5,30], count={len(starting):,}")

    for K_tail in [3, 5, 10, 15, 20, 30]:
        counts = {r: 0 for r in odd_r32}
        for n in starting:
            orbit = forward_orbit_residues(n)
            tail = orbit[-K_tail:] if len(orbit) >= K_tail else orbit
            for r in tail:
                if r in counts:
                    counts[r] += 1
        total = sum(counts.values())
        if total == 0:
            continue
        rho_tail = {r: counts[r] / total for r in odd_r32}
        D = {r: rho_tail[r] / pi_32[r] for r in odd_r32}
        total_dev = sum(abs(D_avg[r] - D[r]) for r in odd_r32)
        log(f"\n  K_tail={K_tail}: total |D - D_avg| = {total_dev:.4f}")
        line = f"  {'r':>3}: "
        for r in odd_r32:
            line += f" {D[r]:.2f}"
        log(line)
        line = f"  D_avg: "
        for r in odd_r32:
            line += f" {D_avg[r]:.2f}"
        log(line)

    (OUT / "inverse_tree_vs_davg_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
