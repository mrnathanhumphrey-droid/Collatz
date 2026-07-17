"""
independence_audit_compute.py — empirical checks for Validation Task 2.

Two key cross-checks:
  Claim 3 — Pearson(R58 prediction, R60 prediction)
  Claim 7 — algebraic D_1 vs 2*log(lambda_max)/log(2)
"""
import math
import csv
from pathlib import Path
import numpy as np

OUT = Path(r"C:\Collatz")


# ---------------- Claim 3: Pearson(R58, R60) ---------------------------
def load_r58():
    """R58 predictions from inverse-tree variants."""
    rows = []
    with open(OUT / "experiments_output" / "inverse_tree_predicted_D.csv") as f:
        rd = csv.DictReader(f)
        for row in rd:
            rows.append({k: float(v) if k != "r" else int(v)
                         for k, v in row.items()})
    return rows


def load_r60():
    rows = []
    with open(OUT / "size_stratified_residue_marginal.csv") as f:
        rd = csv.DictReader(f)
        for row in rd:
            rows.append({k: float(v) if k != "r" else int(v)
                         for k, v in row.items()})
    return rows


def pearson(xs, ys):
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    num = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
    den = math.sqrt(sum((x-mx)**2 for x in xs) *
                    sum((y-my)**2 for y in ys))
    return num/den if den > 0 else 0.0


def claim3():
    print("=" * 78)
    print("Claim 3 — Pearson(R58 prediction, R60 prediction)")
    print("=" * 78)
    r58 = load_r58()
    r60 = load_r60()
    r58_by_r = {row["r"]: row for row in r58}
    r60_by_r = {row["r"]: row for row in r60}
    rs = sorted(set(r58_by_r) & set(r60_by_r))
    print(f"\n  Common residues: {len(rs)}")

    # R58: variant_a is the primary subtree-size weighting
    # R60: D_pred from size-stratified Markov
    # R50/R59 D_avg empirical column for sanity
    r58_var_a = [r58_by_r[r]["variant_a"] for r in rs]
    r58_emp = [r58_by_r[r]["empirical"] for r in rs]
    r60_pred = [r60_by_r[r]["D_pred"] for r in rs]
    r60_emp = [r60_by_r[r]["D_avg"] for r in rs]

    print(f"\n  Pearson(R58 variant_a, R60 D_pred)        = "
          f"{pearson(r58_var_a, r60_pred):+.4f}")
    print(f"  Pearson(R58 variant_a, R58 empirical D)   = "
          f"{pearson(r58_var_a, r58_emp):+.4f}")
    print(f"  Pearson(R60 D_pred,   R60 empirical D_avg)= "
          f"{pearson(r60_pred, r60_emp):+.4f}")
    print(f"  Pearson(R58 variant_a, R60 empirical D_avg)= "
          f"{pearson(r58_var_a, r60_emp):+.4f}")

    # Are R58 empirical and R60 empirical the same data source?
    print(f"  Pearson(R58 empirical, R60 empirical)     = "
          f"{pearson(r58_emp, r60_emp):+.4f}")
    print(f"\n  R58 empirical (D at t=90): " + " ".join(f"{v:.3f}" for v in r58_emp))
    print(f"  R60 empirical (D_avg late-t): " + " ".join(f"{v:.3f}" for v in r60_emp))

    # Per-residue MAE
    diff58 = [abs(p - e) for p, e in zip(r58_var_a, r60_emp)]
    diff60 = [abs(p - e) for p, e in zip(r60_pred, r60_emp)]
    cross = [abs(a - b) for a, b in zip(r58_var_a, r60_pred)]
    print(f"\n  MAE(R58 - D_avg) = {sum(diff58)/len(diff58):.4f}")
    print(f"  MAE(R60 - D_avg) = {sum(diff60)/len(diff60):.4f}")
    print(f"  MAE(R58 - R60)   = {sum(cross)/len(cross):.4f}")


# ---------------- Claim 7: algebraic dim relations -------------
def claim7():
    print("\n" + "=" * 78)
    print("Claim 7 — Are the various dim values algebraically related?")
    print("=" * 78)

    LOG2 = math.log(2.0)

    # R23 lambda_max
    lam_max = 1.263763
    branching_dim = math.log(lam_max) / LOG2
    twice = 2 * branching_dim
    print(f"\n  R23 lambda_max                 = {lam_max:.6f}")
    print(f"  log(lam_max)/log(2)            = {branching_dim:.6f}  (Furstenberg branching dim)")
    print(f"  2 * log(lam_max)/log(2)        = {twice:.6f}  (claimed = 0.6755)")

    # Chang exact
    phi = (1 + math.sqrt(5)) / 2
    chang_exact = math.log(phi) / LOG2
    print(f"  Chang exact = log(phi)/log(2)  = {chang_exact:.6f}  (= log_2(phi))")

    # R61 spatial info dim D_1 from saved CSV
    D_1 = None
    try:
        with open(OUT / "mfdfa_spatial_tau_q.csv") as f:
            rd = csv.DictReader(f)
            for row in rd:
                if abs(float(row["q"]) - 1.0) < 1e-6:
                    D_1 = float(row["D_q"])
    except FileNotFoundError:
        pass
    print(f"\n  R61 spatial info dim D_1       = {D_1:.6f}")

    print(f"\n  Algebraic identity tests:")
    print(f"    D_1 == 2 * log(lam_max)/log(2)?  {D_1:.4f} vs {twice:.4f}  diff = {abs(D_1 - twice):.4f}")
    print(f"      -> {'YES (within 0.01)' if abs(D_1 - twice) < 0.01 else 'NO — coincidence'}")
    print(f"    D_1 == log(phi)/log(2)?           {D_1:.4f} vs {chang_exact:.4f}  diff = {abs(D_1 - chang_exact):.4f}")
    print(f"      -> {'YES (within 0.01)' if abs(D_1 - chang_exact) < 0.01 else 'NO — distinct'}")
    print(f"    2*log(lam_max) == log(phi)?       {twice:.4f} vs {chang_exact:.4f}  diff = {abs(twice - chang_exact):.4f}")
    print(f"      -> {'YES (within 0.01)' if abs(twice - chang_exact) < 0.01 else 'NO — distinct'}")
    print(f"    Chang heuristic 0.68 vs D_1?      0.6800 vs {D_1:.4f}  diff = {abs(0.68 - D_1):.4f}")
    print(f"      -> {'within 0.10' if abs(0.68 - D_1) < 0.10 else 'distinct'}")

    print(f"\n  Conclusion:")
    print(f"    Furstenberg branching dim (0.34), info dim D_1 (~0.61),")
    print(f"    Chang exact (0.69), Chang heuristic (0.68), R23 derived 2x (0.68)")
    print(f"    are five DIFFERENT numbers near 0.34/0.68 — no exact algebraic")
    print(f"    identity holds among them at 4-decimal precision. The 'multiple")
    print(f"    methods agree' claim is at the 0.07-wide 0.6-0.7 region, not")
    print(f"    a sharp identification.")


if __name__ == "__main__":
    claim3()
    claim7()
