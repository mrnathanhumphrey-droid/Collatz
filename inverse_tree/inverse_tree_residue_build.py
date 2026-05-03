"""
inverse_tree_residue_build.py

Brief: Characterize the inverse-Collatz tree's residue mod 2^k stationary
distribution and compare to the forward-orbit trajectory measure.

Method (extending Phase 5 Q3 from a★_6 partition to full mod 2^k partition):
  1. Load existing tree_d50.parquet (379,600 nodes, d=0..50)
  2. For each modulus M = 2^k, k in {5, 6, 8, 10, 11}:
     a. Tabulate residue density P_d(r mod M) at each depth d
     b. Build empirical transition matrix M_emp[r1, r2] = E[# child-r2 per parent-r1]
     c. Build closed-form transition matrix from inverse-Collatz rules:
          parent r --> child 2r mod M (always)
          parent r --> child (2r-1)/3 mod M (when 2r ≡ 1 mod 3 AND (2r-1)/3 odd, etc.)
          (Equivalently: from a forward perspective, parent in tree is forward-T(child),
           so backward step is the inverse: from p, doubling-child is 2p; (n-1)/3-child
           is (p-1)/3 when p ≡ 1 mod 3 AND (p-1)/3 odd AND (p-1)/3 > 1)
     d. Take leading left-eigenvector v of M_closed. Normalize to probability vector.
     e. Compare v to empirical P_d=50(r mod M).
  3. Compare inverse-tree stationary v to forward-orbit trajectory measure:
     - mod 32: spikes at r=5 (1.232), 17 (1.082), 29 (1.081); dips at 21, 19, 25
     - mod 2048: spike at r=341 (1.237)
     Test: does v match these spike ratios?

Outputs:
  inverse_tree_residue.csv -- P_d(r mod 2^k) tables and eigenvector predictions
  inverse_tree_residue_log.txt -- diagnostic log

Cap: <300 lines. Compute: should run in <2 minutes.
"""

import sys
import time
import numpy as np
import polars as pl
from pathlib import Path

OUT = Path(r"C:\Collatz\inverse_tree")
TREE_PARQUET = OUT / "tree_d50.parquet"
LOG = OUT / "inverse_tree_residue_log.txt"

sys.stdout.reconfigure(encoding="utf-8")
log_lines = []


def log(s):
    print(s)
    log_lines.append(s)


def is_eligible(p):
    """Inverse-3 step (p-1)/3 is valid iff p ≡ 1 mod 3 AND (p-1)/3 odd AND (p-1)/3 > 1."""
    if p % 3 != 1:
        return False
    m = (p - 1) // 3
    return m > 1 and m % 2 == 1


# ---------------------------------------------------------------------
# Load tree
# ---------------------------------------------------------------------
t0 = time.perf_counter()
df = pl.read_parquet(TREE_PARQUET)
log(f"[load] tree_d50.parquet -- {df.height:,} nodes, {df.width} columns, {time.perf_counter()-t0:.2f}s")
log(f"[load] depth range: {df['depth'].min()}..{df['depth'].max()}")
log(f"[load] schema: {df.columns}")


# ---------------------------------------------------------------------
# Compute n mod 2^k for each k of interest
# ---------------------------------------------------------------------
KS = [5, 6, 8, 10, 11]
for k in KS:
    M = 1 << k
    if f"r_{k}" not in df.columns:
        df = df.with_columns((pl.col("n") % M).alias(f"r_{k}"))


# ---------------------------------------------------------------------
# Empirical residue density at d=50 (deepest layer)
# ---------------------------------------------------------------------
DMAX = int(df["depth"].max())
deep = df.filter(pl.col("depth") == DMAX)
log(f"\n[empirical] deepest layer d={DMAX} has {deep.height:,} nodes")


# ---------------------------------------------------------------------
# Closed-form transition matrix on residues mod 2^k
# ---------------------------------------------------------------------
def build_closed_transition(k):
    """
    Build M_closed[r_from, r_to] over residues 0..2^k-1, natural-density approximation.

    Each parent integer n has:
      - one "left" child at 2n (always, T^-1 doubling step)
      - one "right" child at (n-1)/3, IFF n ≡ 4 mod 6 (the eligibility condition)

    At the residue level mod 2^k (k>=3), eligibility "n ≡ 4 mod 6" is NOT
    deterministic from n mod 2^k alone — it depends on n mod 6, which is
    independent of n mod 2^k for k>=2 (since gcd(2^k, 6) = 2).

    Within each EVEN residue class r mod 2^k, exactly 1/3 of the lifts
    satisfy n ≡ 4 mod 6 (specifically, those lifts r' with r' mod 6 = 4).
    Within each ODD residue class r mod 2^k, ZERO lifts are eligible
    (since n ≡ 4 mod 6 requires n even).

    Therefore the natural-density transition matrix is:
      M[r, 2r mod 2^k]   += 1                 (always, doubling child)
      M[r, child_r]      += 1/3   if r EVEN   (conditional, with prob 1/3)
    where child_r = ((r' - 1) * inv3) mod 2^k for any lift r' ≡ 4 mod 6
    of r mod 2^k. (All such qualifying lifts give the same child mod 2^k,
    since the mod-2^k child is invariant under r' -> r' + lcm(2^k, 6).)
    """
    Mod = 1 << k
    inv3 = pow(3, -1, Mod)  # inverse of 3 mod 2^k
    M_mat = np.zeros((Mod, Mod), dtype=np.float64)
    for r in range(Mod):
        # Left (doubling) child -- always
        r_left = (2 * r) % Mod
        M_mat[r, r_left] += 1.0
        # Right ((n-1)/3) child -- only EVEN residues, with weight 1/3
        if r % 2 == 0:
            # Find the qualifying lift r' ≡ 4 mod 6 within {r, r+2^k, r+2*2^k, ...}.
            # Since 2^k mod 6 ∈ {2, 4} for k>=1 (specifically 2^k mod 6 cycles 2,4,2,4,...),
            # we just iterate up to 3 lifts to find one with mod 6 == 4.
            for j in range(3):
                lift = r + j * Mod
                if lift % 6 == 4:
                    r_right = ((lift - 1) * inv3) % Mod
                    M_mat[r, r_right] += 1.0 / 3.0
                    break
    return M_mat


def leading_left_eigenvector(M_mat):
    """Leading left eigenvector of M (probability-normalized, |.| if needed)."""
    eigvals, eigvecs = np.linalg.eig(M_mat.T)
    order = np.argsort(-np.abs(eigvals.real))
    lambda_max = float(eigvals.real[order[0]])
    v = eigvecs[:, order[0]].real
    v = np.abs(v)
    v = v / v.sum()
    return lambda_max, v


# ---------------------------------------------------------------------
# Test stationarity at multiple depths (Step 3 of brief)
# ---------------------------------------------------------------------
log(f"\n=== Step 1+3: Residue density across depths (testing stationarity) ===")
log(f"")
log(f"For each (depth d, modulus k), compute chi^2/n_total relative to UNIFORM")
log(f"and KL divergence relative to the closed-form transition-matrix eigenvector.")
log(f"")

# Build closed-form predictions
predictions = {}
log(f"[closed-form] Building transition matrices and leading eigenvectors:")
for k in KS:
    Mod = 1 << k
    t0 = time.perf_counter()
    M_mat = build_closed_transition(k)
    lam, v = leading_left_eigenvector(M_mat)
    predictions[k] = (lam, v)
    log(f"  k={k:>2} (mod {Mod:>5}): lambda_max={lam:.6f} (= exp({np.log(lam):.4f})), "
        f"build+eig {time.perf_counter()-t0:.2f}s")

# Empirical at multiple depths
log(f"\n[empirical at depths] chi^2/n vs uniform AND KL vs closed-form prediction:")
log(f"  {'depth':>5} {'k':>2} {'n_total':>8} {'chi2/n_unif':>12} {'KL_vs_pred':>12}")

dist_rows = []
for d in range(20, DMAX + 1, 5):
    sub = df.filter(pl.col("depth") == d)
    n_total = sub.height
    if n_total < 200:
        continue
    for k in KS:
        Mod = 1 << k
        rk = f"r_{k}"
        # empirical density
        counts = np.zeros(Mod, dtype=np.float64)
        for r, c in zip(sub.group_by(rk).agg(pl.len().alias("c"))[rk].to_list(),
                        sub.group_by(rk).agg(pl.len().alias("c"))["c"].to_list()):
            counts[r] = c
        p_emp = counts / counts.sum()
        # chi^2 vs uniform
        expected_unif = n_total / Mod
        chi2_unif = ((counts - expected_unif) ** 2 / expected_unif).sum()
        chi2_per_n = chi2_unif / n_total
        # KL divergence vs closed-form prediction
        _, v_pred = predictions[k]
        eps = 1e-12
        kl = float(np.sum(p_emp * np.log((p_emp + eps) / (v_pred + eps))))
        log(f"  {d:>5} {k:>2} {n_total:>8} {chi2_per_n:>12.4f} {kl:>12.6f}")
        dist_rows.append({
            "depth": d, "k": k, "n_total": n_total,
            "chi2_per_n_unif": float(chi2_per_n),
            "KL_vs_predicted": float(kl),
        })

# ---------------------------------------------------------------------
# Step 4: Compare to forward-orbit trajectory measure
# ---------------------------------------------------------------------
log(f"\n=== Step 4: Compare to forward-orbit trajectory measure spikes ===")
log(f"")
log(f"From agent2_findings.md (3x+1 trajectory measure on iterates m at q=3):")
log(f"  m mod 32:")
log(f"    Over: 5 (1.232x), 17 (1.082x), 29 (1.081x), 23 (1.052x), 7 (1.051x),")
log(f"          15 (1.045x), 27 (1.040x), 11 (0.996x), 9 (0.981x), 3 (0.981x)")
log(f"    Under: 21 (0.887x), 19 (0.890x), 25 (0.896x), 13 (0.911x), 1 (0.930x), 31 (0.946x)")
log(f"  m mod 2048: 341 (1.237x)")
log(f"")

# Forward-orbit trajectory measure ratios from agent2_findings.md
fwd_mod32 = {5: 1.232, 17: 1.082, 29: 1.081, 23: 1.052, 7: 1.051,
             15: 1.045, 27: 1.040, 11: 0.996, 9: 0.981, 3: 0.981,
             21: 0.887, 19: 0.890, 25: 0.896, 13: 0.911, 1: 0.930, 31: 0.946}

# Inverse-tree closed-form prediction at mod 32
lam5, v5 = predictions[5]
# Restrict to ODD residues (forward Collatz Syracuse iterates m are always odd)
v5_odd = v5.copy()
for r in range(32):
    if r % 2 == 0:
        v5_odd[r] = 0.0
v5_odd = v5_odd / v5_odd.sum() if v5_odd.sum() > 0 else v5_odd
inv_density = {r: v5_odd[r] / (1.0 / 16) for r in range(32) if r % 2 == 1}  # ratio vs uniform-on-odd

log(f"[mod 32 odd residues only]")
log(f"  Inverse-tree closed-form predicted ratio (vs uniform-on-odd-residues):")
log(f"  {'r':>4} {'fwd_orbit':>10} {'inv_pred':>10} {'diff':>10}")
mod32_rows = []
for r in sorted(fwd_mod32.keys()):
    fwd = fwd_mod32[r]
    inv = inv_density[r]
    log(f"  {r:>4} {fwd:>10.4f} {inv:>10.4f} {inv-fwd:>+10.4f}")
    mod32_rows.append({"residue_mod_32": r, "fwd_orbit_ratio": fwd, "inv_tree_pred_ratio": inv,
                       "diff_inv_minus_fwd": inv - fwd})

# Empirical inverse-tree density at d=50, mod 32, ODD residues
deep_mod32 = deep["n"] % 32
empirical_mod32 = np.zeros(32)
for r in deep_mod32.to_list():
    empirical_mod32[r] += 1
empirical_mod32 /= empirical_mod32.sum()
empirical_mod32_odd = empirical_mod32.copy()
for r in range(32):
    if r % 2 == 0:
        empirical_mod32_odd[r] = 0.0
if empirical_mod32_odd.sum() > 0:
    empirical_mod32_odd /= empirical_mod32_odd.sum()
emp_density = {r: empirical_mod32_odd[r] / (1.0 / 16) for r in range(32) if r % 2 == 1}

log(f"")
log(f"[mod 32 odd residues] Inverse-tree EMPIRICAL at d=50 vs forward-orbit ratio:")
log(f"  {'r':>4} {'fwd_orbit':>10} {'inv_emp':>10} {'inv_pred':>10}")
for r in sorted(fwd_mod32.keys()):
    fwd = fwd_mod32[r]
    e = emp_density[r]
    p = inv_density[r]
    log(f"  {r:>4} {fwd:>10.4f} {e:>10.4f} {p:>10.4f}")

# Quantitative match: Pearson r between fwd_orbit and inv_pred ratios
fwd_arr = np.array([fwd_mod32[r] for r in sorted(fwd_mod32.keys())])
inv_pred_arr = np.array([inv_density[r] for r in sorted(fwd_mod32.keys())])
inv_emp_arr = np.array([emp_density[r] for r in sorted(fwd_mod32.keys())])
r_pred_vs_fwd = float(np.corrcoef(fwd_arr, inv_pred_arr)[0, 1])
r_emp_vs_fwd = float(np.corrcoef(fwd_arr, inv_emp_arr)[0, 1])
r_pred_vs_emp = float(np.corrcoef(inv_pred_arr, inv_emp_arr)[0, 1])
log(f"\n[Pearson correlations across 16 odd residues mod 32]")
log(f"  inv_predicted vs forward_orbit:  r = {r_pred_vs_fwd:+.4f}")
log(f"  inv_empirical vs forward_orbit:  r = {r_emp_vs_fwd:+.4f}")
log(f"  inv_predicted vs inv_empirical:  r = {r_pred_vs_emp:+.4f}")

# Also test the m ≡ 341 mod 2048 spike (forward = 1.237x)
lam11, v11 = predictions[11]
v11_odd = v11.copy()
for r in range(1 << 11):
    if r % 2 == 0:
        v11_odd[r] = 0.0
v11_odd_sum = v11_odd.sum()
if v11_odd_sum > 0:
    v11_odd /= v11_odd_sum

# Empirical mod 2048 from inverse tree at d=50
deep_mod2048 = deep["n"] % 2048
emp_mod2048 = np.zeros(2048)
for r in deep_mod2048.to_list():
    emp_mod2048[r] += 1
emp_mod2048 /= emp_mod2048.sum()
emp_mod2048_odd = emp_mod2048.copy()
for r in range(2048):
    if r % 2 == 0:
        emp_mod2048_odd[r] = 0.0
emp_mod2048_odd_sum = emp_mod2048_odd.sum()
if emp_mod2048_odd_sum > 0:
    emp_mod2048_odd /= emp_mod2048_odd_sum

ratio_341_pred = v11_odd[341] / (1.0 / 1024)
ratio_341_emp = emp_mod2048_odd[341] / (1.0 / 1024)
log(f"\n[mod 2048, residue 341 (the 'v=10 spike' from forward trajectory)]")
log(f"  Forward-orbit trajectory ratio: 1.237x")
log(f"  Inverse-tree closed-form pred:  {ratio_341_pred:.4f}x")
log(f"  Inverse-tree empirical (d=50):  {ratio_341_emp:.4f}x")
log(f"  (Note: empirical at d=50 has only {deep.height:,} nodes; mod 2048 sub-cells small)")

# ---------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------
pl.DataFrame(dist_rows).write_csv(OUT / "inverse_tree_residue.csv")
pl.DataFrame(mod32_rows).write_csv(OUT / "inverse_tree_mod32_compare.csv")

# Save full closed-form eigenvector for k=5 (mod 32) and k=6 (mod 64) for the doc
v5_full_rows = [{"residue_mod_32": r, "predicted_density": float(v5[r]),
                 "is_odd": (r % 2 == 1)} for r in range(32)]
pl.DataFrame(v5_full_rows).write_csv(OUT / "inverse_tree_eigvec_mod32.csv")

LOG.write_text("\n".join(log_lines), encoding="utf-8")
log(f"\n[wrote] inverse_tree_residue.csv ({len(dist_rows)} rows)")
log(f"[wrote] inverse_tree_mod32_compare.csv ({len(mod32_rows)} rows)")
log(f"[wrote] inverse_tree_eigvec_mod32.csv (32 rows)")
log(f"[wrote] inverse_tree_residue_log.txt")
LOG.write_text("\n".join(log_lines), encoding="utf-8")
