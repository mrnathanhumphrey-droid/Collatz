"""
higher_q_partition_test.py — Result 65 candidate.

Test whether the R63 q=3 Fourier framework generalizes to higher q ∈ {5, 7, 9,
11, 13, 15}. Key empirical questions:

1. The "predicted vs empirical" comparison is trivially exact: the mass-fraction
   formula and direct-Fourier sum compute the same thing (algebraic identity).
   So we report resonance VALUES at each (a, q), not fit residuals.

2. Are |μ̂(a/q)|^2 at small q substantially > irrational-ξ baseline (~1/N)?
   That's the real "is there a resonance" test.

3. How does mass-fraction asymmetry scale with q? At q=3 we have (a,b,c) =
   (0.007, 0.347, 0.646) — extreme. At higher q we expect the asymmetry to
   shrink as residue classes proliferate.

4. Mechanism: how does inverse-Syracuse predecessor structure depend on m mod q?
   For q=3 the mechanism is clear (m mod 3 controls v selection in pred = (m·2^v−1)/3).
   For other q, document via predecessor-size analysis.
"""
import math
import sys
import time
from collections import deque, defaultdict
from pathlib import Path
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(r"C:\Collatz")
log_lines = []
def log(s):
    print(s, flush=True)
    log_lines.append(s)


# ----------------------------------------------------------------
# Inverse tree builder (consistent with R58/R61)
# ----------------------------------------------------------------
def build_inverse_tree(max_value):
    """BFS from m=1 backward. Returns dict m -> {parent, depth}.
    Includes BOTH even and odd m. Subtree-size weights computed only for odd m."""
    tree = {1: {'parent': None, 'depth': 0}}
    q = deque([1])
    while q:
        m = q.popleft()
        d = tree[m]['depth']
        if m % 3 == 0:
            continue
        v_start = 2 if (m % 3 == 1) else 1
        for v in range(v_start, 64, 2):
            num = m * (1 << v) - 1
            if num <= 0:
                continue
            pred = num // 3
            if pred > max_value:
                break
            if pred & 1 == 0 or pred == m:
                continue
            if pred not in tree:
                tree[pred] = {'parent': m, 'depth': d + 1}
                q.append(pred)
    return tree


def compute_subtree_sizes(tree):
    children = defaultdict(list)
    for m, info in tree.items():
        if info['parent'] is not None:
            children[info['parent']].append(m)
    order = []
    stack = [1]
    visited = set()
    while stack:
        m = stack[-1]
        if m in visited:
            stack.pop()
            order.append(m)
            continue
        visited.add(m)
        for c in children.get(m, []):
            stack.append(c)
    sizes = {}
    for m in order:
        sizes[m] = 1 + sum(sizes[c] for c in children.get(m, []))
    return sizes


# ----------------------------------------------------------------
# Mass-fraction formula and direct Fourier
# ----------------------------------------------------------------
def mass_fractions(m_arr, w_arr, q):
    f = np.zeros(q, dtype=np.float64)
    for r in range(q):
        f[r] = w_arr[(m_arr % q) == r].sum()
    return f / f.sum()


def closed_form_mu_hat_sq(f, a, q):
    """|mu_hat(a/q)|^2 via mass-fraction formula (predicted)."""
    omega = np.exp(2j * np.pi * a / q)
    mu_hat = sum(f[r] * omega**r for r in range(q))
    return float(abs(mu_hat) ** 2)


def empirical_mu_hat_sq(m_arr, w_arr, a, q):
    """Direct |mu_hat(a/q)|^2 from sum over integers."""
    phases = (a * m_arr) % q
    omega = np.exp(2j * np.pi * phases / q)
    mu_hat = (w_arr * omega).sum() / w_arr.sum()
    return float(abs(mu_hat) ** 2)


def empirical_mu_hat_sq_irrational(m_arr, w_arr, xi):
    """Direct |mu_hat(xi)|^2 at arbitrary real xi."""
    mu_hat = (w_arr * np.exp(2j * np.pi * xi * m_arr)).sum() / w_arr.sum()
    return float(abs(mu_hat) ** 2)


def main():
    log("=" * 78)
    log("Higher-q partition test — does R63 q=3 framework generalize?")
    log("=" * 78)

    # ------------------ Build tree ----------------
    max_value = 1 << 22  # match R58 scale
    log(f"\nStep 1: build inverse tree at max_value = 2^22 = {max_value}")
    t0 = time.time()
    tree = build_inverse_tree(max_value)
    log(f"  tree size: {len(tree):,} nodes, {time.time()-t0:.1f}s")
    sizes = compute_subtree_sizes(tree)

    odd_m = np.array([m for m in tree if m & 1], dtype=np.int64)
    odd_w = np.array([float(sizes[m]) for m in odd_m], dtype=np.float64)
    log(f"  odd nodes: {len(odd_m):,}, total weight: {odd_w.sum():.4e}")

    # ------------------ Mass fractions per q ------------
    qs_test = [3, 5, 7, 9, 11, 13, 15]
    log("\n" + "=" * 78)
    log("Step 2: mass-fractions f_r at each q")
    log("=" * 78)

    fractions_by_q = {}
    for q in qs_test:
        f = mass_fractions(odd_m, odd_w, q)
        fractions_by_q[q] = f
        # max-min spread
        spread = f.max() - f.min()
        # Effective number of "active" residues (entropy-based)
        H = -np.sum(f * np.log(f + 1e-30)) / math.log(q)  # normalized entropy
        log(f"\n  q = {q}: mass-fractions f_r = " +
            ", ".join(f"{f[r]:.4f}" for r in range(q)))
        log(f"    spread max−min = {spread:.4f}, H_q (normalized) = {H:.4f} (1=uniform)")

        # Identify near-zero residues (mass < 1%)
        zero_residues = [r for r in range(q) if f[r] < 0.01]
        if zero_residues:
            log(f"    near-zero residues (f<1%): {zero_residues}")

    # ------------------ Resonances |mu_hat(a/q)|^2 ----
    log("\n" + "=" * 78)
    log("Step 3: |mu_hat(a/q)|^2 at all (a, q) with gcd(a, q) = 1, a ≤ q/2")
    log("=" * 78)

    resonances = []  # (q, a, |mu|^2 closed-form, |mu|^2 direct)
    for q in qs_test:
        f = fractions_by_q[q]
        log(f"\n  q = {q}:")
        for a in range(1, (q + 1) // 2 + 1):
            if math.gcd(a, q) != 1:
                continue
            mu_sq_pred = closed_form_mu_hat_sq(f, a, q)
            mu_sq_emp = empirical_mu_hat_sq(odd_m, odd_w, a, q)
            ratio = mu_sq_pred / mu_sq_emp if mu_sq_emp > 0 else np.nan
            log(f"    a={a}/{q}: pred={mu_sq_pred:.6e}, emp={mu_sq_emp:.6e}, "
                f"ratio={ratio:.6f}")
            resonances.append((q, a, mu_sq_pred, mu_sq_emp))

    # ------------------ Irrational baseline -----------
    log("\n" + "=" * 78)
    log("Step 4: irrational-xi baseline (random non-rational ξ)")
    log("=" * 78)

    # Pick 8 irrational xi values not close to any small rationals
    rng = np.random.default_rng(20260503)
    irrationals = []
    while len(irrationals) < 8:
        xi = rng.uniform(0.05, 0.45)
        # Reject if close to a/q for any q ≤ 16
        too_close = False
        for q in range(2, 17):
            for a in range(1, q):
                if abs(xi - a/q) < 0.01:
                    too_close = True
                    break
            if too_close:
                break
        if not too_close:
            irrationals.append(xi)

    log(f"\n  {'xi':>10}  {'|mu_hat|^2':>14}")
    irr_results = []
    for xi in irrationals:
        v = empirical_mu_hat_sq_irrational(odd_m, odd_w, xi)
        irr_results.append((xi, v))
        log(f"  {xi:>10.6f}  {v:>14.4e}")
    irr_baseline = float(np.median([v for _, v in irr_results]))
    log(f"\n  median irrational baseline |mu_hat|^2 = {irr_baseline:.4e}")
    log(f"  (~ 1/N noise floor expected at N = {len(odd_m):,}; "
        f"1/N = {1/len(odd_m):.4e})")

    # ------------------ Summary by q ---------------
    log("\n" + "=" * 78)
    log("Step 5: resonance summary by q vs baseline")
    log("=" * 78)

    log(f"\n  Median irrational baseline = {irr_baseline:.4e}")
    log(f"\n  {'q':>3}  {'a':>3}  {'|mu_hat|^2':>14}  {'/baseline':>10}  signif?")
    for q in qs_test:
        f = fractions_by_q[q]
        for a in range(1, (q + 1) // 2 + 1):
            if math.gcd(a, q) != 1:
                continue
            mu_sq = closed_form_mu_hat_sq(f, a, q)
            ratio_to_baseline = mu_sq / irr_baseline
            sig = ("STRONG" if ratio_to_baseline > 1000 else
                   "moderate" if ratio_to_baseline > 100 else
                   "weak" if ratio_to_baseline > 10 else
                   "noise")
            log(f"  {q:>3}  {a:>3}  {mu_sq:>14.4e}  {ratio_to_baseline:>10.1e}  {sig}")

    # Headline values
    log(f"\n  Headline |mu_hat(a/q)|^2 at smallest a:")
    log(f"  {'q':>3}  {'a':>3}  {'|mu|^2':>10}")
    for q in qs_test:
        f = fractions_by_q[q]
        for a in range(1, q):
            if math.gcd(a, q) == 1:
                mu_sq = closed_form_mu_hat_sq(f, a, q)
                log(f"  {q:>3}  {a:>3}  {mu_sq:>10.4f}")
                break

    # ------------------ Step 6: q-dependence analysis -------
    log("\n" + "=" * 78)
    log("Step 6: how does mass-asymmetry scale with q?")
    log("=" * 78)
    log(f"\n  {'q':>3}  {'spread max-min':>15}  {'normalized H_q':>15}  "
        f"{'effective uniform?':>20}  q ≡ 0 mod 3?")
    for q in qs_test:
        f = fractions_by_q[q]
        spread = f.max() - f.min()
        H = -np.sum(f * np.log(f + 1e-30)) / math.log(q)
        uniform_flag = "ALMOST" if H > 0.99 else "asymmetric" if H > 0.7 else "STRONG"
        is_3div = "YES" if q % 3 == 0 else "no"
        log(f"  {q:>3}  {spread:>15.4f}  {H:>15.4f}  {uniform_flag:>20}  {is_3div}")

    # ------------------ Step 7: predecessor mechanism ----
    log("\n" + "=" * 78)
    log("Step 7: predecessor mechanism analysis per q")
    log("=" * 78)
    log("\n  For each q, examine: does inverse-Syracuse pred(m) systematically")
    log("  depend on m mod q? Specifically: for each r mod q, what is the average")
    log("  v in the smallest-pred relation pred = (m·2^v − 1)/3, restricted to m ≡ r?")
    log("  And what fraction of m ≡ r mod q have NO valid predecessors (leaves)?")
    log("")
    log(f"  {'q':>3}  {'r':>3}  {'#odd-m (≡r mod q)':>18}  {'#leaves':>9}  "
        f"{'leaf-frac':>10}  {'avg subtree size':>17}")
    for q in qs_test:
        if q > 9:  # skip detail for largest q
            continue
        for r in range(q):
            mask = (odd_m % q) == r
            if mask.sum() == 0:
                continue
            sub_m = odd_m[mask]
            sub_w = odd_w[mask]
            # leaf = subtree_size == 1
            n_leaves = int((sub_w == 1).sum())
            leaf_frac = n_leaves / len(sub_m)
            avg_size = float(sub_w.mean())
            log(f"  {q:>3}  {r:>3}  {len(sub_m):>18,}  {n_leaves:>9,}  "
                f"{leaf_frac:>10.4f}  {avg_size:>17.2f}")

    # ------------------ Outputs ----------------
    log("\n" + "=" * 78)
    log("Saving outputs")
    log("=" * 78)

    with open(OUT / "mass_fractions_by_q.csv", "w") as ff:
        ff.write("q,r,f_r\n")
        for q in qs_test:
            f = fractions_by_q[q]
            for r in range(q):
                ff.write(f"{q},{r},{f[r]:.8f}\n")
    log("  [wrote] mass_fractions_by_q.csv")

    with open(OUT / "fourier_predictions_by_q.csv", "w") as ff:
        ff.write("q,a,mu_hat_sq_predicted,mu_hat_sq_empirical,ratio\n")
        for q, a, p, e in resonances:
            r = p / e if e > 0 else float('nan')
            ff.write(f"{q},{a},{p:.8e},{e:.8e},{r:.8f}\n")
    log("  [wrote] fourier_predictions_by_q.csv")

    with open(OUT / "irrational_xi_baseline.csv", "w") as ff:
        ff.write("xi,mu_hat_sq\n")
        for xi, v in irr_results:
            ff.write(f"{xi:.8f},{v:.8e}\n")
    log("  [wrote] irrational_xi_baseline.csv")

    # ------------------ Verdict ----------------
    log("\n" + "=" * 78)
    log("VERDICT")
    log("=" * 78)

    # Determine outcome based on which q have resonances above 100x baseline
    res_above_100 = defaultdict(list)
    for q, a, p, _ in resonances:
        if p / irr_baseline > 100:
            res_above_100[q].append(a)

    qs_with_resonance = sorted(res_above_100.keys())
    qs_3_divisible = [q for q in qs_with_resonance if q % 3 == 0]
    qs_not_3_div = [q for q in qs_with_resonance if q % 3 != 0]
    log(f"\n  q with resonance > 100x baseline: {qs_with_resonance}")
    log(f"  of those, q divisible by 3: {qs_3_divisible}")
    log(f"  of those, q NOT divisible by 3: {qs_not_3_div}")

    if set(qs_with_resonance) == set(qs_test):
        log(f"\n  Outcome (α): framework generalizes to all tested q.")
    elif set(qs_with_resonance) == set([q for q in qs_test if q % 3 == 0]):
        log(f"\n  Outcome (β): framework works for q ∈ multiples of 3 only.")
    elif len(qs_with_resonance) == 1 and 3 in qs_with_resonance:
        log(f"\n  Outcome (δ): framework q=3 specific.")
    else:
        log(f"\n  Outcome (γ): framework works at some q, not others — see table.")

    (OUT / "higher_q_partition_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )
    log("\n  [wrote] higher_q_partition_log.txt")


if __name__ == "__main__":
    main()
