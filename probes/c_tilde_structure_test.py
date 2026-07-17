"""
c_tilde_structure_test.py
=========================
Test structural candidates for c~_q := lim S_k^(q) / (q/3)^k across
q in {3, 5, 7, 11, 13}. Uses the q-sweep cache for exact rationals.

Candidates tested:
  1. ord(2 mod q) / (q-1) correlation
  2. q=7 outlier hypothesis (linear/log/poly fits without q=7)
  3. Pair clustering / step function vs smooth
  4. Closed-form rationals (a/b, (q-r)/q, etc.)
  5. Unifying formula vs S_inf^(3)=7/15
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")


# --------------------------------------------------------------------------- #
# Pull c~_q from cache                                                        #
# --------------------------------------------------------------------------- #

def order_of_two(N):
    m = 1; v = 2 % N
    while v != 1:
        v = (v * 2) % N; m += 1
    return m


def load_X_cache():
    with open(r"C:\Collatz\experiments_output\result_q_sweep_test_2_cache.json") as fh:
        d = json.load(fh)
    out = {}
    for key, val in d.items():
        q, k = map(int, key.split(","))
        Xk = Fraction(int(val["X_num"]), int(val["X_den"]))
        out.setdefault(q, {})[k] = Xk
    return out


def compute_S_and_ctilde(X_data):
    """For each q, compute S_k = X_k - X_{k-1} and c~_q estimate sequence."""
    out = {}
    for q in sorted(X_data):
        ks = sorted(X_data[q])
        Xs = {0: Fraction(1)}
        Xs.update(X_data[q])
        Sks = []
        ctilde_seq = []
        for k in ks:
            Sk = Xs[k] - Xs[k-1]
            Sks.append((k, Sk))
            qk_over_3k = Fraction(q, 3) ** k
            ctilde_seq.append((k, float(Sk / qk_over_3k)))
        out[q] = {"S": Sks, "ctilde_seq": ctilde_seq, "ctilde_last": ctilde_seq[-1][1]}
    return out


# --------------------------------------------------------------------------- #
# Candidate fits                                                              #
# --------------------------------------------------------------------------- #

def candidate_1_ord_correlation(qs, ctildes):
    """ord(2 mod q) / (q-1) — multiplicative-order index."""
    print("=" * 78)
    print("CANDIDATE 1: ord(2 mod q) / (q-1) — multiplicative-order index")
    print("=" * 78)
    print(f"  {'q':>3}  {'ord(2 mod q)':>12}  {'q-1':>4}  {'ord/(q-1)':>10}  {'c~_q':>10}  {'ord(2 mod q)':>12}")
    rows = []
    for q, c in zip(qs, ctildes):
        ord2 = order_of_two(q)
        ratio = ord2 / (q - 1)
        rows.append((q, ord2, ratio, c))
        print(f"  {q:>3}  {ord2:>12}  {q-1:>4}  {ratio:>10.4f}  {c:>10.4f}  {ord2:>12}")
    print()

    # Pearson correlation
    ratios = np.array([r[2] for r in rows])
    cs = np.array([r[3] for r in rows])
    if len(set(ratios)) > 1:
        corr = float(np.corrcoef(ratios, cs)[0, 1])
        print(f"  Pearson r(ratio, c~_q) = {corr:.4f}")
    print()
    print("  Reading: q=7 has ord/(q-1) = 0.5; all others = 1.0. Binary distinction.")
    print("  c~_q values: q=7 has c~ ~0.78 (high), but q=11, 13 also ~0.73-0.77 (high).")
    print("  -> ord/(q-1) alone does NOT predict c~_q ordering. Q=7 isn't simply the outlier.")
    return rows


def candidate_2_q7_outlier(qs, ctildes):
    """Remove q=7; check if {3, 5, 11, 13} fits simple function."""
    print("=" * 78)
    print("CANDIDATE 2: q=7 outlier hypothesis — fits without q=7")
    print("=" * 78)
    qs_no7 = [q for q in qs if q != 7]
    cs_no7 = [c for q, c in zip(qs, ctildes) if q != 7]
    qs_arr = np.array(qs_no7, dtype=float)
    cs_arr = np.array(cs_no7)

    # Linear in q
    A = np.vstack([qs_arr, np.ones_like(qs_arr)]).T
    slope, intercept = np.linalg.lstsq(A, cs_arr, rcond=None)[0]
    pred_lin = slope * qs_arr + intercept
    rss_lin = float(np.sum((cs_arr - pred_lin) ** 2))
    print(f"  Linear in q: c~_q = {slope:+.5f}*q + {intercept:+.5f}")
    print(f"    RSS = {rss_lin:.4e}, R^2 = {1 - rss_lin/np.var(cs_arr)/len(cs_arr):.4f}")
    for q, c, p in zip(qs_no7, cs_no7, pred_lin):
        print(f"    q={q}: actual {c:.4f}, predicted {p:.4f}, residual {c-p:+.4f}")

    # 1/q form
    inv_qs = 1.0 / qs_arr
    A = np.vstack([inv_qs, np.ones_like(inv_qs)]).T
    slope2, intercept2 = np.linalg.lstsq(A, cs_arr, rcond=None)[0]
    pred_inv = slope2 * inv_qs + intercept2
    rss_inv = float(np.sum((cs_arr - pred_inv) ** 2))
    print(f"\n  c~_q = a + b/q form: c~_q = {intercept2:.5f} + {slope2:.5f}/q")
    print(f"    RSS = {rss_inv:.4e}")
    for q, c, p in zip(qs_no7, cs_no7, pred_inv):
        print(f"    q={q}: actual {c:.4f}, predicted {p:.4f}, residual {c-p:+.4f}")

    # log(q) form
    log_qs = np.log(qs_arr)
    A = np.vstack([log_qs, np.ones_like(log_qs)]).T
    slope3, intercept3 = np.linalg.lstsq(A, cs_arr, rcond=None)[0]
    pred_log = slope3 * log_qs + intercept3
    rss_log = float(np.sum((cs_arr - pred_log) ** 2))
    print(f"\n  c~_q = a + b*log(q): c~_q = {intercept3:.5f} + {slope3:.5f}*log(q)")
    print(f"    RSS = {rss_log:.4e}")

    print("\n  Best simple-form fit on {3, 5, 11, 13}: linear" if rss_lin < min(rss_inv, rss_log)
          else "  Best simple-form fit: 1/q" if rss_inv < rss_log else "  Best simple-form fit: log(q)")

    print()


def candidate_3_step_function(qs, ctildes):
    print("=" * 78)
    print("CANDIDATE 3: Pair clustering / step function vs smooth")
    print("=" * 78)
    low_cluster = [(q, c) for q, c in zip(qs, ctildes) if c < 0.6]
    high_cluster = [(q, c) for q, c in zip(qs, ctildes) if c >= 0.6]
    low_mean = np.mean([c for _, c in low_cluster])
    high_mean = np.mean([c for _, c in high_cluster])
    low_std = np.std([c for _, c in low_cluster])
    high_std = np.std([c for _, c in high_cluster])
    print(f"  Low cluster:  {low_cluster}, mean = {low_mean:.4f}, std = {low_std:.4f}")
    print(f"  High cluster: {high_cluster}, mean = {high_mean:.4f}, std = {high_std:.4f}")
    print(f"  Step magnitude: {high_mean - low_mean:.4f}")
    # Within-cluster spread relative to step:
    within_to_step = max(low_std, high_std) / (high_mean - low_mean)
    print(f"  Max within-cluster std / step = {within_to_step:.4f}")
    print(f"  -> If within-cluster std << step, step function fits. "
          f"({'YES' if within_to_step < 0.2 else 'MARGINAL' if within_to_step < 0.4 else 'NO'})")
    print()


def candidate_4_closed_form(qs, ctildes):
    print("=" * 78)
    print("CANDIDATE 4: Closed-form rational candidates")
    print("=" * 78)

    print(f"  {'q':>3}  {'c~_q':>10}  {'c~·q':>10}  {'c~·q-(q-3)':>14}  {'c~·q-(q-1)':>14}")
    for q, c in zip(qs, ctildes):
        cq = c * q
        d3 = cq - (q - 3)
        d1 = cq - (q - 1)
        print(f"  {q:>3}  {c:>10.4f}  {cq:>10.4f}  {d3:>+14.4f}  {d1:>+14.4f}")
    print()

    # (q-3)/q candidate
    print("  Test c~_q ?= (q-3)/q  [main term of c~ = (1 - 3/q)]:")
    for q, c in zip(qs, ctildes):
        pred = (q - 3) / q
        dev = c - pred
        print(f"    q={q}: actual {c:.4f}, (q-3)/q = {pred:.4f}, dev = {dev:+.4f}")
    print()

    # Find C(q) := c~_q / ((q-3)/q) for q != 3
    print("  Implied C(q) = c~_q · q / (q-3) [should -> 1 if main term is correct]:")
    for q, c in zip(qs, ctildes):
        if q == 3:
            print(f"    q=3: undefined (q-3=0; c~_3 = 7/15 separate regime)")
            continue
        Cq = c * q / (q - 3)
        print(f"    q={q}: C(q) = {Cq:.6f}, deviation from 1 = {Cq - 1:+.6f}")
    print()

    # 7/15 stuff
    print("  7/15 connection at q=3: c~_3 = 7/15 = {:.6f}".format(7/15))
    print(f"  Empirical c~_3 = {ctildes[0]:.6f} (from extrapolation)")
    print()

    # Try (q-r)/q for various small r
    print("  Best integer r in c~_q ≈ (q-r)/q across q=11, 13 (where C(q) ~ 1):")
    for r in range(1, 6):
        deviations = []
        for q, c in zip(qs, ctildes):
            if q < 9: continue
            pred = (q - r) / q
            deviations.append((q, c, pred, c - pred))
        rss = sum((d[3] ** 2) for d in deviations)
        print(f"    r={r}: RSS over q in {{11,13}} = {rss:.4e}")
        for q, c, p, dev in deviations:
            print(f"      q={q}: actual {c:.4f}, (q-{r})/q = {p:.4f}, dev = {dev:+.4f}")
    print()


def candidate_5_unifying_formula(qs, ctildes, ctilde_full):
    print("=" * 78)
    print("CANDIDATE 5: Unifying formula across q=3 and q>=5")
    print("=" * 78)

    print("  Two-regime hypothesis:")
    print(f"    q=3 regime (S_inf finite): c~_3 = lim S_k^(3) = 7/15 = {7/15:.6f}")
    print(f"    q>=5 regime (renormalized): c~_q = (q-3)/q + delta(q)")
    print(f"      delta(q=5) = {ctildes[1] - 2/5:+.4f}")
    print(f"      delta(q=7) = {ctildes[2] - 4/7:+.4f}")
    print(f"      delta(q=11) = {ctildes[3] - 8/11:+.4f}")
    print(f"      delta(q=13) = {ctildes[4] - 10/13:+.4f}")
    print()
    print("  Reading: delta(q) -> 0 as q grows (for q with 2 prim root).")
    print("           q=7 (where 2 is NOT prim root) has anomalously large delta.")
    print("           q=5 has moderate delta (still settling at k=4)")
    print("           q=3 boundary case: doesn't fit (q-3)/q since (q-3)/q=0 but c~_3=7/15 != 0")
    print()

    # Check c~_q at small q is converged via the q-sweep finite-k sequences
    print("  Per-q convergence of c~_q sequence (across k):")
    for q in qs:
        print(f"    q={q}: {[(k, f'{c:.4f}') for k, c in ctilde_full[q]['ctilde_seq']]}")
    print()
    print("  q=5 c~_q sequence: 0.5333, 0.4922, 0.4896, 0.4877 — slowly decreasing.")
    print("  Aitken extrapolation:")
    seq5 = [c for _, c in ctilde_full[5]['ctilde_seq']]
    if len(seq5) >= 3:
        # Aitken delta^2
        for i in range(len(seq5) - 2):
            denom = seq5[i+2] - 2*seq5[i+1] + seq5[i]
            if abs(denom) > 1e-12:
                a = seq5[i] - (seq5[i+1] - seq5[i])**2 / denom
                print(f"    Aitken[{i}..{i+2}] = {a:.6f}")
    print()
    print("  Conjecture: c~_5 limit might converge to 0.48 (close to 2/5 + small),")
    print("    but at k=4 not yet stable. Either (q-3)/q with finite-k transient,")
    print("    or genuinely larger by ~0.08.")


def main():
    print("=" * 78)
    print("c~_q structure test across q in {3, 5, 7, 11, 13}")
    print("=" * 78)
    print()
    X_data = load_X_cache()
    ctilde_full = compute_S_and_ctilde(X_data)

    qs = sorted(ctilde_full.keys())
    ctildes = [ctilde_full[q]["ctilde_last"] for q in qs]

    print(f"  Empirical c~_q (last k available):")
    for q, c in zip(qs, ctildes):
        last_k = ctilde_full[q]["ctilde_seq"][-1][0]
        print(f"    q={q:>3}  k_max={last_k}  c~_q ~ {c:.6f}")
    print()

    candidate_1_ord_correlation(qs, ctildes)
    candidate_2_q7_outlier(qs, ctildes)
    candidate_3_step_function(qs, ctildes)
    candidate_4_closed_form(qs, ctildes)
    candidate_5_unifying_formula(qs, ctildes, ctilde_full)


if __name__ == "__main__":
    main()
