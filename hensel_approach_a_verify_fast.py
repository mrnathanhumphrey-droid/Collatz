"""
hensel_approach_a_verify_fast.py — FAST version of Hensel Approach A verification.

Replaces the mpmath direct-sum verification (which OOM'd at 4 GB) with vectorized
numpy float64. The Hensel claim is exact at the algebraic level; we only need
~1e-10 precision to confirm and float64 gives ~1e-15 floor, well sufficient.

CLAIM (Approach A): for prime p >= 3, r >= 2, c=1:
    G_p(a) = p^((r+1)/2) * eta_p(r) * e_q(P_a(s*(r)))
where:
    s*(r) = (C_a - 1) // p  mod p^{r-1}
    P_a(s*(r)) = p*c*s*(r) - C_a * L_p(1+p*s*(r))  (mod p^{r+1})
    eta_p(r) = (1/sqrt(p)) * Sum_{h=0}^{p-1} e_p(h^2 / 2)    at even r
             = 1                                              at odd r

G_p(a) "actual" = Sum_{s=0}^{p^r - 1} e_q(p*c*s - C_a*L_p(1+p*s))
                  (vectorized via numpy)

Tests:
    (p, r) in {(3, 4), (3, 5), (3, 6), (5, 4), (5, 5), (7, 4), (7, 5), (7, 6),
               (11, 4), (11, 5)}.

Pass criteria (pre-reg):
    < ~1e-10 max rel dev across all 12 cells: structure CONFIRMED (H_HENSEL_CLOSES)
    < ~1e-3 (good approx, not exact):         needs cross-term audit
    >  1% miss at any cell:                   structure wrong, identify cell
    Match at r=4, degrades at r=6:            r=4-specific, not family-level
"""

import sys
import os
import math
import csv
from fractions import Fraction
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

OUTDIR = r"C:\Collatz"
CSV_PATH = os.path.join(OUTDIR, "HENSEL_APPROACH_A_VERIFICATION.csv")


def J_for_p(p: int, m: int) -> int:
    """J = largest j such that for all k <= j, k - v_p(k) < m."""
    j = 1
    while True:
        x = j + 1
        v = 0
        while x % p == 0:
            x //= p
            v += 1
        if (j + 1) - v >= m:
            return j
        j += 1


def truncated_padic_log(p: int, s: int, J: int) -> Fraction:
    """L_p(1+p*s) = Sum_{j=1}^J (-1)^{j-1}/j * (p*s)^j  as exact Fraction."""
    if s == 0:
        return Fraction(0)
    L = Fraction(0)
    ps = Fraction(p * s)
    ps_pow = Fraction(1)
    for j in range(1, J + 1):
        ps_pow *= ps
        L += Fraction((-1) ** (j - 1), j) * ps_pow
    return L


def L_mod_q(p: int, s: int, J: int, q: int) -> int:
    """Reduce L_p(1+p*s) mod q. Returns integer in [0, q)."""
    if s == 0:
        return 0
    L_frac = truncated_padic_log(p, s, J)
    num = L_frac.numerator
    den = L_frac.denominator
    # Strip p-factors from den (Fraction is already in lowest terms but den may have p-content)
    v = 0
    d = den
    while d % p == 0:
        d //= p
        v += 1
    if v > 0:
        # num must be divisible by p^v for L_p(1+ps) to be a p-adic integer
        if num % (p ** v) != 0:
            raise ValueError(
                f"L_p(1+{p}*{s}) not p-adic integer: num={num}, den={den}, v={v}"
            )
        num //= p ** v
    # Now den is coprime to p, hence coprime to q = p^{r+1}.
    return (num * pow(d, -1, q)) % q


def precompute_L_vals(p: int, r: int) -> np.ndarray:
    """L_vals[s] = L_p(1+p*s) mod q  for s = 0 .. p^r - 1.  int64 array."""
    q = p ** (r + 1)
    period = p ** r
    J = J_for_p(p, r + 1)
    L_vals = np.zeros(period, dtype=np.int64)
    for s in range(period):
        L_vals[s] = L_mod_q(p, s, J, q)
    return L_vals


def G_actual_vectorized(p: int, r: int, C_a: int, L_vals: np.ndarray) -> complex:
    """Compute G_p(a) = Sum_s e_q(p*s - C_a * L_vals[s]) using numpy vectorization."""
    q = p ** (r + 1)
    period = p ** r
    s_arr = np.arange(period, dtype=np.int64)
    # Use Python int arithmetic for the modular reduction to avoid int64 overflow
    # when p*s - C_a * L_vals[s] could exceed 2^63 at large r. We can compute mod q
    # via numpy if we're careful with intermediate types.
    # Safe approach: keep as Python int via np.array dtype=object for the mod step.
    # Actually for q up to ~10^7 (p=11, r=5 gives q = 11^6 ~ 1.77M; p=11,r=6 gives 19.5M),
    # int64 is fine: p*s up to p^{r+1} ~ q; C_a * L_vals up to (p^r)*(q-1) ~ p^{2r+1}.
    # At p=11, r=5: p^{2r+1} = 11^11 = 2.85e11, fits in int64 (max ~9.2e18). Good.
    # At p=11, r=6: 11^13 = 3.45e13, fits. Good.
    phases_int = (p * s_arr - C_a * L_vals) % q
    angles = 2.0 * np.pi * phases_int.astype(np.float64) / q
    z = np.exp(1j * angles)
    return complex(np.sum(z))


def predict_G_p_a(p: int, r: int, C_a: int, L_vals: np.ndarray) -> complex:
    """Approach A prediction:  G_p(a) = p^{(r+1)/2} * eta_p(r) * e_q(P_a(s*(r)))."""
    q = p ** (r + 1)
    p_rm1 = p ** (r - 1)
    # Hensel-lifted saddle, exact: s*(r) = (C_a - 1) // p mod p^{r-1}
    C_a_canonical = C_a % (p ** r)
    s_star = ((C_a_canonical - 1) // p) % p_rm1

    # P_a(s*(r)) = p*c*s_star - C_a * L_p(1+p*s_star) mod q
    L_val_sstar = int(L_vals[s_star])  # we have it pre-computed
    phase_int = (p * s_star - C_a * L_val_sstar) % q

    mag = float(p) ** ((r + 1) / 2.0)
    phase_factor = complex(math.cos(2.0 * math.pi * phase_int / q),
                           math.sin(2.0 * math.pi * phase_int / q))

    if r % 2 == 0:
        # eta_p = (1/sqrt(p)) sum_h e_p(coeff * h^2) with coeff = 1/2 mod p
        coeff_int = pow(2, -1, p)
        h_arr = np.arange(p, dtype=np.int64)
        ph = (coeff_int * h_arr * h_arr) % p
        s_sum = complex(np.sum(np.exp(2j * np.pi * ph.astype(np.float64) / p)))
        eta = s_sum / math.sqrt(p)
    else:
        eta = 1.0 + 0.0j

    return mag * eta * phase_factor


def verify_cell(p: int, r: int) -> dict:
    """Run the verification at (p, r). Returns summary dict."""
    q = p ** (r + 1)
    period = p ** r
    J = J_for_p(p, r + 1)

    print(f"\n## p={p}, r={r}, q={q}, period={period}, J_p={J}")
    t0 = time.time()

    # Precompute L_vals
    L_vals = precompute_L_vals(p, r)
    t_L = time.time() - t0

    # L_tilde from L_vals[1]
    L1_mod = int(L_vals[1])
    if L1_mod % p != 0:
        return {"p": p, "r": r, "error": f"L_p(1+p) not divisible by p: {L1_mod}"}
    L_tilde = (L1_mod // p) % (p ** r)
    L_tilde_inv = pow(L_tilde, -1, p ** r)

    # Support: {a in Z/p^r : a ≡ 1 mod p}
    support = list(range(1, p ** r, p))
    supp_size = len(support)
    print(f"  L_tilde={L_tilde}, L_tilde_inv={L_tilde_inv}, |support|={supp_size}")
    print(f"  precompute L_vals: {t_L:.2f}s")

    t1 = time.time()
    max_rel_dev = 0.0
    max_abs_dev = 0.0
    worst_a = None
    sample_rows = []  # capture a few rows for the CSV
    for idx, a in enumerate(support):
        C_a = (a * L_tilde_inv) % (p ** r)
        G_actual = G_actual_vectorized(p, r, C_a, L_vals)
        G_pred = predict_G_p_a(p, r, C_a, L_vals)
        abs_dev = abs(G_actual - G_pred)
        abs_actual = abs(G_actual)
        rel_dev = abs_dev / abs_actual if abs_actual > 1e-300 else abs_dev
        if rel_dev > max_rel_dev:
            max_rel_dev = rel_dev
            worst_a = a
        max_abs_dev = max(max_abs_dev, abs_dev)
        if idx < 3 or idx == supp_size - 1:
            sample_rows.append((a, C_a, abs(G_actual), abs(G_pred), rel_dev))

    t_verify = time.time() - t1

    # Diagnostic: pass categories per pre-reg
    if max_rel_dev < 1e-10:
        status = "PASS_1e-10"
    elif max_rel_dev < 1e-3:
        status = "APPROX_1e-3"
    elif max_rel_dev < 1e-2:
        status = "WEAK_APPROX_1e-2"
    else:
        status = "FAIL"

    print(f"  verify all a: {t_verify:.2f}s")
    print(f"  max_abs_dev = {max_abs_dev:.6e}")
    print(f"  max_rel_dev = {max_rel_dev:.6e}  ({status})")
    print(f"  worst_a = {worst_a}")
    for (a, C_a, mA, mP, rd) in sample_rows:
        print(f"    a={a:>6} C_a={C_a:>6} |G_actual|={mA:.6f} |G_pred|={mP:.6f} rel_dev={rd:.2e}")

    return {
        "p": p, "r": r, "q": q, "period": period, "J": J,
        "supp_size": supp_size,
        "max_abs_dev": max_abs_dev,
        "max_rel_dev": max_rel_dev,
        "status": status,
        "worst_a": worst_a,
        "elapsed_s": t_verify + t_L,
    }


def main():
    cells = [
        (3, 4), (3, 5), (3, 6),
        (5, 4), (5, 5),
        (7, 4), (7, 5), (7, 6),
        (11, 4), (11, 5),
    ]
    print("# HENSEL Approach A — fast numerical verification (float64 + numpy)")
    print(f"# {len(cells)} cells")
    print()

    rows_csv = [["p", "r", "q", "period", "J_p", "supp_size",
                 "max_abs_dev", "max_rel_dev", "status", "worst_a", "elapsed_s"]]
    overall_t0 = time.time()
    summaries = []
    for (p, r) in cells:
        try:
            res = verify_cell(p, r)
            summaries.append(res)
            if "error" in res:
                rows_csv.append([p, r, "", "", "", "", "ERROR", res["error"], "", "", ""])
            else:
                rows_csv.append([
                    res["p"], res["r"], res["q"], res["period"], res["J"], res["supp_size"],
                    f"{res['max_abs_dev']:.6e}", f"{res['max_rel_dev']:.6e}",
                    res["status"], res["worst_a"], f"{res['elapsed_s']:.2f}",
                ])
        except Exception as e:
            print(f"  EXCEPTION at p={p}, r={r}: {e}")
            rows_csv.append([p, r, "", "", "", "", "EXCEPTION", str(e), "", "", ""])

    overall_t = time.time() - overall_t0
    print(f"\n# Total elapsed: {overall_t:.1f}s")

    with open(CSV_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerows(rows_csv)
    print(f"# CSV written to {CSV_PATH}")

    print("\n# Summary by cell:")
    print(f"  {'p':>3} {'r':>2} {'max_rel_dev':>14} {'status':>16} {'time(s)':>8}")
    for res in summaries:
        if "error" in res:
            print(f"  {res['p']:>3} {res['r']:>2} {'ERROR':>14} {'':>16} {'':>8}")
        else:
            print(f"  {res['p']:>3} {res['r']:>2} {res['max_rel_dev']:>14.6e} "
                  f"{res['status']:>16} {res['elapsed_s']:>8.2f}")

    # Final disposition
    statuses = [r.get("status", "ERROR") for r in summaries]
    all_pass = all(s == "PASS_1e-10" for s in statuses)
    any_fail = any(s in ("FAIL", "WEAK_APPROX_1e-2") for s in statuses)
    if all_pass:
        print("\n# DISPOSITION: APPROACH_A_VALIDATED (all cells pass at 1e-10)")
    elif any_fail:
        print("\n# DISPOSITION: APPROACH_A_INVALID (at least one cell fails — see worst_a)")
    else:
        print("\n# DISPOSITION: APPROACH_A_APPROXIMATE (all cells in approximate range)")


if __name__ == "__main__":
    main()
