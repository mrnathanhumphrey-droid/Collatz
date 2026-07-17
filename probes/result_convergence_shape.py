"""
result_convergence_shape.py
===========================
Characterize the convergence shape of S_k -> 7/15 for the qx+1 framework
at q=3, using the exact rational S_k for k=1..5 cached in
result_q_sweep_test_1_envelope.csv.

Steps:
  1. Extract eps_k = S_k - 7/15 as exact Fractions.
  2. Ratio analysis eps_{k+1}/eps_k.
  3. Sign pattern (data-corrected: +, +, -, -, -; ONE sign flip, no oscillation).
  4. Functional-form fits:
       (a) eps_k = A * r^k   (pure geometric)
       (c) eps_k = A * r1^k + B * r2^k   (two-mode)
       (d) eps_k = A * r^k * k^alpha   (algebraic correction)
     Damped-osc form (b) skipped: data shows monotone-from-below convergence
     for k in {3, 4, 5}, no oscillation pattern to fit.
  5. Algebraic-form search on rational structure of eps_1, eps_2, eps_3.
  6. SKIPPED: Hateley-style block-average comparison (no lift data on hand).

Outputs:
  result_convergence_shape.md
  result_convergence_shape_data.csv
  result_convergence_shape_fits.csv
"""
from __future__ import annotations
import csv
import math
import os
import sys
from fractions import Fraction

import numpy as np
from scipy.optimize import curve_fit

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz"
ENV_CSV = r"C:\Collatz\result_q_sweep_test_1_envelope.csv"

S_INF = Fraction(7, 15)


# ---------------- Load exact S_k ----------------

def load_S():
    S = {}
    with open(ENV_CSV) as f:
        for row in csv.DictReader(f):
            if int(row["q"]) != 3:
                continue
            k = int(row["n"])
            S[k] = Fraction(int(row["S_n_num"]), int(row["S_n_den"]))
    return S


# ---------------- Algebraic-form helpers ----------------

def factor_int(n, trial_limit=10**7):
    """Return prime factorization as list of (p, e) for |n| using trial division
    up to trial_limit. If a residual > 1 remains, returns it as ('?', residual)
    in the last slot."""
    n = abs(n)
    if n <= 1:
        return []
    out = []
    p = 2
    while p * p <= n and p <= trial_limit:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            out.append((p, e))
        p += 1
    if n > 1:
        if p * p > n:
            out.append((n, 1))   # n itself is prime
        else:
            out.append(("?", n)) # residual too large to fully factor
    return out


def factor_str(n, trial_limit=10**7):
    if n == 0:
        return "0"
    if abs(n) > 10**30:
        return f"({len(str(abs(n)))}-digit, factorization skipped)"
    sign = "-" if n < 0 else ""
    parts = factor_int(n, trial_limit=trial_limit)
    if not parts:
        return f"{sign}1"
    return sign + " * ".join(
        f"{p}^{e}" if isinstance(p, int) and e > 1
        else (f"residual({e})" if p == "?" else f"{p}")
        for p, e in parts)


# ---------------- Fit candidates ----------------

def fit_pure_geom(ks, y_signed):
    """eps_k = A * r^k  (allow A negative). Fit log |eps_k| = log|A| + k log|r|."""
    log_y = np.log(np.abs(y_signed))
    # OLS fit log|y| = c + m*k
    n = len(ks)
    sx = ks.sum()
    sy = log_y.sum()
    sxx = (ks * ks).sum()
    sxy = (ks * log_y).sum()
    m = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    c = (sy - m * sx) / n
    r = math.exp(m)
    A = math.exp(c) * np.sign(y_signed[-1])  # sign from latest data
    return {"A": A, "r": r, "params": (A, r)}


def model_pure_geom(k, A, r):
    return A * r ** k


def model_two_mode(k, A, B, r1, r2):
    return A * r1 ** k + B * r2 ** k


def model_alg_correction(k, A, r, alpha):
    return A * (r ** k) * (k ** alpha)


def safe_fit(model, ks, y_signed, p0, bounds=None):
    try:
        if bounds is None:
            popt, _ = curve_fit(model, ks, y_signed, p0=p0, maxfev=20000)
        else:
            popt, _ = curve_fit(model, ks, y_signed, p0=p0, bounds=bounds,
                                maxfev=20000)
        y_pred = np.array([model(k, *popt) for k in ks])
        ss_res = float(((y_signed - y_pred) ** 2).sum())
        return {"ok": True, "params": popt, "ss_res": ss_res, "y_pred": y_pred}
    except Exception as e:
        return {"ok": False, "error": str(e)}


def aicc(ss_res, n, p):
    """AICc for least-squares fit assuming Gaussian errors."""
    if ss_res <= 0 or n - p - 1 <= 0:
        return float("inf")
    aic = n * math.log(ss_res / n) + 2 * p
    return aic + 2 * p * (p + 1) / (n - p - 1)


# ---------------- main ----------------

def main():
    S = load_S()
    print("=" * 78)
    print("Convergence shape of S_k -> 7/15 (q=3)")
    print("=" * 78)
    print()
    print("Loaded exact S_k (k=1..5):")
    for k in sorted(S):
        print(f"  S_{k} = {float(S[k]):.10f}")
    print()

    # Step 1: eps_k as Fractions
    eps = {k: S[k] - S_INF for k in S}
    print("Step 1: eps_k = S_k - 7/15 (exact rational)")
    for k in sorted(eps):
        e = eps[k]
        print(f"  eps_{k} = {e.numerator}/{e.denominator}  ~ {float(e):+.10e}")
    print()
    print("  Verification of brief-stated values:")
    print(f"    eps_1 = 1/5? -> {eps[1] == Fraction(1, 5)}")
    print(f"    eps_2 = 1/105? -> {eps[2] == Fraction(1, 105)}")
    print()

    # Step 2: ratios
    print("Step 2: eps_{k+1} / eps_k")
    ratios_exact = {}
    for k in range(1, max(eps)):
        r = eps[k + 1] / eps[k]
        ratios_exact[k] = r
        print(f"  eps_{k+1}/eps_{k} = {r.numerator}/{r.denominator}  "
              f"~ {float(r):+.6f}    (rate-1/2 reference: ~0.5 with sign)")
    print()
    abs_ratios_after_transient = [abs(float(ratios_exact[k]))
                                  for k in [2, 3, 4]]
    print("  |ratio| at k=2,3,4 (after k=1->2 transient):",
          [f"{r:.4f}" for r in abs_ratios_after_transient])
    print(f"  -> approaches 1/2 from above: {abs_ratios_after_transient[0]:.3f}, "
          f"{abs_ratios_after_transient[1]:.3f}, "
          f"{abs_ratios_after_transient[2]:.3f}")
    print()

    # Step 3: sign pattern
    signs = ["+" if eps[k] > 0 else ("-" if eps[k] < 0 else "0")
             for k in sorted(eps)]
    print(f"Step 3: sign pattern (k=1..5): ({', '.join(signs)})")
    print(f"  number of sign flips: "
          f"{sum(1 for i in range(1, len(signs)) if signs[i] != signs[i-1])}")
    print(f"  -> ONE sign flip between k=2 and k=3, monotone-from-below for k>=3")
    print(f"  -> damped-oscillation hypothesis (model b) NOT supported by data")
    print()

    # Step 4: numeric fits
    ks = np.array(sorted(eps), dtype=np.float64)
    y_signed = np.array([float(eps[int(k)]) for k in ks])
    abs_y = np.abs(y_signed)
    log_abs_y = np.log(abs_y)
    n = len(ks)
    fits = {}

    # ENVELOPE fits (operate on log|eps_k|, capture |.|-decay only).
    # These are the right object to test "rate-1/2" claims since the model
    # A * r^k has fixed sign and cannot fit signed data with a sign flip.

    # (a-env) Pure-geometric envelope: log|eps_k| = log|A| + k*log(r)
    #   2 params (log|A|, log(r)). Fit by OLS on (k, log|eps_k|).
    sx = ks.sum(); sy = log_abs_y.sum()
    sxx = (ks * ks).sum(); sxy = (ks * log_abs_y).sum()
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    intercept = (sy - slope * sx) / n
    r_geom = math.exp(slope); A_geom = math.exp(intercept)
    pred_log = intercept + slope * ks
    ss_res_a = float(((log_abs_y - pred_log) ** 2).sum())
    fits["(a-env) pure_geom_envelope"] = {
        "ok": True,
        "params": (A_geom, r_geom),
        "ss_res": ss_res_a,
        "aicc": aicc(ss_res_a, n, 2),
        "display": f"|eps_k| = ({A_geom:.4e}) * ({r_geom:.4f})^k",
    }

    # (a-env-late) Pure-geometric envelope on k=2..5 only (drop transient)
    ks_late = ks[1:]; log_late = log_abs_y[1:]; n_late = len(ks_late)
    sx2 = ks_late.sum(); sy2 = log_late.sum()
    sxx2 = (ks_late * ks_late).sum(); sxy2 = (ks_late * log_late).sum()
    slope2 = (n_late * sxy2 - sx2 * sy2) / (n_late * sxx2 - sx2 * sx2)
    intercept2 = (sy2 - slope2 * sx2) / n_late
    r_late = math.exp(slope2); A_late = math.exp(intercept2)
    pred_late = intercept2 + slope2 * ks_late
    ss_res_a_late = float(((log_late - pred_late) ** 2).sum())
    fits["(a-env-late) pure_geom k>=2"] = {
        "ok": True,
        "params": (A_late, r_late),
        "ss_res": ss_res_a_late,
        "aicc": aicc(ss_res_a_late, n_late, 2),
        "display": (f"|eps_k| = ({A_late:.4e}) * ({r_late:.4f})^k  "
                    f"(fit on k=2..5)"),
    }

    # (d-env) Algebraic correction envelope: log|eps_k| = log|A| + k*log(r) + alpha*log(k)
    #   3 params. Closed-form OLS on (k, log(k)) -> log|y|.
    X = np.column_stack([np.ones(n), ks, np.log(ks)])
    coef, *_ = np.linalg.lstsq(X, log_abs_y, rcond=None)
    pred_log_d = X @ coef
    ss_res_d = float(((log_abs_y - pred_log_d) ** 2).sum())
    A_d = math.exp(coef[0]); r_d = math.exp(coef[1]); alpha_d = coef[2]
    fits["(d-env) alg_correction"] = {
        "ok": True,
        "params": (A_d, r_d, alpha_d),
        "ss_res": ss_res_d,
        "aicc": aicc(ss_res_d, n, 3),
        "display": (f"|eps_k| = ({A_d:.4e}) * ({r_d:.4f})^k * "
                    f"k^({alpha_d:+.4f})"),
    }

    # (c-fixed) Two-mode geometric with rates FIXED at 1/2 and 1/3.
    #   eps_k = A*(1/2)^k + B*(1/3)^k. 2 params (signed fit).
    M = np.column_stack([(0.5) ** ks, (1/3) ** ks])
    coef_c, *_ = np.linalg.lstsq(M, y_signed, rcond=None)
    A_c, B_c = float(coef_c[0]), float(coef_c[1])
    pred_c = M @ coef_c
    ss_res_c = float(((y_signed - pred_c) ** 2).sum())
    fits["(c-fixed) two_mode 1/2 + 1/3"] = {
        "ok": True,
        "params": (A_c, B_c),
        "ss_res": ss_res_c,
        "aicc": aicc(ss_res_c, n, 2),
        "display": (f"eps_k = ({A_c:+.4e})*(1/2)^k + "
                    f"({B_c:+.4e})*(1/3)^k"),
    }

    # (c-free) Two-mode with free rates: 4 params, 5 data => AICc undefined
    # (n - p - 1 = 0). Fit anyway and FLAG as overparameterized.
    f_c_free = safe_fit(model_two_mode, ks, y_signed,
                        p0=(0.5, -0.1, 0.4, 0.3),
                        bounds=([-10, -10, 0.01, 0.01],
                                [10, 10, 0.99, 0.99]))
    if f_c_free["ok"]:
        A2, B2, r1, r2 = f_c_free["params"]
        # AICc undefined here (n-p-1 = 0); use BIC instead, or just report SS
        f_c_free["aicc"] = float("inf")  # flagged
        f_c_free["display"] = (f"eps_k = ({A2:+.4e})*({r1:.4f})^k + "
                               f"({B2:+.4e})*({r2:.4f})^k  "
                               f"[OVERPARAMETERIZED: 4 params, 5 points]")
    fits["(c-free) two_mode free"] = f_c_free

    print("Step 4: functional-form fits (envelope and signed)")
    print()
    print(f"  {'model':>32} | {'AICc':>10} | {'SS_res':>14} | params / formula")
    print("  " + "-" * 84)
    rows_fit = []
    for name, f in fits.items():
        if not f.get("ok"):
            print(f"  {name:>32} | {'FAIL':>10} | {'-':>14} | "
                  f"{f.get('error', 'N/A')}")
            continue
        aicc_str = ("inf" if math.isinf(f["aicc"]) else f"{f['aicc']:.3f}")
        print(f"  {name:>32} | {aicc_str:>10} | {f['ss_res']:>14.4e} | "
              f"{f['display']}")
        rows_fit.append({"model": name,
                         "aicc": (f["aicc"] if not math.isinf(f["aicc"])
                                  else "inf"),
                         "ss_res": f["ss_res"],
                         "formula": f["display"],
                         "params": str(list(f["params"]))})
    print()

    # Best fit by AICc among comparable (envelope) models
    env_models = ["(a-env) pure_geom_envelope", "(d-env) alg_correction"]
    best_env = min(env_models, key=lambda n_: fits[n_]["aicc"])
    print(f"  best envelope fit (full k=1..5): {best_env}")
    print(f"  pure-geom rate (k=1..5): {r_geom:.6f}  (vs reference 0.5)")
    print(f"  pure-geom rate (k=2..5 only): {r_late:.6f}  "
          f"<- expected ~0.5 in asymptotic regime")
    print()

    # Step 5: algebraic form search
    print("Step 5: algebraic structure of eps_k")
    print("  numerator and denominator factorizations:")
    for k in sorted(eps):
        e = eps[k]
        nf = factor_str(e.numerator)
        df = factor_str(e.denominator)
        print(f"    eps_{k} = {e.numerator}/{e.denominator}")
        print(f"            num = {nf}")
        print(f"            den = {df}")
    print()

    # Probe specific patterns
    print("  Pattern probes:")
    # eps_1 = 1/5, eps_2 = 1/(3*5*7).  Note: 5, 3*5*7.
    # Check if eps_k = 1 / (product over odd primes from 3 to 2k+3)?
    # k=1: primes 3..5 -> 3*5=15, not 5. Not match.
    # Check 1/( (2k+3) * (2k+1)! / ?) — likely a dead end.
    # Check eps_k / eps_{k-1} structure:
    print(f"    eps_2 / eps_1 = {ratios_exact[1]} = {float(ratios_exact[1]):.4f} = 1/21 = 1/(3*7)")
    if ratios_exact[1] == Fraction(1, 21):
        print("       -> matches 1/(3*7) exactly")
    # What about eps_3 / eps_2?
    r2 = ratios_exact[2]
    print(f"    eps_3 / eps_2 = {r2.numerator}/{r2.denominator} = {float(r2):.4f}")
    # Factor numerator and denominator of this ratio for clarity
    print(f"           num factor = {factor_str(r2.numerator)}")
    print(f"           den factor = {factor_str(r2.denominator)}")
    r3 = ratios_exact[3]
    print(f"    eps_4 / eps_3 = ratio decimal {float(r3):.6f}")
    print(f"           (denominator factorization is large; reporting decimal only)")
    r4 = ratios_exact[4]
    print(f"    eps_5 / eps_4 = {float(r4):.6f}")
    print()

    # OEIS-like sequence test: small numerator sequence
    # ε_1=1/5, ε_2=1/105: numerators are both 1.
    # ε_3 has num |5191| = 29 * 179.
    print("  Numerator sequence |num(eps_k)|: "
          f"{[abs(eps[k].numerator) for k in sorted(eps)]}")
    print("  Denominator sequence den(eps_k): "
          f"{[eps[k].denominator for k in sorted(eps)]}")
    print()
    print("  No obvious OEIS-style closed form jumps out.")
    print("  - eps_1 = 1/5 = 1/(2*1+3)? -> 2+3=5 ✓ but only 1 term")
    print("  - eps_2 = 1/105 = 1/(3*5*7) = 1/(prod of 3 consecutive odd primes)")
    print("  - eps_3 numerator = -29*179, denominator = 3*5*7^2*19*73 — no clean")
    print("    structural form; suggests eps_k ceases to be a simple rational")
    print("    pattern past k=2.")
    print()

    # Step 6 SKIPPED
    print("Step 6: SKIPPED (no Hateley-lift data on hand)")
    print()

    # Save data CSV
    out_data = os.path.join(OUT_DIR, "result_convergence_shape_data.csv")
    with open(out_data, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "S_k_num", "S_k_den", "S_k_decimal",
                    "eps_k_num", "eps_k_den", "eps_k_decimal", "sign",
                    "abs_ratio_to_prev"])
        for k in sorted(eps):
            e = eps[k]
            sign = "+" if e > 0 else ("-" if e < 0 else "0")
            ratio_prev = ""
            if k > 1:
                rr = float(eps[k] / eps[k-1])
                ratio_prev = f"{rr:+.10f}"
            w.writerow([k, S[k].numerator, S[k].denominator, float(S[k]),
                        e.numerator, e.denominator, float(e), sign,
                        ratio_prev])
    print(f"saved {out_data}")

    # Save fits CSV
    out_fits = os.path.join(OUT_DIR, "result_convergence_shape_fits.csv")
    with open(out_fits, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["model", "aicc", "ss_res",
                                          "formula", "params"])
        w.writeheader()
        for r in rows_fit:
            w.writerow(r)
    print(f"saved {out_fits}")

    # Markdown writeup
    out_md = os.path.join(OUT_DIR, "result_convergence_shape.md")
    md = []
    md.append("# Result: convergence shape of S_k -> 7/15 (q=3)")
    md.append("")
    md.append("**Date:** 2026-05-05.  Characterizes the functional shape "
              "of S_k - 7/15 for k=1..5 using exact rational data cached in "
              "`result_q_sweep_test_1_envelope.csv`.")
    md.append("")
    md.append("## Verdict")
    md.append("")
    md.append(f"- Best envelope fit by AICc (full k=1..5): **{best_env}**.")
    if fits[best_env]["ok"]:
        md.append(f"  - Formula: `{fits[best_env]['display']}`")
        md.append(f"  - AICc = {fits[best_env]['aicc']:.3f}, SS_res = "
                  f"{fits[best_env]['ss_res']:.3e}")
    md.append("- Sign pattern (k=1..5): (+, +, -, -, -). One sign flip "
              "between k=2 and k=3; monotone-from-below for k>=3.")
    md.append("- The damped-oscillation hypothesis (model b in the brief) "
              "is **not supported** by the data — there is no oscillation, "
              "only a single sign flip.")
    md.append("- |eps_{k+1}/eps_k| for k=2,3,4 is "
              f"{abs_ratios_after_transient[0]:.4f}, "
              f"{abs_ratios_after_transient[1]:.4f}, "
              f"{abs_ratios_after_transient[2]:.4f} — approaches 1/2 from "
              "above. Consistent with rate-1/2 envelope being approximately "
              "tight in the limit.")
    md.append("")
    md.append("## eps_k table")
    md.append("")
    md.append("| k | S_k (decimal) | eps_k (rational) | eps_k (decimal) | sign | |eps_{k+1}/eps_k| |")
    md.append("|---|---|---|---|---|---|")
    for k in sorted(eps):
        e = eps[k]
        sign = "+" if e > 0 else "-"
        ratio = ""
        if k < max(eps):
            ratio = f"{abs(float(eps[k+1]/e)):.4f}"
        if e.numerator < 0:
            num_str = f"-{-e.numerator}/{e.denominator}"
        else:
            num_str = f"{e.numerator}/{e.denominator}"
        md.append(f"| {k} | {float(S[k]):.10f} | {num_str} | "
                  f"{float(e):+.6e} | {sign} | {ratio} |")
    md.append("")
    md.append("## Functional-form fits")
    md.append("")
    md.append("| model | AICc | SS_res | formula |")
    md.append("|---|---|---|---|")
    for name, f in fits.items():
        if not f["ok"]:
            md.append(f"| {name} | FAIL | — | {f['error']} |")
        else:
            md.append(f"| {name} | {f['aicc']:.3f} | "
                      f"{f['ss_res']:.4e} | `{f['display']}` |")
    md.append("")
    md.append(f"**Pure-geom envelope rate (k=1..5):** {r_geom:.6f}  "
              f"(reference 1/2 = 0.500000)")
    md.append(f"**Pure-geom envelope rate (k=2..5 only, drops k=1 transient):** "
              f"{r_late:.6f}")
    md.append("")
    md.append("The k=2..5 fit is the asymptotic-regime estimate; the k=1 "
              "transient is a different scale and inflates the rate when "
              "included.")
    md.append("")
    md.append("## Algebraic structure")
    md.append("")
    md.append(f"- eps_1 = 1/5")
    md.append(f"- eps_2 = 1/105 = 1/(3·5·7)")
    md.append(f"- eps_3 = -5191/1019445 = -(29·179)/(3·5·7²·19·73)")
    md.append(f"- eps_2/eps_1 = 1/21 = 1/(3·7) — clean rational factor")
    md.append(f"- eps_3/eps_2 ≈ -0.535 — does not factor cleanly into small primes")
    md.append("")
    md.append("Pattern observation: eps_1 and eps_2 have numerator 1 and "
              "denominators built from small odd primes. eps_3 onward has "
              "compound numerators with no clean structural factor. Suggests "
              "the small-k behavior (k = 1, 2) may have a closed-form origin "
              "while k >= 3 is genuinely algebraic-non-elementary at this "
              "computed precision.")
    md.append("")
    md.append("## What was skipped")
    md.append("")
    md.append("- Step 4(b) damped-osc fit: skipped (no oscillation to fit)")
    md.append("- Step 6 Hateley-style block-average comparison: skipped "
              "(no lift data in workspace)")
    md.append("")
    md.append("## Files")
    md.append("")
    md.append("- `result_convergence_shape.py` — script")
    md.append("- `result_convergence_shape_data.csv` — eps_k as exact rationals + decimals")
    md.append("- `result_convergence_shape_fits.csv` — model AICc / params table")
    md.append("- `result_convergence_shape.md` — this writeup")

    with open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"saved {out_md}")


if __name__ == "__main__":
    main()
