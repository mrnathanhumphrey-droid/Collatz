"""
delta_diagnostic.py — Leading-vs-subleading diagnostic for eps_n.

Computes:
  (1) |eps_n| * 2^n  (leading-order-removed coefficient)
  (2) delta_n := |eps_n| * 2^n - 1/30  (subleading correction)
  (3) sign(delta_n)
  (4) delta_n / delta_{n-1}  (consecutive ratios, n=3..6)
  (5) delta_n / eps_n  (relative size of correction)

Then fits five ansatze against delta_n on n=2..5, holding n=6 as falsifier:
  (a) delta_n = c * rho^n           (geometric)
  (b) delta_n = c * n^(-alpha) * (1/2)^n  (power-law correction)
  (c) delta_n = c * log(n) * (1/2)^n     (log correction)
  (d) delta_n = c1*rho1^n + c2*rho2^n   (two-term)
  (e) delta_n = c * cos(omega*n + phi) * (1/2)^n  (oscillating)

Pre-registered residual threshold: ansatz "fits" if
    |predicted_delta_6 - actual_delta_6| / |actual_delta_6| < 0.20.

Inputs: experiments_output/result_77_7_eps_exact_through_k7.json
Outputs: prints quantities table + ansatz fits + disposition.
"""

import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


HERE = Path(__file__).parent
EPS_JSON = HERE / "experiments_output" / "result_77_7_eps_exact_through_k7.json"


def load_eps():
    """Load eps_n as exact Fractions from cached JSON."""
    with open(EPS_JSON, "r") as f:
        raw = json.load(f)
    return {int(k): Fraction(int(v["num"]), int(v["den"])) for k, v in raw.items()}


def compute_quantities(eps):
    """Compute the five diagnostic quantities at each n in {2,3,4,5,6}."""
    ONE_30 = Fraction(1, 30)
    out = []
    prev_delta = None
    for n in range(2, 7):
        e = eps[n]
        abs_e = abs(e)
        leading = abs_e * (2 ** n)             # |eps_n| * 2^n  (exact Fraction)
        delta = leading - ONE_30               # delta_n        (exact Fraction)
        ratio = (delta / prev_delta) if prev_delta is not None else None
        rel = delta / e                        # delta_n / eps_n
        out.append({
            "n": n,
            "eps": e,
            "abs_eps_2n": leading,
            "delta": delta,
            "sign_delta": 1 if delta > 0 else (-1 if delta < 0 else 0),
            "delta_ratio": ratio,
            "delta_rel_to_eps": rel,
        })
        prev_delta = delta
    return out


def fmt_frac(x, digits=10):
    return f"{float(x):.{digits}g}"


def print_quantities(q):
    print("=" * 78)
    print("QUANTITIES")
    print("=" * 78)
    hdr = f"{'n':>2} {'|eps_n|*2^n':>16} {'delta_n':>16} {'sign':>5} {'delta_n/delta_{n-1}':>22} {'delta_n/eps_n':>16}"
    print(hdr)
    print("-" * len(hdr))
    for r in q:
        ratio_str = "         -" if r["delta_ratio"] is None else f"{float(r['delta_ratio']):>22.10g}"
        print(f"{r['n']:>2} {float(r['abs_eps_2n']):>16.10g} {float(r['delta']):>16.10g} "
              f"{r['sign_delta']:>5} {ratio_str} {float(r['delta_rel_to_eps']):>16.10g}")


# -------------------- ANSATZ FITS --------------------

def fit_geometric(ns, deltas):
    """delta_n = c * rho^n. Fit on log|delta| vs n (when sign stable)."""
    # Sign of delta_n flips, so a pure geometric c*rho^n with REAL rho cannot match.
    # We'll still fit by least squares (allowing negative rho is meaningless for real
    # geometric in the standard sense; if rho<0 it would alternate every step, not
    # match a single flip at n=5->6). Provide best-fit by LSQ over n=2..5 ignoring
    # sign flip, then predict n=6 honestly.
    ns = np.asarray(ns, dtype=float)
    deltas = np.asarray(deltas, dtype=float)

    def residual(params):
        c, rho = params
        return c * (rho ** ns) - deltas

    # Heuristic initial: ratio of consecutive |delta|
    x0 = [deltas[0] / (1.5 ** ns[0]), 0.7]
    try:
        res = least_squares(residual, x0, method="lm")
        c, rho = res.x
        return {"c": c, "rho": rho}, lambda n: c * (rho ** n)
    except Exception as exc:
        return {"error": str(exc)}, lambda n: float("nan")


def fit_power_law(ns, deltas):
    """delta_n = c * n^(-alpha) * (1/2)^n."""
    ns = np.asarray(ns, dtype=float)
    deltas = np.asarray(deltas, dtype=float)
    # Equivalent: delta_n * 2^n = c * n^(-alpha). So log(|delta*2^n|) = log|c| - alpha*log(n).
    # But delta changes sign — drop the n=6 (held-out) from the fit anyway.
    # We'll fit by nonlinear LSQ on signed delta.
    def residual(params):
        c, alpha = params
        return c * (ns ** (-alpha)) * (0.5 ** ns) - deltas
    x0 = [1.0, 1.0]
    try:
        res = least_squares(residual, x0, method="lm")
        c, alpha = res.x
        return {"c": c, "alpha": alpha}, lambda n: c * (n ** (-alpha)) * (0.5 ** n)
    except Exception as exc:
        return {"error": str(exc)}, lambda n: float("nan")


def fit_log(ns, deltas):
    """delta_n = c * log(n) * (1/2)^n."""
    ns = np.asarray(ns, dtype=float)
    deltas = np.asarray(deltas, dtype=float)
    # Closed-form linear fit: delta_n = c * log(n) * (1/2)^n  --> c = sum(delta * log(n) * (1/2)^n) / sum((log(n)*(1/2)^n)^2)
    g = np.log(ns) * (0.5 ** ns)
    c = float(np.dot(g, deltas) / np.dot(g, g))
    return {"c": c}, lambda n: c * math.log(n) * (0.5 ** n)


def fit_two_term(ns, deltas):
    """delta_n = c1 * rho1^n + c2 * rho2^n.  4 parameters, fit 4 points (n=2..5)."""
    ns = np.asarray(ns, dtype=float)
    deltas = np.asarray(deltas, dtype=float)

    def residual(params):
        c1, rho1, c2, rho2 = params
        return c1 * (rho1 ** ns) + c2 * (rho2 ** ns) - deltas

    # Multiple starting points to escape local minima:
    starts = [
        [0.05, 0.5, -0.05, 0.7],
        [0.1, 0.6, -0.1, 0.4],
        [0.02, 0.5, -0.02, 0.5],   # near-degenerate
        [1.0, 0.5, -0.5, 0.7],
    ]
    best = None
    for x0 in starts:
        try:
            res = least_squares(residual, x0, method="lm", max_nfev=10000)
            cost = res.cost
            if best is None or cost < best[0]:
                best = (cost, res.x)
        except Exception:
            continue
    if best is None:
        return {"error": "no convergence"}, lambda n: float("nan")
    c1, rho1, c2, rho2 = best[1]
    return ({"c1": c1, "rho1": rho1, "c2": c2, "rho2": rho2, "cost": best[0]},
            lambda n: c1 * (rho1 ** n) + c2 * (rho2 ** n))


def fit_oscillating(ns, deltas):
    """delta_n = c * cos(omega*n + phi) * (1/2)^n."""
    ns = np.asarray(ns, dtype=float)
    deltas = np.asarray(deltas, dtype=float)

    def residual(params):
        c, omega, phi = params
        return c * np.cos(omega * ns + phi) * (0.5 ** ns) - deltas

    # Sign pattern: + + + + - (n=2,3,4,5,6).  delta_n positive for n=2..5, flips at n=6.
    # cos(omega*n + phi) needs to be positive across n=2..5 then negative at n=6.
    # Try several seeds:
    starts = [
        [1.0, 0.5, 0.0],
        [1.0, 1.0, -1.0],
        [1.0, math.pi / 4, math.pi / 4],
        [1.0, 0.3, -0.5],
        [-1.0, 1.5, 0.5],
    ]
    best = None
    for x0 in starts:
        try:
            res = least_squares(residual, x0, method="lm", max_nfev=10000)
            cost = res.cost
            if best is None or cost < best[0]:
                best = (cost, res.x)
        except Exception:
            continue
    if best is None:
        return {"error": "no convergence"}, lambda n: float("nan")
    c, omega, phi = best[1]
    return ({"c": c, "omega": omega, "phi": phi, "cost": best[0]},
            lambda n: c * math.cos(omega * n + phi) * (0.5 ** n))


def run_ansatze(q):
    # n=2..5 fit; n=6 held out as falsifier.
    fit_ns = [r["n"] for r in q if r["n"] <= 5]
    fit_deltas = [float(r["delta"]) for r in q if r["n"] <= 5]
    actual_delta_6 = float(q[-1]["delta"])  # n=6 row

    THRESHOLD = 0.20

    results = []
    for name, fitter in [
        ("(a) geometric          delta=c*rho^n",                  fit_geometric),
        ("(b) power-law          delta=c*n^-alpha*(1/2)^n",       fit_power_law),
        ("(c) log                delta=c*log(n)*(1/2)^n",         fit_log),
        ("(d) two-term           delta=c1*rho1^n + c2*rho2^n",    fit_two_term),
        ("(e) oscillating        delta=c*cos(om*n+phi)*(1/2)^n",  fit_oscillating),
    ]:
        params, pred = fitter(fit_ns, fit_deltas)
        d6_pred = pred(6)
        try:
            rel_residual = abs(d6_pred - actual_delta_6) / abs(actual_delta_6)
        except Exception:
            rel_residual = float("nan")
        fits = (rel_residual < THRESHOLD) if (rel_residual == rel_residual) else False
        results.append({
            "name": name,
            "params": params,
            "pred_delta_6": d6_pred,
            "actual_delta_6": actual_delta_6,
            "rel_residual_at_6": rel_residual,
            "fits": fits,
        })
    return results


def print_ansatze(results):
    print()
    print("=" * 78)
    print("ANSATZ FITS  (fit on n=2..5, n=6 held out; threshold = 20% relative residual)")
    print("=" * 78)
    for r in results:
        print(f"\n{r['name']}")
        if isinstance(r["params"], dict) and "error" not in r["params"]:
            for k, v in r["params"].items():
                print(f"    {k:>8} = {v: .8g}")
        else:
            print(f"    {r['params']}")
        print(f"    predicted delta_6 = {r['pred_delta_6']: .6g}")
        print(f"    actual    delta_6 = {r['actual_delta_6']: .6g}")
        print(f"    rel residual      = {r['rel_residual_at_6']: .4g}")
        print(f"    FITS (<20%)       = {r['fits']}")


def main():
    eps = load_eps()
    q = compute_quantities(eps)
    print_quantities(q)
    results = run_ansatze(q)
    print_ansatze(results)

    fitting = [r for r in results if r["fits"]]
    print()
    print("=" * 78)
    print("DISPOSITION")
    print("=" * 78)
    if not fitting:
        print("H_DELTA_IRREGULAR  — no ansatz fits within 20% at the held-out n=6.")
    elif len(fitting) == 1:
        print(f"H_DELTA_* — single ansatz fits: {fitting[0]['name']}")
    else:
        names = [f["name"] for f in fitting]
        print("INCONCLUSIVE — multiple ansatze fit at similar quality:")
        for n in names:
            print(f"  - {n}")


if __name__ == "__main__":
    main()
