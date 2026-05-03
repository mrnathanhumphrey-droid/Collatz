"""
Experiment 38 (TA.3) — parametric fit of gap vs Delta_log.

Use the 40-cell data from experiment 35 (5 observables × 4 modular
resolutions × 2 N-values; gaps essentially k-invariant so only the 10
observable × N pairs matter for the parametric fit).

Fit candidate forms:
  M1: gap = a + b * Delta_log
  M2: gap = a + b * log(Delta_log)
  M3: gap = a + b / Delta_log
  M4: gap = a + b * Delta_log + c * Delta_log^2
  M5: gap = a + b * log(N) + c * log(threshold)  (= a + (b+c) log N - c Delta_log)
  M6: gap = a + b * Delta_log + c * log(N)   (mixed)
  M7: gap = a + b * log(log N)   (only N dependence)
  M8: gap = a + b * Delta_log + c * log(log N)

Report each model's residual SSE and parameters. Identify the lowest-residual
form and interpret.

Data (raw mean) hardcoded from exp 35 output for N = 2^25 and 2^27 at k=8.
"""
import sys

import numpy as np
from scipy.optimize import curve_fit

sys.stdout.reconfigure(encoding="utf-8")


def main():
    # 10 data points: (observable, N, log_N, Delta_log, gap)
    data_raw = [
        # N = 2^25, log N = 16.3287
        ("sigma",          1 << 25, 16.3287, 16.3287, -2.422),
        ("s @ N^(2/3)",    1 << 25, 16.3287,  5.4429, +2.917),
        ("s @ sqrt*log",   1 << 25, 16.3287,  5.3760, +2.962),  # log N / 2 - log log N
        ("s @ sqrt(N)",    1 << 25, 16.3287,  8.1644, +1.938),  # log N / 2
        ("s @ sqrt/log",   1 << 25, 16.3287, 10.9527, +0.355),  # log N / 2 + log log N
        # N = 2^27, log N = 17.7150
        ("sigma",          1 << 27, 17.7150, 17.7150, -2.445),
        ("s @ N^(2/3)",    1 << 27, 17.7150,  5.9050, +3.009),
        ("s @ sqrt*log",   1 << 27, 17.7150,  5.9831, +3.028),
        ("s @ sqrt(N)",    1 << 27, 17.7150,  8.8575, +2.161),
        ("s @ sqrt/log",   1 << 27, 17.7150, 11.7319, +1.152),
    ]
    # Also include 2^28, 2^30, 2^32 for sigma (from TA.1)
    data_raw.extend([
        ("sigma",          1 << 28, 18.4081, 18.4081, -2.4492),
        ("sigma",          1 << 30, 19.7944, 19.7944, -2.4526),
        ("sigma",          1 << 32, 21.1807, 21.1807, -2.4574),
    ])

    log_N = np.array([d[2] for d in data_raw])
    delta_log = np.array([d[3] for d in data_raw])
    gap = np.array([d[4] for d in data_raw])
    log_log_N = np.log(log_N)

    n = len(gap)
    print(f"Data: {n} (observable, N) cells")
    print(f"  Delta_log range: [{delta_log.min():.3f}, {delta_log.max():.3f}]")
    print(f"  gap range:       [{gap.min():.3f}, {gap.max():.3f}]")
    print()

    def fit_and_score(name, fn, p0, x_data):
        try:
            popt, _ = curve_fit(fn, x_data, gap, p0=p0)
            pred = fn(x_data, *popt)
            resid = gap - pred
            sse = (resid ** 2).sum()
            r2 = 1.0 - sse / ((gap - gap.mean()) ** 2).sum()
            params_str = ", ".join(f"{p:>+.4f}" for p in popt)
            print(f"  {name:>40}    SSE={sse:>9.4f}    R^2={r2:>7.4f}    "
                  f"params=[{params_str}]")
            return name, popt, sse, r2, resid
        except Exception as e:
            print(f"  {name:>40}    fit failed: {e}")
            return name, None, None, None, None

    print("=== Parametric fits ===")
    results = []

    # M1: gap = a + b * Delta_log
    results.append(fit_and_score(
        "gap = a + b*Δlog",
        lambda x, a, b: a + b * x, [0, 0], delta_log
    ))

    # M2: gap = a + b * log(Delta_log)
    results.append(fit_and_score(
        "gap = a + b*log(Δlog)",
        lambda x, a, b: a + b * np.log(x), [0, 0], delta_log
    ))

    # M3: gap = a + b / Delta_log
    results.append(fit_and_score(
        "gap = a + b/Δlog",
        lambda x, a, b: a + b / x, [0, 0], delta_log
    ))

    # M4: gap = a + b*Δlog + c*Δlog^2
    results.append(fit_and_score(
        "gap = a + b*Δlog + c*Δlog^2",
        lambda x, a, b, c: a + b * x + c * x ** 2, [0, 0, 0], delta_log
    ))

    # M5: gap = a + (b+c) log N - c Delta_log = (a, log_N, Delta_log, slopes)
    def m5(X, a, b_log_N, c_delta):
        log_N_x, delta_x = X
        return a + b_log_N * log_N_x + c_delta * delta_x
    results.append(fit_and_score(
        "gap = a + b*log(N) + c*Δlog",
        m5, [0, 0, 0], (log_N, delta_log)
    ))

    # M6: gap = a + b*Δlog + c*log(N)
    results.append(fit_and_score(
        "gap = a + b*Δlog + c*log(N)",
        lambda X, a, b, c: a + b * X[0] + c * X[1], [0, 0, 0], (delta_log, log_N)
    ))

    # M7: gap = a + b*log log N
    results.append(fit_and_score(
        "gap = a + b*log(log N)",
        lambda x, a, b: a + b * x, [0, 0], log_log_N
    ))

    # M8: gap = a + b*Δlog + c*log(log N)
    results.append(fit_and_score(
        "gap = a + b*Δlog + c*log(log N)",
        lambda X, a, b, c: a + b * X[0] + c * X[1], [0, 0, 0], (delta_log, log_log_N)
    ))

    # M9: gap = a + b*Δlog + c*log(Δlog)
    results.append(fit_and_score(
        "gap = a + b*Δlog + c*log(Δlog)",
        lambda x, a, b, c: a + b * x + c * np.log(x), [0, 0, 0], delta_log
    ))

    # M10: gap = a + b * (log N - Δlog) = a + b * log(threshold)
    log_threshold = log_N - delta_log
    results.append(fit_and_score(
        "gap = a + b*log(threshold)",
        lambda x, a, b: a + b * x, [0, 0], log_threshold
    ))

    # M11: gap = a + b/Δlog + c*Δlog
    results.append(fit_and_score(
        "gap = a + b/Δlog + c*Δlog",
        lambda x, a, b, c: a + b / x + c * x, [0, 0, 0], delta_log
    ))

    # Best model
    best = min((r for r in results if r[2] is not None), key=lambda r: r[2])
    print()
    print(f"Best (lowest SSE): {best[0]}    SSE={best[2]:.4f}    R^2={best[3]:.4f}")

    # Print residuals for best model
    print()
    print(f"Per-cell residuals under best model:")
    print(f"  {'observable':>14} {'log N':>9} {'Δlog':>9} {'gap':>9} {'pred':>9} {'resid':>+9}")
    pred_best_resid = best[4]
    # Recompute pred from gap - resid
    for d, r in zip(data_raw, pred_best_resid):
        obs, N, lN, dl, g = d
        pred = g - r
        print(f"  {obs:>14} {lN:>9.4f} {dl:>9.4f} {g:>+9.4f} {pred:>+9.4f} {r:>+9.4f}")


if __name__ == "__main__":
    main()
