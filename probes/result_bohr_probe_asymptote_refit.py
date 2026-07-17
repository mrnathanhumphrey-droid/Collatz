"""
Refit asymptote model using data from k=40..70 (this run) + k=80, 100 from
prior check probe. The k=40..70 window alone is still in rapid growth and
underconstrains A.
"""
import numpy as np
from scipy.optimize import curve_fit

ks = np.array([40, 45, 50, 55, 60, 65, 70, 80, 100], dtype=float)
cpdfs = np.array([50.118, 108.036, 197.417, 320.724, 459.853, 601.643,
                  774.032, 1072.141, 1147.505], dtype=float)


def model(k, A, B, C):
    return A - B * np.exp(-C * k)


print("Fitting chi²/df(k) = A - B·exp(-C·k) on k = 40..100\n")

# Try multiple initial guesses
best = None
for A0 in [1200, 1300, 1500, 2000, 5000]:
    B0 = A0 - cpdfs[0]
    for C0 in [0.03, 0.05, 0.08, 0.10]:
        try:
            popt, pcov = curve_fit(model, ks, cpdfs, p0=[A0, B0, C0],
                                    maxfev=100000)
            A, B, C = popt
            pred = model(ks, *popt)
            ss_res = float(np.sum((cpdfs - pred) ** 2))
            if best is None or ss_res < best[0]:
                best = (ss_res, A, B, C)
        except Exception:
            pass

ss_res, A, B, C = best
print(f"  Best fit: A = {A:.2f}, B = {B:.2f}, C = {C:.5f}")
print(f"  ss_res = {ss_res:.2f}")
print()
print("  k    actual    model     residual")
for k, y in zip(ks, cpdfs):
    pred = model(k, A, B, C)
    print(f"  {int(k):>3}  {y:>9.3f}  {pred:>9.3f}  {y-pred:>+8.3f}")

# Also report per-step rate from the model
print()
print(f"  At k=70: predicted = {model(70, A, B, C):.2f}")
print(f"  At k=100: predicted = {model(100, A, B, C):.2f}")
print(f"  At k=200: predicted = {model(200, A, B, C):.2f}  "
      f"(should be ~A = {A:.2f} if saturating)")
print(f"  At k=500: predicted = {model(500, A, B, C):.2f}")
