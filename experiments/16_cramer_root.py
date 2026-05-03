"""
Experiment 16 — Cramer root (exact gambler's ruin) prediction for qx+1 slopes.

Trajectory step distribution is Geom(1/2) on v_2 (verified empirically in 15b).
Step in log-coords: X = log(q) - v*log(2).
MGF: E[exp(-theta*X)] = q^(-theta) * sum_{v=1}^inf 0.5^v * 2^(theta*v)
                     = q^(-theta) * (2^theta/2) / (1 - 2^theta/2)   for 2^theta < 2

Setting MGF = 1:
    q^(-theta) = (2 - 2^theta) / 2^theta = 2^(1-theta) - 1

Cramer: P(reach 0 from V) = exp(-theta * V)

Starting log-value at j-group: V ~ j*log(q) + const.
Slope of log(conv_rate) vs j: -theta * log(q)

Predicted slopes:
  q=5: theta ~= 0.358, slope ~= -0.576
  q=7: theta ~= 0.720, slope ~= -1.401
  q=9: theta ~= 0.741, slope ~= -1.628
  q=11: theta ~= 0.805, slope ~= -1.930

Compare to empirical slopes from experiment 12.

Usage:
    python 16_cramer_root.py
"""
import numpy as np
from scipy.optimize import brentq

EMPIRICAL_SLOPES = {
    5: -0.5619,   # N=10^8, k=8, 32785 orbits
    7: -1.3685,   # N=10^8, k=6, 258 orbits
    9: -2.0529,   # N=10^8, k=6, 104 orbits
    11: -1.6458,  # N=10^9, k=6, 36 orbits
}


def cramer_root(q):
    """Solve q^(-theta) = 2^(1-theta) - 1 for theta in (0, 1)."""
    def f(theta):
        return q ** (-theta) - (2 ** (1 - theta) - 1)
    # f(0) = 1 - (2 - 1) = 0 (trivial)
    # f(1) = 1/q - (1 - 1) = 1/q > 0
    # We want the non-trivial positive root. f decreases then increases through 0.
    # Search in (0.001, 0.99)
    try:
        return brentq(f, 0.001, 0.99)
    except ValueError:
        return float("nan")


print(f"=== Cramer-root prediction for qx+1 convergence slope (Geom(1/2) v_2) ===")
print(f"{'q':>3} {'theta_Cramer':>14} {'log(q)':>10} {'slope_pred':>12} {'slope_emp':>12} {'ratio_emp/pred':>16}")
for q in [3, 5, 7, 9, 11, 13]:
    theta = cramer_root(q)
    log_q = np.log(q)
    slope_pred = -theta * log_q
    slope_emp = EMPIRICAL_SLOPES.get(q, None)
    if slope_emp is not None:
        ratio = slope_emp / slope_pred
        print(f"{q:>3} {theta:>14.5f} {log_q:>10.4f} {slope_pred:>12.4f} {slope_emp:>12.4f} {ratio:>16.4f}")
    else:
        print(f"{q:>3} {theta:>14.5f} {log_q:>10.4f} {slope_pred:>12.4f} {'--':>12} {'--':>16}")

print()
print(f"=== Coincidental near-constancy of slope/log(q/4) at q in {{5,7,9}} ===")
print(f"{'q':>3} {'slope/log(q/4)':>16} {'theta*log(q)/log(q/4)':>22}")
for q in [5, 7, 9, 11]:
    theta = cramer_root(q)
    pred = theta * np.log(q) / np.log(q / 4.0)
    emp_slope = EMPIRICAL_SLOPES.get(q)
    if emp_slope is not None:
        emp_C = -emp_slope / np.log(q / 4.0)
        print(f"{q:>3} {emp_C:>16.4f} {pred:>22.4f}")
print()
print(f"For q=3, log(q/4)<0 so the formula doesn't apply (orbits mean-reverting).")
print(f"Cramer-pred 'C' (= theta*log(q)/log(q/4)) varies slowly with q in [5,9],")
print(f"crossing 2.5 around q=6-7. The empirical ~2.5 cluster at q=5,7,9 is coincidental")
print(f"in this q range, not a universal constant.")
