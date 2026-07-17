"""
result_80b_extended_r.py — push r = 7, 8 to disambiguate δ → 0 vs δ → positive limit.

Uses faster Gauss-sum factorization F̂(3a) = 3·Σ_{v ∈ principal units} e_q(c·v − 3a·log_4 v)
with period 3^r instead of q = 3^{r+1} (3× speedup).

We run r = 3..8 to get 6 data points on log_q|Σ|.
"""
import sys
import math
import cmath

sys.stdout.reconfigure(encoding="utf-8")


def precompute_principal_units(r):
    """Return (principal_units, log_v_table) where principal_units[u] = 4^u mod q."""
    q = 3 ** (r + 1)
    period = 3 ** r
    pu = [0] * period
    pw = 1
    for u in range(period):
        pu[u] = pw
        pw = (pw * 4) % q
    return pu


def F_hat_via_gauss(r, a, pu, c=1):
    """F̂(3a) = 3 · Σ_{u=0}^{period−1} e_q(c·4^u − 3a·u),  faster than direct."""
    q = 3 ** (r + 1)
    period = 3 ** r
    total = complex(0, 0)
    for u in range(period):
        phase = (c * pu[u] - 3 * a * u) % q
        total += cmath.exp(2j * cmath.pi * phase / q)
    return 3 * total


def ind_hat(r, a):
    q = 3 ** (r + 1)
    N = 3 ** (r - 1)
    if a == 0:
        return complex(N, 0)
    z = cmath.exp(2j * cmath.pi * 3 * a / q)
    z_N = cmath.exp(2j * cmath.pi * 3 * a * N / q)
    return (z_N - 1) / (z - 1)


def supp_iter(r):
    period = 3 ** r
    return [a for a in range(period) if a % 3 == 1]


def bilinear_total(r, pu):
    q = 3 ** (r + 1)
    sqrt_q = math.sqrt(q)
    total = complex(0, 0)
    for a in supp_iter(r):
        ind = ind_hat(r, a)
        F = F_hat_via_gauss(r, a, pu)
        psi = F / (3 * sqrt_q)
        total += ind * psi
    return total


def main():
    print("# Extended bilinear-sum scaling at r = 3..8")
    print(f"  {'r':>2}  {'q':>8}  {'|supp|':>7}  {'|Σ|':>14}  {'log_q|Σ|':>10}  {'δ_emp':>8}")
    rows = []
    for r in range(3, 9):
        q = 3 ** (r + 1)
        pu = precompute_principal_units(r)
        S = bilinear_total(r, pu)
        mag = abs(S)
        lq = math.log(mag) / math.log(q) if mag > 0 else -math.inf
        delta = 1.0 - lq
        rows.append((r, q, mag, lq, delta))
        print(f"  {r:>2}  {q:>8}  {3**(r-1):>7}  {mag:>14.4f}  {lq:>10.4f}  {delta:>8.4f}")
        sys.stdout.flush()

    print()
    print("# Trend analysis")
    print("# If δ → 0:  expect δ values to drop monotonically toward 0")
    print("# If δ → δ∞ > 0: expect δ values to settle toward a positive limit")
    print()
    for i in range(1, len(rows)):
        d_prev = rows[i-1][4]
        d_now = rows[i][4]
        diff = d_now - d_prev
        print(f"  r = {rows[i-1][0]} → {rows[i][0]}: δ {d_prev:.4f} → {d_now:.4f}  (Δ = {diff:+.4f})")

    print()
    print("# log-log fit: |Σ| ≈ A · q^θ")
    import numpy as np
    rs = np.array([row[0] for row in rows], dtype=float)
    log_q = np.array([math.log(row[1]) for row in rows])
    log_S = np.array([math.log(row[2]) for row in rows])
    theta, log_A = np.polyfit(log_q, log_S, 1)
    print(f"  θ_fit = {theta:.4f},  A_fit = {math.exp(log_A):.4f}")
    print(f"  ⟹ extrapolated δ_∞ = {1 - theta:.4f}  (NEGATIVE means δ → 0; POSITIVE means closure)")

    # Power-law decay model: δ(r) = δ_∞ + C/r^α
    print()
    print("# Power-law fit: δ(r) = δ_∞ + C/r^α  using the 4 rightmost points (r=5..8)")
    rs_f = np.array([row[0] for row in rows[2:]], dtype=float)
    delta_f = np.array([row[4] for row in rows[2:]])
    # Try several α values
    for alpha in [0.5, 1.0, 1.5, 2.0]:
        A_mat = np.column_stack([np.ones_like(rs_f), rs_f ** (-alpha)])
        sol, *_ = np.linalg.lstsq(A_mat, delta_f, rcond=None)
        delta_inf, C = sol
        residual = np.linalg.norm(A_mat @ sol - delta_f)
        print(f"  α = {alpha}: δ_∞ = {delta_inf:+.4f},  C = {C:+.4f},  residual = {residual:.4f}")


if __name__ == "__main__":
    main()
