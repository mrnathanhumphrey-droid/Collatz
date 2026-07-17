"""
phase_routeB_class_bprior_geom.py — Route B first probe.

Build the 18-state Markov transition matrix on (class, b_prior mod 9) under
the Geom(1/2) v-distribution (basic DWM trajectory model). Compute spectrum.

State: (c ∈ {+, -}, m ∈ Z/9), index = c*9 + m.
Class evolution (Syracuse): starting from c ∈ {+, -}, ending class depends ONLY on v parity:
  v even → class+
  v odd  → class−
b_prior: m → (m + v) mod 9.

So under v ~ Geom(1/2), the transition probability:
  P((c, m) → (c_v, m+v mod 9)) = 2^{-v}   for v ≥ 1, c_v = +/- by v parity.

Test: does this 18-state chain have a CC pair near (0.984, e^{±2πi/9.2})? If yes, it's
the structural source of empirical PADE period-9 CC pair. If no, route B requires
real-Tao v distribution (NOT Geom(1/2)) or a different operator structure.

Also do M = 18, 27 for resolution.
"""
import sys, os, json
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def build_T(M, V_MAX=64):
    n_states = 2 * M
    T = np.zeros((n_states, n_states), dtype=float)
    for c_in in (0, 1):
        for m_in in range(M):
            idx_in = c_in * M + m_in
            for v in range(1, V_MAX + 1):
                m_out = (m_in + v) % M
                c_out = 0 if v % 2 == 0 else 1
                idx_out = c_out * M + m_out
                T[idx_out, idx_in] += 2.0 ** (-v)
    return T


def analytic_eigvals(M, k_max=None):
    """λ_k = 0.5·e^{2πik/M} / (1 − 0.5·e^{2πik/M}), k = 0..M-1. Plus 0 (degenerate from class rank-1)."""
    if k_max is None:
        k_max = M
    return [0.5 * np.exp(2j * np.pi * k / M) / (1 - 0.5 * np.exp(2j * np.pi * k / M)) for k in range(k_max)]


def analyze(T, M, label):
    col_sums = T.sum(axis=0)
    print(f"  col sums: min={col_sums.min():.6f}, max={col_sums.max():.6f}, mean={col_sums.mean():.6f}")
    eigvals = np.linalg.eigvals(T)
    eigvals_sorted = sorted(eigvals, key=lambda x: -abs(x))
    print(f"  Top 12 eigenvalues of T_{label} (M={M}):")
    for i, e in enumerate(eigvals_sorted[:12]):
        mag = abs(e)
        arg = np.angle(e)
        period = 2 * np.pi / abs(arg) if abs(arg) > 1e-10 else float('inf')
        print(f"    [{i:2d}] |λ|={mag:.6f}, arg={arg:+.4f}rad, period={period:7.3f} steps")
    return [[float(e.real), float(e.imag)] for e in eigvals_sorted]


def main():
    out = {}
    for M in (9, 18, 27):
        print(f"\n=== M = {M} (state space dim {2*M}) ===")
        T = build_T(M, V_MAX=64)
        out[f"M={M}"] = {"eigenvalues_numerical": analyze(T, M, "geom")}

        print(f"  Analytic prediction (closed form): λ_k = 0.5·e^{{2πik/{M}}} / (1 − 0.5·e^{{2πik/{M}}})")
        analytic = analytic_eigvals(M)
        analytic_sorted = sorted(analytic, key=lambda x: -abs(x))
        for i, e in enumerate(analytic_sorted[:6]):
            mag = abs(e)
            arg = np.angle(e)
            period = 2 * np.pi / abs(arg) if abs(arg) > 1e-10 else float('inf')
            print(f"    [{i:2d}] |λ|={mag:.6f}, arg={arg:+.4f}rad, period={period:7.3f} steps")
        out[f"M={M}"]["eigenvalues_analytic_k_top6"] = [[float(e.real), float(e.imag)] for e in analytic_sorted[:6]]

    out["compare"] = {"43/45": 43/45, "empirical_0.984": 0.984, "period_9.2": 9.2}
    with open(os.path.join(OUTDIR, "phase_routeB_class_bprior_geom.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote phase_routeB_class_bprior_geom.json")


if __name__ == "__main__":
    main()
