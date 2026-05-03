"""
chang_qsd_test.py — Door 3 of Chang connection arc.

Compute survivor-conditioned residue distribution rho_t(r mod 32) for fixed
iterate t across many orbits at N=2^32. Compare to Chang's stationary
pi_chang (already computed exactly). The ratio D(r) = rho/pi is the
trajectory-measure deviation: a Lagarias-class observable.

Question: as t grows, does D(r) stabilize? If so, that limit IS the
quasi-stationary distribution under integer absorption.
"""
import sys
import time
from fractions import Fraction
import numpy as np
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_record_residues(starts, T_MAX=120, max_steps=400000):
    n = len(starts)
    # residues[i, t] = m_t mod 32 (or -1 if absorbed)
    residues = np.full((n, T_MAX), -1, dtype=np.int8)
    for i in prange(n):
        m = starts[i]
        residues[i, 0] = m & 31
        steps = 0
        while m != 1 and steps < max_steps:
            if m > MAX_VAL // 3:
                break
            three_m = 3 * m + 1
            while (three_m & 1) == 0:
                three_m >>= 1
            m = three_m
            steps += 1
            if steps < T_MAX:
                residues[i, steps] = m & 31
    return residues


def chang_pi_mod32():
    """Recompute Chang's stationary projected to mod 32. Returns dict r -> Fraction."""
    odd_residues = list(range(1, 64, 2))
    idx = {r: i for i, r in enumerate(odd_residues)}
    N_LIFTS = 128
    P = [[Fraction(0) for _ in range(32)] for _ in range(32)]
    for i, r in enumerate(odd_residues):
        counts = [0] * 32
        for k in range(N_LIFTS):
            n = r + 64 * k
            mm = 3 * n + 1
            while mm % 2 == 0:
                mm //= 2
            counts[idx[mm % 64]] += 1
        for j in range(32):
            P[i][j] = Fraction(counts[j], N_LIFTS)

    # Solve pi P = pi
    n = 32
    A = [[P[j][i] - (Fraction(1) if i == j else Fraction(0)) for j in range(n)] for i in range(n)]
    A[n - 1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n - 1] = Fraction(1)
    for col in range(n):
        pivot = next(row for row in range(col, n) if A[row][col] != 0)
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            b[col], b[pivot] = b[pivot], b[col]
        piv = A[col][col]
        for j in range(col, n):
            A[col][j] /= piv
        b[col] /= piv
        for row in range(n):
            if row != col and A[row][col] != 0:
                factor = A[row][col]
                for j in range(col, n):
                    A[row][j] -= factor * A[col][j]
                b[row] -= factor * b[col]
    pi = b
    # Project to mod 32
    pi_32 = {}
    for r32 in range(1, 32, 2):
        pi_32[r32] = pi[idx[r32]] + pi[idx[r32 + 32]]
    return pi_32


def main():
    log2N = 32
    N = 1 << log2N
    n_orbits = 2_000_000
    seed = 137
    T_MAX = 120

    print(f"# Door 3: trajectory measure deviation from Chang's pi")
    print(f"# N=2^{log2N}, {n_orbits:,} orbits, T_MAX={T_MAX}")

    # Compute Chang's pi_32
    print("\n# Computing Chang's pi_32 (exact rationals)...")
    pi_chang = chang_pi_mod32()

    # Walk
    print(f"\n# Walking {n_orbits:,} orbits...")
    rng = np.random.default_rng(seed)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_orbits, dtype=np.int64) + 1
    t0 = time.perf_counter()
    residues = walk_record_residues(starts, T_MAX=T_MAX)
    print(f"# Walk: {time.perf_counter() - t0:.1f}s")

    # Compute rho_t(r) and deviation D_t(r) = rho_t / pi_chang
    print(f"\n# Survivor-conditioned residue distribution rho_t(r) and ratio D_t(r) = rho_t / pi_chang")
    print(f"# pi_chang(r) is essentially uniform (deviation max 1.6%, see chang_pi.csv)")
    print(f"# So D_t(r) ~ rho_t(r) * 16, with target 1.0 = matches stationary")
    print()

    test_t = [0, 10, 30, 50, 70, 90, 110]
    odd_r32 = list(range(1, 32, 2))

    print("# Format: D_t(r) for each r mod 32")
    header = f"#   {'t':>3}  {'n_alive':>9}"
    for r in odd_r32:
        header += f"  r={r:>2}"
    print(header)

    rows = []
    for t in test_t:
        if t >= T_MAX:
            continue
        r_t = residues[:, t]
        valid = r_t > 0
        r_alive = r_t[valid]
        n_alive = len(r_alive)
        if n_alive < 1000:
            continue
        # Histogram on odd residues mod 32
        odd_idx = (r_alive & 1).astype(bool)
        r_alive_odd = r_alive[odd_idx]
        counts = np.bincount(r_alive_odd, minlength=32)
        # rho_t(r) = counts[r] / sum
        total = counts[1::2].sum()
        rho = {r: counts[r] / total for r in odd_r32}
        # D_t(r) = rho_t(r) / pi_chang(r)
        D = {r: rho[r] / float(pi_chang[r]) for r in odd_r32}
        line = f"#   {t:>3}  {n_alive:>9,}"
        for r in odd_r32:
            line += f"  {D[r]:>5.3f}"
        print(line)
        rows.append((t, n_alive, rho, D))

    # Stability check: how much does D change between t=70, t=90, t=110?
    print("\n# Stability of deviation D(r) at large t:")
    if len(rows) >= 3:
        t_late = [r for r in rows if r[0] >= 50]
        if len(t_late) >= 2:
            print(f"#   {'r':>3}  {'D(t=50)':>9}  {'D(t=70)':>9}  {'D(t=90)':>9}  {'D(t=110)':>10}  {'flag':>10}")
            t50 = next((r for r in rows if r[0] == 50), None)
            t70 = next((r for r in rows if r[0] == 70), None)
            t90 = next((r for r in rows if r[0] == 90), None)
            t110 = next((r for r in rows if r[0] == 110), None)
            for r in odd_r32:
                d50 = t50[3][r] if t50 else float("nan")
                d70 = t70[3][r] if t70 else float("nan")
                d90 = t90[3][r] if t90 else float("nan")
                d110 = t110[3][r] if t110 else float("nan")
                flag = ""
                if d70 < 0.5:
                    flag = "DEPLETED"
                elif d70 > 1.5:
                    flag = "ENHANCED"
                print(f"#   {r:>3}  {d50:>9.4f}  {d70:>9.4f}  {d90:>9.4f}  {d110:>10.4f}  {flag:>10}")

    # Save CSV
    with open("experiments_output/chang_qsd_test.csv", "w") as f:
        f.write("t,n_alive," + ",".join(f"D_r{r}" for r in odd_r32) + "\n")
        for t, n_alive, rho, D in rows:
            f.write(f"{t},{n_alive}," + ",".join(f"{D[r]:.6f}" for r in odd_r32) + "\n")
    print("\n[save] experiments_output/chang_qsd_test.csv")


if __name__ == "__main__":
    main()
