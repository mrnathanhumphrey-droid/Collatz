"""
phase_routeB_class_bprior_mc.py — Monte Carlo of real Tao iterations.

Test whether the 9% magnitude residual (Geom(1/2) prediction 0.898 vs empirical 0.984)
closes when v is sampled from REAL Tao iterations instead of Geom(1/2). The corrections
come from Tao's C_A constants — non-iid behavior of (v_1, v_2, ...).

Setup:
- Sample large random odd n (60, 100, 200 bits).
- Apply Tao iterations n → (3n+1)/2^v.
- Track (c, m) where c = n mod 3, m = (accumulated v) mod M.
- Within each trajectory, count (c_prev, m_prev) → (c_new, m_new) transitions.
- Restart trajectory when n < threshold (m resets to 0 for fresh trajectory).
- Build empirical transition matrix.
- Compare spectrum to Geom(1/2) baseline.

If top |λ| → 0.984 with real-Tao: that's the Tao C_A closure of the magnitude residual.
If top |λ| stays ~0.898: C_A corrections don't close it; need other mechanism.
"""
import sys, os, json, time, random
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def tao_step(n):
    x = 3 * n + 1
    v = 0
    while x & 1 == 0:
        x >>= 1
        v += 1
    return x, v


def run_mc(n_start_bits, n_transitions, M, seed=42, warmup=50, restart_threshold=1 << 30):
    random.seed(seed)
    n_states = 2 * M
    T = np.zeros((n_states, n_states), dtype=np.float64)
    state_idx = lambda c, m: (0 if c == 1 else 1) * M + m

    def fresh_start():
        # Pick a random large odd integer; warm up for `warmup` steps
        n = random.getrandbits(n_start_bits) | 1
        for _ in range(warmup):
            if n <= 3 or n % 3 == 0:
                n = random.getrandbits(n_start_bits) | 1
                continue
            n, _ = tao_step(n)
        return n

    n_traj = 0
    n_done = 0
    while n_done < n_transitions:
        n = fresh_start()
        # ensure n coprime to 3
        if n % 3 == 0:
            continue
        c_prev = n % 3
        m_prev = 0
        n_traj += 1
        # Iterate Tao steps until trajectory drops below threshold
        while n > restart_threshold and n_done < n_transitions:
            if n % 3 == 0 or n == 1:
                break
            n_new, v = tao_step(n)
            c_new = n_new % 3
            if c_new == 0:
                # n_new ≡ 0 mod 3: not in (Z/3^?)*, skip (rare; trajectory ends or skips)
                n = n_new
                continue
            m_new = (m_prev + v) % M
            T[state_idx(c_new, m_new), state_idx(c_prev, m_prev)] += 1
            n = n_new
            c_prev = c_new
            m_prev = m_new
            n_done += 1

    # Normalize columns
    col_sums = T.sum(axis=0)
    col_sums_safe = np.where(col_sums > 0, col_sums, 1)
    T_norm = T / col_sums_safe
    return T, T_norm, n_traj


def analyze(T, M, label):
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
    n_transitions = 2_000_000
    out = {}
    for n_start_bits in (60, 200):
        for M in (9, 18, 27, 36):
            label = f"bits={n_start_bits},M={M}"
            print(f"\n=== {label}, n_iters={n_transitions:,} ===")
            t0 = time.time()
            T_count, T_norm, n_traj = run_mc(
                n_start_bits, n_transitions, M,
                seed=42 + n_start_bits + M,
                warmup=20,
                restart_threshold=1 << (n_start_bits // 4),
            )
            elapsed = time.time() - t0
            actual_transitions = int(T_count.sum())
            print(f"  MC: {elapsed:.1f}s, {n_traj} trajectories, {actual_transitions:,} transitions")
            out[label] = {
                "eigenvalues": analyze(T_norm, M, label),
                "n_transitions": actual_transitions,
                "n_trajectories": n_traj,
            }

    out["compare"] = {
        "geom_M18_top_lambda": 0.8976,
        "geom_M27_top_lambda": 0.9501,
        "empirical_PADE_0.984": 0.984,
        "empirical_period_9.2": 9.2,
    }
    with open(os.path.join(OUTDIR, "phase_routeB_class_bprior_mc.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote phase_routeB_class_bprior_mc.json")


if __name__ == "__main__":
    main()
