"""
phase_routeB_R1_full_chain_mc.py — Full Markov chain on (Z/3^n)* under real Tao iteration.

Per ROUTEB_R2_NEGATIVE.md disposition: the (c, m mod M) chain is a scalar reduction that
sees only class+b_prior aggregates. The FULL chain on (Z/3^n)* (dim 2·3^{n-1}) may have
additional non-trivial eigenvalues that the scalar reduction misses.

R1 hypothesis: at higher n (4, 5, 6), the full chain's top |λ| approaches empirical 0.984
even though the scalar (c, m) reduction caps out at predicted values.

Method:
- MC trajectories from random 200-bit odd integers.
- Track r = trajectory value mod 3^n.
- Build empirical transition matrix on (Z/3^n)* (dim 2·3^{n-1}).
- Eigendecompose, report top eigenvalues + dominant CC pair.

For n=2..6, dim ∈ {6, 18, 54, 162, 486} — all tractable.
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


def run_full_chain_mc(n_lvl, n_transitions, n_start_bits=200, seed=42, restart_threshold=None):
    """Build empirical transition matrix on (Z/3^n_lvl)*."""
    random.seed(seed)
    mod = 3 ** n_lvl
    coprime = [r for r in range(mod) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(coprime)}
    n_states = len(coprime)
    if restart_threshold is None:
        restart_threshold = 1 << (n_start_bits // 4)

    T = np.zeros((n_states, n_states), dtype=np.float64)

    def fresh_start():
        for _ in range(100):
            n = random.getrandbits(n_start_bits) | 1
            if n % 3 != 0:
                # warm up
                for _ in range(50):
                    if n % 3 == 0 or n <= 3:
                        break
                    n, _ = tao_step(n)
                if n % 3 != 0 and n > restart_threshold:
                    return n
        return None

    n_done = 0
    n_traj = 0
    while n_done < n_transitions:
        n = fresh_start()
        if n is None:
            print("    failed to find fresh start; aborting")
            break
        n_traj += 1
        r_prev = n % mod
        while n > restart_threshold and n_done < n_transitions:
            if n % 3 == 0:
                break
            n_new, v = tao_step(n)
            if n_new % 3 == 0:
                n = n_new
                continue
            r_new = n_new % mod
            if r_prev in idx and r_new in idx:
                T[idx[r_new], idx[r_prev]] += 1
                n_done += 1
            n = n_new
            r_prev = r_new

    col_sums = T.sum(axis=0)
    col_sums_safe = np.where(col_sums > 0, col_sums, 1)
    T_norm = T / col_sums_safe
    return T, T_norm, n_traj


def analyze(T, n_lvl, n_print=15):
    eigvals = np.linalg.eigvals(T)
    eigvals_sorted = sorted(eigvals, key=lambda x: -abs(x))
    print(f"  Top {n_print} eigenvalues at n={n_lvl}, dim={T.shape[0]}:")
    for i, e in enumerate(eigvals_sorted[:n_print]):
        mag = abs(e)
        arg = np.angle(e)
        period = 2 * np.pi / abs(arg) if abs(arg) > 1e-10 else float('inf')
        print(f"    [{i:2d}] |λ|={mag:.6f}, arg={arg:+.4f}rad, period={period:7.3f} steps")
    return [[float(e.real), float(e.imag)] for e in eigvals_sorted]


def main():
    n_transitions = 20_000_000
    out = {}
    for n_lvl in (3, 4, 5, 6, 7, 8):
        print(f"\n=== n_lvl={n_lvl} (dim {2*3**(n_lvl-1)}, mod 3^{n_lvl}={3**n_lvl}) ===")
        t0 = time.time()
        T_count, T_norm, n_traj = run_full_chain_mc(
            n_lvl, n_transitions, n_start_bits=200, seed=42 + n_lvl
        )
        actual_t = int(T_count.sum())
        print(f"  MC: {time.time()-t0:.1f}s, {n_traj} trajectories, {actual_t:,} transitions")
        out[f"n={n_lvl}"] = {
            "eigenvalues": analyze(T_norm, n_lvl),
            "n_transitions": actual_t,
            "n_trajectories": n_traj,
        }
        # Print the Phase 4 closed-form prediction for comparison
        angle = np.pi / 3 ** (n_lvl - 1)
        lam_below_phase4 = 0.5 / abs(1 - 0.5 * np.exp(1j * angle))
        period_phase4 = 2 * np.pi / angle
        print(f"  Phase 4 prediction: |λ_below|={lam_below_phase4:.6f}, period={period_phase4:.2f}")

    out["compare"] = {
        "empirical_PADE_0.984": 0.984,
        "empirical_period_9.2": 9.2,
        "T_lead_43/45": 43/45,
    }
    with open(os.path.join(OUTDIR, "phase_routeB_R1_full_chain_mc.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote phase_routeB_R1_full_chain_mc.json")


if __name__ == "__main__":
    main()
