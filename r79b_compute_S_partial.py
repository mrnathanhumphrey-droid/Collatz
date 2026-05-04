"""
Result 79-extension (R79b): empirical computation of |S_partial(r)| for r = 8..max,
with explicit cubic-phase verification and δ-saving fit.

S_partial(r, c, m) := Σ_{u=0}^{N-1} e_q(c·4^u − 9m·u),  N = 3^{r-1}, q = 3^{r+1}

Per Theorem 78.4 / 78.6 of v3.7 framework, S_partial relates by Plancherel to the dual:
    S_partial = (3·e_q(1)/√q) · Σ_{a ∈ supp} ψ(a) · D̄(3a)
where ψ(a) = e_q(P_a(s*(C_a))) (saddle-point closed form, exact at r=3) and
D(3a) = Σ_{u=0}^{N-1} e_q(3au) is the Dirichlet kernel.

We compute BOTH:
  K(r, c, m) directly via Σ_u e_q(c·4^u − 9mu) — O(N) per (c, m), Numba @njit parallel
  Saddle-prediction phase ψ_lead(a) for cross-verification at small r
and report:
  ρ(r) := log|K(r)| / log(q^{1/2}) − 1
  δ_emp := −ρ(r)/2  (empirical saving exponent)
"""
import sys
import os
import csv
import time
import math
import argparse
from fractions import Fraction
import numpy as np
from numba import njit, prange, types
from numba.typed import Dict

sys.stdout.reconfigure(encoding="utf-8")

# ----------------------------------------------------------------------
# Direct Kalafatelis sum (Numba-accelerated)
# ----------------------------------------------------------------------

@njit(cache=True, parallel=True, fastmath=False)
def kalafatelis_sum_direct(r: int, c: int, m: int) -> complex:
    """K(r, c, m) = Σ_{u=0}^{N-1} e_q(c·4^u − 9m·u), via parallel block reduction.

    Implementation: split u into blocks, accumulate each block's complex sum,
    then add. Within a block, run sequentially (4^u recurrence is sequential).
    Block partition: each thread starts from a precomputed 4^{u_block_start}.
    """
    q = 3 ** (r + 1)
    N = 3 ** (r - 1)
    nine_m = 9 * m
    inv_q = 2.0 * math.pi / q

    # Block size: split work across threads
    n_blocks = 32
    block_size = (N + n_blocks - 1) // n_blocks

    # Precompute block-start values 4^{block * block_size} mod q
    starts = np.empty(n_blocks, dtype=np.int64)
    cur = 1
    for b in range(n_blocks):
        starts[b] = cur
        # Advance by block_size
        for _ in range(block_size):
            cur = (cur * 4) % q

    block_sums_re = np.zeros(n_blocks, dtype=np.float64)
    block_sums_im = np.zeros(n_blocks, dtype=np.float64)

    for b in prange(n_blocks):
        u_start = b * block_size
        u_end = min(u_start + block_size, N)
        if u_start >= N:
            continue
        x = starts[b]
        sre = 0.0
        sim = 0.0
        for u in range(u_start, u_end):
            phase = (c * x - nine_m * u) % q
            angle = inv_q * phase
            sre += math.cos(angle)
            sim += math.sin(angle)
            x = (x * 4) % q
        block_sums_re[b] = sre
        block_sums_im[b] = sim

    return complex(block_sums_re.sum(), block_sums_im.sum())


@njit(cache=True)
def kalafatelis_sum_serial(r: int, c: int, m: int) -> complex:
    """Reference single-threaded version — slower but obviously correct."""
    q = 3 ** (r + 1)
    N = 3 ** (r - 1)
    nine_m = 9 * m
    inv_q = 2.0 * math.pi / q
    sre = 0.0
    sim = 0.0
    x = 1
    for u in range(N):
        phase = (c * x - nine_m * u) % q
        angle = inv_q * phase
        sre += math.cos(angle)
        sim += math.sin(angle)
        x = (x * 4) % q
    return complex(sre, sim)


# ----------------------------------------------------------------------
# Saddle-point cubic phase (closed form per Theorem 78.6, exact at r=3)
# ----------------------------------------------------------------------

def J_for_p3(m: int) -> int:
    """Truncation J such that (3s)^j / j has v_3 ≥ m for all j > J."""
    j = 1
    while True:
        x = j + 1
        v = 0
        while x % 3 == 0:
            x //= 3
            v += 1
        if (j + 1) - v >= m:
            return j
        j += 1


def truncated_3adic_log_mod_q(s: int, J: int, q: int) -> int:
    """L(1+3s) = Σ_{j=1}^J (-1)^{j-1}/j · (3s)^j, reduced mod q.

    For j coprime to 3, j has inverse mod q. For j divisible by 3,
    write j = 3^v_3(j) · j', then 3^j / j has v_3 = j - v_3(j) ≥ 1.
    The factor 1/j' has inverse mod q. So contribution is
    (3^j / 3^{v_3(j)}) · (1/j') · (-1)^{j-1} · s^j mod q.
    """
    if s == 0:
        return 0
    val = 0
    for j in range(1, J + 1):
        # Compute v_3(j)
        jp = j
        v3j = 0
        while jp % 3 == 0:
            jp //= 3
            v3j += 1
        # j' = j / 3^{v_3(j)}; sign = (-1)^{j-1}
        sign = 1 if (j - 1) % 2 == 0 else -1
        # Term: (-1)^{j-1} / j · (3s)^j = sign · 3^j · s^j / j
        # = sign · 3^{j - v3j} · s^j / j'
        # mod q, use pow(j', -1, q)
        if 3 ** (j - v3j) >= q:
            # contribution is 0 mod q
            continue
        coeff = (sign * 3 ** (j - v3j)) % q
        # multiply by inverse of j'
        inv_jp = pow(jp, -1, q)
        s_pow = pow(s, j, q)
        val = (val + coeff * inv_jp % q * s_pow) % q
    return val


def compute_saddle_phase_lead(r: int, a: int) -> int:
    """ψ_lead(a) phase: P_a(s*(C_a)) mod q where s* = (C_a − 1)/3 mod 3.

    Returns the phase as an integer mod q (so ψ_lead(a) = e^{2πi·phase/q}).
    """
    q = 3 ** (r + 1)
    p_mm1 = 3 ** r
    J = J_for_p3(r + 1)
    # L(4) = L(1+3·1)
    L4 = truncated_3adic_log_mod_q(1, J, q)
    L_tilde = L4 // 3
    L_tilde_inv = pow(L_tilde, -1, p_mm1)
    C_a = (a * L_tilde_inv) % p_mm1
    s_star = ((C_a - 1) // 3) % 3
    # P_a(s*) = 3 s* − C_a · L(1+3 s*)
    L_at_sstar = truncated_3adic_log_mod_q(s_star, J, q)
    phase = (3 * s_star - C_a * L_at_sstar) % q
    return phase


# ----------------------------------------------------------------------
# Random-phase baseline
# ----------------------------------------------------------------------

@njit(cache=True, parallel=True)
def random_phase_baseline_sq(r: int, seed: int) -> float:
    """Baseline: Σ_{u=0}^{N-1} e^{i·θ_u} with θ_u uniform random.
    Returns |sum|² (E ≈ N for square-root-cancellation random model).
    """
    N = 3 ** (r - 1)
    n_blocks = 32
    block_size = (N + n_blocks - 1) // n_blocks
    block_re = np.zeros(n_blocks, dtype=np.float64)
    block_im = np.zeros(n_blocks, dtype=np.float64)
    for b in prange(n_blocks):
        u_start = b * block_size
        u_end = min(u_start + block_size, N)
        # Linear-congruential RNG seeded by (seed, b)
        rng_state = (seed * 1103515245 + 12345 + b * 6364136223846793005) & ((1 << 63) - 1)
        sre = 0.0
        sim = 0.0
        for _ in range(u_start, u_end):
            rng_state = (rng_state * 6364136223846793005 + 1442695040888963407) & ((1 << 63) - 1)
            theta = (rng_state >> 16) * (2.0 * math.pi / (1 << 47))
            sre += math.cos(theta)
            sim += math.sin(theta)
        block_re[b] = sre
        block_im[b] = sim
    re = block_re.sum()
    im = block_im.sum()
    return re * re + im * im


# ----------------------------------------------------------------------
# Sweep over r
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--r-min", type=int, default=8)
    parser.add_argument("--r-max", type=int, default=20)
    parser.add_argument("--n-c", type=int, default=20, help="number of c-units to sample for max")
    parser.add_argument("--n-m", type=int, default=5, help="number of m to sample (m=0..n_m-1)")
    parser.add_argument("--csv", type=str, default=r"C:\Collatz\r79b_S_partial_data.csv")
    parser.add_argument("--use-serial", action="store_true", help="use serial implementation")
    args = parser.parse_args()

    csv_path = args.csv
    write_header = not os.path.exists(csv_path)
    fieldnames = [
        "r", "q", "N", "K_c1m0_abs", "K_max_abs", "K_max_c", "K_max_m",
        "log_K_max", "log_sqrt_q", "rho", "delta_emp",
        "random_baseline_abs", "trivial_N", "elapsed_s",
    ]

    # Warm up Numba
    print("[warmup] compiling Numba kernels...")
    t0 = time.time()
    _ = kalafatelis_sum_direct(3, 1, 0)
    _ = kalafatelis_sum_serial(3, 1, 0)
    _ = random_phase_baseline_sq(3, 1)
    print(f"[warmup] done in {time.time()-t0:.1f}s")

    with open(csv_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
            f.flush()

        for r in range(args.r_min, args.r_max + 1):
            t_r = time.time()
            q = 3 ** (r + 1)
            N = 3 ** (r - 1)

            # Sample c-units (gcd(c, 3) = 1) and m
            rng = np.random.default_rng(r * 1009 + 17)
            cs_to_try = []
            cs_to_try.append(1)
            while len(cs_to_try) < args.n_c:
                c = int(rng.integers(2, q))
                if c % 3 != 0 and c not in cs_to_try:
                    cs_to_try.append(c)

            ms_to_try = list(range(args.n_m))

            # K at canonical (c=1, m=0)
            K_c1m0 = kalafatelis_sum_direct(r, 1, 0)
            K_c1m0_abs = abs(K_c1m0)

            # Max over (c, m) sample
            K_max_abs = 0.0
            K_max_c = 1
            K_max_m = 0
            for c in cs_to_try:
                for m in ms_to_try:
                    K = kalafatelis_sum_direct(r, c, m) if not args.use_serial else \
                        kalafatelis_sum_serial(r, c, m)
                    if abs(K) > K_max_abs:
                        K_max_abs = abs(K)
                        K_max_c = c
                        K_max_m = m

            # Random baseline (just one sample for shape)
            base_sq = random_phase_baseline_sq(r, seed=r + 1)
            base_abs = base_sq ** 0.5

            log_K_max = math.log(K_max_abs) if K_max_abs > 0 else float("-inf")
            log_sqrt_q = 0.5 * math.log(q)
            rho = log_K_max / log_sqrt_q - 1.0
            delta_emp = -rho / 2.0

            elapsed = time.time() - t_r

            row = {
                "r": r, "q": q, "N": N,
                "K_c1m0_abs": f"{K_c1m0_abs:.6f}",
                "K_max_abs": f"{K_max_abs:.6f}",
                "K_max_c": K_max_c, "K_max_m": K_max_m,
                "log_K_max": f"{log_K_max:.6f}",
                "log_sqrt_q": f"{log_sqrt_q:.6f}",
                "rho": f"{rho:.6f}",
                "delta_emp": f"{delta_emp:.6f}",
                "random_baseline_abs": f"{base_abs:.6f}",
                "trivial_N": N,
                "elapsed_s": f"{elapsed:.3f}",
            }
            writer.writerow(row)
            f.flush()
            print(f"r={r:>3} q=3^{r+1}={q:.3e} N={N:.3e}  |K|_c1m0={K_c1m0_abs:>10.3f}  "
                  f"|K|_max={K_max_abs:>10.3f} (c={K_max_c}, m={K_max_m})  "
                  f"ρ={rho:+.4f} δ={delta_emp:+.4f}  base={base_abs:>8.2f}  "
                  f"[{elapsed:.1f}s]")

    print(f"\n[done] csv = {csv_path}")


if __name__ == "__main__":
    main()
