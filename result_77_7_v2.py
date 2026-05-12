"""
result_77_7_v2.py
=================

R77.7 V2 — New solver for exact-rational eps_n via CRT + modular sparse-dense
Gauss elimination. Replaces the original Fraction-Gauss approach that was killed
at 8.5hr wall-time on k=7.

Pipeline:
  1) Build the integer matrix K_int (with common denominator 2^M - 1).
  2) For each of ~30 primes p in [1e8, 2e9], reduce mod p, Gauss-eliminate
     over F_p with numpy int64 arithmetic. Get pi_p in F_p^N.
  3) CRT-combine residues across primes to recover pi as integers mod P.
  4) Rational-reconstruct each component to recover pi as exact Fraction.
  5) Compute eps_k = (3^k * sum pi^2) - X_{k-1} - 7/15.
  6) Verify against cache for k=2..6; only then run k=7.

Output:
  C:/Collatz/experiments_output/result_77_7_eps_exact_through_k7_v2.json
  C:/Collatz/R77_7_V2_RUN_LOG.txt (stdout tee, written on run)

Runnable from main thread:
  python C:/Collatz/result_77_7_v2.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction
from math import gcd

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)
EPS_CACHE_OLD = os.path.join(OUTDIR, "result_77_7_eps_exact_through_k7.json")
EPS_CACHE_V2 = os.path.join(OUTDIR, "result_77_7_eps_exact_through_k7_v2.json")


# =========================================================================== #
# Markov-chain construction (integer version)                                 #
# =========================================================================== #

def build_K_exponents(k: int):
    """
    Build the Markov kernel as a list of exponent-lists. For each (i, j) where
    state i can transition to state j with weight 2^(M - r_v), record the list of
    exponents (M - r_v) for all r_v that produce this i->j transition.

    K(i, j) = (sum over recorded exponents e of 2^e) / (2^M - 1).

    Returns:
        K_exp : list of N rows, each a dict {j: list of int exponents}
        M     : int, geometric truncation parameter
        coprime_states : list of int, residues coprime to 3 in Z/3^k
        N     : int, number of coprime states
    """
    N_modulus = 3 ** k
    M = 2 * 3 ** (k - 1)
    inv2 = pow(2, -1, N_modulus)
    powers_inv2 = [pow(inv2, v, N_modulus) for v in range(1, M + 1)]
    coprime_states = [r for r in range(N_modulus) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    N = len(coprime_states)

    K_exp = [dict() for _ in range(N)]
    for r in coprime_states:
        i = state_idx[r]
        three_r_plus_1 = 3 * r + 1
        row_dict = K_exp[i]
        for r_v in range(1, M + 1):
            target = (three_r_plus_1 * powers_inv2[r_v - 1]) % N_modulus
            j = state_idx[target]
            e = M - r_v
            if j in row_dict:
                row_dict[j].append(e)
            else:
                row_dict[j] = [e]
    return K_exp, M, coprime_states, N


# =========================================================================== #
# Modular dense Gauss elimination                                             #
# =========================================================================== #

def gauss_solve_mod_p(A: np.ndarray, b: np.ndarray, p: int) -> np.ndarray | None:
    """
    Solve A x ≡ b (mod p) for dense A (N x N int64) using Gauss elimination.

    Returns x as int64 ndarray of length N, or None if A is singular mod p.

    Constraint: p must satisfy 2*(p-1)^2 < 2^63 to avoid int64 overflow in the
    row-subtract step. Practically: p < 2^31 is fully safe; we'll use p < ~2*10^9.
    """
    N = A.shape[0]
    # Stack [A | b] as (N x (N+1)) for joint elimination
    M = np.zeros((N, N + 1), dtype=np.int64)
    M[:, :N] = A % p
    M[:, N] = b % p

    for col in range(N):
        # Find pivot
        pivot = -1
        for row in range(col, N):
            if M[row, col] != 0:
                pivot = row
                break
        if pivot == -1:
            return None  # singular
        if pivot != col:
            tmp = M[col].copy()
            M[col] = M[pivot]
            M[pivot] = tmp

        # Normalize pivot row: divide by pivot value (multiply by inverse)
        piv = int(M[col, col])
        inv_piv = pow(piv, p - 2, p)
        # row scale
        M[col] = (M[col].astype(np.int64) * inv_piv) % p

        # Eliminate from all other rows
        for row in range(N):
            if row == col:
                continue
            factor = int(M[row, col])
            if factor == 0:
                continue
            # M[row] -= factor * M[col]   (mod p)
            # Use the vectorized form, careful with overflow:
            # factor < p and each M[col,j] < p, so factor * M[col,j] < p^2 < 2^63.
            M[row] = (M[row] - factor * M[col]) % p

    return M[:, N].copy()


# =========================================================================== #
# Build mod-p problem and solve                                               #
# =========================================================================== #

def solve_pi_mod_p(K_exp, M_param: int, N: int, p: int, pow2_mod=None):
    """
    Solve A π = e_N mod p, where A = K^T - (2^M - 1) I, with last row replaced
    by all-ones and RHS = e_N (last position = 1, others = 0).

    K_exp is a list of dicts: K_exp[i][j] = list of exponents [e_1, e_2, ...].
    The (i, j) entry of K_int = sum(2^e for e in exponents); mod p it's
    sum(pow2_mod[e] for e in exponents) % p, where pow2_mod[e] = 2^e mod p.

    Returns: numpy int64 vector of length N, or None if singular mod p.
    """
    if pow2_mod is None:
        pow2_mod = [pow(2, e, p) for e in range(M_param + 1)]

    D_mod = (pow2_mod[M_param] - 1) % p  # 2^M - 1 mod p

    # Build A_p as numpy int64 (N x N), starting from K^T mod p
    A_p = np.zeros((N, N), dtype=np.int64)
    for i in range(N):
        row_dict = K_exp[i]
        for j, exponents in row_dict.items():
            w_mod = 0
            for e in exponents:
                w_mod += pow2_mod[e]
            A_p[j, i] = (A_p[j, i] + w_mod) % p

    # Subtract D_mod on the diagonal: A = K^T - (2^M - 1) I
    for i in range(N):
        A_p[i, i] = (A_p[i, i] - D_mod) % p

    # Replace last row with all-ones, RHS = e_N
    A_p[N - 1, :] = 1
    b_p = np.zeros(N, dtype=np.int64)
    b_p[N - 1] = 1

    return gauss_solve_mod_p(A_p, b_p, p)


# =========================================================================== #
# CRT + rational reconstruction                                                #
# =========================================================================== #

def crt_combine(residues_list, moduli_list):
    """
    CRT: given residues[i] and moduli[i] coprime, return (x, P) where
    x ≡ residues[i] mod moduli[i] for all i, and P = product of moduli.

    residues_list and moduli_list: lists of Python ints.
    Returns: (x, P) Python ints, with 0 <= x < P.
    """
    x = 0
    P = 1
    for r, p in zip(residues_list, moduli_list):
        # Combine x mod P with r mod p
        # New solution: x + P * t, where t = (r - x) * P^{-1} mod p
        g, P_inv_p, _ = ext_gcd(P % p, p)
        assert g == 1, f"CRT: moduli not coprime, gcd(P mod {p}, {p}) = {g}"
        t = ((r - x) % p) * (P_inv_p % p) % p
        x = (x + P * t) % (P * p)
        P = P * p
    return x, P


def crt_combine_vectorized(residue_vectors, moduli_list):
    """
    Vectorized CRT across components of a vector. Each residue_vectors[i] is a
    1D numpy int64 array (or list) of length N. Returns a list of Python ints
    of length N representing the CRT-combined integer for each component, and
    the product P of moduli.
    """
    N = len(residue_vectors[0])
    K = len(moduli_list)
    # Accumulate component-wise
    x_arr = [0] * N
    P = 1
    for k_idx in range(K):
        p = moduli_list[k_idx]
        rv = residue_vectors[k_idx]
        if k_idx == 0:
            for i in range(N):
                x_arr[i] = int(rv[i]) % p
            P = p
        else:
            # P_inv_p = (P mod p) ^{-1} mod p
            P_inv_p = pow(P % p, p - 2, p)
            P_new = P * p
            for i in range(N):
                r = int(rv[i]) % p
                t = ((r - x_arr[i]) % p) * P_inv_p % p
                x_arr[i] = (x_arr[i] + P * t) % P_new
            P = P_new
    return x_arr, P


def ext_gcd(a, b):
    """Extended GCD: returns (g, x, y) with a*x + b*y = g = gcd(a, b)."""
    if b == 0:
        return a, 1, 0
    g, x1, y1 = ext_gcd(b, a % b)
    return g, y1, x1 - (a // b) * y1


def rational_reconstruct(u: int, m: int):
    """
    Given u mod m, find rational n/d with |n|, d < sqrt(m/2) such that n ≡ u*d (mod m),
    using the half-extended Euclidean algorithm.

    Returns (n, d) as Python ints with d > 0, or None if no such reconstruction exists.
    """
    if u == 0:
        return (0, 1)
    # Initialize
    r0, r1 = m, u % m
    s0, s1 = 0, 1
    # We want the first r_k with r_k <= sqrt(m/2); then n = r_k, d = s_k (sign-adjusted).
    # Bound: sqrt(m/2). We use integer square root via math.isqrt
    from math import isqrt
    bound = isqrt(m // 2)
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    # Now r1 <= bound. Need |s1| <= bound too for uniqueness.
    n = r1
    d = s1
    if d < 0:
        n = -n
        d = -d
    if d == 0:
        return None
    if abs(d) > bound:
        return None
    if gcd(abs(n), d) != 1:
        # not reduced, reduce
        g = gcd(abs(n), d)
        n //= g
        d //= g
    return (n, d)


# =========================================================================== #
# Choose primes                                                                #
# =========================================================================== #

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False
    # Miller-Rabin deterministic for n < 3.3e24 with these witnesses
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def choose_primes(n_primes: int, start_above: int = 10 ** 9 - 1, bound: int = 2 * 10 ** 9):
    """Pick n_primes primes in (start_above, bound], all < 2*10^9 to keep
    p^2 < 2^63 (well, 4e18 < 9.2e18, safe with margin)."""
    primes = []
    p = start_above + 1
    while len(primes) < n_primes:
        if is_prime(p):
            primes.append(p)
        p += 1
        if p > bound:
            raise RuntimeError(f"Ran out of primes below {bound}")
    return primes


# =========================================================================== #
# Top-level: compute eps_k                                                     #
# =========================================================================== #

def compute_pi_k_exact(k: int, n_primes_initial: int = 30,
                       n_primes_batch_extra: int = 10,
                       max_primes: int = 200,
                       verbose: bool = True):
    """
    Compute the exact-rational stationary distribution pi_k as a list of
    Fraction objects of length N_k.

    Strategy:
      1. Build K_exp once.
      2. Solve mod p for an initial batch of primes.
      3. Try CRT + rational reconstruction.
      4. Verify sum-to-1 + at least one extra prime not used in reconstruction
         agrees with the reconstructed rational.
      5. If verification fails OR reconstruction has failures, add more primes
         and retry.
    """
    if verbose:
        print(f"\n  [k={k}] Building K_exp...")
    t_b = time.time()
    K_exp, M_param, coprime_states, N = build_K_exponents(k)
    if verbose:
        print(f"  [k={k}] N={N}, M={M_param}, build time {time.time() - t_b:.2f}s")

    residue_vectors = []
    moduli_used = []
    t_solves = []
    next_prime_search_above = 10 ** 9 - 1

    def add_primes(n_to_add):
        nonlocal next_prime_search_above
        new_primes = choose_primes(n_to_add, start_above=next_prime_search_above)
        next_prime_search_above = new_primes[-1]
        for idx_p, pp in enumerate(new_primes):
            t0 = time.time()
            pi_p = solve_pi_mod_p(K_exp, M_param, N, pp)
            elapsed = time.time() - t0
            t_solves.append(elapsed)
            if pi_p is None:
                if verbose:
                    print(f"  [k={k}] prime {pp}: SINGULAR mod p, skipping")
                continue
            residue_vectors.append(pi_p)
            moduli_used.append(pp)
            if verbose:
                print(f"  [k={k}] +prime {len(moduli_used)} = {pp}: solved in {elapsed:.2f}s")

    add_primes(n_primes_initial)

    while True:
        # CRT + reconstruct using ALL primes except the last one (held out for verify)
        # Actually: use ALL primes for reconstruction, then test reconstruction by
        # evaluating each Fraction component mod a held-out prime and comparing.
        # We'll do this by adding one extra "witness" prime separately.
        if verbose:
            print(f"  [k={k}] CRT-combining {len(moduli_used)} primes...")
        t_c = time.time()
        x_arr, P = crt_combine_vectorized(residue_vectors, moduli_used)
        if verbose:
            print(f"  [k={k}] CRT done in {time.time() - t_c:.2f}s, "
                  f"P has {len(bin(P)) - 2} bits")

        if verbose:
            print(f"  [k={k}] Reconstructing {N} components...")
        t_r = time.time()
        pi_frac = []
        fail_count = 0
        for i in range(N):
            recon = rational_reconstruct(x_arr[i], P)
            if recon is None:
                fail_count += 1
                pi_frac.append(None)
            else:
                pi_frac.append(Fraction(recon[0], recon[1]))
        if verbose:
            print(f"  [k={k}] reconstruction done in {time.time() - t_r:.2f}s, "
                  f"{fail_count}/{N} failed components")

        verified = (fail_count == 0)

        if verified:
            # Sum-to-1 check
            s = sum(pi_frac)
            sum_ok = (s == 1)
            if verbose:
                print(f"  [k={k}] sum(pi) == 1: {sum_ok}")
            if not sum_ok:
                verified = False

        if verified:
            # Witness-prime check: add ONE more prime, solve mod it, and verify
            # that each reconstructed Fraction matches mod the witness.
            if verbose:
                print(f"  [k={k}] Adding 1 witness prime for verification...")
            wit_primes = choose_primes(1, start_above=next_prime_search_above)
            next_prime_search_above = wit_primes[-1]
            w = wit_primes[0]
            t0 = time.time()
            pi_w = solve_pi_mod_p(K_exp, M_param, N, w)
            t_solves.append(time.time() - t0)
            if pi_w is None:
                if verbose:
                    print(f"  [k={k}] witness prime {w} singular, adding another")
                continue
            # Check: each Fraction.num/Fraction.den mod w should equal pi_w[i]
            all_match = True
            witness_invalid = False
            for i in range(N):
                f = pi_frac[i]
                num_mod_w = f.numerator % w
                den_mod_w = f.denominator % w
                if den_mod_w == 0:
                    # w divides the denominator of one of our reconstructed
                    # fractions. The witness is "blind" at this component;
                    # we can't use it. Get another witness.
                    witness_invalid = True
                    if verbose:
                        print(f"  [k={k}] WARN: witness {w} divides denominator "
                              f"at comp {i}; picking another witness")
                    break
                expected = num_mod_w * pow(den_mod_w, w - 2, w) % w
                if expected != int(pi_w[i]) % w:
                    all_match = False
                    if verbose:
                        print(f"  [k={k}] witness MISMATCH at component {i}: "
                              f"reconstructed={expected}, actual={int(pi_w[i]) % w}")
                    break
            if witness_invalid:
                # Add this prime to the working set and try a new witness next iter
                residue_vectors.append(pi_w)
                moduli_used.append(w)
                continue
            if all_match:
                if verbose:
                    print(f"  [k={k}] Witness verification PASSED (prime {w})")
                # Add the witness to residue_vectors for future use
                residue_vectors.append(pi_w)
                moduli_used.append(w)
                break
            # Otherwise: witness failed; add the witness to set, add more primes, retry
            residue_vectors.append(pi_w)
            moduli_used.append(w)

        if len(moduli_used) >= max_primes:
            raise RuntimeError(f"Hit max_primes={max_primes} without convergence at k={k}")

        if verbose:
            print(f"  [k={k}] Adding {n_primes_batch_extra} more primes...")
        add_primes(n_primes_batch_extra)

    if verbose:
        print(f"  [k={k}] DONE — avg per-prime solve {np.mean(t_solves):.2f}s, "
              f"total primes used {len(moduli_used)}")

    return pi_frac, coprime_states, t_solves, len(moduli_used)


# =========================================================================== #
# Eps from pi                                                                  #
# =========================================================================== #

def eps_from_pi_sequence(pi_dict: dict):
    """
    pi_dict: dict {k: list of Fractions of length N_k} for k = 1, 2, ..., max_k.
    Returns: dict {k: Fraction(eps_k)} for k = 1..max_k.
    """
    target = Fraction(7, 15)
    eps = {}
    for k in sorted(pi_dict.keys()):
        pi = pi_dict[k]
        X_k = Fraction(3 ** k) * sum(p * p for p in pi)
        if k == 1:
            X_km1 = Fraction(1)
        else:
            sum_eps_prior = sum(eps[j] for j in range(1, k))
            X_km1 = Fraction(1) + Fraction(k - 1) * target + sum_eps_prior
        S_k = X_k - X_km1
        eps[k] = S_k - target
    return eps


# =========================================================================== #
# Verification against cached eps                                              #
# =========================================================================== #

def load_eps_cache():
    if not os.path.exists(EPS_CACHE_OLD):
        return {}
    with open(EPS_CACHE_OLD, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    return {int(k): Fraction(int(d["num"]), int(d["den"])) for k, d in data.items()}


def save_eps_cache_v2(eps_dict):
    out = {str(k): {"num": str(v.numerator), "den": str(v.denominator)}
           for k, v in eps_dict.items()}
    with open(EPS_CACHE_V2, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)


# =========================================================================== #
# Main                                                                         #
# =========================================================================== #

def main():
    print("=" * 78)
    print("R77.7 V2 — CRT + modular solver for exact eps_n through k=7")
    print("=" * 78)
    print(f"Output cache: {EPS_CACHE_V2}")
    print()

    cached_eps = load_eps_cache()
    print(f"Loaded {len(cached_eps)} cached eps from original solver:")
    for k in sorted(cached_eps.keys()):
        print(f"  eps_{k} ~ {float(cached_eps[k]):+.10e}")
    print()

    # Phase 4 verification: solve k=2..6, verify against cache
    pi_dict = {}
    timings = {}
    primes_used = {}

    # k=1 trivial: only 2 states (1, 2), pi = (1/2, 1/2) up to symmetry,
    # but actually we should compute it. N=2 is fine.
    print("\nPhase 4: verification k=1..6 against cache")
    print("-" * 78)
    # Empirical: ε_k denominator length grows ~3.3x per step. π_k component
    # denominators are roughly sqrt of that. With 30-bit primes (~9 decimal
    # digits each), we need n_primes ≈ 2·digits(den) / 9.
    # ε denominator digit counts: k=1: 1, k=2: 3, k=3: 7, k=4: 19, k=5: 65, k=6: 200
    # So π denominator digits estimate: 1, 2, 4, 10, 33, 100, 300+ (k=1..7)
    # Required primes ≈ 4·digits(π)/9 (4x for safety + sum-to-1 ratio etc):
    # Adaptive: initial batch is an under-estimate; the compute_pi_k_exact loop
    # adds 10 more at a time if reconstruction or witness verification fails.
    initial_primes_by_k = {1: 6, 2: 6, 3: 8, 4: 12, 5: 24, 6: 50, 7: 150}
    for k in range(1, 7):
        t_total = time.time()
        n_primes = initial_primes_by_k.get(k, 20)
        pi_k, coprime, t_solves, np_used = compute_pi_k_exact(k, n_primes_initial=n_primes)
        wall = time.time() - t_total
        pi_dict[k] = pi_k
        timings[k] = wall
        primes_used[k] = np_used
        print(f"  [k={k}] TOTAL wall time {wall:.2f}s, primes used {np_used}")
        # Incremental save: persist partial results after each k for crash safety
        try:
            partial_eps = eps_from_pi_sequence(pi_dict)
            save_eps_cache_v2(partial_eps)
        except Exception as e:
            print(f"  [k={k}] WARN: incremental save failed: {e}")

    eps_dict = eps_from_pi_sequence(pi_dict)
    print("\n  Computed eps_k from new solver:")
    for k in sorted(eps_dict.keys()):
        ek = eps_dict[k]
        print(f"    eps_{k} = num/den lens ({len(str(abs(ek.numerator)))},"
              f"{len(str(ek.denominator))})  ~ {float(ek):+.10e}")

    print("\n  Verifying against cache:")
    all_match = True
    for k in range(1, 7):
        if k not in cached_eps:
            print(f"    eps_{k}: not in cache, skipping verify")
            continue
        match = eps_dict[k] == cached_eps[k]
        symbol = "OK" if match else "FAIL"
        print(f"    eps_{k}: {symbol}  (computed={float(eps_dict[k]):+.10e}, "
              f"cached={float(cached_eps[k]):+.10e})")
        if not match:
            all_match = False

    if not all_match:
        print("\n  VERIFICATION FAILED at one or more of k=1..6. NOT running k=7.")
        # Still save what we have
        save_eps_cache_v2(eps_dict)
        return
    print("\n  All k=1..6 match cache. Verification PASSED.")

    # Estimate k=7 wall time from k=5, k=6 scaling
    if 5 in timings and 6 in timings:
        ratio = timings[6] / timings[5]
        est_k7 = timings[6] * ratio  # roughly N^3 scales by (1458/486)^3 = 27
        print(f"\n  Estimated k=7 wall time from k=6 scaling: ~{est_k7:.0f}s = ~{est_k7/60:.1f} min")
        # Sanity-cap: theoretical is 27x slower per prime, possibly with more primes too
        theory_est = timings[6] * 27 * 1.5  # 27 for N^3, 1.5 for more primes needed
        print(f"  Theoretical N^3 estimate from k=6: ~{theory_est:.0f}s = ~{theory_est/60:.1f} min")

    print("\nPhase 4 verification complete. Proceeding to k=7.\n")
    print("=" * 78)
    print("Phase 5: compute eps_7")
    print("=" * 78)

    t7_start = time.time()
    # k=7: use ~200 primes initially (denominators very big; safer to over-shoot
    # the initial batch than to retry the full pipeline). Each prime is ~10-15 sec
    # so initial budget is ~30-50 min wall.
    pi_7, coprime_7, t_solves_7, np_used_7 = compute_pi_k_exact(
        7, n_primes_initial=initial_primes_by_k[7])
    t7 = time.time() - t7_start
    pi_dict[7] = pi_7
    timings[7] = t7
    primes_used[7] = np_used_7

    eps_dict = eps_from_pi_sequence(pi_dict)

    print(f"\n  [k=7] TOTAL wall time {t7:.1f}s = {t7/60:.2f} min, "
          f"primes used {np_used_7}, mean per-prime solve {np.mean(t_solves_7):.2f}s")
    print(f"\n  eps_7 = {eps_dict[7].numerator} / {eps_dict[7].denominator}")
    print(f"  eps_7 (float) ~ {float(eps_dict[7]):+.15e}")
    print(f"  |eps_7| * 2^7 = {abs(float(eps_dict[7])) * 128:+.10f}")

    save_eps_cache_v2(eps_dict)
    print(f"\n  Saved eps_1..eps_7 to {EPS_CACHE_V2}")

    # Compare to prior numerical estimate from STATE.md (~-1.18e-3)
    e7 = float(eps_dict[7])
    print(f"\n  Comparison to project numerical (power-iteration) prior: ~-1.18e-3")
    print(f"  Computed exact (float) eps_7 = {e7:+.6e}")
    print(f"  Ratio |eps_7|/|eps_6| = {e7/float(eps_dict[6]):.4f}")
    print(f"  |eps_n| * 2^n series:")
    for kk in sorted(eps_dict.keys()):
        v = float(eps_dict[kk])
        print(f"    n={kk}: |eps_{kk}|*2^{kk} = {abs(v) * (2 ** kk):+.6f}")

    print("\n" + "=" * 78)
    print("Wall time summary")
    print("=" * 78)
    for k in sorted(timings.keys()):
        print(f"  k={k}: {timings[k]:.2f}s, {primes_used[k]} primes")


if __name__ == "__main__":
    main()
