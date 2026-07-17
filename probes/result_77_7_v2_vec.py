"""
result_77_7_v2_vec.py
=====================
Vectorized R77.7 v2 solver. Drop-in replacement for the inner functions of v2.

Key changes vs v2:
  - K_exp is built ONCE as flat numpy arrays (i_arr, j_arr, exp_arr) instead of
    N-dict-of-list-of-int. Memory: ~460 MB at k=8 vs 8 GB for the dict form.
  - Per-prime build A_p uses numpy gather + np.add.at instead of pure-Python
    triple loop. Speedup ~100×.
  - gauss_solve_mod_p uses outer-product elimination instead of per-pivot
    Python row loop. Speedup ~60×.

Expected at k=8 single-prime:
  - K_exp build (one-time): ~50s, ~460 MB
  - C1 build A_p (per prime): ~250s → ~2s
  - C2 gauss (per prime):    ~600s → ~10s
  - Per-prime: 14 min → ~12s
  - Full run @ 500 primes: 118 hr → ~100 min

Usage
-----
Verification (cheap, validates against v2 cache):
  python result_77_7_v2_vec.py --verify

Audit (single-prime k=8, ~1-2 min):
  python result_77_7_v2_vec.py --audit-k 8

Full production at k=8:
  python result_77_7_v2_vec.py --k 8

The CRT / rational reconstruction / eps-from-pi logic is reused verbatim from
result_77_7_v2.py via import. Only the inner per-prime hot paths are
re-implemented.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

# Reuse v2's CRT, rational reconstruction, eps-from-pi, prime selection
sys.path.insert(0, r"C:\Collatz")
from result_77_7_v2 import (
    crt_combine_vectorized,
    rational_reconstruct,
    is_prime,
    choose_primes,
    eps_from_pi_sequence,
    EPS_CACHE_V2,
    OUTDIR,
)

sys.stdout.reconfigure(encoding="utf-8")


# =========================================================================== #
# Vectorized K-flat construction                                              #
# =========================================================================== #

def build_K_flat(k: int):
    """
    Build K as flat int64 arrays. One entry per (i, j, exponent) triple.
    Equivalent to v2's build_K_exponents but with flat-array output.

    Returns:
        i_arr  : int64[E] — source row indices (i in K_exp)
        j_arr  : int64[E] — target column indices (where it gets added in A_p[j, i])
        exp_arr: int64[E] — exponents (each w_mod term is 2^exp_arr[k] mod p)
        M_param: int — geometric truncation
        coprime_states: list[int]
        N: int — matrix dimension
    """
    N_modulus = 3 ** k
    M_param = 2 * 3 ** (k - 1)
    inv2 = pow(2, -1, N_modulus)

    # Precompute (inv2)^v mod 3^k for v = 1..M
    powers_inv2 = [pow(inv2, v, N_modulus) for v in range(1, M_param + 1)]

    coprime_states = [r for r in range(N_modulus) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    N = len(coprime_states)

    # Total entries: N * M_param (each (i, r_v) gives one entry)
    total = N * M_param
    i_arr = np.empty(total, dtype=np.int64)
    j_arr = np.empty(total, dtype=np.int64)
    exp_arr = np.empty(total, dtype=np.int64)

    idx = 0
    for r in coprime_states:
        i = state_idx[r]
        three_r_plus_1 = 3 * r + 1
        for r_v in range(1, M_param + 1):
            target = (three_r_plus_1 * powers_inv2[r_v - 1]) % N_modulus
            j = state_idx[target]
            i_arr[idx] = i
            j_arr[idx] = j
            exp_arr[idx] = M_param - r_v
            idx += 1
    assert idx == total

    return i_arr, j_arr, exp_arr, M_param, coprime_states, N


# =========================================================================== #
# Vectorized per-prime A_p assembly                                           #
# =========================================================================== #

def build_A_p_vec(i_arr, j_arr, exp_arr, M_param, N, p, pow2_mod_arr=None):
    """
    Build the dense N x N matrix A_p = K^T - (2^M - 1) I mod p, last row replaced
    by all-ones, using vectorized numpy.

    Returns: A_p (int64 N x N) with last row = 1s, suitable for solving A_p x = e_N.
    """
    if pow2_mod_arr is None:
        pow2_mod_arr = np.array(
            [pow(2, e, p) for e in range(M_param + 1)], dtype=np.int64
        )

    D_mod = (int(pow2_mod_arr[M_param]) - 1) % p

    A_p = np.zeros((N, N), dtype=np.int64)
    # w[k] = 2^exp_arr[k] mod p
    w = pow2_mod_arr[exp_arr]
    # Scatter-add: A_p[j_arr[k], i_arr[k]] += w[k]
    # Note: np.add.at is the only way to handle duplicate (j, i) pairs correctly.
    # In our case duplicates can occur because multiple (i, r_v) pairs can map to
    # the same (j, i) when different r_v values give the same target.
    np.add.at(A_p, (j_arr, i_arr), w)
    A_p %= p

    # Subtract D_mod on diagonal: A = K^T - (2^M - 1) I
    diag_idx = np.arange(N)
    A_p[diag_idx, diag_idx] = (A_p[diag_idx, diag_idx] - D_mod) % p

    # Replace last row with all-ones (uniqueness constraint / normalization)
    A_p[N - 1, :] = 1

    return A_p


# =========================================================================== #
# Vectorized Gauss elimination mod p                                          #
# =========================================================================== #

def gauss_solve_mod_p_vec(A: np.ndarray, b: np.ndarray, p: int) -> np.ndarray | None:
    """
    Solve A x ≡ b (mod p) using Gauss elimination, with vectorized
    outer-product elimination (replaces v2's per-pivot Python row loop).

    Constraint: p^2 < 2^63 / 2 → p < ~2*10^9 (same as v2). Outer product
    factors[r] * M[col][j] stays in int64 range.

    Returns x as int64 ndarray of length N, or None if A is singular mod p.
    """
    N = A.shape[0]
    # Stack [A | b] as (N x (N+1)) for joint elimination
    M = np.empty((N, N + 1), dtype=np.int64)
    M[:, :N] = A % p
    M[:, N] = b % p

    for col in range(N):
        # Find pivot — first nonzero in column col, rows col..N-1
        col_segment = M[col:, col]
        nz = np.flatnonzero(col_segment)
        if len(nz) == 0:
            return None  # singular
        pivot = col + int(nz[0])
        if pivot != col:
            # Swap rows col and pivot
            tmp = M[col].copy()
            M[col] = M[pivot]
            M[pivot] = tmp

        # Normalize pivot row: M[col] *= inv(M[col, col]) mod p
        piv = int(M[col, col])
        inv_piv = pow(piv, p - 2, p)
        M[col] = (M[col] * inv_piv) % p

        # Eliminate ALL other rows via outer product:
        # M -= factors[:, None] * M[col][None, :]
        # where factors[r] = M[r, col] for r != col, else 0.
        factors = M[:, col].copy()
        factors[col] = 0
        # Outer product is (N, N+1) int64. Each entry < p^2 < 2^63 — safe.
        M -= np.outer(factors, M[col])
        M %= p

    return M[:, N].copy()


# =========================================================================== #
# Vectorized solve_pi_mod_p (wraps build_A_p_vec + gauss_solve_mod_p_vec)     #
# =========================================================================== #

def solve_pi_mod_p_vec(i_arr, j_arr, exp_arr, M_param, N, p, pow2_mod_arr=None):
    """Combined build A_p + solve. Returns pi_p as int64 ndarray of length N, or None."""
    A_p = build_A_p_vec(i_arr, j_arr, exp_arr, M_param, N, p, pow2_mod_arr)
    b_p = np.zeros(N, dtype=np.int64)
    b_p[N - 1] = 1
    return gauss_solve_mod_p_vec(A_p, b_p, p)


# =========================================================================== #
# Vectorized compute_pi_k_exact                                                #
# =========================================================================== #

def compute_pi_k_exact_vec(k: int, n_primes_initial: int = 30,
                           n_primes_batch_extra: int = 10,
                           max_primes: int = 1000,
                           verbose: bool = True,
                           log_func=None):
    """
    Same logic as v2.compute_pi_k_exact, but uses flat K-arrays + vectorized
    per-prime solver.
    """
    def L(msg):
        if log_func is not None:
            log_func(msg)
        elif verbose:
            print(msg, flush=True)

    L(f"  [k={k}] Building K_flat...")
    t_b = time.time()
    i_arr, j_arr, exp_arr, M_param, coprime_states, N = build_K_flat(k)
    L(f"  [k={k}] N={N}, M={M_param}, {len(i_arr):,} triples, "
      f"K_flat built in {time.time() - t_b:.2f}s, "
      f"~{(i_arr.nbytes + j_arr.nbytes + exp_arr.nbytes) / 1024**2:.0f} MB")

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
            pi_p = solve_pi_mod_p_vec(i_arr, j_arr, exp_arr, M_param, N, pp)
            elapsed = time.time() - t0
            t_solves.append(elapsed)
            if pi_p is None:
                L(f"  [k={k}] prime {pp}: SINGULAR mod p, skipping")
                continue
            residue_vectors.append(pi_p)
            moduli_used.append(pp)
            L(f"  [k={k}] +prime {len(moduli_used)} = {pp}: solved in {elapsed:.2f}s")

    add_primes(n_primes_initial)

    while True:
        L(f"  [k={k}] CRT-combining {len(moduli_used)} primes...")
        t_c = time.time()
        x_arr, P = crt_combine_vectorized(residue_vectors, moduli_used)
        L(f"  [k={k}] CRT done in {time.time() - t_c:.2f}s, P has {len(bin(P)) - 2} bits")

        L(f"  [k={k}] Reconstructing {N} components...")
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
        L(f"  [k={k}] reconstruction done in {time.time() - t_r:.2f}s, "
          f"{fail_count}/{N} failed components")

        verified = (fail_count == 0)

        if verified:
            s = sum(pi_frac)
            sum_ok = (s == 1)
            L(f"  [k={k}] sum(pi) == 1: {sum_ok}")
            if not sum_ok:
                verified = False

        if verified:
            L(f"  [k={k}] Adding 1 witness prime for verification...")
            wit_primes = choose_primes(1, start_above=next_prime_search_above)
            next_prime_search_above = wit_primes[-1]
            w = wit_primes[0]
            t0 = time.time()
            pi_w = solve_pi_mod_p_vec(i_arr, j_arr, exp_arr, M_param, N, w)
            t_solves.append(time.time() - t0)
            if pi_w is None:
                L(f"  [k={k}] witness prime {w} singular, adding another")
                continue
            all_match = True
            witness_invalid = False
            for i in range(N):
                f = pi_frac[i]
                num_mod_w = f.numerator % w
                den_mod_w = f.denominator % w
                if den_mod_w == 0:
                    witness_invalid = True
                    L(f"  [k={k}] WARN: witness {w} divides denominator at comp {i}; picking another")
                    break
                expected = num_mod_w * pow(den_mod_w, w - 2, w) % w
                if expected != int(pi_w[i]) % w:
                    all_match = False
                    L(f"  [k={k}] witness MISMATCH at component {i}: "
                      f"reconstructed={expected}, actual={int(pi_w[i]) % w}")
                    break
            if witness_invalid:
                residue_vectors.append(pi_w)
                moduli_used.append(w)
                continue
            if all_match:
                L(f"  [k={k}] Witness verification PASSED (prime {w})")
                residue_vectors.append(pi_w)
                moduli_used.append(w)
                break
            residue_vectors.append(pi_w)
            moduli_used.append(w)

        if len(moduli_used) >= max_primes:
            raise RuntimeError(f"Hit max_primes={max_primes} without convergence at k={k}")

        L(f"  [k={k}] Adding {n_primes_batch_extra} more primes...")
        add_primes(n_primes_batch_extra)

    L(f"  [k={k}] DONE — avg per-prime solve {np.mean(t_solves):.2f}s, "
      f"total primes used {len(moduli_used)}")

    return pi_frac, coprime_states, t_solves, len(moduli_used)


# =========================================================================== #
# Verification against v2 cache                                                #
# =========================================================================== #

def load_cached_eps():
    if not os.path.exists(EPS_CACHE_V2):
        raise RuntimeError(f"Cache not found: {EPS_CACHE_V2}")
    with open(EPS_CACHE_V2, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    return {int(k): Fraction(int(d["num"]), int(d["den"])) for k, d in data.items()}


def verify_against_cache(verify_k: int = 5):
    """Run vec at k=verify_k and check eps_k matches v2 cache."""
    print("=" * 70, flush=True)
    print(f"VERIFICATION: vec vs v2 cache at k={verify_k}", flush=True)
    print("=" * 70, flush=True)

    cached = load_cached_eps()
    if verify_k not in cached:
        raise RuntimeError(f"verify_k={verify_k} not in cache; pick another")

    expected_eps = cached[verify_k]
    print(f"Expected eps_{verify_k} from cache:", flush=True)
    print(f"  num = {expected_eps.numerator}", flush=True)
    print(f"  den = {expected_eps.denominator}", flush=True)
    print(f"  ~ {float(expected_eps):+.15e}", flush=True)
    print(flush=True)

    # Compute all pi_1..pi_verify_k via vec, then derive eps
    pi_dict = {}
    for kk in range(1, verify_k + 1):
        t0 = time.time()
        pi_k, _coprime, _t, _n = compute_pi_k_exact_vec(kk, n_primes_initial=20,
                                                       verbose=False)
        elapsed = time.time() - t0
        pi_dict[kk] = pi_k
        print(f"  k={kk} solved in {elapsed:.2f}s", flush=True)

    eps_vec = eps_from_pi_sequence(pi_dict)
    actual_eps = eps_vec[verify_k]

    print(flush=True)
    print(f"Vec eps_{verify_k}:", flush=True)
    print(f"  num = {actual_eps.numerator}", flush=True)
    print(f"  den = {actual_eps.denominator}", flush=True)
    print(f"  ~ {float(actual_eps):+.15e}", flush=True)
    print(flush=True)

    match = (actual_eps == expected_eps)
    print(f"MATCH: {match}", flush=True)
    if not match:
        diff = actual_eps - expected_eps
        print(f"  diff = {diff} ~ {float(diff):+.6e}", flush=True)
        # Check all k for diagnosis
        for kk in range(1, verify_k + 1):
            print(f"  k={kk}: vec={eps_vec[kk]} cache={cached.get(kk, 'N/A')} "
                  f"match={eps_vec[kk] == cached.get(kk)}", flush=True)

    return match


# =========================================================================== #
# Multiprocessing pool — parallelize over primes                              #
# =========================================================================== #
#
# K_flat (i_arr, j_arr, exp_arr) is immutable across primes (depends only on k),
# so it's the natural shared-state for the worker pool. Workers receive it once
# at initialization (Pool initargs); each worker keeps its own copy in module
# globals. Per-prime work allocates ~440 MB working memory (A_p + Gauss M +
# outer-product temp). With 16 workers on 64 GB: ~7.4 GB K_flat copies +
# ~7 GB working memory = ~15 GB. Comfortable headroom.
#
# Windows spawn semantics: each worker must be able to import this module and
# attach to global state. The init function sets globals; the worker function
# reads them.

# Module-level globals for worker state (set by _init_worker)
_g_i_arr = None
_g_j_arr = None
_g_exp_arr = None
_g_M_param = None
_g_N = None


def _init_worker(i_arr, j_arr, exp_arr, M_param, N):
    """Pool initializer — runs once per worker process, pickle args populate globals."""
    global _g_i_arr, _g_j_arr, _g_exp_arr, _g_M_param, _g_N
    _g_i_arr = i_arr
    _g_j_arr = j_arr
    _g_exp_arr = exp_arr
    _g_M_param = M_param
    _g_N = N


def _worker_solve(prime):
    """Pool worker function — solve mod a single prime, return (prime, pi_p, elapsed)."""
    t0 = time.time()
    try:
        pi_p = solve_pi_mod_p_vec(_g_i_arr, _g_j_arr, _g_exp_arr,
                                  _g_M_param, _g_N, prime)
    except Exception as e:
        return (prime, None, time.time() - t0, f"EXCEPTION: {type(e).__name__}: {e}")
    elapsed = time.time() - t0
    return (prime, pi_p, elapsed, None)


def compute_pi_k_exact_vec_pool(k: int, n_workers: int = 16,
                                n_primes_initial: int = 30,
                                n_primes_batch_extra: int = 16,
                                max_primes: int = 2000,
                                verbose: bool = True):
    """
    Parallelized version of compute_pi_k_exact_vec. Dispatches per-prime solves
    to a multiprocessing Pool.

    n_workers: pool size (typically 16 = physical cores on 9950X3D)
    n_primes_initial: initial batch (rounded up to multiple of n_workers)
    n_primes_batch_extra: per-retry batch if reconstruction fails
    """
    from multiprocessing import Pool, cpu_count
    n_workers = min(n_workers, cpu_count())

    def L(msg=""):
        if verbose:
            print(msg, flush=True)

    L(f"  [k={k}] Building K_flat...")
    t_b = time.time()
    i_arr, j_arr, exp_arr, M_param, coprime_states, N = build_K_flat(k)
    L(f"  [k={k}] N={N}, M={M_param}, {len(i_arr):,} triples, "
      f"K_flat built in {time.time() - t_b:.2f}s, "
      f"~{(i_arr.nbytes + j_arr.nbytes + exp_arr.nbytes) / 1024**2:.0f} MB")

    residue_vectors = []
    moduli_used = []
    t_solves = []
    next_prime_search_above = 10 ** 9 - 1

    L(f"  [k={k}] Spawning Pool with {n_workers} workers...")
    t_pool_start = time.time()
    pool = Pool(processes=n_workers,
                initializer=_init_worker,
                initargs=(i_arr, j_arr, exp_arr, M_param, N))
    L(f"  [k={k}] Pool ready in {time.time() - t_pool_start:.1f}s")

    def dispatch_primes(n_to_add):
        nonlocal next_prime_search_above
        new_primes = choose_primes(n_to_add, start_above=next_prime_search_above)
        next_prime_search_above = new_primes[-1]

        L(f"  [k={k}] Dispatching {len(new_primes)} primes to pool...")
        t_dispatch = time.time()
        # imap_unordered for live progress as workers complete
        n_done = 0
        for prime, pi_p, elapsed, err in pool.imap_unordered(_worker_solve, new_primes):
            n_done += 1
            t_solves.append(elapsed)
            if err is not None:
                L(f"  [k={k}] prime {prime}: ERROR {err}")
                continue
            if pi_p is None:
                L(f"  [k={k}] prime {prime}: SINGULAR mod p, skipping")
                continue
            residue_vectors.append(pi_p)
            moduli_used.append(prime)
            # Progress every ~1 worker-batch or every 30s
            if n_done % n_workers == 0 or n_done == len(new_primes):
                L(f"  [k={k}] progress {n_done}/{len(new_primes)} "
                  f"({len(moduli_used)} usable) | "
                  f"batch wall {time.time() - t_dispatch:.1f}s | "
                  f"last solve {elapsed:.1f}s")
        L(f"  [k={k}] Batch complete: {len(new_primes)} dispatched, "
          f"{len([_ for _ in new_primes])} processed, "
          f"wall {time.time() - t_dispatch:.1f}s, "
          f"effective per-prime {(time.time() - t_dispatch) / len(new_primes):.2f}s "
          f"(serial would be {sum(t_solves[-len(new_primes):]) / len(new_primes):.2f}s/prime, "
          f"speedup {sum(t_solves[-len(new_primes):]) / (time.time() - t_dispatch):.1f}×)")

    try:
        # Round initial up to multiple of n_workers
        n_init_round = ((n_primes_initial + n_workers - 1) // n_workers) * n_workers
        dispatch_primes(n_init_round)

        while True:
            L(f"  [k={k}] CRT-combining {len(moduli_used)} primes...")
            t_c = time.time()
            x_arr, P = crt_combine_vectorized(residue_vectors, moduli_used)
            L(f"  [k={k}] CRT done in {time.time() - t_c:.2f}s, P has {len(bin(P)) - 2} bits")

            L(f"  [k={k}] Reconstructing {N} components...")
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
            L(f"  [k={k}] reconstruction done in {time.time() - t_r:.2f}s, "
              f"{fail_count}/{N} failed components")

            verified = (fail_count == 0)

            if verified:
                s = sum(pi_frac)
                sum_ok = (s == 1)
                L(f"  [k={k}] sum(pi) == 1: {sum_ok}")
                if not sum_ok:
                    verified = False

            if verified:
                L(f"  [k={k}] Adding 1 witness prime for verification...")
                wit_primes = choose_primes(1, start_above=next_prime_search_above)
                next_prime_search_above = wit_primes[-1]
                w = wit_primes[0]
                # Witness via pool (single worker but uses same dispatch path)
                results_w = pool.map(_worker_solve, [w])
                prime, pi_w, elapsed, err = results_w[0]
                t_solves.append(elapsed)
                if err is not None or pi_w is None:
                    L(f"  [k={k}] witness prime {w} failed: {err or 'singular'}, retry")
                    continue
                all_match = True
                witness_invalid = False
                for i in range(N):
                    f = pi_frac[i]
                    num_mod_w = f.numerator % w
                    den_mod_w = f.denominator % w
                    if den_mod_w == 0:
                        witness_invalid = True
                        L(f"  [k={k}] WARN: witness {w} divides denominator at comp {i}; new witness next iter")
                        break
                    expected = num_mod_w * pow(den_mod_w, w - 2, w) % w
                    if expected != int(pi_w[i]) % w:
                        all_match = False
                        L(f"  [k={k}] witness MISMATCH at component {i}: "
                          f"reconstructed={expected}, actual={int(pi_w[i]) % w}")
                        break
                if witness_invalid:
                    residue_vectors.append(pi_w)
                    moduli_used.append(w)
                    continue
                if all_match:
                    L(f"  [k={k}] Witness verification PASSED (prime {w})")
                    residue_vectors.append(pi_w)
                    moduli_used.append(w)
                    break
                residue_vectors.append(pi_w)
                moduli_used.append(w)

            if len(moduli_used) >= max_primes:
                raise RuntimeError(f"Hit max_primes={max_primes} without convergence at k={k}")

            L(f"  [k={k}] Adding {n_primes_batch_extra} more primes...")
            dispatch_primes(n_primes_batch_extra)

    finally:
        pool.close()
        pool.join()

    L(f"  [k={k}] DONE — primes used {len(moduli_used)}, "
      f"avg per-prime solve {np.mean(t_solves):.2f}s")
    return pi_frac, coprime_states, t_solves, len(moduli_used)


def production_run_pool(k: int, n_workers: int = 16,
                        n_primes_initial: int = 500):
    """Pool-parallelized production run at k. Saves to v2_vec_pool cache."""
    cache_path = os.path.join(OUTDIR,
        f"result_77_7_eps_exact_through_k{k}_v2_vec_pool.json")

    print("=" * 70, flush=True)
    print(f"R77.7 v2 VEC POOL production run at k={k}, n_workers={n_workers}", flush=True)
    print("=" * 70, flush=True)
    print(f"Output cache: {cache_path}", flush=True)
    print(flush=True)

    cached = load_cached_eps()
    print(f"Loaded {len(cached)} eps from v2 cache: keys {sorted(cached.keys())}", flush=True)

    pi_dict = {}
    for kk in range(1, k + 1):
        print(flush=True)
        print(f"--- k = {kk} ---", flush=True)
        t0 = time.time()
        if kk <= 5:
            npi = 30
        elif kk == 6:
            npi = 60
        elif kk == 7:
            npi = 200
        else:
            npi = n_primes_initial

        # Use pool only for expensive k; serial is fine for small k
        if kk >= 7:
            pi_k, _, _t, _n = compute_pi_k_exact_vec_pool(kk, n_workers=n_workers,
                                                         n_primes_initial=npi)
        else:
            pi_k, _, _t, _n = compute_pi_k_exact_vec(kk, n_primes_initial=npi,
                                                    verbose=False)
        elapsed = time.time() - t0
        pi_dict[kk] = pi_k
        print(f"  k={kk} done in {elapsed:.2f}s = {elapsed/60:.2f} min", flush=True)

    eps = eps_from_pi_sequence(pi_dict)
    out = {str(kk): {"num": str(v.numerator), "den": str(v.denominator)}
           for kk, v in eps.items()}
    with open(cache_path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)

    print(flush=True)
    print("=" * 70, flush=True)
    print(f"DONE. eps values:", flush=True)
    for kk in sorted(eps.keys()):
        e = eps[kk]
        print(f"  eps_{kk} = {e.numerator}/{e.denominator}  "
              f"~ {float(e):+.15e}  |eps|*2^k = {abs(float(e)) * 2**kk:+.6f}",
              flush=True)
    print(flush=True)
    print(f"Saved to: {cache_path}", flush=True)


# =========================================================================== #
# Audit (single-prime instrumentation at given k)                              #
# =========================================================================== #

def audit_single_prime(k: int = 8, prime: int | None = None,
                       log_path: str | None = None):
    """Run one prime at k with phase timings — compare to non-vec audit."""
    import psutil

    if log_path is None:
        log_path = rf"C:\Collatz\experiments_output\audit_r77_v2_vec_k{k}.log"
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)

    fh = open(log_path, "w", encoding="utf-8", buffering=1)

    def L(msg=""):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        fh.write(line + "\n")
        fh.flush()
        os.fsync(fh.fileno())
        print(line, flush=True)

    L("=" * 70)
    L(f"R77.7 v2 VEC audit — single-prime k={k}")
    L("=" * 70)
    L(f"PID = {os.getpid()}")
    L(f"Python {sys.version.split()[0]}, numpy {np.__version__}")

    proc = psutil.Process(os.getpid())
    L(f"Initial RSS: {proc.memory_info().rss / 1024**2:.1f} MB")
    L("")

    # Phase A: build K_flat
    L("Phase A: build_K_flat")
    t0 = time.time()
    rss0 = proc.memory_info().rss / 1024**2
    i_arr, j_arr, exp_arr, M_param, coprime_states, N = build_K_flat(k)
    t_A = time.time() - t0
    rss1 = proc.memory_info().rss / 1024**2
    K_flat_mb = (i_arr.nbytes + j_arr.nbytes + exp_arr.nbytes) / 1024**2
    L(f"  wall: {t_A:.2f}s")
    L(f"  RSS: {rss0:.1f} -> {rss1:.1f} MB (delta {rss1 - rss0:+.1f})")
    L(f"  N={N}, M_param={M_param}, {len(i_arr):,} triples, K_flat ~{K_flat_mb:.0f} MB")
    L("")

    if prime is None:
        prime = 1_000_000_007
        while not is_prime(prime):
            prime += 2

    # Phase B: pow2_mod table
    L(f"Phase B: pow2_mod table p={prime}")
    t0 = time.time()
    rss0 = proc.memory_info().rss / 1024**2
    pow2_mod_arr = np.array(
        [pow(2, e, prime) for e in range(M_param + 1)], dtype=np.int64
    )
    t_B = time.time() - t0
    rss1 = proc.memory_info().rss / 1024**2
    L(f"  wall: {t_B:.2f}s, RSS delta {rss1 - rss0:+.1f} MB")
    L("")

    # Phase C1: build A_p (vectorized)
    L("Phase C1: build_A_p_vec")
    t0 = time.time()
    rss0 = proc.memory_info().rss / 1024**2
    A_p = build_A_p_vec(i_arr, j_arr, exp_arr, M_param, N, prime, pow2_mod_arr)
    t_C1 = time.time() - t0
    rss1 = proc.memory_info().rss / 1024**2
    L(f"  wall: {t_C1:.2f}s")
    L(f"  RSS: {rss0:.1f} -> {rss1:.1f} MB (delta {rss1 - rss0:+.1f})")
    L(f"  A_p size: {A_p.nbytes / 1024**2:.1f} MB")
    L("")

    # Phase C2: vectorized Gauss
    L("Phase C2: gauss_solve_mod_p_vec")
    t0 = time.time()
    rss0 = proc.memory_info().rss / 1024**2
    b_p = np.zeros(N, dtype=np.int64)
    b_p[N - 1] = 1
    x_p = gauss_solve_mod_p_vec(A_p, b_p, prime)
    t_C2 = time.time() - t0
    rss1 = proc.memory_info().rss / 1024**2
    L(f"  wall: {t_C2:.2f}s")
    L(f"  RSS: {rss0:.1f} -> {rss1:.1f} MB (delta {rss1 - rss0:+.1f})")
    L(f"  solution: {'OK' if x_p is not None else 'SINGULAR'}")
    L("")

    # Extrapolation
    per_prime_s = t_C1 + t_C2
    n_primes = 500
    total_s = per_prime_s * n_primes
    L(f"Per-prime: C1 {t_C1:.2f}s + C2 {t_C2:.2f}s = {per_prime_s:.2f}s")
    L(f"Full k={k} @ {n_primes} primes: {total_s:.0f}s = {total_s/60:.1f} min = {total_s/3600:.2f} hr")
    L(f"Memory ceiling: {rss1:.1f} MB")
    L("")

    # Speedup vs non-vec audit (if we have it)
    nonvec_log = r"C:\Collatz\experiments_output\audit_r77_v2_k8.log"
    if os.path.exists(nonvec_log):
        L("Comparison to non-vec audit (read from prior log)")
        L(f"  See: {nonvec_log}")
    L("")

    fh.close()
    return t_A, t_C1, t_C2


# =========================================================================== #
# Production main                                                              #
# =========================================================================== #

def production_run(k: int, n_primes_initial: int = 500):
    """Run full production at k, save cache."""
    cache_path = os.path.join(OUTDIR, f"result_77_7_eps_exact_through_k{k}_v2_vec.json")

    print("=" * 70, flush=True)
    print(f"R77.7 v2 VEC production run at k={k}", flush=True)
    print("=" * 70, flush=True)
    print(f"Output cache: {cache_path}", flush=True)
    print(flush=True)

    # Load existing cache to skip already-computed levels
    cached = load_cached_eps()
    print(f"Loaded {len(cached)} eps from v2 cache: keys {sorted(cached.keys())}", flush=True)

    # Compute all pi_1..pi_k via vec
    pi_dict = {}
    for kk in range(1, k + 1):
        print(flush=True)
        print(f"--- k = {kk} ---", flush=True)
        t0 = time.time()
        # Use smaller n_primes for small k (k<=5 fast)
        if kk <= 5:
            npi = 30
        elif kk == 6:
            npi = 60
        elif kk == 7:
            npi = 200
        else:
            npi = n_primes_initial
        pi_k, _coprime, _t, _n = compute_pi_k_exact_vec(kk, n_primes_initial=npi)
        elapsed = time.time() - t0
        pi_dict[kk] = pi_k
        print(f"  k={kk} done in {elapsed:.2f}s = {elapsed/60:.2f} min", flush=True)

    # Derive eps
    eps = eps_from_pi_sequence(pi_dict)

    # Save cache
    out = {str(kk): {"num": str(v.numerator), "den": str(v.denominator)}
           for kk, v in eps.items()}
    with open(cache_path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)

    print(flush=True)
    print("=" * 70, flush=True)
    print(f"DONE. eps values:", flush=True)
    for kk in sorted(eps.keys()):
        print(f"  eps_{kk} = {eps[kk].numerator}/{eps[kk].denominator}  "
              f"~ {float(eps[kk]):+.15e}  |eps|*2^k = {abs(float(eps[kk])) * 2**kk:+.6f}",
              flush=True)
    print(flush=True)
    print(f"Saved to: {cache_path}", flush=True)


# =========================================================================== #
# CLI                                                                          #
# =========================================================================== #

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true",
                        help="Verify against v2 cache at k=5")
    parser.add_argument("--verify-k", type=int, default=5,
                        help="k for verification (default 5)")
    parser.add_argument("--audit-k", type=int, default=None,
                        help="Run single-prime audit at this k")
    parser.add_argument("--k", type=int, default=None,
                        help="Production run at this k")
    parser.add_argument("--n-primes", type=int, default=500,
                        help="Initial primes for k=8 (default 500)")
    parser.add_argument("--pool", action="store_true",
                        help="Use multiprocessing pool for parallel primes")
    parser.add_argument("--n-workers", type=int, default=16,
                        help="Number of pool workers (default 16)")
    args = parser.parse_args()

    if args.verify:
        ok = verify_against_cache(verify_k=args.verify_k)
        sys.exit(0 if ok else 1)

    if args.audit_k is not None:
        audit_single_prime(k=args.audit_k)
        return

    if args.k is not None:
        if args.pool:
            production_run_pool(k=args.k, n_workers=args.n_workers,
                                n_primes_initial=args.n_primes)
        else:
            production_run(k=args.k, n_primes_initial=args.n_primes)
        return

    parser.print_help()


if __name__ == "__main__":
    main()
