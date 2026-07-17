"""
result_q_sweep_test_2.py
========================
Q-sweep test 2: compute c_q = S_inf^(q) / q across odd primes q in {3, 5, 7, 11, 13}.

Stage 1: For each q, build the qx+1 Markov chain at level k via the GENERALIZED
         build_markov_q (verified against q=3 in preflight).

Stage 2: Compute S_k^(q) = X_k^(q) - X_{k-1}^(q) for k = 1..k_max(q):
         k_max = 4 for q in {3, 5}
         k_max = 3 for q in {7, 11, 13}

Stage 3: Extrapolate S_inf^(q) per q via the convergence pattern.

Stage 4: Build (q, c_q = S_inf^(q) / q) table; test for closed form.

Stage 5: Outcome classification.

Cache: result_q_sweep_test_2_cache.json stores all computed S_k^(q) as exact
       Fractions; runs are resumable.
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)
CACHE = os.path.join(OUTDIR, "result_q_sweep_test_2_cache.json")
OUT_CSV = os.path.join(OUTDIR, "result_q_sweep_test_2_table.csv")


# --------------------------------------------------------------------------- #
# Generalized qx+1 Markov chain (verified in preflight)                       #
# --------------------------------------------------------------------------- #

def order_of_two(N: int) -> int:
    assert N % 2 == 1
    m = 1
    v = 2 % N
    while v != 1:
        v = (v * 2) % N
        m += 1
    return m


def build_markov_q(q: int, k: int):
    N = q ** k
    M = order_of_two(N)
    inv2 = pow(2, -1, N)
    powers_inv2 = [pow(inv2, v, N) for v in range(1, M + 1)]
    coprime_states = [r for r in range(N) if r % q != 0]
    state_idx = {r: i for i, r in enumerate(coprime_states)}
    n = len(coprime_states)
    K = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    Z_v = Fraction(2 ** M - 1, 2 ** M)
    for r in coprime_states:
        for r_v in range(1, M + 1):
            p = Fraction(1, 2 ** r_v) / Z_v
            target = ((q * r + 1) * powers_inv2[r_v - 1]) % N
            K[state_idx[r]][state_idx[target]] += p
    return K, coprime_states, M


def stationary_rational(K):
    n = len(K)
    A = [[K[j][i] - (Fraction(1) if i == j else Fraction(0)) for j in range(n)] for i in range(n)]
    A[n - 1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n - 1] = Fraction(1)
    for col in range(n):
        pivot = -1
        for row in range(col, n):
            if A[row][col] != 0:
                pivot = row
                break
        if pivot == -1:
            raise ValueError(f"Singular at col {col}")
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
    return b


# --------------------------------------------------------------------------- #
# Cache I/O                                                                   #
# --------------------------------------------------------------------------- #

def load_cache() -> dict:
    """Return {(q, k): (X_k as Fraction, n_states)}. Stored on disk as JSON."""
    if not os.path.exists(CACHE):
        return {}
    with open(CACHE, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    out = {}
    for key, val in data.items():
        q, k = map(int, key.split(","))
        Xk = Fraction(int(val["X_num"]), int(val["X_den"]))
        n_states = int(val["n_states"])
        out[(q, k)] = (Xk, n_states)
    return out


def save_cache(cache: dict):
    data = {}
    for (q, k), (Xk, n_states) in cache.items():
        data[f"{q},{k}"] = {
            "X_num": str(Xk.numerator),
            "X_den": str(Xk.denominator),
            "n_states": n_states,
        }
    with open(CACHE, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)


# --------------------------------------------------------------------------- #
# Stage 2: compute X_k^(q) for all desired (q, k)                              #
# --------------------------------------------------------------------------- #

def compute_X(q: int, k: int, cache: dict, log_prefix: str = "  "):
    """Return X_k^(q) as Fraction. Caches in `cache` and on disk."""
    if (q, k) in cache:
        Xk, n_states = cache[(q, k)]
        print(f"{log_prefix}q={q} k={k}: cached  X_{k} = {Xk}  (~ {float(Xk):.10f})  "
              f"({n_states} states)", flush=True)
        return Xk

    t0 = time.time()
    K, coprime, M = build_markov_q(q, k)
    t_build = time.time() - t0
    n = len(coprime)
    t0b = time.time()
    pi = stationary_rational(K)
    t_solve = time.time() - t0b
    Xk = Fraction(q ** k) * sum(p * p for p in pi)
    elapsed = time.time() - t0
    print(f"{log_prefix}q={q} k={k}: build {t_build:.1f}s  solve {t_solve:.1f}s  "
          f"M={M}  {n} states  X_{k} ~ {float(Xk):.10f}", flush=True)
    cache[(q, k)] = (Xk, n)
    save_cache(cache)
    return Xk


def compute_S_sequence(q: int, k_max: int, cache: dict):
    """Return [(k, S_k^(q) as Fraction) for k=1..k_max]."""
    Xs = {0: Fraction(1)}
    print(f"--- q = {q}, k_max = {k_max} ---")
    for k in range(1, k_max + 1):
        Xs[k] = compute_X(q, k, cache, log_prefix="  ")
    Ss = []
    for k in range(1, k_max + 1):
        Ss.append((k, Xs[k] - Xs[k - 1]))
    print(f"  S_k^({q}) sequence:")
    for k, S in Ss:
        print(f"    S_{k}^({q}) = {S}  ~ {float(S):.10f}")
    return Ss


# --------------------------------------------------------------------------- #
# Stage 3: extrapolate S_inf^(q)                                              #
# --------------------------------------------------------------------------- #

def extrapolate_S_inf(Ss, q: int):
    """Estimate S_inf^(q) from a finite S_k sequence.

    Strategy: fit log|S_k - L| = a + b*k for various candidate L; pick the L
    giving cleanest geometric decay (highest R^2 on log scale). At small k_max
    this is heuristic; report uncertainty.
    """
    ks = np.array([k for k, _ in Ss])
    Sf = np.array([float(S) for _, S in Ss])

    if len(Sf) < 2:
        return {"S_inf_est": float("nan"), "method": "insufficient", "uncertainty": float("nan")}

    if len(Sf) == 2:
        # No way to extrapolate; take last value as estimate
        return {
            "S_inf_est": Sf[-1],
            "method": "last_value (only 2 points)",
            "uncertainty": abs(Sf[-1] - Sf[-2]),
        }

    # Aitken Delta^2 acceleration
    aitken = []
    for i in range(len(Sf) - 2):
        denom = Sf[i + 2] - 2 * Sf[i + 1] + Sf[i]
        if abs(denom) < 1e-300:
            continue
        ai = Sf[i] - (Sf[i + 1] - Sf[i]) ** 2 / denom
        aitken.append(ai)

    # Direct geometric: assume S_k = L + c*r^k; estimate r from successive
    # differences and L from limit
    diffs = np.diff(Sf)
    if len(diffs) >= 2 and abs(diffs[0]) > 1e-300:
        ratios = diffs[1:] / diffs[:-1]
        r_est = float(np.mean(ratios))
        # If fit is good: S_inf = S_{k_max} + diffs[-1] * r / (1 - r)
        if abs(1 - r_est) > 1e-8:
            geom_est = Sf[-1] + diffs[-1] * r_est / (1 - r_est)
        else:
            geom_est = float("nan")
    else:
        r_est = float("nan")
        geom_est = float("nan")

    # Best estimate: average of Aitken-last and geometric, if both available
    candidates = []
    if aitken:
        candidates.append(("aitken_last", aitken[-1]))
    if not np.isnan(geom_est):
        candidates.append(("geometric", geom_est))
    candidates.append(("last_value", Sf[-1]))

    if len(candidates) >= 2:
        vals = [v for _, v in candidates]
        best = float(np.mean(vals))
        uncert = float(np.std(vals))
    else:
        best = candidates[0][1]
        uncert = float("nan")

    return {
        "S_inf_est": best,
        "method": "; ".join(c[0] for c in candidates),
        "uncertainty": uncert,
        "r_est": r_est,
        "aitken_history": aitken,
        "geom_est": geom_est,
        "candidates": candidates,
    }


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #

def main():
    print("=" * 78)
    print("Q-sweep test 2: c_q = S_inf^(q) / q for q in {3, 5, 7, 11, 13}")
    print("=" * 78)
    print()

    cache = load_cache()
    print(f"Cache: {len(cache)} entries already computed.")
    if cache:
        for key in sorted(cache.keys()):
            q, k = key
            Xk, n = cache[key]
            print(f"  (q={q}, k={k}): {n} states, X_{k} ~ {float(Xk):.6f}")
    print()

    # Resource plan
    plan = {
        3: 4,
        5: 4,
        7: 3,
        11: 3,
        13: 3,
    }

    # Order: cheap first (q=3, q=7) then heavier (q=5, q=11, q=13).
    # q=3 verifies path; q=7 confirms generalization; q=5 fills out the family;
    # q=11, q=13 extend the trend.
    run_order = [3, 7, 5, 11, 13]

    all_S = {}
    for q in run_order:
        k_max = plan[q]
        print(f"\n{'#' * 78}\n# q = {q}  (k_max = {k_max})\n{'#' * 78}")
        Ss = compute_S_sequence(q, k_max, cache)
        all_S[q] = Ss
        print()

    # Stage 3 extrapolation
    print("=" * 78)
    print("Stage 3: extrapolation of S_inf^(q)")
    print("=" * 78)
    print()

    extrap = {}
    for q in run_order:
        Ss = all_S[q]
        ext = extrapolate_S_inf(Ss, q)
        extrap[q] = ext
        print(f"q = {q}:")
        print(f"  S_k sequence:        {[f'{float(S):.6f}' for _, S in Ss]}")
        print(f"  S_inf estimate:      {ext['S_inf_est']:.10f}")
        print(f"  uncertainty:         {ext.get('uncertainty', float('nan'))}")
        print(f"  r_est (decay ratio): {ext.get('r_est', float('nan'))}")
        print(f"  candidates:          {ext.get('candidates', [])}")
        print()

    # Stage 4: c_q = S_inf^(q) / q
    print("=" * 78)
    print("Stage 4: c_q = S_inf^(q) / q")
    print("=" * 78)
    print()
    print(f"  {'q':>3}  {'S_inf est':>15}  {'c_q est':>15}  {'c_q float':>15}")
    rows_csv = []
    for q in [3, 5, 7, 11, 13]:
        if q not in extrap:
            continue
        S_inf = extrap[q]["S_inf_est"]
        c_q = S_inf / q
        print(f"  {q:>3}  {S_inf:>15.10f}  {c_q:>15.10f}  {c_q:>15.6e}")
        rows_csv.append({
            "q": q,
            "S_inf_estimate": S_inf,
            "c_q_estimate": c_q,
            "k_max": plan[q],
            "uncertainty": extrap[q].get("uncertainty", float("nan")),
            "r_est": extrap[q].get("r_est", float("nan")),
        })
        # Save the raw S_k sequence too
        for k, S in all_S[q]:
            rows_csv.append({
                "q": q,
                "k": k,
                "S_k_num": S.numerator,
                "S_k_den": S.denominator,
                "S_k_float": float(S),
            })

    # Reference: q=3 ground truth
    print()
    print(f"  q=3 reference: S_inf = 7/15 = {7/15:.10f}; c_3 = 7/45 = {7/45:.10f}")
    if 3 in extrap:
        print(f"  q=3 estimated S_inf = {extrap[3]['S_inf_est']:.10f}  "
              f"(error from 7/15 = {extrap[3]['S_inf_est'] - 7/15:+.6e})")

    # CSV write (heterogeneous columns; use full union)
    fields = ["q", "k", "S_k_num", "S_k_den", "S_k_float",
              "S_inf_estimate", "c_q_estimate", "k_max", "uncertainty", "r_est"]
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows_csv)
    print(f"\n[csv: {OUT_CSV}]")

    # Stage 5: closed-form test (basic)
    print()
    print("=" * 78)
    print("Stage 5: closed-form test")
    print("=" * 78)
    print()
    qs_with_data = [q for q in [3, 5, 7, 11, 13] if q in extrap]
    cqs = [extrap[q]["S_inf_est"] / q for q in qs_with_data]
    print(f"  (q, c_q) data:")
    for q, c in zip(qs_with_data, cqs):
        print(f"    q={q:>3}  c_q = {c:.10f}")

    # Try simple rational families:
    #   c_q = a / (b * q),   c_q = a / (b * q^2),  c_q = (q-1)/(b*q^2),
    #   c_q = a / (q^2 - 1), c_q = (q-1)/(b*q*(q+1))
    # Fit by linear regression in transformed coordinates; check R^2.
    print()
    print("  Closed-form candidate fits:")
    qs_arr = np.array(qs_with_data, dtype=float)
    cqs_arr = np.array(cqs)

    candidates = [
        ("c_q = A / (B * q)",        lambda q: 1.0 / q),
        ("c_q = A / (B * q^2)",      lambda q: 1.0 / q ** 2),
        ("c_q = A / (q^2 - 1)",      lambda q: 1.0 / (q * q - 1)),
        ("c_q = A * (q-1)/(B * q^2)", lambda q: (q - 1) / q ** 2),
        ("c_q = A / (B * q * (q+1))", lambda q: 1.0 / (q * (q + 1))),
        ("c_q = (q-1) / (B * q * (q+1))", lambda q: (q - 1) / (q * (q + 1))),
        ("c_q = (q-1) / (B * q^3)",  lambda q: (q - 1) / q ** 3),
    ]
    for name, f in candidates:
        vals = np.array([f(q) for q in qs_arr])
        # Fit cqs ≈ A * vals; A = mean(cqs / vals); residual is std
        if np.any(vals == 0):
            continue
        ratios = cqs_arr / vals
        A = float(np.mean(ratios))
        std = float(np.std(ratios))
        rel_std = std / A if A != 0 else float("nan")
        print(f"    {name}:  A = {A:.6f},  rel-std of A across q = {rel_std:.4e}")

    # Try expressing c_q rationally over Q assuming exact extrapolation; this
    # works only at high precision. Print decimal forms with q=3 anchor 7/45
    # for cross-check.
    print()
    print(f"  q=3 reference: c_3 = 7/45 = {7/45:.10f}; "
          f"estimated c_3 = {cqs[qs_with_data.index(3)]:.10f}")


if __name__ == "__main__":
    main()
