"""
q-sweep: empirical K-sum |K_p(r, c, m)| across primes p in {3, 5, 7, 11, 13}.

  K_p(r, c, m) = sum_{u=0}^{N_p - 1} e_M( c * (1+p)^u  -  p^2 * m * u )
  M = p^(r+1),  N_p = p^(r-1) = M / p^2,  e_M(x) = exp(2*pi*i*x/M)

For p=3 this reduces to the existing K(r, c, m) at r79b_S_partial_data.csv.
For p>=5 we compute fresh, append-with-fsync to C:/Collatz/q_sweep_data.csv.

Usage:
  python q_sweep.py [--primes 5 7 11 13] [--reload-p3]

Output schema (q_sweep_data.csv):
  p, r, M, N, K_c1m0_abs, K_max_abs, K_max_c, K_max_m,
  log_K_max, log_sqrt_M, ratio, log_ratio, beta_running, elapsed_s

Note on field semantics (matches the brief):
  ratio       = K_max_abs / sqrt(N)            ("prefactor" C)
  log_ratio   = log(K_max_abs) / log(N)        ("saturation exponent" beta, point-est)
  beta_running= OLS slope of log(K_max) vs log(N) for this p over rows so far
"""
import os
import sys
import csv
import time
import math
import argparse
import numpy as np
from numba import njit, prange


sys.stdout.reconfigure(encoding="utf-8")


# ----------------------------------------------------------------------
# Direct K-sum for general prime p (Numba-parallel)
# ----------------------------------------------------------------------

@njit(cache=True, parallel=True, fastmath=False)
def K_direct(p: int, r: int, c: int, m: int, n_blocks: int) -> complex:
    """K_p(r, c, m) = sum_{u=0}^{N-1} e_M(c * (1+p)^u - p^2 * m * u),
    where M = p^(r+1) and N = p^(r-1).

    Block-parallel: blocks share precomputed (1+p)^(b*block_size) mod M starts.
    """
    # Compute M and N as int64 (will fit for our targeted sizes).
    M = 1
    for _ in range(r + 1):
        M = M * p
    N = M // (p * p)
    base = 1 + p
    p2_m = (p * p) * m
    inv_M = 2.0 * math.pi / M

    block_size = (N + n_blocks - 1) // n_blocks

    # Precompute starts: (1+p)^(b*block_size) mod M
    starts = np.empty(n_blocks, dtype=np.int64)
    cur = 1
    for b in range(n_blocks):
        starts[b] = cur
        for _ in range(block_size):
            cur = (cur * base) % M

    block_re = np.zeros(n_blocks, dtype=np.float64)
    block_im = np.zeros(n_blocks, dtype=np.float64)

    for b in prange(n_blocks):
        u_start = b * block_size
        u_end = u_start + block_size
        if u_end > N:
            u_end = N
        if u_start >= N:
            continue
        x = starts[b]
        sre = 0.0
        sim = 0.0
        for u in range(u_start, u_end):
            phase = (c * x - p2_m * u) % M
            angle = inv_M * phase
            sre += math.cos(angle)
            sim += math.sin(angle)
            x = (x * base) % M
        block_re[b] = sre
        block_im[b] = sim

    return complex(block_re.sum(), block_im.sum())


@njit(cache=True)
def K_serial(p: int, r: int, c: int, m: int) -> complex:
    """Reference single-threaded version for verification."""
    M = 1
    for _ in range(r + 1):
        M = M * p
    N = M // (p * p)
    base = 1 + p
    p2_m = (p * p) * m
    inv_M = 2.0 * math.pi / M
    sre = 0.0
    sim = 0.0
    x = 1
    for u in range(N):
        phase = (c * x - p2_m * u) % M
        angle = inv_M * phase
        sre += math.cos(angle)
        sim += math.sin(angle)
        x = (x * base) % M
    return complex(sre, sim)


# ----------------------------------------------------------------------
# Per-(p, r) sweep: K_max over (c, m) in {1..p-1} x {0..p-1}
# ----------------------------------------------------------------------

def compute_K_max_for_p_r(p: int, r: int, n_blocks: int):
    """Returns (K_c1m0_abs, K_max_abs, K_max_c, K_max_m, elapsed_s)."""
    t0 = time.time()
    K_c1m0 = K_direct(p, r, 1, 0, n_blocks)
    K_c1m0_abs = abs(K_c1m0)

    K_max_abs = K_c1m0_abs
    K_max_c = 1
    K_max_m = 0
    for c in range(1, p):
        for m in range(0, p):
            if c == 1 and m == 0:
                continue
            K = K_direct(p, r, c, m, n_blocks)
            a = abs(K)
            if a > K_max_abs:
                K_max_abs = a
                K_max_c = c
                K_max_m = m
    elapsed = time.time() - t0
    return K_c1m0_abs, K_max_abs, K_max_c, K_max_m, elapsed


# ----------------------------------------------------------------------
# Verify p=3 implementation matches existing CSV
# ----------------------------------------------------------------------

def verify_p3_match(existing_csv: str, n_blocks: int) -> bool:
    """Re-compute K for p=3 at r=8 and a few small r; check vs existing CSV
    K_c1m0_abs to 1e-3 absolute tolerance."""
    print("[verify] re-running p=3 at r=8,9 to cross-check existing CSV...")
    targets = {}
    with open(existing_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            r = int(row["r"])
            if r in (8, 9, 10):
                targets[r] = float(row["K_c1m0_abs"])
    ok = True
    for r in sorted(targets):
        K = K_direct(3, r, 1, 0, n_blocks)
        a = abs(K)
        ref = targets[r]
        rel = abs(a - ref) / max(ref, 1.0)
        status = "OK" if rel < 1e-4 else "MISMATCH"
        print(f"  p=3 r={r}: |K|_c1m0 = {a:.6f}  ref={ref:.6f}  rel_err={rel:.2e}  [{status}]")
        if rel > 1e-4:
            ok = False
    return ok


# ----------------------------------------------------------------------
# Load p=3 rows from existing CSV (every 2nd r per brief: 8,10,12,14,16,18,20)
# ----------------------------------------------------------------------

def load_p3_rows(existing_csv: str, r_targets):
    rows = []
    with open(existing_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            r = int(row["r"])
            if r in r_targets:
                rows.append({
                    "p": 3,
                    "r": r,
                    "M": int(row["q"]),
                    "N": int(row["N"]),
                    "K_c1m0_abs": float(row["K_c1m0_abs"]),
                    "K_max_abs": float(row["K_max_abs"]),
                    "K_max_c": int(row["K_max_c"]),
                    "K_max_m": int(row["K_max_m"]),
                    "elapsed_s": float(row["elapsed_s"]),
                })
    return rows


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

OUTPUT_CSV = r"C:\Collatz\q_sweep_data.csv"
P3_CSV     = r"C:\Collatz\r79b_S_partial_data.csv"

PRIME_R_TARGETS = {
    3: [8, 10, 12, 14, 16, 18, 20],
    5: [6, 8, 10, 12],
    7: [6, 8, 10],
    11: [5, 6, 7, 8],
    13: [4, 5, 6, 7],
}


def append_row(csv_path: str, fieldnames, row, write_header: bool):
    with open(csv_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(row)
        f.flush()
        os.fsync(f.fileno())


def beta_running_for_p(rows_for_p):
    """OLS slope of log K_max_abs on log N over rows_for_p so far."""
    if len(rows_for_p) < 2:
        return float("nan")
    xs = np.array([math.log(r["N"]) for r in rows_for_p], dtype=np.float64)
    ys = np.array([math.log(r["K_max_abs"]) for r in rows_for_p], dtype=np.float64)
    if len(xs) < 2:
        return float("nan")
    slope, _ = np.polyfit(xs, ys, 1)
    return float(slope)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", type=int, nargs="+", default=[3, 5, 7, 11, 13])
    parser.add_argument("--n-blocks", type=int, default=32)
    parser.add_argument("--fresh", action="store_true", help="overwrite output CSV")
    args = parser.parse_args()

    fieldnames = [
        "p", "r", "M", "N",
        "K_c1m0_abs", "K_max_abs", "K_max_c", "K_max_m",
        "log_K_max", "log_sqrt_M", "ratio", "log_ratio",
        "beta_running", "elapsed_s",
    ]

    if args.fresh and os.path.exists(OUTPUT_CSV):
        os.remove(OUTPUT_CSV)

    write_header = not os.path.exists(OUTPUT_CSV)
    if write_header:
        with open(OUTPUT_CSV, "a", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            f.flush()
            os.fsync(f.fileno())

    # Warm up Numba kernels
    print("[warmup] compiling Numba kernels...")
    t0 = time.time()
    _ = K_direct(3, 3, 1, 0, args.n_blocks)
    _ = K_direct(5, 3, 1, 0, args.n_blocks)
    _ = K_serial(3, 3, 1, 0)
    print(f"[warmup] done in {time.time()-t0:.1f}s")

    # Verify p=3 implementation if p=3 in primes
    if 3 in args.primes:
        ok = verify_p3_match(P3_CSV, args.n_blocks)
        if not ok:
            print("[verify] WARNING: p=3 mismatch — formulas may differ. Continuing anyway.")

    rows_by_p = {}

    for p in args.primes:
        rows_by_p[p] = []
        r_list = PRIME_R_TARGETS.get(p, [])
        print(f"\n=== p = {p}, r in {r_list} ===")
        for r in r_list:
            M = p ** (r + 1)
            N = p ** (r - 1)

            if p == 3:
                # Reuse existing CSV for p=3
                p3_rows = load_p3_rows(P3_CSV, {r})
                if not p3_rows:
                    print(f"  p=3 r={r}: NOT FOUND in existing CSV — skipping")
                    continue
                src = p3_rows[0]
                K_c1m0_abs = src["K_c1m0_abs"]
                K_max_abs = src["K_max_abs"]
                K_max_c = src["K_max_c"]
                K_max_m = src["K_max_m"]
                elapsed = src["elapsed_s"]
                print(f"  p=3 r={r}: REUSED |K|_c1m0={K_c1m0_abs:.3f}  |K|_max={K_max_abs:.3f}")
            else:
                # Compute fresh
                K_c1m0_abs, K_max_abs, K_max_c, K_max_m, elapsed = \
                    compute_K_max_for_p_r(p, r, args.n_blocks)
                print(f"  p={p} r={r}: M={M:.3e} N={N:.3e}  |K|_c1m0={K_c1m0_abs:.3f}  "
                      f"|K|_max={K_max_abs:.3f} (c={K_max_c}, m={K_max_m})  [{elapsed:.1f}s]")

            log_K_max = math.log(K_max_abs) if K_max_abs > 0 else float("-inf")
            log_sqrt_M = 0.5 * math.log(M)
            ratio = K_max_abs / math.sqrt(N) if N > 0 else float("nan")
            log_ratio = log_K_max / math.log(N) if N > 1 else float("nan")

            row_dict = {
                "p": p, "r": r, "M": M, "N": N,
                "K_c1m0_abs": K_c1m0_abs,
                "K_max_abs": K_max_abs,
                "K_max_c": K_max_c, "K_max_m": K_max_m,
            }
            rows_by_p[p].append(row_dict)

            beta = beta_running_for_p(rows_by_p[p])

            out_row = {
                "p": p, "r": r, "M": M, "N": N,
                "K_c1m0_abs": f"{K_c1m0_abs:.6f}",
                "K_max_abs": f"{K_max_abs:.6f}",
                "K_max_c": K_max_c, "K_max_m": K_max_m,
                "log_K_max": f"{log_K_max:.6f}",
                "log_sqrt_M": f"{log_sqrt_M:.6f}",
                "ratio": f"{ratio:.6f}",
                "log_ratio": f"{log_ratio:.6f}",
                "beta_running": f"{beta:.6f}" if not math.isnan(beta) else "",
                "elapsed_s": f"{elapsed:.3f}",
            }
            append_row(OUTPUT_CSV, fieldnames, out_row, write_header=False)

    print(f"\n[done] csv = {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
