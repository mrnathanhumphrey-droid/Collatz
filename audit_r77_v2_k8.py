"""
audit_r77_v2_k8.py — diagnostic audit harness for the k=8 R77.7 v2 attempt.

Purpose
-------
The k=8 run died at 2.8 GB memory with empty stdout (buffer deadlock).
This audit runs ONE prime end-to-end at k=8 with full instrumentation
so we can decide whether the prior failure was:

  (a) Real memory pressure → matrix grows too big, needs architecture upgrade
      (parallel primes / GPU / sparse / Wiedemann)
  (b) Stdout-buffer-deadlock hiding progress → script was working, just
      invisible; fix logging and run serial as-is

What it measures
----------------
For ONE prime at k=8:
  - Phase A: build_K_exponents (one-time setup, but measures baseline)
  - Phase B: build pow2_mod table (cheap, sanity)
  - Phase C: solve_pi_mod_p — split into:
      C1: build A_p from K_exp + pow2_mod (pure-Python triple loop, expected bottleneck)
      C2: gauss_solve_mod_p (numpy dense Gauss)
  - Phase D: extrapolate full-run estimate (500 primes × per-prime cost)

For each phase:
  - Wall-clock time
  - Peak RSS (resident set size) via psutil
  - Python allocations via tracemalloc snapshot
  - Matrix dimensions / dict sizes / sanity numbers

Output
------
  C:/Collatz/experiments_output/audit_r77_v2_k8.log  (line-by-line log,
                                                      flushed after every write)
  C:/Collatz/experiments_output/audit_r77_v2_k8.json (structured summary)

Run
---
  python C:/Collatz/audit_r77_v2_k8.py

  Or with explicit log path:
  python C:/Collatz/audit_r77_v2_k8.py --log path/to/audit.log

Expected wall-clock for the audit itself: 4-10 minutes at k=8 single-prime.
Memory ceiling expected: ~3 GB. If audit itself dies, log up to the death
point still tells us where the failure is.

Decision rule
-------------
After audit completes (or dies):
  - If Phase A peak ≈ 1.5-2 GB (K_exp baseline) and Phase C1 dominates wall-clock
    → case (b): fix logging, vectorize Phase C1, run serial.
  - If Phase A peak > 3 GB → case (a): K_exp too big for single-process,
    need parallel primes with shared K_exp (or rebuild per-worker).
  - If Phase C1 OK but Phase C2 explodes memory → temporary arrays in Gauss,
    rewrite with in-place ops.
  - If audit completes < 5 min single-prime → the prior k=8 run was just
    slow, not broken. Fix logging and run 500 primes (~24 hr wall).

This audit does NOT modify the solver. Pure observation.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import sys
import time
import tracemalloc
from pathlib import Path

import numpy as np
import psutil

# Reuse v2 solver internals
sys.path.insert(0, r"C:\Collatz")
from result_77_7_v2 import (
    build_K_exponents,
    solve_pi_mod_p,
    is_prime,
)


# ---------------------------------------------------------------------------
# Logging — write directly to file, flush after every line. NEVER use stdout
# piped to tee (that's exactly what caused the buffer deadlock).
# ---------------------------------------------------------------------------

class FlushLogger:
    def __init__(self, path: Path):
        self.path = path
        self.fh = open(path, "w", encoding="utf-8", buffering=1)  # line-buffered

    def log(self, msg: str = ""):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        self.fh.write(line + "\n")
        self.fh.flush()
        os.fsync(self.fh.fileno())  # force OS flush too
        # Mirror to stdout for live tail (but write to file is authoritative)
        print(line, flush=True)

    def close(self):
        self.fh.close()


def get_rss_mb() -> float:
    """Current process RSS in MB."""
    return psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)


def get_tracemalloc_mb() -> float:
    """Total tracked Python allocations in MB."""
    if tracemalloc.is_tracing():
        current, _peak = tracemalloc.get_traced_memory()
        return current / (1024 * 1024)
    return 0.0


# ---------------------------------------------------------------------------
# Audit main
# ---------------------------------------------------------------------------

def audit_k8(log: FlushLogger, k: int = 8, prime: int | None = None):
    """Run one full single-prime pipeline at k with heavy instrumentation."""

    log.log("=" * 70)
    log.log(f"R77.7 v2 audit — single-prime end-to-end at k={k}")
    log.log("=" * 70)
    log.log(f"PID = {os.getpid()}")
    log.log(f"Python {sys.version.split()[0]}, numpy {np.__version__}")

    # System info
    vm = psutil.virtual_memory()
    log.log(f"System RAM: {vm.total / 1024**3:.1f} GB total, "
            f"{vm.available / 1024**3:.1f} GB available")
    log.log(f"Initial RSS: {get_rss_mb():.1f} MB")
    log.log("")

    tracemalloc.start()

    summary = {
        "k": k,
        "phases": {},
    }

    # =====================================================================
    # Phase A: build_K_exponents
    # =====================================================================
    log.log("-" * 70)
    log.log("Phase A: build_K_exponents")
    log.log("-" * 70)
    t0 = time.time()
    rss0 = get_rss_mb()
    tm0 = get_tracemalloc_mb()

    K_exp, M_param, coprime_states, N = build_K_exponents(k)

    t_A = time.time() - t0
    rss1 = get_rss_mb()
    tm1 = get_tracemalloc_mb()

    # K_exp statistics
    total_entries = sum(len(row) for row in K_exp)
    total_exponent_count = sum(sum(len(e) for e in row.values()) for row in K_exp)
    avg_row_density = total_entries / N
    avg_exponents_per_entry = total_exponent_count / max(total_entries, 1)

    log.log(f"  wall time: {t_A:.2f} s")
    log.log(f"  RSS delta: {rss1 - rss0:+.1f} MB (now {rss1:.1f} MB)")
    log.log(f"  Python alloc delta: {tm1 - tm0:+.1f} MB")
    log.log(f"  N (matrix dim) = {N}")
    log.log(f"  M_param (geom truncation) = {M_param}")
    log.log(f"  K_exp: {total_entries:,} total (i,j) entries, "
            f"{total_exponent_count:,} total exponents")
    log.log(f"  avg entries per row: {avg_row_density:.1f}")
    log.log(f"  avg exponents per (i,j): {avg_exponents_per_entry:.2f}")
    log.log("")

    summary["phases"]["A_build_K_exp"] = {
        "wall_s": t_A,
        "rss_delta_mb": rss1 - rss0,
        "rss_after_mb": rss1,
        "tracemalloc_delta_mb": tm1 - tm0,
        "N": N,
        "M_param": M_param,
        "K_exp_total_entries": total_entries,
        "K_exp_total_exponents": total_exponent_count,
        "K_exp_avg_row_density": avg_row_density,
        "K_exp_avg_exponents_per_entry": avg_exponents_per_entry,
    }

    # =====================================================================
    # Phase B: pow2_mod table
    # =====================================================================
    if prime is None:
        # Pick a known-safe prime near the start of v2's range
        prime = 1_000_000_007
        while not is_prime(prime):
            prime += 2

    log.log("-" * 70)
    log.log(f"Phase B: pow2_mod table for prime p = {prime}")
    log.log("-" * 70)
    t0 = time.time()
    rss0 = get_rss_mb()

    pow2_mod = [pow(2, e, prime) for e in range(M_param + 1)]

    t_B = time.time() - t0
    rss1 = get_rss_mb()
    log.log(f"  wall time: {t_B:.2f} s")
    log.log(f"  RSS delta: {rss1 - rss0:+.1f} MB (now {rss1:.1f} MB)")
    log.log(f"  table size: {len(pow2_mod)}")
    log.log("")

    summary["phases"]["B_pow2_mod"] = {
        "wall_s": t_B,
        "rss_delta_mb": rss1 - rss0,
        "rss_after_mb": rss1,
        "table_size": len(pow2_mod),
    }

    # =====================================================================
    # Phase C: solve_pi_mod_p — manually split into C1 (build A_p) and C2 (Gauss)
    # =====================================================================
    log.log("-" * 70)
    log.log("Phase C: solve_pi_mod_p (split C1 build A_p, C2 Gauss)")
    log.log("-" * 70)

    # C1: build A_p (the pure-Python triple-loop that's the suspected bottleneck)
    t0 = time.time()
    rss0 = get_rss_mb()
    tm0 = get_tracemalloc_mb()

    D_mod = (pow2_mod[M_param] - 1) % prime

    A_p = np.zeros((N, N), dtype=np.int64)
    for i in range(N):
        row_dict = K_exp[i]
        for j, exponents in row_dict.items():
            w_mod = 0
            for e in exponents:
                w_mod += pow2_mod[e]
            A_p[j, i] = (A_p[j, i] + w_mod) % prime

    for i in range(N):
        A_p[i, i] = (A_p[i, i] - D_mod) % prime

    A_p[N - 1, :] = 1
    b_p = np.zeros(N, dtype=np.int64)
    b_p[N - 1] = 1

    t_C1 = time.time() - t0
    rss1 = get_rss_mb()
    tm1 = get_tracemalloc_mb()

    A_p_size_mb = A_p.nbytes / (1024 * 1024)
    log.log(f"  C1 build A_p:")
    log.log(f"    wall time: {t_C1:.2f} s")
    log.log(f"    RSS delta: {rss1 - rss0:+.1f} MB (now {rss1:.1f} MB)")
    log.log(f"    Python alloc delta: {tm1 - tm0:+.1f} MB")
    log.log(f"    A_p dense N×N int64: {A_p_size_mb:.1f} MB")
    log.log("")

    summary["phases"]["C1_build_A_p"] = {
        "wall_s": t_C1,
        "rss_delta_mb": rss1 - rss0,
        "rss_after_mb": rss1,
        "tracemalloc_delta_mb": tm1 - tm0,
        "A_p_size_mb": A_p_size_mb,
    }

    # C2: gauss_solve_mod_p
    from result_77_7_v2 import gauss_solve_mod_p
    t0 = time.time()
    rss0 = get_rss_mb()
    tm0 = get_tracemalloc_mb()
    rss_peak_during_gauss = rss0

    # Sample RSS during a 5-second window mid-Gauss to catch the temporary blowup
    # (we can't sample concurrently from same thread, but the peak should reveal
    # in the post-Gauss RSS minus the working-set release)
    x_p = gauss_solve_mod_p(A_p, b_p, prime)

    t_C2 = time.time() - t0
    rss_after_gauss = get_rss_mb()
    tm1 = get_tracemalloc_mb()

    log.log(f"  C2 gauss_solve_mod_p:")
    log.log(f"    wall time: {t_C2:.2f} s")
    log.log(f"    RSS delta from C2 start: {rss_after_gauss - rss0:+.1f} MB "
            f"(now {rss_after_gauss:.1f} MB)")
    log.log(f"    Python alloc delta: {tm1 - tm0:+.1f} MB")
    log.log(f"    solution: {'OK' if x_p is not None else 'SINGULAR'}")
    log.log("")

    summary["phases"]["C2_gauss"] = {
        "wall_s": t_C2,
        "rss_delta_mb": rss_after_gauss - rss0,
        "rss_after_mb": rss_after_gauss,
        "tracemalloc_delta_mb": tm1 - tm0,
        "singular": x_p is None,
    }

    # =====================================================================
    # Phase D: extrapolate
    # =====================================================================
    log.log("-" * 70)
    log.log("Phase D: full-run extrapolation")
    log.log("-" * 70)

    n_primes_estimate = 500
    per_prime_s = t_C1 + t_C2  # K_exp built once, pow2_mod is cheap
    total_s = per_prime_s * n_primes_estimate
    log.log(f"  Per-prime cost: C1 {t_C1:.1f}s + C2 {t_C2:.1f}s = {per_prime_s:.1f}s")
    log.log(f"  Estimated full k={k} run @ {n_primes_estimate} primes:")
    log.log(f"    total wall: {total_s:.0f} s = {total_s/3600:.1f} hr")
    log.log("")
    log.log(f"  Memory ceiling observed: {max(rss_after_gauss, rss1):.1f} MB")
    log.log(f"  Headroom on this machine: ~{vm.total / 1024**2 - rss_after_gauss:.0f} MB free")
    log.log("")

    summary["phases"]["D_extrapolation"] = {
        "per_prime_wall_s": per_prime_s,
        "n_primes_estimate": n_primes_estimate,
        "total_wall_s": total_s,
        "total_wall_hr": total_s / 3600,
        "memory_ceiling_mb": max(rss_after_gauss, rss1),
        "system_ram_total_gb": vm.total / 1024**3,
    }

    # =====================================================================
    # Tracemalloc top allocations (where the memory actually went)
    # =====================================================================
    snapshot = tracemalloc.take_snapshot()
    top_stats = snapshot.statistics("filename")[:10]

    log.log("-" * 70)
    log.log("Top Python allocations by file:")
    log.log("-" * 70)
    for stat in top_stats:
        log.log(f"  {stat.size / 1024 / 1024:7.1f} MB  {stat.traceback.format()[0]}")
    log.log("")

    summary["top_allocations"] = [
        {
            "file": stat.traceback.format()[0],
            "size_mb": stat.size / 1024 / 1024,
        }
        for stat in top_stats
    ]

    # =====================================================================
    # Diagnostic verdict
    # =====================================================================
    log.log("=" * 70)
    log.log("VERDICT")
    log.log("=" * 70)

    if t_C1 > 5 * t_C2:
        log.log("  C1 (build A_p) is the dominant phase (>5x C2 Gauss).")
        log.log("  Fix: vectorize the pure-Python triple loop with numpy assemble.")
    elif t_C2 > 2 * t_C1:
        log.log("  C2 (Gauss) is the dominant phase. Architecture upgrade more "
                "useful than vectorizing C1.")
    else:
        log.log("  C1 and C2 comparable. Both vectorizing + parallelism help.")
    log.log("")

    if max(rss_after_gauss, rss1) > 3 * 1024:
        log.log("  Memory ceiling >3 GB: parallel primes constrained, "
                "GPU or sparse becomes load-bearing for k≥9.")
    else:
        log.log(f"  Memory ceiling {max(rss_after_gauss, rss1):.0f} MB: "
                f"parallel primes feasible at this k.")
    log.log("")

    if total_s / 3600 < 4:
        log.log(f"  Full-run estimate {total_s/3600:.1f}hr: serial is acceptable. "
                f"Fix logging and re-run as-is.")
    elif total_s / 3600 < 24:
        log.log(f"  Full-run estimate {total_s/3600:.1f}hr: parallel primes "
                f"reduces to manageable.")
    else:
        log.log(f"  Full-run estimate {total_s/3600:.1f}hr: architecture upgrade "
                f"(parallel + vectorized C1 + possibly GPU) needed.")

    summary["verdict"] = {
        "dominant_phase": "C1" if t_C1 > 5 * t_C2 else ("C2" if t_C2 > 2 * t_C1 else "balanced"),
        "memory_ceiling_gb": max(rss_after_gauss, rss1) / 1024,
        "full_run_estimate_hr": total_s / 3600,
    }

    tracemalloc.stop()
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=8, help="Markov chain level (default 8)")
    parser.add_argument("--prime", type=int, default=None, help="Test prime (default ~1e9)")
    parser.add_argument(
        "--log",
        type=str,
        default=r"C:\Collatz\experiments_output\audit_r77_v2_k8.log",
    )
    parser.add_argument(
        "--json",
        type=str,
        default=r"C:\Collatz\experiments_output\audit_r77_v2_k8.json",
    )
    args = parser.parse_args()

    log_path = Path(args.log)
    json_path = Path(args.json)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    log = FlushLogger(log_path)
    try:
        summary = audit_k8(log, k=args.k, prime=args.prime)
        with open(json_path, "w", encoding="utf-8") as fh:
            json.dump(summary, fh, indent=2, default=str)
        log.log(f"Audit summary JSON: {json_path}")
        log.log("Audit complete.")
    except Exception as e:
        log.log(f"AUDIT CRASHED: {type(e).__name__}: {e}")
        import traceback
        log.log(traceback.format_exc())
        raise
    finally:
        log.close()


if __name__ == "__main__":
    main()
