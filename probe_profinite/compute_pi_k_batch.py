"""
compute_pi_k_batch.py
=====================
Recompute π_k stationary measure for k=8, 9, 10, 11 and save as .npz files
into probe_profinite/. Uses the same matrix-free MatVecK power iteration
as compute_epsilon_12.py.

π_5, π_6, π_7 already cached at C:/Collatz/probe_mode_amplitudes/pi_k{N}.npy
π_12 already cached at C:/Collatz/probe_epsilon_12/pi_12.npz

Usage: python compute_pi_k_batch.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

# Reuse MatVecK from probe_epsilon_12
sys.path.insert(0, r"C:\Collatz\probe_epsilon_12")
from compute_epsilon_12 import MatVecK, power_iteration  # noqa

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = Path(r"C:\Collatz\probe_profinite")
OUT_DIR.mkdir(exist_ok=True)


def compute_and_save(k, chunk):
    print(f"\n=== k = {k}, chunk = {chunk} ===", flush=True)
    out_path = OUT_DIR / f"pi_{k}.npz"
    if out_path.exists():
        print(f"  exists at {out_path}, skipping")
        return
    t0 = time.time()
    K_op = MatVecK(3, k, chunk=chunk)
    print(f"  init: n = {K_op.n:,}, M = {K_op.M:,}, "
          f"{time.time()-t0:.2f}s")
    t0 = time.time()
    pi_k, iters, res, _ = power_iteration(
        K_op, max_iter=200, tol=1e-13, aitken=False, verbose_every=10
    )
    print(f"  converged: {iters} iters, residual {res:.2e}, "
          f"{time.time()-t0:.1f}s")
    np.savez_compressed(out_path, pi=pi_k, coprime=K_op.coprime, k=k)
    print(f"  saved {out_path}")


def main():
    # Order: cheap first, k=11 is the slow one (~24 min).
    for k, chunk in [(8, 512), (9, 512), (10, 512), (11, 1024)]:
        compute_and_save(k, chunk)
    print("\nAll π_k computed.")


if __name__ == "__main__":
    main()
