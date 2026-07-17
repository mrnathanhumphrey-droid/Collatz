"""
fhat_verification.py — adversarial verification of the F̂_p theorem candidate.

CORRECTED METHOD (A4 deviation from pre-reg; see RESULTS doc §Pre-reg deviation log):

Pre-reg called for zero-padded length-M FFT of f_p, claiming F̂_short^padded vanishes
off the principal-unit sub-support. This is WRONG: zero-padding produces spectral
leakage. The theorem candidate is about F̂_full (periodic-extension of f_p to
length M), which DOES vanish off the predicted support.

Cleanest computation: length-period DFT G[a] := Σ_{s=0}^{p^r - 1} f_p(s) e^{-2πi a s / p^r}
for a ∈ Z/p^r. Algebraically:
    F̂_full(p·a) = p · G[a]   (from u = q·p^r + s splitting; q-sum = p when a integer)
    F̂_full(ξ) = 0    for ξ ∉ p·Z/M (orthogonality on cosets of p in Z/M)

Theorem candidate in this representation:
    |G[a]| = p^{(r+1)/2}    iff a ≡ 1 (mod p)
    |G[a]| = 0              otherwise
equivalently  |F̂_full(p·a)| = p^{(r+3)/2}  for a ≡ 1 (mod p)  and  0 otherwise.

ON-support magnitudes are unchanged from the prior script (FFT artifacts only affect
the zero-padded version's OFF-support values).

Phases (per pre-reg FHAT_THEOREM_VERIFICATION_PRE_REGISTRATION.md):
    Phase 1: p ∈ {11, 13, 17, 19, 23, 29, 31} × r ∈ {2, 3}
    Phase 2: p ∈ {3, 5, 7} × r ∈ {4, 5, 6}
    Phase 4: boundaries (p=2, r ∈ {2,3,4}; r=1, p ∈ {3,5,7,11})
    Phase 5: support characterization (uses Phase 1/2 G[a] outputs)
"""
from __future__ import annotations

import math
import csv
import time
from pathlib import Path

import numpy as np


OUTPATH = Path("C:/Collatz/fhat_verification_results.csv")


def f_p_period(p: int, r: int, c: int = 1) -> np.ndarray:
    """Return f_p(u) = exp(2πi · c · (1+p)^u / M) for u = 0..p^r - 1, complex128."""
    M = p ** (r + 1)
    period = p ** r
    out = np.empty(period, dtype=np.complex128)
    pow_val = 1
    for u in range(period):
        out[u] = np.exp(2j * np.pi * c * pow_val / M)
        pow_val = (pow_val * (1 + p)) % M
    return out


def G_length_period(p: int, r: int, c: int = 1) -> np.ndarray:
    """G[a] := Σ_{s=0}^{p^r - 1} f_p(s) exp(-2πi a s / p^r) for a ∈ Z/p^r.

    Computed via length-period FFT. Equals F̂_full(p·a) / p.
    """
    f = f_p_period(p, r, c)
    return np.fft.fft(f)


def predicted_supp_a_indices(p: int, r: int) -> np.ndarray:
    """Predicted G[a] support: a ∈ {1, 1+p, 1+2p, ..., 1+(p^{r-1}-1)p}, size p^{r-1}.
    For r=1: just {a=1} of size 1."""
    pr = p ** r
    return np.arange(1, pr, p, dtype=np.int64)


def cell_verify_G(p: int, r: int, c: int = 1, off_supp_threshold: float = 1e-10) -> dict:
    """Compute G[a], verify magnitudes match p^{(r+1)/2} on supp {a ≡ 1 mod p}
    and zero elsewhere."""
    t0 = time.time()
    M = p ** (r + 1)
    period = p ** r
    G = G_length_period(p, r, c=c)
    mags = np.abs(G)

    supp = predicted_supp_a_indices(p, r)
    predicted_supp_size = len(supp)
    predicted_mag_short = p ** ((r + 1) / 2)
    predicted_mag_full = p ** ((r + 3) / 2)

    # On-support magnitudes
    mags_supp = mags[supp]
    max_supp = float(mags_supp.max())
    min_supp = float(mags_supp.min())
    median_supp = float(np.median(mags_supp))
    rel_dev_supp = np.abs(mags_supp - predicted_mag_short) / predicted_mag_short
    max_rel_dev = float(rel_dev_supp.max())

    # Off-support
    mask = np.ones(period, dtype=bool)
    mask[supp] = False
    mags_off = mags[mask]
    max_off = float(mags_off.max()) if mags_off.size else 0.0

    # Numerical support
    numerical_supp = np.where(mags > off_supp_threshold)[0]
    numerical_supp_size = int(len(numerical_supp))
    sym_diff = int(len(set(numerical_supp.tolist()) ^ set(supp.tolist())))

    elapsed = time.time() - t0

    return {
        "p": p, "r": r, "M": M, "period": period,
        "predicted_supp_size": predicted_supp_size,
        "numerical_supp_size": numerical_supp_size,
        "supp_symmetric_diff": sym_diff,
        "predicted_mag_short": predicted_mag_short,
        "predicted_mag_full": predicted_mag_full,
        "max_mag_on_supp": max_supp,
        "min_mag_on_supp": min_supp,
        "median_mag_on_supp": median_supp,
        "max_rel_dev_from_pred_short": max_rel_dev,
        "max_mag_off_supp": max_off,
        "off_supp_threshold": off_supp_threshold,
        "elapsed_s": elapsed,
    }


def boundary_p2(r: int, c: int = 1) -> dict:
    """Boundary at p=2: (1+p)=3 has order 2^{r-1} (not 2^r) in (Z/2^{r+1})^×
    for r ≥ 2, since principal units 1 + 2Z_2 mod 2^{r+1} are
    Z/2 × Z/2^{r-1}, NOT cyclic. Test the F̂_full structure."""
    p = 2
    t0 = time.time()
    M = p ** (r + 1)
    # Actual order of (1+p)=3 mod M
    pow_val = 1
    period_actual = None
    for k in range(1, M + 1):
        pow_val = (pow_val * 3) % M
        if pow_val == 1:
            period_actual = k
            break

    # Build f over the actual period, then test G via period-length FFT.
    f_vals = np.empty(period_actual, dtype=np.complex128)
    pow_val = 1
    for u in range(period_actual):
        f_vals[u] = np.exp(2j * np.pi * c * pow_val / M)
        pow_val = (pow_val * 3) % M

    G = np.fft.fft(f_vals)  # length-period_actual FFT
    mags = np.abs(G)

    # If we test the "naive" theorem (supp = {a≡1 mod p=2} in Z/period_actual)
    # the naive predicted magnitude would be 2^{(r+1)/2}.
    # Naive support: a ∈ {1, 3, 5, ...} ⊆ Z/period_actual
    naive_supp = np.arange(1, period_actual, p)
    naive_pred = 2 ** ((r + 1) / 2)

    max_supp_naive = float(mags[naive_supp].max()) if len(naive_supp) else 0.0
    mask = np.ones(period_actual, dtype=bool)
    mask[naive_supp] = False
    max_off_naive = float(mags[mask].max()) if mask.any() else 0.0

    return {
        "p": p, "r": r, "M": M,
        "period_expected_p^r": p ** r,
        "period_actual_of_(1+p)_mod_M": period_actual,
        "naive_pred_mag": naive_pred,
        "max_G_on_naive_supp": max_supp_naive,
        "max_G_off_naive_supp": max_off_naive,
        "max_G_overall": float(mags.max()),
        "G_nonzero_count": int((mags > 1e-10).sum()),
        "elapsed_s": time.time() - t0,
    }


def main():
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    rows = []

    # ============================================================
    # Phase 1: extension to higher primes
    # ============================================================
    print("=" * 88)
    print("PHASE 1: extension to higher primes — p ∈ {11..31}, r ∈ {2, 3}")
    print("  Object: G[a] = length-period FFT of f_p. F̂_full(p·a) = p · G[a].")
    print("=" * 88)
    print(f"  {'p':>3} {'r':>2} {'M':>12} {'period':>8} {'|supp|':>8} {'sym_diff':>9} "
          f"{'max|G|/p^((r+1)/2)':>22} {'min/pred':>14} {'max_off':>12} {'t(s)':>6}")

    phase1_primes = [11, 13, 17, 19, 23, 29, 31]
    phase1_r = [2, 3]
    phase1_passes = 0
    phase1_total = 0

    for p in phase1_primes:
        for r in phase1_r:
            res = cell_verify_G(p, r)
            phase1_total += 1
            pass_cell = (
                res["max_rel_dev_from_pred_short"] < 1e-12
                and res["max_mag_off_supp"] < 1e-10
                and res["supp_symmetric_diff"] == 0
                and res["predicted_supp_size"] == res["numerical_supp_size"]
            )
            if pass_cell:
                phase1_passes += 1
            rows.append({"phase": "1", **res, "pass_cell": pass_cell})
            print(f"  {p:>3} {r:>2} {res['M']:>12} {res['period']:>8} {res['predicted_supp_size']:>8} "
                  f"{res['supp_symmetric_diff']:>9} "
                  f"{res['max_mag_on_supp'] / res['predicted_mag_short']:>22.16f} "
                  f"{res['min_mag_on_supp'] / res['predicted_mag_short']:>14.10f} "
                  f"{res['max_mag_off_supp']:>12.2e} "
                  f"{res['elapsed_s']:>6.2f}")

    print()
    print(f"PHASE 1 RESULT: {phase1_passes}/{phase1_total} cells pass tight threshold (rel_dev<1e-12, off<1e-10, sym_diff=0)")
    print()

    # ============================================================
    # Phase 2: extension to higher r
    # ============================================================
    print("=" * 88)
    print("PHASE 2: extension to higher r — p ∈ {3, 5, 7}, r ∈ {4, 5, 6}")
    print("=" * 88)
    print(f"  {'p':>3} {'r':>2} {'M':>12} {'period':>8} {'|supp|':>8} {'sym_diff':>9} "
          f"{'max|G|/p^((r+1)/2)':>22} {'min/pred':>14} {'max_off':>12} {'t(s)':>6}")

    phase2_primes = [3, 5, 7]
    phase2_r = [4, 5, 6]
    phase2_passes = 0
    phase2_total = 0

    for p in phase2_primes:
        for r in phase2_r:
            res = cell_verify_G(p, r)
            phase2_total += 1
            pass_cell = (
                res["max_rel_dev_from_pred_short"] < 1e-12
                and res["max_mag_off_supp"] < 1e-10
                and res["supp_symmetric_diff"] == 0
                and res["predicted_supp_size"] == res["numerical_supp_size"]
            )
            if pass_cell:
                phase2_passes += 1
            rows.append({"phase": "2", **res, "pass_cell": pass_cell})
            print(f"  {p:>3} {r:>2} {res['M']:>12} {res['period']:>8} {res['predicted_supp_size']:>8} "
                  f"{res['supp_symmetric_diff']:>9} "
                  f"{res['max_mag_on_supp'] / res['predicted_mag_short']:>22.16f} "
                  f"{res['min_mag_on_supp'] / res['predicted_mag_short']:>14.10f} "
                  f"{res['max_mag_off_supp']:>12.2e} "
                  f"{res['elapsed_s']:>6.2f}")

    print()
    print(f"PHASE 2 RESULT: {phase2_passes}/{phase2_total} cells pass tight threshold")
    print()

    # ============================================================
    # Phase 4: boundary cases
    # ============================================================
    print("=" * 88)
    print("PHASE 4: boundary cases — p=2 (excluded by theorem), r=1 (theorem says r≥2)")
    print("=" * 88)
    print()
    print("Boundary 1: p = 2 at r ∈ {2, 3, 4} — testing whether theorem fails (it should).")
    print(f"  {'r':>2} {'M':>8} {'p^r expected':>14} {'(1+p) order':>14} "
          f"{'naive_pred':>12} {'max|G| on naive_supp':>22} "
          f"{'max|G| off naive_supp':>22}")

    for r in [2, 3, 4]:
        res = boundary_p2(r)
        rows.append({"phase": "4_p2", **res})
        print(f"  {r:>2} {res['M']:>8} {res['period_expected_p^r']:>14} "
              f"{res['period_actual_of_(1+p)_mod_M']:>14} "
              f"{res['naive_pred_mag']:>12.4f} {res['max_G_on_naive_supp']:>22.6f} "
              f"{res['max_G_off_naive_supp']:>22.6f}")

    print()
    print("Boundary 2: r = 1 at p ∈ {3, 5, 7, 11} — test whether theorem ACTUALLY excludes r=1.")
    print(f"  {'p':>3} {'r':>2} {'M':>6} {'|supp|':>8} "
          f"{'pred(short)':>12} {'max|G|':>14} {'min|G|':>14} "
          f"{'max_off':>12} {'rel_dev':>12} {'sym_diff':>10} {'PASS':>6}")

    boundary_r1_passes = 0
    for p in [3, 5, 7, 11]:
        res = cell_verify_G(p, 1)
        pass_cell = (
            res["max_rel_dev_from_pred_short"] < 1e-12
            and res["max_mag_off_supp"] < 1e-10
            and res["supp_symmetric_diff"] == 0
        )
        if pass_cell:
            boundary_r1_passes += 1
        rows.append({"phase": "4_r1", **res, "pass_cell": pass_cell})
        print(f"  {p:>3} {res['r']:>2} {res['M']:>6} {res['predicted_supp_size']:>8} "
              f"{res['predicted_mag_short']:>12.6f} {res['max_mag_on_supp']:>14.10f} "
              f"{res['min_mag_on_supp']:>14.10f} {res['max_mag_off_supp']:>12.2e} "
              f"{res['max_rel_dev_from_pred_short']:>12.2e} "
              f"{res['supp_symmetric_diff']:>10} "
              f"{'YES' if pass_cell else 'no':>6}")

    print()
    print(f"  r=1 boundary verdict: {boundary_r1_passes}/4 cells satisfy the theorem.")
    print(f"  → If 4/4 pass, theorem actually holds at r=1 too; pre-reg's r≥2 was conservative.")
    print()

    # ============================================================
    # Phase 5: support characterization summary
    # ============================================================
    print("=" * 88)
    print("PHASE 5: support characterization summary")
    print("=" * 88)
    print()
    print("Predicted G support: {a ∈ Z/p^r : a ≡ 1 (mod p)}, |supp| = p^{r-1}.")
    print("Equivalent: F̂_full has support {p·a (mod M) : a ≡ 1 (mod p)}.")
    print()
    s1 = sum(1 for r_ in rows if r_.get("phase") == "1" and r_.get("supp_symmetric_diff") == 0)
    s2 = sum(1 for r_ in rows if r_.get("phase") == "2" and r_.get("supp_symmetric_diff") == 0)
    sr1 = sum(1 for r_ in rows if r_.get("phase") == "4_r1" and r_.get("supp_symmetric_diff") == 0)
    print(f"  Phase 1: {s1}/{phase1_total} cells have symmetric_diff = 0")
    print(f"  Phase 2: {s2}/{phase2_total} cells have symmetric_diff = 0")
    print(f"  Phase 4 r=1: {sr1}/4 cells have symmetric_diff = 0")
    print()

    max_off_overall = max(r_.get("max_mag_off_supp", 0) for r_ in rows
                          if r_.get("phase") in ("1", "2", "4_r1"))
    print(f"  Max off-support |G| across Phases 1, 2, 4_r1: {max_off_overall:.2e}")
    print()

    # ============================================================
    # Write CSV
    # ============================================================
    all_keys = sorted({k for row in rows for k in row.keys()})
    with open(OUTPATH, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["phase"] + [k for k in all_keys if k != "phase"])
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print(f"[write] {OUTPATH}")

    # ============================================================
    # Summary
    # ============================================================
    print()
    print("=" * 88)
    print("OVERALL SUMMARY")
    print("=" * 88)
    print(f"  Phase 1: {phase1_passes}/{phase1_total} cells passed (extension to higher primes)")
    print(f"  Phase 2: {phase2_passes}/{phase2_total} cells passed (extension to higher r)")
    print(f"  Phase 4: r=1 at {boundary_r1_passes}/4 primes; p=2 documented")
    print(f"  Phase 5: support symmetric_diff = 0 across {s1 + s2 + sr1}/{phase1_total + phase2_total + 4} cells")
    print()
    full_pass = (phase1_passes == phase1_total) and (phase2_passes == phase2_total) and (boundary_r1_passes == 4)
    print(f"  → Phases 1+2 full pass: {(phase1_passes == phase1_total) and (phase2_passes == phase2_total)}")
    print(f"  → Phases 1+2+r=1 full pass: {full_pass}")


if __name__ == "__main__":
    main()
