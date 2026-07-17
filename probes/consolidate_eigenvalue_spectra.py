"""
Consolidate K_k eigenvalue spectra from existing result files for audio
rendering. NO new compute. Pure data extraction + reformatting.

Sources:
  - C:/Collatz/result_eigenvalue_spectrum.csv (top-10 per k=5,6,7; has arg)
  - C:/Collatz/experiments_output/result_77_4_K_spectrum_data.csv
    (all eigenvalues for k=3..6)

Strategy:
  - k=5, k=6: prefer supplementary (full spectrum), take top 20 by |λ|.
  - k=7: only primary file has it, max 10 eigenvalues available.
"""
from __future__ import annotations
import csv
import os
import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = Path(r"C:\Collatz\audio_data")
OUTDIR.mkdir(exist_ok=True)
OUT_CSV = OUTDIR / "eigenvalue_spectra_consolidated.csv"
OUT_TXT = OUTDIR / "eigenvalue_summary.txt"

PRIMARY = Path(r"C:\Collatz\result_eigenvalue_spectrum.csv")
SUPPLEMENTARY = Path(
    r"C:\Collatz\experiments_output\result_77_4_K_spectrum_data.csv"
)
COMPLEX_TOL = 1e-12
TOP_N = 20


def load_primary():
    """Returns dict k -> list of (real, imag, abs_lambda, arg)."""
    by_k = {}
    with open(PRIMARY) as f:
        reader = csv.DictReader(f)
        for row in reader:
            k = int(row["k"])
            re_ = float(row["lambda_real"])
            im_ = float(row["lambda_imag"])
            abs_ = float(row["abs_lambda"])
            arg_ = float(row["arg_lambda"])
            by_k.setdefault(k, []).append((re_, im_, abs_, arg_))
    return by_k


def load_supplementary():
    """Returns dict k -> list of (real, imag)."""
    by_k = {}
    with open(SUPPLEMENTARY) as f:
        reader = csv.DictReader(f)
        for row in reader:
            k = int(row["k"])
            re_ = float(row["lambda_real"])
            im_ = float(row["lambda_imag"])
            by_k.setdefault(k, []).append((re_, im_))
    return by_k


def consolidate():
    primary = load_primary()
    supp = load_supplementary()

    sources = {}  # k -> source label
    eigs_by_k = {}  # k -> list of (real, imag, abs, arg)

    for k in [5, 6]:
        if k in supp and len(supp[k]) >= TOP_N:
            entries = []
            for re_, im_ in supp[k]:
                ab = float(np.hypot(re_, im_))
                arg = float(np.arctan2(im_, re_))
                entries.append((re_, im_, ab, arg))
            entries.sort(key=lambda e: -e[2])
            eigs_by_k[k] = entries[:TOP_N]
            sources[k] = (
                f"supplementary (full spectrum, "
                f"{len(supp[k])} eigvals, top-{TOP_N} by |λ|)"
            )
        elif k in primary:
            entries = list(primary[k])
            entries.sort(key=lambda e: -e[2])
            eigs_by_k[k] = entries
            sources[k] = (
                f"primary fallback (top-{len(primary[k])} only, "
                f"requested top-{TOP_N})"
            )
        else:
            sources[k] = "MISSING"

    if 7 in primary:
        entries = list(primary[7])
        entries.sort(key=lambda e: -e[2])
        eigs_by_k[7] = entries
        n7 = len(entries)
        sources[7] = (
            f"primary only (top-{n7}, requested top-{TOP_N}; supplementary "
            f"file does not contain k=7)"
        )
    else:
        sources[7] = "MISSING"

    return eigs_by_k, sources


def categorize_args(args, mags):
    """Look for clustering near roots of unity for q = 3, 4, 5, 6, 8, 9, 27."""
    args = np.asarray(args)
    mags = np.asarray(mags)
    # Exclude the trivial leading eigenvalue at |λ| ≈ 1; otherwise keep all.
    nontrivial = mags < 0.99
    if not nontrivial.any():
        return "no nontrivial eigenvalues"

    sub = args[nontrivial]
    candidates = [
        ("3rd", 3), ("4th", 4), ("5th", 5), ("6th", 6),
        ("8th", 8), ("9th", 9), ("27th", 27),
    ]
    matches = []
    for name, q in candidates:
        spacing = 2 * np.pi / q
        scaled = sub / spacing
        nearest = np.round(scaled)
        dev_rad = np.abs(scaled - nearest) * spacing
        match_count = int((dev_rad < 0.01).sum())  # within ~0.6°
        if match_count > 0:
            matches.append(f"{match_count}/{len(sub)} within 0.6° of {name} roots")
    if matches:
        return "; ".join(matches)
    return "no alignment with q-th roots for q ∈ {3,4,5,6,8,9,27}"


def main():
    eigs_by_k, sources = consolidate()

    # Write consolidated CSV
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "rank", "magnitude", "argument_rad", "is_complex",
                    "real_part", "imag_part"])
        n_rows = 0
        per_k_count = {}
        for k in sorted(eigs_by_k):
            entries = eigs_by_k[k]
            for rank, (re_, im_, ab, arg) in enumerate(entries, 1):
                is_complex = bool(abs(im_) > COMPLEX_TOL)
                w.writerow([k, rank, f"{ab:.12e}", f"{arg:.12e}",
                            int(is_complex),
                            f"{re_:.12e}", f"{im_:.12e}"])
                n_rows += 1
            per_k_count[k] = len(entries)

    print(f"Wrote {OUT_CSV}")
    print(f"Total rows: {n_rows}")
    for k in sorted(per_k_count):
        print(f"  k = {k}: {per_k_count[k]} eigenvalues   "
              f"[source: {sources[k]}]")

    # Summary text
    lines = []
    lines.append("Eigenvalue Spectrum Consolidation — Summary")
    lines.append("=" * 60)
    lines.append("")
    lines.append("SOURCES")
    lines.append("-" * 60)
    lines.append(f"  primary:       {PRIMARY}")
    lines.append(f"  supplementary: {SUPPLEMENTARY}")
    lines.append("")
    lines.append("PER-k EXTRACTION")
    lines.append("-" * 60)
    for k in sorted(eigs_by_k):
        lines.append(f"  k = {k}: source = {sources[k]}")
        lines.append(f"          rows = {per_k_count[k]}")
    lines.append("")
    lines.append("PER-k STATISTICS")
    lines.append("-" * 60)
    for k in sorted(eigs_by_k):
        entries = eigs_by_k[k]
        mags = np.array([e[2] for e in entries])
        args = np.array([e[3] for e in entries])
        is_complex = np.array([abs(e[1]) > COMPLEX_TOL for e in entries])
        n_complex = int(is_complex.sum())
        # Magnitude range, geomean of top 10
        mag_min = float(mags.min())
        mag_max = float(mags.max())
        top10_mags = mags[:10]
        # Geometric mean of top-10 (skip exactly-zero magnitudes)
        nz = top10_mags[top10_mags > 0]
        gmean10 = float(np.exp(np.log(nz).mean())) if len(nz) > 0 else 0.0
        # Arg distribution
        arg_min = float(args.min())
        arg_max = float(args.max())
        arg_pattern = categorize_args(args, mags)

        lines.append(f"")
        lines.append(f"  k = {k}  ({per_k_count[k]} eigenvalues)")
        lines.append(f"    complex eigvals: {n_complex} / "
                     f"{per_k_count[k]}")
        lines.append(f"    |λ| range:       [{mag_min:.4e}, "
                     f"{mag_max:.4e}]")
        lines.append(f"    |λ| geom mean (top 10): {gmean10:.4e}")
        lines.append(f"    arg range:       [{arg_min:+.4f}, "
                     f"{arg_max:+.4f}] rad")
        lines.append(f"    arg pattern:     {arg_pattern}")

        # Inspect spread/clustering qualitatively
        nontrivial_mask = (mags > 0.001) & (mags < 0.99)
        if nontrivial_mask.any():
            sub_args = args[nontrivial_mask]
            uniq = np.unique(np.round(sub_args, 3))
            n_unique = len(uniq)
            # Are they all near 0? (real eigenvalues clustering at angle 0 or π)
            near_real = ((np.abs(sub_args) < 0.01)
                         | (np.abs(sub_args - np.pi) < 0.01)
                         | (np.abs(sub_args + np.pi) < 0.01))
            n_real_like = int(near_real.sum())
            distribution = (
                "clustered near real axis"
                if n_real_like / max(len(sub_args), 1) > 0.6
                else f"{n_unique} distinct values among {len(sub_args)}"
            )
            lines.append(f"    nontrivial-mode arg distribution: "
                         f"{distribution}")

    lines.append("")
    lines.append("ANOMALIES / CAVEATS")
    lines.append("-" * 60)
    if 7 in eigs_by_k and per_k_count[7] < TOP_N:
        lines.append(f"  k=7: only {per_k_count[7]} eigenvalues available "
                     f"(requested {TOP_N}). Supplementary file does not "
                     "extend to k=7. To get the full top-20 at k=7 would "
                     "require a new compute on K_7 (1458 states).")
    for k in [5, 6]:
        if k in eigs_by_k and per_k_count[k] < TOP_N:
            lines.append(f"  k={k}: only {per_k_count[k]} available; "
                         "fewer than requested.")

    # Look for degenerate eigenvalues (multiple at same |λ| within tolerance)
    for k in sorted(eigs_by_k):
        entries = eigs_by_k[k]
        mags = sorted([e[2] for e in entries], reverse=True)
        # Cluster mags within 1e-9 relative
        for i in range(len(mags) - 1):
            if mags[i] > 1e-12:
                rel_diff = abs(mags[i] - mags[i + 1]) / mags[i]
                if rel_diff < 1e-9 and mags[i] > 1e-6:
                    lines.append(
                        f"  k={k}: ranks {i+1},{i+2} have same |λ| "
                        f"= {mags[i]:.6e} (degenerate / complex pair)"
                    )
                    break

    OUT_TXT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {OUT_TXT}")
    print()
    print("=" * 60)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
