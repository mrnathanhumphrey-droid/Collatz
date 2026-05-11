"""
01_signatures.py — compute prefix-decomposition signatures at k=6 and k=12.

The signature S(r) at modular resolution k is the 4-tuple:
  (a_final, c_final, j, prefix_steps)

produced by symbolic Collatz iteration on state (a, c) starting at (a=2^k, c=r),
with the rules from writeup.md Result 3:
  - (a even, c even) -> (a/2, c/2)
  - (a even, c odd)  -> (3a, 3c + 1)
  - (a odd)          -> terminate

Writes:
  prop_test/prefix_signatures_k6.csv   (32 rows, odd residues 1..63)
  prop_test/prefix_signatures_k12.csv  (2048 rows, odd residues 1..4095)
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT


def prefix_decompose(r: int, k: int) -> tuple[int, int, int, int]:
    """Symbolic prefix decomposition.

    Returns (a_final, c_final, j, prefix_steps) where a_final = 3**j.
    """
    a, c = 2 ** k, r
    steps = 0
    j = 0
    while a % 2 == 0:
        if c % 2 == 0:
            a, c = a // 2, c // 2
        else:
            a, c = 3 * a, 3 * c + 1
            j += 1
        steps += 1
        if steps > 200:
            raise RuntimeError(f"Prefix decomposition exceeded 200 steps for r={r}, k={k}")
    if a != 3 ** j:
        raise AssertionError(
            f"Expected a_final = 3^j = {3 ** j}, got {a} (r={r}, k={k}, j={j})"
        )
    return a, c, j, steps


def write_signatures(k: int, outpath: Path) -> dict[int, tuple[int, int, int, int]]:
    """Compute signatures for all odd residues mod 2^k, write CSV, return dict."""
    sigs: dict[int, tuple[int, int, int, int]] = {}
    rows = []
    for r in range(1, 2 ** k, 2):  # odd residues only
        sig = prefix_decompose(r, k)
        sigs[r] = sig
        rows.append((r, *sig))

    with open(outpath, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["r", "a_final", "c_final", "j", "prefix_steps"])
        for row in rows:
            w.writerow(row)
    return sigs


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    print("# Computing prefix signatures")

    sigs_k6 = write_signatures(6, OUTDIR / "prefix_signatures_k6.csv")
    print(f"  k=6: {len(sigs_k6)} odd residues, signatures written")
    a_finals_k6 = sorted({s[0] for s in sigs_k6.values()})
    print(f"        distinct a_final values: {a_finals_k6}")
    j_dist_k6 = {}
    for r, sig in sigs_k6.items():
        j_dist_k6[sig[2]] = j_dist_k6.get(sig[2], 0) + 1
    print(f"        j distribution: {dict(sorted(j_dist_k6.items()))}")

    sigs_k12 = write_signatures(12, OUTDIR / "prefix_signatures_k12.csv")
    print(f"  k=12: {len(sigs_k12)} odd residues, signatures written")
    a_finals_k12 = sorted({s[0] for s in sigs_k12.values()})
    print(f"         distinct a_final values: {a_finals_k12}")
    j_dist_k12 = {}
    for r, sig in sigs_k12.items():
        j_dist_k12[sig[2]] = j_dist_k12.get(sig[2], 0) + 1
    print(f"         j distribution: {dict(sorted(j_dist_k12.items()))}")

    print()
    print("# Adversarial safeguard A1 (determinism)")
    for r in (1, 27):
        s1 = prefix_decompose(r, 6)
        s2 = prefix_decompose(r, 6)
        match = s1 == s2
        print(f"  r={r}, k=6: run1={s1}  run2={s2}  match={match}")
        if not match:
            raise RuntimeError("Determinism violated!")

    print()
    print("# Adversarial safeguard A2 (n=27 spot check)")
    # n = 27 has stopping time 111, well-known long trajectory.
    # Apply prefix decomposition to r = 27 mod 64; verify that n = 27 itself
    # (m = 0 in the parameterization n = 64m + 27) reaches a_final*0 + c_final
    # at the predicted prefix step count.
    r = 27
    a_final, c_final, j, steps = prefix_decompose(r, 6)
    print(f"  r=27 mod 64: a_final={a_final}, c_final={c_final}, j={j}, steps={steps}")
    # Direct Collatz on n=27, recording the state after `steps` total odd+even steps:
    n = 27
    for i in range(steps):
        if n % 2 == 0:
            n = n // 2
        else:
            n = 3 * n + 1
    # After `steps` Collatz steps starting from n=27, the value should equal
    # a_final * m + c_final where m = (27 - 27) / 64 = 0 ... wait this is wrong.
    # n=27 lifts as n = 64*0 + 27, so m=0. Post-prefix value = a_final*0 + c_final = c_final.
    expected = a_final * 0 + c_final
    match = (n == expected)
    print(f"  After {steps} direct Collatz steps from n=27: n={n}, expected={expected}, match={match}")
    if not match:
        raise RuntimeError(
            f"Prefix decomposition for r=27 mod 64 inconsistent with direct Collatz: "
            f"expected {expected}, got {n}"
        )

    # Also verify for r=27 with m=1: n=91 -> should reach a_final*1 + c_final after steps
    n = 64 * 1 + 27   # = 91
    for i in range(steps):
        if n % 2 == 0:
            n = n // 2
        else:
            n = 3 * n + 1
    expected = a_final * 1 + c_final
    match = (n == expected)
    print(f"  After {steps} direct Collatz steps from n=91 (m=1): n={n}, expected={expected}, match={match}")
    if not match:
        raise RuntimeError(
            f"Prefix decomposition for r=27 mod 64, m=1 inconsistent: "
            f"expected {expected}, got {n}"
        )

    print()
    print("[done] Signature files written:")
    print(f"  {OUTDIR / 'prefix_signatures_k6.csv'}")
    print(f"  {OUTDIR / 'prefix_signatures_k12.csv'}")


if __name__ == "__main__":
    main()
