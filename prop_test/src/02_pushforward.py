"""
02_pushforward.py — compute post-prefix re-entry signatures at k_target=12
for all (r, m) pairs and emit pushforward distributions per starting r.

For each odd r mod 2^k_base, compute the prefix signature S(r) = (a_final, c_final, j, prefix_steps).
For each m in {0, ..., 2^k_target - 1}, compute:
    n'(r, m) = a_final(r) * m + c_final(r)
Then reduce n' mod 2^k_target and look up its prefix signature at k_target.

The pushforward distribution at k_target for starting r is the empirical distribution
of S(n' mod 2^k_target) over m. We summarize each pushforward by:
- mean j', std j'
- mean prefix_steps', std prefix_steps'
- distribution over distinct a_final' values
- mean log(a_final')
- count of (r, m) pairs whose n' mod 2^k_target falls in each a_final' bucket

If the pushforward is uniform over odd residues mod 2^k_target (as the arithmetic
preview predicts), these summaries will be identical across r modulo finite-sample
variation. F1 quantifies whether starting r predicts any of these summaries.

Writes:
  prop_test/post_prefix_signatures_k12.csv  — full (r, m) pairs with signatures
                                                 (kept compressed; subsampled if large)
  prop_test/correspondence_summary.csv      — per-r pushforward summary statistics
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT


def load_signatures(path: Path) -> dict[int, tuple[int, int, int, int]]:
    sigs: dict[int, tuple[int, int, int, int]] = {}
    with open(path, "r", encoding="utf-8") as fh:
        rd = csv.reader(fh)
        next(rd)  # header
        for row in rd:
            r = int(row[0])
            sigs[r] = (int(row[1]), int(row[2]), int(row[3]), int(row[4]))
    return sigs


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    K_BASE = 6
    K_TARGET = 12
    MOD_TARGET = 2 ** K_TARGET

    print(f"# Pushforward computation: k_base={K_BASE} -> k_target={K_TARGET}")

    sigs_base = load_signatures(OUTDIR / "prefix_signatures_k6.csv")
    sigs_target = load_signatures(OUTDIR / "prefix_signatures_k12.csv")

    print(f"  loaded {len(sigs_base)} base signatures, {len(sigs_target)} target signatures")

    # Full enumeration of (r, m) pairs
    # m ranges over {0, ..., 2^K_TARGET - 1}. For each r, all 4096 values of m.
    # Total: 32 * 4096 = 131072 (r, m) pairs.

    full_path = OUTDIR / "post_prefix_signatures_k12.csv"
    summary_rows = []

    with open(full_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "r", "m", "a_final_r", "c_final_r", "j_r", "prefix_steps_r",
            "n_prime_mod_2k_target", "a_final_prime", "c_final_prime",
            "j_prime", "prefix_steps_prime"
        ])

        for r, (a_final, c_final, j, prefix_steps) in sorted(sigs_base.items()):
            # Collect pushforward statistics
            j_primes = []
            prefix_steps_primes = []
            a_final_primes = []
            n_prime_residues = []

            for m in range(MOD_TARGET):
                n_prime = a_final * m + c_final
                n_prime_mod = n_prime % MOD_TARGET
                # n' must be odd (the original Collatz lifts are odd integers).
                # If n_prime_mod is even, skip — but with c_final and a_final both
                # known, parity of n' = parity of (a_final*m + c_final) = parity of
                # (m + c_final) (since a_final is odd). So we only consider m with
                # m + c_final odd.
                if n_prime_mod % 2 == 0:
                    continue
                sig_prime = sigs_target[n_prime_mod]
                af_p, cf_p, j_p, ps_p = sig_prime

                w.writerow([
                    r, m, a_final, c_final, j, prefix_steps,
                    n_prime_mod, af_p, cf_p, j_p, ps_p
                ])

                j_primes.append(j_p)
                prefix_steps_primes.append(ps_p)
                a_final_primes.append(af_p)
                n_prime_residues.append(n_prime_mod)

            n_valid = len(j_primes)
            mean_j = sum(j_primes) / n_valid
            var_j = sum((x - mean_j) ** 2 for x in j_primes) / n_valid
            std_j = math.sqrt(var_j)

            mean_ps = sum(prefix_steps_primes) / n_valid
            std_ps = math.sqrt(sum((x - mean_ps) ** 2 for x in prefix_steps_primes) / n_valid)

            mean_log_af = sum(math.log(x) for x in a_final_primes) / n_valid

            # Distribution over the 12 distinct a_final' values at k=12
            af_counts = {}
            for af in a_final_primes:
                af_counts[af] = af_counts.get(af, 0) + 1

            summary_rows.append({
                "r": r,
                "a_final_r": a_final,
                "c_final_r": c_final,
                "j_r": j,
                "prefix_steps_r": prefix_steps,
                "n_valid": n_valid,
                "mean_j_prime": mean_j,
                "std_j_prime": std_j,
                "mean_prefix_steps_prime": mean_ps,
                "std_prefix_steps_prime": std_ps,
                "mean_log_a_final_prime": mean_log_af,
                "af_counts": af_counts,
            })

            print(
                f"  r={r:>3}  j={j}  steps={prefix_steps:>3}  "
                f"n_valid={n_valid:>5}  "
                f"mean_j'={mean_j:>6.3f}  std={std_j:.3f}  "
                f"mean_log_af'={mean_log_af:>5.3f}"
            )

    # Write summary
    summary_path = OUTDIR / "correspondence_summary.csv"
    with open(summary_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "r", "a_final_r", "c_final_r", "j_r", "prefix_steps_r",
            "n_valid", "mean_j_prime", "std_j_prime",
            "mean_prefix_steps_prime", "std_prefix_steps_prime",
            "mean_log_a_final_prime",
            "af_count_3", "af_count_9", "af_count_27", "af_count_81",
            "af_count_243", "af_count_729", "af_count_2187", "af_count_6561",
            "af_count_19683", "af_count_59049", "af_count_177147", "af_count_531441",
        ])
        af_keys = [3, 9, 27, 81, 243, 729, 2187, 6561, 19683, 59049, 177147, 531441]
        for row in summary_rows:
            af_counts = row["af_counts"]
            w.writerow([
                row["r"], row["a_final_r"], row["c_final_r"],
                row["j_r"], row["prefix_steps_r"], row["n_valid"],
                f"{row['mean_j_prime']:.6f}",
                f"{row['std_j_prime']:.6f}",
                f"{row['mean_prefix_steps_prime']:.6f}",
                f"{row['std_prefix_steps_prime']:.6f}",
                f"{row['mean_log_a_final_prime']:.6f}",
                *[af_counts.get(af, 0) for af in af_keys]
            ])

    print()
    print(f"[done] Pushforward files written:")
    print(f"  {full_path}")
    print(f"  {summary_path}")

    # Quick check: are pushforward summaries similar across r?
    print()
    print("# Summary across all 32 starting r (pushforward at k=12)")
    mean_js = [row["mean_j_prime"] for row in summary_rows]
    mean_ps = [row["mean_prefix_steps_prime"] for row in summary_rows]
    mean_log_afs = [row["mean_log_a_final_prime"] for row in summary_rows]
    print(f"  mean_j': overall mean={sum(mean_js)/len(mean_js):.4f}, "
          f"min={min(mean_js):.4f}, max={max(mean_js):.4f}, "
          f"spread={max(mean_js)-min(mean_js):.4f}")
    print(f"  mean_prefix_steps': overall mean={sum(mean_ps)/len(mean_ps):.4f}, "
          f"min={min(mean_ps):.4f}, max={max(mean_ps):.4f}, "
          f"spread={max(mean_ps)-min(mean_ps):.4f}")
    print(f"  mean_log_a_final': overall mean={sum(mean_log_afs)/len(mean_log_afs):.4f}, "
          f"min={min(mean_log_afs):.4f}, max={max(mean_log_afs):.4f}, "
          f"spread={max(mean_log_afs)-min(mean_log_afs):.4f}")


if __name__ == "__main__":
    main()
