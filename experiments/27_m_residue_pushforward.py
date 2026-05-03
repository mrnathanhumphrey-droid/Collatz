"""
Experiment 27 — pushforward of the trajectory measure onto m mod 2^k.

Mechanism hypothesis test for the v=4 / v=10 / v=16 deviations identified in
exp 25. v_2(3m+1) = k *exactly* iff 3m+1 ≡ 2^k (mod 2^(k+1)), i.e. m is
determined modulo 2^(k+1):
    v=4  ↔  m ≡ 5     (mod 32)        [2^5  = 32]
    v=10 ↔  m ≡ 341   (mod 2048)      [2^11 = 2048]
    v=16 ↔  m ≡ 21845 (mod 131072)    [2^17 = 131072]

If the trajectory measure on iterates m_t over-weights these specific residues,
the spike pattern in the v-distribution is mechanistically explained.

Sample N_start uniform odd starts, walk T Syracuse steps. At each step,
*before* the 3m+1 update, record m mod 32, m mod 2048, m mod 131072. Compute
empirical density on each modulus and compare to uniform-on-odd-residues.

Usage:
    python 27_m_residue_pushforward.py --N_start 100000000 --T 200
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

MOD_SMALL = 32          # 2^5  → v=4 exactly
MOD_MID   = 2048        # 2^11 → v=10 exactly
MOD_LARGE = 131072      # 2^17 → v=16 exactly


@njit(parallel=True, cache=True)
def collect_m_residues(starts, T, mod_s, mod_m, mod_l):
    """Chunk-parallel: each prange iteration owns its own three histogram rows.
    Records m mod {mod_s, mod_m, mod_l} at every Syracuse step before update.
    """
    n = len(starts)
    # Coarser chunking to keep memory budget modest:
    # 32 chunks * 131072 cells * 8 bytes = 32 MB for the large histogram.
    n_chunks = 32
    chunk = (n + n_chunks - 1) // n_chunks
    hist_s = np.zeros((n_chunks, mod_s), dtype=np.int64)
    hist_m = np.zeros((n_chunks, mod_m), dtype=np.int64)
    hist_l = np.zeros((n_chunks, mod_l), dtype=np.int64)
    cnt = np.zeros(n_chunks, dtype=np.int64)
    for c in prange(n_chunks):
        start_i = c * chunk
        end_i = min(start_i + chunk, n)
        for i in range(start_i, end_i):
            m = starts[i]
            for _ in range(T):
                if m == 1:
                    break
                # Record m mod 2^k BEFORE the Syracuse step
                hist_s[c, m & (mod_s - 1)] += 1
                hist_m[c, m & (mod_m - 1)] += 1
                hist_l[c, m & (mod_l - 1)] += 1
                cnt[c] += 1
                tmp = 3 * m + 1
                while tmp % 2 == 0:
                    tmp //= 2
                m = tmp
    return hist_s.sum(axis=0), hist_m.sum(axis=0), hist_l.sum(axis=0), cnt.sum()


def report_residue(name, hist, mod, total, focus_residue=None, top_n=10):
    odd_idx = np.arange(1, mod, 2)
    odd_counts = hist[odd_idx].astype(np.float64)
    n_odd = len(odd_idx)
    expected_uniform = total / n_odd  # uniform over n_odd odd residues
    ratio = odd_counts / expected_uniform
    se_ratio = 1.0 / np.sqrt(odd_counts)  # rough SE of the ratio
    se_ratio[odd_counts == 0] = float("inf")

    # Even residues: m is always odd along Syracuse trajectory, sanity check
    even_idx = np.arange(0, mod, 2)
    even_total = int(hist[even_idx].sum())

    print(f"\n=== {name} (mod {mod}, {n_odd} odd residues) ===")
    print(f"  total samples: {int(odd_counts.sum()):,}    even-residue counts (should be 0): {even_total}")
    print(f"  expected uniform per odd residue: {expected_uniform:,.1f}")
    print(f"  mean ratio: {ratio.mean():.6f}  (should be 1.0)")
    print(f"  SD ratio across {n_odd} odd residues: {ratio.std():.6f}")
    print(f"  min ratio: {ratio.min():.4f} at m≡{odd_idx[ratio.argmin()]}")
    print(f"  max ratio: {ratio.max():.4f} at m≡{odd_idx[ratio.argmax()]}")

    # Top over- and under-represented
    order_desc = np.argsort(-ratio)
    order_asc  = np.argsort(ratio)
    print(f"  top {top_n} over-represented: ", end="")
    print(", ".join(f"m≡{odd_idx[i]}:{ratio[i]:.3f}" for i in order_desc[:top_n]))
    print(f"  top {top_n} under-represented: ", end="")
    print(", ".join(f"m≡{odd_idx[i]}:{ratio[i]:.3f}" for i in order_asc[:top_n]))

    if focus_residue is not None:
        # Find this residue in odd_idx
        fi = np.searchsorted(odd_idx, focus_residue)
        if fi < n_odd and odd_idx[fi] == focus_residue:
            r = ratio[fi]
            se = se_ratio[fi]
            rank = int(np.sum(ratio > r) + 1)
            print(f"  >> focus m≡{focus_residue}: ratio={r:.4f}  SE={se:.4f}  "
                  f"rank={rank}/{n_odd}  (z={(r-1)/se:+.2f})")
        else:
            print(f"  !! focus residue {focus_residue} not in odd-residue set?")
    return ratio, odd_idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N_start", type=int, default=100_000_000)
    ap.add_argument("--T", type=int, default=200)
    ap.add_argument("--N_max", type=int, default=10**9)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    out_dir = here.parent / "experiments_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    starts = 2 * rng.integers(1, args.N_max // 2, size=args.N_start, dtype=np.int64) + 1

    print(f"[setup] N_start={args.N_start:,}  T={args.T}  N_max={args.N_max:,}", flush=True)
    print(f"[run] collecting m mod {{32, 2048, 131072}} along trajectories ...", flush=True)
    t0 = time.perf_counter()
    hist_s, hist_m, hist_l, total = collect_m_residues(starts, args.T, MOD_SMALL, MOD_MID, MOD_LARGE)
    t = time.perf_counter() - t0
    print(f"        done in {t:.1f}s  ({total:,} samples)", flush=True)

    ratio_s, idx_s = report_residue("m mod 32  (v=4 mechanism)",
                                    hist_s, MOD_SMALL, total, focus_residue=5, top_n=10)
    ratio_m, idx_m = report_residue("m mod 2048 (v=10 mechanism)",
                                    hist_m, MOD_MID, total, focus_residue=341, top_n=10)
    ratio_l, idx_l = report_residue("m mod 131072 (v=16 mechanism)",
                                    hist_l, MOD_LARGE, total, focus_residue=21845, top_n=10)

    # Save full histograms for downstream analysis
    import polars as pl
    out_s = pl.DataFrame({"residue": idx_s, "count": hist_s[idx_s], "ratio": ratio_s})
    out_m = pl.DataFrame({"residue": idx_m, "count": hist_m[idx_m], "ratio": ratio_m})
    # mod 131072 has 65536 odd residues; CSV would be a few MB, save it
    out_l = pl.DataFrame({"residue": idx_l, "count": hist_l[idx_l], "ratio": ratio_l})
    out_s.write_csv(out_dir / f"27_m_mod32_Nstart{args.N_start}_T{args.T}.csv")
    out_m.write_csv(out_dir / f"27_m_mod2048_Nstart{args.N_start}_T{args.T}.csv")
    out_l.write_csv(out_dir / f"27_m_mod131072_Nstart{args.N_start}_T{args.T}.csv")
    print(f"\n[save] {out_dir}/27_m_mod{{32,2048,131072}}_Nstart{args.N_start}_T{args.T}.csv")


if __name__ == "__main__":
    main()
