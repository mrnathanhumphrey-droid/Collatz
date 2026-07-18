"""
build_epsilon_consolidated.py
Consolidate eps_k for audio rendering. No new compute — reads the existing
result_epsilon_11.csv (most complete: k=1..11 with signs) and emits:

  audio_data/epsilon_k_consolidated.csv   (k, S_k, |S_k - 7/15|, eps_k, sign)
  audio_data/epsilon_k_summary.txt         (sign / magnitude / ratio sequences,
                                            on-file recurrence coefficients)
"""
from __future__ import annotations

import os
import sys

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")

SRC_EPS = r"C:\Collatz\result_epsilon_11.csv"
SRC_RECUR = r"C:\Collatz\result_renormalization_recurrence_fits.csv"
OUT_DIR = r"C:\Collatz\audio_data"
OUT_CSV = os.path.join(OUT_DIR, "epsilon_k_consolidated.csv")
OUT_TXT = os.path.join(OUT_DIR, "epsilon_k_summary.txt")

SEVEN_FIFTEENTHS = 7.0 / 15.0  # 0.466666...

K_RANGE_LO = 2
K_RANGE_HI = 11

def main():
    eps = pd.read_csv(SRC_EPS)
    eps = eps[(eps["k"] >= K_RANGE_LO) & (eps["k"] <= K_RANGE_HI)].copy()
    eps["S_k"] = eps["eps_k"] + SEVEN_FIFTEENTHS
    eps["S_k_minus_seven_fifteenths"] = eps["abs_eps_k"]
    eps["epsilon_k"] = eps["eps_k"]
    out_csv = eps[["k", "S_k", "S_k_minus_seven_fifteenths", "epsilon_k", "sign"]].copy()
    out_csv.to_csv(OUT_CSV, index=False, float_format="%.15e")
    print(f"[csv] {OUT_CSV}")
    print(out_csv.to_string(index=False))
    print()

    sign_seq = " ".join(eps["sign"].tolist())
    mag_seq = ", ".join(f"{m:.6e}" for m in eps["abs_eps_k"])
    ratio_signed = ", ".join(
        ("--" if pd.isna(r) else f"{r:+.6f}") for r in eps["ratio_to_prev_signed"]
    )
    ratio_abs = ", ".join(
        ("--" if pd.isna(r) else f"{r:.6f}") for r in eps["abs_ratio_to_prev"]
    )

    recur = pd.read_csv(SRC_RECUR)

    with open(OUT_TXT, "w", encoding="utf-8") as fh:
        fh.write("Epsilon_k convergence sequence summary\n")
        fh.write("=======================================\n\n")
        fh.write(f"Source: {os.path.basename(SRC_EPS)} (k = 1..11; consolidating k = {K_RANGE_LO}..{K_RANGE_HI})\n")
        fh.write(f"Limit: 7/15 = {SEVEN_FIFTEENTHS:.15f}\n")
        fh.write(f"epsilon_k = S_k - 7/15  (S_k = stationary Plancherel mass at level k)\n\n")

        fh.write("k     eps_k                  |eps_k|              sign\n")
        for _, r in eps.iterrows():
            fh.write(f"{int(r['k']):>2d}    "
                     f"{r['eps_k']:+.15e}  {r['abs_eps_k']:.15e}  {r['sign']}\n")
        fh.write("\n")

        fh.write(f"Sign sequence (k={K_RANGE_LO}..{K_RANGE_HI}):  {sign_seq}\n\n")

        fh.write(f"Magnitude sequence |eps_k|:\n  {mag_seq}\n\n")

        fh.write(f"Ratio eps_{{k+1}} / eps_k (signed):\n  {ratio_signed}\n\n")

        fh.write(f"Ratio |eps_{{k+1}}/eps_k| (absolute):\n  {ratio_abs}\n\n")

        fh.write("On-file recurrence fits (result_renormalization_recurrence_fits.csv)\n")
        fh.write("--------------------------------------------------------------------\n")
        for _, r in recur.iterrows():
            fh.write(f"\nOrder {int(r['order'])}: n_eq={int(r['n_eq'])}, "
                     f"R^2={r['r2']:.6f}, ss_res={r['ss_res']:.4e}\n")
            fh.write(f"  alphas: {r['alphas (comma-sep)']}\n")
            fh.write(f"  roots:  {r['roots (comma-sep)']}\n")
        fh.write("\n")

        order3 = recur[recur["order"] == 3].iloc[0]
        fh.write("Order-3 highlight (per STATE.md):\n")
        fh.write(f"  Real root rho_slow = 0.826934  (the slow-mode rate)\n")
        fh.write(f"  Complex pair        = 0.154404 +/- 0.114255 i  "
                 f"(magnitude {(0.154404**2 + 0.114255**2)**0.5:.6f}, "
                 f"period 2*pi/atan2(0.114255, 0.154404) ~= "
                 f"{__import__('math').tau / __import__('math').atan2(0.114255, 0.154404):.4f} k-steps)\n")
        fh.write(f"  R^2 (on file)       = {order3['r2']:.6f}  "
                 f"(STATE.md cites ~1.00; on-file CSV has 0.797 — discrepancy worth verifying)\n\n")

        fh.write("Walk-back log (relevant to audio interpretation)\n")
        fh.write("-----------------------------------------------\n")
        fh.write("- Rate-1/2 envelope on |eps_n|*2^n looked stable ~0.04 at k=2..6;\n")
        fh.write("  jumped to 0.150 at k=7 (4x).  WALKED BACK 2026-05-05.\n")
        fh.write("- rho ~ 0.984 single complex-pair model fitted from k=2..10;\n")
        fh.write("  falsified at k=11 (eps_11 grew further to +1.50e-3).\n")
        fh.write("- Current best: order-3 linear recurrence with rho_slow ~ 0.827 (real)\n")
        fh.write("  plus damped complex pair near origin.\n\n")

        fh.write("Notes for audio mapping\n")
        fh.write("-----------------------\n")
        fh.write("- Sign flips (+/-) at k=2->3 and k=9->10 are the major structural events.\n")
        fh.write("- Magnitude minimum at k=9 (|eps_9|=7.5e-6, three orders below typical scale).\n")
        fh.write("- Magnitude bounce at k=10 (95.8x increase) is the largest ratio in the sequence.\n")

    print(f"[txt] {OUT_TXT}")


if __name__ == "__main__":
    main()
