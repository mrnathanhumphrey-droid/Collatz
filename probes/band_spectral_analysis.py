"""C3-extension: spectral decomposition of ĥ_{r,ℓ} restricted to dangerous band.

For each (r, ℓ, t, η):
  1. Compute ĥ_{r,ℓ} via FFT
  2. Build dangerous band D_{r,t}(η) per Kalafatelis Prop 22
  3. Restrict ĥ to J = D, treat as 1D sequence f := ĥ|_J in m-sorted order
  4. DFT of restriction: g(k) := Σ_j f(j) · e^{-2πi kj/|J|}
  5. Compute spectral metrics:
     - max_concentration = max|g|² / Σ|g|²
     - entropy H = -Σ p log p,  p = |g|²/Σ|g|²,  H_norm = H / log|J|
     - r_eff = #modes for 90% energy
     - low-freq mass = mass in |k| ≤ |J|/8 (folded for negative freqs)

Decision rule:
  α (smooth, exploitable): r_eff < |J|/4 AND low_freq_mass > 0.6
  γ (white-noise rough):    r_eff ≥ 0.8|J| AND H_norm > 0.95
  β (mixed):                in between
  δ (varies):               outcome differs across cells at same r
"""
import sys, os, math, time, csv
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"C:\Collatz")

from band_l1_analysis import compute_h, compute_h_hat, multiplier_modulus

PI = math.pi


def spectral_metrics(h_hat: np.ndarray, D_mask: np.ndarray) -> dict:
    """Compute spectral metrics for ĥ restricted to band D."""
    J_indices = np.where(D_mask)[0]
    if len(J_indices) < 2:
        return None
    f = h_hat[J_indices]
    g = np.fft.fft(f)
    g_sq = np.abs(g) ** 2
    energy = g_sq.sum()
    if energy == 0:
        return None

    p = g_sq / energy
    max_p = float(p.max())

    # Entropy (Shannon)
    p_nz = p[p > 1e-15]
    H = float(-(p_nz * np.log(p_nz)).sum())
    log_J = math.log(len(f))
    H_norm = H / log_J if log_J > 0 else 1.0

    # Effective rank (modes for 90% energy)
    sorted_p = np.sort(p)[::-1]
    cum = np.cumsum(sorted_p)
    r_eff = int((cum < 0.9).sum()) + 1

    # Low-frequency mass: |k| ≤ |J|/8 (FFT folds negatives at end)
    cutoff = max(1, len(f) // 8)
    lf_mass = float(g_sq[:cutoff].sum() + g_sq[-cutoff:].sum()) / energy
    if cutoff * 2 >= len(f):
        # Avoid double-counting if cutoff covers full range
        lf_mass = 1.0

    return {
        "J_size": len(f),
        "max_concentration": max_p,
        "uniform_baseline": 1.0 / len(f),
        "concentration_ratio": max_p * len(f),  # 1 = uniform, large = concentrated
        "entropy": H,
        "entropy_normalized": H_norm,
        "r_eff": r_eff,
        "r_eff_ratio": r_eff / len(f),
        "low_freq_mass": lf_mass,
    }


def classify_cell(metrics: dict) -> str:
    """Classify cell as α/β/γ based on metrics."""
    if metrics is None:
        return "skip"
    r_eff_r = metrics["r_eff_ratio"]
    lf = metrics["low_freq_mass"]
    H_norm = metrics["entropy_normalized"]
    is_alpha = (r_eff_r < 0.25) and (lf > 0.6)
    is_gamma = (r_eff_r >= 0.8) and (H_norm > 0.95)
    if is_alpha:
        return "α"
    if is_gamma:
        return "γ"
    return "β"


def main():
    out_csv = r"C:\Collatz\band_spectral_data.csv"
    fieldnames = ["r", "N_r", "ell", "t_label", "t_val", "eta",
                  "J_size", "max_concentration", "concentration_ratio",
                  "entropy", "entropy_normalized", "r_eff", "r_eff_ratio",
                  "low_freq_mass", "classification"]

    print("="*90)
    print("# Spectral decomposition of ĥ_{r,ℓ} restricted to dangerous band D_{r,t}(η)")
    print("# Q: does inter-m structure admit smooth-weight cancellation that ‖·‖_{ℓ¹} doesn't capture?")
    print("="*90)
    print()

    rows = []
    rs_to_run = [6, 8, 10, 12, 14]
    ells_to_run = [0, 1, 2, 3]
    ts = [(1.0, "1"), (2.0, "2"), (3.0, "3")]
    etas = [0.5, 0.25]

    for r in rs_to_run:
        N_r = 2 * 3 ** (r - 1)
        t0 = time.time()
        # Compute ĥ for each ℓ
        h_hats = {}
        for ell in ells_to_run:
            h_hats[ell] = compute_h_hat(r, ell)

        m_arr = np.arange(N_r, dtype=np.int64)

        for (t_val, t_label) in ts:
            mod = multiplier_modulus(r, t_val, m_arr)
            for eta in etas:
                threshold = 1 - eta
                D_mask = mod > threshold
                D_size = int(D_mask.sum())
                if D_size < 4:
                    continue
                for ell in ells_to_run:
                    h_hat = h_hats[ell]
                    metrics = spectral_metrics(h_hat, D_mask)
                    if metrics is None:
                        continue
                    cls = classify_cell(metrics)
                    rows.append({
                        "r": r, "N_r": N_r, "ell": ell,
                        "t_label": t_label, "t_val": t_val, "eta": eta,
                        **metrics, "classification": cls,
                    })

        elapsed = time.time() - t0
        print(f"r={r:>3}  N_r={N_r:>8}  elapsed {elapsed:.2f}s  ({len(ells_to_run)*len(ts)*len(etas)} cells)")

    print()

    # Print per-(r, η) summary
    print("="*90)
    print("# Per-(r, η) summary: averaged over (ℓ, t)")
    print("="*90)
    print(f"  {'r':>3} {'eta':>6} {'|J|':>6} {'r_eff':>6} {'r_eff_r':>8} {'H_norm':>8} "
          f"{'lf_mass':>8} {'conc_ratio':>10} {'class_dist':>20}")

    by_r_eta = {}
    for row in rows:
        key = (row["r"], row["eta"])
        by_r_eta.setdefault(key, []).append(row)

    for key in sorted(by_r_eta.keys()):
        r, eta = key
        cells = by_r_eta[key]
        Jsize = np.mean([c["J_size"] for c in cells])
        r_eff = np.mean([c["r_eff"] for c in cells])
        r_eff_r = np.mean([c["r_eff_ratio"] for c in cells])
        H_norm = np.mean([c["entropy_normalized"] for c in cells])
        lf = np.mean([c["low_freq_mass"] for c in cells])
        conc = np.mean([c["concentration_ratio"] for c in cells])
        from collections import Counter
        class_counts = Counter(c["classification"] for c in cells)
        class_str = ", ".join(f"{k}:{v}" for k, v in sorted(class_counts.items()))
        print(f"  {r:>3} {eta:>6.3f} {Jsize:>6.1f} {r_eff:>6.1f} {r_eff_r:>8.3f} "
              f"{H_norm:>8.3f} {lf:>8.3f} {conc:>10.2f}  {class_str:>20}")

    print()

    # Aggregate verdict
    print("="*90)
    print("# Aggregate verdict")
    print("="*90)
    from collections import Counter
    cls_counts = Counter(row["classification"] for row in rows)
    total = sum(cls_counts.values())
    print(f"  Total cells: {total}")
    for cls in ["α", "β", "γ"]:
        cnt = cls_counts.get(cls, 0)
        pct = 100 * cnt / total if total > 0 else 0
        print(f"    {cls}: {cnt} ({pct:.0f}%)")

    if cls_counts.get("γ", 0) / total > 0.7:
        outcome = "γ (white-noise rough)"
        next_action = "Smooth completion (R78 path 2) loses motivation. Genuinely-novel direction needed."
    elif cls_counts.get("α", 0) / total > 0.7:
        outcome = "α (smooth, exploitable)"
        next_action = "Smooth completion (R78 path 2) becomes a concrete attack."
    elif cls_counts.get("β", 0) / total > 0.5:
        outcome = "β (mixed, partially smooth)"
        next_action = "Mode structure documented; could combine with Esscher tilt or other measure-changing tools."
    else:
        outcome = "δ (varies by cell)"
        next_action = "Smooth subset might admit partial closure; characterize."
    print(f"\n  OUTCOME: {outcome}")
    print(f"  NEXT: {next_action}")

    # CSV
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            # Filter to fieldnames only
            out_row = {k: row.get(k, "") for k in fieldnames}
            writer.writerow(out_row)
    print(f"\n[csv: {out_csv}]")

    # Scaling fit: r_eff_ratio vs r at fixed η — does it shrink as r grows (smoothing) or stay flat?
    print()
    print("="*90)
    print("# Scaling: how do r_eff_ratio and lf_mass evolve with r at fixed η?")
    print("="*90)
    for eta in etas:
        print(f"\n  η = {eta}:")
        print(f"    {'r':>3} {'mean r_eff_ratio':>17} {'mean lf_mass':>14} {'mean H_norm':>13}")
        for r in rs_to_run:
            cells = [row for row in rows if row["r"] == r and row["eta"] == eta]
            if not cells:
                continue
            r_eff_r = np.mean([c["r_eff_ratio"] for c in cells])
            lf = np.mean([c["low_freq_mass"] for c in cells])
            H = np.mean([c["entropy_normalized"] for c in cells])
            print(f"    {r:>3} {r_eff_r:>17.4f} {lf:>14.4f} {H:>13.4f}")


if __name__ == "__main__":
    main()
