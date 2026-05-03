"""
Experiment 46 — test K_eff(N) residuals from K_h against algebraic indexings.

Discriminates monotone-asymptotic / oscillatory / arithmetic / limit-cycle
hypotheses for the K_eff(N) residual structure by checking whether residuals
fall on a smooth curve / cluster by class when indexed by various algebraic
properties of N.

Six data points from exp 45:
  N=2²⁰ K_eff=11.33  resid=-0.904
  N=2²² K_eff=10.33  resid=+0.094
  N=2²³ K_eff=10.05  resid=+0.376
  N=2²⁴ K_eff= 9.97  resid=+0.453
  N=2²⁵ K_eff= 9.96  resid=+0.463
  N=2²⁷ K_eff=10.13  resid=+0.298

Indexings (per prompt):
  1.  log₂(N) mod log₂(3)
  2.  log₂(N) mod log₂(4/3)
  3.  log₂(N) mod 1                — sanity (=0 for all)
  4.  log(N)  mod 1
  5.  k mod 3
  6.  k mod 6
  7.  N mod 3                       — N=2^k mod 3 = 2 if k odd else 1
  8.  N mod 5
  9.  N mod 7
  10. k itself                       — sanity (monotone if asymptotic)
  11. k mod log₂(3)                  — same as #1 since N=2^k
  12. gcd(k, 12)

Quantitative discriminators:
  Continuous (1, 2, 4, 11): cubic / sinusoid fit, R²
  Discrete  (5, 6, 7, 8, 9, 12): between-class / within-class variance ratio

Decisive criteria (per prompt):
  - All R² < 0.5 and all variance ratios < 2 → no algebraic structure
  - log₂(N) mod log₂(3) lights up → limit cycle (Diophantine)
  - Residue-class indexing lights up → arithmetic-function-like
  - Multiple partial → mixed / ambiguous

Six data points is THIN. If anything > threshold, flag for follow-up at
N=2²¹, 2²⁶, 2²⁸.

Usage:
    python 46_K_eff_residual_indexings.py
"""
import sys
from math import gcd
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

K_H = 3.0 / (np.log(4.0) - np.log(3.0))


def main():
    here = Path(__file__).resolve().parent
    out_dir = here.parent / "outputs" / "exp46"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Data
    ks = np.array([20, 22, 23, 24, 25, 27])
    K_eff = np.array([11.3319, 10.3337, 10.0523, 9.9749, 9.9649, 10.1301])
    resid = K_H - K_eff
    Ns = 2 ** ks
    log_N = ks * np.log(2)

    print(f"K_h = {K_H:.6f}")
    print(f"\nData:")
    print(f"  {'k':>3} {'N':>14} {'K_eff':>9} {'resid':>9}")
    for i in range(len(ks)):
        print(f"  {ks[i]:>3} {Ns[i]:>14,d} {K_eff[i]:>9.4f} {resid[i]:>+9.4f}")

    # Indexings
    LOG2_3 = np.log2(3)
    LOG2_4_3 = 2 - LOG2_3  # log₂(4/3)

    indexings = {
        "1_phase_log23":   ("continuous", ks % LOG2_3,
                            "log₂(N) mod log₂(3) ≈ 1.585"),
        "2_phase_log43":   ("continuous", ks % LOG2_4_3,
                            "log₂(N) mod log₂(4/3) ≈ 0.415"),
        "3_log2N_mod1":    ("sanity", ks.astype(float) % 1,
                            "log₂(N) mod 1 (sanity, =0)"),
        "4_logN_mod1":     ("continuous", log_N % 1,
                            "log(N) mod 1"),
        "5_k_mod3":        ("discrete", ks % 3,
                            "k mod 3"),
        "6_k_mod6":        ("discrete", ks % 6,
                            "k mod 6"),
        "7_N_mod3":        ("discrete", np.array([(pow(2, int(k), 3)) for k in ks]),
                            "N mod 3 (= 1 if k even, 2 if k odd)"),
        "8_N_mod5":        ("discrete", np.array([(pow(2, int(k), 5)) for k in ks]),
                            "N mod 5"),
        "9_N_mod7":        ("discrete", np.array([(pow(2, int(k), 7)) for k in ks]),
                            "N mod 7"),
        "10_k":            ("continuous", ks.astype(float),
                            "k itself (monotone if asymptotic)"),
        "11_k_mod_log23":  ("continuous", ks % LOG2_3,
                            "k mod log₂(3) (same as #1)"),
        "12_gcd_k_12":     ("discrete", np.array([gcd(int(k), 12) for k in ks]),
                            "gcd(k, 12)"),
    }

    def cubic_R2(x, y):
        if len(x) < 4:
            return float("nan"), None
        if (x.max() - x.min()) < 1e-12:
            return float("nan"), None
        coef = np.polyfit(x, y, 3)
        pred = np.polyval(coef, x)
        ss_res = ((y - pred) ** 2).sum()
        ss_tot = ((y - y.mean()) ** 2).sum()
        if ss_tot < 1e-20:
            return 1.0, coef
        r2 = 1 - ss_res / ss_tot
        return r2, coef

    def linear_R2(x, y):
        if len(x) < 2:
            return float("nan"), None
        if (x.max() - x.min()) < 1e-12:
            return float("nan"), None
        coef = np.polyfit(x, y, 1)
        pred = np.polyval(coef, x)
        ss_res = ((y - pred) ** 2).sum()
        ss_tot = ((y - y.mean()) ** 2).sum()
        if ss_tot < 1e-20:
            return 1.0, coef
        r2 = 1 - ss_res / ss_tot
        return r2, coef

    def sinusoid_R2(x, y, period_candidates):
        """Try sinusoidal fits at given periods, return best R² and period."""
        best_r2 = -np.inf
        best_period = None
        best_params = None
        for period in period_candidates:
            omega = 2 * np.pi / period
            X = np.column_stack([np.cos(omega * x), np.sin(omega * x), np.ones_like(x)])
            try:
                coef, *_ = np.linalg.lstsq(X, y, rcond=None)
                pred = X @ coef
                ss_res = ((y - pred) ** 2).sum()
                ss_tot = ((y - y.mean()) ** 2).sum()
                if ss_tot < 1e-20:
                    r2 = 1.0
                else:
                    r2 = 1 - ss_res / ss_tot
                if r2 > best_r2:
                    best_r2 = r2
                    best_period = period
                    best_params = coef
            except:
                pass
        return best_r2, best_period, best_params

    def variance_ratio(classes, values):
        unique = np.unique(classes)
        if len(unique) < 2:
            return float("nan"), float("nan"), float("nan")
        within_var_acc = []
        class_means = []
        class_n = []
        for c in unique:
            mask = classes == c
            n_c = int(mask.sum())
            if n_c == 0:
                continue
            v = values[mask]
            class_means.append(v.mean())
            class_n.append(n_c)
            if n_c >= 2:
                within_var_acc.append(((v - v.mean()) ** 2).sum())
        # Within-class variance (pooled)
        n_total = sum(class_n)
        if n_total - len(unique) < 1:
            within_var = float("nan")
        else:
            within_var = sum(within_var_acc) / max(n_total - len(unique), 1)
        # Between-class variance
        grand_mean = sum(m * n for m, n in zip(class_means, class_n)) / n_total
        between_var = sum(n * (m - grand_mean) ** 2 for m, n in zip(class_means, class_n)) / max(len(unique) - 1, 1)
        if within_var > 0:
            ratio = between_var / within_var
        else:
            ratio = float("inf")
        return between_var, within_var, ratio

    print(f"\n{'=' * 70}")
    print(f"FULL DATA (n=6)")
    print(f"{'=' * 70}")
    print(f"\n{'idx':>5} {'kind':>11} {'description':>40} {'metric':>8}  notes")
    print(f"{'-' * 90}")

    full_results = []
    for name, (kind, idx_vals, desc) in indexings.items():
        if kind in ("continuous", "sanity"):
            r2_lin, _ = linear_R2(idx_vals, resid)
            r2_cubic, _ = cubic_R2(idx_vals, resid)
            # Try sinusoidal too
            if "phase_log23" in name or "k_mod_log23" in name:
                periods = [LOG2_3]
            elif "phase_log43" in name:
                periods = [LOG2_4_3]
            elif "logN_mod1" in name:
                periods = [1.0]
            else:
                periods = [LOG2_3, LOG2_4_3, 2.0, 3.0, 6.0]
            r2_sin, period_sin, params_sin = sinusoid_R2(idx_vals, resid, periods)
            print(f"{name:>5} {kind:>11} {desc:>40} "
                  f"R²_lin={r2_lin:>+5.2f} R²_cub={r2_cubic:>+5.2f} R²_sin={r2_sin:>+5.2f}  "
                  f"sin@period={period_sin:.3f}" if period_sin else
                  f"{name:>5} {kind:>11} {desc:>40} "
                  f"R²_lin={r2_lin:>+5.2f} R²_cub={r2_cubic:>+5.2f} R²_sin={r2_sin:>+5.2f}")
            full_results.append((name, kind, desc, r2_lin, r2_cubic, r2_sin))
        else:  # discrete
            bv, wv, ratio = variance_ratio(idx_vals, resid)
            print(f"{name:>5} {kind:>11} {desc:>40} "
                  f"between={bv:.4f}  within={wv:.4f}  ratio={ratio:.2f}  "
                  f"classes={sorted(set(idx_vals.tolist()))}")
            full_results.append((name, kind, desc, bv, wv, ratio))

    # --- Repeat WITHOUT 2²⁰ since it's a clear boundary outlier ---
    print(f"\n{'=' * 70}")
    print(f"WITHOUT 2²⁰ (excluding boundary-regime outlier; n=5)")
    print(f"{'=' * 70}")
    mask = ks != 20
    ks2 = ks[mask]
    resid2 = resid[mask]
    log_N2 = log_N[mask]

    indexings2 = {
        "1_phase_log23":   ("continuous", ks2 % LOG2_3, "log₂(N) mod log₂(3)"),
        "2_phase_log43":   ("continuous", ks2 % LOG2_4_3, "log₂(N) mod log₂(4/3)"),
        "4_logN_mod1":     ("continuous", log_N2 % 1, "log(N) mod 1"),
        "5_k_mod3":        ("discrete", ks2 % 3, "k mod 3"),
        "7_N_mod3":        ("discrete", np.array([(pow(2, int(k), 3)) for k in ks2]),
                            "N mod 3"),
        "8_N_mod5":        ("discrete", np.array([(pow(2, int(k), 5)) for k in ks2]),
                            "N mod 5"),
        "9_N_mod7":        ("discrete", np.array([(pow(2, int(k), 7)) for k in ks2]),
                            "N mod 7"),
        "10_k":            ("continuous", ks2.astype(float), "k itself"),
        "12_gcd_k_12":     ("discrete", np.array([gcd(int(k), 12) for k in ks2]),
                            "gcd(k, 12)"),
    }

    no2to20_results = []
    for name, (kind, idx_vals, desc) in indexings2.items():
        if kind == "continuous":
            r2_lin, _ = linear_R2(idx_vals, resid2)
            r2_cubic, _ = cubic_R2(idx_vals, resid2)
            if "log23" in name:
                periods = [LOG2_3]
            elif "log43" in name:
                periods = [LOG2_4_3]
            elif "logN_mod1" in name:
                periods = [1.0]
            else:
                periods = [2.0, 3.0, 4.0, 5.0, 6.0]
            r2_sin, period_sin, _ = sinusoid_R2(idx_vals, resid2, periods)
            print(f"  {name:>5}: {desc}  "
                  f"R²_lin={r2_lin:>+5.2f} R²_cub={r2_cubic:>+5.2f} R²_sin={r2_sin:>+5.2f}")
            no2to20_results.append((name, kind, desc, r2_lin, r2_cubic, r2_sin))
        else:
            bv, wv, ratio = variance_ratio(idx_vals, resid2)
            print(f"  {name:>5}: {desc}  "
                  f"between={bv:.4f}  within={wv:.4f}  ratio={ratio:.2f}  "
                  f"classes={sorted(set(idx_vals.tolist()))}")
            no2to20_results.append((name, kind, desc, bv, wv, ratio))

    # --- Print residuals sorted by each indexing for visual inspection ---
    print(f"\n{'=' * 70}")
    print(f"RESIDUALS SORTED BY INDEXING (full data)")
    print(f"{'=' * 70}")
    for name, (kind, idx_vals, desc) in list(indexings.items())[:12]:
        if kind == "sanity":
            continue
        order = np.argsort(idx_vals)
        print(f"\n{name}: {desc}")
        print(f"  {'idx':>10} {'k':>4} {'resid':>9}")
        for o in order:
            print(f"  {idx_vals[o]:>10.4f} {ks[o]:>4} {resid[o]:>+9.4f}")

    # --- Plots: 12-panel grid ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(3, 4, figsize=(16, 11))
        axes = axes.flatten()

        for ax, (name, (kind, idx_vals, desc)) in zip(axes, indexings.items()):
            ax.scatter(idx_vals, resid, s=80, c="tab:blue", edgecolor="black", zorder=3)
            for x, y, k in zip(idx_vals, resid, ks):
                ax.annotate(f"k={k}", (x, y), textcoords="offset points",
                            xytext=(5, 5), fontsize=8)
            ax.axhline(0, color="black", linewidth=0.5)
            ax.set_title(f"{name}\n{desc}", fontsize=9)
            ax.set_xlabel("indexing", fontsize=8)
            ax.set_ylabel("K_h - K_eff", fontsize=8)
            ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / "plot_grid_12.png", dpi=110)
        plt.close()
        print(f"\n[plot] grid saved to {out_dir}")
    except Exception as e:
        print(f"  [warn] plotting failed: {e}")

    # --- Summary verdict ---
    print(f"\n{'=' * 70}")
    print(f"VERDICT")
    print(f"{'=' * 70}")
    flagged_full = []
    for name, kind, desc, *metrics in full_results:
        if kind in ("continuous",):
            r2_lin, r2_cubic, r2_sin = metrics
            best = max(r2_lin, r2_cubic, r2_sin)
            if best > 0.7:
                flagged_full.append((name, "STRONG", best, desc))
            elif best > 0.5:
                flagged_full.append((name, "moderate", best, desc))
        elif kind == "discrete":
            _, _, ratio = metrics
            if ratio > 2:
                flagged_full.append((name, "STRONG", ratio, desc))
            elif ratio > 1.5:
                flagged_full.append((name, "moderate", ratio, desc))

    if not flagged_full:
        print(f"  Full data: NO indexing exceeds R²>0.5 or var-ratio>1.5")
    else:
        for name, level, val, desc in flagged_full:
            print(f"  Full data: {name} ({level}, val={val:.2f}): {desc}")

    flagged_no20 = []
    for name, kind, desc, *metrics in no2to20_results:
        if kind == "continuous":
            r2_lin, r2_cubic, r2_sin = metrics
            best = max(r2_lin, r2_cubic, r2_sin)
            if best > 0.7:
                flagged_no20.append((name, "STRONG", best, desc))
            elif best > 0.5:
                flagged_no20.append((name, "moderate", best, desc))
        elif kind == "discrete":
            _, _, ratio = metrics
            if ratio > 2:
                flagged_no20.append((name, "STRONG", ratio, desc))
            elif ratio > 1.5:
                flagged_no20.append((name, "moderate", ratio, desc))

    if not flagged_no20:
        print(f"  Without 2²⁰: NO indexing exceeds R²>0.5 or var-ratio>1.5")
    else:
        for name, level, val, desc in flagged_no20:
            print(f"  Without 2²⁰: {name} ({level}, val={val:.2f}): {desc}")

    print(f"\nNote: 6 (or 5) data points is genuinely thin. Any 'strong' signal here")
    print(f"      should be confirmed at additional N values (2²¹, 2²⁶, 2²⁸).")
    print(f"\n[done]")


if __name__ == "__main__":
    main()
