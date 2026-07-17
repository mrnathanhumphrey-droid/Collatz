"""C3: direct band-l¹ analysis on the dangerous band D_{r,t}(η).

Per Kalafatelis 2026 §4.5:
  q := 3^{r+1}, N_r := 2·3^{r-1}, ω_r := 2^{N_r} mod q (primitive cube root of unity)
  h_{r,ℓ}(j) := e_q(ω_r^ℓ · 2^j),  j ∈ Z/N_r,  ℓ ∈ {0, 1, 2}
  ĥ_{r,ℓ}(m) := Σ_{j=0}^{N_r-1} h_{r,ℓ}(j) e^{-2πi mj/N_r}     (forward DFT)
  m_r,t(θ) := Σ_{b=1}^{N_r} 2^{-b} e^{ib(t-θ)} = (1/2)e^{i(t-θ)} (1−(1/2 e^{i(t-θ)})^{N_r})/(1 − 1/2 e^{i(t-θ)})
  D_{r,t}(η) := {m ∈ Z/N_r : |m_r,t(2πm/N_r)| > 1−η}    (dangerous band)

Eq 190 conjecture: N_r^{-1/2} · ‖ĥ_{r,ℓ}‖_{ℓ¹(D_{r,t}(η))} ≪ η^{1/2+δ} for some δ > 0.

This script measures the LHS at r ∈ {6, 8, 10, 12, 14}, ℓ ∈ {0, 1, 2}, t ∈ {0, π/2, π},
η ∈ {1/16, 1/8, 1/4, 1/2}, and fits its scaling vs N_r at fixed η to determine whether
the LHS is bounded uniformly in r (eq 190 holds at observed r) or grows like √N_r
(saturated, eq 190 unprovable via this route).

Reference scalings:
  - Trivial (no inter-m cancellation): ‖ĥ‖_{ℓ¹(D)} ≈ |D| · sup|ĥ| ≈ √η·N_r · √N_r = √η·N_r^{3/2}
    → N_r^{-1/2}·this = √η·N_r  (slope vs log N_r: 1.0)
  - Cauchy-Schwarz: ‖ĥ‖_{ℓ¹(D)} ≤ √|D| · ‖ĥ‖_{ℓ²(D)} ≤ √|D|·N_r ≈ η^{1/4}·N_r^{3/2}
    → N_r^{-1/2}·this = η^{1/4}·N_r  (slope 1.0)
  - Square-root cancellation in m within band: ‖ĥ‖_{ℓ¹(D)} ≈ √|D| · √N_r·constant
    → N_r^{-1/2}·this = √(√η·N_r) · √N_r/√N_r·constant ∝ η^{1/4}·√N_r·constant. Slope 0.5.
  - Eq 190 (bounded in r): N_r^{-1/2}·‖ĥ‖_{ℓ¹(D)} = O(η^{1/2+δ}). Slope 0.0.

So the measured slope distinguishes:
  slope ≈ 1.0  →  no inter-m cancellation in band (path closed)
  slope ≈ 0.5  →  some cancellation but not enough (intermediate)
  slope ≈ 0.0  →  full cancellation (eq 190 holds, path open)
"""
import sys, os, math, time, csv
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

PI = math.pi


def compute_h(r: int, ell: int) -> np.ndarray:
    """h_{r,ℓ}(j) = e_q(ω_r^ℓ · 2^j) for j = 0..N_r-1, returned as complex128."""
    q = 3 ** (r + 1)
    N_r = 2 * 3 ** (r - 1)
    omega_r = pow(2, N_r, q)
    omega_r_ell = pow(omega_r, ell, q)
    out = np.empty(N_r, dtype=np.complex128)
    inv_q = 2 * PI / q
    pow2 = 1
    for j in range(N_r):
        phase = (omega_r_ell * pow2) % q
        out[j] = math.cos(inv_q * phase) + 1j * math.sin(inv_q * phase)
        pow2 = (pow2 * 2) % q
    return out


def compute_h_hat(r: int, ell: int) -> np.ndarray:
    """ĥ_{r,ℓ}(m) for m = 0..N_r-1 via forward DFT (numpy convention)."""
    h = compute_h(r, ell)
    return np.fft.fft(h)  # numpy fft uses e^{-2πi mj/N} which matches Kalafatelis (170)


def multiplier_modulus(r: int, t: float, m_arr: np.ndarray) -> np.ndarray:
    """|m_r,t(2πm/N_r)| for m in m_arr.

    m_r,t(θ) = (1/2)e^{i(t-θ)} (1 − (z)^{N_r}) / (1 − z), z = (1/2)e^{i(t-θ)}.
    """
    N_r = 2 * 3 ** (r - 1)
    theta = 2 * PI * m_arr / N_r
    phi = t - theta  # length |m_arr|
    # z = (1/2) e^{iφ}
    z = 0.5 * np.exp(1j * phi)
    # numerator: (1/2) e^{iφ} (1 − z^{N_r})
    z_pow = z ** N_r
    num = z * (1 - z_pow)
    den = 1 - z
    # avoid division by zero (when φ = 0, z = 1/2, den = 1/2, OK)
    return np.abs(num / den)


def main():
    out_csv = r"C:\Collatz\band_l1_data.csv"
    fieldnames = ["r", "N_r", "ell", "t_label", "t_val", "eta", "band_size",
                  "sum_l1", "metric", "metric_log_eta", "metric_div_sqrt_Nr",
                  "elapsed_s"]
    print("="*88)
    print("# C3: Band-l¹ analysis of ĥ_{r,ℓ} on dangerous band D_{r,t}(η)")
    print("# metric := N_r^{-1/2} · ‖ĥ‖_{ℓ¹(D)};  eq 190 needs metric ≪ η^{1/2+δ}")
    print("="*88)
    print()

    rows = []
    rs_to_run = [6, 8, 10, 12, 14]
    ells_to_run = [0, 1, 2]
    ts = [(0.0, "0"), (PI/2, "π/2"), (PI, "π")]
    etas = [1/2, 1/4, 1/8, 1/16]

    for r in rs_to_run:
        N_r = 2 * 3 ** (r - 1)
        sqrt_Nr = math.sqrt(N_r)
        t0 = time.time()
        # Per-ℓ FFT once
        h_hats = {}
        for ell in ells_to_run:
            h_hats[ell] = compute_h_hat(r, ell)

        # m_array for multiplier evaluation
        m_arr = np.arange(N_r, dtype=np.int64)

        for (t_val, t_label) in ts:
            mod = multiplier_modulus(r, t_val, m_arr)
            for eta in etas:
                threshold = 1 - eta
                D_mask = mod > threshold
                D_size = int(D_mask.sum())
                if D_size == 0:
                    continue
                for ell in ells_to_run:
                    h_hat = h_hats[ell]
                    abs_hh = np.abs(h_hat)
                    sum_l1 = float(abs_hh[D_mask].sum())
                    metric = sum_l1 / sqrt_Nr
                    # log_eta-form: how does metric scale with η?
                    metric_log_eta = math.log(metric) / math.log(eta) if eta < 1 and metric > 0 else 0.0
                    rows.append({
                        "r": r, "N_r": N_r, "ell": ell,
                        "t_label": t_label, "t_val": t_val,
                        "eta": eta, "band_size": D_size,
                        "sum_l1": sum_l1, "metric": metric,
                        "metric_log_eta": metric_log_eta,
                        "metric_div_sqrt_Nr": metric / sqrt_Nr,
                        "elapsed_s": "",
                    })
        elapsed = time.time() - t0
        print(f"r={r:>3}  N_r={N_r:>8}  elapsed {elapsed:.2f}s  ({len(ts)*len(etas)*len(ells_to_run)} sweeps)")

    # Write CSV
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"\n[csv written: {out_csv}]")
    print()

    # ----------------------------------------------------------------------
    # Per-η scaling fit: log(metric) vs log(N_r), at fixed (ℓ, t, η)
    # ----------------------------------------------------------------------
    print("="*88)
    print("# Per-η scaling: log(metric) = α + β · log(N_r), fixed (ℓ, t)")
    print("#   β ≈ 0.0 → eq 190 holds (bounded uniformly in r)")
    print("#   β ≈ 0.5 → partial cancellation (intermediate)")
    print("#   β ≈ 1.0 → no cancellation (saturated, path closed)")
    print("="*88)
    print()

    # Group by (ell, t_label, eta), fit slope across r
    by_key = {}
    for row in rows:
        key = (row["ell"], row["t_label"], row["eta"])
        by_key.setdefault(key, []).append((row["r"], row["N_r"], row["metric"]))

    print(f"  {'ℓ':>2} {'t':>6} {'η':>8} {'β (slope)':>11} {'α (intercept)':>14} "
          f"{'R²':>7} {'classification':>40}")

    for key in sorted(by_key.keys()):
        ell, t_label, eta = key
        pts = by_key[key]
        if len(pts) < 3:
            continue
        xs = np.array([math.log(p[1]) for p in pts])
        ys = np.array([math.log(p[2]) for p in pts])
        x_m = xs.mean(); y_m = ys.mean()
        ssxy = ((xs - x_m) * (ys - y_m)).sum()
        ssxx = ((xs - x_m) ** 2).sum()
        beta = ssxy / ssxx if ssxx > 0 else 0
        alpha = y_m - beta * x_m
        y_pred = alpha + beta * xs
        ss_res = ((ys - y_pred) ** 2).sum()
        ss_tot = ((ys - y_m) ** 2).sum()
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

        if abs(beta) < 0.15:
            cls = "BOUNDED — eq 190 holds at observed r"
        elif beta < 0.4:
            cls = "intermediate (partial cancellation)"
        elif beta < 0.7:
            cls = "≈ √N_r partial saturation"
        elif beta <= 1.05:
            cls = "saturated (no inter-m cancellation in band)"
        else:
            cls = f"FAST GROWTH (β > 1)"
        print(f"  {ell:>2} {t_label:>6} {eta:>8.4f} {beta:>+11.4f} {alpha:>+14.4f} "
              f"{r2:>7.4f}  {cls}")

    print()
    print("# Summary: aggregate β over all (ℓ, t) per η:")
    by_eta = {}
    for key, pts in by_key.items():
        eta = key[2]
        if len(pts) < 3:
            continue
        xs = np.array([math.log(p[1]) for p in pts])
        ys = np.array([math.log(p[2]) for p in pts])
        x_m = xs.mean(); y_m = ys.mean()
        ssxy = ((xs - x_m) * (ys - y_m)).sum()
        ssxx = ((xs - x_m) ** 2).sum()
        beta = ssxy / ssxx if ssxx > 0 else 0
        by_eta.setdefault(eta, []).append(beta)

    print(f"  {'η':>10} {'mean β':>10} {'std β':>10} {'min β':>10} {'max β':>10} {'n':>4}")
    for eta in sorted(by_eta.keys(), reverse=True):
        bs = by_eta[eta]
        print(f"  {eta:>10.4f} {np.mean(bs):>10.4f} {np.std(bs):>10.4f} "
              f"{np.min(bs):>10.4f} {np.max(bs):>10.4f} {len(bs):>4}")


if __name__ == "__main__":
    main()
