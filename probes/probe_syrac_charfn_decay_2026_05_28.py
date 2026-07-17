"""
probe_syrac_charfn_decay_2026_05_28.py

FIRST HONEST LOOK at Tao's Prop 1.14 / 1.17 object directly.

Computes the EXACT Syracuse characteristic function
    S_chi(n)(xi) = E[ exp(-2*pi*i * xi * Syrac(Z/3^n) / 3^n) ]
for the actual n-Syracuse OFFSET random variable (Tao 1909.03562, Prop 7.1):

    Syrac(Z/3^n) = X_n = sum_{j=1}^{n} 3^{j-1} * 2^{-(a_1+...+a_j)}   mod 3^n,
    (a_1,...,a_n) iid Geom(2):  P(a=k) = 2^{-k}, k>=1.

NOTE: this is DIFFERENT from the project's pi_n (stationary distribution of the
ONE-STEP Syracuse Markov kernel K_n). The eps_k tower studies pi_n; Prop 1.14 is
about the n-fold offset X_n. This probe builds X_n directly.

Method: exact distribution of X_n on Z/3^n via the Horner backward recursion
    S_{n}   = U_n
    S_{j}   = U_j * (1 + 3*S_{j+1})        (j = n-1 .. 1)
    X_n     = S_1
Each S_j is supported on the units (Z/3^n)* (avoids multiples of 3).
U is a single Geom(2) step: P_U[2^{-a} mod 3^n] += 2^{-a}, truncated at A_MAX.

Then mu_hat = np.fft.fft(P_X) gives S_chi(n)(xi) for all xi at once.

Outputs:
  - table of max / mean / median |mu_hat| over xi not div by 3, vs n
  - comparison to n^{-A} (A=1,2,3) and geometric rho^n
  - the slowest-decaying xi (argmax) and its 3-adic valuation
  - L2 energy Sum_{xi != 0 mod 3} |mu_hat|^2  (analog of project's S_n, for X_n)
  - plots: probe_syrac_charfn_decay_2026_05_28.png
  - json dump
"""
from __future__ import annotations
import json, os, sys, time
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)
A_MAX = 100  # 2^-100 ~ 1e-30, effectively exact in float64

def v3(x):
    """3-adic valuation of a positive integer."""
    if x == 0:
        return None
    v = 0
    while x % 3 == 0:
        x //= 3
        v += 1
    return v

def syrac_offset_distribution(n):
    """Exact distribution of X_n = Syrac(Z/3^n) as a length-3^n probability vector."""
    N = 3 ** n
    inv2 = pow(2, -1, N)
    arange = np.arange(N)

    # single Geom(2) step distribution of U = 2^{-a} mod N (a>=1), as value-vector
    P_U = np.zeros(N, dtype=np.float64)
    p = inv2  # 2^{-1}
    for a in range(1, A_MAX + 1):
        P_U[p] += 2.0 ** (-a)
        p = (p * inv2) % N
    P_U /= P_U.sum()
    # support of U (units 2^{-a}); typically <= A_MAX distinct values
    u_support = np.nonzero(P_U)[0]
    u_weight = P_U[u_support]

    # S_n = U
    P_S = P_U.copy()
    for j in range(n - 1, 0, -1):
        # V = 1 + 3*S  (many-to-one mod 3^n: 3*s has period 3^{n-1})
        v_idx = (1 + 3 * arange) % N
        P_V = np.zeros(N, dtype=np.float64)
        np.add.at(P_V, v_idx, P_S)
        # S_j = U * V  (multiplicative convolution; loop over small U support)
        P_new = np.zeros(N, dtype=np.float64)
        for u, w in zip(u_support, u_weight):
            perm = (u * arange) % N   # bijection of Z/N (u is a unit)
            P_new[perm] += w * P_V
        P_S = P_new
    return P_S  # = P_{X_n}

def analyze(n):
    t0 = time.time()
    N = 3 ** n
    P = syrac_offset_distribution(n)
    # sanity: probability mass, support on units only
    mass = P.sum()
    mass_on_mult3 = P[::3].sum()
    mu = np.fft.fft(P)                  # mu[xi] = sum_y P[y] e^{-2pi i xi y / N}
    xi = np.arange(N)
    nontriv = (xi % 3 != 0)
    absmu = np.abs(mu[nontriv])
    xis = xi[nontriv]
    imax = int(np.argmax(absmu))
    xi_star = int(xis[imax])
    out = {
        "n": n, "N": N,
        "mass": float(mass),
        "mass_on_multiples_of_3": float(mass_on_mult3),
        "max_abs_muhat": float(absmu.max()),
        "mean_abs_muhat": float(absmu.mean()),
        "median_abs_muhat": float(np.median(absmu)),
        "L2_energy_nontrivial": float(np.sum(absmu ** 2)),
        "argmax_xi": xi_star,
        "argmax_xi_v3": v3(xi_star),
        "muhat_at_xi1": float(np.abs(mu[1])),
        "n_nontrivial_freqs": int(nontriv.sum()),
        "time_sec": time.time() - t0,
    }
    return out, absmu

def main():
    NS = list(range(1, 12))  # n = 1..11 ; N=3^11=177147 (fast)
    rows = []
    spectra = {}
    print("=" * 88)
    print("Tao Prop 1.17 object: |S_chi(n)(xi)| = |E exp(-2pi i xi Syrac(Z/3^n)/3^n)|, xi not div 3")
    print("=" * 88)
    print(f"{'n':>3} {'N=3^n':>8} {'max|mu|':>11} {'mean|mu|':>11} {'median|mu|':>11} "
          f"{'L2(xi!=0)':>10} {'argmax xi':>9} {'v3':>3} {'t(s)':>6}")
    for n in NS:
        o, absmu = analyze(n)
        rows.append(o)
        spectra[n] = absmu
        # sanity check
        assert abs(o["mass"] - 1.0) < 1e-9, o["mass"]
        assert o["mass_on_multiples_of_3"] < 1e-12, o["mass_on_multiples_of_3"]
        print(f"{o['n']:>3} {o['N']:>8} {o['max_abs_muhat']:>11.6f} {o['mean_abs_muhat']:>11.6f} "
              f"{o['median_abs_muhat']:>11.6f} {o['L2_energy_nontrivial']:>10.5f} "
              f"{o['argmax_xi']:>9} {str(o['argmax_xi_v3']):>3} {o['time_sec']:>6.2f}")

    # ---- decay-rate fits on max|mu| vs n ----
    ns = np.array([r["n"] for r in rows], dtype=float)
    mx = np.array([r["max_abs_muhat"] for r in rows], dtype=float)
    mean = np.array([r["mean_abs_muhat"] for r in rows], dtype=float)

    print("\n" + "=" * 88)
    print("DECAY DIAGNOSIS on max_xi |mu_hat(xi)|")
    print("=" * 88)
    # power-law fit log(max) = log C - A log n  (use n>=2)
    m = ns >= 2
    A_pow, logC = np.polyfit(np.log(ns[m]), np.log(mx[m]), 1)
    print(f"  power-law fit  max|mu| ~ C * n^(-A):   A = {-A_pow:.4f},  C = {np.exp(logC):.4f}")
    # geometric fit log(max) = log C + n log rho
    lr, logC2 = np.polyfit(ns[m], np.log(mx[m]), 1)
    print(f"  geometric fit  max|mu| ~ C * rho^n:    rho = {np.exp(lr):.5f},  C = {np.exp(logC2):.4f}")
    # successive ratios + successive power-law exponents
    print(f"\n  {'n':>3} {'max|mu|':>11} {'ratio':>9} {'local A (n^-A)':>15} {'local rho':>10}")
    for i in range(len(ns)):
        ratio = mx[i] / mx[i-1] if i > 0 else float('nan')
        locA = (-np.log(mx[i]/mx[i-1]) / np.log(ns[i]/ns[i-1])) if i > 0 else float('nan')
        locrho = ratio
        print(f"  {int(ns[i]):>3} {mx[i]:>11.6f} {ratio:>9.4f} {locA:>15.4f} {locrho:>10.5f}")

    # ---- plot ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # (a) semilog max & mean vs n with geometric overlay
        ax = axes[0, 0]
        ax.semilogy(ns, mx, "o-", label="max_xi |mu_hat|", color="C0")
        ax.semilogy(ns, mean, "s-", label="mean_xi |mu_hat|", color="C1")
        ax.semilogy(ns, np.exp(logC2) * np.exp(lr) ** ns, "--",
                    label=f"geom rho={np.exp(lr):.4f}", color="C3")
        ax.set_xlabel("n"); ax.set_ylabel("|S_chi(n)(xi)|  (log)")
        ax.set_title("(a) char-fn magnitude vs n  [semilog]"); ax.legend(); ax.grid(alpha=0.3)

        # (b) log-log max vs n with n^-A overlays
        ax = axes[0, 1]
        ax.loglog(ns, mx, "o-", label="max_xi |mu_hat|", color="C0")
        for A in (1, 2, 3):
            ax.loglog(ns[m], mx[m][0] * (ns[m] / ns[m][0]) ** (-A), ":",
                      label=f"n^(-{A})")
        ax.set_xlabel("n (log)"); ax.set_ylabel("max|mu_hat| (log)")
        ax.set_title("(b) power-law check  [log-log]"); ax.legend(); ax.grid(alpha=0.3, which="both")

        # (c) full spectrum of |mu_hat(xi)| for a few n (sorted descending)
        ax = axes[1, 0]
        for n in [3, 5, 7, 9, 11]:
            if n in spectra:
                s = np.sort(spectra[n])[::-1]
                ax.semilogy(np.linspace(0, 1, len(s)), s, label=f"n={n}")
        ax.set_xlabel("rank fraction of xi (xi not div 3)")
        ax.set_ylabel("|mu_hat(xi)| (log)")
        ax.set_title("(c) full sorted spectrum per n"); ax.legend(); ax.grid(alpha=0.3)

        # (d) local power-law exponent A and local ratio vs n
        ax = axes[1, 1]
        locA = [(-np.log(mx[i]/mx[i-1]) / np.log(ns[i]/ns[i-1])) for i in range(1, len(ns))]
        locrho = [mx[i]/mx[i-1] for i in range(1, len(ns))]
        ax.plot(ns[1:], locA, "o-", label="local A in n^(-A)", color="C2")
        ax.set_xlabel("n"); ax.set_ylabel("local power-law A", color="C2")
        ax.tick_params(axis="y", labelcolor="C2")
        ax2 = ax.twinx()
        ax2.plot(ns[1:], locrho, "s--", label="local ratio (rho)", color="C3")
        ax2.set_ylabel("local ratio max_n/max_{n-1}", color="C3")
        ax2.tick_params(axis="y", labelcolor="C3")
        ax.set_title("(d) is decay polynomial or geometric?"); ax.grid(alpha=0.3)

        fig.suptitle("Tao Prop 1.14/1.17: actual Syracuse char-fn decay (n-fold offset X_n)",
                     fontsize=13)
        fig.tight_layout()
        png = os.path.join(r"C:\Collatz", "probe_syrac_charfn_decay_2026_05_28.png")
        fig.savefig(png, dpi=110)
        print(f"\n  plot saved: {png}")
    except Exception as e:
        print(f"  [plot skipped: {e}]")

    out = {
        "params": {"A_MAX": A_MAX, "ns": NS},
        "rows": rows,
        "fits": {
            "powerlaw_A": float(-A_pow), "powerlaw_C": float(np.exp(logC)),
            "geometric_rho": float(np.exp(lr)), "geometric_C": float(np.exp(logC2)),
        },
    }
    jpath = os.path.join(OUTDIR, "probe_syrac_charfn_decay_2026_05_28.json")
    with open(jpath, "w") as f:
        json.dump(out, f, indent=2)
    print(f"  json saved: {jpath}")

if __name__ == "__main__":
    main()
