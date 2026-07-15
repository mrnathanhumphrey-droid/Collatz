"""
Collatz <-> turbulent-cascade bridge: generator facts + the codimension tension.

Two jobs:
  (1) BANK the durable positive: the Syracuse multiplier W = q/2^v (v ~ Geom(1/2)) is a
      legitimate log-geometric compound-Poisson cascade generator. Verify the closed form,
      energy conservation, Kahane-Peyriere non-degeneracy, log-infinite-divisibility.
  (2) SANITY-CHECK option 2: the lambda tension. The scale ratio lambda that reproduces the
      measured zeta(p) is NOT the one that reproduces the codimension of the most intense
      structures. Test whether that tension is real and Collatz-specific.

Turbulence comparison data is BANKED, zero JHU:
  D:/Turbulence/data/processed/probe11A_analysis.json  (CLM-246, R_lambda = 433/613/1280)

Cascade formalism (b = child count, lambda = LINEAR scale ratio; conflated in 1D where b=lambda):
  density multiplier W, E[W] = 1
  zeta(p) = p/3 - log_lambda E[W^(p/3)]
  check: zeta(3) = 1 - log_lambda E[W] = 1 exactly  (the 4/5 law)
"""
import json
import numpy as np
from pathlib import Path

BANKED = Path("D:/Turbulence/data/processed/probe11A_analysis.json")
OUT = Path("C:/Collatz/experiments_output/cascade_bridge_2026_07_15.json")
OUT.parent.mkdir(parents=True, exist_ok=True)

LN2, LN3 = np.log(2.0), np.log(3.0)


# ---------------------------------------------------------------- generator
def EWq(q, qq=3.0):
    """E[W^q] for W = qq/2^v, v ~ Geom(1/2) on v>=1.  Closed form; converges iff q > -1."""
    return qq ** q / (2.0 ** (q + 1) - 1.0)


def EWq_brute(q, qq=3.0, K=400):
    v = np.arange(1, K)
    return float(np.sum(2.0 ** -v * (qq / 2.0 ** v) ** q))


def log_EWq(q, qq=3.0):
    """ln E[W^q], overflow-safe.  ln(2^(q+1)-1) = (q+1)ln2 + ln1p(-2^-(q+1))."""
    return q * np.log(qq) - ((q + 1.0) * LN2 + np.log1p(-(2.0 ** -(q + 1.0))))


def zeta(p, lam, qq=3.0):
    return p / 3.0 - log_EWq(p / 3.0, qq) / np.log(lam)


res = {"probe": "cascade_bridge_2026_07_15", "generator": {}, "codimension": {}, "fit": {}}

# closed form vs brute force
cf = []
for q in [-0.5, 2 / 3, 1.0, 4 / 3, 2.0, 3.0]:
    a, b = EWq(q), EWq_brute(q)
    cf.append({"q": q, "closed": a, "brute": b, "match": bool(np.isclose(a, b))})
res["generator"]["closed_form_check"] = cf
assert all(c["match"] for c in cf), "closed form disagrees with brute force"

# energy conservation: E[W] = 3 * E[2^-v] = 3 * (1/3) = 1  -- EXACTLY, and circularly
res["generator"]["E_W"] = EWq(1.0)
res["generator"]["E_2mv"] = float(sum(4.0 ** -v for v in range(1, 400)))
res["generator"]["sum_pv_squared"] = float(sum((2.0 ** -v) ** 2 for v in range(1, 400)))

# Kahane-Peyriere non-degeneracy: E[W ln W] < ln(lambda)
EWlnW = LN3 - (4.0 / 3.0) * LN2
res["generator"]["E_W_lnW"] = EWlnW
res["generator"]["KP_pass_lam2"] = bool(EWlnW < LN2)
res["generator"]["KP_pass_lam3"] = bool(EWlnW < LN3)

# mean-critical but a.s.-decaying
res["generator"]["E_lnW"] = LN3 - 2.0 * LN2  # ln(3/4) < 0

# skewness -- the lambda-INDEPENDENT kill (lambda is only a log base; never touches W's law)
m1, m2, m3 = EWq(1.0), EWq(2.0), EWq(3.0)
var = m2 - m1 ** 2
mu3 = m3 - 3 * m1 * m2 + 2 * m1 ** 3
res["generator"]["skew_W"] = float(mu3 / var ** 1.5)
res["generator"]["skew_lnW"] = float(-(2.0 - 0.5) / np.sqrt(0.5))  # -skew(Geom(1/2))

# log-ID: ln W = ln3 - v*ln2, v ~ Geom(1/2) = NegBinom(1,1/2) is ID => compound Poisson
res["generator"]["log_infinitely_divisible"] = True


# ------------------------------------------------------- banked turbulence
bk = json.load(open(BANKED))["by_reynolds"]["1280"]
ess = {int(k): v["val"] for k, v in bk["ess_zeta_over_zeta3"].items()}
res["banked"] = {"R_lambda": 1280, "ess_zeta_over_zeta3": ess, "ess_c1": bk["ess_c1"]}

PS = [2, 3, 4, 5, 6]


def maxdev(lam):
    # Collatz zeta(3) = 1 exactly, so its ESS == its zeta
    return max(abs(zeta(p, lam) - ess[p]) for p in PS)


for lam, name in [(2.0, "lam2_2adic"), (3.0, "lam3_qadic_contraction"), (np.sqrt(2), "lam_sqrt2_codim2")]:
    res["fit"][name] = {
        "lambda": float(lam),
        "zeta": {p: float(zeta(p, lam)) for p in PS},
        "delta_vs_banked_ess": {p: float(zeta(p, lam) - ess[p]) for p in PS},
        "max_abs_dev": float(maxdev(lam)),
    }

# best-fit lambda over the measurable window
grid = np.linspace(1.2, 8.0, 200000)
devs = np.array([maxdev(l) for l in grid])
lam_star = float(grid[np.argmin(devs)])
res["fit"]["best_fit"] = {"lambda_star": lam_star, "max_abs_dev": float(devs.min())}


# -------------------------------------------------- OPTION 2: codimension
# As p -> inf, E[W^q] -> (1/2)(3/2)^q  (the v=1 atom dominates), so
#   zeta(p) -> (p/3)(1 - log_lam(3/2)) + log_lam(2)
# => h_min = (1/3)(1 - log_lam(3/2));  C_inf = log_lam(2) = codim of most intense structures.
def h_min(lam):
    return (1.0 / 3.0) * (1.0 - np.log(1.5) / np.log(lam))


def C_inf(lam):
    return LN2 / np.log(lam)


# numerical confirmation of the asymptote
pbig = np.array([2000.0, 4000.0])
for lam, name in [(2.0, "lam2_2adic"), (3.0, "lam3_qadic_contraction")]:
    z = zeta(pbig, lam)
    slope = float((z[1] - z[0]) / (pbig[1] - pbig[0]))
    res["codimension"][name] = {
        "lambda": float(lam),
        "h_min_analytic": float(h_min(lam)),
        "h_min_numeric": slope,
        "C_inf_analytic": float(C_inf(lam)),
        "C_inf_numeric": float(z[0] - slope * pbig[0]),
    }

# She-Leveque reference: zeta(p) = p/9 + 2(1 - (2/3)^(p/3)) => h_min = 1/9, C_inf = 2 (filaments)
res["codimension"]["she_leveque_ref"] = {"h_min": 1.0 / 9.0, "C_inf": 2.0, "geometry": "filaments (codim 2)"}
# lambda required to force C_inf = 2
res["codimension"]["lambda_forcing_C_inf_2"] = float(2.0 ** 0.5)
res["codimension"]["maxdev_at_that_lambda"] = float(maxdev(np.sqrt(2)))

json.dump(res, open(OUT, "w"), indent=2, default=str)

# ------------------------------------------------------------------ report
g = res["generator"]
print("=" * 74)
print("(1) GENERATOR -- the durable positive")
print("=" * 74)
print(f"  closed form E[W^q] = 3^q/(2^(q+1)-1)   verified vs brute force: {all(c['match'] for c in cf)}")
print(f"  E[W]        = {g['E_W']:.10f}   (energy-conserving -- but circular: restates 3 = 1/E[2^-v])")
print(f"  E[2^-v]     = {g['E_2mv']:.6f}   (annealed multiplier / breakdown coefficient)")
print(f"  sum p_v^2   = {g['sum_pv_squared']:.6f}   (R8 participation ratio -> D2 = log3/log q)")
print(f"    ^ equal ONLY because p_v = 2^-v: the geometric weights ARE the averaged values.")
print(f"  E[W ln W]   = {g['E_W_lnW']:.4f}  <  ln2 = {LN2:.4f}  -> Kahane-Peyriere: {g['KP_pass_lam2']}")
print(f"  E[ln W]     = {g['E_lnW']:.4f}  < 0  -> mean-critical but a.s. decaying")
print(f"  log-ID      = {g['log_infinitely_divisible']} (Geom(1/2) is ID => compound Poisson)")
print(f"  VERDICT: legitimate log-geometric compound-Poisson cascade generator.")
print()
print(f"  skew(W)     = {g['skew_W']:+.4f}")
print(f"  skew(ln W)  = {g['skew_lnW']:+.4f}")
print(f"    ^ NEGATIVE in both conventions. lambda-INDEPENDENT (lambda is only a log base).")
print(f"    ^ JSG 1999: only POSITIVELY skewed laws reproduce measured multipliers. THE KILL.")
print()
print("=" * 74)
print("(2) THE LAMBDA TENSION -- option 2 sanity check")
print("=" * 74)
print("  fit to banked ESS (R_lambda=1280, p=2..6):")
for name in ["lam2_2adic", "lam3_qadic_contraction", "lam_sqrt2_codim2"]:
    f = res["fit"][name]
    print(f"    lambda={f['lambda']:.4f}  ({name:24s})  max|dev| = {f['max_abs_dev']:.4f}")
print(f"    best-fit lambda* = {lam_star:.4f}   max|dev| = {devs.min():.4f}")
print()
print("  codimension of most intense structures  (C_inf = log_lambda 2):")
for name in ["lam2_2adic", "lam3_qadic_contraction"]:
    c = res["codimension"][name]
    print(f"    lambda={c['lambda']:.4f}: h_min={c['h_min_analytic']:.4f} (num {c['h_min_numeric']:.4f})"
          f"   C_inf={c['C_inf_analytic']:.4f} (num {c['C_inf_numeric']:.4f})")
print(f"    She-Leveque:  h_min={1/9:.4f}   C_inf=2.0000  -> filaments")
print()
print(f"  To force C_inf = 2 you need lambda = sqrt(2) = {np.sqrt(2):.4f},")
print(f"  and there max|dev| on zeta(p) explodes to {maxdev(np.sqrt(2)):.4f}.")
print()
print(f"  THE TENSION: lambda={lam_star:.3f} fits zeta(2..6) to {devs.min():.4f} but gives")
print(f"  C_inf = {C_inf(lam_star):.3f} (too space-filling) where turbulence wants 2 (filaments).")
print(f"  No single lambda satisfies both. -> {OUT}")
