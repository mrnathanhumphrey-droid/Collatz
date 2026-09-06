"""
PROBE PEEL-ONE-BRICK (Wilson spec) -- is S_inf's transcendence reducible to ONE generator log2(3)?
Tests the DIMENSION of the unipotent kernel via the frozen-layer contributions Lambda_r, NOT the value.
Kill-first. De-rotate the period-9.06 (omega=ln2) modulation BEFORE any algebraicity test (the trap:
log2(3) is baked into the rotation, so a raw ratio test finds it circularly). High-precision PSLQ with
mandatory negative controls; report the FLOOR. Pre-registered A (one brick) / B (inf-dim) / C (null).

DO NOT evaluate S_inf. Inputs are banked (copy, not recompute).
"""
import os, json
from fractions import Fraction as F
import numpy as np
import mpmath as mp
mp.mp.dps = 60

HERE = os.path.dirname(os.path.abspath(__file__))
R = lambda *p: os.path.join(HERE, '..', *p)

# ---------- inputs: EXACT Lambda_1..7 (arbitrary precision) ----------
EPS = {int(k): F(int(v['num']), int(v['den']))
       for k, v in json.load(open(R('experiments_output',
             'result_77_7_eps_exact_through_k8_v2_vec_pool.json'))).items()}
S = {k: F(7, 15) + EPS[k] for k in EPS}                 # exact S_k, k=1..8
T = {i: S[i + 1] / 2 for i in range(0, 8) if (i + 1) in S}   # T_i = S_{i+1}/2, i=0..7
Lam_exact = {r: T[r] - T[r - 1] for r in range(1, 8) if r in T and (r - 1) in T}  # r=1..7 EXACT
Rexact = max(Lam_exact)

# ---------- inputs: FLOAT tower T_15..20 -> Lambda_16..20 (~12 digits) ----------
Tf = {}
for fn in ['lambda_out/solstice_T.json', 'results/solstice_T17_T18.json', 'lambda_out/solstice_T20.json']:
    p = R(*fn.split('/'))
    if os.path.exists(p):
        try:
            d = json.load(open(p))
            for k, v in (d.items() if isinstance(d, dict) else []):
                if str(k).lstrip('-').isdigit() and isinstance(v, (int, float)):
                    Tf[int(k)] = float(v)
        except Exception as e:
            print('  (skip', fn, e, ')')
Lam_float = {r: Tf[r] - Tf[r - 1] for r in sorted(Tf) if (r - 1) in Tf}   # ~ r=16..20

print("# PROBE PEEL-ONE-BRICK\n")
print("## Inputs")
print(f"   EXACT layers Lambda_r, r=1..{Rexact} (arbitrary precision):")
for r in sorted(Lam_exact):
    print(f"     Lambda_{r} = {Lam_exact[r]}  (~{mp.nstr(mp.mpf(Lam_exact[r].numerator)/Lam_exact[r].denominator, 8)})")
print(f"   FLOAT tower T_i present for i = {sorted(Tf)}  -> Lambda_r float for r = {sorted(Lam_float)}")
for r in sorted(Lam_float):
    print(f"     Lambda_{r} ~ {Lam_float[r]:.6e}  (from T_{r}-T_{r-1}, ~12 sig digits after differencing)")

# ---------- generators (known to arbitrary precision; the LIMIT is the problem, not these) ----------
log2_3 = mp.log(3) / mp.log(2)
log2_5 = mp.log(5) / mp.log(2)
omega = mp.log(2)                 # period 2*pi/ln2 ~ 9.06 in r
period = float(2 * mp.pi / omega)
print(f"\n   omega_lp = ln2 = {mp.nstr(omega,8)}  => modulation PERIOD in r = {period:.3f} levels")
print(f"   log2(3) = {mp.nstr(log2_3,20)}   log2(5) [neg control] = {mp.nstr(log2_5,20)}")

# ---------- STEP 0: de-rotation feasibility on the EXACT contiguous block ----------
print("\n## Step 0 -- de-rotation (MANDATORY, must precede any algebraicity test)")
rs = sorted(Lam_exact); y = [mp.log(abs(mp.mpf(Lam_exact[r].numerator)/Lam_exact[r].denominator)) for r in rs]
span_periods = (rs[-1] - rs[0]) / period
print(f"   exact block spans r={rs[0]}..{rs[-1]} = {span_periods:.2f} PERIODS of the {period:.2f}-level modulation")
# model log|Lambda_r| = c0 + c1*r + A*cos(omega r) + B*sin(omega r)  (4 params)
Amat = np.array([[1.0, float(r), float(mp.cos(omega*r)), float(mp.sin(omega*r))] for r in rs])
yv = np.array([float(v) for v in y])
coef, res, rank, sv = np.linalg.lstsq(Amat, yv, rcond=None)
dof = len(rs) - Amat.shape[1]
cond = sv[0]/sv[-1]
resid = yv - Amat @ coef
print(f"   fit [1, r, cos(omega r), sin(omega r)] on {len(rs)} points: dof = {dof}, design cond# = {cond:.1f}")
print(f"   periodic amplitude sqrt(A^2+B^2) = {np.hypot(coef[2],coef[3]):.4f}  (geom rate exp(c1) = {np.exp(coef[1]):.4f})")
print(f"   fit residual (de-rotated log|Lambda|): max |r~| = {np.abs(resid).max():.2e}")
print(f"   => {dof} degrees of freedom over {span_periods:.2f} periods: the periodic (rotation) component is")
print(f"      NOT SEPARABLE from the trend at <1 period. The de-rotation is UNDERDETERMINED on the exact block;")
print(f"      the residual is fit-noise, not a signal that can be tested for a limit. [de-rotation floor]")

# ---------- STEP 1: the residual limit + PSLQ, at ACHIEVABLE precision ----------
print("\n## Step 1 -- residual limit L and PSLQ (guardrail: >=50 digits stable 50->100)")
if len(Lam_float) >= 3:
    fr = sorted(Lam_float)
    rho = [Lam_float[fr[i+1]]/Lam_float[fr[i]] for i in range(len(fr)-1)]   # float ratios ~0.87
    print(f"   float ratios rho_r (r={fr[0]}..): {[f'{x:.6f}' for x in rho]}")
    # Aitken/Richardson on the ratios -> a limit proxy L (the best 'de-rotated' limit the data allows)
    a = np.array(rho)
    if len(a) >= 3:
        aitken = a[:-2] - (a[1:-1]-a[:-2])**2/((a[2:]-2*a[1:-1]+a[:-2]) + 1e-300)
        L = float(aitken[-1])
    else:
        L = float(a[-1])
    # precision of L: dominated by float differencing (~1e-12 abs on Lambda~1e-4 => ~1e-8 on ratios,
    # amplified by extrapolation). Be generous: ~8 good digits.
    Ldigits = 8
    print(f"   Aitken limit proxy L ~ {L:.10f}  (ACHIEVABLE precision ~ {Ldigits} digits, float-capped)")
    print(f"   NOTE: L is NOT de-rotated (Step 0 floor) and is float-limited; this PSLQ is a FLOOR DEMO, not a test.")
    def pslq_try(Lval, gen, name, digits):
        # working precision fixed >=53 bits; the DATA floor is encoded in tol=10^-digits
        with mp.workdps(60):
            g = gen; Lm = mp.mpf(Lval)
            basis = [mp.mpf(1), Lm, g, Lm*g, g**2, Lm**2]
            rel = mp.pslq(basis, maxcoeff=10**6, maxsteps=10**4, tol=mp.mpf(10)**(-digits))
        return rel
    for digits in (Ldigits, 20, 50):
        r3 = pslq_try(L, log2_3, 'log2_3', digits)
        r5 = pslq_try(L, log2_5, 'log2_5', digits)
        tag = "" if digits <= Ldigits else "  <-- BEYOND DATA PRECISION (spurious either way)"
        print(f"   @ {digits:>3} digits: log2_3 relation = {r3} ; log2_5 [control] = {r5}{tag}")
    print("   => negative control log2(5) is treated identically to log2(3); at the achievable ~8 digits any")
    print("      'relation' cannot be verified STABLE across 50->100 (the data has < 50 digits). Guardrail UNMET.")
else:
    print("   too few float layers to even attempt an extrapolated limit.")

# ---------- STEP 2: lobe-ledger cross-check floor ----------
print("\n## Step 2 -- lobe-ledger (L_r, M_r) cross-check")
print("   The lobe restratification Lambda_r = L_r - M_r (result_P6J) is EXACT only where Lambda_r is exact")
print("   (r<=7); so the lobe cross-check inherits the SAME 7-exact-layer / <1-period floor. No independent lift.")

# ---------- FLOOR + VERDICT ----------
print("\n## FLOOR and VERDICT")
need_periods = 2.5; need_digits = 50
print(f"   The A-vs-B distinction requires: (i) de-rotation over >= ~{need_periods} periods (>= ~{need_periods*period:.0f} contiguous")
print(f"   layers) AND (ii) a de-rotated residual limit to >= {need_digits} digits, stable 50->100.")
print(f"   Achievable: (i) {span_periods:.2f} periods exact ({Rexact} layers; more is blocked -- exact S_k needs the")
print(f"   doubly-exponential stationary chain, infeasible past k~8) ; (ii) ~8 digits (float tower differencing).")
print(f"   Deficit: ~{need_periods-span_periods:.1f} periods of length AND ~{need_digits-8} digits of precision.")
print(f"   ==> OUTCOME C (NULL): the reducible-vs-independent (1-dim vs inf-dim kernel) distinction is BELOW FLOOR")
print(f"       at achievable k and precision. NOT Outcome A, NOT Outcome B -- undecidable with banked data.")
print(f"   The capstone stands UNCHANGED: this probe cannot sharpen it (A) nor is it needed to (B). Value untouched.")
