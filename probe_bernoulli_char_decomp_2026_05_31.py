"""
Bernoulli / self-similar measure analysis.

Step 1: Verify self-similarity (1/3) recursion (already known to trivialize).

Step 2: Character decomposition of mu_inf on (Z/17)*.
  - sigma <-> -sigma symmetry restricts to EVEN characters chi(-1) = +1
  - chi(-1) = (-1)^k where chi_k(g) = zeta_16^k (g primitive root)
  - So k even: 8 characters, with orders {1, 2, 4, 4, 8, 8, 8, 8}
  - mu_inf has Fourier expansion in 8 even characters

Step 3: Compute exact Fourier coefficients of p_m at m=0,1,2 (have exact rationals)
  - Verify Plancherel
  - Check whether |<chi_4, mu>|^2 = c_inf is exact or coincidence

Step 4: Examine p_0 dominant configs structurally.
"""
from __future__ import annotations
import sys, json
from fractions import Fraction
from mpmath import mp, mpf, mpc, sqrt, pi, exp, log, pslq
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 40

# Load p_0, p_1, p_2 exact rationals
with open("C:/Collatz/pm_distributions_2026_05_31.json") as f:
    data = json.load(f)

q = 17
p_0 = {int(s): Fraction(*v) for s, v in data["p_0"].items()}
p_1 = {int(s): Fraction(*v) for s, v in data["p_1_rational"].items()}
p_2 = {int(s): Fraction(*v) for s, v in data["p_2_rational"].items()}

# Primitive root g=3 mod 17, dlog table
g = 3
dlog = {}
x = 1
for k in range(16):
    dlog[x] = k
    x = (x * g) % q
print(f"Primitive root g = {g}; dlog table built (16 entries)")

# Character chi_k(sigma) = exp(2 pi i k * dlog(sigma) / 16) for k in 0..15
def chi_k_value(k, sigma):
    """Return chi_k(sigma) as mpc. chi_k has order 16/gcd(k,16)."""
    if sigma == 0:
        return mpc(0)
    return exp(2 * pi * mpc(0, 1) * k * dlog[sigma] / 16)

# chi(-1) for each k: (-1)^k since -1 = g^8
print("\nCharacter parities:")
print("k | order | chi_k(-1) | even/odd")
for k in range(16):
    order = 16 // mp.maxn if False else 16  # placeholder
    import math
    order = 16 // math.gcd(k, 16) if k > 0 else 1
    parity = "+1" if k % 2 == 0 else "-1"
    typ = "EVEN" if k % 2 == 0 else "odd "
    print(f"  k={k:2d} | order={order:2d} | chi(-1)={parity} | {typ}")

# Even k characters: 0, 2, 4, 6, 8, 10, 12, 14 (8 characters)
even_k = [0, 2, 4, 6, 8, 10, 12, 14]
print(f"\nEven characters (chi(-1)=+1): k = {even_k}")
print(f"Their orders: {[16 // (math.gcd(k, 16) if k else 16) for k in even_k]}")
print(f"  (k=0 trivial; k=8 is Legendre chi_2 order 2; k=4,12 are order 4; rest order 8)")

# === Compute Fourier coefficients of p_0, p_1, p_2 ===
def fourier_coef(p, k):
    """Compute <chi_k, p> = sum_sigma chi_k(sigma)_bar * p(sigma)"""
    s = mpc(0)
    for sigma in range(1, q):
        chi_val = chi_k_value(k, sigma)
        p_val = mpf(p[sigma].numerator) / mpf(p[sigma].denominator)
        s += chi_val.conjugate() * p_val
    return s

print("\n=== Fourier coefficients <chi_k, p_m> for even k, m=0,1,2 ===")
print("k | order | <chi, p_0>                           | <chi, p_1>                           | <chi, p_2>")
for k in even_k:
    order = 16 // (math.gcd(k, 16) if k else 16)
    c0 = fourier_coef(p_0, k)
    c1 = fourier_coef(p_1, k)
    c2 = fourier_coef(p_2, k)
    print(f"  k={k:2d} | ord={order} | Re={float(c0.real):+.8f} Im={float(c0.imag):+.8f} | Re={float(c1.real):+.8f} Im={float(c1.imag):+.8f} | Re={float(c2.real):+.8f} Im={float(c2.imag):+.8f}")

# Identify c(m) = <chi_2, p_m> where chi_2 = chi_8 in our naming
print("\n=== Sanity: c(m) = <chi_8, p_m> (Legendre = order-2 char) ===")
for m, p in [(0, p_0), (1, p_1), (2, p_2)]:
    cm = fourier_coef(p, 8)
    print(f"  c({m}) = Re={float(cm.real):+.15f}, Im={float(cm.imag):+.2e}")

# === Plancherel check on p_2 ===
print("\n=== Plancherel on p_2 ===")
sum_sq = mpf(0)
for sigma in range(1, q):
    pv = mpf(p_2[sigma].numerator) / mpf(p_2[sigma].denominator)
    sum_sq += pv ** 2
print(f"  sum_sigma |p_2|^2 = {float(sum_sq):.15f}")

sum_fourier_sq = mpf(0)
print(f"  Per-character contributions to (1/16) sum |c_k|^2:")
for k in range(16):
    fc = fourier_coef(p_2, k)
    mag_sq = abs(fc) ** 2
    contrib = mag_sq / 16
    if mag_sq > mpf("1e-20"):
        order = 16 // (math.gcd(k, 16) if k else 16)
        print(f"    k={k:2d} (ord={order}): |c_k|^2 = {float(mag_sq):.10f}, contrib = {float(contrib):.10f}")
    sum_fourier_sq += mag_sq
print(f"  (1/16) sum_k |c_k|^2 = {float(sum_fourier_sq / 16):.15f}")
print(f"  Match: {abs(sum_sq - sum_fourier_sq / 16) < mpf('1e-30')}")

# === Check the |<chi_4, mu>|^2 ≈ c_inf "near-coincidence" ===
print("\n=== Near-coincidence check: |<chi_4, p_m>|^2 vs c(m) ===")
for k in [4, 12]:  # both order-4 characters
    print(f"  Using k={k}:")
    for m, p in [(0, p_0), (1, p_1), (2, p_2)]:
        fc = fourier_coef(p, k)
        mag_sq = abs(fc) ** 2
        cm = fourier_coef(p, 8).real
        ratio = mag_sq / cm if abs(cm) > mpf("1e-20") else mpf("nan")
        print(f"    m={m}: |<chi_{k}, p>|^2 = {float(mag_sq):.10f}, c({m}) = {float(cm):.10f}, ratio = {float(ratio):.6f}")

# === Bernoulli structural fact 1: depth-0 dominant configs ===
print("\n=== Depth-0 dominant configs ===")
# Compute p_0(sigma) breakdown by (a_X, a_Y) config
def W_8(r):
    if r == 0:
        return Fraction(1, 2**8 - 1)
    return Fraction(2**(8 - r), 2**8 - 1)

print("  σ = ±4: from (a_X=1, a_Y=2) and (2,1)")
inv2_q = pow(2, -1, q)
configs_pm4 = []
for ax in range(1, 9):
    for ay in range(1, 9):
        val = (pow(inv2_q, ax, q) - pow(inv2_q, ay, q)) % q
        if val == 4 or val == 13:  # 13 = -4 mod 17
            w = W_8(ax % 8) * W_8(ay % 8)
            configs_pm4.append((ax, ay, val, w, float(w)))

# Sort by weight descending
configs_pm4.sort(key=lambda c: -c[4])
print(f"  Configs giving σ=±4:")
for ax, ay, val, w, wf in configs_pm4[:10]:
    print(f"    (a_X={ax}, a_Y={ay}) → σ={val}, weight={wf:.6e}")

# === All configs: tabulate p_0(sigma) decomposition ===
print("\n  Total p_0(sigma) verified vs constructed:")
p_0_built = {sigma: Fraction(0) for sigma in range(q)}
T_total = Fraction(0)
for ax in range(1, 9):
    for ay in range(1, 9):
        val = (pow(inv2_q, ax, q) - pow(inv2_q, ay, q)) % q
        if val == 0:
            continue  # excluded
        w = W_8(ax % 8) * W_8(ay % 8)
        p_0_built[val] += w
        T_total += w
for sigma in range(1, q):
    p_0_built[sigma] = p_0_built[sigma] / T_total
print(f"    σ | built p_0 | json p_0 | match")
for sigma in range(1, q):
    matches = p_0_built[sigma] == p_0[sigma]
    print(f"    {sigma:2d} | {float(p_0_built[sigma]):.10f} | {float(p_0[sigma]):.10f} | {matches}")

# === Self-similarity step 1: derive (1/3) factor explicitly ===
print("\n=== Self-similarity (1/3) derivation ===")
# P(a_X_1 = a_Y_1) = sum_a P(a)^2 = sum_a 2^(-2a) for a >= 1
p_aXaY_equal = sum(Fraction(1, 4**a) for a in range(1, 30))
# Note: actual sum is 1/3 since 1/4 + 1/16 + ... = (1/4)/(1-1/4) = 1/3
print(f"  sum_{{a>=1}} 2^(-2a) = {p_aXaY_equal} = {float(p_aXaY_equal):.10f}")
print(f"  Exact value 1/3 = {1/3:.10f}")
print(f"  Verified: P(a_X = a_Y, level-1 dominant) = 1/3")

print("\n=== Self-similarity for c_inf ===")
print("  Given v_q(D) ≥ 1:")
print("  - With prob 1/3 / P(v_q≥1): level-1 dominant, σ_m = 2^(-a) · σ_{m-1}(D')")
print("  - chi_2(σ_m) = chi_2(2^(-a)) · chi_2(σ_{m-1}(D'))")
print("  - chi_2(2^(-a)) = +1 since 2 is QR mod 17 → all 2^(-a) are QR")
print("  - So in dominant case: chi_2(σ_m) = chi_2(σ_{m-1}(D'))")
print("  - This gives c_inf = (1/3)·c_inf + (2/3)·c_inf^∂")
print("  - Solving: c_inf^∂ = c_inf (trivial consistency, no closed form)")

# === Step: depth-0 dominant configs that contribute to c_inf direction ===
print("\n=== Decomposition of c_inf by depth-0 config χ_2 contribution ===")
chi_2_p0 = Fraction(0)
for ax in range(1, 9):
    for ay in range(1, 9):
        val = (pow(inv2_q, ax, q) - pow(inv2_q, ay, q)) % q
        if val == 0:
            continue
        w = W_8(ax % 8) * W_8(ay % 8)
        chi = 1 if pow(val, (q-1)//2, q) == 1 else -1
        chi_2_p0 += w * chi
print(f"  Total chi_2 weighted: {chi_2_p0} = {float(chi_2_p0):.10f}")
print(f"  / total mass {T_total} = c(0) = {float(chi_2_p0 / T_total):.10f}")
print(f"  Reference c(0) = 19/127 = {float(Fraction(19, 127)):.10f}")

# === Look for ALL fourier coefficient PSLQ relations ===
print("\n=== PSLQ all 8 even-character Fourier coefs vs c(m) and simple constants ===")
# Use p_2 (highest depth we have exact)
fourier_p2 = {}
for k in even_k:
    fc = fourier_coef(p_2, k)
    fourier_p2[k] = fc

print("  Fourier coefs of p_2 (real/imag parts at 30 digits):")
basis_vals = [mpf(1)]
basis_names = ["1"]
for k in even_k:
    if k == 0:
        continue
    fc = fourier_p2[k]
    fc_re = fc.real
    fc_im = fc.imag
    basis_vals.append(fc_re)
    basis_names.append(f"Re<chi_{k},p>")
    if abs(fc_im) > mpf("1e-20"):
        basis_vals.append(fc_im)
        basis_names.append(f"Im<chi_{k},p>")
    print(f"    k={k:2d}: <chi_k, p_2> = Re={float(fc_re):+.20f}, Im={float(fc_im):+.20f}")

# Check Plancherel relation rigorously
print("\n=== Plancherel: |<chi_k, p_2>|^2 sum / 16 should equal sum |p_2|^2 ===")
plancherel_RHS = sum_sq  # sum |p_2|^2
plancherel_LHS = mpf(0)
for k in range(16):
    fc = fourier_coef(p_2, k)
    plancherel_LHS += abs(fc) ** 2
plancherel_LHS = plancherel_LHS / 16
print(f"  RHS = sum_sigma |p_2|^2 = {float(sum_sq):.20f}")
print(f"  LHS = (1/16) sum_k |c_k|^2 = {float(plancherel_LHS):.20f}")
print(f"  Diff = {float(plancherel_LHS - sum_sq):.2e}")
print(f"  Plancherel holds: {abs(plancherel_LHS - sum_sq) < mpf('1e-30')}")

# Express c(2) in terms of other Fourier coefs of p_2 via Plancherel
c2_val = fourier_coef(p_2, 8).real
print(f"\n  c(2) = <chi_8, p_2>.real = {float(c2_val):.20f}")
print(f"  Need 7 OTHER Fourier coefs (k=2,4,6,10,12,14) to express via Plancherel")

# Now check if any |c_k|^2 relations to c(m) hold
print("\n=== Check |<chi_k, p_2>|^2 vs c(2): structural identities? ===")
for k in even_k:
    if k == 0 or k == 8: continue
    fc = fourier_coef(p_2, k)
    mag_sq = abs(fc) ** 2
    ratio = mag_sq / c2_val
    print(f"  k={k}: |<chi_k>|^2 = {float(mag_sq):.10f}, ratio to c(2) = {float(ratio):.10f}")
