"""
probe_lvalue_search_2026_05_30.py

Comprehensive L-value search per user's expanded plan:

TARGETS (test each):
  T1 = c_∞ - c(1) ≈ -1.89e-4  (depth-1 residual)
  T2 = c_∞ - 19/127 ≈ +3.38e-3  (depth-0 residual = Δ_∞)
  T3 = c_∞ ≈ 0.15299  (asymptote itself)

CANDIDATES:
  L(1, χ) for χ in {
    χ_4 mod 17 (two primitive quartic), and their real/imag/sum/diff
    χ_4 mod 5  (two primitive quartic), and same
    χ_2 mod 17 (Legendre)
    χ_2 mod 5  (Legendre, = even nontrivial real char mod 5? actually order 2)
    χ_4 mod 13 (since 13 appears as Gaussian prime norm)
    χ_8 mod 17 (order 8, since ord_17(2)=8 and we suspect octic structure)
  }

NORMALIZATIONS:
  N ∈ {1, q, ϕ(q), √q, Gauss_sum_magnitude}

Then PSLQ each target against the full candidate vector.
"""
from __future__ import annotations
import sys
from mpmath import mp, mpf, mpc, digamma, sqrt, log, pi, exp, pslq, identify
sys.stdout.reconfigure(encoding="utf-8")

mp.dps = 50

def build_chi(q, g_primitive_root, char_indices, order_denom):
    """Build chi(g^k) = exp(2 pi i * char_indices * k / order_denom).
    char_indices is the 'r' index of the character, order_denom is the cycle (e.g. ϕ(q))."""
    table = [mpc(0)] * q
    x = 1
    phi_q = q - 1  # for prime q
    for k in range(phi_q):
        table[x] = exp(2 * pi * mpc(0, 1) * char_indices * k / order_denom)
        x = (x * g_primitive_root) % q
    return table

def L1_chi(chi_table, q):
    """L(1, chi) for primitive non-trivial chi mod q via digamma formula:
    L(1, chi) = -(1/q) sum_{a=1}^{q-1} chi(a) psi(a/q)."""
    total = mpc(0)
    for a in range(1, q):
        total += chi_table[a] * digamma(mpf(a)/mpf(q))
    return -total / mpf(q)

# Primitive roots
g_17 = 3
g_5 = 2
g_13 = 2

# === Build characters ===
chars = {}
# mod 17 (phi = 16): order 4 characters at r=4 and r=12 (conjugates)
chars['chi_4_17_a'] = (build_chi(17, g_17, 4, 16), 17)   # chi(g)=i
chars['chi_4_17_b'] = (build_chi(17, g_17, 12, 16), 17)  # chi(g)=-i (conj)
# Quadratic (Legendre) mod 17
chars['chi_2_17']   = (build_chi(17, g_17, 8, 16), 17)   # chi(g)=-1
# Order 8 mod 17
chars['chi_8_17_a'] = (build_chi(17, g_17, 2, 16), 17)
chars['chi_8_17_b'] = (build_chi(17, g_17, 6, 16), 17)
# Order 16 (faithful) mod 17
chars['chi_16_17_a'] = (build_chi(17, g_17, 1, 16), 17)
chars['chi_16_17_b'] = (build_chi(17, g_17, 5, 16), 17)  # the "Mathar r=5"

# mod 5 (phi = 4): order 4 characters at r=1 and r=3
chars['chi_4_5_a']  = (build_chi(5, g_5, 1, 4), 5)
chars['chi_4_5_b']  = (build_chi(5, g_5, 3, 4), 5)
# Legendre mod 5
chars['chi_2_5']    = (build_chi(5, g_5, 2, 4), 5)

# mod 13 (phi = 12): order 4 at r=3 and r=9
chars['chi_4_13_a'] = (build_chi(13, g_13, 3, 12), 13)
chars['chi_4_13_b'] = (build_chi(13, g_13, 9, 12), 13)
chars['chi_2_13']   = (build_chi(13, g_13, 6, 12), 13)

# === Compute L(1, chi) for each ===
print(f"=== L(1, chi) values (high precision) ===")
L1_vals = {}
for name, (chi, q) in chars.items():
    L = L1_chi(chi, q)
    L1_vals[name] = L
    print(f"  L(1, {name}) = {complex(L.real, L.imag)}  |L|={float(abs(L)):.6e}")

# Verify Legendre mod 17 via closed form
predicted_legendre_17 = log(mpf(4) + sqrt(mpf(17))) / sqrt(mpf(17))
print(f"\n  predicted L(1, χ_2 mod 17) = ln(4+√17)/√17 = {predicted_legendre_17}")
print(f"  (digamma value / 2 = {L1_vals['chi_2_17'].real / 2})  match: {abs(L1_vals['chi_2_17'].real/2 - predicted_legendre_17) < 1e-30}")
# the digamma formula gives 2*(Hurwitz-zeta-style) for even chars; standard is the /2 version

# === Targets ===
# Use exact c(1) and best c_inf estimate
c0 = mpf(19) / mpf(127)
c1_exact_num = 265011804960406635465672455997699
c1_exact_den = 1730087916969634762193659498034425
c1 = mpf(c1_exact_num) / mpf(c1_exact_den)
# c_inf from damped oscillation model (Δ_5 limited precision)
c_inf_model = mpf("0.1529891206058851")
# also use Shanks 3,4,5
c_inf_shanks = mpf("0.1529889944284355")

T1 = c_inf_model - c1
T2 = c_inf_model - c0
T3 = c_inf_model

print(f"\n=== Targets ===")
print(f"  T1 = c_∞ - c(1)     = {T1}")
print(f"  T2 = c_∞ - 19/127   = {T2}")
print(f"  T3 = c_∞            = {T3}")
print(f"\nNote: c_∞ has ~7-digit precision (Δ_5 was float64 FFT). PSLQ tolerance ~1e-6.")

# === Build candidate vector with multiple character L-values and structural forms ===
# For each L-value: real, imag, sum-with-conjugate, diff-with-conjugate
candidates = {}
def add_candidate(name, val):
    candidates[name] = val

# Quartic mod 17 (a, b are conjugates)
La = L1_vals['chi_4_17_a']; Lb = L1_vals['chi_4_17_b']
add_candidate('Re L(1,chi_4_17)', mpf(La.real))
add_candidate('Im L(1,chi_4_17)', mpf(La.imag))
add_candidate('|L(1,chi_4_17)|', abs(La))
add_candidate('La + Lb (real)', mpf((La + Lb).real))
add_candidate('La - Lb (imag·2i)', mpf((La - Lb).imag))

# Quartic mod 5
La5 = L1_vals['chi_4_5_a']; Lb5 = L1_vals['chi_4_5_b']
add_candidate('Re L(1,chi_4_5)', mpf(La5.real))
add_candidate('Im L(1,chi_4_5)', mpf(La5.imag))
add_candidate('|L(1,chi_4_5)|', abs(La5))

# Quartic mod 13
La13 = L1_vals['chi_4_13_a']; Lb13 = L1_vals['chi_4_13_b']
add_candidate('Re L(1,chi_4_13)', mpf(La13.real))
add_candidate('Im L(1,chi_4_13)', mpf(La13.imag))

# Quadratic
add_candidate('L(1, chi_2_17)', mpf(L1_vals['chi_2_17'].real))
add_candidate('L(1, chi_2_5)', mpf(L1_vals['chi_2_5'].real))
add_candidate('L(1, chi_2_13)', mpf(L1_vals['chi_2_13'].real))

# Octic mod 17
L8a = L1_vals['chi_8_17_a']
add_candidate('Re L(1,chi_8_17)', mpf(L8a.real))
add_candidate('Im L(1,chi_8_17)', mpf(L8a.imag))
add_candidate('|L(1,chi_8_17)|', abs(L8a))

# Order-16 mod 17 (Mathar r=5)
L16 = L1_vals['chi_16_17_b']
add_candidate('Re L(1,chi_16_17_r5)', mpf(L16.real))
add_candidate('Im L(1,chi_16_17_r5)', mpf(L16.imag))

# Structural numbers
add_candidate('1/√17', mpf(1)/sqrt(mpf(17)))
add_candidate('1/√5', mpf(1)/sqrt(mpf(5)))
add_candidate('log(4+√17)/√17', log(mpf(4) + sqrt(mpf(17))) / sqrt(mpf(17)))
add_candidate('π/√17', pi/sqrt(mpf(17)))
add_candidate('1', mpf(1))

# === Direct Hecke L over Z[i] for prime π ===
# For π = 1+2i (above 5), N(π) = 5, residue field F_5.
# Quartic residue symbol (α/π)_4 ≡ α^{(N(π)-1)/4} mod π, lifted to {1, i, -1, -i}.
# For π = 1+2i: (α/π)_4 = α mod π, mapped to {1, i, -1, -i} via Z[i]/π ≅ F_5.
# We computed: i ↔ 2 in F_5 (from 2i ≡ -1 mod π, i ≡ -3 ≡ 2 mod 5).
# So unit images: 1→1, i→2, -1→4, -i→3 in F_5.
# Inverse map (F_5* → {1,i,-1,-i}): 1→1, 2→i, 3→-i, 4→-1.
def hecke_quartic_at_1plus2i(alpha_re, alpha_im):
    """Quartic residue symbol (α/π)_4 for π = 1+2i, α coprime to π.
    Computes α mod π in F_5, then lifts back to {1, i, -1, -i}."""
    # α mod π: multiply by 1/(1+2i)? In Z[i], α mod (1+2i) is reduced via norm-form.
    # Use the integer rep: (a + b·i) mod (1+2i) = (a + 2b) mod 5  (since i ≡ 2 mod π).
    r = (alpha_re + 2 * alpha_im) % 5
    if r == 0:
        return None  # not coprime
    # Lift F_5* -> {1, i, -1, -i}
    if r == 1: return mpc(1, 0)
    if r == 2: return mpc(0, 1)
    if r == 3: return mpc(0, -1)
    if r == 4: return mpc(-1, 0)

# L(s, ψ_4) at s=1 via direct sum over Z[i] in a quadrant, with units cancellation handling.
# Since ψ_4 is non-trivial on units (ψ_4(i)=i etc.), the naive sum vanishes.
# Use the "additive" Hecke L: sum over first-quadrant α (one rep per ideal), no unit summation.
# Convergence: slow but mpmath summation handles it.
print(f"\n=== Direct Hecke L(1, ψ_4) on Z[i] mod (1+2i) ===")
N_max = 2000
total = mpc(0)
for a in range(-N_max, N_max + 1):
    for b in range(1, N_max + 1):  # first quadrant: Im > 0
        if a == 0 and b == 0:
            continue
        nrm = a*a + b*b
        if nrm == 0 or nrm > N_max:
            continue
        sym = hecke_quartic_at_1plus2i(a, b)
        if sym is None:
            continue
        total += sym / mpf(nrm)
# Also include positive real axis: b=0, a > 0
for a in range(1, N_max + 1):
    sym = hecke_quartic_at_1plus2i(a, 0)
    if sym is None:
        continue
    total += sym / mpf(a*a)
print(f"  L(1, ψ_4 mod 1+2i) (truncated, N_max={N_max}): {total}")
add_candidate('Re Hecke L(1, ψ_4 / 1+2i)', mpf(total.real))
add_candidate('Im Hecke L(1, ψ_4 / 1+2i)', mpf(total.imag))

print(f"\n=== Candidate vector (size {len(candidates)}) ===")
for name, val in candidates.items():
    print(f"  {name:30} = {float(val):+.12e}")

# === PSLQ tests ===
print(f"\n=== PSLQ search (3 targets × full candidate basis) ===")
cand_names = list(candidates.keys())
cand_vals = list(candidates.values())

# Use a permissive but not absurd tolerance reflecting c_∞ precision
tol = mpf(10) ** (-6)
maxcoeff = 10**8

for target_name, target_val in [('T1 = c_∞ - c(1)', T1), ('T2 = Δ_∞ = c_∞ - 19/127', T2), ('T3 = c_∞', T3)]:
    print(f"\n--- {target_name} = {target_val} ---")
    # PSLQ with target prepended; relation [c, k1, k2, ...] means c·target + sum k_i · cand_i = 0
    vec = [target_val] + cand_vals
    rel = pslq(vec, tol=tol, maxcoeff=maxcoeff)
    if rel is None:
        print(f"  No relation found at tol={tol}, maxcoeff={maxcoeff}")
    else:
        terms = []
        if rel[0] != 0:
            terms.append(f"({rel[0]})·{target_name}")
        for c, name in zip(rel[1:], cand_names):
            if c != 0:
                terms.append(f"({c:+d})·{name}")
        print(f"  Relation: {' + '.join(terms)} ≈ 0")
        # Verify
        residual = sum(c * v for c, v in zip(rel, vec))
        print(f"  Residual: {residual}")

# === Per-candidate single-element PSLQ ===
print(f"\n=== Single-candidate ratio tests (target / candidate) ===")
print(f"   Looking for: target ≈ rational · candidate")
for target_name, target_val in [('T1', T1), ('T2', T2), ('T3', T3)]:
    print(f"\n--- target = {target_name} = {float(target_val):.10e} ---")
    for name, val in candidates.items():
        if abs(val) < mpf(10)**(-30):
            continue
        ratio = target_val / val
        # check if ratio is close to a rational with small denominator
        from fractions import Fraction
        f = Fraction(float(ratio)).limit_denominator(1000)
        approx = mpf(f.numerator) / mpf(f.denominator)
        diff = ratio - approx
        if abs(diff) < mpf(10)**(-5) and f.denominator <= 1000:
            print(f"  {name:30} ratio = {float(ratio):+.10e} ≈ {f}  (diff {float(diff):+.2e})")
