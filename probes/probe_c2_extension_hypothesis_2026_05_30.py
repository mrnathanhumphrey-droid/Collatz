"""
probe_c2_extension_hypothesis_2026_05_30.py

Test the natural extension of the c(1) machinery to c(2):
  c(m) = sum_{k_1,...,k_m} (prod_i W(k_i)) N(2*(sum k_i) mod q)
        / sum_{...} (prod_i W(k_i)) T(2*(sum k_i) mod q)

where W(k) = 2^{-8|k|}, W(0) = 1, and N(s), T(s) are the depth-0 character moments.

If c(2)_hypothesis matches the FFT-measured c(2) ≈ 0.153247920779, the extension is valid
and we can compute c(2), c(3), c(4) exactly.

Equivalently: c(m) = sum_s W_m(s) N(s) / sum_s W_m(s) T(s)
where W_m = (W_1)^{*m} convolved on Z/q under the s = 2k mod q parameterization.
"""
from fractions import Fraction
import sys
sys.stdout.reconfigure(encoding="utf-8")

q = 17
ord_2 = 8

def chi(x, q):
    x %= q
    if x == 0: return 0
    return 1 if pow(x, (q - 1) // 2, q) == 1 else -1

inv2_q = pow(2, q - 2, q)
pow_inv2 = [pow(inv2_q, a, q) for a in range(ord_2)]   # 2^-a mod q

def W_geom(a_res):
    """Sum of 2^-a over a >= 1 with a ≡ a_res mod ord_2."""
    if a_res == 0:
        return Fraction(1, 2**ord_2 - 1)
    return Fraction(2**(ord_2 - a_res), 2**ord_2 - 1)

# N(s), T(s) — depth-0 moments with shift s
N = {}; T = {}
for s in range(q):
    Ns = Fraction(0); Ts = Fraction(0)
    for ar in range(ord_2):
        for br in range(ord_2):
            v = (s + pow_inv2[ar] - pow_inv2[br]) % q
            w = W_geom(ar) * W_geom(br)
            if v != 0:
                Ts += w
                Ns += w * chi(v, q)
    N[s] = Ns; T[s] = Ts

c0 = N[0] / T[0]
print(f"c(0) = {c0} = {float(c0):.15f}  (match 19/127? {c0 == Fraction(19,127)})")

# W_1(s) = sum over k with 2k = s mod q of 2^{-8|k|}
# k ranges over integers; per residue s, k = k_0(s) + 17j where k_0(s) = 9s mod 17 mapped to {-8,...,8}.
def k0_of_s(s):
    """Smallest |k| with 2k ≡ s mod q. Returns int in {-8,...,8}."""
    k = (9 * s) % q  # 9 = 2^-1 mod 17
    if k > q // 2: k -= q
    return k

# For each shift s in Z/q, weights from k ∈ {k_0(s) + 17j : j ∈ Z}
# Truncate |j| <= JMAX. Weight 2^{-8|k|}.
JMAX = 30  # so smallest weight ~ 2^{-8*8*30} = 2^{-1920}, way below any precision needed

def W1(s):
    k0 = k0_of_s(s)
    w = Fraction(0)
    for j in range(-JMAX, JMAX + 1):
        k = k0 + 17 * j
        w += Fraction(1, 2 ** (8 * abs(k))) if k != 0 else Fraction(1)
    return w

W1_vec = {s: W1(s) for s in range(q)}

# c(1)_hypothesis = sum_s W1(s) N(s) / sum_s W1(s) T(s)
num1 = sum(W1_vec[s] * N[s] for s in range(q))
den1 = sum(W1_vec[s] * T[s] for s in range(q))
c1 = num1 / den1
print(f"\nc(1) (hypothesis, JMAX={JMAX}) = {float(c1):.15f}")
print(f"  FFT-measured c(1)            = 0.153178230055")
print(f"  diff = {float(c1) - 0.153178230055:.2e}")

# Now c(2): W_2(s) = sum over (k_1, k_2) with 2(k_1+k_2) ≡ s mod q of 2^{-8(|k_1|+|k_2|)}
# Equivalently: W_2 = W_1 * W_1 (convolution on Z/q via integer-k accumulation).
# But careful: the convolution is over INTEGERS k_1, k_2, and the shift accumulation is 2(k_1+k_2) mod q.
# So W_2(s) = sum_{s_1} W_1(s_1) * W_1((s - s_1) mod q) where the addition is mod q.
# This is the convolution on Z/q of W_1 with itself.
def conv_Zq(f, g):
    h = {s: Fraction(0) for s in range(q)}
    for s1 in range(q):
        for s2 in range(q):
            h[(s1 + s2) % q] += f[s1] * g[s2]
    return h

W2_vec = conv_Zq(W1_vec, W1_vec)
num2 = sum(W2_vec[s] * N[s] for s in range(q))
den2 = sum(W2_vec[s] * T[s] for s in range(q))
c2_hyp = num2 / den2
print(f"\nc(2) (hypothesis, JMAX={JMAX}) = {float(c2_hyp):.15f}")
print(f"  FFT-measured c(2)             = 0.153247920779")
print(f"  diff = {float(c2_hyp) - 0.153247920779:.2e}")

# c(3)
W3_vec = conv_Zq(W2_vec, W1_vec)
num3 = sum(W3_vec[s] * N[s] for s in range(q))
den3 = sum(W3_vec[s] * T[s] for s in range(q))
c3_hyp = num3 / den3
print(f"\nc(3) (hypothesis) = {float(c3_hyp):.15f}")
print(f"  FFT-measured c(3) = 0.153005331692")
print(f"  diff = {float(c3_hyp) - 0.153005331692:.2e}")

# c(4)
W4_vec = conv_Zq(W3_vec, W1_vec)
num4 = sum(W4_vec[s] * N[s] for s in range(q))
den4 = sum(W4_vec[s] * T[s] for s in range(q))
c4_hyp = num4 / den4
print(f"\nc(4) (hypothesis) = {float(c4_hyp):.15f}")
print(f"  FFT-measured c(4) = 0.152988709047")
print(f"  diff = {float(c4_hyp) - 0.152988709047:.2e}")

print(f"\nIf the c(2), c(3), c(4) hypothesis-vs-FFT diffs are all ~1e-12 or smaller,")
print(f"the extension is structurally correct and we have exact c(m) for any m.")
print(f"If they DIVERGE, the shift-accumulation hypothesis is wrong.")
