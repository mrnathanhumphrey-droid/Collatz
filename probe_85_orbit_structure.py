"""
Probe 85 — orbit-structure diagnostic (decides which evidential test is valid).

Before wiring an F-hat substitution, answer the load-bearing structural question
the placeholder skipped: does the ACTUAL DWM j=2 chirp argument

    Delta = 2^{-b} * (2^{-v} - 2^{-v'})   (mod 3^{r+1}, r=n-3)

live on the multiplicative <4> orbit (where R81's F-hat chirp e_q(c*4^u) lives),
or does it walk the full group / the odd coset? A pointwise "DWM chirp vs F-hat
chirp" comparison is only VALID if Delta in <4>. If Delta hits both cosets, the
only valid bridge test is the moment-assembly (Fourier) one, not a pointwise swap.

Method: mod q=3^{r+1} the unit group is cyclic of order 2*3^r, generator g=2
(2 is a primitive root mod 3^k). x in <4> iff dlog_2(x) is EVEN. Build the dlog
table, enumerate every Delta that actually occurs in the n=5,6 moment sums, and
report the dlog-parity distribution. Fast, no moments, decisive.
"""
import sys, math
sys.stdout.reconfigure(encoding="utf-8")

def dlog_table(q):
    """discrete log base 2 for every unit mod q; 2 is a primitive root mod 3^k."""
    order = 2 * (q // 3)  # phi(3^{k}) = 2*3^{k-1}; here q=3^{k}, order=2*3^{k-1}=2*(q/3)
    t = {}
    x = 1
    for e in range(order):
        t[x] = e
        x = (x * 2) % q
    return t, order

def analyze(n):
    r = n - 3
    q = 3 ** (r + 1)          # reduced modulus for j=2 (9 * z mod 3^n = z mod 3^{n-2})
    dl, order = dlog_table(q)
    inv2 = pow(2, -1, q)
    V_MAX = 16
    # b ranges over b1 = v1+vp1 (v1!=vp1 in [1,V_MAX]); (v2,vp2) the j=2 realization
    b_vals = sorted({v1 + vp1 for v1 in range(1, V_MAX + 1)
                             for vp1 in range(1, V_MAX + 1) if v1 != vp1})
    deltas = set()
    for b in b_vals:
        ib = pow(inv2, b, q)
        for v in range(1, V_MAX + 1):
            for vp in range(1, V_MAX + 1):
                if v == vp:
                    continue
                d = (ib * ((pow(inv2, v, q) - pow(inv2, vp, q)) % q)) % q
                deltas.add(d)
    # classify by dlog parity (units only; Delta could be 0 or a non-unit? check)
    even = odd = nonunit = zero = 0
    for d in deltas:
        if d == 0:
            zero += 1
        elif d in dl:
            if dl[d] % 2 == 0:
                even += 1
            else:
                odd += 1
        else:
            nonunit += 1
    tot = len(deltas)
    print(f"\n n={n}  (r={r}, modulus q=3^{r+1}={q}, group order {order}, |<4>|={order//2})")
    print(f"   distinct Delta values occurring: {tot}")
    print(f"     in <4>      (dlog even) : {even:4d}  ({100*even/tot:5.1f}%)")
    print(f"     in 2*<4>    (dlog odd)  : {odd:4d}  ({100*odd/tot:5.1f}%)")
    print(f"     non-unit (v3>0)         : {nonunit:4d}  ({100*nonunit/tot:5.1f}%)")
    print(f"     zero                    : {zero:4d}")
    # verdict
    if odd == 0 and nonunit == 0:
        print(f"   -> Delta subset of <4>: pointwise F-hat comparison is STRUCTURALLY VALID.")
    else:
        print(f"   -> Delta spans BOTH cosets / non-units: pointwise F-hat swap is INVALID;")
        print(f"      only the moment-assembly (Fourier) test is legitimate.")
    # also: how many distinct non-unit valuations? (the 3-adic layer structure)
    if nonunit:
        vals = {}
        for d in deltas:
            if d != 0 and d not in dl:
                k = 0; x = d
                while x % 3 == 0:
                    x //= 3; k += 1
                vals[k] = vals.get(k, 0) + 1
        print(f"      non-unit v3 histogram: {dict(sorted(vals.items()))}")

if __name__ == "__main__":
    print("# PROBE 85 — DWM j=2 chirp-argument orbit structure")
    print("# question: does Delta = 2^{-b}(2^{-v}-2^{-v'}) live on <4>? (decides valid test)")
    for n in (5, 6):
        analyze(n)
    print("\n(<4> = squares mod 3^{r+1} = where R81's F-hat chirp e_q(c*4^u) is supported.)")
