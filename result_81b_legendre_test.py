"""
R81b-Legendre — does v3(c_k) of the arg F-hat Mahler profile close to a Kummer/Legendre form?

Candidate (user): v3(c_k) ~ k*v3(log_3 4) - v3(k!),  v3(k!) = (k - s3(k))/2 (Legendre, p=3).
Test that + variants + integer regression. Extend r high enough to STABILIZE c_k.
"""
import numpy as np

def s3(k):
    s = 0
    while k:
        s += k % 3; k //= 3
    return s

def v3_fact(k):            # Legendre: v3(k!) = (k - s3(k))/2
    return (k - s3(k)) // 2

def v3(n, cap):
    if n % (3**cap) == 0:
        return cap         # not yet stabilized (>= cap)
    k = 0
    while n % 3 == 0:
        n //= 3; k += 1
    return k

def mahler_v3(r, ell=1, eps=0):
    q = 3**(r+1); d = 3**r; a0 = 1 if eps == 0 else 2
    c = (pow(2, eps, q) * pow((1 + 3**r) % q, ell, q)) % q
    pow4 = np.empty(d, dtype=np.int64); acc = 1
    for j in range(d):
        pow4[j] = acc; acc = (acc*4) % q
    chirp = np.exp(2j*np.pi*((c*pow4) % q)/q)
    ghat = np.fft.fft(chirp)
    KHEAD = 40                                   # only need first ~KHEAD support points
    a = np.arange(a0, a0 + 3*KHEAD, 3)
    z = ghat[a]**2 / q
    J2 = np.rint((np.angle(z) % (2*np.pi))*(2*q)/(2*np.pi)).astype(np.int64) % (2*q)
    s = [int(x) % q for x in J2]
    # finite differences on the HEAD only, mod q  (Delta^k s(0) needs s(0..k))
    diffs = [list(s)]
    while len(diffs[-1]) > 1:
        prev = diffs[-1]
        diffs.append([(prev[i+1] - prev[i]) % q for i in range(len(prev)-1)])
    ck = [diffs[k][0] % q for k in range(len(diffs))]
    return q, ck

if __name__ == "__main__":
    R = 15                      # 3^15 FFT ~ stabilizes c_k up to v3 < 15
    print(f"Computing Mahler coeffs at r={R} (q=3^{R+1}) ...")
    q, ck = mahler_v3(R)
    cap = R + 1
    K = 13
    print(f"\n{'k':>2} {'v3(c_k)':>7} {'s3(k)':>5} {'v3(k!)':>6} "
          f"{'k-v3(k!)':>8} {'2k-v3(k!)':>9} {'(3k+s3)/2':>9} {'v3-k':>5}")
    obs = []
    for k in range(K):
        vk = v3(ck[k], cap)
        obs.append(vk)
        f_user = k*1 - v3_fact(k)          # v3(log_3 4)=1
        f_2    = 2*k - v3_fact(k)
        f_alt  = (3*k + s3(k))//2
        stab = "" if vk < cap else "  (NOT STABLE, >=cap)"
        print(f"{k:>2} {vk:>7} {s3(k):>5} {v3_fact(k):>6} {f_user:>8} {f_2:>9} {f_alt:>9} {vk-k:>5}{stab}")

    # integer regression: v3(c_k) = A*k + B*s3(k) + C*v3(k!) + D  over stable k
    stable = [k for k in range(K) if obs[k] < cap and k <= (K-1)]
    M = np.array([[k, s3(k), v3_fact(k), 1] for k in stable], dtype=float)
    y = np.array([obs[k] for k in stable], dtype=float)
    coef, res, *_ = np.linalg.lstsq(M, y, rcond=None)
    pred = M @ coef
    print(f"\nLeast-squares  v3 = A*k + B*s3(k) + C*v3(k!) + D  over stable k={stable}")
    print(f"  coef [A,B,C,D] = {np.round(coef,3)}   max|resid| = {np.max(np.abs(pred-y)):.3f}")
    print(f"  (exact integer closed form requires max|resid|=0 with integer/half-integer coef)")

    # explicit check of the user candidate
    du = [obs[k] - (k*1 - v3_fact(k)) for k in stable]
    print(f"\nUser formula k*v3(log3 4=1) - v3(k!) residuals (obs-pred): {du}")
    print(f"  -> {'MATCHES' if all(x==0 for x in du) else 'DOES NOT MATCH'}")
