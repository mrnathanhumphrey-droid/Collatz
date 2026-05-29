"""
probe_coset_mechanism_2026_05_29.py

Test the closed-form mechanism for the coset-mass asymmetry (q=17, index 2).

The coset masses are governed by A_chi = sum_xi chi(xi)|mu_hat(xi)|^2 for chi a Dirichlet
character trivial on <2> (here chi = Legendre symbol mod 17, the unique nontrivial one).
E_<2> = (1 + A_chi/A_triv)/2,  E_dark = (1 - A_chi/A_triv)/2  => asymmetry = A_chi/A_triv.

Self-similar prediction: X = 2^-a (1+qX'), chi(2)=1, so the chi-value of the unit part of the
difference X-Y reduces (level by level) to chi(1 - 2^Delta), Delta = b-a (iid Geom(2) gaps).
=> a candidate closed-form V = sum_{Delta != 0} w(Delta) chi(1-2^Delta),  w(Delta)=2^-|Delta|/3,
with chi(0)=0 on the deepening residues (ord_2(q) | Delta). chi(1-2^Delta) is periodic in
Delta mod ord_2(q) => V is an explicit finite combination of geometric series (closed form).

Compare A_chi/A_triv (from FFT) to V (and to V times the Gauss factor) to find the exact relation.
"""
from __future__ import annotations
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
from probe_coset_closedform_2026_05_29 import offset_distribution, subgroup_pow2

def legendre(x, q):
    x %= q
    if x == 0: return 0
    return 1 if pow(x, (q - 1) // 2, q) == 1 else -1

def A_ratio(q, n, chi_q):
    """A_chi/A_triv from FFT, chi a Dirichlet char mod q given as dict on units."""
    N = q ** n
    P = offset_distribution(q, n, 100)
    mu = np.fft.fft(P); a2 = np.abs(mu) ** 2
    xi = np.arange(N)
    units = (xi % q != 0)
    chi_vals = np.zeros(N)
    chi_vals[units] = np.array([chi_q[x % q] for x in xi[units]])
    A_chi = np.sum(chi_vals * a2)
    A_triv = np.sum(a2[units])
    return A_chi / A_triv, A_triv

def main():
    q = 17
    ord2 = 1; t = 2 % q
    while t != 1:
        t = (t * 2) % q; ord2 += 1
    print(f"q={q}, ord_2(q)={ord2}, <2> index = {(q-1)//ord2}")
    chi_q = {x: legendre(x, q) for x in range(1, q)}
    print(f"  chi = Legendre mod {q}; <2> = {{x: chi=+1}} = {sorted(x for x in range(1,q) if chi_q[x]==1)}")

    print("\n A_chi/A_triv (= coset asymmetry E_<2> - E_dark) vs n:")
    for n in (2, 3, 4, 5):
        r, _ = A_ratio(q, n, chi_q)
        print(f"   n={n}: asymmetry = {r:.10f}   => E_<2>={(1+r)/2:.10f}")

    # candidate closed-form sum V = sum_{Delta!=0} 2^-|Delta|/3 * chi(1-2^Delta mod q)
    # chi(1-2^Delta) periodic in Delta mod ord2; sum each residue class as a geometric series.
    # For Delta>0: 2^Delta mod q.  For Delta<0: 2^Delta = inv2^|Delta| mod q.
    inv2 = pow(2, -1, q)
    def chi_1m2(Delta):
        if Delta >= 0:
            val = (1 - pow(2, Delta, q)) % q
        else:
            val = (1 - pow(inv2, -Delta, q)) % q
        return chi_q.get(val % q, 0) if val % q != 0 else 0
    # sum over Delta in [-K,K]\{0}
    K = 4000
    V = sum((2.0 ** (-abs(D)) / 3.0) * chi_1m2(D) for D in range(-K, K + 1) if D != 0)
    print(f"\n candidate V = sum_{{Delta!=0}} 2^-|Delta|/3 * chi(1-2^Delta) = {V:.10f}")
    # Gauss sum g(chi) = sum_x chi(x) e(x/q)
    g = sum(chi_q[x] * np.exp(2j * np.pi * x / q) for x in range(1, q))
    print(f" Gauss sum g(chi) = {g.real:.6f} + {g.imag:.6f}i   (|g|=sqrt(q)={np.sqrt(q):.6f})")
    r5, _ = A_ratio(q, 5, chi_q)
    print(f"\n compare: asymmetry(n=5) = {r5:.10f}")
    print(f"          V                = {V:.10f}   ratio asym/V = {r5/V:.6f}")
    print(f"          V * sqrt(q)       = {V*np.sqrt(q):.10f}")
    print(f"          V * q             = {V*q:.10f}")

if __name__ == "__main__":
    main()
