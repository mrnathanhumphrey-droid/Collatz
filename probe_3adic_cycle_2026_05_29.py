"""
3-adic complement of the 2-adic prefix structure, applied to the cycle equation.
(2^K - 3^L) m0 = sum_{i=0}^{L-1} 3^{L-1-i} 2^{S_i},  S_i = k_0+...+k_{i-1}, S_0=0.
Mod 3^ell the RHS is a base-3-style expansion; term i sits at 3-adic level (L-1-i).
We compute the ACTUAL 3-adic structure: digit ladder, the m0 mod 3^ell congruence in terms
of {k_i}, and the distribution of S_i mod ord_2(3^ell). No assumptions about what comes out.
"""
import numpy as np
from collections import Counter

def syr_step(m):
    x=3*m+1; k=(x&-x).bit_length()-1; return x>>k,k

def traj(m0,L):
    m=m0; ks=[]
    for _ in range(L): m,k=syr_step(m); ks.append(k)
    return ks, m  # halving counts, final odd

def cycle_eq(m0,ks,L):
    K=sum(ks); S=[sum(ks[:i]) for i in range(L)]
    RHS=sum(3**(L-1-i)*2**S[i] for i in range(L))
    return K,S,RHS

print("="*72)
print("(1) 3-adic digit ladder of RHS for the trivial cycle and trajectories")
print("="*72)
# trivial cycle
ks=[2]; L=1; K,S,RHS=cycle_eq(1,ks,L)
print(f"trivial {{1,4,2,1}}: K={K} S={S} RHS={RHS}  (2^K-3^L)={2**K-3**L}  m0={RHS//(2**K-3**L)}")
print(f"  RHS mod 3^ell: " + ", ".join(f"3^{e}:{RHS%3**e}" for e in range(1,4)))

import random; random.seed(0)
print("\nReal Syracuse trajectory segments (not cycles, but the identity holds): S_i mod ord, 2^S_i mod 3:")
for m0 in [27, 703, 871]:
    L=12; ks,mL=traj(m0,L); K,S,RHS=cycle_eq(m0,ks,L)
    # 2^{S_i} mod 3 = (-1)^{S_i}; the 3-adic digit pattern
    twoS_mod3=[pow(2,s,3) for s in S]
    print(f" m0={m0:>4} ks={ks}")
    print(f"   S_i        ={S}")
    print(f"   2^S_i mod3 ={twoS_mod3}  (=(-1)^S_i; parity of S_i drives mod-3 digit)")

print("\n"+"="*72)
print("(2) derived congruence: m0 mod 3 in terms of last halving k_{L-1}")
print("="*72)
print("Claim from algebra: m0 = (-1)^{k_{L-1}} mod 3.  Check on the actual identity")
print("  2^K m0 = 3^L m0 + RHS  (mod 3): 2^K m0 == RHS == 2^{S_{L-1}} = (-1)^{K-k_{L-1}}")
for m0 in [27,703,871,6171]:
    L=10; ks,mL=traj(m0,L); K,S,RHS=cycle_eq(m0,ks,L)
    lhs=(pow(2,K,3)*(m0%3))%3
    rhs3=RHS%3
    pred=pow(-1,ks[-1],3)%3 if False else (1 if ks[-1]%2==0 else 2)
    print(f" m0={m0:>4} mod3={m0%3}  k_last={ks[-1]}  2^K m0 mod3={lhs}  RHS mod3={rhs3}  match={lhs==rhs3}")

print("\n"+"="*72)
print("(3) 3-adic 'prefix' count: is the Syracuse map even well-defined mod 3^L? "
      "(the asymmetry vs the 2-adic prefix)")
print("="*72)
# For each residue r mod 3^L, does (3r+1)/2^v land in a single residue mod 3^L independent of lift?
for L in [2,3]:
    N=3**L; ambiguous=0; total=0
    for r in range(N):
        if r%3==0: continue  # coprime classes
        total+=1
        # collect (3r+1)*2^{-v} mod N over v=1..VMAX -- these are the POSSIBLE images
        imgs=set()
        base=(3*r+1)%N
        inv2=pow(2,-1,N)
        for v in range(1,2*3**(L-1)+1):
            imgs.add((base*pow(inv2,v,N))%N)
        if len(imgs)>1: ambiguous+=1
    print(f"  mod 3^{L} (N={N}): {total} coprime residues; {ambiguous} have >1 possible image "
          f"=> map NOT single-valued mod 3^L (v not determined by residue). frac={ambiguous/total:.3f}")
