def syr_step(m):
    x = 3*m + 1
    k = (x & -x).bit_length() - 1   # v_2(3m+1)
    return x >> k, k

def telescope_check(m0, L):
    m = m0; ks = []
    for _ in range(L):
        m, k = syr_step(m); ks.append(k)
    mL = m; K = sum(ks)
    # CORRECT convention: S_i = sum of first i k's (halvings BEFORE step i), S_0=0
    RHS_correct = sum(3**(L-1-i) * 2**(sum(ks[:i])) for i in range(L))
    lhs = 2**K * mL
    rhs = 3**L * m0 + RHS_correct
    # DESKTOP convention: s_i = k_0+...+k_i (through step i)
    RHS_desktop = sum(3**(L-1-i) * 2**(sum(ks[:i+1])) for i in range(L))
    return mL, K, ks, lhs, rhs, (lhs == rhs), RHS_desktop

print("General telescoping identity  2^K*m_L == 3^L*m_0 + sum 3^{L-1-i} 2^{S_i}  (S_i = first i k's):")
print(f"{'m0':>6} {'L':>3} {'k_i':>20} {'2^K m_L':>14} {'3^L m0+RHS':>14} {'match':>6}")
import random
random.seed(1)
for (m0,L) in [(1,1),(3,5),(7,6),(27,10),(31,8),(2**20+1,12)] + [(random.randrange(1,10**6)*2+1, random.randint(3,15)) for _ in range(6)]:
    mL,K,ks,lhs,rhs,ok,rd = telescope_check(m0,L)
    print(f"{m0:>6} {L:>3} {str(ks):>20} {lhs:>14} {rhs:>14} {str(ok):>6}")

print("\nTRIVIAL CYCLE {1,4,2,1}: L=1, m0=1")
mL,K,ks,lhs,rhs,ok,rd = telescope_check(1,1)
RHS_corr = sum(3**(1-1-i)*2**(sum(ks[:i])) for i in range(1))
print(f"  k_i={ks}, K={K}, m_L={mL} (==m0 so cycle)")
print(f"  CORRECT  : (2^K - 3^L) m0 = {(2**K-3**1)*1}   RHS = sum 3^{{L-1-i}} 2^{{S_i}} = {RHS_corr}   -> {'MATCH' if (2**K-3**1)*1==RHS_corr else 'NO'}")
print(f"  DESKTOP  : RHS with s_i (through step i) = {rd}   vs LHS {(2**K-3**1)*1}   -> {'MATCH' if (2**K-3**1)*1==rd else 'MISMATCH (the bug)'}")
