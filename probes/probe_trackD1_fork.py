"""
PROBE D-1 -- THE FORK: H1 (finite Jordan cluster) vs H2 (essential-curve finite sections) + crossing rider.
F1 mode census N(delta,L); F2 phase arithmetic; F3 doublet anatomy; F4 crossing rider (lambda=0.4,0.6).
Direct methods. L=2,3 spectra from D1-A dump; L=4 top modes from D1-C block-SpMV log; F4 = small local power iter.
No fits anywhere. Report discriminators, not conclusions, where they conflict.
"""
import numpy as np, json
from probe_phase2a_q2b_q6 import build_M_gen, subgroup

THIRD = 1.0 / 3.0

# ---- spectra (measured; L2,3 = D1-A dense top-8, L4 = D1-C block-SpMV converged) ----
SPEC = {
    2: [0.3468266586, 0.022320651+0.235129293j, 0.022320651-0.235129293j,
        -0.161688749+0.124496698j, -0.161688749-0.124496698j,
        -0.025202143+0.172497226j, -0.025202143-0.172497226j, -0.151926512+0.025786248j],
    3: [0.3332363004, 0.237639959+0.183030417j, 0.237639959-0.183030417j,
        0.234998610+0.183154983j, 0.234998610-0.183154983j, 0.2731618323,
        0.066874799+0.258213742j, 0.066874799-0.258213742j],
    4: [0.333500000, 0.320423+0.075242j, 0.320423-0.075242j,
        0.320223+0.075252j, 0.320223-0.075252j],
}

def census():
    print("\n## F1  MODE CENSUS  N(delta,L) = #{lambda : |lambda - 1/3| < delta}")
    print("   (also: distance-to-1/3 of the leading complex pair -- the accumulation trend)")
    for d in (0.05, 0.02, 0.01):
        row = []
        for L in (2, 3, 4):
            n = sum(1 for lam in SPEC[L] if abs(complex(lam) - THIRD) < d)
            row.append(n)
        print(f"   delta={d}:  L2={row[0]}  L3={row[1]}  L4={row[2]}")
    print("   leading-pair distance to 1/3:  " +
          "  ".join(f"L{L}={abs(complex(SPEC[L][1])-THIRD):.4f}" for L in (2, 3, 4)) +
          f"   (ratios {abs(complex(SPEC[3][1])-THIRD)/abs(complex(SPEC[2][1])-THIRD):.3f}, "
          f"{abs(complex(SPEC[4][1])-THIRD)/abs(complex(SPEC[3][1])-THIRD):.3f}  -> toward 1/3?)")

def phases():
    print("\n## F2  PHASE ARITHMETIC  (leading complex pair; raw, no fit)")
    ph = [np.angle(complex(SPEC[L][1])) for L in (2, 3, 4)]
    print(f"   arg(pair):  L2={ph[0]:.5f}  L3={ph[1]:.5f}  L4={ph[2]:.5f} rad")
    print(f"   contraction ratios:  L3/L2={ph[1]/ph[0]:.4f}   L4/L3={ph[2]/ph[1]:.4f}   "
          f"(H2 signature = -> 1/3 = 0.3333; H1 = fixed-fan/irregular)")

def doublet():
    print("\n## F3  DOUBLET ANATOMY  (the two adjacent complex pairs)")
    for L in (3, 4):
        p1, p2 = complex(SPEC[L][1]), complex(SPEC[L][3])
        split = abs(p1 - p2)
        print(f"   L={L}: pair1={p1:.6f} (|.|={abs(p1):.5f})  pair2={p2:.6f} (|.|={abs(p2):.5f})  "
              f"splitting={split:.3e}")
    s3 = abs(complex(SPEC[3][1]) - complex(SPEC[3][3])); s4 = abs(complex(SPEC[4][1]) - complex(SPEC[4][3]))
    print(f"   L=2: leading pairs |.|=0.2362, 0.2041 (split 0.032, NO tight doublet yet)")
    print(f"   splitting L3={s3:.3e} -> L4={s4:.3e}  ratio {s4/s3:.4f}  "
          f"(H2 = shrinks with L = adjacent section modes; H1 = stable broken degeneracy)")

def build_tower_rho(L, lam):
    q = 3; qL = q ** L
    sub = subgroup(2 % qL, qL); D = len(sub)
    raw = [lam ** d for d in range(1, D + 1)]
    M, idx, n = build_M_gen(q, L, 2, raw)
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    gam = np.array([s[2] for s in states]); tw = np.where(gam != 0)[0]
    Mt = M[tw][:, tw].tocsr().astype(np.float64)
    # spectral radius via power iteration (gapped off-resonance -> fast)
    rng = np.random.default_rng(0); v = np.abs(rng.standard_normal(Mt.shape[0])); v /= np.linalg.norm(v)
    rho = 0.0
    for _ in range(4000):
        w = Mt.dot(v); rho_new = float(v @ w); v = w / np.linalg.norm(w)
        if abs(rho_new - rho) < 1e-13: rho = rho_new; break
        rho = rho_new
    wn = np.array(raw); wn = wn / wn.sum(); c0 = float(np.sum(wn ** 2))
    return rho, c0

def crossing():
    print("\n## F4  CROSSING RIDER  (c0 = Sum w_norm^2 vs partner = rho(M_tower); lambda=0.4,0.5,0.6)")
    print("   PRE-REG (SHAPE): c0 ~ (1-lambda)/(1+lambda)+fold; partner stays within fluctuation of 1/3.")
    for L in (2, 3):
        print(f"   L={L}:")
        for lam in (0.4, 0.5, 0.6):
            rho, c0 = build_tower_rho(L, lam)
            cf = (1 - lam) / (1 + lam)
            print(f"      lambda={lam}: c0={c0:.6f} ((1-l)/(1+l)={cf:.6f})  partner rho={rho:.6f}  "
                  f"|rho-1/3|={abs(rho-THIRD):.4f}  |c0-1/3|={abs(c0-THIRD):.4f}")

def main():
    print("# PROBE D-1 -- THE FORK (H1 vs H2) + crossing rider. Direct/from-data. No fits.")
    census(); phases(); doublet(); crossing()

if __name__ == "__main__":
    main()
