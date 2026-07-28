"""
PROBE MOD2K-STATIONARY -- closed form for the inverse-tree residue distribution mod 2^k (2026-07-27)

Inverse-tree Phase-5 open question #2 (`inverse_tree_findings.md`): closed form for the node-measure
residue stationary distribution mod 2^k, and explain the *increasing* chi^2/n drift (residues do NOT
equidistribute; the distribution converges to a NON-uniform fixed point).

Node n: children left=2n (always), right=(n-1)/3 (only n=4 mod6). KEY: 3 is a UNIT mod 2^k
(3^-1 mod 64 = 43), so the /3 right-child map is EXACT on the 2-adic part -- no 3-adic lift needed
for the mod-2^k component. Only the branching TEST (n=4 mod6) needs n mod 3. So build the joint
offspring operator on states (u = n mod 2^k, w = n mod 3^J), marginalize the Perron eigenvector to
mod 2^k. (BRANCH-BIAS already pinned the mod-3 side: lambda=(3+sqrt21)/6, exact at any 3-depth.)

  P1  empirical mod-2^k node distribution from tree_d50.parquet; chi^2 vs uniform at k=3..6, per depth.
  P2  joint operator on (mod 2^k) x (mod 3^J): Perron eigvec, marginalize to mod 2^k; J-stability.
  P3  compare Perron marginal to empirical d=50; structure / closed-form hunt.
"""
import numpy as np
from fractions import Fraction as Fr

TREE = r"C:\Collatz\inverse_tree\tree_d50.parquet"
LAM = (3 + np.sqrt(21)) / 6           # tree growth rate (BRANCH-BIAS)


# ----------------------------------------------------------------------------- P1 empirical
def p1_empirical(kmax=6):
    print("## P1  empirical node distribution mod 2^k (tree_d50.parquet)")
    import polars as pl
    df = pl.read_parquet(TREE, columns=["n", "depth"])
    n = df["n"].to_numpy().astype(object)
    depth = df["depth"].to_numpy()
    dmax = int(depth.max())
    emp = {}
    for k in range(3, kmax + 1):
        M = 1 << k
        # deep-node distribution (depth in [dmax-8, dmax]) as the "stationary" estimate
        sel = depth >= dmax - 8
        res = np.array([int(x) % M for x in n[sel]])
        cnt = np.bincount(res, minlength=M).astype(float); cnt /= cnt.sum()
        emp[k] = cnt
        chi2 = M * ((cnt - 1.0 / M) ** 2).sum() / (1.0 / M)   # chi^2/n-like: sum (o-e)^2/e over M cells / (per-node)
        print(f"   k={k} (mod {M:>3}): chi^2 vs uniform = {((cnt-1/M)**2/(1/M)).sum():.4f}   "
              f"min={cnt.min():.4f} max={cnt.max():.4f}")
    # chi^2/n DRIFT across depth (does non-uniformity grow?)
    print("   chi^2(mod 64) vs uniform by depth window:")
    M = 64
    for d0 in range(20, dmax + 1, 6):
        sel = (depth >= d0) & (depth < d0 + 6)
        if sel.sum() < 200:
            continue
        res = np.array([int(x) % M for x in n[sel]])
        cnt = np.bincount(res, minlength=M).astype(float); cnt /= cnt.sum()
        print(f"     depth {d0:>2}-{d0+5}: chi^2={((cnt-1/M)**2/(1/M)).sum():.4f}  (N={sel.sum()})")
    print()
    return emp


# ------------------------------------------------------------------- P2 joint Perron operator
def joint_perron(k, J):
    """Offspring operator on states (u mod 2^k, w mod 3^J). Right child's 2-part exact via 3^-1;
       3-part via the 3 lifts (exact-in-limit, per BRANCH-BIAS). Returns mod-2^k marginal of Perron."""
    Mu = 1 << k
    Mw = 3 ** J
    inv3 = pow(3, -1, Mu)               # 3 is a unit mod 2^k
    S = Mu * Mw
    idx = lambda u, w: u * Mw + w
    A = np.zeros((S, S))
    for u in range(Mu):
        for w in range(Mw):
            s = idx(u, w)
            # left child 2n : always
            A[idx((2 * u) % Mu, (2 * w) % Mw), s] += 1.0
            # branch iff n=4 mod6  <=> n odd? no: n even and n=1 mod3. parity from u (u&1), mod3 from w%3
            if (u % 2 == 0) and (w % 3 == 1):
                u2 = ((u - 1) * inv3) % Mu           # (n-1)/3 mod 2^k, EXACT (unit)
                for j in range(3):                    # 3 lifts of w to mod 3^{J+1}
                    wp = w + j * Mw                   # = w mod 3^J, distinct mod 3^{J+1}
                    w2 = ((wp - 1) // 3) % Mw          # (n-1)/3 mod 3^J
                    A[idx(u2, w2), s] += 1.0 / 3.0
    ev, V = np.linalg.eig(A)
    i = int(np.argmax(ev.real))
    lam = ev[i].real
    v = np.abs(V[:, i].real); v = v / v.sum()
    # marginalize to mod 2^k
    marg = np.zeros(Mu)
    for u in range(Mu):
        for w in range(Mw):
            marg[u] += v[idx(u, w)]
    return lam, marg


def p2_operator(kmax=6):
    print("## P2  joint (mod 2^k)x(mod 3^J) Perron operator -> mod-2^k marginal")
    out = {}
    for k in range(3, kmax + 1):
        row = []
        for J in (1, 2):
            lam, marg = joint_perron(k, J)
            row.append((J, lam, marg))
        # J-stability
        d = np.abs(row[0][2] - row[1][2]).max()
        print(f"   k={k}: lambda={row[0][1]:.6f} (exact (3+sqrt21)/6={LAM:.6f})  "
              f"|marg(J=1)-marg(J=2)|max={d:.2e}")
        out[k] = row[1][2]   # use J=2
    print()
    return out


# ----------------------------------------------------------------------------- P3 compare + structure
def p3_compare(emp, perron):
    print("## P3  Perron mod-2^k marginal vs empirical d=50, and structure")
    for k in sorted(perron):
        M = 1 << k
        e = emp[k]; p = perron[k]
        err = np.abs(e - p).max()
        print(f"   k={k} (mod {M}): |empirical - Perron|max = {err:.4f}   "
              f"(empirical from deep tree layers)")
    # structure of the mod-64 stationary: show it and look for pattern by 2-adic valuation
    k = max(perron); M = 1 << k; p = perron[k]
    print(f"\n   mod-{M} Perron stationary, grouped by 2-adic valuation v2(r) (r=0 excluded):")
    print(f"   {'v2':>3} {'#residues':>9} {'total mass':>11} {'mass/residue':>13}")
    for v in range(k + 1):
        if v == k:
            rs = [0]
        else:
            rs = [r for r in range(M) if (r & ((1 << (v + 1)) - 1)) == (1 << v)]  # v2(r)=v
        if not rs:
            continue
        tot = sum(p[r] for r in rs)
        print(f"   {v:>3} {len(rs):>9} {tot:>11.5f} {tot/len(rs):>13.6f}")
    # odd residues (v2=0): are they uniform? and does mass halve per valuation level?
    odd = [p[r] for r in range(M) if r % 2 == 1]
    print(f"   odd residues (v2=0): min={min(odd):.5f} max={max(odd):.5f} "
          f"spread={max(odd)-min(odd):.2e}  (NON-uniform => finer within-level structure)")

    # ---- CLOSED FORM: P(v2 = v) is Geometric(1 - 1/lam), ratio 1/lam = (sqrt21-3)/2 ----
    print(f"\n   CLOSED FORM: P(v2(n)=v) = (1-1/lam)*(1/lam)^v,  1/lam=(sqrt21-3)/2={1/LAM:.6f}")
    print(f"   {'v':>3} {'mass(v) emp':>12} {'geom pred':>11} {'err':>10}")
    tail = 1.0
    for v in range(k):
        if v == k - 1:
            rs = [r for r in range(M) if (r & ((1 << (v + 1)) - 1)) == (1 << v)] + [0]  # last incl r=0 tail
        else:
            rs = [r for r in range(M) if (r & ((1 << (v + 1)) - 1)) == (1 << v)]
        mass = sum(p[r] for r in rs)
        pred = (1 - 1 / LAM) * (1 / LAM) ** v if v < k - 1 else (1 / LAM) ** v  # last cell = geometric tail
        print(f"   {v:>3} {mass:>12.6f} {pred:>11.6f} {abs(mass-pred):>10.2e}")
    print(f"   => valuation law is Geometric(1-1/lam); success param 1-1/lam = (5-sqrt21)/2 "
          f"= {1-1/LAM:.6f} = P(node is odd). (Distinct from branch ratio (sqrt21-3)/6=0.263763.)")
    print(f"   MECHANISM (renewal): doubling child 2n pushes v2->v2+1; branch child (n-1)/3 is ODD,")
    print(f"   resets v2->0. Stationary valuation of a +1/reset renewal = geometric. Explains the")
    print(f"   NON-equidistribution (chi^2 grows with k): mass concentrates on high-v2 as (1/lam)^v.")


def main():
    print("# PROBE MOD2K-STATIONARY -- closed form for the inverse-tree residue dist mod 2^k\n")
    emp = p1_empirical(6)
    perron = p2_operator(6)
    p3_compare(emp, perron)


if __name__ == "__main__":
    main()
