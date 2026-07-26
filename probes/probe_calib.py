"""
PROBE CALIB (T1) -- calibrate sup|pi-hat| decay against Tao's n^-A (2026-07-26).

Data = the CERTIFIED sup|pi-hat| (LINEAR) values from result_LAMBDA.md (k=3..15). No recompute (uses data we have).

Wilson's model: Tao Prop 1.17 => sup(n) ~ C n^-A. Per-step ratio Srate = sup(n+1)/sup(n) = (n/(n+1))^A, so
    A_eff(n) = ln(Srate) / ln(n/(n+1)).
Predictions to test:
  - A_eff ~ 2.3-2.5 CREEPING UP  => superpolynomial signature (a FIXED power gives CONSTANT fitted A).
  - exp(-c*n) [pure exponential, Remark 1.15 exp(-cm)] gives a CONSTANT Srate; ours climbs => NOT pure exp.
    Fit exp(-c*n^alpha): Wilson expects alpha ~ 0.14 (near-power-law, slowly rising exponent).
Deliverable: A_eff per window across all k; does it converge or keep climbing? (decides T5 / the exp(-cm) reconciliation).
"""
import numpy as np

# sup|pi-hat| (LINEAR), certified in result_LAMBDA.md (k=3..15)
K = np.array([3,4,5,6,7,8,9,10,11,12,13,14,15])
SUP = np.array([0.252237,0.176999,0.129274,0.096106,0.075870,0.060891,
                0.048026,0.038278,0.031944,0.026458,0.022052,0.019128,0.016284])


def main():
    print("# PROBE CALIB (T1) -- sup|pi-hat| vs Tao n^-A ; A_eff per window\n")

    # --- per-step A_eff = ln(Srate)/ln(n/(n+1)) ---
    print("## per-step A_eff(n->n+1) = ln(sup_{n+1}/sup_n) / ln(n/(n+1))")
    print(f"   {'n=k':>4} {'Srate':>8} {'A_eff':>8}")
    A_eff = []
    for i in range(len(K)-1):
        n = K[i]
        srate = SUP[i+1]/SUP[i]
        a = np.log(srate)/np.log(n/(n+1))
        A_eff.append((n, a))
        print(f"   {n:>4} {srate:>8.4f} {a:>8.4f}")
    Aeff_arr = np.array([a for _, a in A_eff])
    print(f"\n   A_eff range {Aeff_arr.min():.3f}..{Aeff_arr.max():.3f}; first->last {Aeff_arr[0]:.3f}->{Aeff_arr[-1]:.3f}; "
          f"trend slope d(A_eff)/dk = {np.polyfit([n for n,_ in A_eff], Aeff_arr, 1)[0]:+.4f}")
    print("   [climbing => superpolynomial (fixed power => constant A_eff). Wilson: ~2.3-2.5 creeping up.]")

    # --- windowed power-law fit: ln sup = c - A ln n, sliding 4-pt windows ---
    print("\n## windowed power-law fit  ln(sup) = c - A*ln(n)  (4-pt sliding windows)")
    print(f"   {'window k':>10} {'A_fit':>8} {'R^2':>8}")
    for i in range(len(K)-3):
        ks = K[i:i+4]; ss = SUP[i:i+4]
        A_fit, c = np.polyfit(np.log(ks), np.log(ss), 1)
        pred = c + A_fit*np.log(ks)
        r2 = 1 - np.sum((np.log(ss)-pred)**2)/np.sum((np.log(ss)-np.log(ss).mean())**2)
        print(f"   {ks[0]:>4}-{ks[-1]:<5} {-A_fit:>8.4f} {r2:>8.5f}")
    print("   [A_fit rising across windows = the same creep, from the fit side.]")

    # --- pure-exponential test: is Srate constant? (exp(-c*n) => constant Srate) ---
    srates = SUP[1:]/SUP[:-1]
    print(f"\n## pure-exp test: Srate should be CONSTANT if sup~exp(-c*n).  "
          f"Srate {srates.min():.3f}..{srates.max():.3f} (climbs {srates[0]:.3f}->{srates[-1]:.3f}) => NOT pure exp.")

    # --- stretched-exp fit ln sup = b - c*n^alpha : grid-search alpha ---
    print("\n## stretched-exp fit  ln(sup) = b - c*n^alpha  (grid over alpha; best R^2)")
    y = np.log(SUP); best = None
    for alpha in np.arange(0.05, 1.61, 0.01):
        x = K.astype(float)**alpha
        slope, intercept = np.polyfit(x, y, 1)
        pred = intercept + slope*x
        r2 = 1 - np.sum((y-pred)**2)/np.sum((y-y.mean())**2)
        if best is None or r2 > best[1]:
            best = (alpha, r2, -slope)
    print(f"   best alpha = {best[0]:.2f}  (R^2={best[1]:.6f}, c={best[2]:.4f})")
    # also the pure power-law (alpha->0 limit is ln n; do it explicitly) and pure-exp alpha=1 for compare
    for al in (1.0,):
        x = K.astype(float)**al; s,ic = np.polyfit(x,y,1); pr=ic+s*x
        r2=1-np.sum((y-pr)**2)/np.sum((y-y.mean())**2)
        print(f"   alpha=1.00 (pure exp) R^2={r2:.6f};  pure power-law ln sup~-A ln n "
              f"R^2={1-np.sum((y-(np.polyfit(np.log(K),y,1)[1]+np.polyfit(np.log(K),y,1)[0]*np.log(K)))**2)/np.sum((y-y.mean())**2):.6f}")
    print(f"   [Wilson predicted alpha~0.14 = near-power-law w/ slowly rising exponent; alpha<<1 corroborates superpoly not exp.]")

    # --- verdict on tower-vs-factorial worry ---
    print(f"\n## USEFUL-A verdict: A_eff sits ~{np.median(Aeff_arr):.1f} at these depths; grind A=2,3,4 and STOP.")
    print("   [tower(2.5) is nothing; blow-up only at large A, which we don't need. Kills the tower-vs-factorial worry.]")


if __name__ == "__main__":
    main()
