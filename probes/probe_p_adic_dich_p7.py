"""
p-ADIC-DICH-p7 (Wilson spec, 2026-09-06) -- confirm the q^2=1+p fusion mechanism at a SECOND prime.
Reuses the validated p=5 machinery from probe_p_adic_dich (change ONLY p=7).  Pre-registered:
biconditional FAILS, forward half (7|k=>enriched) holds, cleanest MISS at k=q^2=4 (miss-LOCATION is
the sharp test), forward fiber-depth growth gamma_inf(7)<gamma_inf(14).
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import probe_p_adic_dich as P

t0 = time.time()
print("# p-ADIC-DICH-p7 -- q^2=1+p fusion mechanism, second-prime confirmation\n")

# ---- GATE 0: white baseline FIRST (spec) ----
N = 7 ** 3
gw = P.gammas(np.ones(N) / N)
bad = max(abs(gw[k] - 1.0) for k in range(1, P.KMAX + 1))
print(f"## GATE 0 white baseline p=7 r=3: max|gamma_white(k)-1| = {bad:.2e}  {'OK' if bad<1e-12 else 'FAIL -- STOP'}")
if bad >= 1e-12:
    sys.exit(1)
print()

# ---- primary p=7 dichotomy (exact r=2..5) ----
v7, tab7, gam7, nus7 = P.dichotomy_table(7, 2, 5, label="PRIMARY p=7")

# ---- the SHARP checks ----
print(f"\n## SHARP CHECKS (pre-registered)")
g4 = tab7[4][1]
print(f"   miss-LOCATION: gamma_inf(4) = {g4:.5f}  ({'ENRICHED' if g4>1 else 'depleted'})   "
      f"[7-nmid 4; PREDICT enriched = the q^2 orbit-axis miss at the SAME k=4 as p=5]")
# is k=4 the cleanest non-multiple miss?  list all enriched non-multiples
enr_nonmult = [k for k in range(1, P.KMAX + 1) if k % 7 != 0 and tab7[k][2]]
print(f"   enriched non-multiples of 7: {enr_nonmult}   (k=4=q^2 present? {4 in enr_nonmult})")
g7, g14 = tab7[7][1], tab7[14][1]
print(f"   forward fiber-depth growth: gamma_inf(7) = {g7:.5f} , gamma_inf(14) = {g14:.5f}   "
      f"[PREDICT gamma_inf(7) < gamma_inf(14)]  -> {'HOLDS' if g7 < g14 else 'FAILS'}")

# ---- tower (engine) exactness at p=7, the dissociation check ----
print(f"\n## ENGINE (tower telescoping) EXACT check p=7, r=3,4:")
for r in (3, 4):
    ok = P.tower_fold_exact(nus7[r], nus7[r - 1], 7, r)
    print(f"   Sum_j rho_{r}(t+j*7^{r-1}) == rho_{r-1}(t) exactly: {ok}")

# ---- verdict ----
n_pk_dep = sum(1 for k in tab7 if tab7[k][0] and not tab7[k][2])
n_np_enr = sum(1 for k in tab7 if not tab7[k][0] and tab7[k][2])
n_pk_enr = sum(1 for k in tab7 if tab7[k][0] and tab7[k][2])
fails = (n_np_enr > 0 or n_pk_dep > 0)
miss_at_4 = (4 in enr_nonmult)
print(f"\n## VERDICT")
if fails and n_pk_dep == 0 and miss_at_4:
    print(f"   CONFIRMS: biconditional FAILS (forward {n_pk_enr}/{n_pk_enr+n_pk_dep} holds, converse breaks "
          f"with {n_np_enr} non-mult enriched), cleanest miss at k=4=q^2 => q^2=1+p mechanism CONFIRMED "
          f"at a SECOND prime.")
elif not fails:
    print(f"   SURPRISE: biconditional HOLDS at p=7 (would need another coincidence at (7,2)). Report.")
elif fails and not miss_at_4:
    print(f"   SURPRISE: biconditional fails but k=4 NOT among enriched non-multiples {enr_nonmult} "
          f"=> q^2=1+p mechanism INCOMPLETE. Flag loudly.")
else:
    print(f"   MIXED: forward half broken (7|k depleted count {n_pk_dep}). Report exactly.")
print(f"\n# ({time.time()-t0:.1f}s)")
