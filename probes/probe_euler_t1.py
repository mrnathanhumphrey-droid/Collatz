"""
PROBE EULER-T1 -- the prize: is the inverse-tree lambda = (3+sqrt21)/6 linked to the Plancherel S_inf?

lambda is a degree-2 algebraic (3 l^2 - 3 l - 1 = 0), the Perron growth rate of the BACKWARD tree
(finite offspring matrix -- an associated-graded / finite-place object).
S_inf is the forward Plancherel ladder LIMIT: infinite Mahler depth, doubly-exp denominators, no closed
form, leading rational 7/15 OVERSHOT (S_16 = 0.471352 > 7/15). A pro-limit / tail / archimedean object.

Sharp algebraic form of a relation: S_inf = a + b*lambda, a,b in Q (since [Q(lambda):Q]=2).

HARD constraints (copied from ledger; re-grep before penning):
  exact floor   2*T_20   = 0.473177     (result_SOLSTICE) -- S_inf STRICTLY ABOVE this
  exact anchor  S_16=2T15 = 0.471352     (result_P6I)
  bracket       S_inf in [0.4714, 0.478] (result_SOLSTICE)
  7/15 = 0.466667 EXCLUDED (overshoot proven)
  best point estimate ~ 0.475

Test A: the two structurally-natural candidates  S_inf = 2 log lambda  and  S_inf = lambda - 1/lambda.
Test B: density -- how many low-complexity a+b*lambda land in the bracket? (if many, numerics cannot decide)
"""
import math
from fractions import Fraction as F

sqrt21 = math.sqrt(21)
lam = (3 + sqrt21) / 6
inv = 1/lam
loglam = math.log(lam)

FLOOR   = 0.473177     # 2*T_20 exact (rounded copy)
S16     = 0.471352     # exact
BR_LO, BR_HI = 0.4714, 0.478
S_pt    = 0.475

print("# PROBE EULER-T1\n")
print(f"lambda = (3+sqrt21)/6 = {lam:.9f}   1/lambda = {inv:.9f}   log lambda = {loglam:.9f}")
print(f"check 3l^2-3l-1 = {3*lam*lam - 3*lam - 1:.2e}")
print(f"S_inf: floor(2*T_20)={FLOOR}  S16={S16}  bracket=[{BR_LO},{BR_HI}]  7/15={7/15:.6f} EXCLUDED\n")

print("## Test A -- the two structurally-natural candidates, against the EXACT FLOOR")
c1 = 2*loglam
c2 = lam - inv
print(f"   candidate  S_inf = 2 log lambda       = {c1:.9f}   (motiv: log lambda ~ lim T_i)")
print(f"      vs floor {FLOOR}: {'ABOVE (survives)' if c1>FLOOR else 'BELOW FLOOR -> REFUTED'}  (gap {c1-FLOOR:+.6f})")
print(f"   candidate  S_inf = lambda - 1/lambda  = {c2:.9f}   (= 1 - 2/(3 lambda), clean algebraic)")
print(f"      vs floor {FLOOR}: {'ABOVE (survives)' if c2>FLOOR else 'BELOW FLOOR -> REFUTED'}  (gap {c2-FLOOR:+.6f})")
print(f"   => both most-natural lambda-candidates fall BELOW the exact rational floor. Refuted, not by")
print(f"      a fit, but by a proven exact lower bound.\n")

print("## Test B -- density: low-complexity  a + b*lambda  inside the bracket (numerics cannot discriminate)")
hits = []
for bd in range(0, 7):
    for bn in range(-6, 7):
        b = F(bn, bd) if bd else F(0)
        if bd and math.gcd(abs(bn), bd) != 1:
            continue
        for ad in range(1, 9):
            for an in range(-9, 19):
                if math.gcd(abs(an), ad) != 1:
                    continue
                a = F(an, ad)
                val = float(a) + float(b)*lam
                if BR_LO <= val <= BR_HI:
                    cx = ad + (bd if bd else 1) + abs(an) + abs(bn)  # crude complexity
                    hits.append((cx, float(a), float(b), val, f"{a}+{b}*lam"))
hits.sort()
print(f"   found {len(hits)} rationals a+b*lambda in [{BR_LO},{BR_HI}] with tiny numerators/denominators.")
print(f"   the 12 lowest-complexity:")
for cx,a,b,val,s in hits[:12]:
    above = "above floor" if val>FLOOR else "BELOW floor"
    print(f"     {s:<16} = {val:.6f}   ({above})")
print(f"   => the bracket (width {BR_HI-BR_LO}) is DENSE in low-complexity lambda-combinations; a numerical")
print(f"      'match' carries no information. Any T1 relation must be STRUCTURAL, not numerical.\n")

print("## Structural note (for the pen)")
print("   lambda: degree-2 Perron root of the finite backward-tree offspring matrix = associated-graded /")
print("           finite-place object (the SYMMETRIC, law-governed side of the structure/value split).")
print("   S_inf : forward Plancherel LIMIT, infinite Mahler depth, value in the tail = the pro-limit /")
print("           archimedean object no framework delivers (the VALUE side of the split).")
print("   They sit on OPPOSITE sides of the corpus's one established split. And duality_S_vs_D_verdict")
print("   already found NO clean forward<->inverse functional bridge (D*S, D+S, D/S all fail; D_n(1)/S_1=1/3")
print("   is 'a trivial coincidence'). So T1 relates the inverse tree's growth constant to the forward")
print("   ladder's limit across a split the corpus has never bridged and a duality already disposed.")
