"""
PROBE EULER-T2  --  is  mu = H1 - H2  (drift = Shannon - collision entropy) our pi^2/6, or trivial?

Express: valuation law P(v) = (1-r) r^{v-1}, v>=1  (Collatz: r=1/2 -> P(v)=2^{-v}).
  H1 (Shannon)   = -sum P log P
  H2 (collision) = -log sum P^2
  mu (drift)     = E[v] log p - log q       (map x -> (q x + 1)/p^v)
Natural coupling for the p-Hydra: r = 1/p (p-adic-valuation-like geometric).

Test 1: Collatz (q,p)=(3,2), r=1/2 -- verify mu = H1 - H2 exactly.
Test 2: GENERALIZE over (q,p) with r=1/p -- find the condition for equality (the cyclotomic-7 test:
        does it hold as a family law, or only at one point?).
Test 3: DECOMPOSE -- show mu = H1 - log q is DEFINITIONAL, and H2 = log q is the known criticality fact.
"""
import math

def law(r, vmax=200000):
    return [ (1-r)*r**(v-1) for v in range(1, vmax) ]

def H1(r):   # Shannon, nats
    P = law(r); return -sum(p*math.log(p) for p in P if p>0)
def H2(r):   # collision (Renyi-2), nats
    P = law(r); return -math.log(sum(p*p for p in P))
def Ev(r):   # mean valuation
    P = law(r); return sum(v*p for v,p in enumerate(P, start=1))
def mu(q, p, r):  # drift
    return Ev(r)*math.log(p) - math.log(q)

# closed forms (derived): with r general,
#   H1 = -log(1-r) - (r/(1-r)) log r
#   H2 = log((1+r)/(1-r))
#   H1 - H2 = -(r/(1-r)) log r - log(1+r)
# with r = 1/p:
#   H1 - H2 = (p/(p-1)) log p - log(p+1)
#   mu      = (p/(p-1)) log p - log q       => EQUAL iff  q = p+1.

def H1_closed(r): return -math.log(1-r) - (r/(1-r))*math.log(r)
def H2_closed(r): return math.log((1+r)/(1-r))

print("# PROBE EULER-T2\n")
print("## Test 1 -- Collatz (q,p)=(3,2), r=1/2")
r=0.5
h1,h2,m = H1(r), H2(r), mu(3,2,r)
print(f"   H1 (Shannon)   = {h1:.12f}   (log 4 = {math.log(4):.12f})")
print(f"   H2 (collision) = {h2:.12f}   (log 3 = {math.log(3):.12f})")
print(f"   mu (drift)     = {m:.12f}   (log(4/3) = {math.log(4/3):.12f})")
print(f"   H1 - H2        = {h1-h2:.12f}")
print(f"   mu - (H1-H2)   = {m-(h1-h2):.2e}   <-- identity holds" )
print(f"   closed forms:  H1={H1_closed(r):.12f}  H2={H2_closed(r):.12f}\n")

print("## Test 2 -- GENERALIZE with r=1/p: does mu = H1-H2 hold as a FAMILY LAW?  (the cyclotomic-7 test)")
print("   row: p, q, r=1/p, H1-H2, mu, equal?, and q vs p+1")
for p in [2,3,4,5,7]:
    r = 1.0/p
    hg = H1(r) - H2(r)
    for q in [p+1, p+2, 2, 3]:
        m = Ev(r)*math.log(p) - math.log(q)
        eq = abs(m-hg) < 1e-9
        tag = "  <== q=p+1" if q==p+1 else ""
        if q in (p+1, p+2) or (p==2 and q==3):
            print(f"   p={p} q={q:<2} r={r:.4f}  H1-H2={hg:.9f}  mu={m:.9f}  equal={eq}  (p+1={p+1}){tag}")
    print()

print("## Test 3 -- DECOMPOSE the identity into known pieces")
r=0.5
print(f"   (a) H1 = E[v]*log2 = mean log-denominator : H1={H1(0.5):.9f}  E[v]*log2={Ev(0.5)*math.log(2):.9f}")
print(f"       so mu = E[v]log2 - log q = H1 - log q  (DEFINITIONAL). H1 - log3 = {H1(0.5)-math.log(3):.9f} = mu")
print(f"   (b) H2 = log q at q=3 (collision-entropy criticality, already banked): H2={H2(0.5):.9f}  log3={math.log(3):.9f}")
print(f"   => mu = H1 - H2  =  (H1 - log q)_definitional  +  (log q - H2)_=0-by-criticality")
print(f"      The identity is TRUE but factors entirely into one definitional + one known-criticality fact.")
print(f"      H1-H2 = (p/(p-1))log p - log(p+1) is a PURE denominator-side function; it equals the")
print(f"      q-containing drift mu ONLY because Collatz has q = p+1 = 3. Not a family law. Leading-layer only.")
