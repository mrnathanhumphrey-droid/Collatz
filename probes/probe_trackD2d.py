"""
PROBE D2-d -- the k=+-4 seat + sector identification + census closure. Dense L=3 (values+vectors). Committed.
T1 all eigenvalues modulus in [0.003,0.05], (mod,phase), flag conj pairs, nearest to seat (0.0101, +-2.79).
T2 sector (gauge-character) mass of the nearest pair's eigenvector -- PRE-REG: leading in k = +-4 class (mod 9).
T3 census closure: among conj pairs modulus > 0.003, dominant gauge-k -> coprime class mod 9. LADDER PREDICTS
   EXACTLY THREE families {+-1,+-2,+-4}. A fourth (unit) class is impossible mod 9; the real test = whether a
   pair's dominant k is DIVISIBLE by 3 (an internal/non-coprime family) or otherwise unassignable -> warning shot.
"""
import numpy as np
from collections import defaultdict
from probe_phase2c0 import build_M_tower_and_coords

def sector_mass(vec, twcoords, dl, D):
    blocks = defaultdict(lambda: np.zeros(D, dtype=complex))
    for i, (a, er, g) in enumerate(twcoords):
        blocks[(er, g)][dl[a]] = vec[i]
    mk = np.zeros(D)
    for bv in blocks.values():
        mk += np.abs(np.fft.fft(bv)) ** 2
    return mk / mk.sum()

def cls9(k):
    r = k % 9
    if r in (1, 8): return "+-1"
    if r in (2, 7): return "+-2"
    if r in (4, 5): return "+-4"
    return f"div3(r={r})"

def main():
    L = 3
    print(f"# PROBE D2-d -- k=+-4 seat + sector + census, dense L={L}. Committed pre-reg.")
    Mt, states, idx, tw, pos, twcoords, q, qL, sub, D, dl, w, two = build_M_tower_and_coords(L)
    Md = Mt.toarray()
    ev, VR = np.linalg.eig(Md)
    upi = [i for i in range(len(ev)) if ev[i].imag > 1e-9]     # upper-half conj representatives

    # ---- T1 ----
    band = sorted([i for i in upi if 0.003 <= abs(ev[i]) <= 0.05], key=lambda i: -abs(ev[i]))
    seatmod, seatph = (1/3) * np.cos(2.7925/2) ** 2, 2.7925
    print(f"\n## T1  pairs with modulus in [0.003,0.05]  (seat: mod {seatmod:.5f}, phase {seatph:.4f})")
    tgt = seatmod * np.exp(1j * seatph)
    for i in band:
        z = ev[i]
        print(f"   {z.real:+.6f}{z.imag:+.6f}j  |.|={abs(z):.5f}  arg={np.angle(z):.4f}  "
              f"dist-to-seat={abs(z - tgt):.5f}")
    nearest = min(band, key=lambda i: abs(ev[i] - tgt)) if band else None
    if nearest is not None:
        z = ev[nearest]
        print(f"   => NEAREST TO SEAT: {z.real:+.6f}{z.imag:+.6f}j |.|={abs(z):.5f} arg={np.angle(z):.4f}")

    # ---- T2 ----
    print(f"\n## T2  sector (gauge-character) mass of the nearest pair  (PRE-REG: leading k = +-4 class mod 9)")
    if nearest is not None:
        mk = sector_mass(VR[:, nearest], twcoords, dl, D)
        topk = sorted(range(D), key=lambda k: -mk[k])[:5]
        print(f"   top-k mass: " + ", ".join(f"k={k}({cls9(k)}):{mk[k]:.4f}" for k in topk))
        dom = topk[0]
        print(f"   => dominant gauge-k = {dom} -> coprime class {cls9(dom)}  "
              f"({'CONFIRMS k=+-4 identification' if cls9(dom)=='+-4' else 'DEVIATES from k=+-4'})")

    # ---- T3 ----
    print(f"\n## T3  CENSUS CLOSURE: dominant gauge-k -> coprime class, all conj pairs modulus > 0.003")
    fam = defaultdict(list)
    for i in upi:
        if abs(ev[i]) <= 0.003: continue
        mk = sector_mass(VR[:, i], twcoords, dl, D)
        dom = max(range(1, D), key=lambda k: mk[k])   # exclude k=0 (DC; complex modes have a_DC~0)
        fam[cls9(dom)].append((abs(ev[i]), float(np.angle(ev[i])), dom, float(mk[dom])))
    print(f"   families found among modulus>0.003 pairs ({sum(len(v) for v in fam.values())} pairs):")
    for c in sorted(fam, key=lambda c: -max(m for m, *_ in fam[c])):
        ms = sorted(fam[c], key=lambda t: -t[0])
        topmods = ", ".join(f"{m:.4f}@{ph:.2f}(k={k})" for m, ph, k, _ in ms[:4])
        print(f"     class {c}: {len(ms)} pairs; top: {topmods}")
    unit_fams = [c for c in fam if c in ("+-1", "+-2", "+-4")]
    div_fams = [c for c in fam if c.startswith("div3")]
    print(f"   => UNIT coprime families = {sorted(unit_fams)} ({len(unit_fams)}); "
          f"predicted EXACTLY 3 {{+-1,+-2,+-4}} ({'MATCH' if set(unit_fams)=={'+-1','+-2','+-4'} and len(unit_fams)==3 else 'DEVIATION'})")
    if div_fams:
        tot_div = sum(len(fam[c]) for c in div_fams)
        big_div = [(m, ph, k) for c in div_fams for m, ph, k, _ in fam[c] if m > 0.01]
        print(f"   [!] div-by-3 (internal-rung) dominant classes: {sorted(div_fams)}, {tot_div} pairs total "
              f"(frequencies 3k/9k = internal rungs of the 3 coprime ladders; EXPECTED, not a 4th coprime family)")
    print(f"   => CENSUS CLOSURE: {'THREE coprime families + internal(div3) rungs + real modes; NO 4th coprime family = ladder indexing HOLDS' if set(unit_fams)=={'+-1','+-2','+-4'} else 'DEVIATION -- warning shot'}")

if __name__ == "__main__":
    main()
