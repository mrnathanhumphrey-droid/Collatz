# LIT DISPOSITION — is S_∞'s value-transcendence over ℚ(log₂3) a compute problem or a named wall? **NAMED FRONTIER (open with current mathematics).** (2026-09-06)

Follow-up to `result_PEEL_BRICK.md` (numerical A-vs-B test = OUTCOME C, below floor). Question: is the
"is `S_∞ = A + B·f(log₂3)` for algebraic `A,B`" (trdeg `ℚ(S_∞, log₂3)/ℚ` = 1 vs 2) question answerable in
principle — a compute problem — or a wall of the same character as the Bernoulli-convolution measure-type
wall the corpus already hit? Literature check (Wilson's call; option 3). Verdict: **it is a named frontier**,
not a compute or missing-lemma problem.

## The tool that WOULD answer it, and exactly why it can't reach S_∞

The machinery for "transcendence degree / algebraic independence of a value" is **Mahler's method**
(Nishioka's theorem, conjectured Mahler 1969, proved Nishioka 1986; general algebraic-independence criteria
by Philippon). It determines `trdeg_ℚ ℚ(f₁(α),…,fₘ(α))` — *exactly our question's shape*. **But it is
restricted to FINITE-order systems:** the method runs on **zero estimates**, and *"proving certain
conjectures remains beyond the reach of current techniques due to the lack of a general zero estimate theory
for nonlinear [/infinite] systems"* (Waldschmidt). Nishioka's transcendence holds for **linear** Mahler
equations (non-rational solution ⟹ natural boundary ⟹ transcendental).

**Our object is INFINITE Mahler order** — the corpus's own `result_MAHLER.md` proved `S` has no finite-order
rational Mahler recurrence at *any* order (denominator rate 3 > 2 > λ_L). So `S_∞` sits **precisely in the
gap** Mahler's method cannot cross. The theorem that walls the value from the analytic side (infinite Mahler
order) is the *same* theorem that puts it outside the algebraic-independence machinery.

## The nearest state-of-the-art result — and it stops exactly at our boundary

**arXiv:2512.14077 (Dec 2025), "Transcendence and algebraic independence of a family of p-adic valuation
generating functions"** — the closest cousin in the literature: it studies `T_p(z) = ∏_j (1−z^{p^j})^{−1/p^j}`
with `log T_p(z) = Σ ν_p(n)/n · z^n`, the **p-adic valuation** generating function (the same *kind* of object
as the Syracuse `÷p^v` valuation structure). Results:
- **Thm 1.1:** `T_p(α)` transcendental for algebraic `α`, `|α|<1`.
- **Thm 1.2:** `T_p(α_1),…,T_p(α_m)` **algebraically independent** at multiplicatively independent algebraic points.

Crucially, this works **because `T_p` satisfies a specific FINITE algebraic Mahler equation**
`T_p(z)^p = T_p(z^p)/(1−z^p)` — "beyond linear, but still a single finite (nonlinear) equation." **`S_∞` has
no finite algebraic Mahler equation at all (infinite order).** So even the cutting edge for
valuation-generating-function transcendence **stops exactly at the finite/infinite boundary that walls us** —
and it is "not connected to Collatz/Syracuse" (its substitution is `z↦z^p`, its coefficient the valuation of
`n`, not the Syracuse stationary mass). The door that would open the question is the one MAHLER closed.

## Independently: trdeg ≥ 2 for two specific constants is open even for e+π

Setting the infinite-order issue aside, the A-vs-B call is a **trdeg-1-vs-2 (algebraic-independence)**
statement about two *specific* constants. That is the hardest class: *"the transcendence of `e + π` is a
longstanding [open] problem."* The only cheap tool — *"two transcendental periods of DIFFERENT transcendence
degree have transcendental sum/ratio"* — doesn't decide the dimension. Under the **Grothendieck–André /
Kontsevich–Zagier period conjecture** the answer is conjecturally pinned (relations among periods are
"geometric in origin"), **but that conjecture is open**.

## Collatz-specific corroboration

**arXiv:1408.6884, "3x+1 inverse orbit generating functions almost always have natural boundaries"** —
establishes natural boundaries (hence transcendence) for the 3x+1 inverse-orbit *generating functions*.
Consistent with our infinite-Mahler / natural-boundary picture, but it concerns the *function*, not the
trdeg of an associated limiting *constant*.

## Verdict

The value-transcendence question — is `S_∞`'s transcendence reducible to the single generator `log₂3` — is a
**NAMED FRONTIER, open with current mathematics**, sitting at the intersection of three separately-open areas:
1. **infinite-order Mahler** (outside the Nishioka/Philippon zero-estimate reach — and MAHLER proved we are
   infinite order);
2. **algebraic independence of specific constants** (open even for `e+π`);
3. **the period conjecture** (Grothendieck–André, open).

So `PEEL_BRICK`'s numerical Outcome C is **matched by a theoretical verdict from the opposite direction**: the
A-vs-B distinction is not merely below the *data* floor — it is beyond the reach of current *mathematics*. Same
character as the Bernoulli-convolution measure-type wall (`result_GARSIA`/`DENOM`): a named-hard problem, not a
shortcut we missed. **Not a compute problem.** The capstone's "∞-dim unipotent kernel, value walled" stands;
the reducibility question is a frontier, not a gap.

**Sources**
- B. Adamczewski & C. Faverjon et al., *An introduction to Mahler's method for transcendence and algebraic independence* — https://arxiv.org/abs/1005.1216
- *Transcendence and algebraic independence of a family of p-adic valuation generating functions* (2025) — https://arxiv.org/abs/2512.14077
- *3x+1 inverse orbit generating functions almost always have natural boundaries* — https://arxiv.org/abs/1408.6884
- M. Waldschmidt, *Conjectures and open problems in transcendence* (AWS Lecture 5) — https://webusers.imj-prg.fr/~michel.waldschmidt/articles/pdf/AWSLecture5.pdf
- *Degrees of periods* — https://arxiv.org/abs/1102.2273
