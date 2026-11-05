# Findings

Every version of the hypothesis that was tested, in the order it was tested, with what
it predicted and what the data said. Scores are Δχ² against the standard model (ΛCDM):
negative is better, Δχ² ≈ −4 is about 2σ.

Data throughout: SPARC (175 galaxies, 3,391 points), Pantheon+SH0ES (1,657 light curves,
1,508 supernovae, full stat+sys covariance), DESI DR2 BAO (13 distances), Planck 2018
compressed CMB distance priors, and 2M++ predicted peculiar velocities for 450 hosts.

---

## Round 0 — the original proposal

> Every cubic metre of empty space grows by a small fixed fraction of itself. Across the
> large empty stretches in a galaxy's sparse outskirts this adds up and changes how fast
> the outer stars move. Between galaxies, expanding voids press inward, so the local
> Hubble rate varies from place to place while averaging to the global value.

Written as an equation, dV/dt = αV, which solves to exponential growth. **This is exactly
how a cosmological constant behaves.** The idea and the standard model's dark energy are
the same idea, so the cosmology half of the proposal was already correct, and already
known.

Fitting Pantheon+SH0ES gives the rate:

| Quantity | Value |
|---|---|
| H₀ | 73.53 ± 1.02 km s⁻¹ Mpc⁻¹ |
| Ωm | 0.332 ± 0.018 |
| w | −0.912 ± 0.150 (consistent with a constant) |
| α | (1.84 ± 0.04) × 10⁻¹⁰ per year |
| in plain units | one cubic metre gains **0.184 mm³ per year** |
| volume doubling time | 3.76 Gyr (distances double every 11.3 Gyr) |
| vacuum density | 6.79 × 10⁻²⁷ kg/m³ ≈ 4 hydrogen atoms per m³ |

The original guess of 0.01 m³ per m³ per year was too large by a factor of 54 million.

---

## Round 1 — expanding space inside galaxies

**Falsified, three independent ways.**

**1. Wrong direction.** A cosmological constant adds an *outward* acceleration
a = H_Λ²r. Outward push on a circular orbit means the orbit needs *less* gravity, so
outer stars orbit slightly *slower*, not faster.

**2. Wrong magnitude.** At the outermost measured radius of each galaxy, compare the
vacuum push with the extra pull the rotation curve actually needs:

| | |
|---|---|
| Median ratio (169 galaxies) | the push is **1 / 13,980** of the required pull |
| Effect of adding it to the fit | changes total χ² by 0.03% |

**3. Wrong shape.** Even flipping the push inward and fitting its strength freely per
galaxy (v² = v_bar² + (hr)²), it fails:

| Model | Free numbers | Median χ²ν | Best fit in N of 171 galaxies |
|---|---|---|---|
| Visible matter only | 0 | 85.7 | — |
| Vacuum push at the measured rate | 0 | 85.7 | — |
| Inward push, one shared strength | 1 | 30.3 | — |
| Inward push, free per galaxy | 171 | 16.3 | 15 |
| Acceleration law (RAR / MOND-like) | 1 | 10.1 | 7 |
| Pseudo-isothermal dark halo | 342 | 0.80 | 149 |

Three further killers for the push:

- **It needs absurd rates.** Median H = 9,487 km s⁻¹ Mpc⁻¹, middle 68% 5,281–14,217,
  against the cosmic 73.5. A galaxy's size would double in ~70 Myr.
- **It is not a property of space.** The required rate spans a factor of 23 across
  galaxies (1,381 to 32,070) and *falls with galaxy mass* (Spearman ρ = −0.54,
  p = 6 × 10⁻¹⁰, n = 115). Empty space cannot know a galaxy's mass.
- **The radial trend is backwards.** Fit the missing acceleration against radius over the
  outer half of each galaxy. The hypothesis needs slope **+1**. Measured: median **−0.64**,
  n = 140; only 8.6% of galaxies have any positive slope and 2.1% reach +0.5.

### What the galaxies demand instead

The discrepancy is not organised by distance, it is organised by acceleration:

| Scatter of the mass discrepancy | |
|---|---|
| against radius | 0.272 dex |
| against the visible matter's own acceleration | **0.147 dex** |

A single number describes all 164 galaxies (2,803 points): extra pull switches on below
**a₀ = (1.141 ± 0.074) × 10⁻¹⁰ m/s²**, intrinsic scatter 0.131 dex. This is the radial
acceleration relation (McGaugh, Lelli & Schombert 2016), recovered independently here.

**The one real hint in the original intuition:** a₀ coincides with cH₀/2π = 1.137 × 10⁻¹⁰
m/s². Galaxy outskirts and cosmic expansion *are* numerically linked. But anything built
from nearby space expanding at rate H gives accelerations of order H²r, and the observed
scale is cH₀ — a ratio of c/(H₀r) ≈ 136,000 at 30 kpc. A working mechanism must couple a
galaxy to the **cosmic horizon**, not to the empty space beside it. That is the road to
MOND and to Verlinde's emergent gravity, not to local void expansion.

Also: the baryonic Tully–Fisher relation holds with slope 3.75 ± 0.10 over 118 galaxies.
Mass ∝ speed⁴. A mechanism keyed to surrounding space has no way to know a galaxy's mass.

---

## Round 1b — does the local Hubble rate vary?

**Confirmed, but already explained by ordinary gravity.**

Per-galaxy Hubble rates from 663 local supernovae scatter widely nearby and converge
with distance, exactly as predicted:

| Distance (Mpc) | Scatter in each host's own H |
|---|---|
| 15–30 | ±21.6% |
| 30–50 | ±14.3% |
| 50–80 | ±9.3% |
| 80–130 | ±7.2% |
| 320–650 | ±6.1% (supernova distance errors dominate) |

One rule explains the funnel: each galaxy carries a random velocity of 334 ± 21 km/s on
top of the expansion, and a fixed velocity matters less the further away the galaxy is.

But gravity from the mapped galaxy distribution (2M++) already predicts the pattern:

| | |
|---|---|
| Observed vs gravity-predicted peculiar velocity, 450 hosts | slope **0.83 ± 0.18** |
| Coherent bulk flow before / after subtracting gravity | 260 ± 50 → 107 ± 37 km/s |
| Residual random motion | 334 ± 21 → 293 ± 25 km/s |

A slope of 1 means gravity accounts for the flows entirely. The measurement sits slightly
*below* 1, so there is no room for an extra void push — and every version of the
hypothesis tested later pushes the predicted slope *up*, away from the data.

---

## Round 2 — dark energy that follows local density

> The lower the density, the more space expands on its own; dark energy is zero inside
> galaxies and strongest in voids.

**The strong form is ruled out at ~8σ.**

Written as ρ_DE(x) = mean × (1 − κ·δ), with κ = 0 the standard model and κ ≈ 1 the strong
form. Lumpy dark energy gravitates (vacuum energy with −2× its density), so voids gain
extra push and dense regions extra pull, changing how structure grows:

| κ | Gravity strength today | Structure clumpiness | Predicted flow slope |
|---|---|---|---|
| 0 (standard) | 1.00 | 1.00 | 1.00 |
| 0.166 (max allowed) | 1.67 | 1.04 | 1.19 |
| 1.0 (your version) | 5.02 | 1.25 | **2.31** |

Measured slope: 0.83 ± 0.18. κ = 1 misses by about 8σ. The data allow **κ < 0.17**: dark
energy at most ~13% denser in a deep void than average.

Two further headwinds:
- **Theory.** Wang, Hui & Khoury (2012) prove that density-screened scalar fields of this
  type have a range of ~1 Mpc and cannot themselves drive cosmic acceleration. That range
  also rules out the 10–100 Mpc pushes this idea needs.
- **Direction of time.** As voids grow, total dark energy would rise; the DESI
  combinations lean the other way.

What the idea got right: switching dark energy off inside galaxies costs nothing
(it was already 14,000× too weak there), and bigger emptier voids do push harder —
though that is already true in plain gravity.

---

## Round 3 — a cap at zero density, and slower compounding

**The cap is fine. Square-root compounding is ruled out at ~27σ.**

Generalising the growth rule to dV/dt ∝ V^p makes the dark energy behave as a fluid with
**w = 1 − 2p**. Fitting Pantheon+ with p free:

| p | Meaning | Result |
|---|---|---|
| 0.96 | best fit (1σ: 0.90–1.02) | — |
| 1.00 | growth in proportion to volume | 0.6σ from best; fine |
| 2/3 | growth through a surface | ruled out at ~10σ |
| 0.50 | square root / diminishing returns | ruled out at ~**27σ** |

Square-root growth gives volume ∝ t², which is a matter-dominated universe: no
acceleration at all. Compounding must stay within a few per cent of proportional.

The **cap** (density cannot go below zero, so dark energy saturates in empty space) is
physically sensible but does not loosen the κ limit, because real voids sit at 10–50% of
mean density, well away from the cap. The strongest version the flows allow:

| Region | Dark energy relative to an average region |
|---|---|
| Completely empty (the cap) | 1.20× |
| Deep void (20% of mean) | 1.15× |
| Average | 1.00× |
| Wall (2× mean) | 0.86× |
| Cluster (200× mean) | 0.03× |
| Inside a galaxy | ≈ 0 |

**Diminishing returns with void size already happens** in plain gravity: edge outflow for
10 / 30 / 100 Mpc voids is roughly 106 / 200 / 266 km/s. Ten times the radius gives about
2.5 times the push, because bigger voids are shallower.

---

## Round 4 — the capped inverse law, over a real cosmic web

The first version with a proper treatment of *when* dark energy switches on. Density is
measured in units of today's mean, so the past matters: the universe was 8× denser at
z = 1 and 27× at z = 2. A law keyed to density is therefore also a law about time.

Shape: S(x) = 1/(1 + (x/x_half)ⁿ) × cutoff(x), averaged over a lognormal cosmic web
(2 Mpc cells, σ₈ = 0.81).

| Result | |
|---|---|
| The law as first specified (x_half = 30, n = 1, cut 200) | **+4.6**, mildly disfavoured |
| Plain inverse (n = 1) needs | x_half ≥ 46× mean |
| Sharp switch-off (large n) works from | x_half ≈ 10–15× mean |
| Largest void excess any allowed law permits | **~3%** over an average patch |
| What the data actually prefer | no density dependence, p ≈ 0.96, Δχ² = −4.8 |

The ~3% limit is far tighter than an earlier rough estimate of 15–20%, because that
estimate ignored the switch-on-over-time constraint. If the law is keyed to density
*relative to the mean at each epoch* instead, the limit relaxes to ~20%.

Galaxies are untouched: the push-to-gravity ratio at galaxy edges falls from ~10⁻⁴ under
uniform dark energy to effectively zero.

---

## Round 5 — head to head with the community, and timescape

Everything fitted on identical data and likelihoods.

| Model | Extra numbers | SN | +BAO | +flows |
|---|---|---|---|---|
| ΛCDM (= dV/dt = αV) | 0 | 0 | 0 | 0 |
| Compounding p free (= wCDM) | 1 | −0.3 | −4.6 | −4.8 |
| DESI evolving dark energy | 2 | −0.3 | −4.8 | −5.0 |
| Capped law as specified | 0 | +0.1 | +2.9 | +4.6 |
| Switch-on tuned to DESI's history | 3 | −0.2 | −4.6 | +4.1 |
| Same history, smooth in space | 3 | −0.2 | −4.6 | **−4.3** |

Supernovae alone cannot separate any of these. **BAO does the work.**

The tuned switch-on law reproduces DESI's preferred history to within 0.5% — but keyed to
*local* density it gives voids 44% more dark energy, pushing the predicted flow slope to
1.36 against the measured 0.83 ± 0.18 (~3σ miss). Keyed to the *average* density, the
same history keeps dark energy smooth in space and the flows come out right.

**Timescape** (Wiltshire), the published model closest to "voids expand faster", fitted
on the same supernovae:

| Treatment | Δχ² vs ΛCDM |
|---|---|
| Raw CMB-frame redshifts | 0.0 (ties) |
| With 2M++ peculiar-velocity corrections | +6.0 |
| Hubble-flow cut z > 0.075 | +2.3 |

The corrections assume standard gravity, which is exactly what timescape disputes, so
this does not settle it. Its own authors report strong evidence in its favour from a
different analysis. Genuinely open.

---

## Round 6 — surface growth, and the CMB

> Voids grow only through their surface, so growth is slower; and the surface itself
> grows as the void grows.

**Surface growth fails. Two other framings survive.**

Scores with the Planck 2018 compressed CMB priors included:

| Model | Extra numbers | SN+BAO+flows | +CMB |
|---|---|---|---|
| ΛCDM (volume compounding, p = 1) | 0 | 0 | 0 |
| Smooth surface (p = 2/3) | 0 | +239 | **+1050** |
| Rough surface, dimension D free | 1 | −4.8 (D = 2.87) | −0.2 (D = 2.98) |
| Void-surface law (Rice formula on the web) | 2 | +8.3 | **+23.8** |
| **Void-share law** | 3 | −4.9 | **−6.0** |
| **Switch-on by average density** | 3 | −5.6 | **−6.2** |
| DESI evolving dark energy (benchmark) | 2 | −5.0 | −5.3 |
| Capped law, local density | 0 | +4.6 | +7.5 |

- **Why surface growth fails.** Surface goes as size², volume as size³, so growth per
  cubic metre falls as 1/size. That is expansion that decelerates forever.
- **The effective growth dimension is D = 2.982, 1σ 2.96–3.00.** Growth must fill the
  volume. A genuine surface (D = 2) is excluded overwhelmingly.
- **Why the void-surface law fails.** Truly empty voids form late, so dark energy would
  switch on too late and too steeply: 17% of today's value at z = 2 where the data want
  70–80%.

### The two survivors

**Switch-on by average density** — the capped inverse law applied to the universe as a
whole: x_half ≈ 24× today's mean density (z ≈ 1.9), n ≈ 1.1 (a plain inverse, as
originally proposed), compounding p ≈ 0.94.

**Void-share law** — dark energy grows with the fraction of the universe below a quarter
of mean density (7% at z = 2, 18% at z = 1, 34% today), to the power 0.59, p ≈ 0.91.

Both give the same history: dark energy rising slightly to z ≈ 0.5, then falling ~6%.
Both keep dark energy smooth in space, so the flow slope stays at 1.00.

| Significance over ΛCDM on this data combination | |
|---|---|
| Switch-on by average density (3 params) | 1.63σ |
| Void-share law (3 params) | 1.59σ |
| DESI evolving dark energy (2 params) | 1.80σ |

Competitive, not winning. DESI's own analysis reaches 2.8σ with these supernovae and
3.8–4.2σ with DES-5yr or Union3, which is the obvious next test.

**What is still missing:** a physical reason for dark energy to track the universe's
*average* density or void share rather than local density. As written these describe the
history, they do not explain it. Buchert averaging and timescape are the nearest existing
frameworks.

---

## The separate galaxy run (`galaxy_run/`)

A cross-validated rotation-curve study with 30 galaxies held out, scored as RMS error in
log₁₀ speed (dex). Lower is better.

| Model | All radii | Outer radii | Stellar mass scale |
|---|---|---|---|
| Acceleration law (RAR) | 0.076 | **0.067** | 1.06 |
| Shared NFW halo | 0.131 | 0.133 | 1.39 |
| Radius gradient (control) | 0.136 | 0.138 | 1.71 |
| **Density gradient** | 0.186 | **0.126** | 1.96 |
| Baryons only | 0.201 | 0.220 | 1.97 |
| Screened Λ | 0.201 | 0.220 | 1.97 |

Paired bootstrap at outer radii, 103 galaxies:

- beats baryons-only by 0.094 dex, in **100%** of resamples (95 of 103 galaxies)
- loses to the acceleration law by 0.059 dex, in **100%** of resamples (wins 21 of 103)
- statistically tied with the radius-gradient control (−0.012 dex, p = 0.83) — **the
  density dependence adds nothing over simply using radius**
- needs stars ~2× heavier than standard (1.96 vs 1.06 for the acceleration law)

Screened Λ scores *identically* to baryons alone, confirming from the galaxy side that
switching dark energy off inside galaxies changes nothing.

Known issues in that run: the field-speed cap (1,500 km/s) was hit in 3 of 5 folds, and
one fold produced an invalid orbit, so "all radii" comparisons are flagged ineligible.

**The 30 reserved test galaxies are still sealed.** Do not touch them until a law is
frozen in phase 3.

---

## Summary

### Correct
- dV/dt = αV is exactly a cosmological constant; α = 1.84 × 10⁻¹⁰ /yr
- the local Hubble rate really does vary, and converges with distance
- switching dark energy off in galaxies costs nothing
- bigger, emptier voids really do push harder
- a cap at zero density is physically sensible
- diminishing returns with void size is real
- a₀ ≈ cH₀/2π: galaxy outskirts and cosmic expansion *are* linked

### Wrong
- expanding space does not cause flat rotation curves (direction, magnitude, shape)
- void push does not improve the Hubble-rate variation — gravity already explains it
- strong local-density dependence is excluded at ~8σ
- square-root compounding at ~27σ; surface growth at Δχ² ≈ +1050
- density-based force laws for rotation curves add nothing over a radius control

### Still open
- switch-on by the universe's average density, or by its void share (both ~1.6σ)
- timescape, where the published literature genuinely disagrees
- the a₀ = cH₀/2π coincidence, which needs a horizon-scale mechanism

---

## Corrections log

Errors found in earlier write-ups, recorded rather than quietly fixed.

| Date | What was wrong | Correction |
|---|---|---|
| 2026-11-05 | The HTML report quoted "+0.3 ± 0.7%" for the expansion step inside vs outside 300 Mpc. That is the figure from *uncorrected* redshifts, quoted inside a paragraph otherwise using gravity-corrected ones. | Corrected figures: raw **+0.33 ± 0.67%**, gravity-corrected **+1.14 ± 0.66%**. Neither detects a large local void; the corrected one is a 1.7σ hint, so "no step" overstated the case. |
| 2026-10-26 | A rough estimate put the allowed void excess at 15–20%. | It ignored the switch-on-over-time constraint. With the cosmic web the limit is **~3%** (local density) or **~20%** (density relative to each epoch's mean). |
| 2026-10-26 | A growth-equation term went 0/0 and overflowed when a law switched dark energy fully off in the early universe, producing spurious +28 scores for sharp laws. | Fixed by computing κ·f directly instead of κ and f separately. |
