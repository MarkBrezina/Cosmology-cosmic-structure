# Pre-registered hypotheses for phase 3

Phase 1 and 2 were exploratory: ideas were proposed and tested in the same breath, which
is how you get a result that looks significant because you went looking for it. Phase 3
fixes that by writing each test down **before** running it.

Format follows the same discipline used in the author's market-research lab: statement,
null, test, decision rule, and **both** implications written out in advance, so that
neither outcome can be re-interpreted after the fact.

## Rules

1. The decision rule is frozen before the data are touched. No rule changes afterwards.
2. Both implications are written before the test. If you cannot say what a rejection
   would mean, the hypothesis is not ready.
3. A predictive claim must survive out of sample, not only in sample.
4. Parameter counts are declared up front; Δχ² is judged against them.
5. The 30 reserved SPARC galaxies stay sealed until a single law is frozen. They get used
   **once**.
6. Results are recorded whichever way they come out, including into FINDINGS.md.

---

## H01 — The void-share law survives stronger supernova samples

**Statement.** Dark energy scaling with the universe's void fraction fits DES-5yr and
Union3 at least as well as DESI's w₀wₐ model does, on the same BAO and CMB.

**Null.** It fits worse than w₀wₐ by more than Δχ² = 2 per extra parameter.

**Test.** Refit both on DESI DR2 BAO + Planck priors + each supernova sample in turn.

**Decision rule.** Survives if, for both samples, its Δχ² vs ΛCDM is within 2 of w₀wₐ's
while using no more than one extra parameter.

**If confirmed.** The law is a genuine competitor, not an artefact of Pantheon+, and
deserves a proper MCMC and a write-up.
**If rejected.** The ~1.6σ on Pantheon+ was sample-specific. Drop it.

---

## H02 — Average-density switch-on is distinguishable from w₀wₐ

**Statement.** The switch-on law and DESI's w₀wₐ predict expansion histories that differ
by more than current data can resolve at some redshift, i.e. there exists a measurement
that separates them.

**Null.** The two histories agree to within observational error everywhere 0 < z < 2.5,
making the switch-on law a reparametrisation with no independent content.

**Test.** Compute the maximum fractional difference in ρ_DE(z) between the two best fits,
and compare it with the DESI DR2 distance errors at each redshift.

**Decision rule.** Distinguishable if the difference exceeds 1σ of the data at any z.

**If confirmed.** Name the redshift and the observable; that is the discriminating test.
**If rejected.** The law adds a story, not physics. Say so plainly and stop fitting it.

---

## H03 — Void galaxies rotate differently from wall galaxies

The most direct test of the original "voids push on galaxies" intuition.

**Statement.** At fixed baryonic mass, rotation curves of galaxies in voids differ
systematically from those of galaxies in walls and filaments.

**Null.** No difference beyond the scatter of the radial acceleration relation.

**Test.** Cross-match SPARC (and, if coverage is thin, SPARC + a larger HI survey) with a
public void catalogue (SDSS / DESI). Compare RAR residuals, matched in baryonic mass and
surface brightness. Permutation test on the environment label.

**Decision rule.** Confirmed if the mean RAR residual differs by more than 0.05 dex with
p < 0.01 under 10,000 permutations.

**If confirmed.** Environment affects galaxy dynamics — a major result, and the first
real support for a void-push effect.
**If rejected.** The void-push idea is dead at galaxy scale, independently of the flow
measurement. Record it as closed.

---

## H04 — No push survives free stellar mass-to-light ratios

**Statement.** With stellar M/L free per galaxy (with a sensible prior), the
density-gradient and expansion-push models still lose to the acceleration law.

**Null.** One of them matches the acceleration law once M/L is free.

**Test.** Refit the SPARC development set with M/L free, log-normal prior centred at 0.5
(disk) / 0.7 (bulge), width 0.1 dex. Same cross-validation folds.

**Decision rule.** The null survives only if a push model comes within 0.02 dex of the
acceleration law at outer radii in at least 95% of paired bootstrap resamples.

**If confirmed.** The galaxy-side push is closed for good; the mass scale of ~1.96 was
the model absorbing stellar mass, not finding new physics.
**If rejected.** Reopen the galaxy side with M/L as a first-class parameter.

---

## H05 — The surviving laws match the growth of structure

**Statement.** The two surviving laws predict fσ₈(z) consistent with redshift-space
distortion measurements over 0 < z < 1.5.

**Null.** They deviate by more than 2σ from the compiled fσ₈ data.

**Test.** Compute fσ₈(z) for each best fit and compare with a published RSD compilation
(BOSS, eBOSS, DESI), using its covariance.

**Decision rule.** Pass if χ² over the compilation is within 4 of the ΛCDM value.

**If confirmed.** The laws are consistent with structure growth as well as distances,
which is a real hurdle cleared.
**If rejected.** They are distance-only fits; say so and stop calling them competitive.

---

## H06 — a₀ = cH₀/2π is a coincidence

**Statement.** The agreement between the fitted a₀ and cH₀/2π is numerological: it does
not hold when H₀ is varied within its uncertainty, or across independent a₀ measurements.

**Null.** The relation holds to better than the combined measurement precision.

**Test.** Propagate the H₀ uncertainty (73.5 ± 1.0 and 67.4 ± 0.5 separately) and compare
with a₀ from SPARC, from the BTFR normalisation, and from published dwarf-galaxy work.

**Decision rule.** Coincidence if the two differ by more than 2σ for either H₀ value.

**If confirmed.** Stop treating the relation as a clue.
**If rejected.** The link is tight enough to justify studying horizon-coupled mechanisms
(emergent gravity, MOND's cosmological-coincidence problem) seriously.

---

## H07 — Full CMB likelihood does not change the verdict

**Statement.** Using the full Planck likelihood instead of three compressed distance
priors leaves the ranking of the surviving laws unchanged.

**Null.** The ranking changes.

**Test.** Re-run the two survivors and w₀wₐ through `cobaya` + `CLASS`/`CAMB` with the
Planck 2018 likelihoods.

**Decision rule.** Unchanged if every Δχ² shifts by less than 2 and no ordering flips.

**If confirmed.** The compressed-prior results in this repo stand.
**If rejected.** Every CMB-including number in FINDINGS.md is retracted and recomputed.

---

## H08 — Timescape, tested without assuming standard gravity

**Statement.** Timescape's standing against ΛCDM depends on whether peculiar-velocity
corrections derived under standard gravity are applied.

**Null.** It loses to ΛCDM under any self-consistent treatment.

**Test.** Fit both models to the same supernovae three ways: raw CMB-frame redshifts, a
Hubble-flow cut at z > 0.075, and with velocities fitted jointly rather than imported.

**Decision rule.** The dependence is real if Δχ² between treatments exceeds 5.

**If confirmed.** The published disagreement is a methods disagreement, which is worth
writing up on its own.
**If rejected.** Timescape's standing is robust either way; report that too.

---

## Candidate hypotheses not yet ready

These need a sharper statement before they can be registered:

- Whether the void-share exponent (β ≈ 0.59) can be derived rather than fitted.
- Whether Buchert averaging produces the average-density switch-on from first principles.
- Whether any horizon-coupled mechanism reproduces a₀ *and* the Tully–Fisher slope of 3.75.
