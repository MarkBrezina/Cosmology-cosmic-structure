# Dark energy and empty space — a hypothesis, tested

An independent test of one idea: **that dark energy depends on how empty space is.**

The idea started as a proposal that empty space expands in proportion to how much empty
space is already there, and that this could explain both galaxy rotation curves and the
acceleration of the universe. Over several rounds it was written down as equations and
run against public data: 175 galaxy rotation curves (SPARC), 1,657 Type Ia supernovae
(Pantheon+SH0ES), 13 baryon acoustic oscillation distances (DESI DR2), the Planck 2018
compressed CMB distance priors, and the peculiar motions of 450 supernova host galaxies.

Most versions of the idea were falsified. Two survive and are competitive with the
current community benchmark. **[FINDINGS.md](FINDINGS.md) is the record of every version
tested, what it predicted, and what the data said.**

## Status

This repository is at **phase 1 of 3**, and is a notes dump by design.

| Phase | What it is | State |
|---|---|---|
| 1 | Upload everything, unstructured, so nothing is lost | you are here |
| 2 | Structure the notes; separate conclusions from scratch work | not started |
| 3 | Rebuild and re-run against pre-registered hypotheses | [HYPOTHESES.md](HYPOTHESES.md) is the seed |

Phase 1 is deliberately messy. It is a record of ideas being run into a wall repeatedly
to see which ones survived. Several did not, and the failures are written down with the
same care as the survivors, because that is the part of this work with the clearest value.

Nothing here is peer reviewed, and nothing here is a claim to have overturned anything.

## What survived, in one table

Scores are Δχ² against the standard model (ΛCDM) on supernovae + BAO + CMB + galaxy
flows. Negative is better. Roughly, Δχ² = −4 is a 2σ improvement.

| Idea | Δχ² | Verdict |
|---|---|---|
| Empty space compounds in proportion to its volume (dV/dt = αV) | 0 | Correct, but this *is* the cosmological constant |
| Dark energy switches on as the **universe's average** density falls | −6.2 | Survives; matches DESI's evolving dark energy |
| Dark energy grows with the **void share** of the universe | −6.0 | Survives; matches DESI's evolving dark energy |
| DESI's evolving dark energy (community benchmark) | −5.3 | For comparison |
| Dark energy keyed to **local** density, as first specified | +7.5 | Disfavoured |
| Empty space grows through its **surface** (dV/dt ∝ V^⅔) | +1050 | Ruled out |
| Expanding space causes flat galaxy rotation curves | — | Ruled out, three independent ways |

The two survivors are at 1.6σ over ΛCDM on this data combination, against 1.8σ for the
community benchmark fitted the same way. That is *competitive*, not *winning*, and the
preference is driven by BAO, not by supernovae.

## What is in here

```
FINDINGS.md          every version tested, with numbers and verdicts   <- start here
HYPOTHESES.md        pre-registered hypotheses for phase 3
DATA.md              data sources, fingerprints, licences
notebooks/           dark_energy_law_lab.ipynb — the main workbench, pre-run
analysis/            the round-1 pipeline and the notebook generator (both standalone)
report/              a standalone HTML write-up of round 1
results/             fitted outputs, scoreboards, figures
galaxy_run/          a separate cross-validated rotation-curve run (30 galaxies held out)
scratch/             working scripts, kept for the record; see the warning below
```

## Running it

```bash
pip install -r requirements.txt
jupyter lab notebooks/dark_energy_law_lab.ipynb
```

Run top to bottom. The first run downloads ~35 MB into `data/` (gitignored) and checks
each file against a SHA-256 fingerprint. A full run takes about 7 minutes.

The notebook is shipped **already executed**, so every figure and number can be read
without running anything.

To change the hypothesis, edit `MY_LAW` in section 2, or write a new dark-energy history
in section 15. Section 9 scans the space of density laws; section 14 compares the
surface-growth laws with and without the CMB.

### Checks that the pipeline is sound

The notebook reproduces published results before it is trusted with anything new:

| Quantity | This repo | Published |
|---|---|---|
| H₀ from Pantheon+SH0ES | 73.53 ± 1.02 | 73.5 ± 1.1 |
| Ωm from DESI DR2 BAO alone | 0.297 | 0.2975 ± 0.0086 |
| ΛCDM from the CMB priors alone | Ωm 0.3153, h 0.6736, ωb 0.02237 | Planck 2018: identical |
| Constant w (SN+BAO+CMB) | −0.988 ± 0.023 | −0.986 ± 0.026 |
| DESI evolving dark energy | w₀ = −0.87, wₐ = −0.47 | w₀ = −0.86, wₐ = −0.58 |
| a₀ from SPARC rotation curves | 1.141 ± 0.074 ×10⁻¹⁰ m/s² | 1.20 ± 0.26 ×10⁻¹⁰ m/s² |

### A warning about `scratch/`

Those scripts produced real numbers in the record, but several of them bootstrap
themselves by `exec`-ing slices of each other's source. They are kept because phase 1 is
about losing nothing, not about being tidy. Do not build on them; rebuild in phase 3.

## How this was produced

The analysis and code were developed in collaboration with Claude (Anthropic) across
several sessions, working from public datasets. The ideas under test are the author's;
the implementation, the literature comparisons and the falsification tests were done
jointly. Anyone checking this work should re-run the notebook rather than trust the prose.

## Licence

Code: MIT. Text and figures: CC BY 4.0. The datasets keep their own terms; see
[DATA.md](DATA.md).
