# Data

No datasets are committed to this repository. The notebook downloads what it needs into
`data/` (gitignored, ~35 MB) and verifies each file against the SHA-256 below. A
"fingerprint differs" message means the provider has updated the file since this work was
done, not necessarily that anything is wrong — but any number in FINDINGS.md was computed
against *these* bytes.

## Supernovae — Pantheon+SH0ES

1,657 light curves of 1,508 Type Ia supernovae with the full statistical + systematic
covariance and SH0ES Cepheid host distances.

| File | SHA-256 |
|---|---|
| `Pantheon+SH0ES.dat` | `1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8` |
| `Pantheon+SH0ES_STAT+SYS.cov` | `abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc` |

Source: <https://github.com/PantheonPlusSH0ES/DataRelease>
Cite: Scolnic et al. (2022), ApJ 938, 113; Brout et al. (2022), ApJ 938, 110.
The `VPEC` column carries peculiar velocities predicted from the 2M++ density field
(Carrick et al. 2015, MNRAS 450, 317), which this work uses for the galaxy-flow test.

## Baryon acoustic oscillations — DESI DR2

13 distance measurements (D_M/r_d, D_H/r_d, D_V/r_d) as a Gaussian summary.

| File | SHA-256 |
|---|---|
| `desi_gaussian_bao_ALL_GCcomb_mean.txt` | `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` |
| `desi_gaussian_bao_ALL_GCcomb_cov.txt` | `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509` |

Source: <https://github.com/CobayaSampler/bao_data> (`desi_bao_dr2/`)
Cite: DESI Collaboration (2025), arXiv:2503.14738.

## Rotation curves — SPARC

175 galaxies, 3,391 velocity points, with gas / disk / bulge mass models at 3.6 µm.

| File | SHA-256 |
|---|---|
| `Rotmod_LTG.zip` | `0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588` |

Source: <http://astroweb.case.edu/SPARC/> (Zenodo mirror:
<https://zenodo.org/records/16284118>)
Cite: Lelli, McGaugh & Schombert (2016), AJ 152, 157.
The BTFR table used for the mass–speed relation is from Lelli, McGaugh & Schombert
(2016), ApJL 816, L14.

Note: during the original sessions the Case Western host was unreachable from the
analysis machine and the files came from a public GitHub mirror of the same release;
distances and velocities were checked against the published database.

## Cosmic microwave background — Planck 2018 compressed priors

Three numbers plus covariance: the shift parameter R, the acoustic scale ℓ_A, and the
baryon density ω_b. Hard-coded in the notebook, not downloaded.

```
R      = 1.74963
l_A    = 301.80845
omega_b= 0.02237
```

Cite: Planck Collaboration (2020), A&A 641, A6; compressed form from Zhai & Wang (2019),
JCAP 07, 005.

**Caveat.** The decoupling redshift uses the Hu & Sugiyama (1996) fitting formula, which
lands at z\* = 1091.9 against Planck's full calculation of 1089.9 — enough on its own to
shift ℓ_A by 5σ. The notebook calibrates this out once, at Planck's best-fit ΛCDM, so
model-to-model differences come from the distance integral, which is exact. This is
reliable for late-time dark energy only. Any model that changes the early universe needs
the full Planck likelihood (see H07 in HYPOTHESES.md).

## The separate galaxy run

`galaxy_run/` holds the outputs of a cross-validated rotation-curve study run separately
by the author. Its own `data_fingerprint.json`, `software_fingerprint.json` and
`data_audit.json` record its provenance. Its 30 held-out test galaxies have not been
looked at.
