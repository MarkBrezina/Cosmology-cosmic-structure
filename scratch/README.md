# scratch/

Working scripts, kept so phase 1 loses nothing. **Do not build on these.**

Several bootstrap themselves by `exec`-ing string slices of each other's source
(`lab3.py` and `lab4.py` read `lab2.py`; `smooth_check.py` reads `followups.py`). That
worked inside a single live session and is not a sane foundation. They also carry absolute
paths from that machine (`/home/claude/lab/...`), so they will not run as-is without editing. They are here as a
record of how the numbers in FINDINGS.md were produced, not as a library.

| File | What it did |
|---|---|
| `lab_core.py` | first engine: data loading, the density law, the cosmic web, SN/BAO likelihoods, growth |
| `lab_analysis.py` | first end-to-end run of that engine |
| `lab_test2.py` | law scan and optimiser timing |
| `lab_debug.py` | traced the 0/0 overflow in the growth term (see the corrections log) |
| `head_to_head.py` | community models on identical data; the timescape distance calculation |
| `followups.py` | timescape with uncorrected redshifts; fitting a law to DESI's history |
| `smooth_check.py` | the same history with dark energy kept smooth in space |
| `lab2.py` | round 6: CMB distance priors, surface-growth laws, the scoreboards |
| `lab3.py` | strict average-density switch-on; the growth-dimension profile |
| `lab4.py` | the void-share law |

The maintained versions of all of this are `analysis/build_notebook.py` (which generates
the notebook) and the notebook itself.
