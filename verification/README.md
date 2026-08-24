# Verification harness

Two layers, kept separate on purpose.

`tests/test_backend.py` is the regression suite: 132 unit tests, a couple of
minutes, run it on every change.

This directory checks the claims rather than guarding against regressions.
Each script verifies one section of the feature catalogue, or one scientific
claim, end to end against real deposited structures. It takes about seven
minutes.

```bash
python verification/run_all.py             # everything
python verification/run_all.py features    # catalogue features 1-121
python verification/run_all.py science     # physics and chemistry only
```

Running an individual script standalone (rather than through `run_all.py`,
which already passes `-u`) or in the background, use `python -u`:

```bash
python -u verification/science_02.py > run.log 2>&1 &
```

Python block-buffers stdout when it is not a terminal. If the process is
killed externally - a timeout, closing the session it was backgrounded
in - before it exits normally, that buffer is never flushed and every
`print()` is lost; only unbuffered logger output survives. `-u` writes each
line as it is produced, so a killed run still leaves a usable partial log.

## Feature scripts

| Script | Features | Module |
|---|---|---|
| `verify_1_11.py` | 1–11 | `io.py`, `cleaner.py` |
| `verify_12_22.py` | 12–22 | `protonator.py` |
| `verify_23_28.py` | 23–28 | `analyzer.py` |
| `verify_29_48.py` | 29–48 | `minimizer.py` |
| `verify_49_69.py` | 49–69 | `site_analyzer.py` |
| `verify_70_80.py` | 70–80 | `exporter.py`, `reporter.py` |
| `verify_81_110.py` | 81–110 | `app.py` |
| `verify_111_121.py` | 111–121 | `cli.py`, frontend |

## Science scripts

| Script | Asks |
|---|---|
| `science_01.py` | Are disulfides left oxidised? Does pH change protonation? Does tier 2 collapse the pocket? |
| `science_02.py` | On real complexes: pocket collapse, rotation invariance, whether "structural" waters are hydrogen bonded |
| `science_03.py` | pKa limits, histidine tautomers, PDBQT charges, whether `converged` means anything |

## Structures

`pdbs/` holds three deposited structures, used because 1CRN is small, capped,
ligand-free and pocket-free and hid every science defect found so far:

- **1STP** streptavidin + biotin — a deeply buried ligand, and an uncapped C-terminus
- **1HSG** HIV-1 protease + MK1 — catalytic aspartate dyad, two chains
- **4INS** insulin + zinc — a metal site

## Three checks deliberately contradict the catalogue

The code is right and the catalogue text is out of date. Each is annotated in
place with the measurement that justified it.

| Feature | Catalogue says | Code does | Why |
|---|---|---|---|
| 48 | `HBonds` constraints | unconstrained | Constraints are for an MD timestep; nothing here runs dynamics. They stop hydrogens relaxing and make the residual force incomparable to the tolerance. 1CRN: −5162.7 kJ/mol at RMS 71.16 constrained, −5166.9 at 5.76 unconstrained. |
| 52 | six axis-aligned rays | 14 rays over 7 axes | Six rays made the answer depend on input orientation — rotating streptavidin 45° changed cavity volume by 47.8%. Now 8.9%. |
| 70 | `-xh` merges non-polar hydrogens | `-xh` not passed | In OpenBabel it *preserves* hydrogens. With it, all 842 hydrogens were typed `HD` — the AutoDock type for a hydrogen on N or O — so carbon-bound hydrogens were presented as hydrogen-bond donors. Without it, 213 genuinely polar ones. |

## Known failing check

`verify_111_121.py` reports one failure: `main.js` still reads the pre-rewrite
report shape (`iterations`, the string `'N/A'`). That is expected — the
frontend is being rewritten, and the current contract is documented at the end
of `BACKEND_AUDIT.md`.
