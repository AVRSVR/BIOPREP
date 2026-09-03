# BioPrep

Prepares protein structures for molecular docking and MD: clean, protonate,
repair, energy-minimise, and export to docking-ready formats — through a web UI,
a CLI, or a batch pipeline.

Everything runs locally. Structures are not uploaded anywhere.

---

## What it does

| Step | Detail |
|---|---|
| **Load** | PDB or mmCIF, uploaded or fetched live by PDB ID from RCSB |
| **Clean** | Chain selection, water removal (with optional structural-water retention), heteroatom removal with per-ligand protection |
| **Repair** | Missing terminal atoms, missing heavy atoms/sidechains, missing loops rebuilt from SEQRES |
| **Protonate** | pH-dependent hydrogen placement via PDBFixer, optionally using PROPKA for structure-specific pKa instead of model-compound defaults |
| **Minimise** | AMBER14 or CHARMM36, optional GBSA implicit solvent, GPU-accelerated when available, with a graded fallback when a residue can't be parameterised |
| **Analyse** | Grid-based binding-pocket detection (LIGSITE-style, seven-axis scan) with pharmacophore points and a druggability breakdown |
| **Export** | PDB, AutoDock/Vina PDBQT (rigid receptor, correct polar-hydrogen typing), GROMACS-named PDB |
| **Scale** | Batch mode for loose files or a zip archive (up to 500 structures), with a combined report |

---

## Quick start

```bash
python -m venv .venv
.venv/Scripts/pip install -r bioprep/requirements.txt
.venv/Scripts/python bioprep/app.py
```

Then open <http://localhost:5000>.

**CLI:**

```bash
bioprep --input structure.pdb --output prepared.pdb --chain A --ph 7.4 --minimize
```

**Docker:**

```bash
docker build -t bioprep .
docker run -p 8000:8000 bioprep
```

`obabel` (OpenBabel) is needed for PDBQT export; everything else installs from
pip. The Docker image installs it for you.

### Configuration

| Variable | Purpose |
|---|---|
| `BIOPREP_DATA_DIR` | Where job history, saved templates and processed structures are written |
| `BIOPREP_PUBLIC_DEMO` | Disables History and Templates. These are shared, unscoped state — fine for one person locally, a cross-visitor leak on a shared deployment |
| `BIOPREP_MINIMIZER_PLATFORM` | Pins OpenMM to one platform (`CPU`, `OpenCL`, `CUDA`) instead of auto-selecting |

---

## Engineering notes

The interesting part of this project isn't the feature list — it's what broke
when it met real deposited structures instead of a tidy test case.
[`BACKEND_AUDIT.md`](BACKEND_AUDIT.md) documents each finding with the evidence
that produced it. A few examples:

- **Chain selection silently defeated ligand protection.** Biopython's
  `Select.accept_chain` is an all-or-nothing gate — rejecting a chain means
  `accept_residue` never runs for anything in it. A ligand routinely sits on a
  different chain than its receptor (HIV-1 protease's inhibitor is chain B of an
  A/B dimer), so `chains=['A']` dropped an explicitly protected ligand and the
  report just said `ligands_preserved: []`. Every isolated test of either
  feature passed; only exercising them together exposed it.

- **A chain end missing backbone atoms crashed protonation.** Real structures
  routinely have no resolved density past a terminal residue's amide nitrogen.
  PDBFixer's terminal-placement code reads that residue's `O` and `CA`
  unconditionally, so discarding them raised `KeyError: 'O'`. Found by sweeping
  89 structures pulled live from RCSB — zinc fingers, protein–DNA/RNA
  complexes, heme/Fe–S/Cu proteins, membrane proteins, NMR ensembles — of which
  4 hit it.

- **`converged` measured the wrong thing.** It reported whether energy had
  fallen, so a 100-iteration run "converged" while 1000 iterations went 900
  kJ/mol lower. It now compares RMS residual force against the tolerance the
  minimiser was actually given.

- **Pocket detection depended on molecular orientation.** Rotating
  streptavidin 45° changed total cavity volume by 47.8% for a physically
  identical cavity. Scanning all seven LIGSITE axes cut worst-case drift to
  8.9%.

- **GPU acceleration, verified before trusting.** Identical starting energy,
  final energies within 2 kJ/mol, and mean atom-position deviation of 0.011 Å
  against a CPU run — for a 3.8× speedup on a small system, and a real job
  (cryo-EM structure, 6,762 atoms, everything enabled) going from *not finished
  after 70 minutes* to **96 seconds**.

Limitations that were found and deliberately **not** papered over — tyrosine is
never deprotonated (no anionic template exists in either force field offered),
the druggability score is an unvalidated heuristic that mis-ranks a known site,
and ligands are never parameterised — are documented with the reasoning in the
audit rather than left for a user to discover.

---

## Testing

```bash
cd tests && python -m unittest test_backend test_webapp
```

152 tests. Each regression test names the finding it guards, and reproduces it
from a synthetic fixture wherever possible rather than depending on one
particular PDB file.

---

## Stack

Flask · Biopython · PDBFixer · OpenMM · PROPKA · OpenBabel · NumPy/SciPy ·
scikit-learn (DBSCAN) · 3Dmol.js · vanilla JS, no framework

## License

MIT — see [LICENSE](LICENSE).
