# BioPrep Backend Audit

Audit date: 2026-08-18. Scope: Python backend only (`bioprep/app.py`, `bioprep/bioprep/**`).
Frontend (`static/`, `templates/`) intentionally excluded — it is being rewritten.

Environment used for verification: `.venv` (Python 3.13, biopython 1.86, OpenMM 8.4.0,
pdbfixer 1.12.0, numpy 2.4.3, scipy 1.17.1, scikit-learn 1.8.0).

Status legend: **CONFIRMED** = reproduced by running code. **INSPECTION** = read from source,
not yet executed.

---

## Progress

### Round 1 — core scientific layer (done, verified)

| ID | Fix | Verification |
|---|---|---|
| A9, A10, A11 | `protonator.py` rewritten. Ligand HETATM records are split out before PDBFixer runs and spliced back with original coordinates, renumbered serials and remapped `CONECT` records. Returns a status dict instead of `None`. | BTN survives protonation with byte-identical coordinates |
| A21, A22, A24, A25, A30, B1 | `minimizer.py` rewritten with the three-tier fallback. Tier 2 strips non-parameterisable residues via `Modeller`, minimises the rest, and merges coordinates back by `(chain, resname, resid, atomname)`. | Ligand structure: `status=partial`, 978.1 → −5155.8 kJ/mol, protein moved, BTN frozen |
| B2 | GBSA force field now also loads `amber14/tip3pfb.xml`, so retained waters no longer abort minimisation. | Water + GBSA: `status=full`, minimised (was: hard failure) |
| B3 | `minimize_structure` always returns `status` (`full`/`partial`/`partial_no_implicit_solvent`/`failed`). A copied-through file is never reported as minimised. | Failure path returns `status=failed` with populated `error` |
| B4 | Starting energies above 1e12 kJ/mol attach a clash warning rather than being silently accepted. | 1.045e13 start now emits the warning |
| B5 | `protect_ligands` is now honoured inside `BioPrepSelect` and outranks `'ALL'`. | Protected BTN survives `remove ALL`; unprotected SO4 removed |
| B6 | Structural waters keyed by `(chain_id, residue.id)`. | Chain-B water 500 Å away no longer preserved |
| B7 | `add_hydrogens` reports `hydrogens_added` truthfully plus a `warnings` list. | Returns `{'hydrogens_added': True, 'ligands_preserved': [...]}` |
| B16 | Energies are `float` or `None`, never the string `'N/A'`. | `json.dumps(..., allow_nan=False)` succeeds |
| B28 | Ligand `CONECT` records preserved and remapped through the pipeline. | — |
| B9, B34 | `cli.py` rewritten: correct kwarg, repeatable `--chain`, plus `--remove-het`, `--protect`, `--minimize`, `--force-field`, `--no-gbsa`. | Full CLI run succeeds end to end |
| — | New `core/residues.py` holds the shared classification tables (amino acids, nucleic acids, waters, ions). | — |

### Round 2 — API layer (done, verified)

New `core/pipeline.py` holds the single clean → protonate → minimise → export
sequence. Precision, batch and high-throughput all call it, which is what stops
the modes drifting apart. `app.py` is now transport only.

| ID | Fix | Verification |
|---|---|---|
| B8 | Failed docking export is recorded as `docking_export.succeeded = false` with the error, and no longer sets `docking_target`. | Report renders an export section with the failure |
| B15 | Each request gets its own `tempfile.mkdtemp` working directory; no fixed filenames in the shared temp dir. | Concurrent requests cannot collide |
| B17, B18 | Working directories are removed in `finally`, and archives are read into memory before the directory is deleted so nothing is left holding a file handle. A startup sweep clears anything orphaned by a crash. | 0 `bioprep_*` directories left after the full API suite |
| B19 | `SessionStore` is a bounded LRU (200 entries) that deletes the stored file on eviction. | — |
| B20 | High-throughput uses the shared pipeline, so chain selection, structural waters and docking export apply there too. | HT run honours `protect_ligands`; settings echoed into the log |
| B21 | `heteroatoms.removed`/`retained` computed from what the selector actually did, including the `ALL` case. | `removed: ['SO4']`, `retained: ['BTN']` |
| B22 | Waters retained by structural-water detection are subtracted from the removed count and reported separately. | — |
| B23 | `atoms_before`/`after` documented as whole-structure counts; `delta` guaranteed `int`. | — |
| B24, B25 | `report_to_text` uses defensive lookups throughout. | Partial and empty reports render instead of raising |
| B29 | `JsonStore` guards history and templates with a lock and writes atomically via `os.replace`. | 25 concurrent template writes all landed; file still valid JSON |
| B30 | The structure is base64-encoded once and reused when the download and viewer files are the same. | — |
| B31 | `debug=False`, bound to `127.0.0.1`. | — |
| B32 | Exceptions are logged server-side; clients get a generic message. | — |
| B33 | ZIP extraction capped at 500 members and 2 GB expanded. | — |
| B32 (path) | History and site routes use Flask's `uuid` converter, so a non-UUID id cannot reach the filesystem. | `..%2f..%2fapp` → 404; valid UUID → 200 |
| B35, B36 | Dead code removed: `allowed_file`, `ALLOWED_EXTENSIONS`, `batch_progress`, and the unused `queue`/`threading`/`Response`/`stream_with_context` imports. | — |
| B40 | `numpy` added to `install_requires`; `python_requires` declared. | — |
| B41, B43 | Stale UTF-16 captures and the ad-hoc scripts deleted; coverage moved into `tests/`. | — |

Regression suite added at `tests/test_backend.py` (16 tests) with `tests/conftest.py`
fixing the namespace-package shadowing that made `bioprep.core` unimportable from
the repository root.

```bash
cd tests && python -m unittest test_backend -v
```

### Round 3 — analysis accuracy (done, verified)

| ID | Fix | Verification |
|---|---|---|
| B10 | `analyze_structure` reads the first model only, so an NMR ensemble no longer multiplies atom counts, water counts and gaps by the model count. | a 3-model file reports the same counts as 1 model |
| B11 | `detect_missing_residues` maps PDBFixer's chain *index* through `topology.chains()`. It had been reading the key as `(model, chain_id)` and reporting an insert position as the chain. | a gap in chain B is labelled B, not 19 |
| B12 | Pocket volume uses the grid resolution actually used, which `_detect_pockets` may coarsen. | doubling the grid multiplies volume by 8 |
| B13 | Site analysis excludes hydrogens and waters; cavities are defined by heavy atoms. | 327 of 642 atoms used, no hydrogen survives the filter |
| B14 | Pharmacophore features are ranked before the 40-cap, so backbone N/O yield to aromatic and sidechain features. | aromatics survive a cap of 6 |
| B26 | `os.path.splitext` instead of `str.replace('.pdb', ...)`, which substituted every occurrence. | a directory named `my.pdb.data` is left intact |
| B27 | `-xr` added, so the receptor is written rigid with no torsion tree. | ROOT/BRANCH absent; without it, 202 BRANCH records |
| B28 | Ligand CONECT records carried through and remapped. | C1-C2, C2-O1, C2-N1 intact in the final output |
| B34 | `--chain` uses `action="append"`, so a list reaches the selector instead of a string being substring-matched. | `--chain A --chain B` keeps exactly A and B |
| B37 | Dead volume and drugability computation removed from `_detect_pockets`; `analyze()` was discarding it. | — |
| B40 | The web app moved into the package as `bioprep.webapp`, with templates and static files as `package_data` and a `bioprep-web` console script. Mutable state moved out of the package to `BIOPREP_DATA_DIR` or the working directory, since site-packages is often read-only and a reinstall would delete the history. `bioprep/app.py` remains a source-checkout shim. | `GET /` returns 200 from the packaged templates; the shim still exposes `app` and its configuration |
| B42 | `BioPrepSelect.accept_model` writes the first model only. | a 3-model ensemble writes one model's worth of atoms |

Every ID in Table B is now either fixed above or recorded under "Known limits".

---

## Table A — Documented features vs. reality

Claims from `bioprep_deep_dive.md` / `bioprep_technical_pitch.md`.

| # | Documented feature | Status | Evidence |
|---|---|---|---|
| A1 | PDB load/save via Biopython, optional Select filter | EXISTS | `core/io.py` |
| A2 | Water removal | EXISTS | `cleaner.py:25-28` |
| A3 | Structural-water detection, 4.0 Å, protein as reference set | EXISTS (buggy — see B6) | `cleaner.py:62-90` |
| A4 | Chain filtering | EXISTS | `cleaner.py:16-19` |
| A5 | Selective heteroatom removal by residue name | EXISTS | `cleaner.py:41-42` |
| A6 | "Remove ALL heteros" mode | EXISTS (buggy — see B5) | `cleaner.py:37-38` |
| A7 | Surgical ligand preservation (keep-by-default) | EXISTS | `cleaner.py:44-45` |
| A8 | pH-dependent hydrogen addition | EXISTS | `protonator.py:44` |
| A9 | **Protein/ligand separation before PDBFixer** | **ABSENT** | No HETATM handling in `protonator.py`. String search across all 3 project copies returns only prose docs. |
| A10 | **Ligand coordinate preservation via re-merge** | **ABSENT** | Same as A9. Ligands survive only because `removeHeterogens()` was deleted (`protonator.py:29`). |
| A11 | **Ligand protonation status report dict** | **ABSENT** | `add_hydrogens()` returns `None`. |
| A12 | Non-standard residue replacement | EXISTS | `protonator.py:26-27` |
| A13 | Optional loop reconstruction | EXISTS | `protonator.py:18-23` |
| A14 | Optional missing-atom addition | EXISTS | `protonator.py:35-37` |
| A15 | Chain ID preservation (`keepIds=True`) | EXISTS | `protonator.py:50` |
| A16 | Chain/water/hetero/atom-count analysis | EXISTS | `analyzer.py` |
| A17 | Sequence gap detection | EXISTS (buggy — see B10) | `analyzer.py:45-53` |
| A18 | PDBFixer missing-residue detection | EXISTS (mislabeled — see B11) | `analyzer.py:67-82` |
| A19 | Force field selection amber14 / charmm36 | EXISTS | `minimizer.py:23-41` |
| A20 | GBSA implicit solvent (OBC2) | EXISTS (breaks waters — see B2) | `minimizer.py:26` |
| A21 | **"Always loads water/ion params to prevent crashes"** | **ABSENT** | GBSA branch loads no water FF. This is the exact crash it claims to prevent. |
| A22 | **Tier 1/2/3 minimization fallback** | **ABSENT** | Single `try` block. No residue stripping, no GBSA retry. |
| A23 | Final fallback: copy input → output | EXISTS | `minimizer.py:94-95` (this is the *only* fallback) |
| A24 | **`_merge_coords_by_name()` coordinate re-merge** | **ABSENT** | Function does not exist anywhere. |
| A25 | **Residue classification sets (26 AA / 6 water / 24 nucleic)** | **ABSENT** | No such sets. |
| A26 | **Auto terminal-group fixing via `Modeller.addHydrogens`** | **ABSENT** | `Modeller` never imported. |
| A27 | Pre-minimization energy validation | PARTIAL | NaN/Inf only. No 10^12 ceiling — see B4. |
| A28 | Energy before/after/delta reporting | EXISTS | `minimizer.py:81-90` |
| A29 | Convergence assessment | EXISTS (misnomer — see B39) | `minimizer.py:79` — rule is `delta < -0.1`, not the documented `> 1.0` |
| A30 | **Partial-minimization warning string** | **ABSENT** | No `warning` key is ever returned. |
| A31 | Grid pocket detection + dynamic scaling | EXISTS | `site_analyzer.py:85-106` |
| A32 | LIGSITE enclosure check, ≥3 of 6 rays, early exit | EXISTS | `site_analyzer.py:123-158` |
| A33 | DBSCAN clustering | EXISTS | `site_analyzer.py:165` |
| A34 | Volume calc / 50 Å³ floor / top-5 cap | EXISTS (buggy — see B12) | `site_analyzer.py:40-41, 200` |
| A35 | 5-factor drugability score with documented weights | EXISTS | `site_analyzer.py:262-268` |
| A36 | Pharmacophore prediction incl. aromatic centroids | EXISTS (buggy — see B14) | `site_analyzer.py:281-347` |
| A37 | PDBQT export via OpenBabel + Gasteiger | EXISTS (see B27) | `exporter.py:41-47` |
| A38 | GROMACS export | EXISTS (cosmetic — renames a copy only) | `exporter.py:21-25` |
| A39 | Structured + text report | EXISTS | `reporter.py` |
| A40 | 11 Flask API routes | EXISTS | `app.py` |
| A41 | Batch mode, high-throughput mode | EXISTS (HT drops settings — see B20) | `app.py:338, 505` |
| A42 | Job history, persistent results, templates | EXISTS | `app.py:308-328, 638-663` |

**Summary:** of the doc's headline engineering claims, **9 are entirely absent** (A9, A10, A11, A21, A22, A24, A25, A26, A30). The two features the documents present as BioPrep's core value — ligand-safe protonation and tiered minimization — were never implemented.

---

## Table B — Bugs, by severity

### CRITICAL — silently produces wrong or unprocessed scientific output

| ID | File:line | Problem | Status |
|---|---|---|---|
| B1 | `minimizer.py:48` | Any retained ligand has no force-field template → `createSystem` raises → input copied to output unchanged. Minimization **never runs on any structure with a preserved ligand**. BioPrep's flagship feature (keep ligands) and its minimizer are mutually exclusive. | **CONFIRMED** — test appended a 2-atom `LIG`; output byte-identical to input, error `No template found for residue 46 (LIG)`. |
| B2 | `minimizer.py:26` | GBSA branch loads `amber14-all.xml` + `implicit/obc2.xml` — **no water parameters**. Any retained water kills minimization. `use_gbsa` defaults to `True`, so "keep structural waters" + minimize always fails. | **CONFIRMED** — same structure: `GBSA=True` → `No template found for residue (HOH)`; `GBSA=False` → succeeded (−4336 kJ/mol). |
| B3 | `minimizer.py:92-109`, `app.py:224` | On failure the function copies input→output and returns a dict with an `error` key, but `/api/process` still returns HTTP 200 `success: true` and the report still contains an `energy_minimization` block. Caller cannot distinguish "minimized" from "not minimized" without inspecting a nested key. | **CONFIRMED** |
| B4 | `minimizer.py:62` | Only NaN/Inf are guarded. A finite but physically absurd starting energy passes through. | **CONFIRMED** — a run started at **1.045 × 10¹³ kJ/mol** and was reported as a successful, converged minimization. |
| B5 | `app.py:199` + `cleaner.py:37` | `hets_actually_removed = [h for h in hets_to_remove if h not in protect_ligands]`. When the user picks "remove ALL", the token `'ALL'` is not in `protect_ligands`, so it survives filtering and `BioPrepSelect` drops **every** hetero — including explicitly protected ligands. Protection is silently ignored. | **CONFIRMED** |
| B6 | `cleaner.py:89` + `:26` | `structural_waters` stores bare `residue.id` tuples (e.g. `('W', 700, ' ')`), which are only unique *within a chain*. A water in chain B is kept because a same-numbered water in chain A was structural. | **CONFIRMED** — chain B water placed 500 Å from all protein atoms was retained. |
| B7 | `protonator.py:43-46`, `app.py:245` | `addMissingHydrogens` failure is caught and logged, then the file is written anyway from possibly half-modified state; the report hardcodes `hydrogens_added=True`. The report claims protonation that did not happen. | INSPECTION |
| B8 | `app.py:257-262` | `export_structure` returns `(False, error_msg)` on OpenBabel failure; the error is discarded. User downloads a PDB while the report says `docking_target: autodock`. | INSPECTION |
| B9 | `cli.py:27` | Calls `clean_structure(structure, target_chain=...)`; parameter is `target_chains`. Every CLI invocation raises `TypeError`. The `bioprep=bioprep.cli:main` console entry point is dead. | **CONFIRMED** — `TypeError: clean_structure() got an unexpected keyword argument 'target_chain'` |

### HIGH — wrong numbers reported, or breaks under normal use

| ID | File:line | Problem | Status |
|---|---|---|---|
| B10 | `analyzer.py:19-55` | Loops over all models. For NMR ensembles, `atoms_total`, `water_count` and `sequence_gaps` are multiplied by the model count (20× for a typical NMR entry). | INSPECTION |
| B11 | `analyzer.py:77` | Unpacks PDBFixer's `missingResidues` keys as `(model_idx, chain_id)`. They are actually `(chain_index, residue_insert_position)`. Reported "chain" values are residue offsets. | INSPECTION (1CRN has no gaps, so not reproducible on the local test file) |
| B12 | `site_analyzer.py:200` | `voxel_volume = 1.5 ** 3` is hardcoded, but `_detect_pockets` raises `grid_res` on large structures. Every reported volume — and the volume-weighted drugability score — is understated exactly when dynamic scaling triggers. | INSPECTION |
| B13 | `site_analyzer.py:12`, `app.py:700` | `self.atoms` includes hydrogens, waters and ligands. Site analysis on a session structure runs on the *protonated* file, so probe distances are measured to hydrogens and retained waters fill the cavities. | INSPECTION |
| B14 | `site_analyzer.py:307-347` | Backbone `N` and `O` match the donor/acceptor lists, so every pocket residue contributes 2 features. `points[:40]` truncates in iteration order and aromatic ring centroids are appended *after* the loop — so the aromatic feature is discarded in essentially every real pocket. | INSPECTION |
| B15 | `app.py:113, 181, 185-187` | Intermediates use fixed names in the shared system temp dir (`raw_{filename}`, `temp_{base}.pdb`, `prot_{base}.pdb`). Two concurrent requests for the same filename overwrite each other mid-pipeline. | INSPECTION |
| B16 | `minimizer.py:83-85` | If minimization yields NaN, `round(nan, 1)` propagates into `jsonify`, which emits a bare `NaN` literal — invalid JSON that throws in any strict parser. | INSPECTION |

### MEDIUM — leaks, dropped settings, inaccurate reporting

| ID | File:line | Problem |
|---|---|---|
| B17 | `app.py:369, 535` | `batch_dir` / `ht_dir` are never removed. **22 stale `bioprep_batch_*` / `bioprep_ht_*` directories currently sit in the temp folder.** |
| B18 | `app.py:181-187` | `/api/process` never deletes `raw_`, `temp_`, `prot_` or final files. |
| B19 | `app.py:40` | `session_pdb_paths` grows without bound for the process lifetime. |
| B20 | `app.py:577-581` | High-throughput never passes `target_chains` or `keep_structural_waters`, and parses `docking_target` (line 533) without ever calling `export_structure`. Three UI settings are silently dropped in this mode. |
| B21 | `app.py:200` | `hets_retained` is computed against `hets_actually_removed`; with `'ALL'` it still lists every ligand as retained after they were all deleted. |
| B22 | `app.py:242` | `waters_removed = water_count_before` ignores structural waters that were kept — overcounts. |
| B23 | `analyzer.py:62` | `atoms_before` is commented "protein-only pre-processing count" but includes waters and heteroatoms; the report's `delta` conflates water removal with hydrogen addition. |
| B24 | `app.py:670-681` | `/api/download-report` passes arbitrary user JSON into `report_to_text`, which indexes `report['chains']` / `['atom_counts']` directly → `KeyError` → 500. |
| B25 | `reporter.py:74-92` | Same direct indexing; no `.get()` defaults. |
| B26 | `exporter.py:22` | `output_path.replace('.pdb', '_gromacs.pdb')` replaces *every* occurrence; use `os.path.splitext`. |
| B27 | `exporter.py:41-47` | Receptor conversion omits OpenBabel's `-xr` rigid-receptor flag, so the PDBQT may carry a torsion tree instead of receptor formatting. Needs verification against AutoDock. |
| B28 | `protonator.py:49-50` | `PDBFile.writeFile` emits no `CONECT` records, so ligand bond orders/connectivity are lost downstream. |
| B29 | `app.py:60-74, 55-57` | `jobs_history.json` and `templates_store.json` are read-modify-written with no locking; concurrent requests corrupt them. History is already 345 KB because full reports are embedded. |
| B30 | `app.py:269-273` | The structure is base64-encoded **twice** (viewer + download) into one JSON response. With the 200 MB upload cap this is a memory blow-up. |
| B31 | `app.py:724` | `app.run(debug=True)` — Werkzeug interactive debugger. |
| B32 | `app.py:130, 305, 498` | Raw `str(e)` returned to clients, leaking absolute filesystem paths. |
| B33 | `app.py:545-546` | `extractall` has no file-count or uncompressed-size limit (zip-bomb). Path traversal itself is mitigated by CPython's sanitization. |
| B34 | `cli.py:12, 27` | Even after fixing B9, `--chain A` passes a *string* where a list is expected; `chain.id not in "AB"` does substring matching. Wrap in a list. |

### LOW — dead code and hygiene

| ID | File:line | Problem |
|---|---|---|
| B35 | `app.py:44, 30` | `allowed_file()` / `ALLOWED_EXTENSIONS` defined, never called; each route re-implements its own check. |
| B36 | `app.py:4, 9, 13, 39` | `queue`, `threading`, `Response`, `stream_with_context` imported but unused; `batch_progress` is only ever `.pop()`ed, never written. Remnants of a removed SSE progress feature. |
| B37 | `site_analyzer.py:176-190` | `_detect_pockets` computes `volume` and `drugability_score` that `analyze()` immediately discards by recomputing in `_analyze_specific_site`. |
| B38 | `minimizer.py:86` | `iterations: 1000` is reported unconditionally regardless of actual iterations used. |
| B39 | `minimizer.py:79` | `converged` really means "energy decreased by > 0.1 kJ/mol" — not convergence. |
| B40 | `setup.py` | `find_packages()` captures only `bioprep/`; `app.py`, `templates/`, `static/` are not installed, so a `pip install`ed copy cannot serve the web app. `numpy` missing from `install_requires`. |
| B41 | repo root | `error.txt` and `test_out.txt` are stale UTF-16 PowerShell captures from the old `Documents/AG PLUGINS` path, describing bugs already fixed (duplicate `delete_template` route, `GBSAOBC2`, `find_binding_sites`). Misleading; delete. |
| B42 | `io.py:31-36` | `PDBIO.save` writes all NMR models; downstream PDBFixer silently uses the first. |
| B43 | `bioprep/tmp/`, root | Ad-hoc test scripts (`test_audit.py`, `test_1stp.py`, `stress_test*.py`) use `sys.exit(1)` and hardcoded relative paths; only `tests/test_basic.py` is a real unittest. No CI-runnable suite. |

---

## Suggested fix order

1. **B9, B34** — restore the CLI (one-line fix, unblocks scripted testing).
2. **B1, B2, B3, B4** — the minimizer cluster. This is where the doc's Tier 2/3 design actually earns its keep: strip unparameterizable residues, minimize, merge coordinates back by `(chain, resname, resid, atomname)`. Fixing B3 first makes the others visible instead of silent.
3. **B5, B6, B7, B8** — correctness of what the report claims happened.
4. **B15, B17, B18, B19** — temp-file and concurrency hygiene, before any multi-user deployment.
5. **B10–B14** — reported-number accuracy.
6. Everything else.

---

## Feature catalogue verification — complete

All 121 features from `bioprep_deep_dive.md` were checked against the running
code, group by group, with an independent verification script per group.

| Features | Module | Checks | Result |
|---|---|---|---|
| 1–4 | `io.py` | 29 (with 5–11) | fixed: empty-file guard, mmCIF, element column, CONECT |
| 5–11 | `cleaner.py` | — | fixed: `OH2` water oxygens, case-insensitive names |
| 12–22 | `protonator.py` | 28 | fixed: ligand status reason |
| 23–28 | `analyzer.py` | 28 | no new defects |
| 29–48 | `minimizer.py` | 52 | added: terminal repair, CPU platform, `converged` |
| 49–69 | `site_analyzer.py` | 51 | no new defects |
| 70–80 | `exporter.py`, `reporter.py` | 63 | no new defects |
| 81–110 | `app.py` | 95 | fixed: SEQRES loss disabling loop reconstruction |
| 111–112 | `cli.py` | 58 | fixed: mmCIF working directory leak |
| 113–121 | frontend | 1 finding | contract drift, see below |

**Total: 404 verification checks, 124 committed regression tests.**

### Defects found and fixed

Six produced wrong output with no error at all, which is the category that does
real damage in a preparation tool:

1. Two-character elements were written one column left of the spec, so a zinc
   cofactor read back as nitrogen and iron as fluorine — in Biopython and
   OpenMM alike.
2. Structural-water detection matched only `O`, `OW`, `O1`, so CHARMM and NAMD
   structures (`OH2`) had no detectable water oxygen and every water was
   deleted despite the user asking to keep them.
3. `--protect lig` failed to match `LIG` and the ligand was deleted by
   `--remove-het ALL` without a word.
4. CONECT records were dropped at the cleaning step, making the protonator's
   ligand-bond preservation dead code.
5. SEQRES was dropped at the same step, so PDBFixer could not tell which
   residues were missing and `reconstruct_loops` rebuilt nothing whether it
   was on or off.
6. Non-PDB uploads blamed the chain-selection settings instead of the file.

Plus: the minimizer's terminal-group repair, CPU platform preference and
`converged` flag were absent; the ligand status report did not say why a ligand
was skipped; and the CLI leaked a working directory on every mmCIF run.

### Frontend contract (for the rewrite)

`static/js/main.js` still reads the pre-rewrite report shape. The current
`/api/process` response is:

```
{
  success, filename, session_id,
  viewer_pdb_b64,          # always a PDB, for the 3D viewer
  pdb_b64,                 # the download; equals viewer_pdb_b64 unless a
                           # docking export produced a different file
  warnings: [str],
  report_text: str,
  report: {
    generated_at, input_file,
    chains: {detected: [], retained: []},
    water_molecules_removed, water_molecules_retained,
    heteroatoms: {removed: [], retained: []},
    protonation: {
      hydrogens_added: bool, ph: float,
      ligands_preserved: [], ligand_status: [{residue, atoms, action, reason}],
      nonstandard_replaced: [], loops_reconstructed: int
    },
    hydrogens_added: bool,          # kept at top level for older readers
    protonation_ph: float,
    missing_residues_detected: [{chain, residue, position}],
    atom_counts: {before_processing, after_processing, delta},
    warnings: [str],
    energy_minimization: {          # present only when minimisation ran
      status: 'full'|'partial'|'partial_no_implicit_solvent'|'failed',
      force_field, gbsa_used,
      energy_before_kJ_mol, energy_after_kJ_mol, delta_energy_kJ_mol,
      energy_decreased: bool, converged: bool,
      excluded_residues: [], terminals_repaired: bool,
      iterations_max, energy_tolerance_kJ_mol_nm,
      warnings: [], error: str|None
    },
    docking_export: {target, succeeded, error?},
    docking_target,                 # only when succeeded
    settings_used: {...},
    processing_time_seconds
  }
}
```

Two changes the old frontend has not caught up with:

- Energies are `float` or `None`, never the string `'N/A'`. `main.js` still
  tests for `'N/A'`.
- `iterations` was replaced by `iterations_max`; the actual iteration count is
  not reported because OpenMM does not expose it.

A returned file is never proof that minimisation ran — check
`energy_minimization.status`.

---

## Science audit

Feature verification asks whether the code does what the catalogue says.
This asks whether the answer is physically defensible. It was run against
deposited structures — 1STP (streptavidin + biotin), 1HSG (HIV-1 protease +
MK1), 4INS (insulin + zinc) — because 1CRN is small, capped, ligand-free and
pocket-free, and hid every defect below.

### Defects found and fixed

| Finding | Evidence | Fix |
|---|---|---|
| Structures with an uncapped C-terminus cannot be minimised at all | 1STP has **zero OXT atoms**; every tier failed. 1HSG likewise. | Terminals are always repaired via PDBFixer's `missingTerminals`, separately from the optional side-chain rebuilding. 1STP now reaches −15149 kJ/mol. |
| Pocket detection depended on molecular orientation | Rotating streptavidin 45° changed total cavity volume by **+47.8%** for a physically identical cavity | Scan all seven LIGSITE axes (three Cartesian + four cubic diagonals) in both senses. Worst drift over five rotations: **8.9%**, the remainder being the axis-aligned grid. |
| Tier 2 let the pocket collapse onto the restored ligand | With biotin deleted, 27 of 67 lining atoms moved inward; closest contact tightened 2.58 → 2.17 Å | Harmonic restraints (5000 kJ/mol/nm²) on atoms within 5 Å of a deleted residue: 7 of 67, closest 2.45 Å |
| "Structural" waters kept with no hydrogen bond | 8 of 78 retained waters on 1STP had no N/O/S within 3.5 Å — the rule counted carbon | Reference set is now hydrogen-bond-capable atoms at 3.5 Å. Zero spurious across all three structures. |
| PDBQT receptors typed every hydrogen as a donor | `-xh` means *preserve* hydrogens, not merge non-polar ones. 842 hydrogens all typed `HD` versus 213 genuinely polar. | Option removed; the default merges non-polar hydrogens, giving the conventional receptor |
| `converged` measured the wrong thing | It reported the energy having fallen, so a 100-iteration run "converged" while 1000 iterations went 900 kJ/mol lower | Compares RMS force against the tolerance `minimizeEnergy` was given. 1CRN: 5 iters → 695, False; 100 → 15.4, False; 1000 → 5.05, True |
| Bond constraints made convergence unmeasurable | `HBonds` put constraint forces into `getForces()`; RMS force never approached the tolerance | Minimisation is unconstrained — nothing here integrates dynamics. 1CRN: −5166.9 kJ/mol and RMS 5.76 versus −5162.7 and 71.16 |
| Chain selection silently defeated ligand protection | Select's `accept_chain` is an all-or-nothing gate: rejecting a chain skips `accept_residue` for everything in it. A ligand routinely sits on a different chain than the receptor — HIV-1 protease's inhibitor MK1 is chain B of the A/B dimer — so `chains=['A']` dropped an explicitly protected ligand with no warning; the report just showed `ligands_preserved: []` | `accept_chain` now always accepts; the chain decision moved into `accept_residue`, checked *after* the protection override so a named ligand survives regardless of which chain it is recorded under. An unprotected heteroatom on an unselected chain is still dropped — confirmed by a dedicated test |
| A chain end missing backbone atoms, not just its sidechain, crashed protonation outright | With `add_missing_atoms` off (the default), `fixer.missingAtoms` was cleared to `{}` entirely on the theory that everything in it is optional sidechain remodelling. It isn't: a residue with no resolved density past its amide N is missing CA/C/O too, and PDBFixer's own terminal-placement code reads `atomPositions['O']` and `['CA']` unconditionally when a `missingTerminals` entry exists for that residue — wiping the dict removed the record that those atoms needed adding, so `addMissingAtoms()` raised `KeyError: 'O'`. Hit 4 of 89 real structures pulled live from RCSB across an "industrial-scale" sweep (zinc fingers, protein–DNA/RNA complexes, heme/Fe–S/Cu proteins, glycoproteins, membrane proteins, NMR ensembles): `1PRC`, `1BGB`, `1K8W`, `2DRP` | Only sidechain atoms are dropped when the flag is off; backbone atoms (N/CA/C/O) are kept regardless, same reasoning as why `missingTerminals` itself is never cleared. All 4 process cleanly now; regression test reproduces it from a synthetic terminal residue rather than depending on any one PDB file |

Found by `verification/stress_combinations.py`, which exercises several
settings together on real complexes rather than one or two at a time — this
particular interaction only shows up with chain selection and ligand
protection combined, and every isolated test of either had passed.

The backbone-atom finding above was found by a from-scratch sweep of 89 real
structures fetched live from RCSB (`search.rcsb.org` for protein–DNA,
protein–RNA and NMR-method structures; hand-picked PDB IDs for zinc fingers,
heme/Fe–S/Cu proteins, glycoproteins, membrane proteins, and small-globular
controls), run through the actual `/api/high-throughput` endpoint with
default settings. 84/89 succeeded on the first pass; the other 5 were the 4
`KeyError('O')` cases above plus one structure whose *download* had been
truncated by a transient network issue during the fetch (not a bioprep bug —
a fresh download of the same entry, PIT-1/DNA complex `1AU7`, processed
correctly). A rerun after the fix: 89/89 succeeded.

### Verified correct

- Disulfide cysteines are left oxidised; no HG is added to a bridged SG.
- Ligand coordinates are byte-identical through the pipeline, and no hydrogens are added to them.
- Pocket detection finds the real site: with the ligand deleted first, 1STP's **top-ranked** pocket sits 6.3 Å from the biotin centroid, lined by the actual binding residues.
- LYS, HIS, ASP and GLU titrate at approximately the right pH.
- Gasteiger charges are chemically sensible: O and N negative, carbon near zero.

### Known limits, not fixed

- **Tyrosine is never deprotonated.** Confirmed by inspecting the template inventories directly, not inferred: `ForceField('amber14-all.xml')` defines `TYR`, `NTYR`, `CTYR` and nothing else for tyrosine — no anionic form, unlike histidine's `HID`/`HIE`/`HIP` or lysine's `LYS`/`LYN`. CHARMM36 has an `STYR` template that looked like a candidate, but it has 16 atoms against `TYR`'s 21 and **no OH or HH atom at all** — it is some other truncated variant, not a tyrosinate. Neither force field this tool offers gives PDBFixer anywhere to select a deprotonated state from, so this cannot be fixed by changing how the protonator calls it; it would need a different hydrogen-placement engine entirely.
- **pKa values are model-compound values, by the library's own design.** `PDBFixer.addMissingHydrogens`'s docstring states it directly: *"No extensive electrostatic analysis is performed; only default residue pKas are used. The pH is only taken into account for standard amino acids."* A buried or salt-bridged residue with a pKa shifted by several units gets the model-compound state regardless — worst exactly where it matters, at a catalytic residue. Fixing this means adding a structure-specific predictor (PROPKA or similar) ahead of PDBFixer, not a change to this tool's own code.
- **The drugability score is an unvalidated heuristic.** On 1HSG the real inhibitor site is found (6.2 Å, largest volume, lined by the catalytic Asp dyad and flap) but ranks **last of five**: the volume term penalises it for exceeding 300–1000 Å³ and concavity for being an open cavity — the very properties that let it bind a peptidomimetic. Two of the five factors, property diversity and pharmacophore density, are 1.00 for every pocket and do no discriminating work. The five sub-scores are returned per pocket so the ranking can be argued with; the weights were deliberately **not** retuned, since fitting them to two structures would be overfitting dressed as improvement.
- **Ligands are never parameterised.** Tier 2 excludes them rather than generating GAFF or CGenFF parameters, so a ligand's internal geometry is never optimised.
- **The element column's case is inconsistent between writers.** RCSB and our own `SpecCompliantPDBIO` write `ZN`; OpenMM's own `PDBFile.writeFile`, used for every post-protonation and post-minimisation file, writes `Zn`. Both identify the correct element to every tool checked (OpenMM, OpenBabel), unlike the earlier column-offset bug that actually changed which element was read. Left alone rather than post-processing every OpenMM-written file to re-case a column that nothing downstream cares about.
