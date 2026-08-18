# BioPrep Backend Audit

Audit date: 2026-08-18. Scope: Python backend only (`bioprep/app.py`, `bioprep/bioprep/**`).
Frontend (`static/`, `templates/`) intentionally excluded — it is being rewritten.

Environment used for verification: `.venv` (Python 3.13, biopython 1.86, OpenMM 8.4.0,
pdbfixer 1.12.0, numpy 2.4.3, scipy 1.17.1, scikit-learn 1.8.0).

Status legend: **CONFIRMED** = reproduced by running code. **INSPECTION** = read from source,
not yet executed.

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
