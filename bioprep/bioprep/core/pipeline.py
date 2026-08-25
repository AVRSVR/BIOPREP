"""
The preparation pipeline, in one place.

Precision, batch and high-throughput modes previously each carried their own
copy of these steps. The copies drifted: high-throughput quietly stopped
honouring chain selection, structural-water preservation and docking export.
Every mode now calls :func:`prepare_structure`, so a setting either works
everywhere or nowhere.
"""

import os
import shutil
import time
from dataclasses import dataclass, field, asdict
from typing import List, Optional

from .io import load_pdb, save_pdb, ensure_pdb
from .cleaner import clean_structure
from .protonator import add_hydrogens
from .analyzer import analyze_structure, detect_missing_residues
from .minimizer import minimize_structure
from .exporter import export_structure
from .reporter import build_report, count_atoms_in_pdb

DEFAULT_PH = 7.4


@dataclass
class PipelineSettings:
    """Every knob the pipeline exposes. One definition, shared by all modes."""

    chains: Optional[List[str]] = None
    remove_water: bool = True
    remove_heteros: List[str] = field(default_factory=list)
    protect_ligands: List[str] = field(default_factory=list)
    ph: float = DEFAULT_PH
    keep_structural_waters: bool = False
    reconstruct_loops: bool = False
    add_missing_atoms: bool = False
    run_minimization: bool = False
    use_gbsa: bool = True
    force_field: str = 'amber14'
    docking_target: Optional[str] = None
    use_propka: bool = False

    @classmethod
    def from_mapping(cls, data):
        """
        Build settings from a dict of already-decoded values (batch/HT config
        JSON). Unknown keys are ignored; malformed values fall back to the
        default rather than raising.
        """
        data = data or {}

        def as_bool(key, default):
            value = data.get(key, default)
            if isinstance(value, str):
                return value.strip().lower() == 'true'
            return bool(value)

        def as_list(key):
            value = data.get(key) or []
            if isinstance(value, str):
                value = [value]
            return [str(v).strip().upper() for v in value if str(v).strip()]

        try:
            ph = float(data.get('ph', DEFAULT_PH))
        except (TypeError, ValueError):
            ph = DEFAULT_PH
        # Outside this range PDBFixer's protonation is not meaningful.
        ph = min(max(ph, 0.0), 14.0)

        chains = data.get('chains') or None
        if isinstance(chains, str):
            chains = [chains]
        if chains:
            chains = [str(c).strip() for c in chains if str(c).strip()]
        if not chains:
            chains = None

        force_field = str(data.get('force_field', 'amber14')).lower()
        if force_field not in ('amber14', 'charmm36'):
            force_field = 'amber14'

        docking = data.get('docking_target') or None
        if docking:
            docking = str(docking).lower()
            if docking not in ('autodock', 'vina', 'gromacs', 'pdb'):
                docking = None

        return cls(
            chains=chains,
            remove_water=as_bool('remove_water', True),
            remove_heteros=as_list('remove_heteros'),
            protect_ligands=as_list('protect_ligands'),
            ph=ph,
            keep_structural_waters=as_bool('keep_structural_waters', False),
            reconstruct_loops=as_bool('reconstruct_loops', False),
            add_missing_atoms=as_bool('add_missing_atoms', False),
            run_minimization=as_bool('run_minimization', False),
            use_gbsa=as_bool('use_gbsa', True),
            force_field=force_field,
            docking_target=docking,
            use_propka=as_bool('use_propka', False),
        )

    def to_dict(self):
        return asdict(self)


def prepare_structure(input_path, workdir, settings, original_filename=None):
    """
    Run clean -> protonate -> (minimise) -> (export) on one structure.

    All intermediates are written inside ``workdir``, which the caller owns and
    is responsible for removing.

    Returns a dict::

        {
          'report':        {...},          # structured preparation report
          'viewer_path':   '<path>.pdb',   # always a PDB, for 3D display
          'download_path': '<path>',       # PDBQT when docking export succeeded
          'warnings':      [str, ...],
        }

    Raises on unrecoverable errors (unreadable PDB, everything filtered out).
    """
    started = time.time()
    warnings = []
    filename = original_filename or os.path.basename(input_path)
    base = os.path.splitext(os.path.basename(filename))[0]

    cleaned_path = os.path.join(workdir, f'{base}__cleaned.pdb')
    protonated_path = os.path.join(workdir, f'{base}__protonated.pdb')
    final_path = os.path.join(workdir, f'{base}__prepared.pdb')

    # ── Pre-processing analysis ──────────────────────────────────────────────
    # mmCIF is converted once, here, so every step below sees PDB.
    # PDBFixer and OpenMM both require it, and detect_missing_residues
    # reads the path directly rather than the parsed structure.
    input_path, converted_from = ensure_pdb(input_path, workdir)
    if converted_from:
        warnings.append(
            f"Input was {converted_from.upper()} and was converted to "
            "PDB for processing; the output is PDB."
        )

    structure = load_pdb(input_path)
    pre = analyze_structure(structure)
    missing_before = detect_missing_residues(input_path)

    # A protected ligand always wins, so it must never appear in the removal
    # list. 'ALL' is handled inside BioPrepSelect, which checks protection first.
    hets_to_remove = [
        het for het in settings.remove_heteros
        if het not in settings.protect_ligands
    ]

    # ── Clean ────────────────────────────────────────────────────────────────
    select_obj = clean_structure(
        structure,
        target_chains=settings.chains,
        remove_water=settings.remove_water,
        remove_heteroatoms=hets_to_remove,
        keep_structural_waters=settings.keep_structural_waters,
        protect_ligands=settings.protect_ligands,
    )
    # source_pdb carries across what Biopython drops on parse: the input's
    # CONECT bonds, so the protonator's ligand-bond preservation has something
    # to preserve, and its SEQRES, without which PDBFixer cannot tell which
    # residues are missing and loop reconstruction silently rebuilds nothing.
    save_pdb(structure, cleaned_path, select=select_obj,
             source_pdb=input_path)

    if count_atoms_in_pdb(cleaned_path) == 0:
        raise ValueError(
            "Cleaning removed every atom. Check the chain selection and "
            "heteroatom removal options."
        )

    waters_kept = len(select_obj.structural_waters)
    waters_removed = (pre['water_count'] - waters_kept) if settings.remove_water else 0

    # ── Protonate ────────────────────────────────────────────────────────────
    protonation = add_hydrogens(
        cleaned_path, protonated_path,
        ph=settings.ph,
        reconstruct_loops=settings.reconstruct_loops,
        add_missing_atoms=settings.add_missing_atoms,
        use_propka=settings.use_propka,
    )
    warnings.extend(protonation['warnings'])

    # ── Minimise (optional) ──────────────────────────────────────────────────
    minimization = None
    if settings.run_minimization:
        minimization = minimize_structure(
            protonated_path, final_path,
            force_field=settings.force_field,
            use_gbsa=settings.use_gbsa,
        )
        warnings.extend(minimization['warnings'])
        if minimization['error']:
            warnings.append(
                "Energy minimization did not run; the protonated structure was "
                "returned unchanged."
            )
    else:
        shutil.copy2(protonated_path, final_path)

    # ── Export (optional) ────────────────────────────────────────────────────
    viewer_path = final_path
    download_path = final_path
    export_status = None
    if settings.docking_target:
        ok, result = export_structure(
            final_path,
            os.path.join(workdir, f'{base}__export.pdb'),
            settings.docking_target,
        )
        if ok:
            download_path = result
            export_status = {'target': settings.docking_target, 'succeeded': True}
        else:
            # Previously this failure was discarded and the report still
            # claimed a docking export had happened.
            export_status = {
                'target': settings.docking_target,
                'succeeded': False,
                'error': str(result),
            }
            warnings.append(f"Docking export failed: {result}")

    # ── Report ───────────────────────────────────────────────────────────────
    all_hets = pre['heteroatoms']
    if 'ALL' in hets_to_remove:
        removed = [h for h in all_hets if h not in settings.protect_ligands]
    else:
        removed = [h for h in hets_to_remove if h in all_hets]
    retained = [h for h in all_hets if h not in removed]

    missing_after = (
        detect_missing_residues(final_path)
        if (settings.reconstruct_loops or settings.add_missing_atoms)
        else missing_before
    )

    report = build_report(
        filename=filename,
        chains_detected=pre['chains'],
        chains_retained=settings.chains if settings.chains else pre['chains'],
        waters_removed=waters_removed,
        waters_retained=waters_kept,
        heteroatoms_removed=removed,
        heteroatoms_retained=retained,
        protonation=protonation,
        ph_used=settings.ph,
        missing_residues=missing_after,
        atoms_before=pre['atoms_total'],
        atoms_after=count_atoms_in_pdb(final_path),
        minimization_stats=minimization,
        docking_export=export_status,
        settings=settings.to_dict(),
        processing_time_s=time.time() - started,
        warnings=warnings,
    )

    return {
        'report': report,
        'viewer_path': viewer_path,
        'download_path': download_path,
        'warnings': warnings,
    }
