"""
Energy minimisation with a graded fallback.

Real PDB entries routinely contain residues no protein force field can
parameterise — novel ligands, unusual cofactors, modified bases. A single
`createSystem` call fails on those, so this module degrades in stages:

  Tier 1  minimise everything as-is
  Tier 2  minimise only force-field-safe residues, then splice the moved
          coordinates back into the full structure by atom identity
  Tier 3  Tier 2 again with implicit solvent switched off

If every tier fails the input is passed through unchanged and ``status`` says
so. Callers must check ``status`` — a returned file is not proof of work.
"""

import logging
import math
import os
import shutil

import numpy as np
import openmm as mm
from openmm import app
from openmm import unit

from .residues import is_force_field_safe
from .io import normalize_element_column_case

logger = logging.getLogger(__name__)

MAX_ITERATIONS = 1000
ENERGY_TOLERANCE = 10.0          # kJ/mol/nm

# No bond constraints. Constraining X-H bonds exists to permit a longer MD
# timestep; nothing here integrates dynamics - the Langevin integrator is
# present only because Simulation requires one. During a pure minimisation
# constraints do two unhelpful things: they stop hydrogens relaxing, and the
# forces reported by getForces() then carry constraint contributions, so the
# residual force cannot be compared against the tolerance at all. Measured on
# 1CRN after 1000 iterations:
#
#                   final energy       RMS force
#   HBonds          -5162.7 kJ/mol     71.16 kJ/mol/nm   (never reaches 10)
#   unconstrained   -5166.9 kJ/mol      5.76 kJ/mol/nm   (converged)
MINIMISATION_CONSTRAINTS = None
SUSPECT_ENERGY = 1.0e12          # kJ/mol — beyond this, geometry is suspect

# Tier 2 deletes the residues no force field can parameterise, which leaves the
# site they occupied empty. These hold the surrounding atoms near their input
# positions so the site cannot relax inward while its occupant is missing.
# Chosen by measurement on streptavidin with biotin deleted, against an input
# closest protein-ligand contact of 2.58 A:
#
#   k        closest   lining atoms      worst      final energy
#            contact   moving inward     inward     (kJ/mol)
#   0        2.17 A    27 of 67          -0.63 A    -15117
#   500      2.37 A    23 of 67          -0.40 A    -15012
#   5000     2.45 A     7 of 67          -0.18 A    -14683
#
# 5000 nearly removes the collapse for about three percent less energy
# reduction, which is the better trade: the point of the tier is to relieve
# strain, not to let the site close on an absent ligand.
POCKET_RESTRAINT_CUTOFF = 5.0    # angstrom around a deleted residue
POCKET_RESTRAINT_K = 5000.0      # kJ/mol/nm^2 applied to those atoms

STATUS_FULL = 'full'
STATUS_PARTIAL = 'partial'
STATUS_PARTIAL_NO_GBSA = 'partial_no_implicit_solvent'
STATUS_FAILED = 'failed'


def _build_forcefield(force_field, use_gbsa):
    """
    Assemble a ForceField. Water parameters are always included — omitting
    them is what made retained crystallographic waters abort minimisation.
    """
    if force_field == 'charmm36':
        return app.ForceField('charmm36.xml', 'charmm36/water.xml')

    if use_gbsa:
        # tip3pfb supplies the water templates that implicit/obc2 lacks.
        return app.ForceField(
            'amber14-all.xml', 'implicit/obc2.xml', 'amber14/tip3pfb.xml'
        )
    return app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')


def _atom_key(atom):
    """Identity of an atom that survives topology rebuilds."""
    residue = atom.residue
    return (residue.chain.id, residue.name, residue.id, atom.name)


def _make_simulation(topology, system, integrator):
    """
    Build a Simulation, preferring a GPU platform when one actually works.

    Measured on this project's own hardware (GTX 1650, OpenCL - no CUDA build
    installed): identical starting energy, final energies within 2 kJ/mol,
    mean atom position deviation 0.011 A against a CPU run of the same
    system - well inside anything that matters structurally - for a 3.8x
    wall-clock speedup on a small (642-atom) system. Larger systems typically
    benefit more, not less, since fixed per-step overhead matters less.

    A platform being *listed* by OpenMM doesn't mean it works - CUDA or
    OpenCL can be present but mis-configured, and that only fails at Simulation
    (i.e. Context) construction, not at Platform.getPlatformByName(). So each
    tier is a real construction attempt, not a name check, falling through to
    the next on any failure: CUDA, then OpenCL, then CPU, then whatever
    OpenMM picks by default. CPU no longer being the fixed, guaranteed-
    reproducible choice is a real trade - the same input can now come out
    minimized to a very slightly different local geometry depending on what
    hardware happened to run it - accepted here because the deviation is far
    below anything the rest of this tool treats as structurally meaningful.

    BIOPREP_MINIMIZER_PLATFORM restricts the cascade to one named platform -
    set by the test suite (see tests/conftest.py) so minimizing crambin
    thirty-odd times doesn't pay GPU context/kernel-compile overhead thirty
    times over. Real usage leaves it unset.

    Returns (simulation, platform_name_used).
    """
    forced = os.environ.get('BIOPREP_MINIMIZER_PLATFORM')
    tiers = (
        ('CUDA', {'Precision': 'mixed'}),
        ('OpenCL', {'Precision': 'mixed'}),
        ('CPU', {}),
    )
    if forced:
        tiers = tuple(t for t in tiers if t[0] == forced) or ((forced, {}),)

    for name, properties in tiers:
        try:
            platform = mm.Platform.getPlatformByName(name)
            sim = app.Simulation(topology, system, integrator, platform, properties)
            return sim, name
        except Exception as exc:
            logger.debug("%s platform unavailable (%s); trying the next tier", name, exc)
    logger.debug("No named platform worked; using OpenMM's own default")
    return app.Simulation(topology, system, integrator), 'default'


def _create_system(topology, positions, forcefield, allow_terminal_fix=True):
    """
    Build the OpenMM system, repairing terminal groups if the first try fails.

    A chain truncated mid-structure often lacks OXT, and a residue missing its
    terminal hydrogens matches no template, so createSystem raises. Running the
    topology through Modeller.addHydrogens caps those ends and usually makes it
    parameterisable. Returns (system, topology, positions), where the topology
    may have gained atoms.
    """
    try:
        system = forcefield.createSystem(
            topology, nonbondedMethod=app.NoCutoff, constraints=MINIMISATION_CONSTRAINTS)
        return system, topology, positions, False
    except Exception as exc:
        if not allow_terminal_fix:
            raise
        logger.debug("createSystem failed (%s); attempting terminal repair", exc)

        modeller = app.Modeller(topology, positions)
        modeller.addHydrogens(forcefield=forcefield)
        system = forcefield.createSystem(
            modeller.topology, nonbondedMethod=app.NoCutoff,
            constraints=MINIMISATION_CONSTRAINTS)
        return system, modeller.topology, modeller.positions, True


def _add_position_restraints(system, positions, indices):
    """
    Pin the given atoms near their input coordinates with a harmonic well.

    Used to hold a pocket open while the residue that filled it is absent.
    """
    if not indices:
        return
    restraint = mm.CustomExternalForce(
        '0.5*k_pocket*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)')
    restraint.addGlobalParameter(
        'k_pocket',
        POCKET_RESTRAINT_K * unit.kilojoule_per_mole / unit.nanometer ** 2)
    for name in ('x0', 'y0', 'z0'):
        restraint.addPerParticleParameter(name)

    reference = positions.value_in_unit(unit.nanometer)
    for index in indices:
        x, y, z = reference[index]
        restraint.addParticle(int(index), [x, y, z])
    system.addForce(restraint)


def _run(topology, positions, forcefield, restrain_indices=()):
    """
    Minimise one system.

    Returns (positions, energy_before, energy_after, topology, repaired,
    rms_force, platform_name). The topology comes back because terminal
    repair may have added atoms.
    """
    system, topology, positions, repaired = _create_system(
        topology, positions, forcefield)

    # Only safe when the topology was not rebuilt; a repair renumbers atoms.
    if restrain_indices and not repaired:
        _add_position_restraints(system, positions, restrain_indices)

    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picosecond
    )
    simulation, platform_name = _make_simulation(topology, system, integrator)
    simulation.context.setPositions(positions)

    before = (simulation.context.getState(getEnergy=True)
              .getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
    if math.isnan(before) or math.isinf(before):
        raise ValueError(
            f"Initial potential energy is {before}. The structure has severe "
            "atomic clashes or overlapping atoms; minimisation cannot start."
        )

    simulation.minimizeEnergy(
        maxIterations=MAX_ITERATIONS, tolerance=ENERGY_TOLERANCE
    )

    state = simulation.context.getState(getEnergy=True, getPositions=True,
                                        getForces=True)
    after = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    if math.isnan(after) or math.isinf(after):
        raise ValueError(
            f"Minimisation produced a non-finite energy ({after}); the result "
            "is unusable."
        )

    # Whether the minimiser actually finished is a question about the residual
    # force, not about the energy having gone down. A run stopped early by the
    # iteration cap still lowers the energy, so an energy-based flag reports
    # success for a structure that is nowhere near a minimum.
    # OpenMM halts when the root-mean-square of all force components reaches
    # the tolerance, so that is the quantity to compare against - not the
    # largest per-atom force, which the HBonds constraints keep high on
    # hydrogens regardless of how well converged the structure is.
    forces = np.asarray(state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer))
    rms_force = float(np.sqrt(np.mean(forces ** 2)))

    return state.getPositions(), before, after, topology, repaired, rms_force, platform_name


def _safe_subset(pdb):
    """
    Build a topology containing only force-field-safe residues.

    Also works out which surviving atoms line the space the deleted residues
    occupied. Those atoms are restrained during the tier-2 minimisation:
    without the ligand present there is nothing holding the pocket open, so
    the walls relax inward and the structure that comes back - ligand restored
    at its crystallographic position - contains contacts that were never in
    the input. Measured on streptavidin, 22 of 67 pocket-lining atoms moved
    toward biotin and the closest contact tightened from 2.58 to 2.33 A.

    Returns (modeller, excluded_resnames, lining_atom_indices), or
    (None, [], []) when nothing was excluded and tier 2 would just repeat
    tier 1.
    """
    excluded = sorted({
        residue.name for residue in pdb.topology.residues()
        if not is_force_field_safe(residue.name)
    })
    if not excluded:
        return None, excluded, []

    positions = np.array(pdb.positions.value_in_unit(unit.nanometer))

    doomed_indices = []
    for residue in pdb.topology.residues():
        if not is_force_field_safe(residue.name):
            doomed_indices.extend(atom.index for atom in residue.atoms())

    modeller = app.Modeller(pdb.topology, pdb.positions)
    doomed = [
        residue for residue in modeller.topology.residues()
        if not is_force_field_safe(residue.name)
    ]
    modeller.delete(doomed)

    # Surviving atoms are in the same order as the original, minus the deleted
    # ones, so the mapping can be rebuilt by walking the kept indices.
    doomed_set = set(doomed_indices)
    kept_indices = [i for i in range(len(positions)) if i not in doomed_set]

    lining = []
    if doomed_indices:
        doomed_xyz = positions[doomed_indices]
        kept_xyz = positions[kept_indices]
        cutoff_nm = POCKET_RESTRAINT_CUTOFF / 10.0
        distances = np.linalg.norm(
            kept_xyz[:, None, :] - doomed_xyz[None, :, :], axis=2).min(axis=1)
        lining = [new for new, d in enumerate(distances) if d <= cutoff_nm]

    return modeller, excluded, lining


def _merge_coords_by_name(pdb, minimised_topology, minimised_positions):
    """
    Write minimised coordinates back into the full structure.

    Matching is by (chain, residue name, residue id, atom name) rather than
    by index: the minimised subset has a different atom count, so positional
    mapping would silently scramble the coordinates.

    Units are stripped to plain Vec3 in nanometres first — mixing bare Vec3
    with unit-carrying Quantity objects in one list produces a sequence that
    OpenMM's PDB writer cannot consume.
    """
    minimised = minimised_positions.value_in_unit(unit.nanometer)
    moved = {
        _atom_key(atom): minimised[i]
        for i, atom in enumerate(minimised_topology.atoms())
    }

    merged = list(pdb.positions.value_in_unit(unit.nanometer))
    for i, atom in enumerate(pdb.topology.atoms()):
        new_position = moved.get(_atom_key(atom))
        if new_position is not None:
            merged[i] = new_position
    return unit.Quantity(merged, unit.nanometer)


def _write(path, topology, positions):
    with open(path, 'w') as fh:
        app.PDBFile.writeFile(topology, positions, fh, keepIds=True)
    normalize_element_column_case(path)


def minimize_structure(input_pdb_path, output_pdb_path,
                       force_field='amber14', use_gbsa=True):
    """
    Minimise a structure, degrading gracefully.

    Always returns a dict containing ``status`` (one of ``full``, ``partial``,
    ``partial_no_implicit_solvent``, ``failed``). Energies are ``None`` when no
    minimisation ran — never the string ``'N/A'``, which is not a number and
    forces callers to type-check.
    """
    result = {
        'force_field': force_field,
        'gbsa_used': use_gbsa,
        'status': STATUS_FAILED,
        'platform': None,
        'energy_before_kJ_mol': None,
        'energy_after_kJ_mol': None,
        'delta_energy_kJ_mol': None,
        'iterations_max': MAX_ITERATIONS,
        'energy_tolerance_kJ_mol_nm': ENERGY_TOLERANCE,
        'energy_decreased': False,
        'converged': False,
        'rms_force_kJ_mol_nm': None,
        'terminals_repaired': False,
        'restrained_atoms': 0,
        'excluded_residues': [],
        'warnings': [],
        'error': None,
        'output_path': output_pdb_path,
    }

    try:
        pdb = app.PDBFile(input_pdb_path)
    except Exception as exc:
        result['error'] = f"Could not read PDB: {exc}"
        shutil.copy2(input_pdb_path, output_pdb_path)
        return result

    if pdb.topology.getNumAtoms() == 0:
        result['error'] = (
            "The structure contains no atoms — every chain may have been "
            "filtered out during cleaning."
        )
        shutil.copy2(input_pdb_path, output_pdb_path)
        return result

    attempts = [
        (STATUS_FULL, use_gbsa, False),
        (STATUS_PARTIAL, use_gbsa, True),
    ]
    if use_gbsa:
        attempts.append((STATUS_PARTIAL_NO_GBSA, False, True))

    errors = []
    for status, gbsa, restrict in attempts:
        try:
            forcefield = _build_forcefield(force_field, gbsa)

            if restrict:
                modeller, excluded, lining = _safe_subset(pdb)
                if modeller is None:
                    # Nothing to strip, so this tier cannot differ from Tier 1.
                    continue
                positions, before, after, sub_topology, repaired, rms_force, platform_name = _run(
                    modeller.topology, modeller.positions, forcefield,
                    restrain_indices=lining
                )
                result['restrained_atoms'] = len(lining) if not repaired else 0
                merged = _merge_coords_by_name(pdb, sub_topology, positions)
                _write(output_pdb_path, pdb.topology, merged)
                result['excluded_residues'] = excluded
                result['warnings'].append(
                    "Minimised only force-field-parameterisable residues. "
                    "These were held at their input coordinates: "
                    + ', '.join(excluded)
                    + f". {len(lining)} atoms lining them were harmonically "
                      "restrained so the site could not relax inward while "
                      "they were absent."
                )
            else:
                positions, before, after, topology, repaired, rms_force, platform_name = _run(
                    pdb.topology, pdb.positions, forcefield
                )
                _write(output_pdb_path, topology, positions)

            if repaired:
                result['terminals_repaired'] = True
                result['warnings'].append(
                    "The force field had no template for the structure as "
                    "given; terminal groups were capped with Modeller before "
                    "minimising."
                )

            if abs(before) > SUSPECT_ENERGY:
                result['warnings'].append(
                    f"Starting energy was {before:.3e} kJ/mol, which indicates "
                    "severe steric clashes in the input. Treat the minimised "
                    "geometry with caution."
                )

            delta = after - before
            result.update({
                'status': status,
                'platform': platform_name,
                'gbsa_used': gbsa,
                'energy_before_kJ_mol': round(before, 1),
                'energy_after_kJ_mol': round(after, 1),
                'delta_energy_kJ_mol': round(delta, 1),
                'energy_decreased': delta < -0.1,
                'rms_force_kJ_mol_nm': round(rms_force, 3),
                # Converged means the minimiser reached the force tolerance it
                # was given, not merely that the energy fell. Capping the
                # iterations stops it early; that is reported as not converged.
                'converged': rms_force <= ENERGY_TOLERANCE,
            })
            return result

        except Exception as exc:
            errors.append(f"[{status}] {exc}")
            logger.warning("Minimisation tier %s failed: %s", status, exc)

    # Every tier failed — pass the input through, but say so plainly.
    shutil.copy2(input_pdb_path, output_pdb_path)
    result['error'] = (
        "Minimisation could not be performed; the structure was passed through "
        "unchanged. Attempts: " + ' | '.join(errors)
    )
    return result
