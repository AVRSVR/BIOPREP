"""
Export a prepared structure into a docking-ready format.
"""

import os
import shutil
import subprocess

OBABEL_TIMEOUT = 120


def export_structure(pdb_path, output_path, format_type):
    """
    Convert a prepared PDB to the requested format.

    Returns ``(True, path)`` on success or ``(False, message)`` on failure.

    - ``autodock`` / ``vina``: PDBQT with Gasteiger charges, written as a rigid
      receptor.
    - ``gromacs``: a suitably named PDB. The topology still has to come from
      ``gmx pdb2gmx``.
    - anything else: a plain PDB copy.
    """
    format_type = (format_type or '').lower()
    base, _ = os.path.splitext(output_path)

    if format_type in ('autodock', 'vina'):
        return convert_to_pdbqt(pdb_path, base + '.pdbqt')

    if format_type == 'gromacs':
        # splitext, not str.replace: replace() substitutes every occurrence of
        # '.pdb' in the path, which mangles any directory containing that text.
        destination = base + '_gromacs.pdb'
        if os.path.abspath(pdb_path) != os.path.abspath(destination):
            shutil.copy2(pdb_path, destination)
        return True, destination

    destination = base + '.pdb'
    if os.path.abspath(pdb_path) != os.path.abspath(destination):
        shutil.copy2(pdb_path, destination)
    return True, destination


def convert_to_pdbqt(pdb_path, pdbqt_path):
    """
    Convert a receptor PDB to PDBQT using OpenBabel.

    ``-xr`` marks the output as a rigid receptor. Without it OpenBabel writes a
    ligand-style torsion tree (ROOT/BRANCH records), which AutoDock and Vina
    will not accept as a receptor.
    """
    command = [
        'obabel', pdb_path,
        '-O', pdbqt_path,
        '-xr',                            # rigid receptor, no torsion tree
        '-xh',                            # merge non-polar hydrogens
        '--partialcharge', 'gasteiger',
    ]

    try:
        result = subprocess.run(command, capture_output=True, text=True,
                                timeout=OBABEL_TIMEOUT)
    except FileNotFoundError:
        return False, ("OpenBabel (obabel) is not installed or not on PATH. "
                       "Install it from https://openbabel.org to enable PDBQT export.")
    except subprocess.TimeoutExpired:
        return False, (f"OpenBabel timed out after {OBABEL_TIMEOUT} seconds. "
                       "The structure may be too large.")
    except OSError as exc:
        return False, f"Could not run OpenBabel: {exc}"

    if result.returncode != 0:
        return False, f"OpenBabel error: {result.stderr.strip() or 'unknown error'}"

    # OpenBabel can exit 0 having written nothing useful.
    if not os.path.exists(pdbqt_path) or os.path.getsize(pdbqt_path) == 0:
        return False, (f"OpenBabel produced no output. "
                       f"{result.stderr.strip() or 'No error message was reported.'}")

    return True, pdbqt_path
