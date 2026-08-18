import os
import subprocess
import shutil

def export_structure(pdb_path, output_path, format_type):
    """
    Exports a prepared PDB to the target format.

    - autodock / vina: Converts to PDBQT using OpenBabel with Gasteiger charges.
    - gromacs: Provides a suitably named PDB. Note: full GROMACS topology (.top,
      .itp) must be generated separately using 'gmx pdb2gmx'.
    - default: Standard PDB copy.
    """
    format_type = format_type.lower()

    if format_type in ['autodock', 'vina']:
        return convert_to_pdbqt(pdb_path, output_path)
    elif format_type == 'gromacs':
        # GROMACS requires a clean PDB. The actual topology (.top) and .gro file
        # must be generated with: gmx pdb2gmx -f input.pdb -o output.gro -water spce
        # We name it _gromacs.pdb to distinguish it clearly.
        gro_pdb_path = output_path.replace('.pdb', '_gromacs.pdb')
        if pdb_path != gro_pdb_path:
            shutil.copy2(pdb_path, gro_pdb_path)
        return True, gro_pdb_path
    else:
        # Default to PDB
        if pdb_path != output_path:
            shutil.copy2(pdb_path, output_path)
        return True, output_path


def convert_to_pdbqt(pdb_path, output_path):
    """
    Uses OpenBabel to convert PDB to PDBQT with Gasteiger partial charges.
    """
    # Robustly handle extension
    base = os.path.splitext(output_path)[0]
    pdbqt_path = base + '.pdbqt'
    try:
        cmd = [
            'obabel', pdb_path,
            '-O', pdbqt_path,
            '-xh',                           # merge non-polar H onto heavy atoms
            '--partialcharge', 'gasteiger',  # FIX: compute real partial charges
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0:
            return True, pdbqt_path
        else:
            return False, f"OpenBabel error: {result.stderr.strip() or 'unknown error'}"
    except FileNotFoundError:
        return False, "OpenBabel (obabel) is not installed or not on PATH. Install it from https://openbabel.org to enable PDBQT export."
    except subprocess.TimeoutExpired:
        return False, "OpenBabel conversion timed out after 120 seconds. The structure may be too large."
    except Exception as e:
        return False, f"Conversion failed: {str(e)}"
