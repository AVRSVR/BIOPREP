"""
Structure loading and saving.

PDB and mmCIF are both accepted on input. mmCIF is converted to PDB at the
edge of the pipeline rather than carried through it: PDBFixer and OpenMM both
consume PDB, so one conversion here is far less risky than teaching every
downstream module a second format.
"""

import os
import tempfile

from Bio.PDB import PDBParser, PDBIO, MMCIFParser

MMCIF_EXTENSIONS = {'.cif', '.mmcif'}
PDB_EXTENSIONS = {'.pdb', '.ent'}

# Hard limits of the fixed-column PDB format.
MAX_PDB_ATOM_SERIAL = 99999
MAX_PDB_CHAIN_ID_LENGTH = 1


def structure_id_for(file_path):
    """
    Derive a Biopython structure id from a path.

    splitext drops only the final suffix, so '1abc.v2.pdb' keeps its '1abc.v2'
    stem instead of being cut at the first dot.
    """
    filename = os.path.basename(file_path.replace("\\", os.sep))
    return os.path.splitext(filename)[0] or filename


def _looks_like_mmcif(file_path):
    """
    Sniff the file contents for mmCIF markers.

    Extension alone is not enough: RCSB serves mmCIF by default now, and a
    downloaded .cif renamed to .pdb would otherwise fail with a confusing
    parse error rather than simply being read.
    """
    try:
        with open(file_path, 'r', errors='ignore') as fh:
            for _ in range(200):
                line = fh.readline()
                if not line:
                    break
                stripped = line.strip()
                if (stripped.startswith('data_')
                        or stripped.startswith('_atom_site.')
                        or stripped.startswith('loop_')):
                    return True
                if stripped.startswith(('ATOM  ', 'HETATM', 'HEADER', 'CRYST1',
                                        'MODEL ', 'SEQRES', 'EXPDTA')):
                    return False
    except OSError:
        pass
    return False


def detect_format(file_path):
    """Return 'mmcif' or 'pdb' for a structure file, preferring content."""
    extension = os.path.splitext(file_path)[1].lower()
    if extension in MMCIF_EXTENSIONS:
        return 'mmcif'
    return 'mmcif' if _looks_like_mmcif(file_path) else 'pdb'


def load_structure(file_path):
    """
    Load a PDB or mmCIF file and return the Biopython Structure.

    Raises ValueError if the file contains no atoms.
    """
    structure_id = structure_id_for(file_path)
    file_format = detect_format(file_path)

    parser = (MMCIFParser(QUIET=True) if file_format == 'mmcif'
              else PDBParser(QUIET=True))

    try:
        structure = parser.get_structure(structure_id, file_path)
    except ValueError:
        # Biopython already raises a clear ValueError for things like an
        # empty file; let those through unchanged.
        raise
    except Exception as exc:
        # A truncated or malformed file surfaces as whatever the parser
        # happened to trip over - MMCIFParser raises KeyError('_atom_site.id')
        # for a mmCIF with no coordinate table. Translate it into something
        # the caller can show a user.
        raise ValueError(
            f"'{os.path.basename(file_path)}' could not be parsed as "
            f"{'mmCIF' if file_format == 'mmcif' else 'PDB'}: "
            f"{type(exc).__name__}: {exc}. The file may be truncated or "
            "malformed."
        ) from exc

    # Both parsers accept arbitrary text and return an empty Structure when
    # they find no coordinate records, so a mis-named or non-structure upload
    # would sail through here and only surface much later as a misleading
    # "cleaning removed every atom" error. Reject it where the cause is known.
    if not any(True for _ in structure.get_atoms()):
        raise ValueError(
            f"No atoms could be read from '{os.path.basename(file_path)}'. "
            "The file does not contain PDB ATOM/HETATM records or an mmCIF "
            "atom_site table - it may be empty, truncated, or in another format."
        )

    return structure


def load_pdb(file_path):
    """
    Backwards-compatible alias for :func:`load_structure`.

    Kept because callers throughout the codebase use this name; it accepts
    mmCIF as well, despite what the name suggests.
    """
    return load_structure(file_path)


def check_pdb_representable(structure):
    """
    Return a list of reasons the structure cannot be written as legacy PDB.

    The format has fixed columns: five digits of atom serial and a single
    character of chain id. Large assemblies routinely exceed both, and writing
    anyway produces a silently corrupt file.
    """
    problems = []

    atom_count = sum(1 for _ in structure.get_atoms())
    if atom_count > MAX_PDB_ATOM_SERIAL:
        problems.append(
            f"it has {atom_count} atoms, over the PDB limit of "
            f"{MAX_PDB_ATOM_SERIAL}"
        )

    long_chains = sorted({
        str(chain.id) for model in structure for chain in model
        if len(str(chain.id).strip()) > MAX_PDB_CHAIN_ID_LENGTH
    })
    if long_chains:
        problems.append(
            "these chain identifiers are longer than the single character PDB "
            f"allows: {', '.join(long_chains)}"
        )

    return problems


def convert_to_pdb(input_path, output_path):
    """
    Write any supported structure file out as PDB.

    Raises ValueError with a specific reason when the structure is too large
    or too complex for the legacy format.
    """
    structure = load_structure(input_path)

    problems = check_pdb_representable(structure)
    if problems:
        raise ValueError(
            f"'{os.path.basename(input_path)}' cannot be converted to PDB "
            "format because " + "; and ".join(problems) + ". Prepare this "
            "structure as individual chains or a smaller sub-assembly."
        )

    save_pdb(structure, output_path)
    return output_path


def ensure_pdb(input_path, workdir=None):
    """
    Return a path to a PDB rendering of ``input_path``.

    PDB input passes straight through. mmCIF is converted into ``workdir``
    (or a temporary directory).

    Returns ``(pdb_path, converted_from)``, where ``converted_from`` is the
    original format name when a conversion happened and None otherwise.
    """
    if detect_format(input_path) == 'pdb':
        return input_path, None

    if workdir is None:
        workdir = tempfile.mkdtemp(prefix='bioprep_cif_')

    output_path = os.path.join(
        workdir, f'{structure_id_for(input_path)}__from_mmcif.pdb')
    convert_to_pdb(input_path, output_path)
    return output_path, 'mmcif'


def save_pdb(structure, output_path, select=None):
    """
    Save a Biopython Structure to a PDB file.

    Args:
        structure (Structure): Biopython Structure object.
        output_path (str): Path to write the PDB file.
        select (Select, optional): Biopython Select used to filter atoms.
    """
    io = PDBIO()
    io.set_structure(structure)
    if select:
        io.save(output_path, select)
    else:
        io.save(output_path)
