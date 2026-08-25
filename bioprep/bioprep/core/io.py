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


def normalize_element_column_case(pdb_path):
    """
    Upper-case the element column (77-78) of every ATOM/HETATM line in place.

    RCSB and this project's own SpecCompliantPDBIO write 'ZN'; OpenMM's
    PDBFile.writeFile - used for every post-protonation and post-minimisation
    write - writes 'Zn'. Both identify the correct element to every tool this
    project checked (OpenMM, OpenBabel), so this was previously left alone
    deliberately (see BACKEND_AUDIT.md); it is fixed here only because doing
    so is now nearly free, not because anything downstream actually needed it.
    """
    with open(pdb_path, 'r') as fh:
        lines = fh.readlines()

    changed = False
    for i, line in enumerate(lines):
        if line[:6] in ('ATOM  ', 'HETATM') and len(line) >= 78:
            upper = line[76:78].upper()
            if upper != line[76:78]:
                lines[i] = line[:76] + upper + line[78:]
                changed = True

    if changed:
        with open(pdb_path, 'w') as fh:
            fh.writelines(lines)


class SpecCompliantPDBIO(PDBIO):
    """
    PDBIO that puts the element symbol where the PDB specification puts it.

    Biopython's ``_ATOM_FORMAT_STRING`` places five spaces between the
    B-factor and the segment identifier where the format requires six, so
    every field from column 67 onward is written one column early. The
    element then straddles the boundary: a reader taking columns 77-78 sees
    the symbol's second character followed by a blank.

    Single-character elements survive by luck. Two-character metals do not,
    and the failure is silent whenever the leftover character is itself a
    valid symbol - ZN reads back as nitrogen, FE as fluorine. Both Biopython
    and OpenMM reproduce that, which matters for any structure containing a
    metal cofactor.
    """

    def _get_atom_line(self, atom, hetfield, segid, atom_number, resname,
                       resseq, icode, chain_id, charge="  "):
        line = super()._get_atom_line(atom, hetfield, segid, atom_number,
                                      resname, resseq, icode, chain_id, charge)

        # Only coordinate records carry an element column.
        if not line.startswith(("ATOM  ", "HETATM")):
            return line

        # Columns 1-66 (x, y, z, occupancy, B-factor) already match the spec;
        # rebuild everything after them.
        element = (atom.element or "").strip().upper()
        return (
            f"{line[:66]}"
            f"{'':6}"        # columns 67-72, blank
            f"{segid:<4}"    # columns 73-76, segment identifier
            f"{element:>2}"  # columns 77-78, element, right-justified
            f"{charge:>2}"   # columns 79-80, formal charge
            "\n"
        )


def _atom_key_from_line(line):
    """
    Identify an atom by its PDB fields rather than its serial number.

    Serials are reassigned on every write, so they cannot be used to match an
    atom across files. Chain, residue and atom name survive.
    """
    return (
        line[21],              # chain id
        line[22:26].strip(),   # residue sequence number
        line[26].strip(),      # insertion code
        line[17:20].strip(),   # residue name
        line[12:16].strip(),   # atom name
        line[16].strip(),      # altloc
    )


def _serial_index(pdb_path):
    """Return (serial -> key, key -> serial) for one PDB file."""
    by_serial, by_key = {}, {}
    try:
        with open(pdb_path, 'r') as fh:
            for line in fh:
                if line[:6] not in ('ATOM  ', 'HETATM'):
                    continue
                try:
                    serial = int(line[6:11])
                except ValueError:
                    continue
                key = _atom_key_from_line(line)
                by_serial[serial] = key
                by_key.setdefault(key, serial)
    except OSError:
        pass
    return by_serial, by_key


def conect_serials(line):
    """Yield the atom serials referenced by one CONECT record."""
    body = line[6:].rstrip('\n')
    for start in range(0, len(body), 5):
        chunk = body[start:start + 5].strip()
        if chunk:
            try:
                yield int(chunk)
            except ValueError:
                continue


def transfer_conect_records(source_path, output_path):
    """
    Copy CONECT records from ``source_path`` onto ``output_path``.

    Biopython drops CONECT entirely when parsing, so a saved structure has no
    bond records at all. Ligand connectivity would otherwise be lost the first
    time a structure was written, before anything downstream could preserve it.

    Bonds are matched by atom identity, because both files number their atoms
    independently. A bond is carried over only when both of its atoms survived
    into the output, so filtering out a ligand also drops its bonds.

    Returns the number of bonds written.
    """
    source_by_serial, _ = _serial_index(source_path)
    _, output_by_key = _serial_index(output_path)
    if not source_by_serial or not output_by_key:
        return 0

    bonds = set()
    try:
        with open(source_path, 'r') as fh:
            for line in fh:
                if not line.startswith('CONECT'):
                    continue
                serials = list(conect_serials(line))
                if len(serials) < 2:
                    continue

                central = output_by_key.get(source_by_serial.get(serials[0]))
                if central is None:
                    continue
                for partner_serial in serials[1:]:
                    partner = output_by_key.get(
                        source_by_serial.get(partner_serial))
                    if partner is not None and partner != central:
                        bonds.add((min(central, partner),
                                   max(central, partner)))
    except OSError:
        return 0

    if not bonds:
        return 0

    partners = {}
    for first, second in bonds:
        partners.setdefault(first, set()).add(second)
        partners.setdefault(second, set()).add(first)

    records = []
    for central in sorted(partners):
        listed = sorted(partners[central])
        # A CONECT record holds at most four partners; spill into extra records.
        for start in range(0, len(listed), 4):
            chunk = listed[start:start + 4]
            records.append('CONECT' + f'{central:>5}'
                           + ''.join(f'{p:>5}' for p in chunk) + '\n')

    with open(output_path, 'r') as fh:
        existing = [l for l in fh if not l.startswith(('CONECT', 'END'))]

    with open(output_path, 'w') as fh:
        fh.writelines(existing)
        fh.writelines(records)
        fh.write('END\n')

    return len(bonds)


def transfer_header_records(source_path, output_path, records=('SEQRES',)):
    """
    Copy header records from ``source_path`` onto ``output_path``.

    Biopython drops these on parse and PDBIO writes none, so a saved structure
    carries no SEQRES. PDBFixer compares SEQRES against the observed residues
    to work out what is missing, so without it loop reconstruction has nothing
    to rebuild and silently does nothing.

    SEQRES is filtered to the chains that survived into the output. Carrying
    the sequence of a chain that was filtered out would make PDBFixer treat
    that whole chain as missing and try to build it from nothing.

    Returns the number of records written.
    """
    chains_present = set()
    try:
        with open(output_path, 'r') as fh:
            for line in fh:
                if line[:6] in ('ATOM  ', 'HETATM') and len(line) > 21:
                    chains_present.add(line[21])
    except OSError:
        return 0
    if not chains_present:
        return 0

    carried = []
    try:
        with open(source_path, 'r') as fh:
            for line in fh:
                name = line[:6].strip()
                if name not in records:
                    continue
                if name == 'SEQRES':
                    # Column 12 (index 11) holds the chain identifier.
                    if len(line) > 11 and line[11] in chains_present:
                        carried.append(line)
                else:
                    carried.append(line)
    except OSError:
        return 0

    if not carried:
        return 0

    with open(output_path, 'r') as fh:
        body = fh.readlines()

    insert_at = len(body)
    for index, line in enumerate(body):
        if line[:6] in ('ATOM  ', 'HETATM', 'MODEL '):
            insert_at = index
            break

    with open(output_path, 'w') as fh:
        fh.writelines(body[:insert_at])
        fh.writelines(carried)
        fh.writelines(body[insert_at:])

    return len(carried)


def save_pdb(structure, output_path, select=None, source_pdb=None,
             conect_source=None):
    """
    Save a Biopython Structure to a PDB file.

    Args:
        structure (Structure): Biopython Structure object.
        output_path (str): Path to write the PDB file.
        select (Select, optional): Biopython Select used to filter atoms.
        source_pdb (str, optional): Path of the file the structure was read
            from. Records Biopython discards on parse are carried onto the
            output: CONECT bonds remapped to the new atom numbering, and
            SEQRES for the chains that survived. Without it, ligand bonds and
            the reference sequence are both lost.
        conect_source (str, optional): Deprecated alias for ``source_pdb``.
    """
    parent = os.path.dirname(os.path.abspath(output_path))
    if parent:
        os.makedirs(parent, exist_ok=True)

    io = SpecCompliantPDBIO()
    io.set_structure(structure)
    if select:
        io.save(output_path, select)
    else:
        io.save(output_path)

    origin = source_pdb or conect_source
    if origin and os.path.exists(origin):
        transfer_conect_records(origin, output_path)
        transfer_header_records(origin, output_path)
