"""
pH-dependent protonation that does not let PDBFixer touch ligands.

What the pH argument actually buys, measured on HIV-1 protease (1HSG):

  LYS   three HZ below pH 11, two above          correct, pKa about 10.5
  HIS   HIP at pH 2, neutral HID from 7.4 up     correct, pKa about 6
  ASP   protonated below about pH 5              correct
  GLU   protonated below about pH 5              correct
  TYR   keeps its HH at every pH, up to 13       NOT handled

Two limits worth knowing about before trusting a protonation state:

* Tyrosine is never deprotonated. OpenMM's Modeller has no tyrosinate
  variant, so above pKa 10.1 the hydroxyl stays put.
* The pKa values are model-compound values applied uniformly. There is no
  structure-specific pKa prediction here - nothing like PROPKA - so a buried
  or salt-bridged residue whose pKa is shifted by several units will be
  assigned the wrong state at physiological pH. That matters most for a
  catalytic residue, which is exactly where the shift tends to be largest.

The histidine tautomer, HID or HIE, is chosen by OpenMM from the local
hydrogen-bonding environment rather than fixed; on 1HSG both histidines come
out HID.

PDBFixer has no templates for arbitrary small molecules. Handed a
protein-ligand complex it will try to template-match the ligand and can
strip atoms or add nonsensical hydrogens to it. The fix is to hold the
ligand records out of the file entirely, protonate the biopolymer, then
splice the untouched ligand coordinates back in.
"""

import os
import tempfile

from pdbfixer import PDBFixer
from openmm.app import PDBFile

from .residues import is_standard, is_water

# A residue can be missing more than its sidechain - real deposited structures
# routinely have a chain terminus with no density past the amide nitrogen, so
# CA/C/O are absent too. Those aren't a modelling choice the way a sidechain
# is: without them PDBFixer has nothing to compute a terminal OXT position
# from, and addMissingAtoms() raises KeyError trying to read them. So only
# sidechain atoms get dropped when add_missing_atoms is off; backbone atoms
# stay regardless, same reasoning as why missingTerminals is never cleared.
BACKBONE_ATOM_NAMES = frozenset(('N', 'CA', 'C', 'O'))


def _split_records(pdb_path):
    """
    Partition a PDB into the part PDBFixer may rewrite and the part it may not.

    Waters stay with the biopolymer: PDBFixer protonates them correctly and
    the force fields have templates for them. Everything else flagged HETATM
    (drugs, cofactors, lipids, unparameterised ions) is held out verbatim.

    Only the first model is read. The MODEL/ENDMDL markers are dropped on
    write, so keeping later models would weld an NMR ensemble into a single
    chimeric protein with every residue repeated once per model.

    Waters are returned separately from the polymer so a TER can be written
    between them. Crystallographic waters usually carry the same chain id as
    the protein; with no TER, PDBFixer reads them as a continuation of that
    chain, the last amino acid is no longer the chain's end, and no terminal
    OXT is added. The structure is then unparameterisable and every
    minimisation tier fails - which is what happened to 1STP whenever
    structural waters were kept.

    Returns (polymer_lines, water_lines, ligand_lines, conect_lines).
    """
    biopolymer, waters, ligand, conect = [], [], [], []
    ligand_serials = set()
    past_first_model = False

    with open(pdb_path, 'r') as fh:
        for line in fh:
            record = line[:6]

            if past_first_model:
                continue

            if record in ('ATOM  ', 'HETATM'):
                resname = line[17:20].strip()
                if is_water(resname):
                    waters.append(line)
                elif record == 'ATOM  ' or is_standard(resname):
                    biopolymer.append(line)
                else:
                    ligand.append(line)
                    try:
                        ligand_serials.add(int(line[6:11]))
                    except ValueError:
                        pass

            elif record == 'CONECT':
                conect.append(line)

            elif record.strip() == 'ENDMDL':
                # Everything past here belongs to a later model.
                past_first_model = True

            elif record.strip() in ('TER', 'END', 'MODEL'):
                # Structural records are regenerated on write; drop them.
                # Compared after stripping: a real file writes 'END' as three
                # characters, which never matched the padded 'END   ' this
                # once tested for. The record then fell through into the
                # biopolymer list, and once waters were split out it landed
                # ahead of them - so the file handed to PDBFixer ended before
                # its own waters and every one of them was silently dropped.
                continue
            else:
                # SEQRES, HELIX, CRYST1 etc. must reach PDBFixer — it uses
                # SEQRES to work out which residues are missing.
                biopolymer.append(line)

    # Keep only the CONECT records that describe held-out ligands.
    ligand_conect = [
        c for c in conect
        if any(s in ligand_serials for s in _conect_serials(c))
    ]
    return biopolymer, waters, ligand, ligand_conect


def _conect_serials(line):
    """Yield the atom serials referenced by a CONECT record."""
    body = line[6:].rstrip('\n')
    for start in range(0, len(body), 5):
        chunk = body[start:start + 5].strip()
        if chunk:
            try:
                yield int(chunk)
            except ValueError:
                continue


def _merge(protein_pdb_path, ligand_lines, ligand_conect, output_path):
    """
    Append the untouched ligand records after the protonated biopolymer,
    renumbering atom serials so they continue on from the protein and
    rewriting CONECT records to match.
    """
    protein_lines = []
    protein_conect = []
    max_serial = 0
    with open(protein_pdb_path, 'r') as fh:
        for line in fh:
            if line.startswith('END'):
                continue
            if line[:6] == 'CONECT':
                # Held back so every CONECT lands after the last atom record,
                # which is where the PDB spec puts them.
                protein_conect.append(line)
                continue
            if line[:6] in ('ATOM  ', 'HETATM'):
                try:
                    max_serial = max(max_serial, int(line[6:11]))
                except ValueError:
                    pass
            protein_lines.append(line)

    serial_map = {}
    renumbered = []
    next_serial = max_serial + 1
    for line in ligand_lines:
        try:
            old = int(line[6:11])
            serial_map[old] = next_serial
        except ValueError:
            pass
        renumbered.append(f"{line[:6]}{next_serial:>5}{line[11:].rstrip()}\n")
        next_serial += 1

    remapped_conect = []
    for line in ligand_conect:
        serials = [serial_map.get(s) for s in _conect_serials(line)]
        if not serials or serials[0] is None:
            continue
        # Drop references to atoms that were not carried over.
        serials = [s for s in serials if s is not None]
        remapped_conect.append('CONECT' + ''.join(f"{s:>5}" for s in serials) + '\n')

    with open(output_path, 'w') as out:
        out.writelines(protein_lines)
        if renumbered:
            out.write('TER\n')
            out.writelines(renumbered)
        out.writelines(protein_conect)
        out.writelines(remapped_conect)
        out.write('END\n')


def add_hydrogens(input_pdb_path, output_pdb_path, ph=7.4,
                  reconstruct_loops=False, add_missing_atoms=False):
    """
    Add hydrogens at the given pH, protecting ligands from PDBFixer.

    Returns a dict describing what actually happened, so callers can report
    the truth instead of assuming success::

        {
          'hydrogens_added': bool,
          'ph': float,
          'ligands_preserved': ['BTN', ...],
          'nonstandard_replaced': ['MSE', ...],
          'loops_reconstructed': int,
          'warnings': [str, ...],
        }
    """
    result = {
        'hydrogens_added': False,
        'ph': ph,
        'ligands_preserved': [],
        'ligand_status': [],
        'nonstandard_replaced': [],
        'terminals_repaired': 0,
        'loops_reconstructed': 0,
        'warnings': [],
    }

    biopolymer, water_lines, ligand_lines, ligand_conect = _split_records(
        input_pdb_path)

    # Record not just which ligands were held back but why, so the preparation
    # report can explain an otherwise surprising result: the ligand comes out
    # with no hydrogens on it, and that is deliberate.
    atoms_per_ligand = {}
    for line in ligand_lines:
        atoms_per_ligand[line[17:20].strip()] = \
            atoms_per_ligand.get(line[17:20].strip(), 0) + 1

    result['ligands_preserved'] = sorted(atoms_per_ligand)
    result['ligand_status'] = [
        {
            'residue': name,
            'atoms': atoms_per_ligand[name],
            'action': 'preserved',
            'reason': ('Held out of PDBFixer, which has no template for it. '
                       'Blind protonation would alter the ligand geometry, so '
                       'its coordinates are carried through unchanged and no '
                       'hydrogens are added to it.'),
        }
        for name in sorted(atoms_per_ligand)
    ]

    if not any(l[:6] in ('ATOM  ', 'HETATM') for l in biopolymer):
        raise ValueError(
            "No biopolymer atoms found to protonate. The structure may consist "
            "entirely of ligands, or every chain was filtered out during cleaning."
        )

    fd, protein_only = tempfile.mkstemp(suffix='_protein.pdb')
    os.close(fd)
    fd, protonated = tempfile.mkstemp(suffix='_protonated.pdb')
    os.close(fd)

    try:
        with open(protein_only, 'w') as fh:
            fh.writelines(biopolymer)
            if water_lines:
                # A TER, or PDBFixer reads the waters as a continuation of the
                # protein chain: the last amino acid stops being the chain's
                # end, no terminal OXT is added, and the structure will not
                # parameterise. That silently broke minimisation for any
                # structure processed with structural waters kept.
                fh.write('TER\n')
                fh.writelines(water_lines)
            fh.write('END\n')

        # PDBFixer(filename=...) closes its handle only on the success path —
        # a parse failure leaves the file open, and on Windows the open handle
        # makes the cleanup below fail with WinError 32, masking the real
        # error. Owning the handle here guarantees it is released either way.
        with open(protein_only, 'r') as fh:
            fixer = PDBFixer(pdbfile=fh)

        fixer.findMissingResidues()
        if reconstruct_loops:
            result['loops_reconstructed'] = sum(
                len(v) for v in fixer.missingResidues.values()
            )
        else:
            fixer.missingResidues = {}

        fixer.findNonstandardResidues()
        if fixer.nonstandardResidues:
            result['nonstandard_replaced'] = sorted(
                {res.name for res, _ in fixer.nonstandardResidues}
            )
            fixer.replaceNonstandardResidues()

        fixer.findMissingAtoms()

        # A terminal OXT is not optional modelling. Without it the last residue
        # matches no force-field template, so createSystem fails and the whole
        # structure falls through every minimisation tier - which is what
        # happens to any crystal structure deposited without a capped terminus.
        # Rebuilding absent side chains is a modelling choice and stays behind
        # the flag, so when the flag is off, sidechain atoms are dropped from
        # the list - but backbone atoms (N/CA/C/O) stay: a residue with no
        # resolved density past its amide nitrogen is missing those too, and
        # PDBFixer's own terminal-placement math needs CA and O to exist.
        terminal_atoms = sum(len(v) for v in fixer.missingTerminals.values())
        if not (add_missing_atoms or reconstruct_loops):
            fixer.missingAtoms = {
                residue: [atom for atom in atoms if atom.name in BACKBONE_ATOM_NAMES]
                for residue, atoms in fixer.missingAtoms.items()
            }
            fixer.missingAtoms = {
                residue: atoms for residue, atoms in fixer.missingAtoms.items() if atoms
            }

        result['terminals_repaired'] = terminal_atoms
        if terminal_atoms:
            result['warnings'].append(
                f"Added {terminal_atoms} missing terminal atom(s); the input "
                "chain ended without one and would not have parameterised."
            )

        fixer.addMissingAtoms()

        try:
            fixer.addMissingHydrogens(ph)
            result['hydrogens_added'] = True
        except Exception as exc:
            result['warnings'].append(f"Hydrogen addition failed: {exc}")

        with open(protonated, 'w') as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

        _merge(protonated, ligand_lines, ligand_conect, output_pdb_path)

    finally:
        # Never let cleanup raise: an error here would replace whatever went
        # wrong above, and the caller would see a temp-file path instead of
        # the real diagnostic.
        for path in (protein_only, protonated):
            try:
                os.remove(path)
            except OSError:
                pass

    return result
