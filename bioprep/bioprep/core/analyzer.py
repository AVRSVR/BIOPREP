"""
Structural metadata extraction.

Counts come from the first model only. Iterating every model silently
multiplied atom counts, water counts and detected gaps by the model count,
so a 20-model NMR ensemble reported twenty times the real numbers.
"""

from .residues import is_water


def _first_model(structure):
    """Return the first model, or None for an empty structure."""
    for model in structure:
        return model
    return None


def analyze_structure(structure):
    """
    Extract metadata for the UI from the first model of a structure.

    Returns chain ids, water and heteroatom inventory, atom counts broken down
    by category, sequence gaps, and the number of models present.
    """
    chains = []
    water_count = 0
    heteroatoms = set()
    protein_atoms = 0
    water_atoms = 0
    hetero_atoms = 0
    sequence_gaps = []

    model_count = sum(1 for _ in structure)
    model = _first_model(structure)

    if model is not None:
        for chain in model:
            chains.append(chain.id)
            prev_resseq = None

            for residue in chain:
                hetfield = residue.id[0]
                resseq = residue.id[1]
                atom_count = len(residue.get_list())

                if hetfield == 'W' or is_water(residue.resname):
                    water_count += 1
                    water_atoms += atom_count
                    continue

                if hetfield != ' ':
                    heteroatoms.add(residue.resname.strip())
                    hetero_atoms += atom_count
                    continue

                protein_atoms += atom_count

                # Gap detection: a jump in residue numbering between two
                # consecutive observed standard residues.
                if prev_resseq is not None and (resseq - prev_resseq) > 1:
                    sequence_gaps.append({
                        "chain": chain.id,
                        "from": prev_resseq,
                        "to": resseq,
                        "missing_count": resseq - prev_resseq - 1,
                    })
                prev_resseq = resseq

    atoms_total = protein_atoms + water_atoms + hetero_atoms

    return {
        'chains': sorted(set(chains)),
        'water_count': water_count,
        'heteroatoms': sorted(heteroatoms),
        # Whole-structure count for the first model, all categories included.
        'atoms_total': atoms_total,
        'atoms_before': atoms_total,
        'atom_breakdown': {
            'protein': protein_atoms,
            'water': water_atoms,
            'heteroatom': hetero_atoms,
        },
        'sequence_gaps': sequence_gaps,
        'model_count': model_count,
    }


def detect_missing_residues(pdb_path):
    """
    Report residues present in SEQRES but absent from the coordinates.

    PDBFixer keys ``missingResidues`` by ``(chain_index, insertion_position)``.
    Reading that tuple as ``(model, chain_id)`` — as this once did — labelled
    every gap with a residue offset in place of the chain.
    """
    try:
        from pdbfixer import PDBFixer

        fixer = PDBFixer(filename=pdb_path)
        fixer.findMissingResidues()

        chains = list(fixer.topology.chains())
        missing = []
        for (chain_index, insert_at), residue_names in fixer.missingResidues.items():
            try:
                chain_id = chains[chain_index].id
            except (IndexError, AttributeError):
                chain_id = str(chain_index)

            for offset, name in enumerate(residue_names):
                missing.append({
                    "chain": chain_id,
                    "residue": name,
                    "position": insert_at + offset,
                })
        return missing
    except Exception:
        return []
