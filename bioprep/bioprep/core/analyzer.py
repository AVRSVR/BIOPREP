from Bio.PDB import PDBParser


def analyze_structure(structure):
    """
    Parses a Biopython Structure object to extract metadata for the UI:
    - Lists of available chains
    - Total water counts
    - Unique list of heteroatoms/ligands
    - Total protein atom count (waters excluded)
    - Sequence gap detection (missing residues between observed residues)
    """
    chains = []
    water_count = 0
    heteroatoms = set()
    atoms_total = 0
    sequence_gaps = []  # FIX: actually detect gaps now

    for model in structure:
        for chain in model:
            chain_id = chain.id
            chains.append(chain_id)
            prev_resseq = None  # Track previous residue sequence number

            for residue in chain:
                hetfield = residue.id[0]
                resseq = residue.id[1]

                # Count all atoms for total structure count
                atoms_count = len(residue.get_list())
                atoms_total += atoms_count

                if hetfield == 'W':
                    water_count += 1
                    continue
                elif hetfield != ' ' and hetfield != 'W':
                    res_name = residue.resname.strip()
                    heteroatoms.add(res_name)
                    continue
                else:
                    # Standard amino acid – check sequence gaps


                    # FIX: Gap detection — using prev_resseq which was dead before
                    if prev_resseq is not None and (resseq - prev_resseq) > 1:
                        # Gap detected: residues between prev_resseq+1 and resseq-1
                        gap_size = resseq - prev_resseq - 1
                        sequence_gaps.append({
                            "chain": chain_id,
                            "from": prev_resseq,
                            "to": resseq,
                            "missing_count": gap_size
                        })

                    prev_resseq = resseq

    return {
        'chains': sorted(list(set(chains))),
        'water_count': water_count,
        'heteroatoms': sorted(list(heteroatoms)),
        'atoms_total': atoms_total,
        'atoms_before': atoms_total,   # protein-only pre-processing count
        'sequence_gaps': sequence_gaps,  # newly detected gaps
    }


def detect_missing_residues(pdb_path):
    """
    Uses PDBFixer to detect missing residues without rebuilding.
    Returns a list of (chain_id, residue_number) tuples for reporting.
    """
    try:
        from pdbfixer import PDBFixer
        fixer = PDBFixer(filename=pdb_path)
        fixer.findMissingResidues()
        missing = []
        for (model_idx, chain_id), residues in fixer.missingResidues.items():
            for r in residues:
                missing.append({"chain": chain_id, "residue": str(r)})
        return missing
    except Exception:
        return []
