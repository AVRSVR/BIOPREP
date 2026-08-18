from pdbfixer import PDBFixer
from openmm.app import PDBFile

def add_hydrogens(input_pdb_path, output_pdb_path, ph=7.4, reconstruct_loops=False, add_missing_atoms=False):
    """
    Uses PDBFixer to add missing hydrogens to a PDB file at a specific pH.
    
    Args:
        input_pdb_path (str): Path to the temporary cleaned PDB.
        output_pdb_path (str): Desired output path for the final protonated PDB.
        ph (float): Target pH for physiological protonation.
        reconstruct_loops (bool): Whether to reconstruct missing loops.
        add_missing_atoms (bool): Whether to add missing heavy atoms.
    """
    fixer = PDBFixer(filename=input_pdb_path)
    
    # 1. Identify and add missing residues (if any loop domains are abruptly broken)
    if reconstruct_loops:
        fixer.findMissingResidues()
    else:
        # We still need to call it but we clear it so it doesn't add them
        fixer.findMissingResidues()
        fixer.missingResidues = {}
    
    # 2. Convert non-standard residues to standard ones where possible
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    
    # FIX: REMOVED removeHeterogens(False). 
    # This was the 'predator' stripping out our ligands/lipids.
    # Since we already cleaned the PDB surgically via Biopython, 
    # we trust the input structure is exactly what we want.
    
    # 3. Fill in missing heavy atoms
    if add_missing_atoms or reconstruct_loops:
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
    
    # 4. Add missing hydrogens appropriate for the pH
    # We wrap this in a try-block because OpenMM can occasionally struggle 
    # with extremely non-standard ligands; we want to preserve the residue 
    # even if we can't perfectly protonate it.
    try:
        fixer.addMissingHydrogens(ph)
    except Exception as e:
        print(f"Warning: Could not perfectly protonate all residues: {e}")
    
    # Write topology and positions to the output file
    with open(output_pdb_path, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
