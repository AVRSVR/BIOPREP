import os

from Bio.PDB import PDBParser, PDBIO

def load_pdb(file_path):
    """
    Loads a PDB file and returns the Biopython Structure object.
    
    Args:
        file_path (str): Path to the PDB file.
    
    Returns:
        Structure: Biopython Structure object.
    """
    parser = PDBParser(QUIET=True)

    # Structure ID is the filename without its extension. splitext drops only
    # the final suffix, so a name like '1abc.v2.pdb' keeps its '1abc.v2' stem
    # instead of being cut at the first dot.
    filename = os.path.basename(file_path.replace("\\", os.sep))
    structure_id = os.path.splitext(filename)[0] or filename

    structure = parser.get_structure(structure_id, file_path)

    # PDBParser accepts any text file and returns an empty Structure when it
    # finds no coordinate records, so a mis-named or non-PDB upload would sail
    # through here and only surface much later as a misleading "cleaning
    # removed every atom" error. Reject it where the cause is still known.
    if not any(True for _ in structure.get_atoms()):
        raise ValueError(
            f"No atoms could be read from '{filename}'. The file does not "
            "contain PDB ATOM or HETATM records - it may be empty, truncated, "
            "or in another format such as mmCIF."
        )

    return structure

def save_pdb(structure, output_path, select=None):
    """
    Saves a Biopython Structure object to a PDB file.
    
    Args:
        structure (Structure): Biopython Structure object.
        output_path (str): Path to write the PDB file.
        select (Select, optional): A Biopython Select object used for filtering out atoms/residues.
    """
    io = PDBIO()
    io.set_structure(structure)
    if select:
        io.save(output_path, select)
    else:
        io.save(output_path)
