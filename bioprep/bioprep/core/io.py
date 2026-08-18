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
    
    # Generate a dummy structure ID based on the filename
    filename = file_path.replace("\\", "/").split("/")[-1]
    structure_id = filename.split(".")[0]
    
    structure = parser.get_structure(structure_id, file_path)
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
