from Bio.PDB import Select, NeighborSearch

from .residues import is_water


def _is_water_residue(residue):
    """
    True for any water, however it is named.

    Biopython only tags HOH and WAT with hetfield 'W'. SOL, TIP3, TIP4, DOD and
    friends arrive as ordinary heteroatoms, so a bare ``hetfield == 'W'`` test
    lets them slip past water removal and into the ligand branch.
    """
    return residue.id[0] == 'W' or is_water(residue.resname)


class BioPrepSelect(Select):
    """
    Biopython Select class to optionally filter out water, specific heteroatoms,
    and restrict to specific chains without altering the base coordinate geometry.
    """
    def __init__(self, target_chains=None, remove_water=True, remove_heteroatoms=None,
                 keep_structural_waters=False, structural_waters=None, protect_ligands=None,
                 first_model_id=0):
        self._first_model_id = first_model_id
        self.target_chains = target_chains if target_chains else []
        self.remove_water = remove_water
        self.remove_heteroatoms = remove_heteroatoms if remove_heteroatoms else []
        self.keep_structural_waters = keep_structural_waters
        # Keys are (chain_id, residue.id) — residue ids repeat across chains.
        self.structural_waters = structural_waters if structural_waters else set()
        self.protect_ligands = protect_ligands if protect_ligands else []

    def accept_model(self, model):
        # Write only the first model. NMR ensembles otherwise emit every model,
        # and the downstream tools silently use the first one anyway — so the
        # extra models only inflate the file and the atom counts.
        return 1 if model.id == self._first_model_id else 0

    def accept_chain(self, chain):
        if self.target_chains and chain.id not in self.target_chains:
            return 0
        return 1

    def accept_residue(self, residue):
        hetfield = residue.id[0]
        chain_id = residue.get_parent().id

        # 1. Handle Water
        if _is_water_residue(residue):
            if self.keep_structural_waters and (chain_id, residue.id) in self.structural_waters:
                return 1
            return 0 if self.remove_water else 1

        # 2. Handle Heteroatoms (Ligands, ions, lipids, etc.)
        # Biopython heteroatoms have id[0] NOT starting with a space ' '
        # (Standard residues are ' ', waters are 'W', heteroatoms are 'H_xxx')
        if hetfield != ' ':
            res_name = residue.resname.strip()

            # An explicitly protected ligand outranks every removal rule,
            # including 'ALL'. Without this, "remove all heteroatoms" silently
            # deleted ligands the user had asked to keep.
            if res_name in self.protect_ligands:
                return 1

            # If user selected "Remove ALL Heteros"
            if 'ALL' in self.remove_heteroatoms:
                return 0

            # If THIS specific residue name is in the removal list, remove it
            if res_name in self.remove_heteroatoms:
                return 0

            # OTHERWISE: Keep it (Surgical preservation)
            return 1

        # 3. Handle Standard Protein Residues
        return 1


def clean_structure(structure, target_chains=None, remove_water=True, remove_heteroatoms=None,
                    keep_structural_waters=False, protect_ligands=None):
    """
    Returns a configured Biopython Select object that filters unwanted
    residues, water molecules, and restricts to targeted chains.

    Structural water detection: Build a NeighborSearch from PROTEIN atoms,
    then query each water atom to find waters within 4.0Å of any protein atom.
    This is the correct direction — protein is the reference set.
    """
    structural_waters = set()
    first_model_id = next((model.id for model in structure), 0)

    if keep_structural_waters and remove_water:
        # Collect protein (standard residue) atoms as the reference set
        protein_atoms = []
        water_residues = []

        for model in structure:
            if model.id != first_model_id:
                break  # only the first model is written out
            for chain in model:
                if target_chains and chain.id not in target_chains:
                    continue
                for residue in chain:
                    if residue.id[0] == ' ':
                        # Standard amino acid — add all atoms to reference
                        protein_atoms.extend(residue.get_atoms())
                    elif _is_water_residue(residue):
                        water_residues.append((chain.id, residue))

        # Build NeighborSearch from protein atoms (the reference), then query
        # each water's oxygen to see if it's close to ANY protein atom.
        if protein_atoms and water_residues:
            ns = NeighborSearch(protein_atoms)
            for chain_id, water_res in water_residues:
                for w_atom in water_res.get_atoms():
                    # Standard water oxygen is usually 'O' or 'OW'
                    if w_atom.get_name() in ('O', 'OW', 'O1'):
                        nearby_protein_atoms = ns.search(w_atom.get_coord(), 4.0)
                        if nearby_protein_atoms:
                            # Qualify by chain: residue ids repeat across chains,
                            # so a bare id would keep unrelated waters elsewhere.
                            structural_waters.add((chain_id, water_res.id))
                            break  # Only need one hit per water molecule

    return BioPrepSelect(
        target_chains=target_chains,
        remove_water=remove_water,
        remove_heteroatoms=remove_heteroatoms,
        keep_structural_waters=keep_structural_waters,
        structural_waters=structural_waters,
        protect_ligands=protect_ligands,
        first_model_id=first_model_id,
    )
