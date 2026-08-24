from Bio.PDB import Select, NeighborSearch

from .residues import is_water


#: A water is kept as "structural" when its oxygen sits this close to a protein
#: atom that can actually hydrogen bond to it. The earlier rule was 4.0 A to
#: ANY protein atom, carbon included, which keeps water packed against a
#: hydrophobic surface with nothing holding it there. On streptavidin that was
#: 8 of 78 retained waters. 3.5 A is the usual upper bound for an O/N-H...O
#: hydrogen bond between heavy atoms.
STRUCTURAL_WATER_CUTOFF = 3.5
HBOND_CAPABLE_ELEMENTS = {'N', 'O', 'S'}


def _can_hydrogen_bond(atom):
    """True for protein atoms able to donate or accept a hydrogen bond."""
    element = (atom.element or '').strip().upper()
    if element:
        return element in HBOND_CAPABLE_ELEMENTS
    return atom.get_name().strip().upper()[:1] in HBOND_CAPABLE_ELEMENTS


def _normalise_names(names):
    """Upper-case and strip a list of residue names for comparison."""
    return [str(name).strip().upper() for name in (names or []) if str(name).strip()]


def _is_water_oxygen(atom):
    """
    True for the oxygen of a water molecule, whatever it is called.

    Crystallographic waters use O; GROMACS uses OW; CHARMM and NAMD use OH2.
    Matching a fixed list of names meant a CHARMM-derived structure had no
    detectable water oxygen at all, so structural-water preservation silently
    kept nothing and every water was deleted.
    """
    if (atom.element or '').strip().upper() == 'O':
        return True
    return atom.get_name().strip().upper().startswith('O')


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
        # Residue names are compared case-insensitively. PDB files use upper
        # case, but a name typed by hand on the CLI may not be, and a lower
        # case entry in protect_ligands used to silently fail to match - so a
        # ligand the user had asked to keep was deleted without a word.
        self.remove_heteroatoms = _normalise_names(remove_heteroatoms)
        self.keep_structural_waters = keep_structural_waters
        # Keys are (chain_id, residue.id) — residue ids repeat across chains.
        self.structural_waters = structural_waters if structural_waters else set()
        self.protect_ligands = _normalise_names(protect_ligands)

    def accept_model(self, model):
        # Write only the first model. NMR ensembles otherwise emit every model,
        # and the downstream tools silently use the first one anyway — so the
        # extra models only inflate the file and the atom counts.
        return 1 if model.id == self._first_model_id else 0

    def accept_chain(self, chain):
        # Always accept at this level. Select's accept_chain is an all-or-
        # nothing gate - returning 0 here means accept_residue is never even
        # called for anything in the chain, which made target_chains outrank
        # protect_ligands rather than the other way round. A co-crystallised
        # ligand is routinely recorded under a different chain letter than
        # the receptor: HIV protease's inhibitor sits on chain B of an A/B
        # dimer. Selecting chains=['A'] silently dropped it even with the
        # ligand explicitly protected, because the chain gate ran first and
        # accept_residue's protection check never got the chance to fire.
        # The actual chain decision now happens per residue, below.
        return 1

    def accept_residue(self, residue):
        hetfield = residue.id[0]
        chain_id = residue.get_parent().id
        is_water = _is_water_residue(residue)

        # An explicitly protected ligand survives regardless of which chain
        # it is recorded under. Chain selection targets the receptor; a bound
        # ligand naming itself distinct from the protein is exactly the case
        # protect_ligands exists for, so requiring target_chains membership
        # as well would make "protect this ligand" an unreliable promise.
        if hetfield != ' ' and not is_water:
            if residue.resname.strip().upper() in self.protect_ligands:
                return 1

        if self.target_chains and chain_id not in self.target_chains:
            return 0

        # 1. Handle Water
        if is_water:
            if self.keep_structural_waters and (chain_id, residue.id) in self.structural_waters:
                return 1
            return 0 if self.remove_water else 1

        # 2. Handle Heteroatoms (Ligands, ions, lipids, etc.)
        # Biopython heteroatoms have id[0] NOT starting with a space ' '
        # (Standard residues are ' ', waters are 'W', heteroatoms are 'H_xxx')
        if hetfield != ' ':
            res_name = residue.resname.strip().upper()

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
                        # Only atoms that can hydrogen bond count as anchors;
                        # proximity to a carbon does not hold a water in place.
                        protein_atoms.extend(
                            a for a in residue.get_atoms() if _can_hydrogen_bond(a))
                    elif _is_water_residue(residue):
                        water_residues.append((chain.id, residue))

        # Build NeighborSearch from protein atoms (the reference), then query
        # each water's oxygen to see if it's close to ANY protein atom.
        if protein_atoms and water_residues:
            ns = NeighborSearch(protein_atoms)
            for chain_id, water_res in water_residues:
                for w_atom in water_res.get_atoms():
                    if _is_water_oxygen(w_atom):
                        nearby_protein_atoms = ns.search(w_atom.get_coord(),
                                                    STRUCTURAL_WATER_CUTOFF)
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
