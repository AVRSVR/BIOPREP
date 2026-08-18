"""
Shared residue classification tables.

Both the protonator and the minimizer need to answer the same question —
"is this residue something a protein force field knows about?" — so the sets
live here rather than being duplicated (and drifting) in each module.
"""

WATER = {
    'HOH', 'WAT', 'H2O', 'TIP', 'TIP3', 'TIP4', 'SOL', 'DOD', 'D2O',
}

# 20 standard amino acids plus the protonation/termini variants that AMBER and
# CHARMM emit. PDBFixer and OpenMM both understand these.
AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    # protonation states / disulfide / CHARMM naming
    'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP',
    'CYX', 'CYM', 'ASH', 'GLH', 'LYN', 'TYM', 'ARN',
    # terminal caps
    'ACE', 'NME', 'NMA', 'NH2', 'FOR',
    # selenomethionine — common in crystal structures
    'MSE',
}

NUCLEIC_ACIDS = {
    'A', 'C', 'G', 'T', 'U', 'I',
    'DA', 'DC', 'DG', 'DT', 'DU', 'DI',
    'RA', 'RC', 'RG', 'RU',
    'A3', 'A5', 'C3', 'C5', 'G3', 'G5', 'T3', 'T5', 'U3', 'U5',
    'DA3', 'DA5', 'DC3', 'DC5', 'DG3', 'DG5', 'DT3', 'DT5',
}

# Monatomic ions that standard force fields do parameterise.
IONS = {
    'NA', 'K', 'LI', 'RB', 'CS', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'CO',
    'NI', 'CD', 'CL', 'BR', 'I', 'F', 'SO4', 'PO4',
}

#: Residues a protein force field can build a template for without extra work.
FORCE_FIELD_SAFE = AMINO_ACIDS | NUCLEIC_ACIDS | WATER


def is_water(resname):
    return resname.strip().upper() in WATER


def is_standard(resname):
    """True if the residue is a normal biopolymer unit (not a ligand)."""
    return resname.strip().upper() in (AMINO_ACIDS | NUCLEIC_ACIDS)


def is_force_field_safe(resname):
    """True if a protein force field is expected to have a template for it."""
    return resname.strip().upper() in FORCE_FIELD_SAFE
