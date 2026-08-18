import os
from Bio.PDB import PDBList

pdbl = PDBList()
# 1CRN is a small, well-known protein structure (Crambin)
pdbl.retrieve_pdb_file('1CRN', pdir='.', file_format='pdb')
os.rename('pdb1crn.ent', '1crn.pdb')
print("Downloaded 1crn.pdb")
