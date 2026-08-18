# BioPrep

BioPrep is a Python command-line tool for preparing protein structures for molecular docking. It performs essential preprocessing steps on PDB files to ensure they are clean, properly protonated, and ready for computational analysis.

## Features
- Removes water molecules (HOH)
- Removes ligands and heteroatoms, keeping only standard protein chains
- Allows selection of a specific protein chain
- Adds missing hydrogens at a requested physiological pH (default: 7.4) using `pdbfixer`

## Installation
You can install this via pip:
```bash
pip install -e .
```
This requires `conda` to install `pdbfixer` and `openmm` successfully in most environments.

## Usage
Basic usage to clean a structure and protonate at pH 7.4:
```bash
bioprep --input raw_structure.pdb --output cleaned_structure.pdb
```

Select a specific chain (e.g., Chain A) and adjust the pH to 7.0:
```bash
bioprep --input raw_structure.pdb --output cleaned_structure.pdb --chain A --ph 7.0
```
