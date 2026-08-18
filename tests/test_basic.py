import unittest
import os
from bioprep.core.analyzer import analyze_structure, detect_missing_residues
from bioprep.core.reporter import build_report, count_atoms_in_pdb
from Bio.PDB import PDBParser

class TestBioPrepCore(unittest.TestCase):
    def setUp(self):
        # We need a small PDB for testing. 1CRN is a standard small one.
        self.test_pdb = os.path.join(os.path.dirname(__file__), '..', 'bioprep', '1crn.pdb')
        self.parser = PDBParser(QUIET=True)

    def test_analyzer(self):
        if not os.path.exists(self.test_pdb):
            self.skipTest("1crn.pdb not found for testing")
        
        structure = self.parser.get_structure("test", self.test_pdb)
        meta = analyze_structure(structure)
        
        self.assertIn('chains', meta)
        self.assertIn('atoms_total', meta)
        self.assertGreater(meta['atoms_total'], 0)
        self.assertIsInstance(meta['chains'], list)

    def test_atom_count(self):
        if not os.path.exists(self.test_pdb):
            self.skipTest("1crn.pdb not found for testing")
            
        count = count_atoms_in_pdb(self.test_pdb)
        self.assertGreater(count, 0)

    def test_report_builder(self):
        report = build_report(
            filename="test.pdb",
            chains_detected=['A'],
            chains_retained=['A'],
            waters_removed=0,
            heteroatoms_removed=[],
            heteroatoms_retained=[],
            hydrogens_added=True,
            ph_used=7.4,
            missing_residues=[],
            atoms_before=100,
            atoms_after=110
        )
        self.assertEqual(report['input_file'], "test.pdb")
        self.assertEqual(report['atom_counts']['delta'], 10)

if __name__ == '__main__':
    unittest.main()
