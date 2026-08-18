"""
Regression tests for the bugs catalogued in BACKEND_AUDIT.md.

Each test names the audit ID it guards so a future change that reintroduces
the bug fails with an obvious label.
"""

import json
import os
import pathlib
import sys
import tempfile
import unittest

import conftest  # noqa: F401  - sets up sys.path

from bioprep.core import io, cleaner, protonator, minimizer, reporter
from bioprep.core.pipeline import PipelineSettings

HERE = os.path.dirname(os.path.abspath(__file__))
CRN = os.path.join(HERE, '..', 'bioprep', '1crn.pdb')


def _write(path, lines):
    with open(path, 'w') as fh:
        fh.writelines(lines)
        fh.write('END\n')
    return path


def _protein_lines():
    with open(CRN) as fh:
        return [l for l in fh if l.startswith('ATOM')]


def _first_ca_xyz(lines):
    ca = next(l for l in lines if l[12:16].strip() == 'CA')
    return float(ca[30:38]), float(ca[38:46]), float(ca[46:54])


class TempDirTest(unittest.TestCase):
    def setUp(self):
        if not os.path.exists(CRN):
            self.skipTest('1crn.pdb not available')
        self._dir = tempfile.TemporaryDirectory(prefix='bioprep_test_')
        self.tmp = self._dir.name

    def tearDown(self):
        self._dir.cleanup()

    def path(self, name):
        return os.path.join(self.tmp, name)


class TestCleaner(TempDirTest):

    def test_b5_protected_ligand_survives_remove_all(self):
        """B5: an explicitly protected ligand outranks 'ALL'."""
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        lines += [
            'HETATM 9001  C1  BTN A 900    %8.3f%8.3f%8.3f  1.00  0.00           C  \n' % (x + 5, y, z),
            'HETATM 9002  S   SO4 A 901    %8.3f%8.3f%8.3f  1.00  0.00           S  \n' % (x + 9, y, z),
        ]
        source = _write(self.path('in.pdb'), lines)

        structure = io.load_pdb(source)
        select = cleaner.clean_structure(
            structure, remove_heteroatoms=['ALL'], protect_ligands=['BTN'])
        io.save_pdb(structure, self.path('out.pdb'), select=select)

        text = pathlib.Path(self.path('out.pdb')).read_text()
        self.assertIn('BTN', text, 'protected ligand was deleted by ALL')
        self.assertNotIn('SO4', text, 'unprotected heteroatom survived ALL')

    def test_b6_structural_water_is_chain_qualified(self):
        """B6: a same-numbered water in another chain must not be kept."""
        protein = _protein_lines()[:40]
        x, y, z = _first_ca_xyz(protein)
        lines = list(protein) + [l[:21] + 'B' + l[22:] for l in protein] + [
            'HETATM 9500  O   HOH A 700    %8.3f%8.3f%8.3f  1.00  0.00           O  \n' % (x + 2.8, y, z),
            'HETATM 9501  O   HOH B 700    %8.3f%8.3f%8.3f  1.00  0.00           O  \n' % (x + 500, y + 500, z + 500),
        ]
        source = _write(self.path('in.pdb'), lines)

        structure = io.load_pdb(source)
        select = cleaner.clean_structure(
            structure, remove_water=True, keep_structural_waters=True)
        io.save_pdb(structure, self.path('out.pdb'), select=select)

        kept = [l for l in pathlib.Path(self.path('out.pdb')).read_text().splitlines()
                if l.startswith('HETATM') and 'HOH' in l]
        chains = sorted({l[21] for l in kept})
        self.assertEqual(chains, ['A'],
                         'a distant water in another chain was preserved')


class TestProtonator(TempDirTest):

    def test_b7_a9_ligand_held_out_and_status_reported(self):
        """A9/A10/B7: ligands bypass PDBFixer and the result is reported."""
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        lines.append(
            'HETATM 9001  C1  BTN A 900    %8.3f%8.3f%8.3f  1.00  0.00           C  \n' % (x + 5, y, z))
        source = _write(self.path('in.pdb'), lines)

        result = protonator.add_hydrogens(source, self.path('out.pdb'), ph=7.4)

        self.assertTrue(result['hydrogens_added'])
        self.assertIn('BTN', result['ligands_preserved'])
        self.assertIsInstance(result['warnings'], list)

        out = [l for l in pathlib.Path(self.path('out.pdb')).read_text().splitlines() if 'BTN' in l]
        self.assertEqual(len(out), 1, 'ligand lost during protonation')
        self.assertEqual(out[0][30:54], lines[-1][30:54],
                         'ligand coordinates were modified by PDBFixer')

    def test_ligand_only_input_raises_clearly(self):
        source = _write(self.path('lig.pdb'), [
            'HETATM 9001  C1  BTN A 900       1.000   1.000   1.000  1.00  0.00           C  \n'])
        with self.assertRaises(ValueError):
            protonator.add_hydrogens(source, self.path('out.pdb'))


class TestMinimizer(TempDirTest):

    def _prepared(self, extra_lines=()):
        lines = _protein_lines() + list(extra_lines)
        source = _write(self.path('raw.pdb'), lines)
        structure = io.load_pdb(source)
        io.save_pdb(structure, self.path('clean.pdb'),
                    select=cleaner.clean_structure(structure, remove_water=False))
        protonator.add_hydrogens(self.path('clean.pdb'), self.path('prot.pdb'))
        return self.path('prot.pdb')

    def test_b1_b3_ligand_falls_back_to_partial_tier(self):
        """B1/B3: a ligand must trigger Tier 2, not a silent pass-through."""
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        source = self._prepared([
            'HETATM 9001  C1  BTN A 900    %8.3f%8.3f%8.3f  1.00  0.00           C  \n' % (x + 5, y, z)])

        stats = minimizer.minimize_structure(source, self.path('min.pdb'))

        self.assertEqual(stats['status'], 'partial')
        self.assertIn('BTN', stats['excluded_residues'])
        self.assertIsNone(stats['error'])
        self.assertTrue(stats['energy_decreased'])

        before = [l for l in pathlib.Path(source).read_text().splitlines() if 'BTN' in l][0][30:54]
        after = [l for l in pathlib.Path(self.path('min.pdb')).read_text().splitlines() if 'BTN' in l][0][30:54]
        self.assertEqual(before, after, 'excluded ligand was moved')

        original = pathlib.Path(source).read_bytes()
        minimised = pathlib.Path(self.path('min.pdb')).read_bytes()
        self.assertNotEqual(original, minimised,
                            'output identical to input - nothing was minimised')

    def test_b2_retained_water_minimises_with_gbsa(self):
        """B2: the GBSA force field must include water parameters."""
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        source = self._prepared([
            'HETATM 9100  O   HOH A 800    %8.3f%8.3f%8.3f  1.00  0.00           O  \n' % (x + 2.8, y, z)])

        stats = minimizer.minimize_structure(source, self.path('min.pdb'),
                                             use_gbsa=True)
        self.assertEqual(stats['status'], 'full')
        self.assertIsNone(stats['error'])

    def test_b16_result_is_strict_json_serialisable(self):
        """B16: energies are numbers or None, never NaN or the string 'N/A'."""
        stats = minimizer.minimize_structure(CRN, self.path('min.pdb'))
        json.dumps(stats, allow_nan=False)
        for key in ('energy_before_kJ_mol', 'energy_after_kJ_mol'):
            self.assertTrue(stats[key] is None or isinstance(stats[key], float))

    def test_failure_reports_status_failed(self):
        """B3: an unusable input must not be reported as minimised."""
        source = _write(self.path('junk.pdb'), ['REMARK nothing here\n'])
        stats = minimizer.minimize_structure(source, self.path('min.pdb'))
        self.assertEqual(stats['status'], 'failed')
        self.assertIsNotNone(stats['error'])
        self.assertTrue(os.path.exists(self.path('min.pdb')))


class TestReporter(unittest.TestCase):

    def test_b24_b25_renders_partial_report(self):
        """B24/B25: a partial report must render, not raise KeyError."""
        text = reporter.report_to_text({'input_file': 'x.pdb'})
        self.assertIn('x.pdb', text)

    def test_b24_renders_empty_report(self):
        self.assertIn('BIOPREP', reporter.report_to_text({}))

    def test_b7_protonation_status_comes_from_protonator(self):
        report = reporter.build_report(
            filename='x.pdb', chains_detected=['A'], chains_retained=['A'],
            waters_removed=0, heteroatoms_removed=[], heteroatoms_retained=[],
            ph_used=7.4, missing_residues=[], atoms_before=100, atoms_after=110,
            protonation={'hydrogens_added': False, 'warnings': ['boom']},
        )
        self.assertFalse(report['hydrogens_added'])
        self.assertFalse(report['protonation']['hydrogens_added'])
        self.assertEqual(report['atom_counts']['delta'], 10)


class TestPipelineSettings(unittest.TestCase):

    def test_defaults(self):
        s = PipelineSettings.from_mapping({})
        self.assertEqual(s.ph, 7.4)
        self.assertTrue(s.remove_water)
        self.assertIsNone(s.chains)

    def test_malformed_values_fall_back(self):
        s = PipelineSettings.from_mapping(
            {'ph': 'not-a-number', 'force_field': 'bogus', 'docking_target': 'evil'})
        self.assertEqual(s.ph, 7.4)
        self.assertEqual(s.force_field, 'amber14')
        self.assertIsNone(s.docking_target)

    def test_ph_is_clamped(self):
        self.assertEqual(PipelineSettings.from_mapping({'ph': 99}).ph, 14.0)
        self.assertEqual(PipelineSettings.from_mapping({'ph': -5}).ph, 0.0)

    def test_string_booleans_from_form(self):
        s = PipelineSettings.from_mapping({'run_minimization': 'true',
                                           'remove_water': 'false'})
        self.assertTrue(s.run_minimization)
        self.assertFalse(s.remove_water)

    def test_residue_names_normalised(self):
        s = PipelineSettings.from_mapping({'protect_ligands': [' btn ', '']})
        self.assertEqual(s.protect_ligands, ['BTN'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
