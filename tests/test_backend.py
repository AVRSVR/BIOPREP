"""
Regression tests for the bugs catalogued in BACKEND_AUDIT.md.

Each test names the audit ID it guards so a future change that reintroduces
the bug fails with an obvious label.
"""

import json
import os
import pathlib
import shutil
import tempfile
import unittest

import conftest  # noqa: F401  - sets up sys.path

from bioprep.core import (analyzer, cleaner, exporter, io, minimizer,
                          protonator, reporter, site_analyzer)
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


class TestIO(TempDirTest):
    """Feature 1: PDB loading via PDBParser(QUIET=True)."""

    ATOM_LINE = ('ATOM      1  N   ALA A   1      11.104   6.134  -6.504'
                 '  1.00  0.00           N\n')

    def test_feature1_loads_pdb_and_returns_structure(self):
        structure = io.load_pdb(CRN)
        self.assertEqual(sum(1 for _ in structure.get_atoms()), 327)
        self.assertEqual(structure.id, '1crn')

    def test_feature1_quiet_suppresses_construction_warnings(self):
        """QUIET=True is the feature; a messy file must not emit warnings."""
        import warnings as _warnings

        messy = _write(self.path('messy.pdb'), [
            self.ATOM_LINE,
            self.ATOM_LINE.replace('  N   ALA', '  CA  ALA'),
            self.ATOM_LINE.replace('  N   ALA', '  CA  ALA'),   # duplicate atom
            self.ATOM_LINE.replace('ALA A   1', 'GLY B   1'),   # chain break
            self.ATOM_LINE.replace('ALA A   1', 'ALA A   2'),
        ])

        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter('always')
            io.load_pdb(messy)

        self.assertEqual(len(caught), 0,
                         'QUIET=True did not suppress parser warnings')

    def test_feature2_structure_id_keeps_the_whole_stem(self):
        for name, expected in (('1abc.pdb', '1abc'), ('1abc.v2.pdb', '1abc.v2')):
            path = _write(self.path(name), [self.ATOM_LINE])
            self.assertEqual(io.load_pdb(path).id, expected)

    def test_feature1_non_pdb_file_is_rejected_at_load(self):
        """
        PDBParser returns an empty Structure for any text file. Without a guard
        the failure surfaced later as a misleading 'cleaning removed every atom'.
        """
        junk = _write(self.path('junk.pdb'), ['this is not a pdb file\n'])
        with self.assertRaises(ValueError) as caught:
            io.load_pdb(junk)
        self.assertIn('No atoms', str(caught.exception))

    def test_feature1_unparseable_mmcif_is_rejected(self):
        cif = _write(self.path('s.cif'), ['data_1CRN\n', '_entry.id 1CRN\n'])
        with self.assertRaises(ValueError):
            io.load_pdb(cif)


class TestSavePdb(TempDirTest):
    """Features 3 and 4: PDBIO saving, and saving through a Select filter."""

    METALS = [
        'HETATM  900 ZN    ZN A 300      12.000  12.000  12.000  1.00 10.00          ZN  \n',
        'HETATM  901 MG    MG A 301      16.000  12.000  12.000  1.00 10.00          MG  \n',
        'HETATM  902 FE    FE A 302      20.000  12.000  12.000  1.00 10.00          FE  \n',
        'HETATM  903 CL    CL A 303      24.000  12.000  12.000  1.00 10.00          CL  \n',
    ]

    def test_feature3_element_column_matches_the_spec(self):
        """
        Biopython writes column 67 onward one place early, which truncates
        two-character element symbols. Compare against a real RCSB line.
        """
        structure = io.load_pdb(CRN)
        io.save_pdb(structure, self.path('out.pdb'))

        rcsb = next(l.rstrip('\n') for l in pathlib.Path(CRN).read_text().splitlines()
                    if l.startswith('ATOM'))
        ours = next(l for l in pathlib.Path(self.path('out.pdb')).read_text().splitlines()
                    if l.startswith('ATOM'))

        self.assertEqual(len(ours), 80, 'line is not 80 columns')
        self.assertEqual(ours[76:78], rcsb[76:78],
                         'element column does not match the reference file')

    def test_feature3_two_character_metals_round_trip(self):
        """ZN must not come back as N, nor FE as F."""
        source = _write(self.path('metals.pdb'), _protein_lines() + self.METALS)

        original = io.load_pdb(source)
        expected = {a.get_parent().get_resname().strip(): a.element
                    for a in original.get_atoms() if a.get_parent().id[0] != ' '}
        io.save_pdb(original, self.path('out.pdb'))

        reloaded = {a.get_parent().get_resname().strip(): a.element
                    for a in io.load_pdb(self.path('out.pdb')).get_atoms()
                    if a.get_parent().id[0] != ' '}

        self.assertEqual(expected, reloaded)
        self.assertEqual(reloaded.get('ZN'), 'ZN')
        self.assertEqual(reloaded.get('FE'), 'FE')

    def test_feature3_openmm_reads_metals_from_our_output(self):
        """The downstream consumer must agree; OpenMM read ZN as N before."""
        from openmm.app import PDBFile

        source = _write(self.path('metals.pdb'), _protein_lines() + self.METALS)
        io.save_pdb(io.load_pdb(source), self.path('out.pdb'))

        symbols = {
            atom.residue.name.strip(): (atom.element.symbol if atom.element else None)
            for atom in PDBFile(self.path('out.pdb')).topology.atoms()
            if atom.residue.name.strip() in ('ZN', 'MG', 'FE', 'CL')
        }
        self.assertEqual(symbols,
                         {'ZN': 'Zn', 'MG': 'Mg', 'FE': 'Fe', 'CL': 'Cl'})

    def test_feature3_preserves_occupancy_bfactor_altloc_icode(self):
        rows = [
            'ATOM      1  N   ALA A   1       1.000   2.000   3.000  0.75 12.34           N  \n',
            'ATOM      2  CA AALA A  50      10.000  10.000  10.000  0.50 20.00           C  \n',
            'ATOM      3  CA BALA A  50      10.500  10.000  10.000  0.50 20.00           C  \n',
            'ATOM      4  N   GLY A  51A     11.000  11.000  11.000  1.00 15.00           N  \n',
        ]
        source = _write(self.path('in.pdb'), rows)
        io.save_pdb(io.load_pdb(source), self.path('out.pdb'))

        written = [l for l in pathlib.Path(self.path('out.pdb')).read_text().splitlines()
                   if l.startswith('ATOM')]
        self.assertEqual(len(written), 4)

        first = next(l for l in written if l[22:26].strip() == '1')
        self.assertEqual(first[54:60].strip(), '0.75')
        self.assertEqual(first[60:66].strip(), '12.34')

        self.assertEqual({l[16] for l in written if l[22:26].strip() == '50'},
                         {'A', 'B'})
        self.assertTrue(any(l[26] == 'A' for l in written
                            if l[22:26].strip() == '51'))

    def test_feature3_creates_missing_parent_directories(self):
        io.save_pdb(io.load_pdb(CRN), self.path('a/b/c/out.pdb'))
        self.assertTrue(os.path.isfile(self.path('a/b/c/out.pdb')))

    def test_feature4_select_filters_what_is_written(self):
        from Bio.PDB import Select

        source = _write(self.path('in.pdb'), _protein_lines() + self.METALS)
        structure = io.load_pdb(source)

        class ProteinOnly(Select):
            def accept_residue(self, residue):
                return 1 if residue.id[0] == ' ' else 0

        io.save_pdb(structure, self.path('out.pdb'), select=ProteinOnly())
        text = pathlib.Path(self.path('out.pdb')).read_text()

        for metal in ('ZN', 'MG', 'FE'):
            self.assertNotIn(metal, text)
        self.assertTrue(any(l.startswith('ATOM') for l in text.splitlines()))

    def test_feature4_select_rejecting_everything_still_writes_a_file(self):
        from Bio.PDB import Select

        class RejectAll(Select):
            def accept_residue(self, residue):
                return 0

        io.save_pdb(io.load_pdb(CRN), self.path('empty.pdb'), select=RejectAll())
        self.assertTrue(os.path.isfile(self.path('empty.pdb')))
        text = pathlib.Path(self.path('empty.pdb')).read_text()
        self.assertEqual([l for l in text.splitlines()
                          if l.startswith(('ATOM', 'HETATM'))], [])


class TestMmcifSupport(TempDirTest):
    """mmCIF input: RCSB serves it by default, so it must be accepted."""

    def _make_cif(self, name='1crn.cif'):
        from Bio.PDB import MMCIFIO
        writer = MMCIFIO()
        writer.set_structure(io.load_pdb(CRN))
        path = self.path(name)
        writer.save(path)
        return path

    def test_detect_format_by_extension(self):
        self.assertEqual(io.detect_format(CRN), 'pdb')
        self.assertEqual(io.detect_format(self._make_cif()), 'mmcif')

    def test_detect_format_prefers_content_over_extension(self):
        """A .cif downloaded from RCSB and renamed .pdb must still work."""
        cif = self._make_cif()
        renamed = self.path('sneaky.pdb')
        shutil.copy(cif, renamed)

        self.assertEqual(io.detect_format(renamed), 'mmcif')
        self.assertEqual(sum(1 for _ in io.load_structure(renamed).get_atoms()), 327)

    def test_mmcif_loads_with_same_atom_count_as_pdb(self):
        cif = io.load_structure(self._make_cif())
        pdb = io.load_structure(CRN)
        self.assertEqual(sum(1 for _ in cif.get_atoms()),
                         sum(1 for _ in pdb.get_atoms()))

    def test_ensure_pdb_passes_pdb_through_untouched(self):
        path, converted = io.ensure_pdb(CRN, self.tmp)
        self.assertEqual(path, CRN)
        self.assertIsNone(converted)

    def test_ensure_pdb_converts_mmcif(self):
        path, converted = io.ensure_pdb(self._make_cif(), self.tmp)
        self.assertEqual(converted, 'mmcif')
        self.assertTrue(path.endswith('.pdb'))
        self.assertEqual(sum(1 for _ in io.load_pdb(path).get_atoms()), 327)

    def test_structures_too_large_for_pdb_are_refused(self):
        """
        PDB has five columns of atom serial and one of chain id. Writing a
        structure that exceeds either produces a silently corrupt file, so
        conversion must refuse instead.
        """
        class Chain:
            def __init__(self, identifier):
                self.id = identifier

        class Model:
            def __init__(self, chains):
                self._chains = chains

            def __iter__(self):
                return iter(self._chains)

        class Structure:
            def __init__(self, chains):
                self._models = [Model(chains)]

            def __iter__(self):
                return iter(self._models)

            def get_atoms(self):
                return iter([])

        problems = io.check_pdb_representable(Structure([Chain('AAA'), Chain('B')]))
        self.assertTrue(problems)
        self.assertIn('AAA', problems[0])

        self.assertEqual(io.check_pdb_representable(io.load_pdb(CRN)), [])


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

    def test_water_removal_covers_non_hoh_names(self):
        """Biopython tags only HOH/WAT as 'W'; SOL and TIP3 arrive as heteroatoms."""
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        lines += [
            'HETATM 9100  O   HOH A 700    %8.3f%8.3f%8.3f  1.00  0.00           O  \n' % (x + 40, y, z),
            'HETATM 9101  O   SOL A 701    %8.3f%8.3f%8.3f  1.00  0.00           O  \n' % (x + 44, y, z),
            'HETATM 9102  O   TIP3A 702    %8.3f%8.3f%8.3f  1.00  0.00           O  \n' % (x + 48, y, z),
            'HETATM 9103  C1  BTN A 900    %8.3f%8.3f%8.3f  1.00  0.00           C  \n' % (x + 5, y, z),
        ]
        source = _write(self.path('in.pdb'), lines)

        structure = io.load_pdb(source)
        io.save_pdb(structure, self.path('out.pdb'),
                    select=cleaner.clean_structure(structure, remove_water=True))

        text = pathlib.Path(self.path('out.pdb')).read_text()
        for name in ('HOH', 'SOL', 'TIP'):
            self.assertNotIn(name, text, '%s water survived remove_water' % name)
        self.assertIn('BTN', text, 'ligand was removed along with the waters')


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

    def test_failed_parse_reports_the_real_error_and_leaks_nothing(self):
        """PDBFixer leaves its handle open on a parse failure.

        On Windows the open handle made the temp-file cleanup raise WinError 32
        from inside the finally block, replacing the real diagnostic with a
        message naming an internal temp path.
        """
        import glob as _glob

        source = _write(self.path('bad.pdb'), [
            'ATOM      1  N   ALA A   1      NOTNUMERIC\n'])

        pattern = os.path.join(tempfile.gettempdir(), '*_prot*.pdb')
        before = set(_glob.glob(pattern))

        with self.assertRaises(Exception) as caught:
            protonator.add_hydrogens(source, self.path('out.pdb'))

        message = str(caught.exception)
        self.assertNotIn('WinError', message, 'cleanup masked the real error')
        self.assertNotIn(tempfile.gettempdir(), message,
                         'an internal temp path leaked into the error')
        self.assertEqual(set(_glob.glob(pattern)), before,
                         'temp files leaked on the failure path')

    def test_multi_model_input_is_not_concatenated(self):
        """MODEL/ENDMDL are dropped on write, so later models must not be read."""
        lines = _protein_lines()
        single = _write(self.path('single.pdb'), lines)

        multi = self.path('multi.pdb')
        with open(multi, 'w') as fh:
            for index in (1, 2, 3):
                fh.write('MODEL     %4d\n' % index)
                fh.writelines(lines)
                fh.write('ENDMDL\n')
            fh.write('END\n')

        def atom_count(path):
            return sum(1 for l in pathlib.Path(path).read_text().splitlines()
                       if l.startswith(('ATOM', 'HETATM')))

        protonator.add_hydrogens(single, self.path('out_single.pdb'))
        protonator.add_hydrogens(multi, self.path('out_multi.pdb'))

        self.assertEqual(atom_count(self.path('out_multi.pdb')),
                         atom_count(self.path('out_single.pdb')),
                         'models were welded into one chimeric protein')

    def test_conect_records_follow_the_last_atom(self):
        """The PDB spec puts CONECT after every coordinate record."""
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        lines += [
            'HETATM 9001  C1  BTN A 900    %8.3f%8.3f%8.3f  1.00  0.00           C  \n' % (x + 5, y, z),
            'HETATM 9002  C2  BTN A 900    %8.3f%8.3f%8.3f  1.00  0.00           C  \n' % (x + 6.5, y, z),
            'CONECT 9001 9002\n',
            'CONECT 9002 9001\n',
        ]
        source = _write(self.path('in.pdb'), lines)
        protonator.add_hydrogens(source, self.path('out.pdb'))

        written = pathlib.Path(self.path('out.pdb')).read_text().splitlines()
        last_atom = max(i for i, l in enumerate(written)
                        if l.startswith(('ATOM', 'HETATM')))
        first_conect = min(i for i, l in enumerate(written) if l.startswith('CONECT'))
        self.assertGreater(first_conect, last_atom,
                           'a CONECT record precedes an atom record')

        # The ligand's own bonds must point at its renumbered serials.
        btn = [l for l in written if l[17:20].strip() == 'BTN']
        serials = {int(l[6:11]) for l in btn}
        bonded = set()
        for line in (l for l in written if l.startswith('CONECT')):
            body = line[6:].rstrip()
            refs = {int(body[i:i + 5]) for i in range(0, len(body), 5)
                    if body[i:i + 5].strip()}
            if refs & serials:
                bonded |= refs
        self.assertEqual(bonded, serials, 'ligand CONECT serials were not remapped')


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


class TestAnalyzer(TempDirTest):

    def _multi_model(self, models=3):
        lines = _protein_lines()
        path = self.path('multi.pdb')
        with open(path, 'w') as fh:
            for index in range(1, models + 1):
                fh.write('MODEL     %4d\n' % index)
                fh.writelines(lines)
                fh.write('ENDMDL\n')
            fh.write('END\n')
        return path, len(lines)

    def test_b10_counts_do_not_scale_with_model_count(self):
        """B10: an NMR ensemble must not multiply atom and water counts."""
        multi, per_model = self._multi_model(3)
        single = _write(self.path('single.pdb'), _protein_lines())

        one = analyzer.analyze_structure(io.load_pdb(single))
        many = analyzer.analyze_structure(io.load_pdb(multi))

        self.assertEqual(one['atoms_total'], many['atoms_total'])
        self.assertEqual(many['atoms_total'], per_model)
        self.assertEqual(many['model_count'], 3)

    def test_b42_only_first_model_is_written(self):
        """B42: saving must not emit every model of an ensemble."""
        multi, per_model = self._multi_model(3)
        structure = io.load_pdb(multi)
        io.save_pdb(structure, self.path('out.pdb'),
                    select=cleaner.clean_structure(structure))

        written = pathlib.Path(self.path('out.pdb')).read_text().splitlines()
        atoms = [l for l in written if l.startswith('ATOM')]
        self.assertEqual(len(atoms), per_model)

    def test_b11_missing_residues_report_real_chain_ids(self):
        """B11: the PDBFixer key is (chain index, position), not (model, chain)."""
        with open(CRN) as fh:
            seqres = [l for l in fh if l.startswith('SEQRES')]
        kept = [l for l in _protein_lines() if not (20 <= int(l[22:26]) <= 24)]
        source = _write(self.path('gap.pdb'), seqres + kept)

        missing = analyzer.detect_missing_residues(source)
        self.assertTrue(missing, 'no missing residues detected')
        self.assertEqual({m['chain'] for m in missing}, {'A'})
        for entry in missing:
            self.assertIsInstance(entry['residue'], str)

    def test_sequence_gap_detected_once(self):
        with open(CRN) as fh:
            seqres = [l for l in fh if l.startswith('SEQRES')]
        kept = [l for l in _protein_lines() if not (20 <= int(l[22:26]) <= 24)]
        source = _write(self.path('gap.pdb'), seqres + kept)

        gaps = analyzer.analyze_structure(io.load_pdb(source))['sequence_gaps']
        self.assertEqual(len(gaps), 1)
        self.assertEqual(gaps[0]['missing_count'], 5)
        self.assertEqual(gaps[0]['chain'], 'A')


class TestSiteAnalyzer(TempDirTest):

    def _protonated(self):
        structure = io.load_pdb(CRN)
        io.save_pdb(structure, self.path('clean.pdb'),
                    select=cleaner.clean_structure(structure))
        protonator.add_hydrogens(self.path('clean.pdb'), self.path('prot.pdb'))
        return self.path('prot.pdb')

    def test_b13_hydrogens_excluded_from_pocket_geometry(self):
        """B13: pocket geometry is defined by heavy atoms."""
        source = self._protonated()
        total = sum(1 for l in pathlib.Path(source).read_text().splitlines()
                    if l.startswith(('ATOM', 'HETATM')))

        analyzer_obj = site_analyzer.BindingSiteAnalyzer(source)
        self.assertLess(len(analyzer_obj.atoms), total)
        self.assertTrue(all(a.element != 'H' for a in analyzer_obj.atoms))

    def test_b12_volume_tracks_grid_resolution(self):
        """B12: volume must use the grid actually used, not a hardcoded 1.5."""
        import numpy as np
        analyzer_obj = site_analyzer.BindingSiteAnalyzer(self._protonated())
        points = np.array([[float(i), 0.0, 0.0] for i in range(100)])

        fine = analyzer_obj._analyze_specific_site(points, 1.5)['volume']
        coarse = analyzer_obj._analyze_specific_site(points, 3.0)['volume']
        self.assertAlmostEqual(coarse / fine, 8.0, places=2)

    def test_b14_aromatics_survive_a_tight_cap(self):
        """B14: backbone N/O must yield to aromatic and sidechain features."""
        original = site_analyzer.MAX_PHARMACOPHORES
        site_analyzer.MAX_PHARMACOPHORES = 8
        try:
            analyzer_obj = site_analyzer.BindingSiteAnalyzer(self._protonated())
            sites = analyzer_obj.analyze()
            self.assertTrue(sites, 'no pockets found in 1CRN')

            aromatic = sum(1 for s in sites for p in s['pharmacophore_points']
                           if p['type'] == 'AROMATIC')
            self.assertGreater(aromatic, 0,
                               'aromatic features dropped by the feature cap')
            for site in sites:
                self.assertLessEqual(len(site['pharmacophore_points']), 8)
        finally:
            site_analyzer.MAX_PHARMACOPHORES = original

    def test_empty_structure_returns_no_pockets(self):
        source = _write(self.path('tiny.pdb'), _protein_lines()[:3])
        self.assertEqual(site_analyzer.BindingSiteAnalyzer(source).analyze(), [])


class TestExporter(TempDirTest):

    def test_b26_gromacs_export_does_not_mangle_paths(self):
        """B26: str.replace('.pdb') corrupts any path containing that text."""
        nested = os.path.join(self.tmp, 'my.pdb.data')
        os.makedirs(nested, exist_ok=True)
        source = _write(self.path('in.pdb'), _protein_lines())

        ok, result = exporter.export_structure(
            source, os.path.join(nested, 'out.pdb'), 'gromacs')

        self.assertTrue(ok)
        self.assertIn('my.pdb.data', result)
        self.assertTrue(os.path.isfile(result))

    def test_b27_pdbqt_receptor_has_no_torsion_tree(self):
        """B27: a receptor PDBQT must not carry ROOT/BRANCH records."""
        if not shutil.which('obabel'):
            self.skipTest('OpenBabel not on PATH')

        structure = io.load_pdb(CRN)
        io.save_pdb(structure, self.path('clean.pdb'),
                    select=cleaner.clean_structure(structure))
        protonator.add_hydrogens(self.path('clean.pdb'), self.path('prot.pdb'))

        ok, result = exporter.export_structure(
            self.path('prot.pdb'), self.path('rec.pdb'), 'vina')
        self.assertTrue(ok, f'export failed: {result}')

        text = pathlib.Path(result).read_text()
        self.assertNotIn('ROOT', text, 'receptor was written as a torsion tree')
        self.assertNotIn('BRANCH', text)
        atoms = [l for l in text.splitlines() if l.startswith(('ATOM', 'HETATM'))]
        self.assertTrue(atoms, 'receptor has no atom records')

    def test_missing_obabel_is_reported_not_raised(self):
        """A missing OpenBabel must return a message, not blow up the request."""
        import subprocess as _subprocess

        def explode(*args, **kwargs):
            raise FileNotFoundError('obabel')

        original = _subprocess.run
        exporter.subprocess.run = explode
        try:
            ok, message = exporter.convert_to_pdbqt('in.pdb', 'out.pdbqt')
        finally:
            exporter.subprocess.run = original

        self.assertFalse(ok)
        self.assertIn('OpenBabel', message)

    def test_obabel_nonzero_exit_is_reported(self):
        import subprocess as _subprocess

        class Result:
            returncode = 1
            stderr = 'something went wrong'

        original = _subprocess.run
        exporter.subprocess.run = lambda *a, **k: Result()
        try:
            ok, message = exporter.convert_to_pdbqt('in.pdb', 'out.pdbqt')
        finally:
            exporter.subprocess.run = original

        self.assertFalse(ok)
        self.assertIn('something went wrong', message)


if __name__ == '__main__':
    unittest.main(verbosity=2)
