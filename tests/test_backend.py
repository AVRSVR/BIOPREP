"""
Regression tests for the bugs catalogued in BACKEND_AUDIT.md.

Each test names the audit ID it guards so a future change that reintroduces
the bug fails with an obvious label.
"""

import json
import os
import pathlib
import shutil
import sys
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


class TestConectPreservation(TempDirTest):
    """Bond records must survive a save; Biopython drops them on parse."""

    LIGAND = [
        'HETATM  900  C1  LIG A 900      12.000  12.000  12.000  1.00 10.00           C  \n',
        'HETATM  901  C2  LIG A 900      13.500  12.000  12.000  1.00 10.00           C  \n',
        'HETATM  902  O1  LIG A 900      14.200  13.100  12.000  1.00 10.00           O  \n',
        'HETATM  903  N1  LIG A 900      14.200  10.900  12.000  1.00 10.00           N  \n',
    ]
    CONECT = [
        'CONECT  900  901\n',
        'CONECT  901  900  902  903\n',
        'CONECT  902  901\n',
        'CONECT  903  901\n',
    ]

    def _source(self):
        return _write(self.path('in.pdb'),
                      _protein_lines() + self.LIGAND + self.CONECT)

    def _bond_names(self, path):
        """Resolve CONECT records to ligand atom-name pairs."""
        by_serial, _ = io._serial_index(path)
        bonds = set()
        for line in pathlib.Path(path).read_text().splitlines():
            if not line.startswith('CONECT'):
                continue
            serials = list(io.conect_serials(line))
            for partner in serials[1:]:
                left, right = by_serial.get(serials[0]), by_serial.get(partner)
                if left and right and left[3] == 'LIG' and right[3] == 'LIG':
                    bonds.add(tuple(sorted([left[4], right[4]])))
        return bonds

    def test_conect_is_carried_across_a_save(self):
        source = self._source()
        structure = io.load_pdb(source)
        io.save_pdb(structure, self.path('out.pdb'), conect_source=source)

        self.assertEqual(self._bond_names(self.path('out.pdb')),
                         {('C1', 'C2'), ('C2', 'O1'), ('C2', 'N1')})

    def test_serials_are_remapped_not_copied(self):
        """Both files number atoms independently, so serials must be rewritten."""
        source = self._source()
        io.save_pdb(io.load_pdb(source), self.path('out.pdb'),
                    conect_source=source)

        _, out_by_key = io._serial_index(self.path('out.pdb'))
        new_serial = out_by_key[('A', '900', '', 'LIG', 'C1', '')]
        self.assertNotEqual(new_serial, 900,
                            'output still uses the input serial numbering')

        written = [l for l in pathlib.Path(self.path('out.pdb')).read_text().splitlines()
                   if l.startswith('CONECT')]
        self.assertTrue(written)
        for line in written:
            for serial in io.conect_serials(line):
                self.assertIn(serial, out_by_key.values(),
                              'CONECT references an atom that is not in the file')

    def test_bonds_to_removed_atoms_are_dropped(self):
        """Filtering a ligand out must not leave dangling bond records."""
        source = self._source()
        structure = io.load_pdb(source)
        select = cleaner.clean_structure(structure, remove_heteroatoms=['ALL'])
        io.save_pdb(structure, self.path('out.pdb'), select=select,
                    conect_source=source)

        text = pathlib.Path(self.path('out.pdb')).read_text()
        self.assertNotIn('LIG', text)
        self.assertEqual([l for l in text.splitlines() if l.startswith('CONECT')], [])

    def test_without_conect_source_nothing_is_written(self):
        source = self._source()
        io.save_pdb(io.load_pdb(source), self.path('out.pdb'))
        self.assertEqual(self._bond_names(self.path('out.pdb')), set())

    def test_conect_survives_the_full_pipeline(self):
        from bioprep.core.pipeline import PipelineSettings, prepare_structure

        source = self._source()
        workdir = self.path('work')
        os.makedirs(workdir, exist_ok=True)

        outcome = prepare_structure(source, workdir,
                                    PipelineSettings.from_mapping({}),
                                    original_filename='in.pdb')

        self.assertEqual(self._bond_names(outcome['viewer_path']),
                         {('C1', 'C2'), ('C2', 'O1'), ('C2', 'N1')},
                         'ligand connectivity lost somewhere in the pipeline')


class TestSeqresPreservation(TempDirTest):
    """
    SEQRES must survive a save, or PDBFixer cannot tell what is missing.

    Biopython drops it on parse and PDBIO writes none, so cleaning stripped it
    before the protonator ran and loop reconstruction silently rebuilt nothing.
    """

    def _gapped(self):
        with open(CRN) as fh:
            seqres = [l for l in fh if l.startswith('SEQRES')]
        kept = [l for l in _protein_lines() if not (20 <= int(l[22:26]) <= 24)]
        return _write(self.path('gap.pdb'), seqres + kept)

    @staticmethod
    def _seqres_count(path):
        return sum(1 for l in pathlib.Path(path).read_text().splitlines()
                   if l.startswith('SEQRES'))

    def test_seqres_carried_across_a_save(self):
        source = self._gapped()
        structure = io.load_pdb(source)
        io.save_pdb(structure, self.path('out.pdb'), source_pdb=source)

        self.assertGreater(self._seqres_count(self.path('out.pdb')), 0)

    def test_without_source_no_seqres_is_written(self):
        source = self._gapped()
        io.save_pdb(io.load_pdb(source), self.path('out.pdb'))
        self.assertEqual(self._seqres_count(self.path('out.pdb')), 0)

    def test_seqres_filtered_to_surviving_chains(self):
        """
        Carrying the sequence of a filtered-out chain would make PDBFixer treat
        that whole chain as missing and try to build it from nothing.
        """
        with open(CRN) as fh:
            seqres = [l for l in fh if l.startswith('SEQRES')]
        seqres_b = [l[:11] + 'B' + l[12:] for l in seqres]
        lines = _protein_lines()
        source = _write(self.path('two.pdb'),
                        seqres + seqres_b + lines
                        + [l[:21] + 'B' + l[22:] for l in lines])

        structure = io.load_pdb(source)
        io.save_pdb(structure, self.path('out.pdb'),
                    select=cleaner.clean_structure(structure, target_chains=['A']),
                    source_pdb=source)

        carried = [l for l in pathlib.Path(self.path('out.pdb')).read_text().splitlines()
                   if l.startswith('SEQRES')]
        self.assertTrue(carried)
        self.assertEqual({l[11] for l in carried}, {'A'})

    def test_loop_reconstruction_works_through_the_pipeline(self):
        from bioprep.core.pipeline import PipelineSettings, prepare_structure

        source = self._gapped()
        workdir = self.path('work')
        os.makedirs(workdir, exist_ok=True)

        outcome = prepare_structure(
            source, workdir,
            PipelineSettings.from_mapping({'reconstruct_loops': True}),
            original_filename='gap.pdb')

        self.assertEqual(
            outcome['report']['protonation']['loops_reconstructed'], 5,
            'loops were not rebuilt - SEQRES probably did not survive cleaning')

        written = {int(l[22:26]) for l in
                   pathlib.Path(outcome['viewer_path']).read_text().splitlines()
                   if l.startswith('ATOM')}
        self.assertLessEqual({20, 21, 22, 23, 24}, written)

    def test_loops_not_rebuilt_when_flag_is_off(self):
        from bioprep.core.pipeline import PipelineSettings, prepare_structure

        source = self._gapped()
        workdir = self.path('work2')
        os.makedirs(workdir, exist_ok=True)

        outcome = prepare_structure(source, workdir,
                                    PipelineSettings.from_mapping({}),
                                    original_filename='gap.pdb')
        self.assertEqual(
            outcome['report']['protonation']['loops_reconstructed'], 0)


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


class TestWaterRemoval(TempDirTest):
    """Feature 5: water removal, whatever the water is called."""

    SPELLINGS = ('HOH', 'WAT', 'H2O', 'TIP', 'SOL', 'DOD')

    def _with_waters(self):
        def het(serial, name, resn, resseq, x, elem):
            return (f'HETATM{serial:>5} {name:<4} {resn:>3} A{resseq:>4}    '
                    f'{x:8.3f}  12.000  12.000  1.00 10.00          {elem:>2}  \n')

        rows = [het(9000 + i, 'O', name, 700 + i, 40.0 + i * 4, 'O')
                for i, name in enumerate(self.SPELLINGS)]
        rows.append(het(9100, 'C1', 'LIG', 800, 12.0, 'C'))
        rows.append(het(9101, 'ZN', ' ZN', 801, 20.0, 'ZN'))
        return _write(self.path('in.pdb'), _protein_lines() + rows)

    def _heteroatom_names(self, path):
        return sorted({l[17:20].strip()
                       for l in pathlib.Path(path).read_text().splitlines()
                       if l.startswith('HETATM')})

    def test_every_water_spelling_is_removed(self):
        """Biopython tags only HOH and WAT as 'W'; the rest arrive as H_xxx."""
        source = self._with_waters()
        structure = io.load_pdb(source)
        io.save_pdb(structure, self.path('out.pdb'),
                    select=cleaner.clean_structure(structure, remove_water=True))

        self.assertEqual(self._heteroatom_names(self.path('out.pdb')),
                         ['LIG', 'ZN'])

    def test_waters_are_kept_when_not_asked_to_remove(self):
        source = self._with_waters()
        structure = io.load_pdb(source)
        io.save_pdb(structure, self.path('out.pdb'),
                    select=cleaner.clean_structure(structure, remove_water=False))

        kept = self._heteroatom_names(self.path('out.pdb'))
        for name in self.SPELLINGS:
            self.assertIn(name, kept)

    def test_analyzer_counts_every_spelling_as_water(self):
        metadata = analyzer.analyze_structure(io.load_pdb(self._with_waters()))

        self.assertEqual(metadata['water_count'], len(self.SPELLINGS))
        self.assertEqual(metadata['heteroatoms'], ['LIG', 'ZN'],
                         'waters were double-counted as ligands')
        self.assertEqual(metadata['atom_breakdown']['water'], len(self.SPELLINGS))

    def test_lowercase_water_name_is_recognised(self):
        row = ('HETATM 9200  O   hoh A 700      40.000  12.000  12.000'
               '  1.00 10.00           O  \n')
        source = _write(self.path('lc.pdb'), _protein_lines() + [row])
        structure = io.load_pdb(source)
        io.save_pdb(structure, self.path('out.pdb'),
                    select=cleaner.clean_structure(structure, remove_water=True))

        self.assertEqual(self._heteroatom_names(self.path('out.pdb')), [])


class TestCleanerRules(TempDirTest):
    """Features 6-11: structural waters, chain filtering, heteroatom rules."""

    @staticmethod
    def _het(serial, name, resn, chain, resseq, x, y, z, elem):
        return (f'HETATM{serial:>5} {name:<4} {resn:>3} {chain:1}{resseq:>4}    '
                f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00          {elem:>2}  \n')

    def _near_far_waters(self, oxygen_name):
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        return _write(self.path('w.pdb'), lines + [
            self._het(9100, oxygen_name, 'HOH', 'A', 750, x + 2.8, y, z, 'O'),
            self._het(9101, oxygen_name, 'HOH', 'A', 751, x + 80, y, z, 'O'),
        ])

    def _kept_resseqs(self, path):
        return sorted(int(l[22:26]) for l in pathlib.Path(path).read_text().splitlines()
                      if l.startswith('HETATM'))

    def _het_names(self, path):
        return sorted({l[17:20].strip()
                       for l in pathlib.Path(path).read_text().splitlines()
                       if l.startswith('HETATM')})

    def test_feature6_water_oxygen_naming_conventions(self):
        """
        Crystallographic files use O, GROMACS uses OW, CHARMM uses OH2. Matching
        a fixed name list meant CHARMM waters had no detectable oxygen, so
        structural-water preservation kept nothing at all.
        """
        for oxygen in ('O', 'OW', 'OH2', 'O1'):
            with self.subTest(oxygen=oxygen):
                source = self._near_far_waters(oxygen)
                structure = io.load_pdb(source)
                io.save_pdb(structure, self.path('out.pdb'),
                            select=cleaner.clean_structure(
                                structure, remove_water=True,
                                keep_structural_waters=True))
                self.assertEqual(self._kept_resseqs(self.path('out.pdb')), [750])

    def test_feature6_distant_water_is_removed(self):
        source = self._near_far_waters('O')
        structure = io.load_pdb(source)
        select = cleaner.clean_structure(structure, remove_water=True,
                                         keep_structural_waters=True)
        self.assertEqual(len(select.structural_waters), 1)

    def test_feature7_chain_selection_keeps_only_that_chain(self):
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        both = lines + [l[:21] + 'B' + l[22:] for l in lines] + [
            self._het(9200, 'C1', 'LIG', 'A', 800, x, y, z, 'C'),
            self._het(9201, 'C1', 'LIG', 'B', 800, x, y, z, 'C'),
        ]
        source = _write(self.path('ch.pdb'), both)

        structure = io.load_pdb(source)
        io.save_pdb(structure, self.path('out.pdb'),
                    select=cleaner.clean_structure(structure, target_chains=['A']))

        text = pathlib.Path(self.path('out.pdb')).read_text().splitlines()
        chains = {l[21] for l in text if l.startswith(('ATOM', 'HETATM'))}
        self.assertEqual(chains, {'A'})

    def _mixed(self):
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        return _write(self.path('m.pdb'), lines + [
            self._het(9300, 'C1', 'LIG', 'A', 800, x + 10, y, z, 'C'),
            self._het(9301, 'S', 'SO4', 'A', 801, x + 14, y, z, 'S'),
            self._het(9302, 'ZN', ' ZN', 'A', 802, x + 18, y, z, 'ZN'),
        ])

    def _clean(self, **kwargs):
        source = self._mixed()
        structure = io.load_pdb(source)
        io.save_pdb(structure, self.path('out.pdb'),
                    select=cleaner.clean_structure(structure, **kwargs))
        return self._het_names(self.path('out.pdb'))

    def test_feature8_selective_removal(self):
        self.assertEqual(self._clean(remove_heteroatoms=['SO4']), ['LIG', 'ZN'])

    def test_feature9_remove_all_heteroatoms(self):
        self.assertEqual(self._clean(remove_heteroatoms=['ALL']), [])

    def test_feature10_unlisted_heteroatoms_are_kept(self):
        self.assertEqual(self._clean(remove_heteroatoms=[]), ['LIG', 'SO4', 'ZN'])

    def test_residue_names_are_matched_case_insensitively(self):
        """A lower-case entry used to silently fail to match."""
        self.assertEqual(self._clean(remove_heteroatoms=['so4']), ['LIG', 'ZN'])
        self.assertEqual(self._clean(remove_heteroatoms=['  Zn ']), ['LIG', 'SO4'])
        self.assertEqual(self._clean(remove_heteroatoms=['all']), [])

    def test_lowercase_protect_ligands_still_protects(self):
        """The dangerous case: the ligand was deleted without a word."""
        self.assertEqual(
            self._clean(remove_heteroatoms=['ALL'], protect_ligands=['lig']),
            ['LIG'])

    def test_feature11_is_a_biopython_select_subclass(self):
        from Bio.PDB import Select

        select = cleaner.clean_structure(io.load_pdb(CRN))
        self.assertIsInstance(select, Select)
        self.assertIsInstance(select, cleaner.BioPrepSelect)
        for method in ('accept_model', 'accept_chain', 'accept_residue', 'accept_atom'):
            self.assertTrue(hasattr(select, method))


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


class TestProtonationFeatures(TempDirTest):
    """Features 12-22: pH, ligand separation, repairs, status report."""

    @staticmethod
    def _het(serial, name, resn, chain, resseq, x, y, z, elem):
        return (f'HETATM{serial:>5} {name:<4} {resn:>3} {chain:1}{resseq:>4}    '
                f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00          {elem:>2}  \n')

    @staticmethod
    def _hydrogens(path):
        return sum(1 for l in pathlib.Path(path).read_text().splitlines()
                   if l.startswith(('ATOM', 'HETATM')) and l[76:78].strip() == 'H')

    def _with_ligand(self):
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        return _write(self.path('in.pdb'), lines + [
            self._het(900, 'C1', 'LIG', 'A', 900, x + 6, y, z, 'C'),
            self._het(901, 'C2', 'LIG', 'A', 900, x + 7.5, y, z, 'C'),
            self._het(902, 'O1', 'LIG', 'A', 900, x + 8.2, y + 1.1, z, 'O'),
        ])

    def test_feature12_hydrogen_count_depends_on_ph(self):
        source = _write(self.path('p.pdb'), _protein_lines())
        counts = {}
        for ph in (1.0, 13.0):
            out = self.path(f'ph{ph}.pdb')
            protonator.add_hydrogens(source, out, ph=ph)
            counts[ph] = self._hydrogens(out)

        self.assertGreater(counts[1.0], 0)
        self.assertGreater(counts[1.0], counts[13.0],
                           'acidic pH should protonate more than basic')

    def test_feature12_default_ph_is_7_4(self):
        source = _write(self.path('p.pdb'), _protein_lines())
        result = protonator.add_hydrogens(source, self.path('o.pdb'))
        self.assertEqual(result['ph'], 7.4)

    def test_feature13_14_ligand_is_untouched(self):
        source = self._with_ligand()
        result = protonator.add_hydrogens(source, self.path('out.pdb'))

        original = [l for l in pathlib.Path(source).read_text().splitlines()
                    if l[17:20].strip() == 'LIG']
        written = [l for l in pathlib.Path(self.path('out.pdb')).read_text().splitlines()
                   if l[17:20].strip() == 'LIG']

        self.assertEqual(result['ligands_preserved'], ['LIG'])
        self.assertEqual(len(written), len(original))
        for before, after in zip(original, written):
            self.assertEqual(before[30:54], after[30:54], 'ligand coordinates moved')
        self.assertTrue(all(l[76:78].strip() != 'H' for l in written),
                        'hydrogens were added to the ligand')

    def test_feature13_waters_are_not_held_out_as_ligands(self):
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        source = _write(self.path('w.pdb'), lines + [
            self._het(950, 'O', 'HOH', 'A', 700, x + 2.8, y, z, 'O')])

        result = protonator.add_hydrogens(source, self.path('out.pdb'))
        self.assertEqual(result['ligands_preserved'], [])

    def test_feature19_merge_produces_unique_ascending_serials(self):
        source = self._with_ligand()
        protonator.add_hydrogens(source, self.path('out.pdb'))

        serials = [int(l[6:11]) for l in
                   pathlib.Path(self.path('out.pdb')).read_text().splitlines()
                   if l[:6] in ('ATOM  ', 'HETATM')]
        self.assertEqual(serials, sorted(serials))
        self.assertEqual(len(serials), len(set(serials)))

    def test_feature16_loops_only_rebuilt_when_asked(self):
        with open(CRN) as fh:
            seqres = [l for l in fh if l.startswith('SEQRES')]
        gapped = [l for l in _protein_lines() if not (20 <= int(l[22:26]) <= 24)]
        source = _write(self.path('gap.pdb'), seqres + gapped)

        off = protonator.add_hydrogens(source, self.path('off.pdb'),
                                       reconstruct_loops=False)
        on = protonator.add_hydrogens(source, self.path('on.pdb'),
                                      reconstruct_loops=True)

        self.assertEqual(off['loops_reconstructed'], 0)
        self.assertEqual(on['loops_reconstructed'], 5)

        def resseqs(path):
            return {int(l[22:26]) for l in
                    pathlib.Path(path).read_text().splitlines()
                    if l.startswith('ATOM')}

        self.assertFalse(resseqs(self.path('off.pdb')) & {20, 21, 22, 23, 24})
        self.assertLessEqual({20, 21, 22, 23, 24}, resseqs(self.path('on.pdb')))

    def test_feature18_hydrogen_failure_is_caught_not_raised(self):
        """The documented try/except: warn, keep the structure, do not crash."""
        import pdbfixer

        source = _write(self.path('p.pdb'), _protein_lines())
        original = pdbfixer.PDBFixer.addMissingHydrogens

        def explode(self, ph=7.0, *args, **kwargs):
            raise RuntimeError('simulated template failure')

        pdbfixer.PDBFixer.addMissingHydrogens = explode
        try:
            result = protonator.add_hydrogens(source, self.path('out.pdb'))
        finally:
            pdbfixer.PDBFixer.addMissingHydrogens = original

        self.assertFalse(result['hydrogens_added'])
        self.assertTrue(any('simulated' in w for w in result['warnings']))
        self.assertTrue(os.path.isfile(self.path('out.pdb')),
                        'structure was lost when protonation failed')

    def test_feature20_status_report_explains_why(self):
        source = self._with_ligand()
        result = protonator.add_hydrogens(source, self.path('out.pdb'))

        self.assertTrue(result['ligand_status'])
        entry = result['ligand_status'][0]
        self.assertEqual(entry['residue'], 'LIG')
        self.assertEqual(entry['atoms'], 3)
        self.assertEqual(entry['action'], 'preserved')
        self.assertIn('protonation', entry['reason'].lower())

    def test_feature20_reason_reaches_the_text_report(self):
        from bioprep.core.pipeline import PipelineSettings, prepare_structure

        source = self._with_ligand()
        workdir = self.path('work')
        os.makedirs(workdir, exist_ok=True)
        outcome = prepare_structure(source, workdir,
                                    PipelineSettings.from_mapping({}),
                                    original_filename='in.pdb')

        text = reporter.report_to_text(outcome['report'])
        self.assertIn('LIG', text)
        self.assertIn('Blind protonation', text)

    def test_feature21_no_temp_files_left_behind(self):
        import glob as _glob

        def leaked():
            root = tempfile.gettempdir()
            return set(_glob.glob(os.path.join(root, '*_protein.pdb'))) | \
                set(_glob.glob(os.path.join(root, '*_protonated.pdb')))

        source = self._with_ligand()
        before = leaked()
        protonator.add_hydrogens(source, self.path('out.pdb'))
        self.assertEqual(leaked(), before)

        # and on the error path
        ligand_only = _write(self.path('lig.pdb'), [
            self._het(900, 'C1', 'LIG', 'A', 900, 1, 1, 1, 'C')])
        with self.assertRaises(ValueError):
            protonator.add_hydrogens(ligand_only, self.path('x.pdb'))
        self.assertEqual(leaked(), before)

    def test_feature22_chain_ids_and_numbering_preserved(self):
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        source = _write(self.path('multi.pdb'),
                        lines + [l[:21] + 'B' + l[22:] for l in lines] +
                        [self._het(960, 'C1', 'LIG', 'C', 900, x + 6, y, z, 'C')])

        protonator.add_hydrogens(source, self.path('out.pdb'))
        written = pathlib.Path(self.path('out.pdb')).read_text().splitlines()

        chains = sorted({l[21] for l in written if l[:6] in ('ATOM  ', 'HETATM')})
        self.assertEqual(chains, ['A', 'B', 'C'])

        numbers = [int(l[22:26]) for l in written
                   if l.startswith('ATOM') and l[21] == 'A']
        self.assertEqual((min(numbers), max(numbers)), (1, 46))


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


class TestMinimizerFeatures(TempDirTest):
    """Features 29-48 that the earlier rewrite did not cover."""

    def _prepared(self, extra=()):
        source = _write(self.path('raw.pdb'), _protein_lines() + list(extra))
        structure = io.load_pdb(source)
        io.save_pdb(structure, self.path('clean.pdb'),
                    select=cleaner.clean_structure(structure, remove_water=False))
        protonator.add_hydrogens(self.path('clean.pdb'), self.path('prot.pdb'))
        return self.path('prot.pdb')

    def test_feature29_charmm36_force_field_runs(self):
        stats = minimizer.minimize_structure(
            self._prepared(), self.path('out.pdb'), force_field='charmm36')
        self.assertEqual(stats['status'], 'full', stats.get('error'))
        self.assertEqual(stats['force_field'], 'charmm36')

    def test_feature33_terminal_repair_is_available(self):
        """createSystem failure should trigger a Modeller capping retry."""
        import inspect
        source = inspect.getsource(minimizer._create_system)
        self.assertIn('addHydrogens(forcefield=', source)
        self.assertIn('Modeller', source)

    def test_feature37_excluded_residue_is_held_fixed(self):
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        ligand = (f'HETATM  900  C1  LIG A 900    {x + 6:8.3f}{y:8.3f}{z:8.3f}'
                  '  1.00 10.00           C  \n')
        source = self._prepared([ligand])

        stats = minimizer.minimize_structure(source, self.path('out.pdb'))
        self.assertEqual(stats['status'], 'partial')

        def ligand_coords(path):
            # ATOM/HETATM only: PDBFile also emits a TER carrying the resname.
            return [l[30:54] for l in pathlib.Path(path).read_text().splitlines()
                    if l[:6] in ('ATOM  ', 'HETATM') and l[17:20].strip() == 'LIG']

        self.assertEqual(ligand_coords(source), ligand_coords(self.path('out.pdb')))

    def test_feature41_convergence_follows_the_documented_rule(self):
        stats = minimizer.minimize_structure(self._prepared(), self.path('out.pdb'))
        delta = stats['delta_energy_kJ_mol']
        expected = (delta < -1.0) or (abs(delta) < 1.0
                                      and stats['energy_after_kJ_mol'] < 0)
        self.assertEqual(stats['converged'], expected)
        self.assertTrue(stats['converged'], f'delta={delta}')

    def test_feature41_failed_run_is_not_converged(self):
        junk = _write(self.path('junk.pdb'), ['REMARK nothing\n'])
        stats = minimizer.minimize_structure(junk, self.path('out.pdb'))
        self.assertEqual(stats['status'], 'failed')
        self.assertFalse(stats['converged'])

    def test_feature43_44_caps_are_reported(self):
        stats = minimizer.minimize_structure(self._prepared(), self.path('out.pdb'))
        self.assertEqual(stats['iterations_max'], 1000)
        self.assertEqual(stats['energy_tolerance_kJ_mol_nm'], 10.0)

    def test_feature46_cpu_platform_is_preferred_with_fallback(self):
        import inspect
        source = inspect.getsource(minimizer._make_simulation)
        self.assertIn("getPlatformByName('CPU')", source)
        self.assertIn('except Exception', source)

    def test_feature38_residue_classification_sets(self):
        from bioprep.core import residues

        self.assertLessEqual({'HID', 'HIE', 'HIP', 'CYX', 'ASH', 'GLH'},
                             residues.AMINO_ACIDS)
        self.assertLessEqual({'DA', 'DC', 'DG', 'DT'}, residues.NUCLEIC_ACIDS)
        self.assertEqual(residues.FORCE_FIELD_SAFE,
                         residues.AMINO_ACIDS | residues.NUCLEIC_ACIDS
                         | residues.WATER)
        self.assertFalse(residues.is_force_field_safe('LIG'))


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


class TestReportContent(unittest.TestCase):
    """Features 74-80: the report is an audit trail, so it must be accurate."""

    PROTONATION = {
        'hydrogens_added': True, 'ph': 7.4, 'ligands_preserved': ['BTN'],
        'ligand_status': [{'residue': 'BTN', 'atoms': 16,
                           'action': 'preserved', 'reason': 'no template'}],
        'nonstandard_replaced': ['MSE'], 'loops_reconstructed': 2,
    }
    MINIMIZATION = {
        'status': 'partial', 'force_field': 'amber14',
        'energy_before_kJ_mol': 900.0, 'energy_after_kJ_mol': -5000.0,
        'delta_energy_kJ_mol': -5900.0, 'converged': True,
        'excluded_residues': ['BTN'], 'error': None,
    }

    def _report(self, **overrides):
        kwargs = dict(
            filename='1crn.pdb', chains_detected=['A', 'B'],
            chains_retained=['A'], waters_removed=12, waters_retained=3,
            heteroatoms_removed=['SO4'], heteroatoms_retained=['BTN'],
            ph_used=7.4, missing_residues=[{'chain': 'A', 'residue': 'GLY'}],
            atoms_before=327, atoms_after=642,
            protonation=self.PROTONATION, minimization_stats=self.MINIMIZATION,
            docking_export={'target': 'vina', 'succeeded': True},
            processing_time_s=3.14159, warnings=['something to note'])
        kwargs.update(overrides)
        return reporter.build_report(**kwargs)

    def test_feature74_structure_and_atom_delta(self):
        report = self._report()
        self.assertEqual(report['chains'],
                         {'detected': ['A', 'B'], 'retained': ['A']})
        self.assertEqual(report['atom_counts']['delta'], 315)
        self.assertTrue(report['generated_at'].endswith('UTC'))

    def test_feature76_failed_export_does_not_claim_a_target(self):
        report = self._report(docking_export={'target': 'vina',
                                              'succeeded': False,
                                              'error': 'obabel missing'})
        self.assertNotIn('docking_target', report)
        self.assertFalse(report['docking_export']['succeeded'])

    def test_feature77_processing_time_rounded(self):
        self.assertEqual(self._report()['processing_time_seconds'], 3.14)

    def test_feature78_ligand_status_reaches_the_report(self):
        report = self._report()
        self.assertEqual(report['protonation']['ligand_status'],
                         self.PROTONATION['ligand_status'])

    def test_feature79_text_report_sections(self):
        text = reporter.report_to_text(self._report())
        for heading in ('PREPARATION REPORT', 'CHAINS', 'WATER MOLECULES',
                        'HETEROATOMS', 'PROTONATION', 'STRUCTURE',
                        'ENERGY MINIMIZATION', 'DOCKING EXPORT', 'WARNINGS'):
            self.assertIn(heading, text)
        self.assertIn('+315', text)
        self.assertIn('no template', text)
        self.assertIn('something to note', text)

    def test_feature79_degenerate_reports_render(self):
        self.assertIn('BIOPREP', reporter.report_to_text({}))
        self.assertIn('BIOPREP', reporter.report_to_text(None))
        self.assertIn('x.pdb', reporter.report_to_text({'input_file': 'x.pdb'}))

    def test_feature80_atom_counter(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'm.pdb')
            with open(path, 'w') as fh:
                fh.write('REMARK header\n')
                fh.write('ATOM      1  N   ALA A   1       1.000   2.000'
                         '   3.000  1.00  0.00           N  \n')
                fh.write('HETATM    2 ZN    ZN A   2       1.000   2.000'
                         '   3.000  1.00  0.00          ZN  \n')
                fh.write('TER\nEND\n')
            self.assertEqual(reporter.count_atoms_in_pdb(path), 2)
            self.assertEqual(
                reporter.count_atoms_in_pdb(os.path.join(tmp, 'nope.pdb')), 0)


class TestExporterFormats(TempDirTest):
    """Features 70-72: the three export targets."""

    def _prepared(self):
        source = _write(self.path('raw.pdb'), _protein_lines())
        structure = io.load_pdb(source)
        io.save_pdb(structure, self.path('clean.pdb'),
                    select=cleaner.clean_structure(structure))
        protonator.add_hydrogens(self.path('clean.pdb'), self.path('prot.pdb'))
        return self.path('prot.pdb')

    def test_feature71_gromacs_suffix_and_content(self):
        source = self._prepared()
        ok, result = exporter.export_structure(source, self.path('g.pdb'), 'gromacs')

        self.assertTrue(ok)
        self.assertTrue(str(result).endswith('_gromacs.pdb'))
        self.assertEqual(reporter.count_atoms_in_pdb(result),
                         reporter.count_atoms_in_pdb(source))

    def test_feature72_unknown_format_falls_back_to_pdb(self):
        source = self._prepared()
        ok, result = exporter.export_structure(source, self.path('u.pdb'), 'nonsense')

        self.assertTrue(ok)
        self.assertTrue(str(result).endswith('.pdb'))
        self.assertEqual(reporter.count_atoms_in_pdb(result),
                         reporter.count_atoms_in_pdb(source))

    def test_feature72_export_onto_own_path_is_a_noop(self):
        source = self._prepared()
        before = reporter.count_atoms_in_pdb(source)
        ok, _ = exporter.export_structure(source, source, 'pdb')

        self.assertTrue(ok)
        self.assertEqual(reporter.count_atoms_in_pdb(source), before)

    def test_feature73_exit_zero_without_output_is_a_failure(self):
        """OpenBabel can return 0 having written nothing."""
        import subprocess as _subprocess

        class Empty:
            returncode = 0
            stderr = ''

        original = _subprocess.run
        exporter.subprocess.run = lambda *a, **k: Empty()
        try:
            ok, message = exporter.convert_to_pdbqt(
                'in.pdb', self.path('never.pdbqt'))
        finally:
            exporter.subprocess.run = original

        self.assertFalse(ok)
        self.assertIn('no output', message.lower())

    def test_feature73_timeout_is_reported(self):
        import subprocess as _subprocess

        original = _subprocess.run

        def timeout(*args, **kwargs):
            raise _subprocess.TimeoutExpired('obabel', 120)

        exporter.subprocess.run = timeout
        try:
            ok, message = exporter.convert_to_pdbqt('in.pdb', self.path('x.pdbqt'))
        finally:
            exporter.subprocess.run = original

        self.assertFalse(ok)
        self.assertIn('timed out', message.lower())


class TestCli(TempDirTest):
    """Features 111-112: the command-line interface."""

    def _run(self, *argv):
        """Call main() in-process with a patched argv; returns (code, stdout)."""
        import contextlib
        import io as _stdio
        from bioprep import cli

        buffer = _stdio.StringIO()
        original = sys.argv
        sys.argv = ['bioprep', *argv]
        try:
            with contextlib.redirect_stdout(buffer):
                code = cli.main()
        finally:
            sys.argv = original
        return code, buffer.getvalue()

    def _rich(self):
        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)
        # Chain B is offset; superimposed chains give an infinite start energy.
        chain_b = [l[:21] + 'B' + l[22:30] + f'{float(l[30:38]) + 60.0:8.3f}' + l[38:]
                   for l in lines]
        het = (f'HETATM  900  C1  LIG A 900    {x + 6:8.3f}{y:8.3f}{z:8.3f}'
               '  1.00 10.00           C  \n')
        water = (f'HETATM  902  O   HOH A 700    {x + 2.8:8.3f}{y:8.3f}{z:8.3f}'
                 '  1.00 10.00           O  \n')
        return _write(self.path('rich.pdb'), lines + chain_b + [het, water])

    def _resnames(self, path):
        return {l[17:20].strip() for l in pathlib.Path(path).read_text().splitlines()
                if l.startswith('HETATM')}

    def test_feature112_basic_run(self):
        code, out = self._run('--input', self._rich(), '--output', self.path('o.pdb'))
        self.assertEqual(code, 0, out)
        self.assertIn('Cleaning structure', out)
        self.assertIn('Adding hydrogens', out)
        self.assertTrue(os.path.isfile(self.path('o.pdb')))

    def test_feature112_water_removed_ligand_kept_by_default(self):
        self._run('--input', self._rich(), '--output', self.path('o.pdb'))
        names = self._resnames(self.path('o.pdb'))
        self.assertNotIn('HOH', names)
        self.assertIn('LIG', names)

    def test_feature112_chain_flag_is_repeatable(self):
        source = self._rich()
        self._run('--input', source, '--output', self.path('a.pdb'), '--chain', 'A')
        self._run('--input', source, '--output', self.path('ab.pdb'),
                  '--chain', 'A', '--chain', 'B')

        def chains(path):
            return {l[21] for l in pathlib.Path(path).read_text().splitlines()
                    if l[:6] in ('ATOM  ', 'HETATM')}

        self.assertEqual(chains(self.path('a.pdb')), {'A'})
        self.assertEqual(chains(self.path('ab.pdb')), {'A', 'B'})

    def test_feature112_protect_outranks_remove_all_case_insensitively(self):
        self._run('--input', self._rich(), '--output', self.path('o.pdb'),
                  '--remove-het', 'ALL', '--protect', 'lig')
        self.assertEqual(self._resnames(self.path('o.pdb')), {'LIG'})

    def test_feature112_minimize_reports_energies(self):
        code, out = self._run('--input', self._rich(), '--output', self.path('o.pdb'),
                              '--minimize')
        self.assertEqual(code, 0, out)
        self.assertIn('kJ/mol', out)
        self.assertIn('status:', out)

    def test_feature111_unreadable_input_exits_nonzero(self):
        code, out = self._run('--input', self.path('nope.pdb'),
                              '--output', self.path('o.pdb'))
        self.assertEqual(code, 1)
        self.assertIn('Error loading', out)

    def test_feature112_leaves_no_working_directory_behind(self):
        import glob as _glob

        def leaked():
            root = tempfile.gettempdir()
            return set(_glob.glob(os.path.join(root, 'bioprep_cli_*'))) | \
                set(_glob.glob(os.path.join(root, 'bioprep_cif_*')))

        before = leaked()
        self._run('--input', self._rich(), '--output', self.path('o.pdb'))
        self.assertEqual(leaked(), before)

        # and on the failure path
        self._run('--input', self.path('nope.pdb'), '--output', self.path('x.pdb'))
        self.assertEqual(leaked(), before)


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


class TestAnalyzerFeatures(TempDirTest):
    """Features 23-28: chains, waters, ligands, atom counts, gaps."""

    @staticmethod
    def _het(serial, name, resn, chain, resseq, x, elem):
        return (f'HETATM{serial:>5} {name:<4} {resn:>3} {chain:1}{resseq:>4}    '
                f'{x:8.3f}  12.000  12.000  1.00 10.00          {elem:>2}  \n')

    def _seqres(self):
        with open(CRN) as fh:
            return [l for l in fh if l.startswith('SEQRES')]

    def test_feature23_lists_every_chain_once(self):
        lines = _protein_lines()
        source = _write(self.path('c.pdb'),
                        lines + [l[:21] + 'B' + l[22:] for l in lines]
                        + [l[:21] + 'C' + l[22:] for l in lines])
        self.assertEqual(analyzer.analyze_structure(io.load_pdb(source))['chains'],
                         ['A', 'B', 'C'])

    def test_feature25_heteroatom_names_are_unique_and_exclude_water(self):
        rows = [self._het(9100, 'C1', 'LIG', 'A', 800, 12.0, 'C'),
                self._het(9101, 'C2', 'LIG', 'A', 800, 13.0, 'C'),
                self._het(9102, 'S', 'SO4', 'A', 801, 16.0, 'S'),
                self._het(9103, 'O', 'HOH', 'A', 700, 40.0, 'O')]
        source = _write(self.path('h.pdb'), _protein_lines() + rows)

        metadata = analyzer.analyze_structure(io.load_pdb(source))
        self.assertEqual(metadata['heteroatoms'], ['LIG', 'SO4'])

    def test_feature26_atom_breakdown_sums_to_total(self):
        rows = [self._het(9100, 'C1', 'LIG', 'A', 800, 12.0, 'C'),
                self._het(9103, 'O', 'HOH', 'A', 700, 40.0, 'O')]
        source = _write(self.path('h.pdb'), _protein_lines() + rows)

        metadata = analyzer.analyze_structure(io.load_pdb(source))
        self.assertEqual(sum(metadata['atom_breakdown'].values()),
                         metadata['atoms_total'])
        self.assertEqual(metadata['atom_breakdown'],
                         {'protein': len(_protein_lines()), 'water': 1,
                          'heteroatom': 1})

    def test_feature27_gap_reports_chain_from_to_and_size(self):
        gapped = [l for l in _protein_lines() if not (20 <= int(l[22:26]) <= 24)]
        source = _write(self.path('g.pdb'), gapped)

        gaps = analyzer.analyze_structure(io.load_pdb(source))['sequence_gaps']
        self.assertEqual(gaps, [{'chain': 'A', 'from': 19, 'to': 25,
                                 'missing_count': 5}])

    def test_feature27_no_false_positives(self):
        """Insertion codes, interleaved ligands and negative numbering."""
        lines = _protein_lines()

        with_insertions = []
        for line in lines:
            with_insertions.append(line)
            if int(line[22:26]) == 5:
                with_insertions.append(line[:26] + 'A' + line[27:])
        source = _write(self.path('ins.pdb'), with_insertions)
        self.assertEqual(
            analyzer.analyze_structure(io.load_pdb(source))['sequence_gaps'], [])

        interleaved = list(lines)
        interleaved.insert(10, self._het(9500, 'C1', 'LIG', 'A', 500, 50.0, 'C'))
        source = _write(self.path('int.pdb'), interleaved)
        self.assertEqual(
            analyzer.analyze_structure(io.load_pdb(source))['sequence_gaps'], [])

        shifted = [l[:22] + f'{int(l[22:26]) - 3:>4}' + l[26:] for l in lines]
        source = _write(self.path('neg.pdb'), shifted)
        self.assertEqual(
            analyzer.analyze_structure(io.load_pdb(source))['sequence_gaps'], [])

    def test_feature27_ensemble_does_not_duplicate_gaps(self):
        gapped = [l for l in _protein_lines() if not (20 <= int(l[22:26]) <= 24)]
        ensemble = []
        for index in range(1, 4):
            ensemble.append(f'MODEL     {index:>4}\n')
            ensemble.extend(gapped)
            ensemble.append('ENDMDL\n')
        source = _write(self.path('e.pdb'), ensemble)

        self.assertEqual(
            len(analyzer.analyze_structure(io.load_pdb(source))['sequence_gaps']), 1)

    def test_feature28_reports_chain_id_not_insert_position(self):
        """
        PDBFixer keys missingResidues by (chain index, insert position).
        Reading that as (model, chain_id) reported an integer offset - for a
        gap in the second chain it said chain 19 instead of chain B.
        """
        lines = _protein_lines()
        seqres = self._seqres()
        seqres_b = [l[:11] + 'B' + l[12:] for l in seqres]
        chain_b = [l[:21] + 'B' + l[22:] for l in lines
                   if not (20 <= int(l[22:26]) <= 24)]

        source = _write(self.path('b.pdb'), seqres + seqres_b + lines + chain_b)
        missing = analyzer.detect_missing_residues(source)

        self.assertEqual(len(missing), 5)
        self.assertEqual({m['chain'] for m in missing}, {'B'})
        for entry in missing:
            self.assertIsInstance(entry['residue'], str)
            self.assertTrue(entry['residue'].isalpha())

    def test_feature28_gaps_in_both_chains_attributed_separately(self):
        lines = _protein_lines()
        seqres = self._seqres()
        seqres_b = [l[:11] + 'B' + l[12:] for l in seqres]
        chain_a = [l for l in lines if not (5 <= int(l[22:26]) <= 7)]
        chain_b = [l[:21] + 'B' + l[22:] for l in lines
                   if not (20 <= int(l[22:26]) <= 24)]

        source = _write(self.path('ab.pdb'),
                        seqres + seqres_b + chain_a + chain_b)
        missing = analyzer.detect_missing_residues(source)

        counts = {}
        for entry in missing:
            counts[entry['chain']] = counts.get(entry['chain'], 0) + 1
        self.assertEqual(counts, {'A': 3, 'B': 5})

    def test_feature28_degrades_quietly(self):
        gapped = [l for l in _protein_lines() if not (20 <= int(l[22:26]) <= 24)]
        # No SEQRES: nothing to compare against.
        self.assertEqual(
            analyzer.detect_missing_residues(_write(self.path('n.pdb'), gapped)), [])
        # Complete structure.
        self.assertEqual(
            analyzer.detect_missing_residues(
                _write(self.path('f.pdb'), self._seqres() + _protein_lines())), [])
        # Unreadable path must not raise.
        self.assertEqual(
            analyzer.detect_missing_residues(self.path('nope.pdb')), [])


class TestScienceCorrectness(TempDirTest):
    """
    Physics and chemistry checks, as distinct from feature checks.

    These ask whether the result is defensible, not whether the code ran.
    """

    def test_disulfide_cysteines_are_not_protonated(self):
        """
        1CRN has three disulfide bridges. A bridged cysteine is oxidised;
        an HG on its SG invents a fourth bond on sulphur and breaks the bridge.
        """
        import numpy as np

        structure = io.load_pdb(CRN)
        io.save_pdb(structure, self.path('clean.pdb'),
                    select=cleaner.clean_structure(structure), source_pdb=CRN)
        protonator.add_hydrogens(self.path('clean.pdb'), self.path('out.pdb'))

        rows = [l for l in pathlib.Path(self.path('out.pdb')).read_text().splitlines()
                if l[:6] in ('ATOM  ', 'HETATM')]
        sg = {int(l[22:26]): np.array([float(l[30:38]), float(l[38:46]),
                                       float(l[46:54])])
              for l in rows if l[12:16].strip() == 'SG'}

        bridged = set()
        keys = sorted(sg)
        for i, a in enumerate(keys):
            for b in keys[i + 1:]:
                if np.linalg.norm(sg[a] - sg[b]) < 2.5:
                    bridged.update({a, b})
        self.assertEqual(len(bridged), 6, 'expected three disulfide bridges')

        with_hg = {int(l[22:26]) for l in rows
                   if l[12:16].strip() == 'HG' and l[17:20].strip() == 'CYS'}
        self.assertFalse(with_hg & bridged,
                         'a bridged cysteine was given a thiol hydrogen')

    def test_terminal_oxt_is_added_when_missing(self):
        """
        A chain deposited without a capped terminus has no OXT, matches no
        force-field template, and fails every minimisation tier.
        """
        lines = [l for l in _protein_lines() if l[12:16].strip() != 'OXT']
        source = _write(self.path('noterm.pdb'), lines)

        result = protonator.add_hydrogens(source, self.path('out.pdb'))
        written = pathlib.Path(self.path('out.pdb')).read_text()

        self.assertIn('OXT', written, 'terminal oxygen was not added')
        self.assertGreater(result['terminals_repaired'], 0)

    def test_structural_water_requires_a_polar_contact(self):
        """
        Proximity to a carbon does not hold a water in place. Every retained
        water must have an N, O or S within the hydrogen-bond cutoff.
        """
        import numpy as np

        lines = _protein_lines()
        x, y, z = _first_ca_xyz(lines)

        # one water hydrogen bonded to a backbone oxygen, one packed against
        # a carbon with no polar partner, one far away
        polar_anchor = next(l for l in lines if l[12:16].strip() == 'O')
        px, py, pz = (float(polar_anchor[30:38]), float(polar_anchor[38:46]),
                      float(polar_anchor[46:54]))

        rows = lines + [
            f'HETATM  900  O   HOH A 700    {px + 2.8:8.3f}{py:8.3f}{pz:8.3f}'
            '  1.00 10.00           O  \n',
            f'HETATM  901  O   HOH A 701    {x + 500:8.3f}{y:8.3f}{z:8.3f}'
            '  1.00 10.00           O  \n',
        ]
        source = _write(self.path('w.pdb'), rows)

        structure = io.load_pdb(source)
        select = cleaner.clean_structure(structure, remove_water=True,
                                         keep_structural_waters=True)
        io.save_pdb(structure, self.path('out.pdb'), select=select,
                    source_pdb=source)

        kept = [l for l in pathlib.Path(self.path('out.pdb')).read_text().splitlines()
                if l.startswith('HETATM') and l[17:20].strip() == 'HOH']
        self.assertEqual(len(kept), 1, 'expected only the hydrogen-bonded water')
        self.assertEqual(int(kept[0][22:26]), 700)

        # and every kept water genuinely has a polar partner
        protein = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                            for l in lines
                            if l[76:78].strip() in ('N', 'O', 'S')])
        for row in kept:
            w = np.array([float(row[30:38]), float(row[38:46]), float(row[46:54])])
            self.assertLessEqual(np.linalg.norm(protein - w, axis=1).min(),
                                 cleaner.STRUCTURAL_WATER_CUTOFF + 0.01)

    def test_pocket_detection_scans_seven_axes(self):
        """
        Six axis-aligned rays make the result depend on input orientation.
        LIGSITE scans the three Cartesian axes plus four cubic diagonals.
        """
        import numpy as np

        directions = site_analyzer.ENCLOSURE_DIRECTIONS
        self.assertEqual(len(directions), 14)
        norms = np.linalg.norm(directions, axis=1)
        self.assertTrue(np.allclose(norms, 1.0), 'directions must be unit vectors')

        # every direction must have its opposite present
        for d in directions:
            self.assertTrue(any(np.allclose(d, -other) for other in directions))

        # and the diagonals must actually be there
        diagonal = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)
        self.assertTrue(any(np.allclose(diagonal, d) for d in directions))

    def test_drugability_factors_are_exposed(self):
        """The ranking is an unvalidated heuristic; it must be auditable."""
        source = self._protonated_crn()
        sites = site_analyzer.BindingSiteAnalyzer(source).analyze()
        self.assertTrue(sites)

        factors = sites[0]['drugability_factors']
        self.assertEqual(
            set(factors), {'volume', 'property_diversity', 'hydrophobic_balance',
                           'concavity', 'pharmacophore_density'})
        self.assertAlmostEqual(sum(f['weight'] for f in factors.values()), 1.0,
                               places=6)

        recomputed = sum(f['score'] * f['weight'] for f in factors.values())
        self.assertAlmostEqual(recomputed, sites[0]['drugability_score'],
                               delta=0.02)

    def test_convergence_is_a_force_test_not_an_energy_test(self):
        """
        A run stopped by the iteration cap still lowers the energy, so an
        energy-based flag calls it converged when it is nowhere near a minimum.
        """
        source = self._protonated_crn()
        original = minimizer.MAX_ITERATIONS
        try:
            minimizer.MAX_ITERATIONS = 5
            short = minimizer.minimize_structure(source, self.path('short.pdb'))
            minimizer.MAX_ITERATIONS = 1000
            full = minimizer.minimize_structure(source, self.path('full.pdb'))
        finally:
            minimizer.MAX_ITERATIONS = original

        # both lower the energy
        self.assertTrue(short['energy_decreased'])
        self.assertTrue(full['energy_decreased'])

        # but only the complete run is converged
        self.assertFalse(short['converged'],
                         'a 5-iteration run must not report convergence')
        self.assertTrue(full['converged'])

        self.assertGreater(short['rms_force_kJ_mol_nm'],
                           full['rms_force_kJ_mol_nm'])
        self.assertLessEqual(full['rms_force_kJ_mol_nm'],
                             minimizer.ENERGY_TOLERANCE)

    def test_minimisation_is_unconstrained(self):
        """
        Constraints exist for an MD timestep; nothing here runs dynamics, and
        they make the residual force incomparable to the tolerance.
        """
        self.assertIsNone(minimizer.MINIMISATION_CONSTRAINTS)

    def test_pdbqt_receptor_has_only_polar_hydrogens(self):
        """
        OpenBabel's -xh means 'preserve hydrogens', not 'merge non-polar
        hydrogens'. With it every hydrogen survived and was typed HD, the
        AutoDock type for a hydrogen on N or O, so carbon-bound hydrogens were
        presented to the scoring function as hydrogen-bond donors.
        """
        if not shutil.which('obabel'):
            self.skipTest('OpenBabel not on PATH')

        source = self._protonated_crn()
        total_h = sum(1 for l in pathlib.Path(source).read_text().splitlines()
                      if l[:6] in ('ATOM  ', 'HETATM') and l[76:78].strip() == 'H')

        ok, result = exporter.export_structure(source, self.path('r.pdb'), 'vina')
        self.assertTrue(ok, result)

        types = {}
        for line in pathlib.Path(result).read_text().splitlines():
            if line[:6] in ('ATOM  ', 'HETATM') and len(line) > 77:
                key = line[77:79].strip()
                types[key] = types.get(key, 0) + 1

        self.assertIn('HD', types)
        self.assertLess(types['HD'], total_h,
                        'every hydrogen was typed HD; non-polar ones were kept')
        self.assertNotIn('ROOT', pathlib.Path(result).read_text())

    def test_waters_survive_protonation_and_get_hydrogens(self):
        """
        Waters are handed to PDBFixer with the polymer so they are protonated.
        A stray END record ahead of them once truncated the file and every
        water was silently dropped.
        """
        lines = _protein_lines()
        anchor = next(l for l in lines if l[12:16].strip() == 'O')
        x, y, z = (float(anchor[30:38]), float(anchor[38:46]), float(anchor[46:54]))

        waters = [
            f'HETATM  9{i:02d}  O   HOH A {700 + i:3d}    '
            f'{x + 2.8 + i * 0.2:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00           O  \n'
            for i in range(3)
        ]
        source = _write(self.path('w.pdb'), lines + waters)

        protonator.add_hydrogens(source, self.path('out.pdb'))
        written = [l for l in pathlib.Path(self.path('out.pdb')).read_text().splitlines()
                   if l[:6] in ('ATOM  ', 'HETATM')]

        residues = {(l[21], l[22:26]) for l in written
                    if l[17:20].strip() == 'HOH'}
        self.assertEqual(len(residues), 3, 'waters were lost during protonation')

        hydrogens = sum(1 for l in written
                        if l[17:20].strip() == 'HOH' and l[76:78].strip() == 'H')
        self.assertEqual(hydrogens, 6, 'expected two hydrogens per water')

    def test_end_record_does_not_truncate_the_input(self):
        """`_write` appends END; it must not cut the file short."""
        from bioprep.core.protonator import _split_records

        lines = _protein_lines()
        water = ('HETATM  900  O   HOH A 700      10.000  10.000  10.000'
                 '  1.00 10.00           O  \n')
        source = _write(self.path('e.pdb'), lines + [water])

        polymer, waters, ligands, _ = _split_records(source)
        self.assertTrue(waters, 'water was not separated out')
        for line in polymer:
            self.assertNotEqual(line[:6].strip(), 'END',
                                'an END record survived into the polymer')

    def test_terminal_is_repaired_even_with_waters_present(self):
        """
        Waters carry the protein's chain id. Without a TER between them,
        PDBFixer stops seeing the last amino acid as the chain's end and never
        adds the terminal OXT, leaving a structure that will not parameterise.
        """
        lines = [l for l in _protein_lines() if l[12:16].strip() != 'OXT']
        water = ('HETATM  900  O   HOH A 700      10.000  10.000  10.000'
                 '  1.00 10.00           O  \n')
        source = _write(self.path('t.pdb'), lines + [water])

        result = protonator.add_hydrogens(source, self.path('out.pdb'))
        text = pathlib.Path(self.path('out.pdb')).read_text()

        self.assertGreater(result['terminals_repaired'], 0,
                           'no terminal added while waters were present')
        self.assertIn('OXT', text)
        self.assertIn('HOH', text)

    def _protonated_crn(self):
        structure = io.load_pdb(CRN)
        io.save_pdb(structure, self.path('c.pdb'),
                    select=cleaner.clean_structure(structure), source_pdb=CRN)
        protonator.add_hydrogens(self.path('c.pdb'), self.path('p.pdb'))
        return self.path('p.pdb')


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

    def test_feature50_large_box_coarsens_the_grid(self):
        """Memory guard: a 1.5 A grid over a virus-sized box is unaffordable."""
        import numpy as np

        class FakeAtom:
            def __init__(self, coord):
                self._coord = np.array(coord, dtype=float)

            def get_coord(self):
                return self._coord

        analyzer_obj = site_analyzer.BindingSiteAnalyzer.__new__(
            site_analyzer.BindingSiteAnalyzer)
        span = np.linspace(0, 200, 10)
        analyzer_obj.coords = np.array(
            [[x, y, z] for x in span for y in span for z in span])
        analyzer_obj.atoms = [FakeAtom(c) for c in analyzer_obj.coords]
        analyzer_obj.RESIDUE_PROPS = site_analyzer.RESIDUE_PROPS

        _, resolution = analyzer_obj._detect_pockets(grid_res=1.5)
        self.assertGreater(resolution, 1.5)

        box = float(np.prod(analyzer_obj.coords.max(axis=0)
                            - analyzer_obj.coords.min(axis=0) + 10))
        self.assertLessEqual(box / resolution ** 3, 250000 * 1.05)

    def test_feature57_pockets_sorted_and_renumbered(self):
        source = self._protonated()
        sites = site_analyzer.BindingSiteAnalyzer(source).analyze()

        self.assertLessEqual(len(sites), site_analyzer.MAX_POCKETS)
        scores = [s['drugability_score'] for s in sites]
        self.assertEqual(scores, sorted(scores, reverse=True))
        self.assertEqual([s['id'] for s in sites], list(range(1, len(sites) + 1)))
        self.assertTrue(all(s['volume'] >= site_analyzer.MIN_POCKET_VOLUME
                            for s in sites))

    def test_feature58_residue_property_table(self):
        props = site_analyzer.RESIDUE_PROPS
        self.assertEqual(len(props), 20)
        self.assertEqual(set(props.values()),
                         {'HYDROPHOBIC', 'POLAR', 'CHARGED_NEG',
                          'CHARGED_POS', 'NEUTRAL'})
        self.assertEqual(props['ASP'], 'CHARGED_NEG')
        self.assertEqual(props['LYS'], 'CHARGED_POS')

    def test_feature69_summary_matches_the_sites(self):
        source = self._protonated()
        analyzer_obj = site_analyzer.BindingSiteAnalyzer(source)
        sites = analyzer_obj.analyze()
        summary = analyzer_obj.get_summary(sites)

        self.assertEqual(summary['site_count'], len(sites))
        self.assertAlmostEqual(summary['total_volume'],
                               sum(s['volume'] for s in sites), places=1)
        if sites:
            self.assertEqual(summary['primary_volume'], sites[0]['volume'])

        empty = analyzer_obj.get_summary([])
        self.assertEqual(empty['site_count'], 0)
        self.assertIn('No significant', empty['text'])

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
