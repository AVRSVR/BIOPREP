"""
Route-level regression tests for webapp.py.

test_backend.py covers bioprep.core directly; this file is for bugs that only
show up at the Flask route layer - request handling, session storage, what
gets served back to the browser - rather than inside the pipeline itself.
"""

import os
import tempfile
import unittest

import conftest  # noqa: F401  - sets up sys.path

# The data directory is read from this env var at import time, so it has to
# be set before bioprep.webapp is imported - otherwise the test suite writes
# into the real jobs_history.json / processed_data next to the app.
_DATA_DIR = tempfile.mkdtemp(prefix='bioprep_webapp_test_')
os.environ['BIOPREP_DATA_DIR'] = _DATA_DIR

from bioprep.webapp import app  # noqa: E402  - see above
from Bio.PDB import PDBParser, MMCIFIO  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CRN = os.path.join(HERE, '..', 'bioprep', '1crn.pdb')


def _crn_as_mmcif(path):
    """Write crambin as mmCIF - a small, real, valid fixture, not hand-rolled."""
    structure = PDBParser(QUIET=True).get_structure('crn', CRN)
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(path)
    return path


class TestAnalyzeRoute(unittest.TestCase):
    def setUp(self):
        if not os.path.exists(CRN):
            self.skipTest('1crn.pdb not available')
        self.client = app.test_client()

    def test_mmcif_upload_is_converted_before_being_stored_for_the_viewer(self):
        """An mmCIF upload's session must serve PDB, not raw mmCIF.

        /api/history/pdb/<session_id> is fed straight to the 3Dmol viewer as
        format "pdb". /api/process already converts mmCIF via ensure_pdb();
        /api/analyze did not, so uploading a .cif file stored (and served)
        the raw mmCIF text under a session id, and the viewer silently
        rendered nothing - it was told to parse mmCIF as PDB. Caught live
        against a real 3OIE.cif upload, not just this synthetic one.
        """
        with tempfile.TemporaryDirectory(prefix='bioprep_test_') as tmp:
            cif_path = _crn_as_mmcif(os.path.join(tmp, 'crn.cif'))

            with open(cif_path, 'rb') as fh:
                resp = self.client.post(
                    '/api/analyze',
                    data={'file': (fh, 'crn.cif')},
                    content_type='multipart/form-data',
                )
            self.assertEqual(resp.status_code, 200, resp.get_json())
            session_id = resp.get_json()['session_id']

            served = self.client.get(f'/api/history/pdb/{session_id}')
            self.assertEqual(served.status_code, 200)
            text = served.get_data(as_text=True)

            self.assertFalse(text.lstrip().startswith('data_'),
                             'session still serves raw mmCIF, not a PDB conversion')
            self.assertIn('ATOM', text)
            atom_lines = [l for l in text.splitlines() if l.startswith('ATOM')]
            self.assertEqual(len(atom_lines), 327,
                             'converted PDB does not carry the same atoms as the source')


if __name__ == '__main__':
    unittest.main()
