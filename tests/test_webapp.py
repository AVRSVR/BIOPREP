"""
Route-level regression tests for webapp.py.

test_backend.py covers bioprep.core directly; this file is for bugs that only
show up at the Flask route layer - request handling, session storage, what
gets served back to the browser - rather than inside the pipeline itself.
"""

import io
import os
import tempfile
import unittest
import unittest.mock as mock
import urllib.error

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


class TestPublicDemoGate(unittest.TestCase):
    """BIOPREP_PUBLIC_DEMO=1 turns off the two routes that leak across visitors.

    /api/history returns every job's filename and full report to anyone, and
    saved templates are global too - fine on one person's own machine, a
    real cross-visitor leak on a shared public deployment. get_history_pdb is
    deliberately NOT gated: the Prepare/Sites viewer depends on it to load
    back the structure from the session id its own /api/analyze or
    /api/process call just returned, and that id isn't discoverable without
    already having it.
    """

    def setUp(self):
        self.client = app.test_client()

    def test_history_and_templates_disabled_under_public_demo(self):
        import bioprep.webapp as webapp
        with mock.patch.object(webapp, 'PUBLIC_DEMO', True):
            self.assertEqual(self.client.get('/api/history').status_code, 404)
            self.assertEqual(self.client.get('/api/templates').status_code, 404)
            self.assertEqual(
                self.client.post('/api/templates', json={'name': 'x', 'settings': {}}).status_code, 404)
            self.assertEqual(self.client.delete('/api/templates/x').status_code, 404)

    def test_history_and_templates_open_by_default(self):
        self.assertEqual(self.client.get('/api/history').status_code, 200)
        self.assertEqual(self.client.get('/api/templates').status_code, 200)

    def test_history_pdb_by_id_stays_open_under_public_demo(self):
        """Gating the list must not gate the viewer's own fetch-by-id call."""
        import bioprep.webapp as webapp
        with mock.patch.object(webapp, 'PUBLIC_DEMO', True):
            with tempfile.TemporaryDirectory(prefix='bioprep_test_') as tmp:
                cif_path = _crn_as_mmcif(os.path.join(tmp, 'crn.cif'))
                with open(cif_path, 'rb') as fh:
                    resp = self.client.post('/api/analyze', data={'file': (fh, 'crn.cif')},
                                            content_type='multipart/form-data')
            session_id = resp.get_json()['session_id']
            self.assertEqual(
                self.client.get(f'/api/history/pdb/{session_id}').status_code, 200)

    def test_nav_and_footer_reflect_public_demo(self):
        import bioprep.webapp as webapp
        with mock.patch.object(webapp, 'PUBLIC_DEMO', True):
            html = self.client.get('/').get_data(as_text=True)
        self.assertNotIn('data-view="history"', html)
        self.assertNotIn('data-view="templates"', html)
        self.assertIn('Hosted demo', html)


class TestFetchByPdbId(unittest.TestCase):
    """/api/analyze with pdb_id instead of a file upload.

    Network calls are mocked throughout - real RCSB access is exercised
    manually (see BACKEND_AUDIT.md), not on every test run, so the suite
    stays fast and doesn't fail in an offline CI environment.
    """

    def setUp(self):
        self.client = app.test_client()

    def test_bad_id_format_is_rejected_without_a_network_call(self):
        with mock.patch('urllib.request.urlopen', side_effect=AssertionError(
                'should not have tried the network for an invalid id')) as m:
            resp = self.client.post('/api/analyze', data={'pdb_id': 'nope!'})
        m.assert_not_called()
        self.assertEqual(resp.status_code, 400)
        self.assertIn('id', resp.get_json()['error'])

    def test_no_legacy_pdb_falls_back_to_mmcif(self):
        """Some entries (large assemblies, cryo-EM) 404 on .pdb but not .cif."""
        cif_path = _crn_as_mmcif(os.path.join(tempfile.mkdtemp(), 'crn.cif'))
        with open(cif_path, 'rb') as fh:
            cif_bytes = fh.read()

        def fake_urlopen(url, timeout=None):
            if url.endswith('.pdb'):
                raise urllib.error.HTTPError(url, 404, 'Not Found', {}, None)
            return io.BytesIO(cif_bytes)

        with mock.patch('urllib.request.urlopen', side_effect=fake_urlopen):
            resp = self.client.post('/api/analyze', data={'pdb_id': '1abc'})

        self.assertEqual(resp.status_code, 200, resp.get_json())
        data = resp.get_json()
        self.assertEqual(data['filename'], '1ABC.cif')
        self.assertEqual(data['metadata']['atoms_total'], 327)

    def test_network_failure_does_not_try_the_second_format(self):
        attempted = []

        def fake_urlopen(url, timeout=None):
            attempted.append(url)
            raise urllib.error.URLError('name resolution failed')

        with mock.patch('urllib.request.urlopen', side_effect=fake_urlopen):
            resp = self.client.post('/api/analyze', data={'pdb_id': '1abc'})

        self.assertEqual(len(attempted), 1,
                         'tried a second format after a network-level failure, '
                         'not just a 404 - the format was never the problem')
        self.assertEqual(resp.status_code, 400)
        self.assertIn('RCSB', resp.get_json()['error'])


if __name__ == '__main__':
    unittest.main()
