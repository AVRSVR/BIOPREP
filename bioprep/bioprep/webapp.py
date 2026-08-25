"""
BioPrep HTTP layer.

This module handles requests and nothing else: parsing, storage and response
shaping. All preparation logic lives in ``bioprep.core.pipeline`` so the three
processing modes cannot drift apart.
"""

import base64
import glob
import io
import json
import logging
import os
import re
import shutil
import tempfile
import threading
import urllib.error
import urllib.request
import uuid
import zipfile
from collections import OrderedDict
from datetime import datetime, timezone

from flask import Flask, render_template, request, send_file, jsonify
from werkzeug.utils import secure_filename

from bioprep.core.pipeline import PipelineSettings, prepare_structure
from bioprep.core.io import load_pdb, ensure_pdb
from bioprep.core.analyzer import analyze_structure, detect_missing_residues
from bioprep.core.reporter import report_to_text
from bioprep.core.site_analyzer import BindingSiteAnalyzer

logger = logging.getLogger(__name__)

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 200 * 1024 * 1024   # 200 MB
# This is a single-user local tool (nothing served templates to more than one
# client), so there's no reason to keep the compiled-template cache around;
# it only made template edits silently invisible until a restart.
app.config['TEMPLATES_AUTO_RELOAD'] = True

HERE = os.path.dirname(os.path.abspath(__file__))

# Saved templates, job history and stored structures are mutable state, so they
# must not live inside the package: site-packages is frequently read-only, and
# reinstalling would delete the history. BIOPREP_DATA_DIR chooses where they go;
# otherwise the working directory. The bioprep/app.py shim points this at the
# project directory so running from source keeps using the existing files.
DATA_DIR = os.path.abspath(
    os.environ.get('BIOPREP_DATA_DIR') or os.getcwd())
os.makedirs(DATA_DIR, exist_ok=True)

TEMPLATES_FILE = os.path.join(DATA_DIR, 'templates_store.json')
JOBS_FILE = os.path.join(DATA_DIR, 'jobs_history.json')
RESULTS_DIR = os.path.join(DATA_DIR, 'processed_data')
os.makedirs(RESULTS_DIR, exist_ok=True)

HISTORY_LIMIT = 50
SESSION_LIMIT = 200              # bounded so long-running servers do not leak
MAX_ZIP_MEMBERS = 500
MAX_ZIP_UNCOMPRESSED = 2 * 1024 ** 3   # 2 GB


# ─────────────────────────────────────────────────────────────────────────────
# Storage helpers
# ─────────────────────────────────────────────────────────────────────────────

class JsonStore:
    """
    A small JSON file guarded by a lock and written atomically.

    Read-modify-write without either of those corrupts the file when two
    requests land at once, which is how a 345 KB history file ends up
    unparseable.
    """

    def __init__(self, path, default):
        self._path = path
        self._default = default
        self._lock = threading.Lock()

    def read(self):
        with self._lock:
            return self._read_unlocked()

    def _read_unlocked(self):
        if not os.path.exists(self._path):
            return json.loads(json.dumps(self._default))
        try:
            with open(self._path, 'r', encoding='utf-8') as fh:
                return json.load(fh)
        except (OSError, json.JSONDecodeError):
            logger.warning("Could not read %s; starting fresh", self._path)
            return json.loads(json.dumps(self._default))

    def update(self, mutator):
        """Apply ``mutator`` to the stored value and persist the result."""
        with self._lock:
            data = self._read_unlocked()
            data = mutator(data)
            tmp_path = f'{self._path}.tmp'
            with open(tmp_path, 'w', encoding='utf-8') as fh:
                json.dump(data, fh, indent=2)
            os.replace(tmp_path, self._path)   # atomic on POSIX and Windows
            return data


templates_store = JsonStore(TEMPLATES_FILE, {})
jobs_store = JsonStore(JOBS_FILE, [])


class SessionStore:
    """Bounded map of session id -> stored structure path."""

    def __init__(self, limit):
        self._limit = limit
        self._lock = threading.Lock()
        self._items = OrderedDict()

    def put(self, session_id, path):
        with self._lock:
            self._items[session_id] = path
            self._items.move_to_end(session_id)
            while len(self._items) > self._limit:
                _, evicted = self._items.popitem(last=False)
                _quietly_remove(evicted)

    def get(self, session_id):
        with self._lock:
            path = self._items.get(session_id)
            if path:
                self._items.move_to_end(session_id)
            return path


sessions = SessionStore(SESSION_LIMIT)


def _quietly_remove(path):
    try:
        if path and os.path.isfile(path):
            os.remove(path)
    except OSError:
        pass


def _store_result(session_id, source_path):
    """Copy a prepared structure into persistent storage; return its path."""
    destination = os.path.join(RESULTS_DIR, f'{session_id}.pdb')
    shutil.copy2(source_path, destination)
    sessions.put(session_id, destination)
    return destination


def _record_job(session_id, filename, report, status='success'):
    minimization = report.get('energy_minimization') or {}
    entry = {
        'id': session_id,
        'filename': filename,
        'timestamp': datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M'),
        'status': status,
        'energy_before': minimization.get('energy_before_kJ_mol'),
        'energy_after': minimization.get('energy_after_kJ_mol'),
        'minimization_status': minimization.get('status'),
        'format': (report.get('docking_target') or 'PDB').upper(),
        'report': report,
    }

    def mutator(history):
        if not isinstance(history, list):
            history = []
        history.insert(0, entry)
        for stale in history[HISTORY_LIMIT:]:
            _quietly_remove(os.path.join(RESULTS_DIR, f"{stale.get('id')}.pdb"))
        return history[:HISTORY_LIMIT]

    jobs_store.update(mutator)


def _send_and_cleanup(path, workdir, download_name, mimetype):
    """
    Send a generated archive and remove its working directory immediately.

    The archive is read into memory first. Streaming straight from disk would
    keep the file open past the end of the view, so the directory could only be
    removed by a ``call_on_close`` hook — and that hook depends on the WSGI
    server closing the iterable, which is not guaranteed. Deferred cleanup that
    quietly does not happen is what left 22 orphaned directories behind before.
    """
    with open(path, 'rb') as fh:
        payload = io.BytesIO(fh.read())
    shutil.rmtree(workdir, ignore_errors=True)
    payload.seek(0)
    return send_file(payload, as_attachment=True,
                     download_name=download_name, mimetype=mimetype)


def _sweep_orphaned_workdirs():
    """
    Remove bioprep working directories left by a previous crash.

    Cleanup inside each request is the primary mechanism; this only catches
    directories orphaned by a hard kill.
    """
    import time as _time
    cutoff = _time.time() - 3600
    pattern = os.path.join(tempfile.gettempdir(), 'bioprep_*')
    for path in glob.glob(pattern):
        try:
            if os.path.isdir(path) and os.path.getmtime(path) < cutoff:
                shutil.rmtree(path, ignore_errors=True)
                logger.info("Removed orphaned working directory %s", path)
        except OSError:
            pass


def _fail(message, status=400, **extra):
    payload = {'error': message}
    payload.update(extra)
    return jsonify(payload), status


def _server_error(message, exc):
    """Log the detail, return a generic message. Tracebacks leak paths."""
    logger.exception(message)
    return jsonify({'error': f'{message}. See server logs for details.'}), 500


STRUCTURE_EXTENSIONS = ('.pdb', '.ent', '.cif', '.mmcif')


def _require_pdb_upload(field='file'):
    if field not in request.files:
        return None, _fail('No file was uploaded.')
    upload = request.files[field]
    if not upload.filename:
        return None, _fail('No file was selected.')
    if not upload.filename.lower().endswith(STRUCTURE_EXTENSIONS):
        return None, _fail(
            'Unsupported file type. Upload a PDB (.pdb, .ent) or '
            'mmCIF (.cif, .mmcif) structure.')
    return upload, None


PDB_ID_PATTERN = re.compile(r'^[0-9][A-Za-z0-9]{3}$')
RCSB_FETCH_TIMEOUT = 20


def _fetch_pdb_by_id(pdb_id, workdir):
    """
    Download a structure straight from RCSB by its 4-character id.

    Tries legacy PDB format first, then mmCIF - some entries (typically
    large assemblies or cryo-EM structures, like 30IE elsewhere in this
    project's own testing) are deposited without a legacy PDB rendering at
    all and only 404 on that format, not on mmCIF.

    Returns (local_path, filename). Raises ValueError with a message safe to
    show the user on a bad id, a 404, or a network failure - the caller
    already knows how to turn that into a clean JSON error.
    """
    pdb_id = pdb_id.strip()
    if not PDB_ID_PATTERN.match(pdb_id):
        raise ValueError(
            f"'{pdb_id}' doesn't look like a PDB id - it should be 4 "
            "characters, starting with a digit (e.g. 1HSG).")
    pdb_id = pdb_id.upper()

    errors = []
    for ext in ('pdb', 'cif'):
        url = f'https://files.rcsb.org/download/{pdb_id}.{ext}'
        dest = os.path.join(workdir, f'{pdb_id}.{ext}')
        try:
            with urllib.request.urlopen(url, timeout=RCSB_FETCH_TIMEOUT) as resp:
                with open(dest, 'wb') as fh:
                    shutil.copyfileobj(resp, fh)
            return dest, f'{pdb_id}.{ext}'
        except urllib.error.HTTPError as exc:
            errors.append(f'{ext.upper()}: HTTP {exc.code}')
        except urllib.error.URLError as exc:
            # Not worth trying the other format too - the network itself is
            # the problem, not the format.
            raise ValueError(f'Could not reach RCSB ({exc.reason}).')

    raise ValueError(
        f"RCSB has no entry '{pdb_id}' in a format this tool reads "
        f"({'; '.join(errors)}). Check the id, or upload the file directly.")


def _settings_from_form():
    """Precision mode posts a flat multipart form."""
    data = {}
    for key in ('remove_water', 'keep_structural_waters', 'reconstruct_loops',
                'add_missing_atoms', 'run_minimization', 'use_gbsa', 'use_propka',
                'ph', 'force_field', 'docking_target'):
        if key in request.form:
            data[key] = request.form.get(key)

    for key in ('chains', 'remove_heteros', 'protect_ligands'):
        raw = request.form.get(key, '[]')
        try:
            parsed = json.loads(raw)
        except json.JSONDecodeError:
            parsed = []
        data[key] = parsed

    return PipelineSettings.from_mapping(data)


def _settings_from_config():
    """Batch and high-throughput modes post a `config` JSON blob."""
    try:
        config = json.loads(request.form.get('config', '{}'))
    except json.JSONDecodeError:
        config = {}
    return PipelineSettings.from_mapping(config)


# ─────────────────────────────────────────────────────────────────────────────
# Pages
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/')
def index():
    return render_template('index.html')


# ─────────────────────────────────────────────────────────────────────────────
# Analysis
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/api/analyze', methods=['POST'])
def analyze():
    """Return structural metadata for a PDB, uploaded or fetched by id, without processing it."""
    pdb_id = request.form.get('pdb_id', '').strip()
    if not pdb_id:
        upload, error = _require_pdb_upload()
        if error:
            return error

    workdir = tempfile.mkdtemp(prefix='bioprep_analyze_')
    try:
        if pdb_id:
            # Raises ValueError on a bad id, a 404, or an unreachable RCSB -
            # caught below and reported the same way a bad upload would be.
            input_path, filename = _fetch_pdb_by_id(pdb_id, workdir)
        else:
            filename = secure_filename(upload.filename)
            input_path = os.path.join(workdir, filename)
            upload.save(input_path)

        metadata = analyze_structure(load_pdb(input_path))
        metadata['missing_residues'] = detect_missing_residues(input_path)

        # The stored copy is served straight to the 3Dmol viewer (and to any
        # /api/analyze-site/<session_id> call) as PDB - an mmCIF upload has to
        # be converted here too, not just at /api/process time, or the viewer
        # gets handed raw mmCIF text told to parse as PDB and silently shows
        # nothing.
        viewer_path, _ = ensure_pdb(input_path, workdir)
        session_id = str(uuid.uuid4())
        _store_result(session_id, viewer_path)

        return jsonify({
            'success': True,
            'filename': filename,
            'metadata': metadata,
            'session_id': session_id,
        })
    except ValueError as exc:
        # Unreadable/non-PDB upload, bad id, or a fetch failure: the user can
        # fix this, so say what is wrong rather than returning a generic 500.
        return _fail(str(exc))
    except Exception as exc:
        return _server_error('Analysis failed', exc)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ─────────────────────────────────────────────────────────────────────────────
# Precision mode
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/api/process', methods=['POST'])
def process():
    """Prepare a single structure and return the result inline."""
    upload, error = _require_pdb_upload()
    if error:
        return error

    settings = _settings_from_form()
    session_id = str(uuid.uuid4())
    workdir = tempfile.mkdtemp(prefix='bioprep_process_')
    try:
        filename = secure_filename(upload.filename)
        input_path = os.path.join(workdir, filename)
        upload.save(input_path)

        outcome = prepare_structure(input_path, workdir, settings,
                                    original_filename=filename)
        report = outcome['report']

        _store_result(session_id, outcome['viewer_path'])
        _record_job(session_id, filename, report)

        with open(outcome['viewer_path'], 'rb') as fh:
            viewer_b64 = base64.b64encode(fh.read()).decode('ascii')

        # Only encode a second copy when the download really is a different
        # file; otherwise the response carried the structure twice.
        if outcome['download_path'] == outcome['viewer_path']:
            download_b64 = viewer_b64
            download_name = filename
        else:
            with open(outcome['download_path'], 'rb') as fh:
                download_b64 = base64.b64encode(fh.read()).decode('ascii')
            download_name = os.path.basename(outcome['download_path'])

        return jsonify({
            'success': True,
            'filename': download_name,
            'report': report,
            'report_text': report_to_text(report),
            'warnings': outcome['warnings'],
            'viewer_pdb_b64': viewer_b64,
            'pdb_b64': download_b64,
            'session_id': session_id,
        })

    except ValueError as exc:
        # Expected, user-correctable problems.
        return _fail(str(exc))
    except Exception as exc:
        return _server_error('Processing failed', exc)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ─────────────────────────────────────────────────────────────────────────────
# Batch mode
# ─────────────────────────────────────────────────────────────────────────────

def _process_many(named_streams, settings, workdir):
    """
    Prepare several structures into ``workdir/output``.

    ``named_streams`` yields (filename, save_callable). Returns
    (output_dir, results).
    """
    output_dir = os.path.join(workdir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    results = []
    used_names = set()

    for filename, save in named_streams:
        safe_name = secure_filename(filename)
        if not safe_name.lower().endswith('.pdb'):
            results.append({'file': filename, 'status': 'skipped',
                            'reason': 'not a .pdb file'})
            continue

        item_dir = os.path.join(workdir, f'item_{len(results)}')
        os.makedirs(item_dir, exist_ok=True)
        try:
            input_path = os.path.join(item_dir, safe_name)
            save(input_path)

            outcome = prepare_structure(input_path, item_dir, settings,
                                        original_filename=safe_name)

            base = os.path.splitext(safe_name)[0]
            extension = os.path.splitext(outcome['download_path'])[1]
            out_name = f'{base}_prepared{extension}'
            suffix = 1
            while out_name.lower() in used_names:
                out_name = f'{base}_{suffix}_prepared{extension}'
                suffix += 1
            used_names.add(out_name.lower())

            shutil.copy2(outcome['download_path'],
                         os.path.join(output_dir, out_name))
            results.append({'file': safe_name, 'status': 'success',
                            'output': out_name, 'report': outcome['report']})
        except Exception as exc:
            logger.exception("Failed to prepare %s", safe_name)
            results.append({'file': safe_name, 'status': 'failed',
                            'reason': str(exc)})
        finally:
            shutil.rmtree(item_dir, ignore_errors=True)

    return output_dir, results


def _summary_text(title, results, settings, elapsed):
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'failed']
    skipped = [r for r in results if r['status'] == 'skipped']

    lines = [
        '=' * 60,
        f'         BIOPREP - {title}',
        '=' * 60,
        f"  Generated  : {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        f'  Total      : {len(results)}',
        f'  Succeeded  : {len(successful)}',
        f'  Failed     : {len(failed)}',
        f'  Skipped    : {len(skipped)}',
        f'  Total time : {round(elapsed, 2)} s',
        '',
        '-- SETTINGS USED --------------------------------------------',
    ]
    for key, value in settings.to_dict().items():
        lines.append(f'  {key:<24}: {value}')

    if failed or skipped:
        lines += ['', '-- NOT PROCESSED --------------------------------------------']
        for r in failed + skipped:
            lines.append(f"  {r['file']}: {r.get('reason', 'unknown')}")

    lines += ['', '-- INDIVIDUAL REPORTS ---------------------------------------']
    for r in successful:
        lines.append(report_to_text(r['report']))

    return '\n'.join(lines)


@app.route('/api/batch', methods=['POST'])
def batch():
    """Prepare several uploaded PDBs and return them as a ZIP."""
    import time
    uploads = request.files.getlist('files')
    if not uploads:
        return _fail('No files were uploaded.')

    settings = _settings_from_config()
    workdir = tempfile.mkdtemp(prefix='bioprep_batch_')
    started = time.time()
    try:
        streams = [(u.filename, u.save) for u in uploads]
        output_dir, results = _process_many(streams, settings, workdir)

        if not any(r['status'] == 'success' for r in results):
            shutil.rmtree(workdir, ignore_errors=True)
            return _fail('No files could be prepared.', details=results)

        with open(os.path.join(output_dir, 'processing_report.txt'),
                  'w', encoding='utf-8') as fh:
            fh.write(_summary_text('BATCH PROCESSING REPORT', results,
                                   settings, time.time() - started))

        archive = shutil.make_archive(
            os.path.join(workdir, 'bioprep_batch_results'), 'zip', output_dir)
        return _send_and_cleanup(archive, workdir,
                                 'bioprep_batch_results.zip', 'application/zip')
    except Exception as exc:
        shutil.rmtree(workdir, ignore_errors=True)
        return _server_error('Batch processing failed', exc)


# ─────────────────────────────────────────────────────────────────────────────
# High-throughput mode
# ─────────────────────────────────────────────────────────────────────────────

def _safe_extract(zip_path, destination):
    """Extract a ZIP with member-count and expanded-size limits."""
    with zipfile.ZipFile(zip_path) as archive:
        members = [m for m in archive.infolist() if not m.is_dir()]
        if len(members) > MAX_ZIP_MEMBERS:
            raise ValueError(
                f'Archive contains {len(members)} files; the limit is '
                f'{MAX_ZIP_MEMBERS}.')
        total = sum(m.file_size for m in members)
        if total > MAX_ZIP_UNCOMPRESSED:
            raise ValueError(
                f'Archive expands to {total // 1024 ** 2} MB; the limit is '
                f'{MAX_ZIP_UNCOMPRESSED // 1024 ** 2} MB.')
        archive.extractall(destination)


@app.route('/api/high-throughput', methods=['POST'])
def high_throughput():
    """Prepare every PDB inside an uploaded ZIP and return a ZIP of results."""
    import time
    if 'file' not in request.files:
        return _fail('No file was uploaded.')
    upload = request.files['file']
    if not upload.filename or not upload.filename.lower().endswith('.zip'):
        return _fail('Please upload a .zip archive containing PDB files.')

    settings = _settings_from_config()
    workdir = tempfile.mkdtemp(prefix='bioprep_ht_')
    started = time.time()
    try:
        zip_path = os.path.join(workdir, secure_filename(upload.filename))
        upload.save(zip_path)

        extract_dir = os.path.join(workdir, 'input')
        os.makedirs(extract_dir, exist_ok=True)
        try:
            _safe_extract(zip_path, extract_dir)
        except (zipfile.BadZipFile, ValueError) as exc:
            shutil.rmtree(workdir, ignore_errors=True)
            return _fail(f'Could not read archive: {exc}')

        discovered = []
        for root, _, names in os.walk(extract_dir):
            for name in names:
                if name.lower().endswith('.pdb'):
                    source = os.path.join(root, name)
                    discovered.append(
                        (name, (lambda src: lambda dst: shutil.copy2(src, dst))(source))
                    )

        if not discovered:
            shutil.rmtree(workdir, ignore_errors=True)
            return _fail('The archive contained no .pdb files.')

        output_dir, results = _process_many(discovered, settings, workdir)

        if not any(r['status'] == 'success' for r in results):
            shutil.rmtree(workdir, ignore_errors=True)
            return _fail('No structures could be prepared.', details=results)

        with open(os.path.join(output_dir, 'processing_logs.txt'),
                  'w', encoding='utf-8') as fh:
            fh.write(_summary_text('HIGH-THROUGHPUT PROCESSING LOG', results,
                                   settings, time.time() - started))

        archive = shutil.make_archive(
            os.path.join(workdir, 'bioprep_ht_results'), 'zip', output_dir)
        return _send_and_cleanup(archive, workdir,
                                 'bioprep_ht_results.zip', 'application/zip')
    except Exception as exc:
        shutil.rmtree(workdir, ignore_errors=True)
        return _server_error('High-throughput processing failed', exc)


# ─────────────────────────────────────────────────────────────────────────────
# History
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/api/history', methods=['GET'])
def list_history():
    return jsonify(jobs_store.read())


@app.route('/api/history/pdb/<uuid:job_id>', methods=['GET'])
def get_history_pdb(job_id):
    """
    Serve a stored structure. The ``uuid`` converter constrains the id to a
    real UUID, so it cannot be used to walk out of the results directory.
    """
    job_id = str(job_id)
    path = sessions.get(job_id) or os.path.join(RESULTS_DIR, f'{job_id}.pdb')
    if not os.path.isfile(path):
        return _fail('That structure is no longer stored on this server.', 404)
    return send_file(path, mimetype='chemical/x-pdb')


# ─────────────────────────────────────────────────────────────────────────────
# Templates
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/api/templates', methods=['GET'])
def get_templates():
    return jsonify(templates_store.read())


@app.route('/api/templates', methods=['POST'])
def save_template():
    data = request.get_json(silent=True) or {}
    name = str(data.get('name', '')).strip()
    if not name or 'settings' not in data:
        return _fail("Both 'name' and 'settings' are required.")

    entry = {
        'settings': data['settings'],
        'created_at': datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC'),
    }

    def mutator(store):
        store[name] = entry
        return store

    return jsonify({'success': True, 'templates': templates_store.update(mutator)})


@app.route('/api/templates/<name>', methods=['DELETE'])
def delete_template(name):
    def mutator(store):
        store.pop(name, None)
        return store

    return jsonify({'success': True, 'templates': templates_store.update(mutator)})


# ─────────────────────────────────────────────────────────────────────────────
# Report download
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/api/download-report', methods=['POST'])
def download_report():
    data = request.get_json(silent=True)
    if not isinstance(data, dict):
        return _fail('Expected a JSON report object.')

    # Small enough to build in memory, which avoids a temp file entirely.
    payload = io.BytesIO(report_to_text(data).encode('utf-8'))
    return send_file(payload, as_attachment=True,
                     download_name='preparation_report.txt',
                     mimetype='text/plain')


# ─────────────────────────────────────────────────────────────────────────────
# Binding site analysis
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/api/analyze-site', methods=['POST'])
@app.route('/api/analyze-site/<uuid:session_id>', methods=['POST'])
def analyze_binding_sites(session_id=None):
    workdir = None
    try:
        if 'file' in request.files and request.files['file'].filename:
            upload, error = _require_pdb_upload()
            if error:
                return error
            workdir = tempfile.mkdtemp(prefix='bioprep_site_')
            pdb_path = os.path.join(workdir, secure_filename(upload.filename))
            upload.save(pdb_path)
        elif session_id is not None:
            pdb_path = sessions.get(str(session_id)) or os.path.join(
                RESULTS_DIR, f'{session_id}.pdb')
        else:
            return _fail('Upload a structure or supply a session id.',
                         code='NO_STRUCTURE')

        if not pdb_path or not os.path.isfile(pdb_path):
            return _fail('No stored structure found for analysis. Process a '
                         'structure first, or upload one directly.',
                         code='NO_STRUCTURE')

        analyzer = BindingSiteAnalyzer(pdb_path)
        sites = analyzer.analyze()
        return jsonify({
            'success': True,
            'sites': sites,
            'summary': analyzer.get_summary(sites),
            'message': f'Found {len(sites)} candidate binding sites.',
        })
    except Exception as exc:
        return _server_error('Binding site analysis failed', exc)
    finally:
        if workdir:
            shutil.rmtree(workdir, ignore_errors=True)


def main(host='127.0.0.1', port=5000):
    """Entry point for the ``bioprep-web`` console script."""
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)s %(name)s: %(message)s')
    logger.info("Storing templates, history and results under %s", DATA_DIR)
    _sweep_orphaned_workdirs()
    # debug=False: the Werkzeug debugger executes arbitrary code from the browser.
    app.run(host=host, port=port, debug=False)


if __name__ == '__main__':
    main()
