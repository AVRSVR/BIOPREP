import os
import json
import time
import queue
import shutil
import base64
import zipfile
import tempfile
import threading
import uuid
from datetime import datetime, timezone

from flask import Flask, render_template, request, send_file, jsonify, Response, stream_with_context
from werkzeug.utils import secure_filename

# Import BioPrep core logic
from bioprep.core.io import load_pdb, save_pdb
from bioprep.core.cleaner import clean_structure
from bioprep.core.protonator import add_hydrogens
from bioprep.core.analyzer import analyze_structure, detect_missing_residues
from bioprep.core.reporter import build_report, report_to_text, count_atoms_in_pdb
from bioprep.core.minimizer import minimize_structure
from bioprep.core.site_analyzer import BindingSiteAnalyzer
from bioprep.core.exporter import export_structure

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = tempfile.gettempdir()
app.config['MAX_CONTENT_LENGTH'] = 200 * 1024 * 1024  # 200 MB max

ALLOWED_EXTENSIONS = {'pdb', 'zip'}
TEMPLATES_FILE = os.path.join(os.path.dirname(__file__), 'templates_store.json')
JOBS_FILE = os.path.join(os.path.dirname(__file__), 'jobs_history.json')

# Persistent results store
RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'processed_data')
os.makedirs(RESULTS_DIR, exist_ok=True)

# Global state
batch_progress = {}
session_pdb_paths = {} # {session_id: {'raw': path, 'current': path}}



def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def load_templates_store():
    if os.path.exists(TEMPLATES_FILE):
        with open(TEMPLATES_FILE, 'r') as f:
            return json.load(f)
    return {}


def save_templates_store(data):
    with open(TEMPLATES_FILE, 'w') as f:
        json.dump(data, f, indent=2)


def save_to_history(job_data):
    """Save an entry to the persistent job history."""
    history = []
    if os.path.exists(JOBS_FILE):
        with open(JOBS_FILE, 'r') as f:
            try:
                history = json.load(f)
            except:
                history = []
    
    history.insert(0, job_data)
    history = history[:50] # Limit to 50 entries
    
    with open(JOBS_FILE, 'w') as f:
        json.dump(history, f, indent=2)


def get_history():
    """Retrieve all history entries."""
    if os.path.exists(JOBS_FILE):
        with open(JOBS_FILE, 'r') as f:
            try:
                return json.load(f)
            except:
                return []
    return []


# ─────────────────────────────────────────────────────────────────────────────
# ROUTES: Pages
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/')
def index():
    return render_template('index.html')


# ─────────────────────────────────────────────────────────────────────────────
# ROUTES: Analysis
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/api/analyze', methods=['POST'])
def analyze():
    """Reads the uploaded PDB and returns structural metadata."""
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400

    file = request.files['file']
    if file.filename == '' or not file.filename.lower().endswith('.pdb'):
        return jsonify({'error': 'Invalid file'}), 400

    try:
        filename = secure_filename(file.filename)
        input_path = os.path.join(app.config['UPLOAD_FOLDER'], f"raw_anlz_{filename}")
        file.save(input_path)

        structure = load_pdb(input_path)
        metadata = analyze_structure(structure)
        metadata['missing_residues'] = detect_missing_residues(input_path)

        session_id = str(uuid.uuid4())
        session_pdb_paths[session_id] = {'raw': input_path, 'current': input_path}

        return jsonify({
            'success': True, 
            'filename': filename, 
            'metadata': metadata,
            'session_id': session_id
        })
    except Exception as e:
        return jsonify({'error': f"Analysis failed: {str(e)}"}), 500


# ─────────────────────────────────────────────────────────────────────────────
# ROUTES: Precision Mode (single PDB)
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/api/process', methods=['POST'])
def process():
    """Precision Mode – single PDB, full control, returns JSON with PDB + report."""
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400

    file = request.files['file']
    if file.filename == '' or not file.filename.lower().endswith('.pdb'):
        return jsonify({'error': 'Invalid file'}), 400

    try:
        t_start = time.time()

        # Generate session_id early for persistent storage
        session_id = str(uuid.uuid4())

        # ── Parse parameters ──────────────────────────────────────────────────
        chains_input = request.form.get('chains', '[]')
        chains_to_keep = json.loads(chains_input) if chains_input != '[]' else None

        remove_water = request.form.get('remove_water', 'true').lower() == 'true'

        hets_input = request.form.get('remove_heteros', '[]')
        hets_to_remove = json.loads(hets_input)

        protect_ligands_input = request.form.get('protect_ligands', '[]')
        protect_ligands = json.loads(protect_ligands_input)

        ph_input = request.form.get('ph', '7.4')
        try:
            target_ph = float(ph_input)
        except ValueError:
            target_ph = 7.4

        keep_structural_waters = request.form.get('keep_structural_waters', 'false').lower() == 'true'
        reconstruct_loops = request.form.get('reconstruct_loops', 'false').lower() == 'true'
        add_missing_atoms_flag = request.form.get('add_missing_atoms', 'false').lower() == 'true'
        run_minimization = request.form.get('run_minimization', 'false').lower() == 'true'
        use_gbsa = request.form.get('use_gbsa', 'true').lower() == 'true'
        force_field = request.form.get('force_field', 'amber14')
        docking_target = request.form.get('docking_target', '') or None

        # ── Save uploaded file ────────────────────────────────────────────────
        filename = secure_filename(file.filename)
        input_path = os.path.join(app.config['UPLOAD_FOLDER'], f"raw_{filename}")
        file.save(input_path)

        base_name = os.path.splitext(filename)[0]
        temp_cleaned_path = os.path.join(app.config['UPLOAD_FOLDER'], f"temp_{base_name}.pdb")
        protonated_path = os.path.join(app.config['UPLOAD_FOLDER'], f"prot_{base_name}.pdb")
        final_output_path = os.path.join(app.config['UPLOAD_FOLDER'], f"{base_name}_clean.pdb")

        # ── Pre-processing analysis ───────────────────────────────────────────
        structure = load_pdb(input_path)
        pre_meta = analyze_structure(structure)
        atoms_before = pre_meta['atoms_total']
        all_chains = pre_meta['chains']
        all_hets = pre_meta['heteroatoms']
        water_count_before = pre_meta['water_count']
        missing_residues = detect_missing_residues(input_path)

        # Determine what will actually be removed for report
        hets_actually_removed = [h for h in hets_to_remove if h not in protect_ligands]
        hets_retained = [h for h in all_hets if h not in hets_actually_removed]

        # ── Step 1: Clean ─────────────────────────────────────────────────────
        select_obj = clean_structure(
            structure,
            target_chains=chains_to_keep,
            remove_water=remove_water,
            remove_heteroatoms=hets_actually_removed,
            keep_structural_waters=keep_structural_waters,
        )
        save_pdb(structure, temp_cleaned_path, select=select_obj)

        # ── Step 2: Protonate ─────────────────────────────────────────────────
        add_hydrogens(
            temp_cleaned_path,
            protonated_path,
            ph=target_ph,
            reconstruct_loops=reconstruct_loops,
            add_missing_atoms=add_missing_atoms_flag,
        )

        # ── Step 3: Optional energy minimization ──────────────────────────────
        minimization_stats = None
        if run_minimization:
            minimization_stats = minimize_structure(protonated_path, final_output_path, force_field=force_field, use_gbsa=use_gbsa)
        else:
            shutil.copy2(protonated_path, final_output_path)

        # ── Count atoms after ─────────────────────────────────────────────────
        atoms_after = count_atoms_in_pdb(final_output_path)

        # ── Recount missing residues after repair ─────────────────────────────
        if reconstruct_loops or add_missing_atoms_flag:
            missing_residues_after = detect_missing_residues(final_output_path)
        else:
            missing_residues_after = missing_residues

        # ── Build report ──────────────────────────────────────────────────────
        report = build_report(
            filename=filename,
            chains_detected=all_chains,
            chains_retained=chains_to_keep if chains_to_keep else all_chains,
            waters_removed=water_count_before if remove_water else 0,
            heteroatoms_removed=hets_actually_removed,
            heteroatoms_retained=hets_retained,
            hydrogens_added=True,
            ph_used=target_ph,
            missing_residues=missing_residues_after,
            atoms_before=atoms_before,
            atoms_after=atoms_after,
            minimization_stats=minimization_stats,
            docking_target=docking_target,
            processing_time_s=time.time() - t_start,
        )

        # ── Step 4: Export Formatting ─────────────────────────────────────────
        viewer_path = final_output_path
        if docking_target:
            success, export_res = export_structure(final_output_path, final_output_path, docking_target)
            if success and str(export_res).endswith('.pdbqt'):
                final_output_path = str(export_res)
                # Filename in response should reflect the change
                filename = os.path.basename(final_output_path)

        # ── Persist structure for history viewing after restart ─────────────
        persistent_path = os.path.join(RESULTS_DIR, f"{session_id}.pdb")
        shutil.copy2(viewer_path, persistent_path) # Use viewer_path for persistent storage

        # ── Encode results as base64 for JSON response ───────────────────
        with open(viewer_path, 'rb') as f_view:
            viewer_b64 = base64.b64encode(f_view.read()).decode('utf-8')
            
        with open(final_output_path, 'rb') as f_out:
            download_b64 = base64.b64encode(f_out.read()).decode('utf-8')

        session_pdb_paths[session_id] = {
            'raw': input_path, 
            'current': viewer_path,  # Use viewer path for history 3D rendering
            'download': final_output_path
        }

        # ── Record job to history ─────────────────────────────────────────────
        job_info = {
            'id': session_id,
            'filename': filename,
            'timestamp': datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M'),
            'energy_before': minimization_stats['energy_before_kJ_mol'] if minimization_stats else 'N/A',
            'energy_after': minimization_stats['energy_after_kJ_mol'] if minimization_stats else 'N/A',
            'status': 'success',
            'format': (docking_target.upper() if docking_target else 'PDB'),
            'report': report
        }
        save_to_history(job_info)

        return jsonify({
            'success': True,
            'filename': filename,
            'report': report,
            'report_text': report_to_text(report),
            'viewer_pdb_b64': viewer_b64,
            'pdb_b64': download_b64,
            'session_id': session_id
        })

    except Exception as e:
        return jsonify({'error': f"Processing failed: {str(e)}"}), 500


@app.route('/api/history', methods=['GET'])
def list_history():
    """Return the list of past jobs."""
    return jsonify(get_history())


@app.route('/api/history/pdb/<job_id>', methods=['GET'])
def get_history_pdb(job_id):
    """Retrieve the PDB content for a historical job if it still exists."""
    # Try memory first
    path_data = session_pdb_paths.get(job_id)
    pdb_path = path_data.get('current') if path_data else None

    # Fallback to persistent storage
    if not pdb_path or not os.path.exists(pdb_path):
        pdb_path = os.path.join(RESULTS_DIR, f"{job_id}.pdb")

    if not os.path.exists(pdb_path):
        return "Structure no longer available on this server.", 404

    return send_file(pdb_path, mimetype='text/plain')


# (Redundant template routes removed, consolidated at bottom of file)


# ─────────────────────────────────────────────────────────────────────────────
# ROUTES: Batch Mode (multiple PDB files)
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/api/batch', methods=['POST'])
def batch_process():
    """
    Batch Mode – accept multiple PDB files + settings JSON.
    Returns a ZIP with cleaned PDBs + processing_report.txt.
    """
    try:
        files = request.files.getlist('files')
        if not files:
            return jsonify({'error': 'No files uploaded'}), 400

        config_input = request.form.get('config', '{}')
        try:
            config = json.loads(config_input)
        except json.JSONDecodeError:
            return jsonify({'error': 'Invalid config JSON'}), 400

        # Parse config
        chains_to_keep = config.get('chains', None)
        remove_water = config.get('remove_water', True)
        hets_to_remove = config.get('remove_heteros', [])
        protect_ligands = config.get('protect_ligands', [])
        target_ph = float(config.get('ph', 7.4))
        keep_structural_waters = config.get('keep_structural_waters', False)
        reconstruct_loops = config.get('reconstruct_loops', False)
        add_missing_atoms_flag = config.get('add_missing_atoms', False)
        run_minimization = config.get('run_minimization', False)
        use_gbsa = config.get('use_gbsa', True)
        force_field = config.get('force_field', 'amber14')
        docking_target = config.get('docking_target') or None

        batch_dir = tempfile.mkdtemp(prefix="bioprep_batch_")
        output_dir = os.path.join(batch_dir, 'output')
        os.makedirs(output_dir, exist_ok=True)

        results = []
        total_waters_removed = 0
        t_batch_start = time.time()

        hets_actually_removed = [h for h in hets_to_remove if h not in protect_ligands]

        for file in files:
            if not file.filename.lower().endswith('.pdb'):
                results.append({'file': file.filename, 'status': 'skipped', 'reason': 'not a PDB'})
                continue

            fname = secure_filename(file.filename)
            base_name = os.path.splitext(fname)[0]
            in_path = os.path.join(batch_dir, fname)
            temp_path = os.path.join(batch_dir, f"temp_{base_name}.pdb")
            prot_path = os.path.join(batch_dir, f"prot_{base_name}.pdb")
            out_path = os.path.join(output_dir, f"{base_name}_clean.pdb")

            try:
                file.save(in_path)
                structure = load_pdb(in_path)
                pre_meta = analyze_structure(structure)
                atoms_before = pre_meta['atoms_total']
                water_count_before = pre_meta['water_count']
                missing_residues = detect_missing_residues(in_path)

                select_obj = clean_structure(
                    structure,
                    target_chains=chains_to_keep,
                    remove_water=remove_water,
                    remove_heteroatoms=hets_actually_removed,
                    keep_structural_waters=keep_structural_waters,
                )
                save_pdb(structure, temp_path, select=select_obj)
                add_hydrogens(temp_path, prot_path, ph=target_ph,
                             reconstruct_loops=reconstruct_loops,
                             add_missing_atoms=add_missing_atoms_flag)

                mini_stats = None
                if run_minimization:
                    mini_stats = minimize_structure(prot_path, out_path, force_field=force_field, use_gbsa=use_gbsa)
                else:
                    shutil.copy2(prot_path, out_path)

                atoms_after = count_atoms_in_pdb(out_path)
                total_waters_removed += int(water_count_before) if remove_water else 0

                chains_retained = chains_to_keep if chains_to_keep else pre_meta['chains']
                all_hets = pre_meta['heteroatoms']
                hets_retained = [h for h in all_hets if h not in hets_actually_removed]

                report = build_report(
                    filename=fname,
                    chains_detected=pre_meta['chains'],
                    chains_retained=chains_retained,
                    waters_removed=water_count_before if remove_water else 0,
                    heteroatoms_removed=hets_actually_removed,
                    heteroatoms_retained=hets_retained,
                    hydrogens_added=True,
                    ph_used=target_ph,
                    missing_residues=missing_residues,
                    atoms_before=atoms_before,
                    atoms_after=atoms_after,
                    minimization_stats=mini_stats,
                    docking_target=docking_target,
                )
                results.append({'file': fname, 'status': 'success', 'report': report})

            except Exception as e:
                results.append({'file': fname, 'status': 'failed', 'reason': str(e)})

        # Build batch report text
        successful = [r for r in results if r.get('status') == 'success']
        failed = [r for r in results if r.get('status') == 'failed']
        avg_waters = (total_waters_removed / len(successful)) if successful else 0
        elapsed = round(time.time() - t_batch_start, 2)

        batch_report_lines = [
            "=" * 60,
            "         BIOPREP – BATCH PROCESSING REPORT",
            "=" * 60,
            f"  Generated       : {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
            f"  Total uploaded  : {len(files)}",
            f"  Successfully processed : {len(successful)}",
            f"  Failed          : {len(failed)}",
            f"  Avg waters removed : {avg_waters:.1f}",
            f"  Total time      : {elapsed} s",
            "",
            "── FAILED FILES ────────────────────────────────────────────",
        ]
        if failed:
            for r in failed:
                batch_report_lines.append(f"  {r['file']}: {r.get('reason', 'unknown error')}")
        else:
            batch_report_lines.append("  None")

        batch_report_lines += [
            "",
            "── INDIVIDUAL REPORTS ──────────────────────────────────────",
        ]
        for r in successful:
            if 'report' in r:
                batch_report_lines.append(report_to_text(r['report']))

        batch_report_text = "\n".join(batch_report_lines)

        # Write individual PDB reports alongside output files
        with open(os.path.join(output_dir, 'processing_report.txt'), 'w', encoding='utf-8') as f:
            f.write(batch_report_text)

        # Package everything into a ZIP
        output_zip_base = os.path.join(batch_dir, "bioprep_batch_results")
        shutil.make_archive(output_zip_base, 'zip', output_dir)
        output_zip_path = output_zip_base + ".zip"

        if not successful:
            return jsonify({'error': 'No files were successfully processed', 'details': results}), 400

        return send_file(
            output_zip_path,
            as_attachment=True,
            download_name="bioprep_batch_results.zip",
            mimetype="application/zip",
        )
    except Exception as e:
        return jsonify({'error': f"Batch processing crashed: {str(e)}"}), 500


# ─────────────────────────────────────────────────────────────────────────────
# ROUTES: High-Throughput Mode (ZIP in / ZIP out)
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/api/high-throughput', methods=['POST'])
def high_throughput():
    """
    High-Throughput Mode – accept a ZIP of PDBs, auto-pipeline, return ZIP + logs.
    Supports large batches (50–200+ structures).
    """
    if 'file' not in request.files:
        return jsonify({'error': 'No file uploaded'}), 400

    file = request.files['file']
    if not file.filename.lower().endswith('.zip'):
        return jsonify({'error': 'Please upload a .zip file containing PDB files'}), 400

    config_input = request.form.get('config', '{}')
    try:
        config = json.loads(config_input)
    except json.JSONDecodeError:
        config = {}

    remove_water = config.get('remove_water', True)
    hets_to_remove = config.get('remove_heteros', [])
    protect_ligands = config.get('protect_ligands', [])
    target_ph = float(config.get('ph', 7.4))
    run_minimization = config.get('run_minimization', False)
    force_field = config.get('force_field', 'amber14')
    reconstruct_loops = config.get('reconstruct_loops', False)
    add_missing_atoms_flag = config.get('add_missing_atoms', False)
    use_gbsa = config.get('use_gbsa', True)
    docking_target = config.get('docking_target') or None

    ht_dir = tempfile.mkdtemp(prefix="bioprep_ht_")
    extract_dir = os.path.join(ht_dir, 'input')
    output_dir = os.path.join(ht_dir, 'output')
    os.makedirs(extract_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    zip_path = os.path.join(ht_dir, secure_filename(file.filename))
    file.save(zip_path)

    try:
        with zipfile.ZipFile(zip_path, 'r') as zp:
            zp.extractall(extract_dir)
    except Exception as e:
        return jsonify({'error': f"Failed to extract ZIP: {str(e)}"}), 400

    logs = []
    processed = 0
    failed_list = []
    t_start = time.time()
    hets_actually_removed = [h for h in hets_to_remove if h not in protect_ligands]

    for root, _, files_list in os.walk(extract_dir):
        for fname in files_list:
            if not fname.lower().endswith('.pdb'):
                continue
            in_path = os.path.join(root, fname)
            base_name = os.path.splitext(fname)[0]
            
            # Prevent naming collisions in bulk output
            final_fname = f"{base_name}_clean.pdb"
            counter = 1
            while os.path.exists(os.path.join(output_dir, final_fname)):
                final_fname = f"{base_name}_{counter}_clean.pdb"
                counter += 1
                
            temp_path = os.path.join(ht_dir, f"temp_{base_name}_{counter}.pdb")
            prot_path = os.path.join(ht_dir, f"prot_{base_name}_{counter}.pdb")
            out_path = os.path.join(output_dir, final_fname)

            try:
                structure = load_pdb(in_path)
                pre_meta = analyze_structure(structure)
                select_obj = clean_structure(
                    structure,
                    remove_water=remove_water,
                    remove_heteroatoms=hets_actually_removed,
                )
                save_pdb(structure, temp_path, select=select_obj)
                add_hydrogens(temp_path, prot_path, ph=target_ph,
                             reconstruct_loops=reconstruct_loops,
                             add_missing_atoms=add_missing_atoms_flag)

                if run_minimization:
                    minimize_structure(prot_path, out_path, force_field=force_field, use_gbsa=use_gbsa)
                else:
                    shutil.copy2(prot_path, out_path)

                atoms_after = count_atoms_in_pdb(out_path)
                logs.append(f"[OK]  {fname}  →  {atoms_after} atoms")
                processed += 1
            except Exception as e:
                logs.append(f"[FAIL] {fname}: {str(e)}")
                failed_list.append(fname)

    elapsed = round(time.time() - t_start, 2)
    log_header = [
        "=" * 60,
        "    BIOPREP – HIGH-THROUGHPUT PROCESSING LOG",
        "=" * 60,
        f"  Generated  : {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        f"  Processed  : {processed}",
        f"  Failed     : {len(failed_list)}",
        f"  Total time : {elapsed} s",
        "=" * 60,
        "",
    ]
    full_log = "\n".join(log_header + logs)
    with open(os.path.join(output_dir, 'processing_logs.txt'), 'w', encoding='utf-8') as f:
        f.write(full_log)

    output_zip_base = os.path.join(ht_dir, "bioprep_ht_results")
    shutil.make_archive(output_zip_base, 'zip', output_dir)
    output_zip_path = output_zip_base + ".zip"

    if processed == 0:
        return jsonify({'error': 'No valid PDB files were processed', 'log': full_log}), 400

    # FIX: Cleanup batch_progress entry to prevent memory leak on long-running server.
    # high_throughput uses ht_dir as the identifier, not session_id.
    batch_progress.pop(ht_dir, None)

    return send_file(
        output_zip_path,
        as_attachment=True,
        download_name="bioprep_ht_results.zip",
        mimetype="application/zip",
    )


# ─────────────────────────────────────────────────────────────────────────────
# ROUTES: Templates
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/api/templates', methods=['GET'])
def get_templates():
    return jsonify(load_templates_store())


@app.route('/api/templates', methods=['POST'])
def save_template():
    data = request.get_json()
    if not data or 'name' not in data or 'settings' not in data:
        return jsonify({'error': 'Need name and settings fields'}), 400
    store = load_templates_store()
    store[data['name']] = {
        'settings': data['settings'],
        'created_at': datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC'),
    }
    save_templates_store(store)
    return jsonify({'success': True, 'templates': store})


@app.route('/api/templates/<name>', methods=['DELETE'])
def delete_template(name):
    store = load_templates_store()
    if name in store:
        del store[name]
        save_templates_store(store)
    return jsonify({'success': True, 'templates': store})


# ─────────────────────────────────────────────────────────────────────────────
# ROUTES: Report download
# ─────────────────────────────────────────────────────────────────────────────

@app.route('/api/download-report', methods=['POST'])
def download_report():
    """Accept a report JSON body and return a downloadable .txt file."""
    data = request.get_json()
    if not data:
        return jsonify({'error': 'No report data'}), 400
    from bioprep.core.reporter import report_to_text
    text = report_to_text(data)
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, encoding='utf-8')
    tmp.write(text)
    tmp.close()
    return send_file(tmp.name, as_attachment=True, download_name='preparation_report.txt', mimetype='text/plain')


@app.route('/api/analyze-site', methods=['POST'])
@app.route('/api/analyze-site/<session_id>', methods=['POST'])
def analyze_binding_sites(session_id=None):
    """Detect pockets and analyze binding sites."""
    try:
        pdb_path = None
        if 'file' in request.files:
            file = request.files['file']
            if file and file.filename.endswith('.pdb'):
                temp_dir = tempfile.mkdtemp()
                pdb_path = os.path.join(temp_dir, secure_filename(file.filename))
                file.save(pdb_path)
        elif session_id:
            # Check if session_id exists in session_pdb_paths
            session_data = session_pdb_paths.get(session_id)
            if session_data:
                pdb_path = session_data.get('current')
            
        if not pdb_path or not os.path.exists(pdb_path):
            return jsonify({
                'error': 'No structure found for analysis. Please upload a file or process a structure first.',
                'code': 'NO_STRUCTURE'
            }), 400

        analyzer = BindingSiteAnalyzer(pdb_path)
        sites = analyzer.analyze()
        summary = analyzer.get_summary(sites)
        
        return jsonify({
            'success': True,
            'sites': sites,
            'summary': summary,
            'message': f"Found {len(sites)} potential binding sites"
        })

    except Exception as e:
        app.logger.error(f"Site analysis failed: {str(e)}")
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5000)
