/* ═══════════════════════════════════════════
   BioPrep Pro – main.js
   All 3 modes + templates + visual diff
   ═══════════════════════════════════════════ */
document.addEventListener('DOMContentLoaded', () => {

    // ─────────────────────────────────────── //
    // SHARED STATE
    // ─────────────────────────────────────── //
    let precisionFile = null;
    let precisionRawPdb = null;
    let precisionCleanPdb = null;
    let precisionReport = null;
    let glviewer = null;
    let diffMode = false;
    let templatesStore = {}; // Global store for user presets

    let batchFiles = [];
    let batchResultBlob = null;

    let htFile = null;
    let htResultBlob = null;
    let currentSessionId = null;
    let siteManualFile = null; // Stored file for manual analysis

    let discoveredSites = [];
    let siteShapes = []; // Array to track 3Dmol shapes for sites
    let historyStore = {}; // Memory-safe store for job data


    let modalGlViewer = null;
    let currentZipFiles = {};

    // ─────────────────────────────────────── //
    // TAB SWITCHING
    // ─────────────────────────────────────── //
    document.querySelectorAll('.tab-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
            document.querySelectorAll('.tab-content').forEach(p => {
                p.classList.remove('active');
                p.classList.add('hidden');
            });
            btn.classList.add('active');
            const panel = document.getElementById(`tab-panel-${btn.dataset.tab}`);
            panel.classList.add('active');
            panel.classList.remove('hidden');

            // Manage right stack panels
            updateRightStackVisibility(btn.dataset.tab);

            // FIX: Ensure viewer resizes after tab switch to avoid black/distorted canvas
            if (glviewer) {
                setTimeout(() => {
                    glviewer.resize();
                    glviewer.render();
                }, 50);
            }
        });
    });

    function updateRightStackVisibility(activeTab) {
        // Hide all result panels first
        document.getElementById('report-panel').classList.add('hidden');
        document.getElementById('batch-results').classList.add('hidden');
        document.getElementById('ht-results').classList.add('hidden');

        if (activeTab === 'precision') {
            if (precisionReport) document.getElementById('report-panel').classList.remove('hidden');
        } else if (activeTab === 'batch') {
            if (batchResultBlob) document.getElementById('batch-results').classList.remove('hidden');
        } else if (activeTab === 'highthroughput') {
            if (htResultBlob) document.getElementById('ht-results').classList.remove('hidden');
        } else if (activeTab === 'site') {
            if (discoveredSites.length > 0) document.getElementById('site-results').classList.remove('hidden');
        }
    }

    // ─────────────────────────────────────── //
    // VIEWER UTILS
    // ─────────────────────────────────────── //
    const viewerEl = document.getElementById('3dmol-viewer');
    const placeholder = document.getElementById('viewer-placeholder');

    function initViewer() {
        if (!glviewer) {
            viewerEl.innerHTML = '';
            // Ensure 3Dmol is available
            if (typeof $3Dmol === 'undefined') {
                console.error("3Dmol.js not loaded. Structural viewer will not be available.");
                showStatus('precision-status', '3D Viewer library failed to load. Please check your internet connection and refresh.', 'error');
                return;
            }
            glviewer = $3Dmol.createViewer(viewerEl, { 
                defaultcolors: $3Dmol.rasmolElementColors, 
                backgroundColor: '#060d1a',
                nomouse: false // Ensure mouse interaction is ON
            });
            glviewer.setClickable(true);

            // Handle global window resize
            window.addEventListener('resize', () => {
                if (glviewer) {
                    glviewer.resize();
                    glviewer.render();
                }
            });
        }
    }

    function loadPdbInViewer(pdbText, colorScheme = 'spectrum', label = 'Raw', viewer = null) {
        if (!viewer && !glviewer) initViewer();
        const targetViewer = viewer || glviewer;
        if (!targetViewer) return;

        // Force hide placeholder
        const ph = document.getElementById('viewer-placeholder');
        if (ph) ph.classList.add('hidden');
        if (placeholder) placeholder.classList.add('hidden');

        targetViewer.clear();
        targetViewer.removeAllModels();

        // PDBQT files from OpenBabel often lack the CONECT or secondary structure 
        // metadata required for 'cartoon' representation. We detect PDBQT and 
        // provide a line fallback to ensure it's always visible.
        let isPdbqt = false;
        const upText = pdbText.toUpperCase();
        if (upText.includes('ROOT') || upText.includes('BRANCH') || 
            upText.includes('REMARK VINA') || upText.includes('REMARK  5')) {
            isPdbqt = true;
        }

        try {
            // Always parse as standard PDB for maximum compatibility
            targetViewer.addModel(pdbText, 'pdb');
            
            if (isPdbqt) {
                // PDBQT Fallback: lines for backbone, sticks for everything else
                targetViewer.setStyle({}, { line: { color: 'white', opacity: 0.5 } });
                targetViewer.setStyle({ cartoon: { color: colorScheme } }, { cartoon: { color: colorScheme } });
            } else {
                // Standard PDB: cartoon for backbone
                targetViewer.setStyle({}, { cartoon: { color: colorScheme } });
            }
            
            // LIGAND RADIUS FIX: Ensure all non-protein residues (ligands/lipids/ions) 
            // are visible as thick sticks. We use a 'not' filter for common amino acids 
            // to catch everything else (like 2RH1's aromatic rings or 3PBR's lipids).
            const standardAmino = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HOH"];
            targetViewer.setStyle({ resn: standardAmino, invert: true }, { stick: { radius: 0.25, colorscheme: 'Jmol' } });
            
            // Explicit sphere fallback for ions
            targetViewer.setStyle({ hetflag: true, atom: 'CA,MG,ZN,FE,NA,CL' }, { sphere: { radius: 0.8 } });
            
            targetViewer.zoomTo();
            targetViewer.render();
            
            // Multi-pass resize to handle DOM layout shifts
            setTimeout(() => { if (targetViewer) { targetViewer.resize(); targetViewer.render(); } }, 100);
            setTimeout(() => { if (targetViewer) { targetViewer.resize(); targetViewer.render(); } }, 300);
            
            targetViewer.setClickable(true); 

            if (!viewer) {
                const sb = document.getElementById('state-badge');
                if (sb) {
                    sb.textContent = label;
                    sb.className = colorScheme === 'spectrum' ? 'badge badge-info' : 'badge badge-ok';
                    sb.classList.remove('hidden');
                }
            }
        } catch (err) {
            console.error("Error loading structure in viewer:", err);
            showStatus('precision-status', 'Failed to render structure in 3D viewer. Check console for details.', 'error');
        }
    }

    // Reset Viewer Button
    document.getElementById('reset-viewer-btn').addEventListener('click', () => {
        if (glviewer) {
            glviewer.zoomTo();
            glviewer.render();
            glviewer.resize();
            // Force a few more delayed resizes
            setTimeout(() => { glviewer.resize(); glviewer.render(); }, 100);
            setTimeout(() => { glviewer.resize(); glviewer.render(); }, 500);
        }
    });

    // Modal Viewer Implementation
    const modalViewerEl = document.getElementById('modal-3d-viewer');
    const viewerModal = document.getElementById('viewer-modal');

    function openModalViewer(pdbText, title = "Structure Preview") {
        document.getElementById('modal-viewer-title').textContent = title;
        viewerModal.classList.remove('hidden');
        document.getElementById('modal-file-sidebar').classList.add('hidden');

        if (!modalGlViewer) {
            modalGlViewer = $3Dmol.createViewer(modalViewerEl, { defaultcolors: $3Dmol.rasmolElementColors });
            modalGlViewer.setBackgroundColor('#000000');
        }
        loadPdbInViewer(pdbText, 'spectrum', '', modalGlViewer);
    }

    function closeModal() {
        viewerModal.classList.add('hidden');
    }

    document.getElementById('close-viewer-btn').onclick = closeModal;
    viewerModal.onclick = (e) => { if (e.target === viewerModal) closeModal(); };

    // ZIP Results Previewer
    async function previewZipResults(blob, title = "Preview Results") {
        const zip = await JSZip.loadAsync(blob);
        const pdbFiles = {};

        // Find all PDB files (excluding potential macOS __MACOSX meta files)
        const entries = Object.keys(zip.files).filter(f => f.endsWith('.pdb') && !f.includes('__MACOSX'));

        if (entries.length === 0) {
            alert("No PDB files found in the results ZIP.");
            return;
        }

        for (const filename of entries) {
            const content = await zip.files[filename].async("text");
            pdbFiles[filename] = content;
        }

        currentZipFiles = pdbFiles;
        renderZipFileList(entries, title);
    }

    function renderZipFileList(filenames, title) {
        // We use the MODAL for the zip list because it has a sidebar we built for it.
        // It's better than trying to cram it into the main UI for now.
        document.getElementById('modal-viewer-title').textContent = title;
        viewerModal.classList.remove('hidden');
        const sidebar = document.getElementById('modal-file-sidebar');
        const listContainer = document.getElementById('modal-file-list');
        sidebar.classList.remove('hidden');
        listContainer.innerHTML = '';

        if (!modalGlViewer) {
            modalGlViewer = $3Dmol.createViewer(modalViewerEl, { defaultcolors: $3Dmol.rasmolElementColors });
            modalGlViewer.setBackgroundColor('#000000');
        }

        filenames.forEach((fname, idx) => {
            const item = document.createElement('div');
            item.className = 'modal-file-item' + (idx === 0 ? ' active' : '');
            item.textContent = fname.split('/').pop(); // Show basename
            item.onclick = () => {
                document.querySelectorAll('.modal-file-item').forEach(el => el.classList.remove('active'));
                item.classList.add('active');
                loadPdbInViewer(currentZipFiles[fname], 'spectrum', '', modalGlViewer);
            };
            listContainer.appendChild(item);
        });

        // Load first file by default
        loadPdbInViewer(currentZipFiles[filenames[0]], 'spectrum', '', modalGlViewer);
    }

    function showDiff() {
        if (!precisionRawPdb || !precisionCleanPdb) return;
        initViewer();
        glviewer.clear();
        glviewer.removeAllModels();

        // 1. Raw Structure as translucent grey cartoon
        glviewer.addModel(precisionRawPdb, 'pdb');
        glviewer.setStyle({ model: 0 }, { cartoon: { color: '#555', opacity: 0.5 } });

        // 2. Processed Structure overlaid in spectrum
        // We now use precisionViewerPdb which is guaranteed to be a standard PDB
        glviewer.addModel(precisionViewerPdb, 'pdb');
        glviewer.setStyle({ model: 1 }, { cartoon: { color: 'spectrum' } });
        glviewer.setStyle({ model: 1, hetflag: true }, { stick: { colorscheme: 'Jmol' } });

        glviewer.zoomTo();
        glviewer.render();

        document.getElementById('state-badge').textContent = '⚖️ Diff Mode';
    }

    // ─────────────────────────────────────── //
    // STATUS MESSAGES
    // ─────────────────────────────────────── //
    function showStatus(elId, msg, type) {
        const el = document.getElementById(elId);
        el.textContent = msg;
        el.className = `status-msg status-${type}`;
        el.classList.remove('hidden');
    }
    function hideStatus(elId) {
        document.getElementById(elId).classList.add('hidden');
    }

    // ─────────────────────────────────────── //
    // ══════ PRECISION MODE ══════════════════
    // ─────────────────────────────────────── //

    // Drop zone precision - now a <label> so clicking works natively via for="file-input-precision"
    const fileInputP = document.getElementById('file-input-precision');
    const dropZoneP = document.getElementById('drop-zone-precision');

    // Drag-and-drop still needs JS
    ['dragenter', 'dragover'].forEach(ev => dropZoneP.addEventListener(ev, e => { e.preventDefault(); dropZoneP.classList.add('drag-active'); }));
    ['dragleave', 'drop'].forEach(ev => dropZoneP.addEventListener(ev, e => { e.preventDefault(); dropZoneP.classList.remove('drag-active'); }));
    dropZoneP.addEventListener('drop', e => {
        const file = e.dataTransfer.files[0];
        if (file && file.name.toLowerCase().endsWith('.pdb')) handlePrecisionFile(file);
    });
    fileInputP.addEventListener('change', function () {
        if (this.files[0]) handlePrecisionFile(this.files[0]);
    });

    function handlePrecisionFile(file) {
        precisionFile = file;
        document.getElementById('file-name-precision').textContent = `📄 ${file.name}`;
        document.getElementById('file-name-precision').classList.remove('hidden');

        // Cross-clear Site Analyzer manual state
        siteManualFile = null;
        const siteTag = document.getElementById('site-file-name');
        if (siteTag) siteTag.classList.add('hidden');

        // FULL RESET on new file (Sync)
        currentSessionId = null;
        precisionCleanPdb = null;
        precisionReport = null;
        diffMode = false;
        
        document.getElementById('report-panel').classList.add('hidden');
        document.getElementById('diff-mode-btn').classList.add('hidden');
        document.getElementById('diff-mode-btn').textContent = '⚖️ Diff Mode';
        document.getElementById('step-2-precision').classList.add('hidden');
        document.getElementById('atom-count-badge').classList.add('hidden');
        
        if (glviewer) {
            glviewer.clear();
            glviewer.render();
        }
        clearSiteViz();

        const reader = new FileReader();
        reader.onload = e => {
            precisionRawPdb = e.target.result;
            loadPdbInViewer(precisionRawPdb, 'spectrum', 'Raw Structure');
            console.log("Structure loaded:", file.name);
        };
        reader.readAsText(file);

        analyzeStructure(file);
    }

    // Analyze
    async function analyzeStructure(file) {
        showStatus('precision-status', 'Analyzing structure…', 'info');
        const fd = new FormData();
        fd.append('file', file);
        try {
            const res = await fetch('/api/analyze', { method: 'POST', body: fd });
            const data = await res.json();
            if (!res.ok) throw new Error(data.error || 'Analysis failed');
            populatePrecisionConfig(data.metadata);
            currentSessionId = data.session_id; // Store for other tools
            document.getElementById('step-2-precision').classList.remove('hidden');
            document.getElementById('atom-count-badge').textContent = `${data.metadata.atoms_total} Atoms`;
            document.getElementById('atom-count-badge').classList.remove('hidden');
            hideStatus('precision-status');
        } catch (e) {
            showStatus('precision-status', e.message, 'error');
        }
    }

    function populatePrecisionConfig(meta) {
        // Chains
        const cc = document.getElementById('chains-container');
        cc.innerHTML = '';
        if (meta.chains.length) {
            meta.chains.forEach(ch => {
                cc.innerHTML += `<div class="toggle-row"><input type="checkbox" id="chain-${ch}" class="chain-cb" value="${ch}" checked><label for="chain-${ch}">✅ Keep Chain ${ch}</label></div>`;
            });
        } else {
            cc.innerHTML = '<p class="empty-msg">No chains detected.</p>';
        }

        // Waters
        document.getElementById('water-count').textContent = meta.water_count;

        // Heteroatoms – single list, KEEP approach:
        // checked = keep this molecule, unchecked = remove it
        const hc = document.getElementById('heteros-container');
        hc.innerHTML = '';
        if (meta.heteroatoms.length) {
            meta.heteroatoms.forEach(h => {
                hc.innerHTML += `
                <div class="het-item">
                    <div class="toggle-row">
                        <input type="checkbox" id="keep-${h}" class="het-keep-cb" value="${h}" checked>
                        <label for="keep-${h}">
                            <span class="het-name">${h}</span>
                            <span class="het-action-badge" id="badge-het-${h}">✅ Will be KEPT</span>
                        </label>
                    </div>
                </div>`;
            });
            // Update badge when toggled
            hc.querySelectorAll('.het-keep-cb').forEach(cb => {
                cb.addEventListener('change', function () {
                    const badge = document.getElementById(`badge-het-${this.value}`);
                    badge.textContent = this.checked ? '✅ Will be KEPT' : '🗑 Will be REMOVED';
                    badge.style.color = this.checked ? 'var(--green)' : 'var(--red)';
                });
            });
        } else {
            hc.innerHTML = '<p class="empty-msg">No heteroatoms detected in this file.</p>';
        }

        // Missing residues notice
        if (meta.missing_residues && meta.missing_residues.length) {
            showStatus('precision-status', `⚠️ ${meta.missing_residues.length} missing residue(s) detected in the structure.`, 'info');
        }
    }

    // Minimization toggle
    document.getElementById('run-minimization').addEventListener('change', function () {
        document.getElementById('mini-settings').style.display = this.checked ? 'block' : 'none';
    });

    // pH Slider value display
    document.getElementById('ph-slider').addEventListener('input', function() {
        document.getElementById('ph-val').textContent = this.value;
    });

    // Remove-waters toggle → show/hide smart water sub-toggle
    document.getElementById('remove-waters').addEventListener('change', function () {
        document.getElementById('structural-water-row').style.opacity = this.checked ? '1' : '.3';
    });

    // Process button
    document.getElementById('process-btn').addEventListener('click', async () => {
        if (!precisionFile) return;

        const chains = Array.from(document.querySelectorAll('.chain-cb:checked')).map(c => c.value);
        // het-keep-cb: checked = KEEP, unchecked = REMOVE → send unchecked as remove list
        const allHets = Array.from(document.querySelectorAll('.het-keep-cb')).map(c => c.value);
        const keptHets = Array.from(document.querySelectorAll('.het-keep-cb:checked')).map(c => c.value);
        const removedHets = allHets.filter(h => !keptHets.includes(h));
        const protectedLig = keptHets; // kept ligands ARE the protected ones
        const removeWater = document.getElementById('remove-waters').checked;
        const keepStructural = document.getElementById('keep-structural-waters').checked;
        const ph = document.getElementById('ph-slider').value;
        const addAtoms = document.getElementById('add-missing-atoms').checked;
        const reconLoops = document.getElementById('reconstruct-loops').checked;
        const runMini = document.getElementById('run-minimization').checked;
        const useGbsa = document.getElementById('use-gbsa').checked;
        const ff = document.getElementById('force-field-select').value;
        const docking = document.querySelector('input[name="docking"]:checked')?.value || '';

        // UI setup
        const btn = document.getElementById('process-btn');
        const btnText = document.getElementById('process-btn-text');
        const spinner = document.getElementById('process-spinner');
        btn.disabled = true;
        btnText.textContent = 'Processing…';
        spinner.classList.remove('hidden');
        hideStatus('precision-status');

        const fd = new FormData();
        fd.append('file', precisionFile);
        fd.append('chains', JSON.stringify(chains));
        fd.append('remove_heteros', JSON.stringify(removedHets));
        fd.append('protect_ligands', JSON.stringify(protectedLig));
        fd.append('remove_water', removeWater);
        fd.append('keep_structural_waters', keepStructural);
        fd.append('ph', ph);
        fd.append('add_missing_atoms', addAtoms);
        fd.append('reconstruct_loops', reconLoops);
        fd.append('run_minimization', runMini);
        fd.append('use_gbsa', useGbsa);
        fd.append('force_field', ff);
        fd.append('docking_target', docking);

        try {
            const res = await fetch('/api/process', { method: 'POST', body: fd });
            const data = await res.json();
            
            if (!res.ok || !data.success) throw new Error(data.error || 'Processing failed');

            // HIGH PERFORMANCE DECODING
            const fastBase64ToText = (b64) => {
                const binary = atob(b64);
                const bytes = new Uint8Array(binary.length);
                for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
                return new TextDecoder().decode(bytes);
            };

            precisionCleanPdb = fastBase64ToText(data.pdb_b64);
            precisionViewerPdb = data.viewer_pdb_b64 ? fastBase64ToText(data.viewer_pdb_b64) : precisionCleanPdb;
            precisionReport = data.report;

            renderReport(data.report, data.report_text);
            updateRightStackVisibility('precision');
            loadPdbInViewer(precisionViewerPdb, 'spectrum', `✅ Processed (${data.report.protonation_ph})`);
            
            if (glviewer) {
                setTimeout(() => {
                    glviewer.zoomTo();
                    glviewer.render();
                    glviewer.resize();
                }, 200);
            }

            const delta = data.report.atom_counts.delta;
            const abc = document.getElementById('atom-count-badge');
            if (abc) {
                abc.textContent = `${data.report.atom_counts.after_processing} Atoms (${delta > 0 ? '+' : ''}${delta})`;
                abc.classList.remove('hidden');
            }

            showStatus('precision-status', '✅ Structure prepared successfully!', 'success');
        } catch (e) {
            console.error(e);
            showStatus('precision-status', `Error: ${e.message}`, 'error');
        } finally {
            btn.disabled = false;
            btnText.textContent = 'Process Structure';
            spinner.classList.add('hidden');
        }
    });





    // Diff mode toggle
    document.getElementById('diff-mode-btn').addEventListener('click', () => {
        diffMode = !diffMode;
        if (diffMode) {
            showDiff();
            document.getElementById('diff-mode-btn').textContent = '⬅ Exit Diff';
        } else {
            loadPdbInViewer(precisionViewerPdb, 'spectrum', '✅ Processed');
            document.getElementById('diff-mode-btn').textContent = '⚖️ Diff Mode';
        }
    });

    // Render Report Panel
    function renderReport(report, reportText) {
        const panel = document.getElementById('report-panel');
        panel.classList.remove('hidden');

        const grid = document.getElementById('report-grid');
        const ac = report.atom_counts;
        const chains = report.chains;
        const hets = report.heteroatoms;

        grid.innerHTML = '';
        // Reported by the protonator rather than assumed. The old markup
        // stated hydrogens were added whether or not they were.
        const prot = report.protonation || {};
        const hydrogensAdded = prot.hydrogens_added !== undefined
            ? prot.hydrogens_added : report.hydrogens_added;

        const items = [
            { label: 'Chains Detected', val: chains.detected.join(', ') || 'None' },
            { label: 'Chains Retained', val: chains.retained.join(', ') || 'None' },
            { label: 'Waters Removed', val: report.water_molecules_removed },
            { label: 'Waters Retained', val: report.water_molecules_retained || 0 },
            { label: 'Heteroatoms Removed', val: hets.removed.join(', ') || 'None' },
            { label: 'Ligands Retained', val: hets.retained.join(', ') || 'None' },
            { label: 'Hydrogens Added', val: hydrogensAdded ? '✅ Yes' : '❌ No' },
            { label: 'Target pH', val: report.protonation_ph },
            { label: 'Missing Residues', val: report.missing_residues_detected.length },
            { label: 'Atoms Before', val: ac.before_processing },
            { label: 'Atoms After', val: ac.after_processing },
            { label: 'Atom Delta', val: (ac.delta >= 0 ? '+' : '') + ac.delta },
            { label: 'Processing Time', val: `${report.processing_time_seconds}s` },
        ];

        if (prot.ligands_preserved && prot.ligands_preserved.length) {
            items.push({
                label: 'Ligands Protected',
                val: prot.ligands_preserved.join(', '),
            });
        }
        if (prot.loops_reconstructed) {
            items.push({ label: 'Loops Rebuilt', val: prot.loops_reconstructed });
        }
        if (prot.terminals_repaired) {
            items.push({ label: 'Terminals Repaired', val: prot.terminals_repaired });
        }
        if (prot.nonstandard_replaced && prot.nonstandard_replaced.length) {
            items.push({
                label: 'Nonstandard Replaced',
                val: prot.nonstandard_replaced.join(', '),
            });
        }
        items.forEach(({ label, val }) => {
            const div = document.createElement('div');
            div.className = 'report-item';
            div.innerHTML = `<div class="ri-label">${label}</div><div class="ri-val">${val}</div>`;
            grid.appendChild(div);
        });

        // Reproducibility tags
        const reproGrid = document.getElementById('repro-grid');
        reproGrid.innerHTML = '';
        const tags = [
            `remove_water: ${report.water_molecules_removed > 0}`,
            `ph: ${report.protonation_ph}`,
            `chains: [${chains.retained.join(',')}]`,
            `removed_hets: [${hets.removed.join(',')}]`,
            `hydrogens: ${!!hydrogensAdded}`,
        ];
        if (report.energy_minimization) tags.push(`force_field: ${report.energy_minimization.force_field}`);
        if (report.docking_target) tags.push(`docking: ${report.docking_target}`);
        tags.forEach(t => {
            const span = document.createElement('span');
            span.className = 'repro-tag';
            span.textContent = t;
            reproGrid.appendChild(span);
        });

        // Energy minimization section
        if (report.energy_minimization) {
            const em = report.energy_minimization;
            const miniPanel = document.getElementById('mini-report-panel');
            miniPanel.classList.remove('hidden');
            const miniGrid = document.getElementById('mini-grid');
            miniGrid.innerHTML = '';
            
            if (em.error) {
                // Display the error gracefully instead of N/A
                const d = document.createElement('div');
                d.className = 'report-item' + (em.error ? ' error-item' : '');
                d.style.gridColumn = "1 / -1"; // Span full width
                d.innerHTML = `<div class="ri-label" style="color: var(--red);">Minimization Failed</div>
                               <div class="ri-val" style="font-size: 0.85rem; white-space: normal;">${em.error}</div>`;
                miniGrid.appendChild(d);
            } else {
                // Energies are numbers or null now, never the string 'N/A'.
                const kJ = (v) => (v === null || v === undefined)
                    ? 'not run' : `${v} kJ/mol`;

                // status is the field that matters: a returned file is not
                // proof that minimization ran. 'partial' means some residues
                // were held at their input coordinates.
                const STATUS_LABEL = {
                    full: 'Full structure',
                    partial: 'Partial — some residues held fixed',
                    partial_no_implicit_solvent: 'Partial, no implicit solvent',
                    failed: 'Did not run',
                };

                const miniItems = [
                    { label: 'Status', val: STATUS_LABEL[em.status] || em.status },
                    { label: 'Force Field', val: em.force_field },
                    { label: 'Energy Before', val: kJ(em.energy_before_kJ_mol) },
                    { label: 'Energy After', val: kJ(em.energy_after_kJ_mol) },
                    { label: 'Max Iterations', val: em.iterations_max },
                    {
                        label: 'Converged',
                        val: em.converged ? '✅ Yes' : '❌ No'
                            + (em.rms_force_kJ_mol_nm !== null
                                && em.rms_force_kJ_mol_nm !== undefined
                                ? ` (RMS force ${em.rms_force_kJ_mol_nm} kJ/mol/nm)` : ''),
                    },
                ];

                if (em.excluded_residues && em.excluded_residues.length) {
                    miniItems.push({
                        label: 'Held at input coordinates',
                        val: em.excluded_residues.join(', '),
                    });
                }
                if (em.restrained_atoms) {
                    miniItems.push({
                        label: 'Restrained pocket atoms',
                        val: em.restrained_atoms,
                    });
                }
                miniItems.forEach(({ label, val }) => {
                    const d = document.createElement('div');
                    d.className = 'report-item';
                    d.innerHTML = `<div class="ri-label">${label}</div><div class="ri-val">${val}</div>`;
                    miniGrid.appendChild(d);
                });
            }
        }

        // Warnings. The pipeline records things the user needs to know about
        // the result - a ligand held at its input coordinates, a repaired
        // terminus, a suspicious starting energy - and nothing displayed them.
        const warningsPanel = document.getElementById('warnings-panel');
        const warningsList = document.getElementById('warnings-list');
        const warnings = report.warnings || [];
        warningsList.innerHTML = '';
        if (warnings.length) {
            warningsPanel.classList.remove('hidden');
            warnings.forEach(text => {
                const li = document.createElement('li');
                li.textContent = text;
                warningsList.appendChild(li);
            });
        } else {
            warningsPanel.classList.add('hidden');
        }

        // Download PDB/PDBQT button
        document.getElementById('download-pdb-btn').onclick = () => {
            const isPdbqt = report.docking_target && ['vina', 'autodock'].includes(report.docking_target.toLowerCase());
            const ext = isPdbqt ? '.pdbqt' : '.pdb';
            const blob = new Blob([precisionCleanPdb], { type: isPdbqt ? 'text/plain' : 'chemical/x-pdb' });
            triggerDownload(blob, report.input_file.replace('.pdb', ext));
        };

        // Download report
        document.getElementById('download-report-btn').onclick = async () => {
            try {
                const res = await fetch('/api/download-report', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify(report),
                });
                const blob = await res.blob();
                triggerDownload(blob, 'preparation_report.txt');
            } catch (e) {
                alert('Failed to download report.');
            }
        };

        // Copy report to clipboard
        document.getElementById('copy-report-btn').onclick = () => {
            navigator.clipboard.writeText(reportText);
            document.getElementById('copy-report-btn').textContent = '✅ Copied';
            setTimeout(() => { document.getElementById('copy-report-btn').textContent = '📋 Copy'; }, 2000);
        };
    }

    function refreshTemplateSelect() {
        const sel = document.getElementById('template-select');
        if (!sel) return;
        sel.innerHTML = '<option value="">-- No template --</option>';
        Object.keys(templatesStore).forEach(name => {
            sel.innerHTML += `<option value="${name}">${name}</option>`;
        });
    }

    window.applyTemplate = function (name) {
        const s = templatesStore[name]?.settings;
        if (!s) return;

        if (s.ph !== undefined) {
            document.getElementById('ph-slider').value = s.ph;
            document.getElementById('ph-val').textContent = s.ph;
        }
        if (s.remove_water !== undefined) document.getElementById('remove-waters').checked = s.remove_water;
        if (s.keep_structural_waters !== undefined) document.getElementById('keep-structural-waters').checked = s.keep_structural_waters;
        if (s.add_missing_atoms !== undefined) document.getElementById('add-missing-atoms').checked = s.add_missing_atoms;
        if (s.reconstruct_loops !== undefined) document.getElementById('reconstruct-loops').checked = s.reconstruct_loops;
        if (s.run_minimization !== undefined) {
            document.getElementById('run-minimization').checked = s.run_minimization;
            document.getElementById('mini-settings').style.display = s.run_minimization ? 'block' : 'none';
        }
        if (s.use_gbsa !== undefined) document.getElementById('use-gbsa').checked = s.use_gbsa;
        if (s.force_field) document.getElementById('force-field-select').value = s.force_field;

        // FIX: Always reset ALL docking radios first, then set the correct one.
        // Previously, if docking_target was '' (none), the old radio stayed checked.
        document.querySelectorAll('input[name="docking"]').forEach(r => r.checked = false);
        if (s.docking_target) {
            const radio = document.querySelector(`input[name="docking"][value="${s.docking_target}"]`);
            if (radio) radio.checked = true;
        } else {
            // Explicitly select "None"
            const noneRadio = document.querySelector('input[name="docking"][value=""]');
            if (noneRadio) noneRadio.checked = true;
        }

        showStatus('precision-status', `⚡ Preset "${name}" applied!`, 'success');
        setTimeout(() => hideStatus('precision-status'), 3000);
    };


    window.deleteTemplate = async function (name) {
        await fetch(`/api/templates/${encodeURIComponent(name)}`, { method: 'DELETE' });
        await loadTemplates();
    };

    document.getElementById('load-template-btn').addEventListener('click', () => {
        const name = document.getElementById('template-select').value;
        if (name) window.applyTemplate(name);
        else alert('Select a preset from the dropdown first.');
    });

    document.getElementById('save-template-btn').addEventListener('click', async () => {
        const name = document.getElementById('new-template-name').value.trim();
        if (!name) {
            alert('Type a name for the preset first (e.g. "Docking Prep")');
            return;
        }
        const settings = getCurrentSettings();
        await fetch('/api/templates', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ name, settings }),
        });
        await loadTemplates();
        document.getElementById('new-template-name').value = '';
        document.getElementById('template-select').value = name;
        showStatus('precision-status', `✅ Preset "${name}" saved!`, 'success');
        setTimeout(() => hideStatus('precision-status'), 3000);
    });

    document.getElementById('export-template-btn').addEventListener('click', () => {
        const name = document.getElementById('template-select').value;
        if (!name) {
            alert('Select a preset to export first.');
            return;
        }
        const settings = templatesStore[name];
        if (!settings) return;
        
        const blob = new Blob([JSON.stringify({ name, settings }, null, 2)], { type: 'application/json' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `BioPrep_Preset_${name.replace(/\s+/g, '_')}.json`;
        a.click();
        URL.revokeObjectURL(url);
    });

    document.getElementById('import-template-btn').addEventListener('click', () => {
        document.getElementById('import-template-input').click();
    });

    document.getElementById('import-template-input').addEventListener('change', async (e) => {
        const file = e.target.files[0];
        if (!file) return;
        
        const reader = new FileReader();
        reader.onload = async (ev) => {
            try {
                const data = JSON.parse(ev.target.result);
                if (!data.name || !data.settings) {
                    throw new Error("Invalid preset file format. Must contain 'name' and 'settings'.");
                }
                
                await fetch('/api/templates', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ name: data.name, settings: data.settings }),
                });
                await loadTemplates();
                showStatus('precision-status', `✅ Preset "${data.name}" imported!`, 'success');
            } catch (err) {
                alert(`Import failed: ${err.message}`);
            }
        };
        reader.readAsText(file);
        e.target.value = ''; // Reset input
    });

    function getCurrentSettings() {
        return {
            ph: parseFloat(document.getElementById('ph-slider').value),
            remove_water: document.getElementById('remove-waters').checked,
            keep_structural_waters: document.getElementById('keep-structural-waters').checked,
            add_missing_atoms: document.getElementById('add-missing-atoms').checked,
            reconstruct_loops: document.getElementById('reconstruct-loops').checked,
            run_minimization: document.getElementById('run-minimization').checked,
            use_gbsa: document.getElementById('use-gbsa').checked,
            force_field: document.getElementById('force-field-select').value,
            docking_target: document.querySelector('input[name="docking"]:checked')?.value || '',
        };
    }

    // ─────────────────────────────────────── //
    // ══════ BATCH MODE ══════════════════════
    // ─────────────────────────────────────── //
    const dropZoneB = document.getElementById('drop-zone-batch');
    const fileInputB = document.getElementById('file-input-batch');

    // batch multi-file drop
    dropZoneB.addEventListener('click', () => fileInputB.click());
    fileInputB.addEventListener('change', () => handleBatchFiles(fileInputB.files));
    ['dragenter', 'dragover'].forEach(e => dropZoneB.addEventListener(e, ev => { ev.preventDefault(); dropZoneB.classList.add('drag-active'); }));
    ['dragleave', 'drop'].forEach(e => dropZoneB.addEventListener(e, ev => { ev.preventDefault(); dropZoneB.classList.remove('drag-active'); }));
    dropZoneB.addEventListener('drop', e => handleBatchFiles(e.dataTransfer.files));

    function handleBatchFiles(files) {
        const newFiles = Array.from(files).filter(f => f.name.toLowerCase().endsWith('.pdb'));
        batchFiles = [...batchFiles, ...newFiles];
        const list = document.getElementById('batch-file-list');
        list.innerHTML = '';
        batchFiles.forEach(f => {
            const pill = document.createElement('div');
            pill.className = 'file-pill';
            pill.innerHTML = `<span>${f.name}</span> <span class="view-file-btn" title="View Structure">👁</span>`;

            pill.querySelector('.view-file-btn').onclick = () => {
                const reader = new FileReader();
                reader.onload = e => {
                    initViewer();
                    loadPdbInViewer(e.target.result, 'spectrum', f.name);
                };
                reader.readAsText(f);
            };

            list.appendChild(pill);
        });

        // Automatically load first file from the batch into the viewer
        if (batchFiles.length > 0) {
            const first = batchFiles[0];
            const reader = new FileReader();
            reader.onload = e => {
                loadPdbInViewer(e.target.result, 'spectrum', first.name);
            };
            reader.readAsText(first);
        }
    }

    document.getElementById('b-ph-slider').addEventListener('input', function () {
        document.getElementById('b-ph-val').textContent = this.value;
    });

    document.getElementById('batch-process-btn').addEventListener('click', async () => {
        if (!batchFiles.length) { showStatus('batch-status', 'No PDB files selected.', 'error'); return; }

        const btn = document.getElementById('batch-process-btn');
        const btnText = document.getElementById('batch-btn-text');
        const spinner = document.getElementById('batch-spinner');
        btn.disabled = true;
        btnText.textContent = 'Processing…';
        spinner.classList.remove('hidden');
        hideStatus('batch-status');

        // Show fake progress animation
        const progressWrapper = document.getElementById('batch-progress-wrapper');
        const progressBar = document.getElementById('batch-progress-bar');
        const progressLabel = document.getElementById('batch-progress-label');
        progressWrapper.classList.remove('hidden');
        progressLabel.textContent = `Processing ${batchFiles.length} files…`;

        let fakeProgress = 0;
        const fakeTimer = setInterval(() => {
            fakeProgress = Math.min(fakeProgress + (100 / batchFiles.length / 3), 90);
            progressBar.style.width = `${fakeProgress}%`;
        }, 400);

        const config = {
            remove_water: document.getElementById('b-remove-waters').checked,
            remove_heteros: document.getElementById('b-remove-hets').checked ? ['ALL'] : [],
            ph: parseFloat(document.getElementById('b-ph-slider').value),
            run_minimization: document.getElementById('b-run-mini').checked,
            use_gbsa: document.getElementById('b-use-gbsa').checked,
            reconstruct_loops: document.getElementById('b-reconstruct').checked,
            add_missing_atoms: document.getElementById('b-add-atoms').checked,
            docking_target: document.getElementById('b-docking').value,
        };

        const fd = new FormData();
        fd.append('config', JSON.stringify(config));
        batchFiles.forEach(f => fd.append('files', f));

        try {
            const res = await fetch('/api/batch', { method: 'POST', body: fd });
            clearInterval(fakeTimer);
            progressBar.style.width = '100%';

            if (!res.ok) {
                const err = await res.json();
                throw new Error(err.error || 'Batch failed');
            }
            batchResultBlob = await res.blob();
            progressLabel.textContent = `✅ ${batchFiles.length} files processed!`;

            // Show summary
            updateRightStackVisibility('batch');
            const summaryGrid = document.getElementById('batch-summary');
            summaryGrid.innerHTML = '';
            const items = [
                { label: 'Files Uploaded', val: batchFiles.length },
                { label: 'Output', val: 'bioprep_batch_results.zip' },
            ];
            items.forEach(({ label, val }) => {
                summaryGrid.innerHTML += `<div class="report-item"><div class="ri-label">${label}</div><div class="ri-val">${val}</div></div>`;
            });

            // Automatically preview results
            previewZipResults(batchResultBlob, "Batch Processing Results");

        } catch (e) {
            clearInterval(fakeTimer);
            showStatus('batch-status', `Error: ${e.message}`, 'error');
            progressWrapper.classList.add('hidden');
        } finally {
            btn.disabled = false;
            btnText.textContent = 'Process All Files';
            spinner.classList.add('hidden');
        }
    });

    document.getElementById('download-batch-btn').addEventListener('click', () => {
        if (batchResultBlob) triggerDownload(batchResultBlob, 'bioprep_batch_results.zip');
    });

    document.getElementById('preview-batch-btn').onclick = () => {
        if (batchResultBlob) previewZipResults(batchResultBlob, "Batch Processing Results");
    };

    // ─────────────────────────────────────── //
    // ══════ HIGH-THROUGHPUT MODE ════════════
    // ─────────────────────────────────────── //
    const dropZoneHT = document.getElementById('drop-zone-ht');
    const fileInputHT = document.getElementById('file-input-ht');

    setupDropZone(dropZoneHT, fileInputHT, '.zip', (file) => {
        htFile = file;
        htResultBlob = null; // Clear old result
        document.getElementById('ht-results').classList.add('hidden');
        document.getElementById('ht-file-tag').textContent = `📦 ${file.name}`;
        document.getElementById('ht-file-tag').classList.remove('hidden');
        hideStatus('ht-status');
    });

    document.getElementById('ht-ph-slider').addEventListener('input', function () {
        document.getElementById('ht-ph-val').textContent = this.value;
    });

    document.getElementById('b-run-mini').addEventListener('change', function () {
        const row = document.getElementById('b-gbsa-row');
        if (row) row.style.display = this.checked ? 'block' : 'none';
    });

    document.getElementById('ht-run-mini').addEventListener('change', function () {
        const row = document.getElementById('ht-gbsa-row');
        if (row) row.style.display = this.checked ? 'block' : 'none';
    });

    document.getElementById('ht-process-btn').addEventListener('click', async () => {
        if (!htFile) { showStatus('ht-status', 'No ZIP file selected.', 'error'); return; }

        const btn = document.getElementById('ht-process-btn');
        const btnText = document.getElementById('ht-btn-text');
        const spinner = document.getElementById('ht-spinner');
        btn.disabled = true;
        btnText.textContent = 'Running Pipeline…';
        spinner.classList.remove('hidden');
        hideStatus('ht-status');

        const progressWrapper = document.getElementById('ht-progress-wrapper');
        progressWrapper.classList.remove('hidden');
        document.getElementById('ht-log-panel').classList.add('hidden');

        const config = {
            remove_water: document.getElementById('ht-remove-waters').checked,
            remove_heteros: document.getElementById('ht-remove-hets').checked ? ['ALL'] : [],
            ph: parseFloat(document.getElementById('ht-ph-slider').value),
            run_minimization: document.getElementById('ht-run-mini').checked,
            use_gbsa: document.getElementById('ht-use-gbsa').checked,
            reconstruct_loops: document.getElementById('ht-reconstruct').checked,
            add_missing_atoms: document.getElementById('ht-add-atoms').checked,
            docking_target: document.getElementById('ht-docking').value,
        };

        const fd = new FormData();
        fd.append('file', htFile);
        fd.append('config', JSON.stringify(config));

        try {
            const res = await fetch('/api/high-throughput', { method: 'POST', body: fd });
            progressWrapper.classList.add('hidden');

            if (!res.ok) {
                const err = await res.json();
                throw new Error(err.error || 'Pipeline failed');
            }
            htResultBlob = await res.blob();

            updateRightStackVisibility('highthroughput');
            showStatus('ht-status', '✅ High-throughput pipeline complete!', 'success');

            // Automatically preview results
            previewZipResults(htResultBlob, "HT Pipeline Results");

        } catch (e) {
            progressWrapper.classList.add('hidden');
            showStatus('ht-status', `Error: ${e.message}`, 'error');
        } finally {
            btn.disabled = false;
            btnText.textContent = 'Start Pipeline';
            spinner.classList.add('hidden');
        }
    });

    document.getElementById('download-ht-btn').addEventListener('click', () => {
        if (htResultBlob) triggerDownload(htResultBlob, 'bioprep_ht_results.zip');
    });

    document.getElementById('preview-ht-btn').onclick = () => {
        if (htResultBlob) previewZipResults(htResultBlob, "HT Pipeline Results");
    };

    // ─────────────────────────────────────── //
    // UTILS
    // ─────────────────────────────────────── //
    function triggerDownload(blob, filename) {
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        setTimeout(() => URL.revokeObjectURL(url), 5000);
    }

    function setupDropZone(zone, input, ext, onFile) {
        zone.addEventListener('click', () => input.click());
        ['dragenter', 'dragover'].forEach(ev => {
            zone.addEventListener(ev, e => { e.preventDefault(); zone.classList.add('drag-active'); });
        });
        ['dragleave', 'drop'].forEach(ev => {
            zone.addEventListener(ev, e => { e.preventDefault(); zone.classList.remove('drag-active'); });
        });
        zone.addEventListener('drop', e => {
            const file = e.dataTransfer.files[0];
            if (file && file.name.toLowerCase().endsWith(ext)) onFile(file);
        });
        input.addEventListener('change', function () {
            if (this.files[0]) onFile(this.files[0]);
        });
    }

    // ── SITE ANALYZER LOGIC ──────────────────────────────────
    window.runSiteAnalysis = async function () {
        const btn = document.getElementById('analyze-site-btn');
        const spinner = document.getElementById('site-spinner');
        const status = document.getElementById('site-status');

        if (!siteManualFile && !currentSessionId) {
            showStatus('site-status', 'Please upload a structure first.', 'error');
            return;
        }

        const fileToUse = siteManualFile;

        btn.disabled = true;
        if (spinner) spinner.classList.remove('hidden');
        showStatus('site-status', 'Analyzing binding sites and pharmacophores...', 'info');

        try {
            const fd = new FormData();
            let url = `/api/analyze-site`;
            if (fileToUse) {
                fd.append('file', fileToUse);
            } else {
                url += `/${currentSessionId}`;
            }

            const response = await fetch(url, { 
                method: 'POST', 
                body: fileToUse ? fd : null 
            });
            const data = await response.json();

            if (data.success) {
                discoveredSites = data.sites;
                renderSiteList(data.sites);
                if (data.summary) renderSiteSummary(data.summary);
                
                document.getElementById('site-layers-panel').classList.remove('hidden');
                document.getElementById('site-results').classList.remove('hidden');
                showStatus('site-status', data.message, 'success');
                updateSiteViz(); // Initial render
            } else {
                throw new Error(data.error || 'Analysis failed');
            }
        } catch (err) {
            showStatus('site-status', err.message, 'error');
        } finally {
            btn.disabled = false;
            if (spinner) spinner.classList.add('hidden');
        }
    };

    function renderSiteSummary(summary) {
        const panel = document.getElementById('site-summary-panel');
        const textEl = document.getElementById('site-summary-text');
        const statsEl = document.getElementById('site-summary-stats');
        
        if (!panel || !summary) return;
        
        panel.classList.remove('hidden');
        textEl.innerHTML = summary.text.replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>');
        
        const siteCount = summary.site_count || 0;
        const primaryVol = (summary.primary_volume || 0).toFixed(1);
        const totalVol = (summary.total_volume || 0).toFixed(1);
        const features = summary.total_features || 0;

        statsEl.innerHTML = `
            <div class="report-item">
                <div class="ri-label">Total Sites</div>
                <div class="ri-val">${siteCount}</div>
            </div>
            <div class="report-item">
                <div class="ri-label">Primary Vol</div>
                <div class="ri-val">${primaryVol} Å³</div>
            </div>
            <div class="report-item">
                <div class="ri-label">Total Space</div>
                <div class="ri-val">${totalVol} Å³</div>
            </div>
            <div class="report-item">
                <div class="ri-label">Pharmacophores</div>
                <div class="ri-val">${features}</div>
            </div>
        `;
    }

    function renderSiteList(sites) {
        const list = document.getElementById('site-list');
        list.innerHTML = '';

        if (sites.length === 0) {
            list.innerHTML = '<p class="muted">No significant pockets detected.</p>';
            return;
        }

        sites.forEach(site => {
            const card = document.createElement('div');
            card.className = 'site-card';
            card.innerHTML = `
                <div class="site-card-header">
                    <strong>Site #${site.id}</strong>
                    <span class="score-pill">Score: ${site.drugability_score}</span>
                </div>
                <div class="site-card-meta">
                    <span>Vol: ${site.volume} Å³</span>
                    <span>Res: ${site.residues.length}</span>
                </div>
                <div class="site-card-coords" style="font-size: 0.75rem; color: var(--cyan-2); margin-top: 5px; font-family: monospace;">
                    [${(site.centroid[0] || 0).toFixed(1)}, ${(site.centroid[1] || 0).toFixed(1)}, ${(site.centroid[2] || 0).toFixed(1)}]
                </div>
            `;
            card.onclick = (e) => focusSite(site, e.currentTarget);
            list.appendChild(card);
        });
    }

    window.focusSite = function (site, element) {
        if (!glviewer) return;
        // Gently center on the site without extreme zoom
        // Use a selection sphere around the centroid so we keep structural context
        const cx = site.centroid[0], cy = site.centroid[1], cz = site.centroid[2];
        glviewer.center({ x: cx, y: cy, z: cz }, 600);
        glviewer.render();
        
        const info = document.getElementById('active-site-info');
        info.innerHTML = `
            <div id="active-site-details">
                <h5>Site #${site.id} Details</h5>
                <div class="site-detail-body">
                    <p class="coords-display">
                        <strong>Coordinates (X, Y, Z):</strong> 
                        <code>${(site.centroid[0] || 0).toFixed(2)}, ${(site.centroid[1] || 0).toFixed(2)}, ${(site.centroid[2] || 0).toFixed(2)}</code>
                    </p>
                    <p><strong>Top Residues:</strong> ${site.residues.join(', ')}</p>
                    <div class="prop-tags">
                        <span class="tag tag-hydro">H-Phob: ${site.properties.HYDROPHOBIC}</span>
                        <span class="tag tag-polar">Polar: ${site.properties.POLAR}</span>
                        <span class="tag tag-charged">Charged: ${site.properties.CHARGED}</span>
                    </div>
                </div>
            </div>
        `;
        
        document.querySelectorAll('.site-card').forEach(c => c.classList.remove('active'));
        if (element) element.classList.add('active');
    };

    window.updateSiteViz = function () {
        if (!glviewer) return;
        siteShapes.forEach(s => glviewer.removeShape(s));
        siteShapes.length = 0;

        const showSpheres = document.getElementById('show-site-spheres').checked;
        const showAcc = document.getElementById('show-site-acc').checked;
        const showDon = document.getElementById('show-site-don').checked;
        const showHphob = document.getElementById('show-site-hphob').checked;

        discoveredSites.forEach(site => {
            if (showSpheres) {
                // Main site bounding sphere - make it a soft glowing neon cyan
                const s = glviewer.addSphere({
                    center: { x: site.centroid[0], y: site.centroid[1], z: site.centroid[2] },
                    radius: 4.5, // Slightly larger for better prominence
                    color: '#00ffff', // Neon Cyan
                    alpha: 0.3 // Increased opacity
                });
                siteShapes.push(s);
            }

            site.pharmacophore_points.forEach(p => {
                let color = '#ffffff'; // Default bright white
                let visible = false;
                
                // Use bright neon colors for pharmacophores
                if (p.type === 'ACCEPTOR') { 
                    color = '#ff0055'; // Vibrant Neon Red/Pink
                    visible = showAcc; 
                }
                else if (p.type === 'DONOR') { 
                    color = '#00bbff'; // Bright Neon Blue
                    visible = showDon; 
                }
                else if (p.type === 'HYDROPHOBIC') { 
                    color = '#ccff00'; // Neon Lime/Yellow
                    visible = showHphob; 
                }

                if (visible) {
                    const sp = glviewer.addSphere({
                        center: { x: p.coords[0], y: p.coords[1], z: p.coords[2] },
                        radius: 0.7, // Slightly larger points
                        color: color,
                        alpha: 1.0, // Full opacity for pharmacophores
                        clickable: true,
                        hoverable: true,
                        hover_callback: function() {
                            if (!glviewer._hoverLabels) glviewer._hoverLabels = [];
                            const lbl = glviewer.addLabel(p.type, {
                                position: { x: p.coords[0], y: p.coords[1] + 1.0, z: p.coords[2] },
                                backgroundColor: 'black',
                                fontColor: 'white',
                                backgroundOpacity: 0.8,
                                inFront: true,
                                fontSize: 12
                            });
                            glviewer._hoverLabels.push(lbl);
                            glviewer.render();
                        },
                        unhover_callback: function() {
                            if (glviewer._hoverLabels) {
                                glviewer._hoverLabels.forEach(l => glviewer.removeLabel(l));
                                glviewer._hoverLabels = [];
                                glviewer.render();
                            }
                        }
                    });
                    siteShapes.push(sp);
                }
            });
        });
        glviewer.render();
    };

    window.clearSiteViz = function () {
        if (!glviewer) return;
        siteShapes.forEach(s => glviewer.removeShape(s));
        siteShapes.length = 0;
        discoveredSites = [];
        document.getElementById('site-list').innerHTML = '';
        document.getElementById('site-results').classList.add('hidden');
        document.getElementById('site-layers-panel').classList.add('hidden');
        const summaryPanel = document.getElementById('site-summary-panel');
        if (summaryPanel) summaryPanel.classList.add('hidden');
        glviewer.render();
    };

    // ── INIT ───────────────────────────────────────
    setupDropZone(
        document.getElementById('drop-zone-site-manual'),
        document.getElementById('file-input-site-manual'),
        '.pdb',
        (file) => {
            // FULL RESET for Site Analyzer results
            clearSiteViz();
            
            // Cross-clear Precision state
            precisionFile = null;
            precisionRawPdb = null;
            precisionCleanPdb = null;
            currentSessionId = null;
            const pTag = document.getElementById('file-name-precision');
            if (pTag) pTag.classList.add('hidden');
            document.getElementById('step-2-precision').classList.add('hidden');

            siteManualFile = file;
            const tag = document.getElementById('site-file-name');
            tag.textContent = `🧬 ${file.name}`;
            tag.classList.remove('hidden');
            // Clear input so same file can be re-selected
            document.getElementById('file-input-site-manual').value = '';
            
            // Load into visualizer immediately
            const reader = new FileReader();
            reader.onload = e => {
                loadPdbInViewer(e.target.result, 'spectrum', `🧬 ${file.name}`);
            };
            reader.readAsText(file);
            
            showStatus('site-status', 'Structure loaded. Ready for discovery.', 'info');
        }
    );

    // ─────────────────────────────────────── //
    // JOB HISTORY
    // ─────────────────────────────────────── //
    async function loadHistory() {
        const list = document.getElementById('history-list');
        try {
            const res = await fetch('/api/history');
            const data = await res.json();
            renderHistory(data);
        } catch (e) {
            console.error("Failed to load history:", e);
            list.innerHTML = '<div class="empty-state">Failed to load history items.</div>';
        }
    }

    function renderHistory(data) {
        const list = document.getElementById('history-list');
        list.innerHTML = '';

        if (!data || data.length === 0) {
            list.innerHTML = `
                <div class="empty-state">
                    <p>No previous jobs found. Start processing proteins to build your history!</p>
                </div>
            `;
            return;
        }

        data.forEach(job => {
            const card = document.createElement('div');
            card.className = 'history-card glass-panel hover-scale';
            
            // Energies are numbers or null. The old check compared against the
            // string 'N/A', so a null pair passed it and rendered "0.00".
            const haveEnergies = typeof job.energy_after === 'number'
                && typeof job.energy_before === 'number';
            const energyDelta = haveEnergies
                ? (job.energy_after - job.energy_before).toFixed(2)
                : 'not minimized';
            
            card.innerHTML = `
                <div class="history-card-header">
                    <span class="history-filename" title="${job.filename}">${job.filename}</span>
                    <span class="chip chip-${job.format && job.format.toLowerCase().includes('pdb') ? 'blue' : 'purple'}">${job.format || 'PDB'}</span>
                </div>
                <div class="history-meta">
                    <div class="meta-item">📅 ${job.timestamp}</div>
                    <div class="meta-item">📉 Δ Energy: <span class="${energyDelta < 0 ? 'energy-good' : ''}">${energyDelta} kJ/mol</span></div>
                </div>
                <div class="history-actions">
                    <button class="btn-ghost small" onclick="viewHistoryReport('${job.id}')">View Report</button>
                    <button class="btn-primary small" onclick="openPdbInViewer('${job.id}')">3D View</button>
                </div>
            `;
            // Store full job data in memory-safe store instead of DOM dataset
            historyStore[job.id] = job;
            list.appendChild(card);
        });
    }

    window.viewHistoryReport = function(jobId) {
        const job = historyStore[jobId];
        if (!job) return;
        
        if (job.report) {
            precisionReport = job.report; // Global for updateRightStackVisibility
            renderReport(job.report, job.report_text || "Report Loaded from History");
            
            // Switch tab to show the report panel
            document.getElementById('tab-precision').click();
            showStatus('precision-status', `Loaded report for ${job.filename}`, 'success');
        } else {
            alert('No report data found for this job.');
        }
    };


    window.openPdbInViewer = async function(jobId) {
        // Show loading on the History tab's status, not Precision, since user is on History
        const historyStatus = document.getElementById('history-status-msg');
        
        try {
            const res = await fetch(`/api/history/pdb/${jobId}`);
            if (!res.ok) {
                const msg = await res.text();
                // FIX: Show clear, friendly error directly on the history tab
                // so user doesn't get confused by navigating to Precision tab for an error.
                if (historyStatus) {
                    historyStatus.textContent = '⚠️ Structure unavailable after server restart. Re-process the file to view it in 3D.';
                    historyStatus.className = 'status-msg status-error';
                    historyStatus.classList.remove('hidden');
                    setTimeout(() => historyStatus.classList.add('hidden'), 6000);
                } else {
                    alert('Structure no longer available in session memory. Re-process the file to view it in 3D.');
                }
                return;
            }

            const pdbText = await res.text();
            loadPdbInViewer(pdbText, 'spectrum', `📜 Historic Structure`);

            // Switch to precision tab to show viewer
            document.getElementById('tab-precision').click();
            showStatus('precision-status', 'Historical structure loaded in viewer.', 'success');
        } catch (e) {
            console.error("Failed to load historical PDB:", e);
            if (historyStatus) {
                historyStatus.textContent = `❌ Could not load structure: ${e.message}`;
                historyStatus.className = 'status-msg status-error';
                historyStatus.classList.remove('hidden');
                setTimeout(() => historyStatus.classList.add('hidden'), 6000);
            }
        }
    };

    document.getElementById('refresh-history-btn').addEventListener('click', loadHistory);
    document.getElementById('tab-history').addEventListener('click', loadHistory);

    // ── INIT ───────────────────────────────────────
    loadTemplates();
    loadHistory();
});
