/* BioPrep frontend. Talks to the Flask API documented in BACKEND_AUDIT.md. */
(() => {
  "use strict";

  // ── tiny utilities ────────────────────────────────────────────────────

  const qs = (sel, root = document) => root.querySelector(sel);
  const qsa = (sel, root = document) => Array.from(root.querySelectorAll(sel));
  const el = (tag, cls, html) => {
    const node = document.createElement(tag);
    if (cls) node.className = cls;
    if (html !== undefined) node.innerHTML = html;
    return node;
  };
  const escapeHtml = (s) =>
    String(s).replace(/[&<>"']/g, (c) => ({
      "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
    }[c]));
  const fmtNum = (n, digits = 1) =>
    typeof n === "number" ? n.toLocaleString(undefined, { maximumFractionDigits: digits }) : "—";

  function toast(message, kind = "") {
    const stack = qs("#toast-stack");
    const node = el("div", `toast${kind ? " toast-" + kind : ""}`, escapeHtml(message));
    stack.appendChild(node);
    setTimeout(() => node.remove(), 4200);
  }

  function base64ToBlob(b64, mime) {
    const bytes = atob(b64);
    const arr = new Uint8Array(bytes.length);
    for (let i = 0; i < bytes.length; i++) arr[i] = bytes.charCodeAt(i);
    return new Blob([arr], { type: mime });
  }
  function base64ToText(b64) {
    return decodeURIComponent(escape(atob(b64)));
  }
  function downloadBlob(blob, filename) {
    const url = URL.createObjectURL(blob);
    const a = el("a");
    a.href = url; a.download = filename;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 4000);
  }

  async function api(path, options = {}) {
    const res = await fetch(path, options);
    const contentType = res.headers.get("content-type") || "";
    if (!res.ok) {
      let message = `Request failed (${res.status})`;
      if (contentType.includes("application/json")) {
        try { message = (await res.json()).error || message; } catch (_) {}
      }
      throw new Error(message);
    }
    return res;
  }
  async function apiJson(path, options = {}) {
    return (await api(path, options)).json();
  }

  // ── icons (inline, single stroke, no fills) ──────────────────────────

  const ICON = {
    check: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M5 12l5 5L19 7"/></svg>',
    warn: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"><path d="M12 3l9 16H3z"/><path d="M12 10v4M12 17.5v.1"/></svg>',
    x: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"><path d="M6 6l12 12M18 6L6 18"/></svg>',
    dna: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round"><path d="M4 5c4 5 12 5 16 0M4 19c4-5 12-5 16 0"/><path d="M8 8l8 8M16 8l-8 8" opacity="0.5"/></svg>',
  };

  // ── navigation ────────────────────────────────────────────────────────

  function initNav() {
    qsa(".nav-item").forEach((btn) => {
      btn.addEventListener("click", () => {
        qsa(".nav-item").forEach((b) => b.classList.remove("active"));
        qsa(".view").forEach((v) => v.classList.remove("active"));
        btn.classList.add("active");
        qs(`#view-${btn.dataset.view}`).classList.add("active");
        if (btn.dataset.view === "history") loadHistory();
        if (btn.dataset.view === "templates") loadTemplates();
      });
    });
  }

  function goToView(name) {
    qs(`.nav-item[data-view="${name}"]`).click();
  }

  // ── dropzone helper ───────────────────────────────────────────────────

  function wireDropzone(zoneEl, inputEl, onFile) {
    zoneEl.addEventListener("click", () => inputEl.click());
    inputEl.addEventListener("change", () => {
      if (inputEl.files.length) onFile(inputEl.multiple ? Array.from(inputEl.files) : inputEl.files[0]);
    });
    ["dragenter", "dragover"].forEach((evt) =>
      zoneEl.addEventListener(evt, (e) => { e.preventDefault(); zoneEl.classList.add("drag-over"); }));
    ["dragleave", "drop"].forEach((evt) =>
      zoneEl.addEventListener(evt, (e) => { e.preventDefault(); zoneEl.classList.remove("drag-over"); }));
    zoneEl.addEventListener("drop", (e) => {
      const files = Array.from(e.dataTransfer.files || []);
      if (!files.length) return;
      onFile(inputEl.multiple ? files : files[0]);
    });
  }

  // ── 3Dmol wrapper (defensive: CDN may be unavailable) ────────────────

  function createViewer(containerEl) {
    if (typeof $3Dmol === "undefined") {
      containerEl.parentElement.innerHTML =
        '<div class="viewer-empty"><span>3D viewer script did not load</span></div>';
      return null;
    }
    return $3Dmol.createViewer(containerEl, { backgroundColor: "white" });
  }

  function renderStructure(viewer, pdbText, style = "cartoon") {
    if (!viewer) return;
    viewer.clear();
    viewer.addModel(pdbText, "pdb");
    if (style === "cartoon") {
      viewer.setStyle({}, { cartoon: { color: "spectrum" }, stick: { radius: 0.12, hidden: true } });
      viewer.setStyle({ hetflag: true }, { stick: { colorscheme: "grayCarbon", radius: 0.18 } });
    } else {
      viewer.setStyle({}, { stick: { radius: 0.15 } });
    }
    viewer.zoomTo();
    viewer.render();
  }

  // ── shared settings panel (used by Batch and High-throughput) ───────

  function buildSettingsPanel(container, idPrefix) {
    container.innerHTML = `
      <div class="stack" style="gap:var(--space-3)">
        <label class="field">
          <span class="field-label">pH — <span class="value-mono" id="${idPrefix}-ph-readout">7.4</span></span>
          <input type="range" id="${idPrefix}-ph" min="0" max="14" step="0.1" value="7.4">
        </label>
        <label class="field">
          <span class="field-label">Chains to keep <span class="optional">(comma-separated, blank = all)</span></span>
          <input type="text" id="${idPrefix}-chains" placeholder="A, B">
        </label>
        <label class="field">
          <span class="field-label">Heteroatoms to remove <span class="optional">(comma-separated, or ALL)</span></span>
          <input type="text" id="${idPrefix}-remove-hets" placeholder="SO4, PO4">
        </label>
        <label class="field">
          <span class="field-label">Ligands to protect <span class="optional">(survives ALL)</span></span>
          <input type="text" id="${idPrefix}-protect" placeholder="LIG">
        </label>
        <label class="check-row"><input type="checkbox" id="${idPrefix}-remove-water" checked><span class="check-text">Remove water</span></label>
        <label class="check-row"><input type="checkbox" id="${idPrefix}-structural-water"><span class="check-text">Keep structural waters</span></label>
        <label class="check-row"><input type="checkbox" id="${idPrefix}-loops"><span class="check-text">Reconstruct missing loops</span></label>
        <label class="check-row"><input type="checkbox" id="${idPrefix}-missing-atoms"><span class="check-text">Add missing heavy atoms</span></label>
        <label class="check-row"><input type="checkbox" id="${idPrefix}-minimize"><span class="check-text">Run energy minimisation</span></label>
        <div class="stack hidden" id="${idPrefix}-minimize-options" style="gap:var(--space-3);padding-left:23px">
          <label class="field">
            <span class="field-label">Force field</span>
            <select id="${idPrefix}-force-field">
              <option value="amber14">AMBER14</option>
              <option value="charmm36">CHARMM36</option>
            </select>
          </label>
          <label class="check-row"><input type="checkbox" id="${idPrefix}-gbsa" checked><span class="check-text">GBSA implicit solvent</span></label>
        </div>
        <label class="field">
          <span class="field-label">Docking target <span class="optional">(optional)</span></span>
          <select id="${idPrefix}-docking">
            <option value="">Standard PDB</option>
            <option value="vina">AutoDock Vina (PDBQT)</option>
            <option value="autodock">AutoDock (PDBQT)</option>
            <option value="gromacs">GROMACS</option>
          </select>
        </label>
      </div>`;

    const ph = qs(`#${idPrefix}-ph`, container);
    ph.addEventListener("input", () => { qs(`#${idPrefix}-ph-readout`, container).textContent = ph.value; });
    const minCheck = qs(`#${idPrefix}-minimize`, container);
    minCheck.addEventListener("change", () =>
      qs(`#${idPrefix}-minimize-options`, container).classList.toggle("hidden", !minCheck.checked));
  }

  function splitList(value) {
    return value.split(",").map((s) => s.trim()).filter(Boolean);
  }

  function readSettingsPanel(idPrefix) {
    const val = (id) => document.getElementById(`${idPrefix}-${id}`);
    return {
      ph: parseFloat(val("ph").value),
      chains: splitList(val("chains").value),
      remove_heteros: splitList(val("remove-hets").value).map((s) => s.toUpperCase()),
      protect_ligands: splitList(val("protect").value).map((s) => s.toUpperCase()),
      remove_water: val("remove-water").checked,
      keep_structural_waters: val("structural-water").checked,
      reconstruct_loops: val("loops").checked,
      add_missing_atoms: val("missing-atoms").checked,
      run_minimization: val("minimize").checked,
      force_field: val("force-field").value,
      use_gbsa: val("gbsa").checked,
      docking_target: val("docking").value || null,
    };
  }

  // ══════════════════════════════════════════════════════ PREPARE ═════

  const Prepare = (() => {
    let sessionId = null;
    let currentFile = null;
    let metadata = null;
    let originalPdbText = null;
    let preparedPdbText = null;
    let currentReport = null;
    const chainState = new Set();          // selected = keep only these; empty = all
    const hetState = new Map();            // name -> 'remove' | 'protect'
    let viewer = null;
    let viewerMode = "prepared";

    function reset() {
      sessionId = null; currentFile = null; metadata = null;
      originalPdbText = null; preparedPdbText = null; currentReport = null;
      chainState.clear(); hetState.clear();
      qs("#prepare-upload").classList.remove("hidden");
      qs("#prepare-workspace").classList.add("hidden");
      qs("#report-section").classList.add("hidden");
      qs("#prepare-file-input").value = "";
    }

    async function onFileChosen(file) {
      currentFile = file;
      qs("#prepare-file-chip").innerHTML =
        `${escapeHtml(file.name)} <button id="prepare-file-clear" title="Remove">${ICON.x}</button>`;
      qs("#prepare-upload").classList.add("hidden");
      qs("#prepare-workspace").classList.remove("hidden");
      qs("#report-section").classList.add("hidden");
      qs("#prepare-file-clear").addEventListener("click", reset);

      const summaryEl = qs("#analysis-summary");
      summaryEl.innerHTML = `<div class="row" style="color:var(--ink-faint);font-size:12.5px"><span class="spinner"></span>Analysing…</div>`;

      try {
        const form = new FormData();
        form.append("file", file);
        const data = await apiJson("/api/analyze", { method: "POST", body: form });
        sessionId = data.session_id;
        metadata = data.metadata;
        renderAnalysis();
        loadOriginalIntoViewer();
      } catch (err) {
        summaryEl.innerHTML = "";
        toast(err.message, "rust");
      }
    }

    function renderAnalysis() {
      // chains
      const chainPills = qs("#chain-pills");
      chainPills.innerHTML = "";
      metadata.chains.forEach((c) => {
        const pill = el("button", "pill", escapeHtml(c));
        pill.addEventListener("click", () => {
          if (chainState.has(c)) chainState.delete(c); else chainState.add(c);
          pill.classList.toggle("selected", chainState.has(c));
        });
        chainPills.appendChild(pill);
      });
      if (!metadata.chains.length) chainPills.innerHTML = '<span class="dropzone-hint">None detected</span>';

      // heteroatoms
      const hetPills = qs("#hetero-pills");
      hetPills.innerHTML = "";
      if (!metadata.heteroatoms.length) {
        hetPills.innerHTML = '<span class="dropzone-hint">None detected</span>';
      }
      metadata.heteroatoms.forEach((name) => {
        const pill = el("button", "pill", escapeHtml(name));
        const cycle = () => {
          const current = hetState.get(name);
          if (!current) hetState.set(name, "remove");
          else if (current === "remove") hetState.set(name, "protect");
          else hetState.delete(name);
          paintHetPill(pill, name);
        };
        pill.addEventListener("click", cycle);
        hetPills.appendChild(pill);
        paintHetPill(pill, name);
      });

      // structural water only matters if water is present
      qs("#opt-structural-water-row").classList.toggle("hidden", !metadata.water_count);

      // summary
      const summary = qs("#analysis-summary");
      summary.innerHTML = "";
      const rows = [
        ["Chains", metadata.chains.join(", ") || "—"],
        ["Total atoms", fmtNum(metadata.atoms_total, 0)],
        ["Waters", fmtNum(metadata.water_count, 0)],
        ["Heteroatoms", metadata.heteroatoms.length ? metadata.heteroatoms.join(", ") : "None"],
        ["Models", fmtNum(metadata.model_count || 1, 0)],
      ];
      rows.forEach(([k, v]) => {
        const row = el("div", "kv-row");
        row.innerHTML = `<span class="kv-key">${k}</span><span class="kv-val">${escapeHtml(String(v))}</span>`;
        summary.appendChild(row);
      });

      const gapEl = qs("#gap-callout");
      gapEl.innerHTML = "";
      if (metadata.sequence_gaps && metadata.sequence_gaps.length) {
        const cal = el("div", "callout callout-amber");
        cal.innerHTML = `${ICON.warn}<span class="callout-body">${metadata.sequence_gaps.length} sequence gap(s) detected. Enable "Reconstruct missing loops" to rebuild them.</span>`;
        gapEl.appendChild(cal);
      }
    }

    function paintHetPill(pillEl, name) {
      const state = hetState.get(name);
      pillEl.classList.toggle("selected", !!state);
      pillEl.innerHTML = escapeHtml(name) +
        (state === "remove" ? ' <span class="count">remove</span>' :
         state === "protect" ? ' <span class="count">protect</span>' : "");
    }

    async function loadOriginalIntoViewer() {
      try {
        const res = await api(`/api/history/pdb/${sessionId}`);
        originalPdbText = await res.text();
        if (!viewer) viewer = createViewer(qs("#prepare-viewer"));
        if (viewer) {
          renderStructure(viewer, originalPdbText);
          qs("#prepare-viewer-toolbar").classList.remove("hidden");
          qs("#viewer-label").textContent = "Original";
          viewerMode = "original";
        }
      } catch (_) { /* viewer is a convenience, not required */ }
    }

    function gatherSettings() {
      const removeHets = [];
      const protectLigands = [];
      hetState.forEach((state, name) => {
        if (state === "remove") removeHets.push(name); else protectLigands.push(name);
      });
      // Removing every named heteroatom individually has the same effect as
      // the surgical default (unlisted heteroatoms are kept); ALL is only
      // needed when the user wants to remove something not explicitly listed.
      return {
        chains: chainState.size ? Array.from(chainState) : null,
        remove_water: qs("#opt-remove-water").checked,
        remove_heteros: removeHets,
        protect_ligands: protectLigands,
        keep_structural_waters: qs("#opt-structural-water").checked,
        reconstruct_loops: qs("#opt-reconstruct-loops").checked,
        add_missing_atoms: qs("#opt-missing-atoms").checked,
        run_minimization: qs("#opt-minimize").checked,
        use_gbsa: qs("#opt-gbsa").checked,
        force_field: qs("#opt-force-field").value,
        docking_target: qs("#opt-docking-target").value || null,
        ph: parseFloat(qs("#ph-slider").value),
      };
    }

    async function process() {
      if (!currentFile) return;
      const btn = qs("#process-btn");
      const label = qs("#process-btn-label");
      btn.disabled = true;
      const original = label.textContent;
      label.innerHTML = '<span class="spinner"></span> Processing…';

      const settings = gatherSettings();
      const form = new FormData();
      form.append("file", currentFile);
      form.append("ph", settings.ph);
      form.append("remove_water", settings.remove_water);
      form.append("keep_structural_waters", settings.keep_structural_waters);
      form.append("reconstruct_loops", settings.reconstruct_loops);
      form.append("add_missing_atoms", settings.add_missing_atoms);
      form.append("run_minimization", settings.run_minimization);
      form.append("use_gbsa", settings.use_gbsa);
      form.append("force_field", settings.force_field);
      if (settings.docking_target) form.append("docking_target", settings.docking_target);
      if (settings.chains) form.append("chains", JSON.stringify(settings.chains));
      form.append("remove_heteros", JSON.stringify(settings.remove_heteros));
      form.append("protect_ligands", JSON.stringify(settings.protect_ligands));

      try {
        const data = await apiJson("/api/process", { method: "POST", body: form });
        currentReport = data.report;
        preparedPdbText = base64ToText(data.viewer_pdb_b64);
        // Binding-site analysis and the viewer's "prepared" mode must read the
        // structure /api/process actually produced (cleaned, protonated,
        // minimised), not the original upload's session — swap to the new id.
        sessionId = data.session_id;
        renderReport(data);
        toast("Structure processed", "moss");
      } catch (err) {
        toast(err.message, "rust");
      } finally {
        btn.disabled = false;
        label.textContent = original;
      }
    }

    function renderReport(data) {
      const report = data.report;
      qs("#report-section").classList.remove("hidden");

      const ac = report.atom_counts || {};
      const stats = qs("#report-stats");
      stats.innerHTML = "";
      const statDefs = [
        ["Atoms before", fmtNum(ac.before_processing, 0)],
        ["Atoms after", fmtNum(ac.after_processing, 0)],
        ["Delta", (ac.delta >= 0 ? "+" : "") + fmtNum(ac.delta, 0)],
        ["Waters removed", fmtNum(report.water_molecules_removed, 0)],
        ["Waters retained", fmtNum(report.water_molecules_retained, 0)],
        ["Time", fmtNum(report.processing_time_seconds, 2) + " s"],
      ];
      statDefs.forEach(([k, v]) => {
        const s = el("div", "stat");
        s.innerHTML = `<div class="stat-label">${k}</div><div class="stat-value">${v}</div>`;
        stats.appendChild(s);
      });

      // warnings
      const warnEl = qs("#report-warnings");
      warnEl.innerHTML = "";
      const warnings = data.warnings || report.warnings || [];
      if (warnings.length) {
        const list = el("ul", "warnings-list");
        warnings.forEach((w) => list.appendChild(el("li", "", escapeHtml(w))));
        warnEl.appendChild(list);
      }

      // protonation
      const prot = report.protonation || {};
      const protEl = qs("#report-protonation");
      protEl.innerHTML = "";
      const protRows = [
        ["Hydrogens added", prot.hydrogens_added ? "Yes" : "No"],
        ["pH", fmtNum(prot.ph ?? report.protonation_ph, 2)],
        ["Ligands protected", (prot.ligands_preserved || []).join(", ") || "None"],
        ["Nonstandard replaced", (prot.nonstandard_replaced || []).join(", ") || "None"],
        ["Loops rebuilt", fmtNum(prot.loops_reconstructed || 0, 0)],
        ["Terminals repaired", fmtNum(prot.terminals_repaired || 0, 0)],
      ];
      protRows.forEach(([k, v]) => {
        const row = el("div", "kv-row");
        row.innerHTML = `<span class="kv-key">${k}</span><span class="kv-val">${escapeHtml(String(v))}</span>`;
        protEl.appendChild(row);
      });

      // minimisation
      const mini = report.energy_minimization;
      const miniPanel = qs("#report-minimization-panel");
      if (mini) {
        miniPanel.classList.remove("hidden");
        const miniEl = qs("#report-minimization");
        miniEl.innerHTML = "";
        const statusBadge = {
          full: '<span class="badge badge-moss">full structure</span>',
          partial: '<span class="badge badge-amber">partial</span>',
          partial_no_implicit_solvent: '<span class="badge badge-amber">partial, no solvent</span>',
          failed: '<span class="badge badge-rust">did not run</span>',
        }[mini.status] || escapeHtml(String(mini.status));
        const miniRows = [
          ["Status", statusBadge],
          ["Force field", escapeHtml(mini.force_field || "—")],
          ["Energy before", mini.energy_before_kJ_mol != null ? fmtNum(mini.energy_before_kJ_mol) + " kJ/mol" : "—"],
          ["Energy after", mini.energy_after_kJ_mol != null ? fmtNum(mini.energy_after_kJ_mol) + " kJ/mol" : "—"],
          ["Converged", mini.converged ? `${ICON.check} yes` : "no"],
          ["RMS residual force", mini.rms_force_kJ_mol_nm != null ? fmtNum(mini.rms_force_kJ_mol_nm, 2) + " kJ/mol/nm" : "—"],
        ];
        if (mini.excluded_residues && mini.excluded_residues.length) {
          miniRows.push(["Held at input coordinates", mini.excluded_residues.join(", ")]);
        }
        if (mini.restrained_atoms) {
          miniRows.push(["Pocket atoms restrained", fmtNum(mini.restrained_atoms, 0)]);
        }
        if (mini.error) miniRows.push(["Error", escapeHtml(mini.error)]);
        miniRows.forEach(([k, v]) => {
          const row = el("div", "kv-row");
          row.innerHTML = `<span class="kv-key">${k}</span><span class="kv-val">${v}</span>`;
          miniEl.appendChild(row);
        });
      } else {
        miniPanel.classList.add("hidden");
      }

      qs("#report-text-block").textContent = data.report_text || "";

      // viewer: switch to prepared
      if (viewer) {
        renderStructure(viewer, preparedPdbText);
        qs("#viewer-label").textContent = "Prepared";
        viewerMode = "prepared";
      }

      // downloads
      qs("#download-pdb-btn").onclick = () => {
        const isPdbqt = (data.filename || "").endsWith(".pdbqt");
        const blob = base64ToBlob(data.pdb_b64, isPdbqt ? "chemical/x-pdbqt" : "chemical/x-pdb");
        downloadBlob(blob, data.filename);
      };
      qs("#report-copy-btn").onclick = () => {
        navigator.clipboard.writeText(data.report_text || "").then(() => toast("Report copied"));
      };
      qs("#report-download-btn").onclick = async () => {
        const res = await api("/api/download-report", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify(report),
        });
        downloadBlob(await res.blob(), "preparation_report.txt");
      };
      qs("#prepare-sites-btn").onclick = () => {
        goToView("sites");
        Sites.loadFromSession(sessionId, currentFile.name);
      };
    }

    function init() {
      wireDropzone(qs("#prepare-dropzone"), qs("#prepare-file-input"), onFileChosen);
      qs("#prepare-reset-btn").addEventListener("click", reset);
      qs("#ph-slider").addEventListener("input", (e) => {
        qs("#ph-readout").textContent = e.target.value;
      });
      qs("#opt-minimize").addEventListener("change", (e) => {
        qs("#minimize-options").classList.toggle("hidden", !e.target.checked);
      });
      qs("#process-btn").addEventListener("click", process);
      qs("#viewer-toggle-btn").addEventListener("click", () => {
        if (!viewer) return;
        if (viewerMode === "original" && preparedPdbText) {
          renderStructure(viewer, preparedPdbText);
          qs("#viewer-label").textContent = "Prepared";
          qs("#viewer-toggle-btn").textContent = "Show original";
          viewerMode = "prepared";
        } else if (originalPdbText) {
          renderStructure(viewer, originalPdbText);
          qs("#viewer-label").textContent = "Original";
          qs("#viewer-toggle-btn").textContent = "Show prepared";
          viewerMode = "original";
        }
      });
    }

    return { init, applySettings: applyTemplateSettings, gatherSettingsForTemplate: () => gatherSettings() };

    function applyTemplateSettings(settings) {
      if (!settings) return;
      if (settings.ph != null) { qs("#ph-slider").value = settings.ph; qs("#ph-readout").textContent = settings.ph; }
      if ("remove_water" in settings) qs("#opt-remove-water").checked = !!settings.remove_water;
      if ("keep_structural_waters" in settings) qs("#opt-structural-water").checked = !!settings.keep_structural_waters;
      if ("reconstruct_loops" in settings) qs("#opt-reconstruct-loops").checked = !!settings.reconstruct_loops;
      if ("add_missing_atoms" in settings) qs("#opt-missing-atoms").checked = !!settings.add_missing_atoms;
      if ("run_minimization" in settings) {
        qs("#opt-minimize").checked = !!settings.run_minimization;
        qs("#minimize-options").classList.toggle("hidden", !settings.run_minimization);
      }
      if (settings.force_field) qs("#opt-force-field").value = settings.force_field;
      if ("use_gbsa" in settings) qs("#opt-gbsa").checked = !!settings.use_gbsa;
      if ("docking_target" in settings) qs("#opt-docking-target").value = settings.docking_target || "";
      toast("Template applied to Prepare panel");
    }
  })();

  // ══════════════════════════════════════════════════════ SITES ═══════

  const Sites = (() => {
    let viewer = null;

    function reset() {
      qs("#sites-dropzone").classList.remove("hidden");
      qs("#sites-workspace").classList.add("hidden");
      qs("#sites-file-input").value = "";
    }

    async function run(payload, filenameForChip) {
      qs("#sites-dropzone").classList.add("hidden");
      qs("#sites-workspace").classList.remove("hidden");
      qs("#sites-file-chip").innerHTML =
        `${escapeHtml(filenameForChip)} <button id="sites-file-clear" title="Remove">${ICON.x}</button>`;
      qs("#sites-file-clear").addEventListener("click", reset);

      const listEl = qs("#sites-list");
      listEl.innerHTML = `<div class="row" style="color:var(--ink-faint);font-size:12.5px"><span class="spinner"></span>Scanning for pockets…</div>`;

      try {
        const data = await apiJson(payload.url, payload.init);
        renderSites(data);
      } catch (err) {
        listEl.innerHTML = "";
        toast(err.message, "rust");
      }
    }

    function loadFromSession(sessionId, filename) {
      run({ url: `/api/analyze-site/${sessionId}`, init: { method: "POST" } }, filename);
      loadStructureForViewer(sessionId);
    }

    async function loadStructureForViewer(sessionId) {
      try {
        const text = await (await api(`/api/history/pdb/${sessionId}`)).text();
        if (!viewer) viewer = createViewer(qs("#sites-viewer"));
        if (viewer) renderStructure(viewer, text, "cartoon");
      } catch (_) {}
    }

    function onFileChosen(file) {
      const form = new FormData();
      form.append("file", file);
      run({ url: "/api/analyze-site", init: { method: "POST", body: form } }, file.name);
      file.text().then((text) => {
        if (!viewer) viewer = createViewer(qs("#sites-viewer"));
        if (viewer) renderStructure(viewer, text, "cartoon");
      });
    }

    function renderSites(data) {
      const listEl = qs("#sites-list");
      listEl.innerHTML = "";

      const summary = el("div", "kv-row");
      summary.style.marginBottom = "8px";
      summary.innerHTML = `<span class="kv-key">${escapeHtml(data.message || "")}</span>`;
      listEl.appendChild(summary);

      if (!data.sites || !data.sites.length) {
        const empty = el("div", "empty-state");
        empty.innerHTML = `<div class="empty-state-title">No pockets found</div><div class="empty-state-desc">Try a structure with a defined binding cavity, or check that hydrogens/waters were removed.</div>`;
        listEl.appendChild(empty);
        return;
      }

      const rankColors = ["#2f6e64", "#5b8fae", "#9a6b1f", "#a1432f", "#6b6d64"];
      if (viewer) {
        qsa(".sphere-shape-marker"); // no-op, placeholder for clarity
        data.sites.forEach((site, i) => {
          viewer.addSphere({
            center: { x: site.centroid[0], y: site.centroid[1], z: site.centroid[2] },
            radius: 1.6,
            color: rankColors[i % rankColors.length],
            opacity: 0.85,
          });
        });
        viewer.render();
      }

      data.sites.forEach((site, i) => {
        const card = el("div", "panel");
        card.style.padding = "var(--space-4)";
        const factors = site.drugability_factors || {};
        const factorRows = Object.entries(factors)
          .map(([name, f]) => `<div class="kv-row"><span class="kv-key">${escapeHtml(name.replace(/_/g, " "))}</span><span class="kv-val">${fmtNum(f.score, 2)} × ${f.weight}</span></div>`)
          .join("");
        card.innerHTML = `
          <div class="row-between" style="margin-bottom:8px">
            <span class="row" style="gap:8px">
              <span style="width:9px;height:9px;border-radius:50%;background:${rankColors[i % rankColors.length]};display:inline-block"></span>
              <strong style="font-size:13px">Pocket ${site.id}</strong>
            </span>
            <span class="badge badge-accent">druggability ${fmtNum(site.drugability_score, 2)}</span>
          </div>
          <div class="kv-row"><span class="kv-key">Volume</span><span class="kv-val">${fmtNum(site.volume, 0)} Å³</span></div>
          <div class="kv-row"><span class="kv-key">Concavity</span><span class="kv-val">${fmtNum(site.concavity, 2)}</span></div>
          <div class="kv-row"><span class="kv-key">Pharmacophore points</span><span class="kv-val">${(site.pharmacophore_points || []).length}</span></div>
          <details style="margin-top:6px">
            <summary style="cursor:pointer;font-size:11.5px;color:var(--ink-faint)">Scoring breakdown</summary>
            <div class="kv-list" style="margin-top:4px">${factorRows}</div>
          </details>
          <div class="eyebrow" style="margin-top:10px">Lining residues</div>
          <p style="font-size:12px;color:var(--ink-soft);line-height:1.6">${(site.residues || []).join(", ") || "—"}</p>
        `;
        listEl.appendChild(card);
      });
    }

    function init() {
      wireDropzone(qs("#sites-dropzone"), qs("#sites-file-input"), onFileChosen);
      qs("#sites-reset-btn").addEventListener("click", reset);
    }

    return { init, loadFromSession };
  })();

  // ══════════════════════════════════════════════════════ BATCH ═══════

  const Batch = (() => {
    let files = [];

    function renderFileList() {
      const listEl = qs("#batch-file-list");
      listEl.innerHTML = "";
      files.forEach((f, i) => {
        const row = el("div", "file-chip");
        row.style.display = "inline-flex";
        row.innerHTML = `${escapeHtml(f.name)} <button data-i="${i}" title="Remove">${ICON.x}</button>`;
        row.querySelector("button").addEventListener("click", () => {
          files.splice(i, 1); renderFileList();
        });
        listEl.appendChild(row);
      });
      qs("#batch-run-btn").disabled = files.length === 0;
    }

    function onFilesChosen(chosen) {
      files = files.concat(chosen);
      renderFileList();
    }

    async function run() {
      const btn = qs("#batch-run-btn");
      btn.disabled = true;
      const original = btn.textContent;
      btn.innerHTML = '<span class="spinner"></span> Running…';

      const settings = readSettingsPanel("batch");
      const config = { ...settings, ph: settings.ph };
      const form = new FormData();
      files.forEach((f) => form.append("files", f));
      form.append("config", JSON.stringify(config));

      try {
        const res = await api("/api/batch", { method: "POST", body: form });
        const blob = await res.blob();
        downloadBlob(blob, "bioprep_batch_results.zip");
        toast(`Batch complete — ${files.length} file(s) submitted`, "moss");
      } catch (err) {
        toast(err.message, "rust");
      } finally {
        btn.disabled = false;
        btn.textContent = original;
      }
    }

    function init() {
      buildSettingsPanel(qs("#batch-settings"), "batch");
      wireDropzone(qs("#batch-dropzone"), qs("#batch-file-input"), onFilesChosen);
      qs("#batch-run-btn").addEventListener("click", run);
    }

    return { init };
  })();

  // ══════════════════════════════════════════════════ HIGH-THROUGHPUT ═

  const Throughput = (() => {
    let zipFile = null;

    function onFileChosen(file) {
      zipFile = file;
      qs("#ht-file-chip").classList.remove("hidden");
      qs("#ht-file-chip").innerHTML =
        `${escapeHtml(file.name)} <button id="ht-file-clear" title="Remove">${ICON.x}</button>`;
      qs("#ht-file-clear").addEventListener("click", () => {
        zipFile = null;
        qs("#ht-file-chip").classList.add("hidden");
        qs("#ht-file-input").value = "";
        qs("#ht-run-btn").disabled = true;
      });
      qs("#ht-run-btn").disabled = false;
    }

    async function run() {
      if (!zipFile) return;
      const btn = qs("#ht-run-btn");
      btn.disabled = true;
      const original = btn.textContent;
      btn.innerHTML = '<span class="spinner"></span> Running…';

      const settings = readSettingsPanel("ht");
      const form = new FormData();
      form.append("file", zipFile);
      form.append("config", JSON.stringify(settings));

      try {
        const res = await api("/api/high-throughput", { method: "POST", body: form });
        const blob = await res.blob();
        downloadBlob(blob, "bioprep_ht_results.zip");
        toast("High-throughput run complete", "moss");
      } catch (err) {
        toast(err.message, "rust");
      } finally {
        btn.disabled = false;
        btn.textContent = original;
      }
    }

    function init() {
      buildSettingsPanel(qs("#ht-settings"), "ht");
      wireDropzone(qs("#ht-dropzone"), qs("#ht-file-input"), onFileChosen);
      qs("#ht-run-btn").addEventListener("click", run);
    }

    return { init };
  })();

  // ══════════════════════════════════════════════════════ HISTORY ═════

  async function loadHistory() {
    const root = qs("#history-content");
    root.innerHTML = `<div class="row" style="color:var(--ink-faint);font-size:12.5px"><span class="spinner"></span>Loading…</div>`;
    try {
      const jobs = await apiJson("/api/history");
      if (!jobs.length) {
        root.innerHTML = `<div class="empty-state"><div class="empty-state-title">No jobs yet</div><div class="empty-state-desc">Processed structures appear here.</div></div>`;
        return;
      }
      const table = el("table", "data-table");
      table.innerHTML = `<thead><tr>
        <th>File</th><th>Time</th><th>Format</th><th>Energy Δ (kJ/mol)</th><th>Status</th><th></th>
      </tr></thead>`;
      const tbody = el("tbody");
      jobs.forEach((job) => {
        const tr = el("tr");
        const hasEnergy = typeof job.energy_before === "number" && typeof job.energy_after === "number";
        const delta = hasEnergy ? fmtNum(job.energy_after - job.energy_before) : "—";
        tr.innerHTML = `
          <td>${escapeHtml(job.filename)}</td>
          <td class="mono">${escapeHtml(job.timestamp)}</td>
          <td>${escapeHtml(job.format || "PDB")}</td>
          <td class="mono">${delta}</td>
          <td><span class="badge ${job.status === "success" ? "badge-moss" : "badge-rust"}">${escapeHtml(job.status)}</span></td>
          <td><button class="btn btn-ghost btn-sm" data-id="${job.id}">View report</button></td>`;
        tr.querySelector("button").addEventListener("click", () => showHistoryModal(job));
        tbody.appendChild(tr);
      });
      table.appendChild(tbody);
      root.innerHTML = "";
      root.appendChild(table);
    } catch (err) {
      root.innerHTML = "";
      toast(err.message, "rust");
    }
  }

  function showHistoryModal(job) {
    const root = qs("#modal-root");
    const backdrop = el("div", "modal-backdrop");
    const reportText = job.report ? undefined : null;
    backdrop.innerHTML = `
      <div class="modal">
        <div class="modal-header">
          <strong style="font-size:13.5px">${escapeHtml(job.filename)}</strong>
          <button class="btn btn-ghost btn-sm" id="modal-close">${ICON.x}</button>
        </div>
        <div class="modal-body">
          <div class="viewer-shell" style="margin-bottom:var(--space-4)">
            <div class="viewer-canvas" id="modal-viewer" style="height:280px"></div>
          </div>
          <div class="log-block">${escapeHtml(JSON.stringify(job.report, null, 2))}</div>
        </div>
      </div>`;
    root.appendChild(backdrop);
    backdrop.addEventListener("click", (e) => { if (e.target === backdrop) backdrop.remove(); });
    qs("#modal-close", backdrop).addEventListener("click", () => backdrop.remove());

    api(`/api/history/pdb/${job.id}`).then((r) => r.text()).then((text) => {
      const v = createViewer(qs("#modal-viewer", backdrop));
      if (v) renderStructure(v, text);
    }).catch(() => {});
  }

  // ══════════════════════════════════════════════════════ TEMPLATES ═══

  async function loadTemplates() {
    const root = qs("#templates-content");
    root.innerHTML = `<div class="row" style="color:var(--ink-faint);font-size:12.5px"><span class="spinner"></span>Loading…</div>`;
    try {
      const templates = await apiJson("/api/templates");
      root.innerHTML = "";

      const savePanel = el("div", "panel");
      savePanel.style.marginBottom = "var(--space-5)";
      savePanel.innerHTML = `
        <div class="panel-title">Save current Prepare settings</div>
        <div class="row" style="gap:var(--space-2)">
          <input type="text" id="template-name-input" placeholder="Template name" style="flex:1">
          <button class="btn btn-primary" id="template-save-btn">Save</button>
        </div>`;
      root.appendChild(savePanel);
      qs("#template-save-btn", savePanel).addEventListener("click", async () => {
        const name = qs("#template-name-input", savePanel).value.trim();
        if (!name) { toast("Name the template first", "rust"); return; }
        try {
          await apiJson("/api/templates", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ name, settings: Prepare.gatherSettingsForTemplate() }),
          });
          toast("Template saved", "moss");
          loadTemplates();
        } catch (err) { toast(err.message, "rust"); }
      });

      const names = Object.keys(templates);
      if (!names.length) {
        root.appendChild(el("div", "empty-state",
          '<div class="empty-state-title">No saved templates</div><div class="empty-state-desc">Configure the Prepare panel, then save it here.</div>'));
        return;
      }

      const grid = el("div", "stack", "");
      names.forEach((name) => {
        const entry = templates[name];
        const row = el("div", "panel");
        row.style.padding = "var(--space-3) var(--space-4)";
        row.innerHTML = `
          <div class="row-between">
            <div>
              <div style="font-weight:600;font-size:13px">${escapeHtml(name)}</div>
              <div class="dropzone-hint">Saved ${escapeHtml(entry.created_at)}</div>
            </div>
            <div class="row" style="gap:6px">
              <button class="btn btn-sm" data-act="apply">Apply</button>
              <button class="btn btn-ghost btn-sm" data-act="delete">Delete</button>
            </div>
          </div>`;
        row.querySelector('[data-act="apply"]').addEventListener("click", () => {
          Prepare.applySettings(entry.settings);
          goToView("prepare");
        });
        row.querySelector('[data-act="delete"]').addEventListener("click", async () => {
          await apiJson(`/api/templates/${encodeURIComponent(name)}`, { method: "DELETE" });
          loadTemplates();
        });
        grid.appendChild(row);
      });
      root.appendChild(grid);
    } catch (err) {
      root.innerHTML = "";
      toast(err.message, "rust");
    }
  }

  // ══════════════════════════════════════════════════════ BOOT ════════

  document.addEventListener("DOMContentLoaded", () => {
    initNav();
    Prepare.init();
    Sites.init();
    Batch.init();
    Throughput.init();
  });
})();
