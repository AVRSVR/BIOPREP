"""
reporter.py - Builds the structured Preparation Report and renders it as text.

The report is an audit trail, so it records what the pipeline actually did.
Earlier versions hardcoded ``hydrogens_added: True`` and claimed a docking
export even when OpenBabel had failed.
"""

from datetime import datetime, timezone


def build_report(
    filename,
    chains_detected,
    chains_retained,
    waters_removed,
    heteroatoms_removed,
    heteroatoms_retained,
    ph_used,
    missing_residues,
    atoms_before,
    atoms_after,
    protonation=None,
    waters_retained=0,
    minimization_stats=None,
    docking_export=None,
    settings=None,
    processing_time_s=None,
    warnings=None,
):
    """Construct the full Preparation Report as a plain dict."""
    protonation = protonation or {}

    report = {
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "input_file": filename,
        "chains": {
            "detected": list(chains_detected or []),
            "retained": list(chains_retained or []),
        },
        "water_molecules_removed": waters_removed,
        "water_molecules_retained": waters_retained,
        "heteroatoms": {
            "removed": list(heteroatoms_removed or []),
            "retained": list(heteroatoms_retained or []),
        },
        "protonation": {
            # Reported by the protonator rather than assumed.
            "hydrogens_added": protonation.get("hydrogens_added", False),
            "ph": protonation.get("ph", ph_used),
            "ligands_preserved": protonation.get("ligands_preserved", []),
            "nonstandard_replaced": protonation.get("nonstandard_replaced", []),
            "loops_reconstructed": protonation.get("loops_reconstructed", 0),
        },
        # Kept at the top level for backwards compatibility with older readers.
        "hydrogens_added": protonation.get("hydrogens_added", False),
        "protonation_ph": ph_used,
        "missing_residues_detected": missing_residues or [],
        "atom_counts": {
            "before_processing": atoms_before,
            "after_processing": atoms_after,
            "delta": (atoms_after or 0) - (atoms_before or 0),
        },
        "warnings": list(warnings or []),
    }

    if minimization_stats:
        report["energy_minimization"] = minimization_stats

    if docking_export:
        report["docking_export"] = docking_export
        if docking_export.get("succeeded"):
            report["docking_target"] = docking_export.get("target")

    if settings:
        report["settings_used"] = settings

    if processing_time_s is not None:
        report["processing_time_seconds"] = round(processing_time_s, 2)

    return report


def _join(values, empty="None"):
    values = [str(v) for v in (values or [])]
    return ", ".join(values) if values else empty


def report_to_text(report):
    """
    Render a report dict as human-readable text.

    Every lookup is defensive: this also renders reports posted back by a
    client, which may be partial or hand-edited.
    """
    report = report or {}
    chains = report.get("chains") or {}
    hets = report.get("heteroatoms") or {}
    atoms = report.get("atom_counts") or {}
    protonation = report.get("protonation") or {}

    before = atoms.get("before_processing")
    after = atoms.get("after_processing")
    delta = atoms.get("delta")
    delta_text = f"{delta:+d}" if isinstance(delta, int) else "N/A"

    lines = [
        "=" * 60,
        "          BIOPREP - PROTEIN PREPARATION REPORT",
        "=" * 60,
        f"  Generated : {report.get('generated_at', 'N/A')}",
        f"  Input File: {report.get('input_file', 'N/A')}",
        "",
        "-- CHAINS ---------------------------------------------------",
        f"  Detected  : {_join(chains.get('detected'))}",
        f"  Retained  : {_join(chains.get('retained'))}",
        "",
        "-- WATER MOLECULES ------------------------------------------",
        f"  Removed   : {report.get('water_molecules_removed', 0)}",
        f"  Retained  : {report.get('water_molecules_retained', 0)}",
        "",
        "-- HETEROATOMS / LIGANDS ------------------------------------",
        f"  Removed   : {_join(hets.get('removed'))}",
        f"  Retained  : {_join(hets.get('retained'))}",
        "",
        "-- PROTONATION ----------------------------------------------",
        f"  Hydrogens Added   : {protonation.get('hydrogens_added', report.get('hydrogens_added', 'N/A'))}",
        f"  Target pH         : {report.get('protonation_ph', 'N/A')}",
        f"  Ligands Protected : {_join(protonation.get('ligands_preserved'))}",
        f"  Nonstandard Fixed : {_join(protonation.get('nonstandard_replaced'))}",
        f"  Loops Rebuilt     : {protonation.get('loops_reconstructed', 0)}",
        "",
        "-- STRUCTURE ------------------------------------------------",
        f"  Missing Residues: {len(report.get('missing_residues_detected') or [])}",
        f"  Atoms Before    : {before if before is not None else 'N/A'}",
        f"  Atoms After     : {after if after is not None else 'N/A'}",
        f"  Delta (change)  : {delta_text}",
    ]

    minimization = report.get("energy_minimization")
    if minimization:
        lines += [
            "",
            "-- ENERGY MINIMIZATION --------------------------------------",
            f"  Status          : {minimization.get('status', 'N/A')}",
            f"  Force Field     : {minimization.get('force_field', 'N/A')}",
            f"  Energy Before   : {minimization.get('energy_before_kJ_mol', 'not run')} kJ/mol",
            f"  Energy After    : {minimization.get('energy_after_kJ_mol', 'not run')} kJ/mol",
            f"  Energy Decreased: {minimization.get('energy_decreased', 'N/A')}",
        ]
        if minimization.get("excluded_residues"):
            lines.append(
                f"  Held Fixed      : {_join(minimization['excluded_residues'])}"
            )
        if minimization.get("error"):
            lines.append(f"  ERROR           : {minimization['error']}")

    export = report.get("docking_export")
    if export:
        lines += [
            "",
            "-- DOCKING EXPORT -------------------------------------------",
            f"  Target Software : {str(export.get('target', 'N/A')).upper()}",
            f"  Succeeded       : {export.get('succeeded', False)}",
        ]
        if export.get("error"):
            lines.append(f"  ERROR           : {export['error']}")

    if report.get("warnings"):
        lines += ["", "-- WARNINGS -------------------------------------------------"]
        lines += [f"  ! {w}" for w in report["warnings"]]

    if report.get("processing_time_seconds") is not None:
        lines += ["", f"  Processing Time : {report['processing_time_seconds']} s"]

    lines += ["", "=" * 60]
    return "\n".join(lines)


def count_atoms_in_pdb(pdb_path):
    """Count ATOM + HETATM lines without Biopython overhead."""
    count = 0
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    count += 1
    except OSError:
        pass
    return count
