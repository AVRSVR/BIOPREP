"""
reporter.py – Generates structured Preparation Reports after processing runs.
"""
import json
from datetime import datetime, timezone



def build_report(
    filename,
    chains_detected,
    chains_retained,
    waters_removed,
    heteroatoms_removed,
    heteroatoms_retained,
    hydrogens_added,
    ph_used,
    missing_residues,
    atoms_before,
    atoms_after,
    minimization_stats=None,
    docking_target=None,
    processing_time_s=None,
):
    """
    Constructs and returns a dictionary representing the full Preparation Report.
    """
    report = {
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "input_file": filename,
        "chains": {
            "detected": chains_detected,
            "retained": chains_retained,
        },
        "water_molecules_removed": waters_removed,
        "heteroatoms": {
            "removed": heteroatoms_removed,
            "retained": heteroatoms_retained,
        },
        "hydrogens_added": hydrogens_added,
        "protonation_ph": ph_used,
        "missing_residues_detected": missing_residues if missing_residues else [],
        "atom_counts": {
            "before_processing": atoms_before,
            "after_processing": atoms_after,
            "delta": atoms_after - atoms_before,
        },
    }

    if minimization_stats:
        report["energy_minimization"] = minimization_stats

    if docking_target:
        report["docking_target"] = docking_target

    if processing_time_s is not None:
        report["processing_time_seconds"] = round(processing_time_s, 2)

    return report


def report_to_text(report):
    """
    Renders a PreparationReport dict as a human-readable plain-text string.
    """
    lines = [
        "=" * 60,
        "          BIOPREP – PROTEIN PREPARATION REPORT",
        "=" * 60,
        f"  Generated : {report.get('generated_at', 'N/A')}",
        f"  Input File: {report.get('input_file', 'N/A')}",
        "",
        "── CHAINS ─────────────────────────────────────────────────",
        f"  Detected  : {', '.join(report['chains']['detected']) or 'None'}",
        f"  Retained  : {', '.join(report['chains']['retained']) or 'None'}",
        "",
        "── WATER MOLECULES ─────────────────────────────────────────",
        f"  Removed   : {report.get('water_molecules_removed', 0)}",
        "",
        "── HETEROATOMS / LIGANDS ───────────────────────────────────",
        f"  Removed   : {', '.join(report['heteroatoms']['removed']) or 'None'}",
        f"  Retained  : {', '.join(report['heteroatoms']['retained']) or 'None'}",
        "",
        "── PROTONATION ─────────────────────────────────────────────",
        f"  Hydrogens Added : {report.get('hydrogens_added', 'N/A')}",
        f"  Target pH       : {report.get('protonation_ph', 'N/A')}",
        "",
        "── STRUCTURE ───────────────────────────────────────────────",
        f"  Missing Residues: {len(report.get('missing_residues_detected', []))}",
        f"  Atoms Before    : {report['atom_counts']['before_processing']}",
        f"  Atoms After     : {report['atom_counts']['after_processing']}",
        f"  Delta (change)  : {report['atom_counts']['delta']:+d}",
    ]

    if "energy_minimization" in report:
        em = report["energy_minimization"]
        lines += [
            "",
            "── ENERGY MINIMIZATION ──────────────────────────────────",
            f"  Force Field     : {em.get('force_field', 'N/A')}",
            f"  Energy Before   : {em.get('energy_before_kJ_mol', 'N/A')} kJ/mol",
            f"  Energy After    : {em.get('energy_after_kJ_mol', 'N/A')} kJ/mol",
            f"  Iterations      : {em.get('iterations', 'N/A')}",
            f"  Converged       : {em.get('converged', 'N/A')}",
        ]

    if "docking_target" in report:
        lines += [
            "",
            "── DOCKING EXPORT ───────────────────────────────────────",
            f"  Target Software : {report['docking_target'].upper()}",
        ]

    if "processing_time_seconds" in report:
        lines += [
            "",
            f"  Processing Time : {report['processing_time_seconds']} s",
        ]

    lines += ["", "=" * 60]
    return "\n".join(lines)


def count_atoms_in_pdb(pdb_path):
    """
    Quick count of ATOM + HETATM lines in a raw PDB file without Biopython overhead.
    """
    count = 0
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    count += 1
    except Exception:
        pass
    return count
