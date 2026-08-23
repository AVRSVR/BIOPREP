"""Independent verification of catalogue features 70-80 (exporter, reporter)."""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)

import inspect
import os
import shutil
import subprocess
import sys
import tempfile
import time
import warnings

warnings.filterwarnings("ignore")
BIOPREP = _os.path.join(_ROOT, "bioprep")
sys.path.insert(0, BIOPREP)
os.chdir(BIOPREP)

from bioprep.core import io as bio, cleaner, protonator      # noqa: E402
from bioprep.core import exporter, reporter                  # noqa: E402

TMP = tempfile.mkdtemp(prefix="exp_")
P = lambda n: os.path.join(TMP, n)      # noqa: E731
RESULTS = []


def check(feature, description, passed, detail=""):
    RESULTS.append((feature, description, passed, detail))


PROT = [l for l in open("1crn.pdb") if l.startswith("ATOM")]
raw = P("raw.pdb")
with open(raw, "w") as fh:
    fh.write("".join(PROT) + "END\n")
st = bio.load_pdb(raw)
bio.save_pdb(st, P("clean.pdb"), select=cleaner.clean_structure(st))
protonator.add_hydrogens(P("clean.pdb"), P("prot.pdb"))
PREPARED = P("prot.pdb")

HAVE_OBABEL = shutil.which("obabel") is not None
exp_src = inspect.getsource(exporter)
rep_src = inspect.getsource(reporter)


# ---- F70: PDBQT export -----------------------------------------------------
check(70, "OpenBabel available in this environment", HAVE_OBABEL, "obabel not on PATH")
check(70, "Gasteiger partial charges requested",
      "'--partialcharge', 'gasteiger'" in exp_src, "")
# DEVIATION FROM THE CATALOGUE, on purpose. Feature 70 says -xh merges
# non-polar hydrogens. In OpenBabel it means "preserve hydrogens": with it
# every hydrogen survived and was typed HD, the AutoDock type for a hydrogen
# on N or O, so carbon-bound hydrogens were presented to the scoring
# function as hydrogen-bond donors - 842 of them on streptavidin against 213
# genuinely polar. The default merges them, which is the conventional receptor.
check(70, "-xh is NOT passed (it preserves hydrogens, not merges them)",
      "'-xh'" not in exp_src, "")
check(70, "rigid receptor flag (-xr)", "'-xr'" in exp_src, "")

if HAVE_OBABEL:
    for target in ("autodock", "vina"):
        ok, result = exporter.export_structure(
            PREPARED, P(f"{target}.pdb"), target)
        check(70, f"{target} produces a .pdbqt",
              ok and str(result).endswith(".pdbqt"), str(result)[:70])
        if ok:
            text = open(result).read()
            check(70, f"{target} output has atom records",
                  any(l.startswith(("ATOM", "HETATM")) for l in text.splitlines()), "")
            check(70, f"{target} receptor carries no torsion tree",
                  "ROOT" not in text and "BRANCH" not in text, "")
            check(70, f"{target} atoms carry partial charges",
                  any(len(l) > 70 and l[:6] in ("ATOM  ", "HETATM")
                      for l in text.splitlines()), "")

# ---- F71: GROMACS export ---------------------------------------------------
ok, result = exporter.export_structure(PREPARED, P("g.pdb"), "gromacs")
check(71, "gromacs export succeeds", ok, str(result)[:60])
check(71, "output named with the _gromacs suffix",
      str(result).endswith("_gromacs.pdb"), os.path.basename(str(result)))
check(71, "gromacs output is a real PDB",
      os.path.isfile(result)
      and any(l.startswith("ATOM") for l in open(result)), "")
check(71, "gromacs output has the same atom count as the input",
      sum(1 for l in open(result) if l[:6] in ("ATOM  ", "HETATM"))
      == sum(1 for l in open(PREPARED) if l[:6] in ("ATOM  ", "HETATM")), "")

# splitext, not str.replace: a directory containing '.pdb' must survive
nested = os.path.join(TMP, "my.pdb.data")
os.makedirs(nested, exist_ok=True)
ok, result = exporter.export_structure(PREPARED, os.path.join(nested, "o.pdb"), "gromacs")
check(71, "a directory containing '.pdb' is not mangled",
      ok and "my.pdb.data" in str(result) and os.path.isfile(result), str(result)[:70])

# ---- F72: plain PDB export -------------------------------------------------
ok, result = exporter.export_structure(PREPARED, P("plain.pdb"), "pdb")
check(72, "default export writes a PDB", ok and str(result).endswith(".pdb"), "")
check(72, "plain export preserves every atom",
      sum(1 for l in open(result) if l[:6] in ("ATOM  ", "HETATM"))
      == sum(1 for l in open(PREPARED) if l[:6] in ("ATOM  ", "HETATM")), "")
ok, result = exporter.export_structure(PREPARED, P("unknown.pdb"), "something-else")
check(72, "an unrecognised format falls back to PDB",
      ok and str(result).endswith(".pdb"), str(result)[:60])
ok, result = exporter.export_structure(PREPARED, PREPARED, "pdb")
check(72, "export onto its own path does not truncate the file",
      ok and sum(1 for l in open(PREPARED) if l[:6] in ("ATOM  ", "HETATM")) > 0, "")

# ---- F73: OpenBabel error handling -----------------------------------------
real_run = subprocess.run

exporter.subprocess.run = lambda *a, **k: (_ for _ in ()).throw(FileNotFoundError("obabel"))
try:
    ok, message = exporter.convert_to_pdbqt("in.pdb", P("x.pdbqt"))
finally:
    exporter.subprocess.run = real_run
check(73, "missing OpenBabel reported, not raised",
      ok is False and "OpenBabel" in message, str(message)[:60])

exporter.subprocess.run = lambda *a, **k: (_ for _ in ()).throw(
    subprocess.TimeoutExpired("obabel", 120))
try:
    ok, message = exporter.convert_to_pdbqt("in.pdb", P("x.pdbqt"))
finally:
    exporter.subprocess.run = real_run
check(73, "timeout reported, not raised",
      ok is False and "timed out" in message.lower(), str(message)[:60])


class Failed:
    returncode = 1
    stderr = "obabel exploded"


exporter.subprocess.run = lambda *a, **k: Failed()
try:
    ok, message = exporter.convert_to_pdbqt("in.pdb", P("x.pdbqt"))
finally:
    exporter.subprocess.run = real_run
check(73, "non-zero exit reported with stderr",
      ok is False and "exploded" in message, str(message)[:60])


class Empty:
    returncode = 0
    stderr = ""


exporter.subprocess.run = lambda *a, **k: Empty()
try:
    ok, message = exporter.convert_to_pdbqt("in.pdb", P("never_written.pdbqt"))
finally:
    exporter.subprocess.run = real_run
check(73, "exit 0 with no output file is still a failure",
      ok is False and "no output" in message.lower(), str(message)[:60])
check(73, "timeout is 120 seconds", exporter.OBABEL_TIMEOUT == 120, "")


# ---- F74-F78: structured report --------------------------------------------
protonation = {"hydrogens_added": True, "ph": 7.4,
               "ligands_preserved": ["BTN"],
               "ligand_status": [{"residue": "BTN", "atoms": 16,
                                  "action": "preserved", "reason": "no template"}],
               "nonstandard_replaced": ["MSE"], "loops_reconstructed": 2}
minimization = {"status": "partial", "force_field": "amber14",
                "energy_before_kJ_mol": 900.0, "energy_after_kJ_mol": -5000.0,
                "delta_energy_kJ_mol": -5900.0, "converged": True,
                "energy_decreased": True, "excluded_residues": ["BTN"],
                "gbsa_used": True, "error": None}
export_status = {"target": "vina", "succeeded": True}

report = reporter.build_report(
    filename="1crn.pdb", chains_detected=["A", "B"], chains_retained=["A"],
    waters_removed=12, waters_retained=3,
    heteroatoms_removed=["SO4"], heteroatoms_retained=["BTN"],
    ph_used=7.4, missing_residues=[{"chain": "A", "residue": "GLY"}],
    atoms_before=327, atoms_after=642,
    protonation=protonation, minimization_stats=minimization,
    docking_export=export_status, settings={"ph": 7.4},
    processing_time_s=3.14159, warnings=["something to note"])

for key in ("generated_at", "input_file", "chains", "water_molecules_removed",
            "heteroatoms", "protonation", "missing_residues_detected",
            "atom_counts"):
    check(74, f"report contains {key}", key in report, "")
check(74, "chains split into detected and retained",
      report["chains"] == {"detected": ["A", "B"], "retained": ["A"]}, "")
check(74, "atom delta computed",
      report["atom_counts"]["delta"] == 642 - 327, str(report["atom_counts"]))
check(74, "timestamp is UTC", report["generated_at"].endswith("UTC"), "")

check(75, "minimisation stats embedded",
      report["energy_minimization"] == minimization, "")
check(76, "docking target recorded when the export succeeded",
      report.get("docking_target") == "vina", str(report.get("docking_target")))

failed_export = {"target": "vina", "succeeded": False, "error": "obabel missing"}
r2 = reporter.build_report(
    filename="x.pdb", chains_detected=[], chains_retained=[], waters_removed=0,
    heteroatoms_removed=[], heteroatoms_retained=[], ph_used=7.4,
    missing_residues=[], atoms_before=1, atoms_after=1,
    docking_export=failed_export)
check(76, "a failed export does not claim a docking target",
      "docking_target" not in r2 and r2["docking_export"]["succeeded"] is False, "")

check(77, "processing time recorded and rounded",
      report["processing_time_seconds"] == 3.14, str(report.get("processing_time_seconds")))

check(78, "ligand protonation status carried into the report",
      report["protonation"]["ligand_status"] == protonation["ligand_status"], "")
check(78, "hydrogens_added comes from the protonator, not assumed",
      report["hydrogens_added"] is True
      and reporter.build_report(
          filename="x", chains_detected=[], chains_retained=[], waters_removed=0,
          heteroatoms_removed=[], heteroatoms_retained=[], ph_used=7,
          missing_residues=[], atoms_before=0, atoms_after=0,
          protonation={"hydrogens_added": False})["hydrogens_added"] is False, "")

# ---- F79: text report ------------------------------------------------------
text = reporter.report_to_text(report)
for heading in ("PREPARATION REPORT", "CHAINS", "WATER MOLECULES",
                "HETEROATOMS", "PROTONATION", "STRUCTURE",
                "ENERGY MINIMIZATION", "DOCKING EXPORT", "WARNINGS"):
    check(79, f"text report has a {heading} section", heading in text, "")
check(79, "delta rendered with a sign", "+315" in text, "")
check(79, "ligand reason rendered", "no template" in text, "")
check(79, "minimisation status rendered", "partial" in text, "")
check(79, "warning rendered", "something to note" in text, "")
check(79, "partial report renders instead of raising",
      "x.pdb" in reporter.report_to_text({"input_file": "x.pdb"}), "")
check(79, "empty report renders", "BIOPREP" in reporter.report_to_text({}), "")
check(79, "None renders", "BIOPREP" in reporter.report_to_text(None), "")

# ---- F80: quick atom counter -----------------------------------------------
by_counter = reporter.count_atoms_in_pdb(PREPARED)
by_hand = sum(1 for l in open(PREPARED) if l.startswith(("ATOM", "HETATM")))
check(80, "counts ATOM and HETATM lines", by_counter == by_hand,
      f"{by_counter} vs {by_hand}")
check(80, "missing file returns zero rather than raising",
      reporter.count_atoms_in_pdb(P("nope.pdb")) == 0, "")

mixed = P("mixed.pdb")
with open(mixed, "w") as fh:
    fh.write("REMARK header\n")
    fh.write("ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N  \n")
    fh.write("HETATM    2 ZN    ZN A   2       1.000   2.000   3.000  1.00  0.00          ZN  \n")
    fh.write("TER\nEND\n")
check(80, "ignores non-coordinate records",
      reporter.count_atoms_in_pdb(mixed) == 2,
      str(reporter.count_atoms_in_pdb(mixed)))

start = time.time()
reporter.count_atoms_in_pdb(PREPARED)
check(80, "counter is fast (no Biopython parse)", (time.time() - start) < 0.5, "")


# ---- report ----------------------------------------------------------------
print("=" * 82)
print(f"{'F':>3}  {'RESULT':<6}  DESCRIPTION")
print("=" * 82)
failed = 0
for feature, description, passed, detail in sorted(RESULTS, key=lambda x: x[0]):
    mark = "ok" if passed else "FAIL"
    if not passed:
        failed += 1
    line = f"{feature:>3}  {mark:<6}  {description}"
    if detail and not passed:
        line += f"   [{detail}]"
    print(line)
print("=" * 82)
print(f"{len(RESULTS)} checks across features 70-80, {failed} failed")
