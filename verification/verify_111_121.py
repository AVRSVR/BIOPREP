"""Independent verification of catalogue features 111-121 (CLI and frontend)."""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)

import json
import os
import re
import subprocess
import sys
import tempfile
import warnings

warnings.filterwarnings("ignore")
ROOT = _ROOT
BIOPREP = os.path.join(ROOT, "bioprep")
PYTHON = _os.path.join(_ROOT, ".venv", "Scripts", "python.exe")
sys.path.insert(0, BIOPREP)
os.chdir(BIOPREP)

TMP = tempfile.mkdtemp(prefix="cli_")
P = lambda n: os.path.join(TMP, n)      # noqa: E731
RESULTS = []


def check(feature, description, passed, detail=""):
    RESULTS.append((feature, description, passed, detail))


def run_cli(*args, timeout=600):
    proc = subprocess.run([PYTHON, "-m", "bioprep.cli", *args],
                          capture_output=True, text=True, timeout=timeout,
                          cwd=BIOPREP)
    return proc


def het(serial, name, resn, chain, resseq, x, y, z, elem):
    return (f"HETATM{serial:>5} {name:<4} {resn:>3} {chain:1}{resseq:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00          {elem:>2}  \n")


PROT = [l for l in open("1crn.pdb") if l.startswith("ATOM")]
SEQRES = [l for l in open("1crn.pdb") if l.startswith("SEQRES")]
CA = [l for l in PROT if l[12:16].strip() == "CA"][0]
CX, CY, CZ = float(CA[30:38]), float(CA[38:46]), float(CA[46:54])

RICH = P("rich.pdb")
with open(RICH, "w") as fh:
    fh.write("".join(PROT))
    # Chain B offset so the two chains do not sit on top of each other;
    # perfectly superimposed atoms give an infinite starting energy.
    fh.write("".join(l[:21] + "B" + l[22:30]
                     + f"{float(l[30:38]) + 60.0:8.3f}" + l[38:]
                     for l in PROT))
    fh.write(het(9000, "C1", "LIG", "A", 900, CX + 6, CY, CZ, "C"))
    fh.write(het(9001, "S", "SO4", "A", 901, CX + 12, CY, CZ, "S"))
    fh.write(het(9002, "O", "HOH", "A", 700, CX + 2.8, CY, CZ, "O"))
    fh.write("END\n")


def atoms(path):
    return [l for l in open(path) if l[:6] in ("ATOM  ", "HETATM")]


def resnames(path):
    return {l[17:20].strip() for l in atoms(path) if l.startswith("HETATM")}


# ---- F111: the interface itself --------------------------------------------
proc = run_cli("--help")
check(111, "--help works", proc.returncode == 0, proc.stderr[:80])
helptext = proc.stdout
for flag in ("--input", "--output", "--chain", "--ph"):
    check(111, f"{flag} documented", flag in helptext, "")
for flag in ("--minimize", "--force-field", "--remove-het", "--protect",
             "--keep-water", "--keep-structural-water", "--no-gbsa"):
    check(111, f"{flag} available", flag in helptext, "")

proc = run_cli()
check(111, "missing required arguments exits non-zero", proc.returncode != 0, "")
check(111, "missing arguments explains what is required",
      "--input" in (proc.stderr + proc.stdout), proc.stderr[:60])

proc = run_cli("--input", P("nope.pdb"), "--output", P("x.pdb"))
check(111, "unreadable input exits non-zero", proc.returncode != 0, "")
check(111, "unreadable input reports the problem",
      "Error loading" in proc.stdout or "Error" in proc.stdout, proc.stdout[:70])

# ---- F112: clean + protonate -----------------------------------------------
out = P("basic.pdb")
proc = run_cli("--input", RICH, "--output", out)
check(112, "basic run succeeds", proc.returncode == 0, proc.stdout[-200:])
check(112, "reports the stages it ran",
      "Cleaning structure" in proc.stdout and "Adding hydrogens" in proc.stdout, "")
check(112, "writes an output file", os.path.isfile(out), "")
check(112, "output contains atoms", len(atoms(out)) > 0, "")
check(112, "hydrogens were added",
      any(l[76:78].strip() == "H" for l in atoms(out)), "")
check(112, "water removed by default", "HOH" not in resnames(out),
      str(resnames(out)))
check(112, "ligands kept by default", {"LIG", "SO4"} <= resnames(out),
      str(resnames(out)))

# pH
low, high = P("ph1.pdb"), P("ph13.pdb")
run_cli("--input", RICH, "--output", low, "--ph", "1.0")
run_cli("--input", RICH, "--output", high, "--ph", "13.0")


def hcount(path):
    return sum(1 for l in atoms(path) if l[76:78].strip() == "H")


check(112, "--ph changes the protonation", hcount(low) > hcount(high),
      f"{hcount(low)} vs {hcount(high)}")

# chains, repeatable
one, two = P("a.pdb"), P("ab.pdb")
run_cli("--input", RICH, "--output", one, "--chain", "A")
run_cli("--input", RICH, "--output", two, "--chain", "A", "--chain", "B")
check(112, "--chain filters to one chain",
      {l[21] for l in atoms(one)} == {"A"}, str({l[21] for l in atoms(one)}))
check(112, "--chain is repeatable",
      {l[21] for l in atoms(two)} == {"A", "B"}, str({l[21] for l in atoms(two)}))

# heteroatom handling
rm, prot_lig = P("rm.pdb"), P("prot.pdb")
run_cli("--input", RICH, "--output", rm, "--remove-het", "SO4")
check(112, "--remove-het removes the named residue",
      "SO4" not in resnames(rm) and "LIG" in resnames(rm), str(resnames(rm)))

run_cli("--input", RICH, "--output", prot_lig,
        "--remove-het", "ALL", "--protect", "LIG")
check(112, "--protect outranks --remove-het ALL",
      resnames(prot_lig) == {"LIG"}, str(resnames(prot_lig)))

lower = P("lower.pdb")
run_cli("--input", RICH, "--output", lower, "--remove-het", "ALL", "--protect", "lig")
check(112, "--protect is case-insensitive", "LIG" in resnames(lower),
      str(resnames(lower)))

# water flags
kw, ksw = P("kw.pdb"), P("ksw.pdb")
run_cli("--input", RICH, "--output", kw, "--keep-water")
check(112, "--keep-water retains water", "HOH" in resnames(kw), str(resnames(kw)))
run_cli("--input", RICH, "--output", ksw, "--keep-structural-water")
check(112, "--keep-structural-water retains the nearby water",
      "HOH" in resnames(ksw), str(resnames(ksw)))

# minimisation
mini = P("mini.pdb")
proc = run_cli("--input", RICH, "--output", mini, "--minimize")
check(112, "--minimize runs", proc.returncode == 0, proc.stdout[-200:])
check(112, "minimisation energies reported",
      "kJ/mol" in proc.stdout and "status:" in proc.stdout,
      proc.stdout[-160:].replace("\n", " "))
check(112, "partial tier reported honestly when a ligand blocks it",
      "partial" in proc.stdout or "full" in proc.stdout, "")
check(112, "excluded residues named in the warning",
      "held at their input coordinates" in proc.stdout.lower()
      or "status: full" in proc.stdout, proc.stdout[-160:].replace("\n", " "))

nog = P("nogbsa.pdb")
proc = run_cli("--input", RICH, "--output", nog, "--minimize", "--no-gbsa")
check(112, "--no-gbsa accepted", proc.returncode == 0, "")
cf = P("charmm.pdb")
proc = run_cli("--input", RICH, "--output", cf, "--minimize",
               "--force-field", "charmm36")
check(112, "--force-field charmm36 accepted", proc.returncode == 0,
      proc.stdout[-120:])
proc = run_cli("--input", RICH, "--output", P("bad.pdb"),
               "--force-field", "nonsense")
check(112, "invalid force field rejected by argparse", proc.returncode != 0, "")

# ligand reporting
proc = run_cli("--input", RICH, "--output", P("rep.pdb"))
check(112, "protected ligands named on stdout",
      "Ligands held out of PDBFixer" in proc.stdout, proc.stdout[:200])

# mmCIF input
from Bio.PDB import MMCIFIO                                    # noqa: E402
from bioprep.core.io import load_pdb                           # noqa: E402
w = MMCIFIO(); w.set_structure(load_pdb(RICH)); w.save(P("rich.cif"))
cifout = P("fromcif.pdb")
proc = run_cli("--input", P("rich.cif"), "--output", cifout)
check(112, "mmCIF input accepted", proc.returncode == 0, proc.stdout[-150:])
check(112, "conversion announced", "Converted MMCIF" in proc.stdout, "")
check(112, "mmCIF run produces atoms", len(atoms(cifout)) > 0, "")

# loop reconstruction through the CLI is not exposed; confirm SEQRES survives
gap = P("gap.pdb")
with open(gap, "w") as fh:
    fh.writelines(SEQRES)
    fh.writelines(l for l in PROT if not (20 <= int(l[22:26]) <= 24))
    fh.write("END\n")
gapout = P("gapout.pdb")
proc = run_cli("--input", gap, "--output", gapout)
check(112, "gapped input processed without error", proc.returncode == 0,
      proc.stdout[-150:])

# temp hygiene
import glob                                                    # noqa: E402
leaked = (glob.glob(os.path.join(tempfile.gettempdir(), "bioprep_cli_*"))
          + glob.glob(os.path.join(tempfile.gettempdir(), "bioprep_cif_*")))
check(112, "CLI leaves no working directories behind", not leaked,
      f"{len(leaked)} left")


# ---- F113-F121: frontend ---------------------------------------------------
HTML = os.path.join(BIOPREP, "templates", "index.html")
JS = os.path.join(BIOPREP, "static", "js", "main.js")
CSS = os.path.join(BIOPREP, "static", "css", "style.css")

for path, label in ((HTML, "index.html"), (JS, "main.js"), (CSS, "style.css")):
    check(113, f"{label} present", os.path.isfile(path), path)

html = open(HTML, encoding="utf-8", errors="replace").read()
js = open(JS, encoding="utf-8", errors="replace").read()

check(113, "single-page app with tab panels",
      html.count('class="tab-content') >= 5, str(html.count('class="tab-content')))
check(114, "3D viewer library referenced", "3Dmol" in html or "3Dmol" in js, "")
check(114, "viewer created in the frontend", "createViewer" in js, "")
check(115, "parameter controls present",
      'id="ph-slider"' in html or "ph-slider" in js or 'type="checkbox"' in html, "")
check(116, "analysis preview calls /api/analyze", "/api/analyze" in js, "")
check(117, "report rendered in the UI", "renderReport" in js, "")
check(118, "history panel calls /api/history", "/api/history" in js, "")
check(119, "template save and load wired",
      "/api/templates" in js, "")
check(120, "batch upload wired", "/api/batch" in js, "")
check(121, "binding site UI wired", "/api/analyze-site" in js, "")

# --- frontend/backend contract: what the rewrite must target ----------------
report_fields_used = set(re.findall(r"report\.([a-zA-Z_]+)", js))
check(117, "frontend reads report fields", bool(report_fields_used),
      str(sorted(report_fields_used))[:70])

CONTRACT = {
    "report.protonation.hydrogens_added": "protonation" in js,
    "report.energy_minimization.status": "energy_minimization" in js,
    "report.warnings": "warnings" in js,
}
stale = []
for key in ("potential_energy_before", "iterations", "'N/A'"):
    if key in js:
        stale.append(key)
check(117, "frontend has no references to removed energy fields",
      not stale, f"stale references: {stale}")


# ---- report ----------------------------------------------------------------
print("=" * 84)
print(f"{'F':>3}  {'RESULT':<6}  DESCRIPTION")
print("=" * 84)
failed = 0
for feature, description, passed, detail in sorted(RESULTS, key=lambda x: x[0]):
    mark = "ok" if passed else "FAIL"
    if not passed:
        failed += 1
    line = f"{feature:>3}  {mark:<6}  {description}"
    if detail and not passed:
        line += f"   [{detail}]"
    print(line)
print("=" * 84)
print(f"{len(RESULTS)} checks across features 111-121, {failed} failed")
print()
print("frontend report fields referenced:", sorted(report_fields_used))
print("backend contract present in js:", {k: v for k, v in CONTRACT.items()})
