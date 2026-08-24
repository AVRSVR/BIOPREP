"""
Stress the pipeline with settings combined, on real deposited structures.

Every check so far exercises one or two options at a time. Defects tend to hide
in the interaction: SEQRES survived on its own and CONECT survived on its own,
but carrying both plus structural waters is what exposed the terminus bug. This
turns several options on at once and asks whether the result is self-consistent,
not merely produced.
"""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)

import base64
import io as _io
import json
import os
import sys
import warnings
import zipfile

warnings.filterwarnings("ignore")
BIOPREP = _os.path.join(_ROOT, "bioprep")
PDBS = _os.path.join(_HERE, "pdbs")
sys.path.insert(0, BIOPREP)
os.chdir(BIOPREP)

from bioprep import webapp                                   # noqa: E402
from bioprep.core import io as bio                           # noqa: E402

RESULTS = []
client = webapp.app.test_client()


def check(label, description, passed, detail=""):
    RESULTS.append((label, description, passed, detail))


def decode(b64):
    return base64.b64decode(b64).decode("utf-8", "replace")


def atom_lines(text, resname=None):
    rows = [l for l in text.splitlines() if l[:6] in ("ATOM  ", "HETATM")]
    if resname:
        rows = [l for l in rows if l[17:20].strip() == resname]
    return rows


def post(route, data):
    return client.post(route, data=data, content_type="multipart/form-data")


# ---------------------------------------------------------------- everything on
print("=" * 78)
print("STRESS 1  Every option at once on HIV-1 protease (1HSG)")
print("=" * 78)
print("Two chains, an inhibitor, crystallographic waters. Chain selection,")
print("structural waters, ligand protection, loop rebuilding, missing atoms,")
print("minimisation and a docking export, all together.")
print()

raw = open(os.path.join(PDBS, "1HSG.pdb"), "rb").read()
response = post("/api/process", {
    "file": (_io.BytesIO(raw), "1HSG.pdb"),
    "chains": json.dumps(["A"]),
    "remove_heteros": json.dumps(["ALL"]),
    "protect_ligands": json.dumps(["MK1"]),
    "keep_structural_waters": "true",
    "reconstruct_loops": "true",
    "add_missing_atoms": "true",
    "run_minimization": "true",
    "use_gbsa": "true",
    "docking_target": "vina",
    "ph": "7.0",
})
check("1", "request succeeds", response.status_code == 200, str(response.status_code))

if response.status_code == 200:
    body = response.get_json()
    report = body["report"]
    viewer = decode(body["viewer_pdb_b64"])
    download = decode(body["pdb_b64"])

    chains = {l[21] for l in atom_lines(viewer)}
    check("1", "only the requested chain survives", chains == {"A"}, str(chains))

    ligand = atom_lines(viewer, "MK1")
    check("1", "protected ligand survives remove-ALL", bool(ligand),
          f"{len(ligand)} atoms")

    waters = {(l[21], l[22:26]) for l in atom_lines(viewer, "HOH")}
    check("1", "structural waters retained", bool(waters), f"{len(waters)}")
    check("1", "report agrees with the file on waters",
          report["water_molecules_retained"] == len(waters),
          f"report {report['water_molecules_retained']} vs file {len(waters)}")

    prot = report["protonation"]
    check("1", "protonation reported truthfully", prot["hydrogens_added"] is True, "")
    check("1", "protected ligand named in the report",
          "MK1" in prot["ligands_preserved"], str(prot["ligands_preserved"]))
    check("1", "pH honoured", report["protonation_ph"] == 7.0,
          str(report["protonation_ph"]))

    mini = report.get("energy_minimization") or {}
    check("1", "minimisation ran", mini.get("status") not in (None, "failed"),
          f"{mini.get('status')}: {str(mini.get('error'))[:70]}")
    if mini.get("status") not in (None, "failed"):
        check("1", "energy decreased",
              mini["energy_after_kJ_mol"] < mini["energy_before_kJ_mol"],
              f"{mini['energy_before_kJ_mol']} -> {mini['energy_after_kJ_mol']}")
        check("1", "excluded residues are the protected ligand",
              mini["excluded_residues"] == ["MK1"] or mini["status"] == "full",
              str(mini["excluded_residues"]))
        check("1", "pocket atoms were restrained",
              mini.get("restrained_atoms", 0) > 0 or mini["status"] == "full",
              str(mini.get("restrained_atoms")))

    export = report.get("docking_export") or {}
    check("1", "docking export succeeded", export.get("succeeded") is True,
          str(export.get("error"))[:70])
    check("1", "download is a PDBQT while the viewer stays PDB",
          body["filename"].endswith(".pdbqt")
          and "ATOM" in viewer and viewer is not download,
          body["filename"])
    check("1", "receptor carries no torsion tree",
          "ROOT" not in download and "BRANCH" not in download, "")

    # the ligand must not have moved, whatever else happened
    src_lig = atom_lines(open(os.path.join(PDBS, "1HSG.pdb")).read(), "MK1")
    check("1", "ligand coordinates untouched end to end",
          [l[30:54] for l in src_lig] == [l[30:54] for l in ligand],
          f"{len(src_lig)} in, {len(ligand)} out")

    print(f"   chains={chains}  waters={len(waters)}  MK1 atoms={len(ligand)}")
    print(f"   minimisation {mini.get('status')} "
          f"{mini.get('energy_before_kJ_mol')} -> {mini.get('energy_after_kJ_mol')}")
    print(f"   restrained atoms={mini.get('restrained_atoms')}  "
          f"loops rebuilt={prot.get('loops_reconstructed')}  "
          f"terminals={prot.get('terminals_repaired')}")
    print(f"   warnings: {len(report.get('warnings', []))}")
    for w in report.get("warnings", []):
        print(f"      - {w[:96]}")

# ---------------------------------------------------------------- mmCIF + all
print()
print("=" * 78)
print("STRESS 2  The same options, from mmCIF input")
print("=" * 78)

from Bio.PDB import MMCIFIO                                   # noqa: E402
import tempfile                                               # noqa: E402
tmp = tempfile.mkdtemp()
cif = os.path.join(tmp, "1hsg.cif")
writer = MMCIFIO()
writer.set_structure(bio.load_pdb(os.path.join(PDBS, "1HSG.pdb")))
writer.save(cif)

response = post("/api/process", {
    "file": (_io.BytesIO(open(cif, "rb").read()), "1hsg.cif"),
    "protect_ligands": json.dumps(["MK1"]),
    "keep_structural_waters": "true",
    "run_minimization": "true",
})
check("2", "mmCIF accepted with everything on", response.status_code == 200,
      str(response.status_code))
if response.status_code == 200:
    report = response.get_json()["report"]
    converted = [w for w in report.get("warnings", []) if "MMCIF" in w.upper()]
    check("2", "conversion recorded in the report", bool(converted),
          str(report.get("warnings"))[:80])
    mini = report.get("energy_minimization") or {}
    check("2", "minimisation ran from mmCIF",
          mini.get("status") not in (None, "failed"), str(mini.get("status")))
    print(f"   status={mini.get('status')} "
          f"{mini.get('energy_before_kJ_mol')} -> {mini.get('energy_after_kJ_mol')}")

# ---------------------------------------------------------------- metal site
print()
print("=" * 78)
print("STRESS 3  Metal site: zinc must not be silently discarded (4INS)")
print("=" * 78)

raw = open(os.path.join(PDBS, "4INS.pdb"), "rb").read()
response = post("/api/process", {
    "file": (_io.BytesIO(raw), "4INS.pdb"),
    "protect_ligands": json.dumps(["ZN"]),
    "run_minimization": "true",
})
check("3", "request succeeds", response.status_code == 200, str(response.status_code))
if response.status_code == 200:
    body = response.get_json()
    viewer = decode(body["viewer_pdb_b64"])
    zinc = atom_lines(viewer, "ZN")
    check("3", "zinc retained", bool(zinc), f"{len(zinc)} atoms")
    # Case-insensitive on purpose. RCSB and our own Biopython-path writer use
    # 'ZN'; OpenMM's own PDBFile.writeFile (used for every post-minimisation
    # file) writes 'Zn'. Both identify the same element correctly to every
    # tool that reads this file - documented as a minor known inconsistency,
    # not chased further, unlike the earlier column-offset bug that actually
    # changed the identified element.
    check("3", "zinc element still identifies as Zn/ZN",
          all(l[76:78].strip().upper() == "ZN" for l in zinc),
          str({l[76:78] for l in zinc}))
    mini = body["report"].get("energy_minimization") or {}
    check("3", "minimisation ran with a metal present",
          mini.get("status") not in (None, "failed"), str(mini.get("status")))
    print(f"   ZN atoms={len(zinc)} elements={ {l[76:78].strip() for l in zinc} }")
    print(f"   status={mini.get('status')}")

# ---------------------------------------------------------------- NMR ensemble
print()
print("=" * 78)
print("STRESS 4  A multi-model file must not be welded into one chimera")
print("=" * 78)

single = [l for l in open("1crn.pdb") if l.startswith("ATOM")]
ensemble = []
for index in range(1, 4):
    ensemble.append(f"MODEL     {index:>4}\n")
    ensemble.extend(single)
    ensemble.append("ENDMDL\n")
ensemble.append("END\n")

response = post("/api/process", {
    "file": (_io.BytesIO("".join(ensemble).encode()), "ens.pdb"),
    "run_minimization": "true",
})
check("4", "request succeeds", response.status_code == 200, str(response.status_code))
if response.status_code == 200:
    body = response.get_json()
    viewer = decode(body["viewer_pdb_b64"])
    heavy = [l for l in atom_lines(viewer) if l[76:78].strip() != "H"]
    residues = {(l[21], l[22:26]) for l in atom_lines(viewer)}
    check("4", "one model's worth of heavy atoms, not three",
          len(heavy) == len(single), f"{len(heavy)} vs {len(single)}")
    check("4", "no residue appears three times", len(residues) == 46,
          f"{len(residues)} residues")
    mini = body["report"].get("energy_minimization") or {}
    check("4", "ensemble minimises", mini.get("status") not in (None, "failed"),
          str(mini.get("status")))
    print(f"   heavy atoms={len(heavy)} (one model has {len(single)})  "
          f"residues={len(residues)}  status={mini.get('status')}")

# ---------------------------------------------------------------- report
print()
print("=" * 78)
print(f"{'#':>3}  {'RESULT':<6}  DESCRIPTION")
print("=" * 78)
failed = 0
for label, description, passed, detail in RESULTS:
    mark = "ok" if passed else "FAIL"
    if not passed:
        failed += 1
    line = f"{label:>3}  {mark:<6}  {description}"
    if detail and not passed:
        line += f"   [{detail}]"
    print(line)
print("=" * 78)
print(f"{len(RESULTS)} checks across features 84-110, {failed} failed")
