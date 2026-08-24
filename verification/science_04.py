"""
Science audit 4: structural waters through the whole pipeline.

Waters are the part most easily lost without anyone noticing - they are
numerous, individually unimportant-looking, and nothing downstream complains
when they vanish. Two separate defects did exactly that:

  * a stray END record ahead of them truncated the file handed to PDBFixer
  * no TER between polymer and water made PDBFixer read the waters as a
    continuation of the protein chain, so the last amino acid stopped being
    the chain's end and no terminal OXT was added - which then failed every
    minimisation tier

So this checks the count in and the count out, on real structures, with
minimisation on.
"""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)

import os
import sys
import tempfile
import warnings

warnings.filterwarnings("ignore")
BIOPREP = _os.path.join(_ROOT, "bioprep")
PDBS = _os.path.join(_HERE, "pdbs")
sys.path.insert(0, BIOPREP)
os.chdir(BIOPREP)

from bioprep.core import io as bio, cleaner, protonator, minimizer  # noqa: E402

TMP = tempfile.mkdtemp(prefix="sci4_")
P = lambda n: os.path.join(TMP, n)      # noqa: E731
RESULTS = []


def check(feature, description, passed, detail=""):
    RESULTS.append((feature, description, passed, detail))


def residues(path, resname):
    return {(l[21], l[22:26]) for l in open(path)
            if l[:6] in ("ATOM  ", "HETATM") and l[17:20].strip() == resname}


def atom_count(path, resname):
    return sum(1 for l in open(path)
               if l[:6] in ("ATOM  ", "HETATM") and l[17:20].strip() == resname)


def hydrogens_on(path, resname):
    return sum(1 for l in open(path)
               if l[:6] in ("ATOM  ", "HETATM") and l[17:20].strip() == resname
               and l[76:78].strip() == "H")


print("=" * 78)
print("SCIENCE 4  Structural waters survive the pipeline, with minimisation on")
print("=" * 78)
print()

for pdb_id, ligand in (("1STP", "BTN"), ("1HSG", "MK1"), ("4INS", "ZN")):
    source = os.path.join(PDBS, f"{pdb_id}.pdb")
    structure = bio.load_pdb(source)
    select = cleaner.clean_structure(
        structure, remove_water=True, keep_structural_waters=True,
        protect_ligands=[ligand])

    cleaned = P(f"{pdb_id}_c.pdb")
    bio.save_pdb(structure, cleaned, select=select, source_pdb=source)
    kept = len(residues(cleaned, "HOH"))

    protonated = P(f"{pdb_id}_p.pdb")
    prot = protonator.add_hydrogens(cleaned, protonated)
    out = len(residues(protonated, "HOH"))
    water_h = hydrogens_on(protonated, "HOH")
    ligand_atoms = atom_count(protonated, ligand)

    minimised = P(f"{pdb_id}_m.pdb")
    stats = minimizer.minimize_structure(protonated, minimised)
    after_min = len(residues(minimised, "HOH"))

    print(f"{pdb_id}  ligand {ligand}")
    print(f"   waters kept by cleaning      : {kept}")
    print(f"   waters after protonation     : {out}")
    print(f"   hydrogens added to waters    : {water_h}")
    print(f"   waters after minimisation    : {after_min}")
    print(f"   ligand atoms preserved       : {ligand_atoms}")
    print(f"   terminals repaired           : {prot['terminals_repaired']}")
    print(f"   minimisation                 : {stats['status']}  "
          f"{stats['energy_before_kJ_mol']} -> {stats['energy_after_kJ_mol']}")
    print()

    check(6, f"{pdb_id}: no water lost in protonation", out == kept,
          f"{kept} -> {out}")
    check(6, f"{pdb_id}: waters protonated (2 H each)", water_h == 2 * out,
          f"{water_h} H for {out} waters")
    check(6, f"{pdb_id}: no water lost in minimisation", after_min == out,
          f"{out} -> {after_min}")
    check(14, f"{pdb_id}: ligand preserved", ligand_atoms > 0,
          f"{ligand_atoms} atoms")
    check(33, f"{pdb_id}: minimisation ran with waters present",
          stats["status"] != "failed", str(stats.get("error"))[:80])


print("=" * 78)
print(f"{'F':>3}  {'RESULT':<6}  DESCRIPTION")
print("=" * 78)
failed = 0
for feature, description, passed, detail in RESULTS:
    mark = "ok" if passed else "FAIL"
    if not passed:
        failed += 1
    line = f"{feature:>3}  {mark:<6}  {description}"
    if detail and not passed:
        line += f"   [{detail}]"
    print(line)
print("=" * 78)
print(f"{len(RESULTS)} checks across features 6-33, {failed} failed")
