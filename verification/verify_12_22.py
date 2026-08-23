"""Independent verification of catalogue features 12-22 (protonator.py)."""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)

import glob
import os
import sys
import tempfile
import warnings

warnings.filterwarnings("ignore")
BIOPREP = _os.path.join(_ROOT, "bioprep")
sys.path.insert(0, BIOPREP)
os.chdir(BIOPREP)

from bioprep.core import io as bio                     # noqa: E402
from bioprep.core import protonator                    # noqa: E402

TMP = tempfile.mkdtemp(prefix="prot_")
P = lambda n: os.path.join(TMP, n)        # noqa: E731
RESULTS = []


def check(feature, description, passed, detail=""):
    RESULTS.append((feature, description, passed, detail))


def het(serial, name, resn, chain, resseq, x, y, z, elem):
    return (f"HETATM{serial:>5} {name:<4} {resn:>3} {chain:1}{resseq:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00          {elem:>2}  \n")


PROT = [l for l in open("1crn.pdb") if l.startswith("ATOM")]
CA = [l for l in PROT if l[12:16].strip() == "CA"][0]
CX, CY, CZ = float(CA[30:38]), float(CA[38:46]), float(CA[46:54])
LIG = [het(900, "C1", "LIG", "A", 900, CX + 6, CY, CZ, "C"),
       het(901, "C2", "LIG", "A", 900, CX + 7.5, CY, CZ, "C"),
       het(902, "O1", "LIG", "A", 900, CX + 8.2, CY + 1.1, CZ, "O")]


def count_h(path):
    return sum(1 for l in open(path)
               if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() == "H")


# ---- F12: pH-dependent hydrogen addition ----------------------------------
plain = P("plain.pdb")
with open(plain, "w") as fh:
    fh.write("".join(PROT) + "END\n")

counts = {}
for ph in (1.0, 7.4, 13.0):
    out = P(f"ph{ph}.pdb")
    protonator.add_hydrogens(plain, out, ph=ph)
    counts[ph] = count_h(out)

check(12, "hydrogens are added at all", counts[7.4] > 0, f"{counts[7.4]} H at pH 7.4")
check(12, "hydrogen count varies with pH", len(set(counts.values())) > 1, str(counts))
check(12, "low pH gives more hydrogens than high pH",
      counts[1.0] > counts[13.0], f"pH1={counts[1.0]} pH13={counts[13.0]}")

res = protonator.add_hydrogens(plain, P("def.pdb"))
check(12, "default pH is 7.4", res["ph"] == 7.4, str(res["ph"]))


# ---- F13/F14/F19: separation, preservation, re-merge ----------------------
withlig = P("withlig.pdb")
with open(withlig, "w") as fh:
    fh.write("".join(PROT) + "".join(LIG) + "END\n")

out = P("lig_out.pdb")
res = protonator.add_hydrogens(withlig, out, ph=7.4)

src_lig = [l for l in open(withlig) if l[17:20].strip() == "LIG"]
out_lig = [l for l in open(out) if l[17:20].strip() == "LIG"]
check(13, "ligand held out of PDBFixer", res["ligands_preserved"] == ["LIG"],
      str(res["ligands_preserved"]))
check(14, "ligand atom count unchanged",
      len(out_lig) == len(src_lig), f"{len(src_lig)} -> {len(out_lig)}")
check(14, "ligand coordinates byte-identical",
      all(s[30:54] == o[30:54] for s, o in zip(src_lig, out_lig)), "")
check(14, "no hydrogens added to the ligand",
      all(l[76:78].strip() != "H" for l in out_lig), "")
check(19, "protein and ligand merged into one file",
      any(l.startswith("ATOM") for l in open(out)) and out_lig, "")
check(19, "merged file re-parses cleanly",
      sum(1 for _ in bio.load_pdb(out).get_atoms()) > 327, "")
check(19, "protein was protonated despite the ligand",
      count_h(out) > 0, f"{count_h(out)} H")

serials = [int(l[6:11]) for l in open(out) if l[:6] in ("ATOM  ", "HETATM")]
check(19, "atom serials are unique and ascending after merge",
      serials == sorted(serials) and len(serials) == len(set(serials)), "")

# waters must stay WITH the protein, not be held out as ligands
wat = P("wat.pdb")
with open(wat, "w") as fh:
    fh.write("".join(PROT)
             + het(950, "O", "HOH", "A", 700, CX + 2.8, CY, CZ, "O") + "END\n")
res_w = protonator.add_hydrogens(wat, P("wat_out.pdb"), ph=7.4)
hoh = [l for l in open(P("wat_out.pdb")) if l[17:20].strip() == "HOH"]
check(13, "waters go through PDBFixer, not held out as ligands",
      res_w["ligands_preserved"] == [] and len(hoh) > 1,
      f"preserved={res_w['ligands_preserved']} HOH atoms={len(hoh)}")


# ---- F15: non-standard residue replacement --------------------------------
mse = []
for line in PROT:
    if line[17:20] == "MET":
        mse.append("HETATM" + line[6:17] + "MSE" + line[20:])
    else:
        mse.append(line)
if not any(l[17:20] == "MSE" for l in mse):     # 1CRN has no MET; convert a THR
    mse = []
    for line in PROT:
        if line[17:20] == "THR" and line[22:26].strip() == "1":
            mse.append("HETATM" + line[6:17] + "MSE" + line[20:])
        else:
            mse.append(line)
msefile = P("mse.pdb")
with open(msefile, "w") as fh:
    fh.write("".join(mse) + "END\n")
res_m = protonator.add_hydrogens(msefile, P("mse_out.pdb"), ph=7.4)
check(15, "non-standard residue detected and replaced",
      bool(res_m["nonstandard_replaced"]), str(res_m["nonstandard_replaced"]))
check(15, "replaced residue no longer present as MSE",
      not any(l[17:20].strip() == "MSE" for l in open(P("mse_out.pdb"))), "")


# ---- F16: loop reconstruction ---------------------------------------------
seqres = [l for l in open("1crn.pdb") if l.startswith("SEQRES")]
gapped = [l for l in PROT if not (20 <= int(l[22:26]) <= 24)]
gapfile = P("gap.pdb")
with open(gapfile, "w") as fh:
    fh.write("".join(seqres + gapped) + "END\n")

res_off = protonator.add_hydrogens(gapfile, P("gap_off.pdb"),
                                   reconstruct_loops=False)
res_on = protonator.add_hydrogens(gapfile, P("gap_on.pdb"),
                                  reconstruct_loops=True)


def resseqs(path):
    return {int(l[22:26]) for l in open(path)
            if l.startswith("ATOM") and l[17:20].strip() != "HOH"}


check(16, "loops NOT rebuilt when the flag is off",
      res_off["loops_reconstructed"] == 0
      and not (resseqs(P("gap_off.pdb")) & {20, 21, 22, 23, 24}),
      f"count={res_off['loops_reconstructed']}")
check(16, "loops rebuilt when the flag is on",
      res_on["loops_reconstructed"] == 5, f"count={res_on['loops_reconstructed']}")
check(16, "rebuilt residues actually appear in the output",
      {20, 21, 22, 23, 24} <= resseqs(P("gap_on.pdb")),
      f"missing {sorted({20,21,22,23,24} - resseqs(P('gap_on.pdb')))}")


# ---- F17: missing atom addition -------------------------------------------
stripped = []
for line in PROT:
    if line[22:26].strip() == "2" and line[12:16].strip() in ("CG1", "CG2", "OG1"):
        continue
    stripped.append(line)
sfile = P("strip.pdb")
with open(sfile, "w") as fh:
    fh.write("".join(stripped) + "END\n")

protonator.add_hydrogens(sfile, P("strip_off.pdb"), add_missing_atoms=False)
protonator.add_hydrogens(sfile, P("strip_on.pdb"), add_missing_atoms=True)


def heavy_of_res2(path):
    return {l[12:16].strip() for l in open(path)
            if l.startswith("ATOM") and l[22:26].strip() == "2"
            and l[76:78].strip() != "H"}


check(17, "missing side-chain atoms added when the flag is on",
      len(heavy_of_res2(P("strip_on.pdb"))) > len(heavy_of_res2(P("strip_off.pdb"))),
      f"off={sorted(heavy_of_res2(P('strip_off.pdb')))} "
      f"on={sorted(heavy_of_res2(P('strip_on.pdb')))}")


# ---- F18: error handling ---------------------------------------------------
try:
    ligonly = P("ligonly.pdb")
    with open(ligonly, "w") as fh:
        fh.write("".join(LIG) + "END\n")
    protonator.add_hydrogens(ligonly, P("ligonly_out.pdb"))
    check(18, "ligand-only input raises a clear error", False, "no error")
except ValueError as exc:
    check(18, "ligand-only input raises a clear error",
          "No biopolymer atoms" in str(exc), "")

check(18, "result carries a warnings list", isinstance(res.get("warnings"), list), "")


# ---- F20: status report ----------------------------------------------------
expected_keys = {"hydrogens_added", "ph", "ligands_preserved",
                 "nonstandard_replaced", "loops_reconstructed", "warnings"}
check(20, "returns a status dict with the documented keys",
      expected_keys <= set(res.keys()), str(sorted(res.keys())))
check(20, "reports which ligands were skipped", res["ligands_preserved"] == ["LIG"], "")
check(20, "reports hydrogens_added truthfully", res["hydrogens_added"] is True, "")


# ---- F21: temp file cleanup ------------------------------------------------
pattern = os.path.join(tempfile.gettempdir(), "*_protein.pdb")
pattern2 = os.path.join(tempfile.gettempdir(), "*_protonated.pdb")
before = set(glob.glob(pattern)) | set(glob.glob(pattern2))
for _ in range(3):
    protonator.add_hydrogens(withlig, P("tmp_out.pdb"))
after = set(glob.glob(pattern)) | set(glob.glob(pattern2))
check(21, "no temp files left behind", after == before,
      f"{len(after - before)} leaked")

try:
    protonator.add_hydrogens(P("ligonly.pdb"), P("x.pdb"))
except ValueError:
    pass
after_err = set(glob.glob(pattern)) | set(glob.glob(pattern2))
check(21, "no temp files left behind on the error path",
      after_err == before, f"{len(after_err - before)} leaked")


# ---- F22: chain ID preservation --------------------------------------------
multi = PROT + [l[:21] + "B" + l[22:] for l in PROT] + [
    het(960, "C1", "LIG", "C", 900, CX + 6, CY, CZ, "C")]
mfile = P("multi.pdb")
with open(mfile, "w") as fh:
    fh.write("".join(multi) + "END\n")
protonator.add_hydrogens(mfile, P("multi_out.pdb"), ph=7.4)
chains = sorted({l[21] for l in open(P("multi_out.pdb"))
                 if l[:6] in ("ATOM  ", "HETATM")})
check(22, "chain identifiers preserved through protonation",
      chains == ["A", "B", "C"], f"{chains}")

resnums = [int(l[22:26]) for l in open(P("multi_out.pdb"))
           if l.startswith("ATOM") and l[21] == "A"]
check(22, "residue numbering preserved",
      min(resnums) == 1 and max(resnums) == 46, f"{min(resnums)}-{max(resnums)}")


# ---- report ----------------------------------------------------------------
print("=" * 80)
print(f"{'F':>3}  {'RESULT':<6}  DESCRIPTION")
print("=" * 80)
failed = 0
for feature, description, passed, detail in RESULTS:
    mark = "ok" if passed else "FAIL"
    if not passed:
        failed += 1
    line = f"{feature:>3}  {mark:<6}  {description}"
    if detail and not passed:
        line += f"   [{detail}]"
    print(line)
print("=" * 80)
print(f"{len(RESULTS)} checks across features 12-22, {failed} failed")
