"""Independent verification of catalogue features 1-11 against the live code."""
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
sys.path.insert(0, BIOPREP)
os.chdir(BIOPREP)

from Bio.PDB import Select, MMCIFIO                      # noqa: E402
from bioprep.core import io as bio                        # noqa: E402
from bioprep.core.cleaner import clean_structure, BioPrepSelect  # noqa: E402
from bioprep.core.analyzer import analyze_structure       # noqa: E402

TMP = tempfile.mkdtemp(prefix="verify_")
P = lambda n: os.path.join(TMP, n)          # noqa: E731
RESULTS = []


def check(feature, description, passed, detail=""):
    RESULTS.append((feature, description, passed, detail))


def het(serial, name, resn, chain, resseq, x, y, z, elem):
    return (f"HETATM{serial:>5} {name:<4} {resn:>3} {chain:1}{resseq:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00          {elem:>2}  \n")


PROT = [l for l in open("1crn.pdb") if l.startswith("ATOM")]
CA = [l for l in PROT if l[12:16].strip() == "CA"][0]
CX, CY, CZ = float(CA[30:38]), float(CA[38:46]), float(CA[46:54])


# ---- F1 -------------------------------------------------------------------
s = bio.load_pdb("1crn.pdb")
check(1, "loads PDB via PDBParser, returns Structure",
      sum(1 for _ in s.get_atoms()) == 327, "327 atoms")

messy = P("messy.pdb")
with open(messy, "w") as _fh:
    _fh.write("".join([PROT[0], PROT[0], PROT[0].replace(" A ", " B ")]) + "END\n")
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    bio.load_pdb(messy)
# A ResourceWarning raised by this harness is not a parser warning.
parser_warnings = [w for w in caught if w.category is not ResourceWarning]
check(1, "QUIET=True suppresses parser warnings", not parser_warnings,
      f"{[w.category.__name__ for w in parser_warnings]}")

try:
    bio.load_pdb(_ := P("junk.pdb")) if open(P("junk.pdb"), "w").write("junk\n") else None
    check(1, "non-structure file rejected", False, "no error raised")
except ValueError as exc:
    check(1, "non-structure file rejected", "No atoms" in str(exc), "clear ValueError")

w = MMCIFIO(); w.set_structure(bio.load_pdb("1crn.pdb")); w.save(P("1crn.cif"))
check(1, "mmCIF input accepted",
      sum(1 for _ in bio.load_structure(P("1crn.cif")).get_atoms()) == 327,
      "327 atoms from .cif")

import shutil                                              # noqa: E402
shutil.copy(P("1crn.cif"), P("sneaky.pdb"))
check(1, "format detected by content, not extension",
      bio.detect_format(P("sneaky.pdb")) == "mmcif", ".cif renamed .pdb")


# ---- F2 -------------------------------------------------------------------
ids = {n: bio.structure_id_for(n) for n in
       ("1crn.pdb", r"C:\d\1xyz.pdb", "1abc.v2.pdb", "/dir.pdb/file.pdb", "noext")}
check(2, "structure id derived from filename stem",
      list(ids.values()) == ["1crn", "1xyz", "1abc.v2", "file", "noext"], str(ids))


# ---- F3 -------------------------------------------------------------------
bio.save_pdb(bio.load_pdb("1crn.pdb"), P("rt.pdb"))
rcsb = next(l.rstrip("\n") for l in open("1crn.pdb") if l.startswith("ATOM"))
ours = next(l.rstrip("\n") for l in open(P("rt.pdb")) if l.startswith("ATOM"))
check(3, "writes 80-column lines matching RCSB layout",
      len(ours) == 80 and ours[76:78] == rcsb[76:78],
      f"len={len(ours)} element={ours[76:78]!r}")

metals = [het(900, "ZN", " ZN", "A", 300, 12, 12, 12, "ZN"),
          het(901, "FE", " FE", "A", 301, 16, 12, 12, "FE")]
open(P("met.pdb"), "w").write("".join(PROT) + "".join(metals) + "END\n")
st = bio.load_pdb(P("met.pdb")); bio.save_pdb(st, P("met_out.pdb"))
after = {a.get_parent().get_resname().strip(): a.element
         for a in bio.load_pdb(P("met_out.pdb")).get_atoms()
         if a.get_parent().id[0] != " "}
check(3, "two-character elements survive a save",
      after.get("ZN") == "ZN" and after.get("FE") == "FE", str(after))

from openmm.app import PDBFile                             # noqa: E402
sym = {a.residue.name.strip(): (a.element.symbol if a.element else None)
       for a in PDBFile(P("met_out.pdb")).topology.atoms()
       if a.residue.name.strip() in ("ZN", "FE")}
check(3, "OpenMM agrees on the elements", sym == {"ZN": "Zn", "FE": "Fe"}, str(sym))

rows = ["ATOM      1  N   ALA A   1       1.000   2.000   3.000  0.75 12.34           N  \n",
        "ATOM      2  CA AALA A  50      10.000  10.000  10.000  0.50 20.00           C  \n",
        "ATOM      3  CA BALA A  50      10.500  10.000  10.000  0.50 20.00           C  \n",
        "ATOM      4  N   GLY A  51A     11.000  11.000  11.000  1.00 15.00           N  \n"]
open(P("f.pdb"), "w").write("".join(rows) + "END\n")
bio.save_pdb(bio.load_pdb(P("f.pdb")), P("f_out.pdb"))
o = [l for l in open(P("f_out.pdb")) if l.startswith("ATOM")]
check(3, "occupancy, B-factor, altloc and insertion code preserved",
      len(o) == 4 and o[0][54:60].strip() == "0.75" and o[0][60:66].strip() == "12.34"
      and {l[16] for l in o if l[22:26].strip() == "50"} == {"A", "B"}
      and any(l[26] == "A" for l in o if l[22:26].strip() == "51"), "")

lig = [het(900, "C1", "LIG", "A", 900, 12.0, 12.0, 12.0, "C"),
       het(901, "C2", "LIG", "A", 900, 13.5, 12.0, 12.0, "C")]
open(P("c.pdb"), "w").write("".join(PROT) + "".join(lig) + "CONECT  900  901\n" + "END\n")
bio.save_pdb(bio.load_pdb(P("c.pdb")), P("c_out.pdb"), conect_source=P("c.pdb"))
check(3, "CONECT records carried across a save",
      sum(1 for l in open(P("c_out.pdb")) if l.startswith("CONECT")) > 0, "")

bio.save_pdb(bio.load_pdb("1crn.pdb"), P("a/b/c.pdb"))
check(3, "missing parent directories created", os.path.isfile(P("a/b/c.pdb")), "")


# ---- F4 -------------------------------------------------------------------
mix = [het(930, "C1", "LIG", "A", 800, CX + 10, CY, CZ, "C"),
       het(931, "S", "SO4", "A", 801, CX + 14, CY, CZ, "S"),
       het(932, "ZN", " ZN", "A", 802, CX + 18, CY, CZ, "ZN"),
       het(933, "O", "HOH", "A", 803, CX + 2.8, CY, CZ, "O")]
open(P("mix.pdb"), "w").write("".join(PROT) + "".join(mix) + "END\n")


class ProteinOnly(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0


stx = bio.load_pdb(P("mix.pdb"))
bio.save_pdb(stx, P("sel.pdb"), select=ProteinOnly())
txt = open(P("sel.pdb")).read()
check(4, "Select filters what reaches disk",
      "LIG" not in txt and "SO4" not in txt and "ATOM" in txt, "")


def hets(path):
    return sorted({l[17:20].strip() for l in open(path) if l.startswith("HETATM")})


def clean(**kw):
    st_ = bio.load_pdb(P("mix.pdb"))
    bio.save_pdb(st_, P("cl.pdb"), select=clean_structure(st_, **kw))
    return hets(P("cl.pdb"))


# ---- F5 -------------------------------------------------------------------
spellings = ["HOH", "WAT", "H2O", "TIP", "SOL", "DOD"]
wat = [het(940 + i, "O", nm, "A", 700 + i, CX + 40 + i * 4, CY, CZ, "O")
       for i, nm in enumerate(spellings)]
open(P("wat.pdb"), "w").write("".join(PROT) + "".join(wat) + het(960, "C1", "LIG", "A", 900, CX + 10, CY, CZ, "C") + "END\n")
stw = bio.load_pdb(P("wat.pdb"))
bio.save_pdb(stw, P("wat_out.pdb"), select=clean_structure(stw, remove_water=True))
check(5, "every water spelling removed", hets(P("wat_out.pdb")) == ["LIG"],
      f"remaining {hets(P('wat_out.pdb'))}")

stw = bio.load_pdb(P("wat.pdb"))
bio.save_pdb(stw, P("wat_keep.pdb"), select=clean_structure(stw, remove_water=False))
check(5, "waters retained when remove_water=False",
      all(n in hets(P("wat_keep.pdb")) for n in spellings), "")

meta = analyze_structure(bio.load_pdb(P("wat.pdb")))
check(5, "analyzer counts waters, does not list them as ligands",
      meta["water_count"] == len(spellings) and meta["heteroatoms"] == ["LIG"],
      f"count={meta['water_count']}")


# ---- F6 -------------------------------------------------------------------
for oxy in ("O", "OW", "OH2", "O1"):
    rows2 = [het(970, oxy, "HOH", "A", 750, CX + 2.8, CY, CZ, "O"),
             het(971, oxy, "HOH", "A", 751, CX + 80, CY, CZ, "O")]
    open(P("sw.pdb"), "w").write("".join(PROT) + "".join(rows2) + "END\n")
    s6 = bio.load_pdb(P("sw.pdb"))
    sel6 = clean_structure(s6, remove_water=True, keep_structural_waters=True)
    bio.save_pdb(s6, P("sw_out.pdb"), select=sel6)
    kept = sorted(int(l[22:26]) for l in open(P("sw_out.pdb")) if l.startswith("HETATM"))
    check(6, f"structural water found via oxygen named {oxy}", kept == [750],
          f"kept {kept}")

prot40 = PROT[:40]
two = prot40 + [l[:21] + "B" + l[22:] for l in prot40] + [
    het(980, "O", "HOH", "A", 700, CX + 2.8, CY, CZ, "O"),
    het(981, "O", "HOH", "B", 700, CX + 500, CY + 500, CZ + 500, "O")]
open(P("cw.pdb"), "w").write("".join(two) + "END\n")
s6b = bio.load_pdb(P("cw.pdb"))
sel6b = clean_structure(s6b, remove_water=True, keep_structural_waters=True)
bio.save_pdb(s6b, P("cw_out.pdb"), select=sel6b)
chains_kept = sorted({l[21] for l in open(P("cw_out.pdb")) if l.startswith("HETATM")})
check(6, "structural waters qualified by chain", chains_kept == ["A"],
      f"chains {chains_kept}")


# ---- F7 -------------------------------------------------------------------
two7 = PROT + [l[:21] + "B" + l[22:] for l in PROT] + [
    het(990, "C1", "LIG", "A", 800, CX, CY, CZ, "C"),
    het(991, "C1", "LIG", "B", 800, CX, CY, CZ, "C")]
open(P("ch.pdb"), "w").write("".join(two7) + "END\n")
s7 = bio.load_pdb(P("ch.pdb"))
bio.save_pdb(s7, P("ch_out.pdb"), select=clean_structure(s7, target_chains=["A"]))
ch = sorted({l[21] for l in open(P("ch_out.pdb")) if l.startswith(("ATOM", "HETATM"))})
check(7, "chain selection retains only the requested chain", ch == ["A"], f"{ch}")


# ---- F8/9/10 --------------------------------------------------------------
check(8, "named heteroatom removed, others kept",
      clean(remove_heteroatoms=["SO4"]) == ["LIG", "ZN"], "")
check(9, "'ALL' strips every heteroatom",
      clean(remove_heteroatoms=["ALL"]) == [], "")
check(10, "unlisted heteroatoms preserved by default",
      clean(remove_heteroatoms=[]) == ["LIG", "SO4", "ZN"], "")
check(10, "protect_ligands outranks 'ALL'",
      clean(remove_heteroatoms=["ALL"], protect_ligands=["LIG"]) == ["LIG"], "")
check(10, "residue names matched case-insensitively",
      clean(remove_heteroatoms=["so4"]) == ["LIG", "ZN"]
      and clean(remove_heteroatoms=["ALL"], protect_ligands=["lig"]) == ["LIG"], "")


# ---- F11 ------------------------------------------------------------------
sel11 = clean_structure(bio.load_pdb("1crn.pdb"))
check(11, "BioPrepSelect subclasses Bio.PDB.Select",
      isinstance(sel11, Select) and isinstance(sel11, BioPrepSelect), "")
check(11, "implements model, chain, residue and atom hooks",
      all(hasattr(sel11, m) for m in
          ("accept_model", "accept_chain", "accept_residue", "accept_atom")), "")


# ---- report ---------------------------------------------------------------
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
print(f"{len(RESULTS)} checks across features 1-11, {failed} failed")
