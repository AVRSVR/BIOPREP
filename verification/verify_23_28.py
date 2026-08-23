"""Independent verification of catalogue features 23-28 (analyzer.py)."""
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

from bioprep.core import io as bio                                  # noqa: E402
from bioprep.core.analyzer import analyze_structure, detect_missing_residues  # noqa: E402

TMP = tempfile.mkdtemp(prefix="anz_")
P = lambda n: os.path.join(TMP, n)          # noqa: E731
RESULTS = []


def check(feature, description, passed, detail=""):
    RESULTS.append((feature, description, passed, detail))


def het(serial, name, resn, chain, resseq, x, y, z, elem):
    return (f"HETATM{serial:>5} {name:<4} {resn:>3} {chain:1}{resseq:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00          {elem:>2}  \n")


def write(name, lines):
    path = P(name)
    with open(path, "w") as fh:
        fh.writelines(lines)
        fh.write("END\n")
    return path


PROT = [l for l in open("1crn.pdb") if l.startswith("ATOM")]
SEQRES = [l for l in open("1crn.pdb") if l.startswith("SEQRES")]
N_PROT = len(PROT)


# ---- F23: chain enumeration ------------------------------------------------
multi = PROT + [l[:21] + "B" + l[22:] for l in PROT] + \
        [l[:21] + "C" + l[22:] for l in PROT]
meta = analyze_structure(bio.load_pdb(write("3ch.pdb", multi)))
check(23, "lists every chain id", meta["chains"] == ["A", "B", "C"], str(meta["chains"]))
check(23, "chain ids are unique and sorted",
      meta["chains"] == sorted(set(meta["chains"])), "")

single = analyze_structure(bio.load_pdb(write("1ch.pdb", PROT)))
check(23, "single-chain structure reports one chain",
      single["chains"] == ["A"], str(single["chains"]))


# ---- F24: water count ------------------------------------------------------
waters = [het(9000 + i, "O", nm, "A", 700 + i, 40.0 + i * 4, 12, 12, "O")
          for i, nm in enumerate(["HOH", "WAT", "H2O", "TIP", "SOL", "DOD"])]
m = analyze_structure(bio.load_pdb(write("wat.pdb", PROT + waters)))
check(24, "counts waters of every spelling", m["water_count"] == 6,
      f"{m['water_count']}")
check(24, "reports zero when there are no waters",
      analyze_structure(bio.load_pdb(write("nw.pdb", PROT)))["water_count"] == 0, "")


# ---- F25: heteroatom discovery ---------------------------------------------
hets = [het(9100, "C1", "LIG", "A", 800, 12, 12, 12, "C"),
        het(9101, "C2", "LIG", "A", 800, 13, 12, 12, "C"),   # same residue twice
        het(9102, "S", "SO4", "A", 801, 16, 12, 12, "S"),
        het(9103, "ZN", " ZN", "A", 802, 20, 12, 12, "ZN")]
m = analyze_structure(bio.load_pdb(write("het.pdb", PROT + hets + waters)))
check(25, "lists unique heteroatom names",
      m["heteroatoms"] == ["LIG", "SO4", "ZN"], str(m["heteroatoms"]))
check(25, "waters excluded from the heteroatom list",
      not any(w in m["heteroatoms"] for w in ("HOH", "SOL", "TIP")), "")
check(25, "duplicate residue names collapsed",
      m["heteroatoms"].count("LIG") == 1, "")


# ---- F26: total atom count -------------------------------------------------
expected = N_PROT + len(hets) + len(waters)
check(26, "counts every atom of every residue type",
      m["atoms_total"] == expected, f"{m['atoms_total']} vs {expected}")
check(26, "breakdown sums to the total",
      sum(m["atom_breakdown"].values()) == m["atoms_total"], str(m["atom_breakdown"]))
check(26, "breakdown splits protein/water/heteroatom correctly",
      m["atom_breakdown"] == {"protein": N_PROT, "water": len(waters),
                             "heteroatom": len(hets)}, str(m["atom_breakdown"]))

# multi-model must not multiply the counts
ens = []
for i in range(1, 4):
    ens.append(f"MODEL     {i:>4}\n")
    ens.extend(PROT)
    ens.append("ENDMDL\n")
me = analyze_structure(bio.load_pdb(write("ens.pdb", ens)))
check(26, "NMR ensemble does not multiply the atom count",
      me["atoms_total"] == N_PROT, f"{me['atoms_total']} vs {N_PROT}")
check(26, "model count still reported", me["model_count"] == 3, str(me["model_count"]))


# ---- F27: sequence gap detection -------------------------------------------
gapped = [l for l in PROT if not (20 <= int(l[22:26]) <= 24)]
g = analyze_structure(bio.load_pdb(write("gap.pdb", gapped)))
check(27, "detects a gap", len(g["sequence_gaps"]) == 1, str(g["sequence_gaps"]))
if g["sequence_gaps"]:
    gap = g["sequence_gaps"][0]
    check(27, "reports chain, from, to and size",
          gap == {"chain": "A", "from": 19, "to": 25, "missing_count": 5}, str(gap))

# no gap in a continuous chain
check(27, "no false positive on a continuous chain",
      analyze_structure(bio.load_pdb(write("cont.pdb", PROT)))["sequence_gaps"] == [], "")

# two gaps in one chain
two_gaps = [l for l in PROT if not (10 <= int(l[22:26]) <= 12)
            and not (30 <= int(l[22:26]) <= 31)]
g2 = analyze_structure(bio.load_pdb(write("g2.pdb", two_gaps)))
check(27, "detects multiple gaps in one chain", len(g2["sequence_gaps"]) == 2,
      str(g2["sequence_gaps"]))

# gaps reported per chain
per_chain = gapped + [l[:21] + "B" + l[22:] for l in PROT]
g3 = analyze_structure(bio.load_pdb(write("g3.pdb", per_chain)))
check(27, "gap attributed to the right chain",
      [x["chain"] for x in g3["sequence_gaps"]] == ["A"], str(g3["sequence_gaps"]))

# insertion codes must not look like a gap
ins = []
for l in PROT:
    if int(l[22:26]) == 5:
        ins.append(l)
        ins.append(l[:26] + "A" + l[27:])
        ins.append(l[:26] + "B" + l[27:])
    else:
        ins.append(l)
gi = analyze_structure(bio.load_pdb(write("ins.pdb", ins)))
check(27, "insertion codes do not create a false gap",
      gi["sequence_gaps"] == [], str(gi["sequence_gaps"]))

# heteroatoms interleaved must not break the run
inter = []
for l in PROT:
    inter.append(l)
    if int(l[22:26]) == 5 and l[12:16].strip() == "CA":
        inter.append(het(9500, "C1", "LIG", "A", 500, 50, 50, 50, "C"))
gh = analyze_structure(bio.load_pdb(write("inter.pdb", inter)))
check(27, "an interleaved ligand does not create a false gap",
      gh["sequence_gaps"] == [], str(gh["sequence_gaps"]))

# ensemble must not duplicate gaps
ens_gap = []
for i in range(1, 4):
    ens_gap.append(f"MODEL     {i:>4}\n")
    ens_gap.extend(gapped)
    ens_gap.append("ENDMDL\n")
ge = analyze_structure(bio.load_pdb(write("ensgap.pdb", ens_gap)))
check(27, "ensemble does not report the same gap three times",
      len(ge["sequence_gaps"]) == 1, str(ge["sequence_gaps"]))


# ---- F28: PDBFixer missing residue detection -------------------------------
with_seqres = write("seq.pdb", SEQRES + gapped)
missing = detect_missing_residues(with_seqres)
check(28, "detects residues present in SEQRES but absent from coordinates",
      len(missing) == 5, f"{len(missing)} found")
check(28, "reports a real chain id, not an index",
      {mm["chain"] for mm in missing} == {"A"}, str({mm["chain"] for mm in missing}))
check(28, "reports residue identity",
      all(isinstance(mm["residue"], str) and mm["residue"].isalpha()
          for mm in missing), str(missing[:2]))

# multi-chain: the chain label must follow the right chain
mc = SEQRES + gapped + [l[:21] + "B" + l[22:] for l in PROT]
missing_mc = detect_missing_residues(write("mc.pdb", mc))
check(28, "multi-chain: gap attributed to chain A only",
      {mm["chain"] for mm in missing_mc} == {"A"},
      str({mm["chain"] for mm in missing_mc}))

# without SEQRES there is nothing to compare against
check(28, "no SEQRES means nothing reported (not a crash)",
      detect_missing_residues(write("noseq.pdb", gapped)) == [], "")

# a complete structure has none
check(28, "complete structure reports no missing residues",
      detect_missing_residues(write("full.pdb", SEQRES + PROT)) == [], "")

# unreadable input must not raise
check(28, "unreadable file returns empty rather than raising",
      detect_missing_residues(P("does_not_exist.pdb")) == [], "")


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
print(f"{len(RESULTS)} checks across features 23-28, {failed} failed")
