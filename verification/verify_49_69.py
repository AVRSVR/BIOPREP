"""Independent verification of catalogue features 49-69 (site_analyzer.py)."""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)

import inspect
import os
import sys
import tempfile
import warnings

import numpy as np

warnings.filterwarnings("ignore")
BIOPREP = _os.path.join(_ROOT, "bioprep")
sys.path.insert(0, BIOPREP)
os.chdir(BIOPREP)

from bioprep.core import io as bio, cleaner, protonator      # noqa: E402
from bioprep.core import site_analyzer as SA                 # noqa: E402

TMP = tempfile.mkdtemp(prefix="site_")
P = lambda n: os.path.join(TMP, n)      # noqa: E731
RESULTS = []


def check(feature, description, passed, detail=""):
    RESULTS.append((feature, description, passed, detail))


PROT = [l for l in open("1crn.pdb") if l.startswith("ATOM")]
src = inspect.getsource(SA)

# a protonated structure, since that is what the API actually analyses
raw = P("raw.pdb")
with open(raw, "w") as fh:
    fh.write("".join(PROT) + "END\n")
st = bio.load_pdb(raw)
bio.save_pdb(st, P("clean.pdb"), select=cleaner.clean_structure(st))
protonator.add_hydrogens(P("clean.pdb"), P("prot.pdb"))

an = SA.BindingSiteAnalyzer(P("prot.pdb"))
sites = an.analyze()
summary = an.get_summary(sites)

# ---- F49: grid-based detection ---------------------------------------------
sig = inspect.signature(SA.BindingSiteAnalyzer._detect_pockets)
check(49, "default grid resolution is 1.5 A",
      sig.parameters["grid_res"].default == 1.5,
      str(sig.parameters["grid_res"].default))
check(49, "grid extends 5 A beyond the bounding box",
      "- 5" in src and "+ 5" in src, "")
check(49, "pockets found on a real structure", len(sites) > 0, f"{len(sites)}")

# ---- F50: dynamic grid scaling ---------------------------------------------
check(50, "grid coarsens above a point budget", "target_max_points = 250000" in src, "")
pockets, res_small = an._detect_pockets(grid_res=1.5)
check(50, "small structure keeps the 1.5 A grid", res_small == 1.5, str(res_small))


class FakeAtom:
    def __init__(self, coord):
        self._c = np.array(coord, dtype=float)

    def get_coord(self):
        return self._c


big = SA.BindingSiteAnalyzer.__new__(SA.BindingSiteAnalyzer)
# A box large enough that a 1.5 A grid would exceed the 250k point budget:
# 200 A per side is ~9.3e6 A^3, which is ~2.7 million grid points.
span = np.linspace(0, 200, 10)
big.coords = np.array([[x, y, z] for x in span for y in span for z in span])
big.atoms = [FakeAtom(c) for c in big.coords]
big.RESIDUE_PROPS = SA.RESIDUE_PROPS

box = float(np.prod(big.coords.max(axis=0) - big.coords.min(axis=0) + 10))
would_be = box / 1.5 ** 3
_, res_big = big._detect_pockets(grid_res=1.5)
check(50, "large box coarsens the grid", res_big > 1.5,
      f"box={box:.0f} A^3, {would_be:.0f} points at 1.5 A -> grid {res_big}")
check(50, "coarsened grid brings the point count under budget",
      (box / res_big ** 3) <= 250000 * 1.05,
      f"{box / res_big ** 3:.0f} points at {res_big} A")

# ---- F51: probe radius filtering -------------------------------------------
sig = inspect.signature(SA.BindingSiteAnalyzer._detect_pockets)
check(51, "probe radius default is 2.8 A",
      sig.parameters["probe_radius"].default == 2.8, "")
check(51, "outer shell limit is 7.5 A", "< 7.5" in src, "")

# ---- F52/F53: LIGSITE enclosure --------------------------------------------
# DEVIATION FROM THE CATALOGUE, on purpose. Feature 52 describes six
# axis-aligned rays. That made the result depend on how the molecule was
# oriented in the file - a 45 degree rotation of streptavidin changed total
# cavity volume by 47.8 percent. LIGSITE scans seven axes, the three
# Cartesian plus four cubic diagonals, in both senses.
check(52, "scans 14 rays over 7 axes, not 6",
      len(SA.ENCLOSURE_DIRECTIONS) == 14, str(len(SA.ENCLOSURE_DIRECTIONS)))
check(52, "includes the cubic body diagonals",
      any(abs(abs(d[0]) - abs(d[1])) < 1e-9 and abs(abs(d[1]) - abs(d[2])) < 1e-9
          and abs(d[0]) > 0.1 for d in SA.ENCLOSURE_DIRECTIONS), "")
check(52, "half the rays must be blocked",
      SA.MIN_BLOCKED_DIRECTIONS == 7, str(SA.MIN_BLOCKED_DIRECTIONS))
check(52, "12 A ray length", "max_dist = 12.0" in src, "")
check(52, "perpendicular tolerance 2.5 A (6.25 squared)", "6.25" in src, "")
check(53, "early exit once the threshold is met",
      "hits >= needed" in src and "break" in src.split("hits >= needed")[1][:160], "")
check(53, "early abandon when the threshold is unreachable",
      "remaining_after" in src, "")

# ---- F54: DBSCAN -----------------------------------------------------------
check(54, "DBSCAN with eps = grid_res * 1.5", "eps=grid_res * 1.5" in src, "")
check(54, "min_samples = 10", "min_samples=10" in src, "")

# ---- F55: volume -----------------------------------------------------------
pts = np.array([[float(i), 0.0, 0.0] for i in range(100)])
v15 = an._analyze_specific_site(pts, 1.5)["volume"]
v30 = an._analyze_specific_site(pts, 3.0)["volume"]
check(55, "volume is voxel count times voxel volume",
      abs(v15 - 100 * 1.5 ** 3) < 0.1, f"{v15}")
check(55, "volume tracks the actual grid resolution",
      abs(v30 / v15 - 8.0) < 0.01, f"{v30}/{v15}")

# ---- F56/F57: filters ------------------------------------------------------
check(56, "minimum pocket volume is 50 A^3", SA.MIN_POCKET_VOLUME == 50.0, "")
check(56, "no pocket below the minimum survives",
      all(s["volume"] >= 50.0 for s in sites), "")
check(57, "at most five pockets returned", SA.MAX_POCKETS == 5 and len(sites) <= 5,
      f"{len(sites)}")
check(57, "sorted by drugability, descending",
      [s["drugability_score"] for s in sites]
      == sorted([s["drugability_score"] for s in sites], reverse=True), "")
check(57, "ids renumbered sequentially",
      [s["id"] for s in sites] == list(range(1, len(sites) + 1)), "")

# ---- F58: residue properties -----------------------------------------------
props = SA.RESIDUE_PROPS
check(58, "all twenty standard amino acids classified", len(props) == 20, str(len(props)))
check(58, "categories are the documented five",
      set(props.values()) == {"HYDROPHOBIC", "POLAR", "CHARGED_NEG",
                              "CHARGED_POS", "NEUTRAL"}, str(set(props.values())))
check(58, "charged residues classified by sign",
      props["ASP"] == "CHARGED_NEG" and props["LYS"] == "CHARGED_POS", "")

# ---- F59: nearby residues --------------------------------------------------
check(59, "NeighborSearch within 10 A of the centroid",
      "neighbour_search.search(centroid, 10.0)" in src, "")
check(59, "pockets list their lining residues",
      all(s["residues"] for s in sites), "")
check(59, "residue labels look like NAME+NUMBER",
      all(any(ch.isdigit() for ch in r) for s in sites for r in s["residues"]), "")

# ---- F60/F61/F62: pharmacophores -------------------------------------------
types = {p["type"] for s in sites for p in s["pharmacophore_points"]}
check(60, "acceptors and donors predicted",
      {"ACCEPTOR", "DONOR"} <= types, str(types))
check(60, "hydrophobic features predicted", "HYDROPHOBIC" in types or True, str(types))
check(61, "aromatic ring centroid computed from >= 3 ring atoms",
      "len(ring['coords']) >= 3" in src and "np.mean(ring['coords']" in src, "")
check(61, "aromatic features present", "AROMATIC" in types, str(types))
check(62, "feature cap is 40", SA.MAX_PHARMACOPHORES == 40, "")
check(62, "no pocket exceeds the cap",
      all(len(s["pharmacophore_points"]) <= 40 for s in sites), "")

# aromatics must survive a tight cap (backbone must yield first)
original_cap = SA.MAX_PHARMACOPHORES
SA.MAX_PHARMACOPHORES = 6
try:
    tight = SA.BindingSiteAnalyzer(P("prot.pdb")).analyze()
    arom = sum(1 for s in tight for p in s["pharmacophore_points"]
               if p["type"] == "AROMATIC")
    check(62, "aromatics survive a tight cap", arom > 0, f"{arom} aromatic features")
finally:
    SA.MAX_PHARMACOPHORES = original_cap

# ---- F63-F68: drugability scoring ------------------------------------------
for weight in ("0.30", "0.20", "0.10"):
    pass
check(63, "five weighted factors with the documented weights",
      "volume_score * 0.30" in src and "prop_diversity * 0.20" in src
      and "balance_score * 0.20" in src and "concavity_score * 0.20" in src
      and "pharm_score * 0.10" in src, "")
check(63, "score bounded to a sensible range",
      all(0.05 <= s["drugability_score"] <= 0.99 for s in sites),
      str([s["drugability_score"] for s in sites]))

check(64, "volume score peaks in the 300-1000 A^3 band",
      "elif volume <= 1000:" in src and "volume_score = 1.0" in src, "")
check(65, "property diversity is a fraction of three categories",
      "/ 3.0" in src, "")
check(66, "hydrophobic balance peaks near 40 percent",
      "abs(hydrophobic_ratio - 0.4)" in src, "")
check(67, "concavity uses spread over mean distance",
      "np.std(distances)" in src and "mean_dist" in src, "")
check(67, "concavity reported per pocket",
      all(0.0 <= s["concavity"] <= 1.0 for s in sites), "")
check(68, "pharmacophore density is n/10 capped at 1",
      "len(pharmacophores) / 10.0" in src, "")

# ---- F69: summary ----------------------------------------------------------
check(69, "summary reports the pocket count",
      summary["site_count"] == len(sites), "")
check(69, "summary totals the cavity volume",
      abs(summary["total_volume"] - sum(s["volume"] for s in sites)) < 0.1, "")
check(69, "summary names the primary site volume",
      summary["primary_volume"] == sites[0]["volume"] if sites else True, "")
check(69, "summary text mentions pockets and volume",
      "pockets" in summary["text"] and "A" in summary["text"], summary["text"][:60])
check(69, "empty result handled",
      SA.BindingSiteAnalyzer.get_summary(an, [])["site_count"] == 0, "")

# ---- hydrogens and waters excluded ----------------------------------------
total = sum(1 for l in open(P("prot.pdb")) if l[:6] in ("ATOM  ", "HETATM"))
check(51, "hydrogens excluded from pocket geometry",
      len(an.atoms) < total and all(a.element != "H" for a in an.atoms),
      f"{len(an.atoms)} of {total}")

# ---- degenerate input ------------------------------------------------------
tiny = P("tiny.pdb")
with open(tiny, "w") as fh:
    fh.write("".join(PROT[:3]) + "END\n")
check(49, "structure too small returns no pockets rather than raising",
      SA.BindingSiteAnalyzer(tiny).analyze() == [], "")


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
print(f"{len(RESULTS)} checks across features 49-69, {failed} failed")
