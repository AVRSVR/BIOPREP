"""
Science audit 2: run the real physics questions against real complexes.

1STP  streptavidin + biotin, a deeply buried ligand
1HSG  HIV-1 protease + MK1 inhibitor, catalytic aspartate dyad
4INS  insulin hexamer + zinc, a metal site
"""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)

import os
import sys
import tempfile
import warnings

import numpy as np

warnings.filterwarnings("ignore")
BIOPREP = _os.path.join(_ROOT, "bioprep")
PDBS = _os.path.join(_HERE, "pdbs")
sys.path.insert(0, BIOPREP)
os.chdir(BIOPREP)

from bioprep.core import io as bio, cleaner, protonator, minimizer  # noqa: E402
from bioprep.core import site_analyzer as SA                         # noqa: E402

TMP = tempfile.mkdtemp(prefix="sci2_")
P = lambda n: os.path.join(TMP, n)      # noqa: E731


def prepare(pdb_id, ligand=None, ph=7.4, keep_water=False):
    src = os.path.join(PDBS, f"{pdb_id}.pdb")
    st = bio.load_pdb(src)
    sel = cleaner.clean_structure(st, remove_water=not keep_water,
                                  protect_ligands=[ligand] if ligand else None)
    clean = P(f"{pdb_id}_clean.pdb")
    bio.save_pdb(st, clean, select=sel, source_pdb=src)
    out = P(f"{pdb_id}_prot.pdb")
    protonator.add_hydrogens(clean, out, ph=ph)
    return out


def atoms_of(path, resname=None, exclude=None):
    rows = [l for l in open(path) if l[:6] in ("ATOM  ", "HETATM")]
    if resname:
        rows = [l for l in rows if l[17:20].strip() == resname]
    if exclude:
        rows = [l for l in rows if l[17:20].strip() != exclude]
    return rows


def xyz(rows):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in rows])


def heavy(rows):
    return [l for l in rows if l[76:78].strip() != "H"]


print("=" * 78)
print("SCIENCE 2a  Does Tier-2 minimisation collapse the pocket onto the ligand?")
print("=" * 78)
print("1STP: biotin sits in a deep streptavidin pocket. Tier 2 deletes BTN,")
print("minimises the protein with the pocket EMPTY, then restores BTN where it")
print("was. If the pocket relaxes inward, the merged complex has clashes that")
print("were never in the crystal structure.")
print()

src = prepare("1STP", ligand="BTN")
lig_in = heavy(atoms_of(src, resname="BTN"))
prot_in = heavy(atoms_of(src, exclude="BTN"))
print(f"prepared: {len(prot_in)} protein heavy atoms, {len(lig_in)} biotin heavy atoms")

stats = minimizer.minimize_structure(src, P("1stp_min.pdb"))
print(f"minimisation status={stats['status']} excluded={stats['excluded_residues']}")
print(f"energy {stats['energy_before_kJ_mol']} -> {stats['energy_after_kJ_mol']}")

lig_out = heavy(atoms_of(P("1stp_min.pdb"), resname="BTN"))
prot_out = heavy(atoms_of(P("1stp_min.pdb"), exclude="BTN"))


def contact_profile(lig, prot):
    d = np.linalg.norm(xyz(lig)[:, None, :] - xyz(prot)[None, :, :], axis=2)
    return d.min(), int((d < 2.2).sum()), int((d < 3.0).sum())


mn_i, c22_i, c30_i = contact_profile(lig_in, prot_in)
mn_o, c22_o, c30_o = contact_profile(lig_out, prot_out)
print()
print(f"{'':22}{'before':>10}{'after':>10}")
print(f"{'closest contact (A)':22}{mn_i:10.2f}{mn_o:10.2f}")
print(f"{'contacts < 2.2 A':22}{c22_i:10d}{c22_o:10d}")
print(f"{'contacts < 3.0 A':22}{c30_i:10d}{c30_o:10d}")

# how far did the pocket-lining atoms move inward?
d_in = np.linalg.norm(xyz(prot_in)[:, None, :] - xyz(lig_in)[None, :, :], axis=2).min(axis=1)
lining = np.where(d_in < 5.0)[0]
key = {}
for i, l in enumerate(prot_in):
    key[(l[21], l[22:26], l[12:16])] = i
moved = []
for l in prot_out:
    k = (l[21], l[22:26], l[12:16])
    if k in key and key[k] in lining:
        i = key[k]
        before = xyz([prot_in[i]])[0]
        after = xyz([l])[0]
        d_before = np.linalg.norm(xyz(lig_in) - before, axis=1).min()
        d_after = np.linalg.norm(xyz(lig_in) - after, axis=1).min()
        moved.append(d_after - d_before)
moved = np.array(moved)
print()
print(f"pocket-lining atoms within 5 A of biotin: {len(moved)}")
if len(moved):
    print(f"   mean change in distance to biotin: {moved.mean():+.3f} A "
          f"(negative = moved toward the ligand)")
    print(f"   atoms that moved INTO the pocket:  {int((moved < -0.1).sum())}")
    print(f"   largest inward move:               {moved.min():+.3f} A")

verdict = []
if mn_o < mn_i - 0.2:
    verdict.append("closest contact tightened")
if c22_o > c22_i:
    verdict.append("new hard clashes")
if len(moved) and moved.mean() < -0.05:
    verdict.append("pocket relaxed inward on average")
print()
print("VERDICT:", "; ".join(verdict) if verdict
      else "pocket did not collapse measurably on this structure")

print()
print("=" * 78)
print("SCIENCE 2b  Is pocket detection rotation-invariant?")
print("=" * 78)
print("The enclosure test casts rays along +-X, +-Y, +-Z only. LIGSITE uses")
print("seven directions including the cubic diagonals precisely so the result")
print("does not depend on how the molecule happens to be oriented. With six")
print("axis-aligned rays, rotating the protein should change the answer -")
print("which would mean the pocket list is an artefact of input orientation.")
print()

base = prepare("1STP", ligand="BTN")
an = SA.BindingSiteAnalyzer(base)
sites0 = an.analyze()
print(f"original orientation: {len(sites0)} pockets, "
      f"volumes {[s['volume'] for s in sites0]}")


def rotate_pdb(path, out, angles):
    ax, ay, az = np.radians(angles)
    Rx = np.array([[1, 0, 0], [0, np.cos(ax), -np.sin(ax)], [0, np.sin(ax), np.cos(ax)]])
    Ry = np.array([[np.cos(ay), 0, np.sin(ay)], [0, 1, 0], [-np.sin(ay), 0, np.cos(ay)]])
    Rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
    R = Rz @ Ry @ Rx
    rows = []
    for l in open(path):
        if l[:6] in ("ATOM  ", "HETATM"):
            v = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            w = R @ v
            rows.append(l[:30] + f"{w[0]:8.3f}{w[1]:8.3f}{w[2]:8.3f}" + l[54:])
        else:
            rows.append(l)
    with open(out, "w") as fh:
        fh.writelines(rows)
    return out


for angles in ((30, 0, 0), (0, 45, 0), (37, 53, 19)):
    rot = rotate_pdb(base, P(f"rot_{angles}.pdb".replace(" ", "")), angles)
    sites = SA.BindingSiteAnalyzer(rot).analyze()
    vols = [s["volume"] for s in sites]
    total0 = sum(s["volume"] for s in sites0)
    total = sum(vols)
    drift = (total - total0) / total0 * 100 if total0 else 0
    print(f"rotated {str(angles):14}: {len(sites)} pockets, total volume "
          f"{total:8.1f} A^3  ({drift:+6.1f}% vs original {total0:.1f})")

print()
print("A rotation is a rigid-body move: the cavity is physically identical.")
print("Any change in pocket count or volume is a defect of the sampling.")

print()
print("=" * 78)
print("SCIENCE 2c  Structural water: proximity or actual hydrogen bonding?")
print("=" * 78)
print("The rule is 'within 4.0 A of ANY protein atom'. Carbon counts. A water")
print("packed against a hydrophobic surface is not a structural water - it has")
print("no hydrogen bond holding it. A defensible criterion is a polar contact")
print("(N, O, S) inside about 3.5 A.")
print()

src_pdb = os.path.join(PDBS, "1STP.pdb")
st = bio.load_pdb(src_pdb)
sel = cleaner.clean_structure(st, remove_water=True, keep_structural_waters=True)
kept = sel.structural_waters
print(f"waters the current rule keeps: {len(kept)}")

model = next(iter(st))
prot_atoms, waters = [], []
for ch in model:
    for res in ch:
        if res.id[0] == " ":
            prot_atoms.extend(res.get_atoms())
        elif res.get_resname().strip() in ("HOH", "WAT"):
            waters.append((ch.id, res))

pa = np.array([a.get_coord() for a in prot_atoms])
polar_mask = np.array([(a.element or a.get_name()[0]).strip().upper() in ("N", "O", "S")
                       for a in prot_atoms])
pa_polar = pa[polar_mask]
print(f"protein atoms: {len(pa)} total, {int(polar_mask.sum())} polar (N/O/S)")

any_atom_4 = 0
polar_35 = 0
kept_but_no_polar = 0
# Measure the set the code ACTUALLY keeps, not the old rule this script
# originally hardcoded - that comparison stopped being meaningful once the
# criterion changed.
for chain_id, res in waters:
    o = [a for a in res.get_atoms()
         if (a.element or "").strip().upper() == "O"
         or a.get_name().strip().upper().startswith("O")]
    if not o:
        continue
    c = o[0].get_coord()
    d_any = np.linalg.norm(pa - c, axis=1).min()
    d_pol = np.linalg.norm(pa_polar - c, axis=1).min() if len(pa_polar) else 99
    if d_any <= 4.0:
        any_atom_4 += 1
    if d_pol <= 3.5:
        polar_35 += 1
    if (chain_id, res.id) in kept and d_pol > 3.5:
        kept_but_no_polar += 1

print(f"waters within 4.0 A of ANY atom      : {any_atom_4}")
print(f"waters within 3.5 A of a POLAR atom  : {polar_35}")
print(f"the OLD rule (4.0 A to any atom) would keep : {any_atom_4}")
print(f"actually kept with NO polar contact <= 3.5 A : {kept_but_no_polar}")
print()
print("VERDICT:", "current rule keeps waters with no hydrogen bond"
      if kept_but_no_polar else "every kept water has a polar contact")
