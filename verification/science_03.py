"""Science audit 3: pKa handling, HIS tautomers, PDBQT charges, convergence."""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)

import os
import subprocess
import sys
import tempfile
import warnings
from collections import Counter

import numpy as np

warnings.filterwarnings("ignore")
BIOPREP = _os.path.join(_ROOT, "bioprep")
PDBS = _os.path.join(_HERE, "pdbs")
sys.path.insert(0, BIOPREP)
os.chdir(BIOPREP)

from bioprep.core import io as bio, cleaner, protonator, minimizer, exporter  # noqa: E402

TMP = tempfile.mkdtemp(prefix="sci3_")
P = lambda n: os.path.join(TMP, n)      # noqa: E731


def prep(pid, ph=7.4, remove_hets=True):
    src = os.path.join(PDBS, f"{pid}.pdb")
    st = bio.load_pdb(src)
    bio.save_pdb(st, P(f"{pid}_c.pdb"),
                 select=cleaner.clean_structure(
                     st, remove_water=True,
                     remove_heteroatoms=["ALL"] if remove_hets else []),
                 source_pdb=src)
    out = P(f"{pid}_ph{ph}.pdb")
    protonator.add_hydrogens(P(f"{pid}_c.pdb"), out, ph=ph)
    return out


def rows(path):
    return [l for l in open(path) if l[:6] in ("ATOM  ", "HETATM")]


print("=" * 78)
print("SCIENCE 3a  Does the pH parameter produce correct titration behaviour?")
print("=" * 78)
print("1HSG has 12 LYS, 8 ASP, 8 GLU, 2 HIS - enough to see real titration.")
print("Model-compound pKa: ASP 3.9, GLU 4.3, HIS 6.0, LYS 10.5, TYR 10.1.")
print()

# atom names that only exist in the protonated form of each side chain
MARKERS = {
    "ASP": ("HD2",), "GLU": ("HE2",),
    "LYS": ("HZ1", "HZ2", "HZ3"),
    "TYR": ("HH",),
    "HIS": ("HD1", "HE2"),
    "CYS": ("HG",),
}

print(f"{'pH':>5} | {'ASP-COOH':>9} {'GLU-COOH':>9} {'LYS-NH3':>8} "
      f"{'TYR-OH':>7} {'HIS-H':>6}")
print("-" * 58)
for ph in (2.0, 4.0, 6.0, 7.4, 9.0, 11.0, 13.0):
    out = prep("1HSG", ph=ph)
    r = rows(out)
    counts = {}
    for res, names in MARKERS.items():
        residues = set()
        for l in r:
            if l[17:20].strip() == res and l[12:16].strip() in names:
                residues.add((l[21], l[22:26]))
        counts[res] = len(residues)
    print(f"{ph:>5} | {counts['ASP']:>9} {counts['GLU']:>9} {counts['LYS']:>8} "
          f"{counts['TYR']:>7} {counts['HIS']:>6}")

print()
print("Expected: ASP/GLU protonated only below ~4, LYS protonated below ~10.5,")
print("TYR above ~10 loses HH, HIS protonated below ~6.")

print()
print("=" * 78)
print("SCIENCE 3b  Which histidine tautomer is produced?")
print("=" * 78)
print("Neutral HIS has one ring NH: HD1 only (HID) or HE2 only (HIE). Which one")
print("is chosen changes the hydrogen-bond pattern in a binding site, and the")
print("choice should follow the local environment, not a fixed default.")
print()

out = prep("1HSG", ph=7.4)
r = rows(out)
his = {}
for l in r:
    if l[17:20].strip() == "HIS":
        key = (l[21], l[22:26].strip())
        his.setdefault(key, set()).add(l[12:16].strip())
for key, names in sorted(his.items()):
    hd1, he2 = "HD1" in names, "HE2" in names
    if hd1 and he2:
        state = "HIP (both, +1 charge)"
    elif hd1:
        state = "HID (delta tautomer)"
    elif he2:
        state = "HIE (epsilon tautomer)"
    else:
        state = "neither - unusual"
    print(f"   HIS {key[0]}{key[1]:>4}: {state}")
tauts = Counter(
    ("HIP" if {"HD1", "HE2"} <= n else "HID" if "HD1" in n else
     "HIE" if "HE2" in n else "?") for n in his.values())
print(f"   tautomer distribution: {dict(tauts)}")

print()
print("=" * 78)
print("SCIENCE 3c  PDBQT partial charges: are they chemically sensible?")
print("=" * 78)

import shutil                                                   # noqa: E402
if not shutil.which("obabel"):
    print("obabel not on PATH - skipped")
else:
    src = prep("1STP", ph=7.4, remove_hets=False)
    ok, pdbqt = exporter.export_structure(src, P("rec.pdb"), "vina")
    print(f"export ok={ok} -> {os.path.basename(str(pdbqt))}")
    if ok:
        charges, by_element = [], {}
        for l in open(pdbqt):
            if l[:6] in ("ATOM  ", "HETATM") and len(l) > 76:
                try:
                    q = float(l[66:76])
                except ValueError:
                    continue
                charges.append(q)
                el = l[77:79].strip() or l[12:16].strip()[:1]
                by_element.setdefault(el, []).append(q)
        charges = np.array(charges)
        print(f"   atoms with a charge: {len(charges)}")
        print(f"   total charge       : {charges.sum():+.3f} e")
        print(f"   range              : {charges.min():+.3f} to {charges.max():+.3f} e")
        print()
        print("   mean charge by atom type (Gasteiger):")
        for el in sorted(by_element, key=lambda e: -len(by_element[e]))[:8]:
            v = np.array(by_element[el])
            print(f"      {el:3} n={len(v):4}  mean {v.mean():+.3f}  "
                  f"range {v.min():+.3f}..{v.max():+.3f}")
        print()
        print("   Sanity: N and O should be net negative, carbon near zero,")
        print("   polar hydrogens positive. A total far from an integer would")
        print("   indicate the charge model failed.")

print()
print("=" * 78)
print("SCIENCE 3d  Is 'converged' a real convergence test?")
print("=" * 78)
print("The flag is derived from the energy change, not from the gradient.")
print("L-BFGS is told to stop at a force tolerance; whether it reached that")
print("tolerance or simply hit the iteration cap is not currently reported.")
print()

src = prep("1STP", ph=7.4, remove_hets=False)
for iters in (10, 100, 1000):
    minimizer.MAX_ITERATIONS = iters
    st = minimizer.minimize_structure(src, P(f"m{iters}.pdb"))
    print(f"   maxIterations={iters:5}: status={st['status']:8} "
          f"E {st['energy_before_kJ_mol']} -> {st['energy_after_kJ_mol']}  "
          f"converged={st['converged']}")
minimizer.MAX_ITERATIONS = 1000
print()
print("   If a 10-iteration run also reports converged=True, the flag is")
print("   measuring 'energy went down', not 'the minimiser finished'.")
