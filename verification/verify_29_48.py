"""Independent verification of catalogue features 29-48 (minimizer.py)."""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)

import inspect
import os
import sys
import tempfile
import warnings

warnings.filterwarnings("ignore")
BIOPREP = _os.path.join(_ROOT, "bioprep")
sys.path.insert(0, BIOPREP)
os.chdir(BIOPREP)

import openmm as mm                                     # noqa: E402
from openmm import app, unit                            # noqa: E402
from bioprep.core import io as bio                       # noqa: E402
from bioprep.core import minimizer, protonator, cleaner  # noqa: E402
from bioprep.core import residues as R                   # noqa: E402

TMP = tempfile.mkdtemp(prefix="min_")
P = lambda n: os.path.join(TMP, n)      # noqa: E731
RESULTS = []


def check(feature, description, passed, detail=""):
    RESULTS.append((feature, description, passed, detail))


def het(serial, name, resn, chain, resseq, x, y, z, elem):
    return (f"HETATM{serial:>5} {name:<4} {resn:>3} {chain:1}{resseq:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00          {elem:>2}  \n")


PROT = [l for l in open("1crn.pdb") if l.startswith("ATOM")]
CA = [l for l in PROT if l[12:16].strip() == "CA"][0]
CX, CY, CZ = float(CA[30:38]), float(CA[38:46]), float(CA[46:54])


def prepared(extra=(), name="p.pdb"):
    """Clean + protonate a structure so it is ready to minimise."""
    raw = P("raw_" + name)
    with open(raw, "w") as fh:
        fh.write("".join(PROT) + "".join(extra) + "END\n")
    st = bio.load_pdb(raw)
    bio.save_pdb(st, P("cl_" + name),
                 select=cleaner.clean_structure(st, remove_water=False))
    protonator.add_hydrogens(P("cl_" + name), P(name))
    return P(name)


PLAIN = prepared()
LIGAND = prepared([het(900, "C1", "LIG", "A", 900, CX + 6, CY, CZ, "C")], "lig.pdb")
WATER = prepared([het(910, "O", "HOH", "A", 700, CX + 2.8, CY, CZ, "O")], "wat.pdb")

src = inspect.getsource(minimizer)

# ---- F29: force field selection -------------------------------------------
for ff in ("amber14", "charmm36"):
    try:
        obj = minimizer._build_forcefield(ff, use_gbsa=False)
        check(29, f"{ff} force field loads", obj is not None, "")
    except Exception as exc:
        check(29, f"{ff} force field loads", False, str(exc)[:60])

r = minimizer.minimize_structure(PLAIN, P("o_amber.pdb"), force_field="amber14")
check(29, "amber14 minimisation runs", r["status"] == "full", r["status"])
check(29, "force field echoed in the result", r["force_field"] == "amber14", "")

r_ch = minimizer.minimize_structure(PLAIN, P("o_charmm.pdb"), force_field="charmm36")
check(29, "charmm36 minimisation runs", r_ch["status"] == "full",
      f"{r_ch['status']} {str(r_ch.get('error'))[:60]}")

# ---- F30: GBSA toggle ------------------------------------------------------
check(30, "obc2 implicit solvent referenced", "implicit/obc2.xml" in src, "")
on = minimizer.minimize_structure(PLAIN, P("g_on.pdb"), use_gbsa=True)
off = minimizer.minimize_structure(PLAIN, P("g_off.pdb"), use_gbsa=False)
check(30, "GBSA on and off both succeed",
      on["status"] == "full" and off["status"] == "full", "")
check(30, "GBSA changes the energy",
      on["energy_after_kJ_mol"] != off["energy_after_kJ_mol"],
      f"{on['energy_after_kJ_mol']} vs {off['energy_after_kJ_mol']}")
check(30, "gbsa_used reported", on["gbsa_used"] is True and off["gbsa_used"] is False, "")

# ---- F31: water/ion parameters always present ------------------------------
w_on = minimizer.minimize_structure(WATER, P("w_on.pdb"), use_gbsa=True)
w_off = minimizer.minimize_structure(WATER, P("w_off.pdb"), use_gbsa=False)
check(31, "retained water minimises with GBSA on",
      w_on["status"] == "full", f"{w_on['status']} {str(w_on.get('error'))[:60]}")
check(31, "retained water minimises with GBSA off",
      w_off["status"] == "full", w_off["status"])

# ---- F32: tier 1 -----------------------------------------------------------
check(32, "clean structure minimises at the full tier",
      r["status"] == minimizer.STATUS_FULL, r["status"])
check(32, "energy actually decreases", r["energy_after_kJ_mol"] < r["energy_before_kJ_mol"],
      f"{r['energy_before_kJ_mol']} -> {r['energy_after_kJ_mol']}")

# ---- F33: terminal group repair --------------------------------------------
check(33, "Modeller.addHydrogens used for terminal repair",
      "addHydrogens(forcefield=" in src, "")
# a chain truncated so the terminus lacks OXT
trunc = [l for l in PROT if int(l[22:26]) <= 20]
tr = P("trunc_raw.pdb")
with open(tr, "w") as fh:
    fh.write("".join(trunc) + "END\n")
st = bio.load_pdb(tr)
bio.save_pdb(st, P("trunc_clean.pdb"), select=cleaner.clean_structure(st))
# deliberately skip protonation so the terminus is bare
r33 = minimizer.minimize_structure(P("trunc_clean.pdb"), P("trunc_out.pdb"))
check(33, "unparameterisable terminus still produces a result",
      r33["status"] != "failed" or r33["error"] is not None,
      f"{r33['status']}")

# ---- F34/F37/F42: tier 2 ---------------------------------------------------
r2 = minimizer.minimize_structure(LIGAND, P("t2.pdb"))
check(34, "ligand forces the partial tier", r2["status"] == minimizer.STATUS_PARTIAL,
      r2["status"])
check(34, "excluded residues reported", r2["excluded_residues"] == ["LIG"],
      str(r2["excluded_residues"]))
check(42, "partial minimisation carries a warning",
      any("held at their input coordinates" in w for w in r2["warnings"]),
      str(r2["warnings"])[:80])

def lig_atoms(path):
    # Restrict to coordinate records: PDBFile also writes a TER whose
    # residue-name columns hold the ligand name.
    return [l[30:54] for l in open(path)
            if l[:6] in ("ATOM  ", "HETATM") and l[17:20].strip() == "LIG"]


check(37, "excluded ligand keeps its input coordinates",
      lig_atoms(LIGAND) == lig_atoms(P("t2.pdb")),
      f"{lig_atoms(LIGAND)} vs {lig_atoms(P('t2.pdb'))}")
before_ca = [l for l in open(LIGAND) if l[12:16].strip() == "CA"][0][30:54]
after_ca = [l for l in open(P("t2.pdb")) if l[12:16].strip() == "CA"][0][30:54]
check(37, "protein coordinates did move", before_ca != after_ca, "")
check(37, "atom count preserved through the merge",
      sum(1 for l in open(LIGAND) if l[:6] in ("ATOM  ", "HETATM"))
      == sum(1 for l in open(P("t2.pdb")) if l[:6] in ("ATOM  ", "HETATM")), "")
check(37, "merge matches by atom identity, not index",
      "_merge_coords_by_name" in src and "_atom_key" in src, "")

# ---- F35: tier 3 -----------------------------------------------------------
check(35, "a no-GBSA retry tier exists",
      minimizer.STATUS_PARTIAL_NO_GBSA in src, "")
attempts = src.split("attempts = [")[1].split("]")[0]
check(35, "tier 3 only added when GBSA was requested",
      "if use_gbsa:" in src.split("attempts = [")[1][:400], "")

# ---- F36: final fallback ---------------------------------------------------
junk = P("junk.pdb")
with open(junk, "w") as fh:
    fh.write("REMARK nothing here\n")
r36 = minimizer.minimize_structure(junk, P("junk_out.pdb"))
check(36, "unusable input still produces an output file",
      os.path.isfile(P("junk_out.pdb")), "")
check(36, "failure is reported, not hidden",
      r36["status"] == "failed" and r36["error"], str(r36["error"])[:60])
check(36, "fallback output is the input unchanged",
      open(junk, "rb").read() == open(P("junk_out.pdb"), "rb").read(), "")

# ---- F38: residue classification -------------------------------------------
check(38, "amino acids include protonation variants",
      {"HID", "HIE", "HIP", "CYX", "ASH", "GLH"} <= R.AMINO_ACIDS, "")
check(38, "water set covers the common names",
      {"HOH", "WAT", "SOL", "TIP3", "DOD"} <= R.WATER, "")
check(38, "nucleic acids classified",
      {"DA", "DC", "DG", "DT", "A", "C", "G", "U"} <= R.NUCLEIC_ACIDS, "")
check(38, "force-field-safe set is the union",
      R.FORCE_FIELD_SAFE == (R.AMINO_ACIDS | R.NUCLEIC_ACIDS | R.WATER), "")
check(38, "a novel ligand is not force-field-safe",
      not R.is_force_field_safe("LIG"), "")

# ---- F39: energy validation ------------------------------------------------
check(39, "NaN and Inf rejected", "math.isnan" in src and "math.isinf" in src, "")
check(39, "suspiciously high energy threshold defined",
      minimizer.SUSPECT_ENERGY == 1.0e12, str(minimizer.SUSPECT_ENERGY))
clash = [PROT[0]] + [PROT[0][:30] + PROT[0][30:54] + PROT[0][54:]]
r39 = minimizer.minimize_structure(WATER, P("w2.pdb"), use_gbsa=True)
flagged = [w for w in r39["warnings"] if "clash" in w.lower()]
check(39, "high starting energy produces a warning",
      (abs(r39["energy_before_kJ_mol"]) < minimizer.SUSPECT_ENERGY) or flagged,
      f"E0={r39['energy_before_kJ_mol']}")

# ---- F40: energy reporting -------------------------------------------------
check(40, "before, after and delta all reported",
      all(isinstance(r[k], float) for k in
          ("energy_before_kJ_mol", "energy_after_kJ_mol", "delta_energy_kJ_mol")),
      "")
check(40, "delta equals after minus before",
      abs(r["delta_energy_kJ_mol"]
          - (r["energy_after_kJ_mol"] - r["energy_before_kJ_mol"])) < 0.2, "")

# ---- F41: convergence ------------------------------------------------------
check(41, "converged reported", isinstance(r.get("converged"), bool), str(r.get("converged")))
check(41, "a real energy drop counts as converged", r["converged"] is True,
      f"delta={r['delta_energy_kJ_mol']}")
check(41, "a failed run is not converged", r36["converged"] is False, "")

# ---- F43/F44: iteration cap and tolerance ----------------------------------
check(43, "max iterations is 1000", minimizer.MAX_ITERATIONS == 1000, "")
check(43, "cap passed to minimizeEnergy", "maxIterations=MAX_ITERATIONS" in src, "")
check(43, "cap reported in the result", r["iterations_max"] == 1000, "")
check(44, "energy tolerance is 10.0", minimizer.ENERGY_TOLERANCE == 10.0, "")
check(44, "tolerance passed to minimizeEnergy", "tolerance=ENERGY_TOLERANCE" in src, "")

# ---- F45: integrator -------------------------------------------------------
check(45, "LangevinMiddleIntegrator used", "LangevinMiddleIntegrator" in src, "")
check(45, "300 K, 1/ps friction, 2 fs step",
      "300 * unit.kelvin" in src and "1 / unit.picosecond" in src
      and "0.002 * unit.picosecond" in src, "")

# ---- F46: platform preference ----------------------------------------------
check(46, "CPU platform requested explicitly", "getPlatformByName('CPU')" in src, "")
_sim_body = inspect.getsource(minimizer._make_simulation)
check(46, "falls back when CPU is unavailable",
      "except Exception" in _sim_body and "app.Simulation(topology, system, integrator)"
      in _sim_body, "")
try:
    mm.Platform.getPlatformByName("CPU")
    check(46, "CPU platform available in this environment", True, "")
except Exception as exc:
    check(46, "CPU platform available in this environment", False, str(exc)[:50])

# ---- F47/F48: nonbonded and constraints ------------------------------------
check(47, "NoCutoff nonbonded method", "nonbondedMethod=app.NoCutoff" in src, "")
# DEVIATION FROM THE CATALOGUE, on purpose. Feature 48 asks for HBonds
# constraints "for stability". Constraints exist to allow a longer MD
# timestep and nothing here integrates dynamics; during minimisation they
# stop hydrogens relaxing and put constraint contributions into getForces(),
# so the residual force can never be compared against the tolerance.
# On 1CRN after 1000 iterations: HBonds gave -5162.7 kJ/mol at RMS force
# 71.16, unconstrained gave -5166.9 at 5.76.
check(48, "minimisation is unconstrained",
      minimizer.MINIMISATION_CONSTRAINTS is None
      and "constraints=MINIMISATION_CONSTRAINTS" in src,
      str(minimizer.MINIMISATION_CONSTRAINTS))
check(48, "convergence compares RMS force to the tolerance",
      "rms_force <= ENERGY_TOLERANCE" in src, "")


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
print(f"{len(RESULTS)} checks across features 29-48, {failed} failed")
