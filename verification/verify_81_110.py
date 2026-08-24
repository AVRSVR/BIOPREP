"""Independent verification of catalogue features 81-110 (app.py / HTTP layer)."""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)

import base64
import io as _io
import json
import os
import sys
import tempfile
import warnings
import zipfile

warnings.filterwarnings("ignore")
BIOPREP = _os.path.join(_ROOT, "bioprep")
sys.path.insert(0, BIOPREP)
os.chdir(BIOPREP)

from bioprep.core import io as bio, cleaner                  # noqa: E402
from bioprep import webapp                                    # noqa: E402
import app as shim                                            # noqa: E402

TMP = tempfile.mkdtemp(prefix="api_")
P = lambda n: os.path.join(TMP, n)      # noqa: E731
RESULTS = []
client = webapp.app.test_client()


def check(feature, description, passed, detail=""):
    RESULTS.append((feature, description, passed, detail))


def het(serial, name, resn, chain, resseq, x, y, z, elem):
    return (f"HETATM{serial:>5} {name:<4} {resn:>3} {chain:1}{resseq:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00          {elem:>2}  \n")


PROT = [l for l in open("1crn.pdb") if l.startswith("ATOM")]
SEQRES = [l for l in open("1crn.pdb") if l.startswith("SEQRES")]
CA = [l for l in PROT if l[12:16].strip() == "CA"][0]
CX, CY, CZ = float(CA[30:38]), float(CA[38:46]), float(CA[46:54])

# a rich structure: two chains, a ligand, a sulphate, a near water, a gap
CHAIN_B = [l[:21] + "B" + l[22:] for l in PROT]
RICH = ("".join(PROT) + "".join(CHAIN_B)
        + het(9000, "C1", "LIG", "A", 900, CX + 6, CY, CZ, "C")
        + het(9001, "S", "SO4", "A", 901, CX + 12, CY, CZ, "S")
        + het(9002, "O", "HOH", "A", 700, CX + 2.8, CY, CZ, "O")
        + "END\n").encode()

PLAIN = ("".join(PROT) + "END\n").encode()
GAPPED = ("".join(SEQRES)
          + "".join(l for l in PROT if not (20 <= int(l[22:26]) <= 24))
          + "END\n").encode()


def upload(data, name="s.pdb"):
    return (_io.BytesIO(data), name)


def post(route, data, **kw):
    return client.post(route, data=data, content_type="multipart/form-data", **kw)


def decode(b64):
    return base64.b64decode(b64).decode("utf-8", "replace")


# ---- F81/F82: server and upload limit --------------------------------------
check(81, "Flask application object exists", webapp.app is not None, "")
check(81, "the bioprep/app.py shim still exposes the app and its config",
      shim.app is webapp.app and shim.RESULTS_DIR == webapp.RESULTS_DIR, "")
check(81, "index route registered",
      any(r.rule == "/" for r in webapp.app.url_map.iter_rules()), "")
import inspect as _inspect                                     # noqa: E402
# The application moved into the package so an installed copy can serve
# it; bioprep/app.py is only a source-checkout shim now.
src = _inspect.getsource(webapp)
check(81, "served on port 5000", "port=5000" in src, "")
check(81, "debug disabled", "debug=False" in src, "")
check(82, "upload limit is 200 MB",
      webapp.app.config["MAX_CONTENT_LENGTH"] == 200 * 1024 * 1024,
      str(webapp.app.config["MAX_CONTENT_LENGTH"]))

routes = {r.rule for r in webapp.app.url_map.iter_rules()}
for feature, rule in ((83, "/api/analyze"), (84, "/api/process"),
                      (99, "/api/batch"), (100, "/api/high-throughput"),
                      (103, "/api/history"), (106, "/api/templates"),
                      (108, "/api/download-report"), (109, "/api/analyze-site")):
    check(feature, f"{rule} registered", rule in routes, str(sorted(routes))[:60])

# ---- F83: analyze ----------------------------------------------------------
r = post("/api/analyze", {"file": upload(RICH)})
check(83, "analyze returns 200", r.status_code == 200, str(r.status_code))
body = r.get_json()
meta = body["metadata"]
check(83, "reports chains", meta["chains"] == ["A", "B"], str(meta["chains"]))
check(83, "reports heteroatoms", set(meta["heteroatoms"]) == {"LIG", "SO4"},
      str(meta["heteroatoms"]))
check(83, "reports water count", meta["water_count"] == 1, str(meta["water_count"]))
check(83, "reports gaps key", "sequence_gaps" in meta, "")
check(83, "returns a session id", bool(body.get("session_id")), "")
ANALYZE_SESSION = body["session_id"]

r = post("/api/analyze", {"file": upload(b"not a pdb\n")})
check(83, "bad upload gives 400 with a reason",
      r.status_code == 400 and "No atoms" in r.get_json().get("error", ""),
      str(r.status_code))
r = client.post("/api/analyze", data={}, content_type="multipart/form-data")
check(83, "missing file gives 400", r.status_code == 400, str(r.status_code))

# ---- F84: process ----------------------------------------------------------
r = post("/api/process", {"file": upload(RICH), "ph": "7.4"})
check(84, "process returns 200", r.status_code == 200, str(r.status_code))
body = r.get_json()
for key in ("report", "report_text", "viewer_pdb_b64", "pdb_b64", "session_id"):
    check(84, f"response contains {key}", key in body, "")
BASE_SESSION = body["session_id"]

# ---- F85: pH ---------------------------------------------------------------
def h_count(b64):
    return sum(1 for l in decode(b64).splitlines()
               if l[:6] in ("ATOM  ", "HETATM") and l[76:78].strip() == "H")


low = post("/api/process", {"file": upload(PLAIN), "ph": "1.0"}).get_json()
high = post("/api/process", {"file": upload(PLAIN), "ph": "13.0"}).get_json()
check(85, "pH reaches the pipeline",
      h_count(low["viewer_pdb_b64"]) > h_count(high["viewer_pdb_b64"]),
      f"{h_count(low['viewer_pdb_b64'])} vs {h_count(high['viewer_pdb_b64'])}")
check(85, "pH echoed in the report", low["report"]["protonation_ph"] == 1.0, "")
bad_ph = post("/api/process", {"file": upload(PLAIN), "ph": "not-a-number"}).get_json()
check(85, "invalid pH falls back to 7.4",
      bad_ph["report"]["protonation_ph"] == 7.4, "")

# ---- F86: chain selection --------------------------------------------------
r = post("/api/process", {"file": upload(RICH), "chains": json.dumps(["A"])}).get_json()
chains_out = {l[21] for l in decode(r["viewer_pdb_b64"]).splitlines()
              if l[:6] in ("ATOM  ", "HETATM")}
check(86, "only the requested chain is kept", chains_out == {"A"}, str(chains_out))
check(86, "report records the retained chain",
      r["report"]["chains"]["retained"] == ["A"], "")

# ---- F87/F88: heteroatom removal and ligand protection ---------------------
r = post("/api/process", {"file": upload(RICH),
                          "remove_heteros": json.dumps(["SO4"])}).get_json()
names = {l[17:20].strip() for l in decode(r["viewer_pdb_b64"]).splitlines()
         if l.startswith("HETATM")}
check(87, "named heteroatom removed", "SO4" not in names, str(names))
check(87, "unlisted ligand kept", "LIG" in names, str(names))
check(87, "report lists what was removed",
      "SO4" in r["report"]["heteroatoms"]["removed"], "")

r = post("/api/process", {"file": upload(RICH),
                          "remove_heteros": json.dumps(["ALL"]),
                          "protect_ligands": json.dumps(["LIG"])}).get_json()
names = {l[17:20].strip() for l in decode(r["viewer_pdb_b64"]).splitlines()
         if l.startswith("HETATM")}
check(88, "protected ligand survives 'ALL'", "LIG" in names, str(names))
check(88, "unprotected heteroatoms removed by 'ALL'", "SO4" not in names, str(names))

# ---- F89: structural water toggle ------------------------------------------
off = post("/api/process", {"file": upload(RICH),
                            "keep_structural_waters": "false"}).get_json()
on = post("/api/process", {"file": upload(RICH),
                           "keep_structural_waters": "true"}).get_json()


def water_lines(b64):
    return [l for l in decode(b64).splitlines()
            if l[:6] in ("ATOM  ", "HETATM") and l[17:20].strip() == "HOH"]


check(89, "structural water kept only when the toggle is on",
      len(water_lines(off["viewer_pdb_b64"])) == 0
      and len(water_lines(on["viewer_pdb_b64"])) > 0,
      f"off={len(water_lines(off['viewer_pdb_b64']))} "
      f"on={len(water_lines(on['viewer_pdb_b64']))}")
check(89, "report counts retained waters",
      on["report"]["water_molecules_retained"] >= 1, "")

# ---- F90: loop reconstruction ----------------------------------------------
off = post("/api/process", {"file": upload(GAPPED),
                            "reconstruct_loops": "false"}).get_json()
on = post("/api/process", {"file": upload(GAPPED),
                           "reconstruct_loops": "true"}).get_json()
check(90, "loops rebuilt only when the toggle is on",
      on["report"]["protonation"]["loops_reconstructed"] == 5
      and off["report"]["protonation"]["loops_reconstructed"] == 0,
      f"on={on['report']['protonation']['loops_reconstructed']} "
      f"off={off['report']['protonation']['loops_reconstructed']}")

# ---- F91: missing atoms ----------------------------------------------------
stripped = ("".join(l for l in PROT
                    if not (l[22:26].strip() == "2"
                            and l[12:16].strip() in ("CG1", "CG2", "OG1")))
            + "END\n").encode()


def heavy_res2(b64):
    return {l[12:16].strip() for l in decode(b64).splitlines()
            if l.startswith("ATOM") and l[22:26].strip() == "2"
            and l[76:78].strip() != "H"}


off = post("/api/process", {"file": upload(stripped),
                            "add_missing_atoms": "false"}).get_json()
on = post("/api/process", {"file": upload(stripped),
                           "add_missing_atoms": "true"}).get_json()
check(91, "missing side-chain atoms added only when the toggle is on",
      len(heavy_res2(on["viewer_pdb_b64"])) > len(heavy_res2(off["viewer_pdb_b64"])),
      f"off={sorted(heavy_res2(off['viewer_pdb_b64']))} "
      f"on={sorted(heavy_res2(on['viewer_pdb_b64']))}")

# ---- F92/F93/F94: minimisation toggles -------------------------------------
r_off = post("/api/process", {"file": upload(PLAIN),
                              "run_minimization": "false"}).get_json()
r_on = post("/api/process", {"file": upload(PLAIN),
                             "run_minimization": "true"}).get_json()
check(92, "minimisation block absent unless requested",
      "energy_minimization" not in r_off["report"], "")
check(92, "minimisation block present when requested",
      "energy_minimization" in r_on["report"], "")
check(92, "minimisation actually ran",
      r_on["report"]["energy_minimization"]["status"] == "full",
      str(r_on["report"]["energy_minimization"]["status"]))

g_off = post("/api/process", {"file": upload(PLAIN), "run_minimization": "true",
                              "use_gbsa": "false"}).get_json()
check(93, "GBSA toggle reaches the minimiser",
      r_on["report"]["energy_minimization"]["gbsa_used"] is True
      and g_off["report"]["energy_minimization"]["gbsa_used"] is False, "")

ff = post("/api/process", {"file": upload(PLAIN), "run_minimization": "true",
                           "force_field": "charmm36"}).get_json()
check(94, "force field selection honoured",
      ff["report"]["energy_minimization"]["force_field"] == "charmm36",
      str(ff["report"]["energy_minimization"]["force_field"]))
bad_ff = post("/api/process", {"file": upload(PLAIN), "run_minimization": "true",
                               "force_field": "bogus"}).get_json()
check(94, "unknown force field falls back to amber14",
      bad_ff["report"]["energy_minimization"]["force_field"] == "amber14", "")

# ---- F95/F97/F98: docking export, base64, viewer vs download ---------------
r = post("/api/process", {"file": upload(PLAIN),
                          "docking_target": "vina"}).get_json()
check(95, "docking target recorded",
      r["report"].get("docking_target") == "vina", str(r["report"].get("docking_target")))
check(98, "download filename is the pdbqt", r["filename"].endswith(".pdbqt"),
      r["filename"])
check(98, "viewer payload is still a PDB",
      decode(r["viewer_pdb_b64"]).lstrip().startswith(("REMARK", "ATOM", "CRYST")), "")
check(98, "viewer and download payloads differ",
      r["viewer_pdb_b64"] != r["pdb_b64"], "")

plain = post("/api/process", {"file": upload(PLAIN)}).get_json()
check(98, "without export both payloads are the same object",
      plain["viewer_pdb_b64"] == plain["pdb_b64"], "")
check(97, "payload is valid base64 of a PDB",
      "ATOM" in decode(plain["viewer_pdb_b64"]), "")

# ---- F96/F105: sessions and persistence ------------------------------------
import uuid as _uuid
try:
    _uuid.UUID(BASE_SESSION)
    check(96, "session id is a UUID", True, "")
except ValueError:
    check(96, "session id is a UUID", False, BASE_SESSION)
check(96, "sessions are distinct per request", BASE_SESSION != ANALYZE_SESSION, "")
stored = os.path.join(webapp.RESULTS_DIR, f"{BASE_SESSION}.pdb")
check(105, "processed structure persisted to disk", os.path.isfile(stored), stored)
check(105, "session store is bounded",
      webapp.SESSION_LIMIT == 200 and hasattr(webapp.sessions, "_limit"), "")

# ---- F99/F102: batch -------------------------------------------------------
r = client.post("/api/batch", data={
    "files": [upload(RICH, "a.pdb"), upload(PLAIN, "b.pdb"),
              upload(b"junk", "c.pdb")],
    "config": json.dumps({"remove_heteros": ["ALL"], "protect_ligands": ["LIG"],
                          "ph": 7.0}),
}, content_type="multipart/form-data")
check(99, "batch returns a zip", r.status_code == 200
      and r.headers["Content-Type"].startswith("application/zip"), str(r.status_code))
zf = zipfile.ZipFile(_io.BytesIO(r.data))
names = zf.namelist()
check(99, "batch prepared both good files",
      len([n for n in names if n.endswith(".pdb")]) == 2, str(names))
check(102, "batch report included", "processing_report.txt" in names, str(names))
report_txt = zf.read("processing_report.txt").decode()
for token in ("BATCH PROCESSING REPORT", "Succeeded", "Failed", "SETTINGS USED"):
    check(102, f"batch report mentions {token}", token in report_txt, "")
check(102, "batch report names the file that failed", "c.pdb" in report_txt, "")
check(99, "batch honours config",
      "LIG" in zf.read([n for n in names if n.startswith("a")][0]).decode(), "")

# ---- F100/F101: high throughput --------------------------------------------
buf = _io.BytesIO()
with zipfile.ZipFile(buf, "w") as z:
    z.writestr("a.pdb", RICH)
    z.writestr("nested/a.pdb", PLAIN)        # same basename, different folder
    z.writestr("readme.txt", b"ignore me")
buf.seek(0)
r = client.post("/api/high-throughput", data={
    "file": (buf, "in.zip"),
    "config": json.dumps({"chains": ["A"], "remove_heteros": ["ALL"],
                          "protect_ligands": ["LIG"]}),
}, content_type="multipart/form-data")
check(100, "high-throughput returns a zip", r.status_code == 200, str(r.status_code))
zf = zipfile.ZipFile(_io.BytesIO(r.data))
names = zf.namelist()
pdbs = [n for n in names if n.endswith(".pdb")]
check(100, "every pdb in the archive processed", len(pdbs) == 2, str(names))
check(100, "processing log included", "processing_logs.txt" in names, str(names))
check(101, "colliding basenames given distinct output names",
      len(set(pdbs)) == 2, str(pdbs))
check(100, "high-throughput honours config",
      any("LIG" in zf.read(n).decode() for n in pdbs), "")

r = client.post("/api/high-throughput", data={"file": (_io.BytesIO(b"x"), "in.zip")},
                content_type="multipart/form-data")
check(100, "a corrupt archive gives 400", r.status_code == 400, str(r.status_code))

# ---- F103/F104: history ----------------------------------------------------
r = client.get("/api/history")
check(103, "history returns a list", r.status_code == 200
      and isinstance(r.get_json(), list), str(r.status_code))
history = r.get_json()
check(103, "history capped at 50", len(history) <= 50 and webapp.HISTORY_LIMIT == 50,
      str(len(history)))
check(103, "entries carry id, filename, timestamp and report",
      all(k in history[0] for k in ("id", "filename", "timestamp", "report")),
      str(sorted(history[0].keys()))[:70])

r = client.get(f"/api/history/pdb/{BASE_SESSION}")
check(104, "history pdb served", r.status_code == 200 and b"ATOM" in r.data,
      str(r.status_code))
r = client.get(f"/api/history/pdb/{_uuid.uuid4()}")
check(104, "unknown job gives 404", r.status_code == 404, str(r.status_code))
r = client.get("/api/history/pdb/not-a-uuid")
check(104, "non-uuid job id rejected by the router", r.status_code == 404,
      str(r.status_code))

# ---- F106/F107: templates --------------------------------------------------
r = client.post("/api/templates", json={"name": "verify-tmpl",
                                        "settings": {"ph": 6.5}})
check(106, "template saved", r.status_code == 200
      and "verify-tmpl" in r.get_json()["templates"], str(r.status_code))
r = client.get("/api/templates")
check(106, "template listed", "verify-tmpl" in r.get_json(), "")
check(106, "template keeps its settings and a timestamp",
      r.get_json()["verify-tmpl"]["settings"] == {"ph": 6.5}
      and "created_at" in r.get_json()["verify-tmpl"], "")
r = client.post("/api/templates", json={"name": ""})
check(106, "template needs a name and settings", r.status_code == 400,
      str(r.status_code))
r = client.delete("/api/templates/verify-tmpl")
check(107, "template deleted", r.status_code == 200
      and "verify-tmpl" not in r.get_json()["templates"], str(r.status_code))
r = client.delete("/api/templates/never-existed")
check(107, "deleting a missing template is not an error", r.status_code == 200,
      str(r.status_code))

# ---- F108: report download -------------------------------------------------
r = client.post("/api/download-report", json=plain["report"])
check(108, "report downloads as text", r.status_code == 200
      and b"BIOPREP" in r.data, str(r.status_code))
check(108, "served as an attachment",
      "attachment" in r.headers.get("Content-Disposition", ""),
      r.headers.get("Content-Disposition", ""))
r = client.post("/api/download-report", json={"input_file": "partial.pdb"})
check(108, "a partial report still renders", r.status_code == 200, str(r.status_code))
r = client.post("/api/download-report", json=[])
check(108, "a non-object body gives 400", r.status_code == 400, str(r.status_code))

# ---- F109/F110: binding site analysis --------------------------------------
r = post("/api/analyze-site", {"file": upload(PLAIN)})
check(109, "site analysis on an upload", r.status_code == 200, str(r.status_code))
body = r.get_json()
check(109, "returns sites and a summary",
      "sites" in body and "summary" in body, str(sorted(body.keys())))
check(109, "sites carry volume and drugability",
      all("volume" in s and "drugability_score" in s for s in body["sites"]), "")

r = client.post(f"/api/analyze-site/{BASE_SESSION}")
check(110, "site analysis from a stored session", r.status_code == 200,
      str(r.status_code))
check(110, "session route finds pockets", len(r.get_json()["sites"]) > 0,
      str(len(r.get_json()["sites"])))
r = client.post(f"/api/analyze-site/{_uuid.uuid4()}")
check(110, "unknown session gives 400 with a code", r.status_code == 400
      and r.get_json().get("code") == "NO_STRUCTURE", str(r.status_code))

# ---- temp hygiene ----------------------------------------------------------
import glob                                                    # noqa: E402
leaked = [d for d in glob.glob(os.path.join(tempfile.gettempdir(), "bioprep_*"))
          if os.path.isdir(d)]
check(84, "no working directories leaked across all requests",
      not leaked, f"{len(leaked)} left")


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
print(f"{len(RESULTS)} checks across features 81-110, {failed} failed")
