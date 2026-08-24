"""
Run every verification and science script and summarise the result.

    python verification/run_all.py            # everything
    python verification/run_all.py features   # catalogue features 1-121
    python verification/run_all.py science    # physics and chemistry only

These are deliberately separate from tests/test_backend.py. The unit tests
guard against regressions and run in a couple of minutes; these scripts check
each catalogue feature and each scientific claim end to end against real
deposited structures, and take longer.
"""

import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PYTHON = os.path.join(ROOT, ".venv", "Scripts", "python.exe")
if not os.path.exists(PYTHON):
    PYTHON = sys.executable

FEATURE_SCRIPTS = [
    "verify_1_11.py",       # io.py loading and saving
    "verify_12_22.py",      # protonator.py
    "verify_23_28.py",      # analyzer.py
    "verify_29_48.py",      # minimizer.py
    "verify_49_69.py",      # site_analyzer.py
    "verify_70_80.py",      # exporter.py and reporter.py
    "verify_81_110.py",     # app.py HTTP layer
    "verify_111_121.py",    # cli.py and frontend
]

SCIENCE_SCRIPTS = [
    "science_01.py",        # disulfides, pH, tier-2 pocket collapse
    "science_02.py",        # real complexes: collapse, rotation, water
    "science_03.py",        # pKa limits, tautomers, charges, convergence
]

SUMMARY = re.compile(r"(\d+) checks across (features [\d\-]+), (\d+) failed")


def run(script):
    path = os.path.join(HERE, script)
    started = time.time()
    # -u: unbuffered stdout. If a run is killed externally (a timeout, a
    # session teardown) before it exits normally, Python never flushes a
    # block-buffered stream, so everything but unbuffered logger output is
    # lost. Unbuffered mode means each line is on disk as it is printed.
    proc = subprocess.run([PYTHON, "-u", path], capture_output=True, text=True)
    elapsed = time.time() - started

    output = proc.stdout
    failures = [l for l in output.splitlines() if "FAIL" in l]
    match = SUMMARY.search(output)

    if match:
        total, label, failed = match.group(1), match.group(2), int(match.group(3))
        status = "ok" if failed == 0 else f"{failed} FAILED"
        print(f"  {script:<22} {label:<18} {total:>3} checks  {elapsed:6.1f}s  {status}")
    else:
        status = "ok" if proc.returncode == 0 else "ERROR"
        print(f"  {script:<22} {'science':<18} {'':>3}          {elapsed:6.1f}s  {status}")
        failed = 0 if proc.returncode == 0 else 1

    for line in failures:
        print(f"      {line.strip()}")
    if proc.returncode != 0 and not match:
        print(f"      {proc.stderr.strip().splitlines()[-1] if proc.stderr else ''}")

    return failed, len(failures)


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    scripts = []
    if which in ("all", "features"):
        scripts += FEATURE_SCRIPTS
    if which in ("all", "science"):
        scripts += SCIENCE_SCRIPTS
    if not scripts:
        print(__doc__)
        return 2

    print(f"Running {len(scripts)} scripts with {PYTHON}")
    print()
    total_failed = 0
    started = time.time()
    for script in scripts:
        failed, _ = run(script)
        total_failed += failed
    print()
    print(f"{len(scripts)} scripts in {time.time() - started:.0f}s, "
          f"{total_failed} failing checks")
    return 1 if total_failed else 0


if __name__ == "__main__":
    sys.exit(main())
