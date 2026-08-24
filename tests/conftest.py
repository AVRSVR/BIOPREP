"""
Put the real package ahead of the project directory on sys.path.

The layout is ``protein prep/bioprep/bioprep/``. Run from the repository root,
``import bioprep`` resolves to the *outer* directory as an implicit namespace
package, which shadows the real one and makes ``bioprep.core`` unimportable.
Inserting the inner path first resolves it.
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT = os.path.join(ROOT, 'bioprep')

for path in (PROJECT,):
    if path in sys.path:
        sys.path.remove(path)
    sys.path.insert(0, path)

# minimizer._make_simulation() prefers a GPU platform (see its docstring for
# why). GPU context creation and kernel compilation cost real fixed time per
# Simulation, though - a win for one real job, but a 5x slowdown across a
# suite that minimizes tiny fixtures thirty-odd times over (measured: 104s ->
# 514s on this machine). Tests care about the science being correct, not
# which hardware ran it, so default the suite to CPU; a test that wants to
# exercise the real cascade clears this itself for its own duration.
os.environ.setdefault('BIOPREP_MINIMIZER_PLATFORM', 'CPU')
