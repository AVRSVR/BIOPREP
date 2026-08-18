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
