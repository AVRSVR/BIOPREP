"""
Run the BioPrep web app from a source checkout.

The application itself lives in ``bioprep.webapp`` so that a ``pip install``ed
copy can serve it - templates and static files ship as package data, and
``bioprep-web`` is installed as a console script. This shim exists so the
long-standing ``python app.py`` workflow keeps working, and so that anything
importing ``app`` still finds the Flask object.

It also pins the data directory to this folder, which is where the existing
templates_store.json, jobs_history.json and processed_data/ already live.
Without that, running from source would start writing them to whatever the
current working directory happened to be.
"""

import os

os.environ.setdefault(
    'BIOPREP_DATA_DIR', os.path.dirname(os.path.abspath(__file__)))

from bioprep import webapp  # noqa: E402
from bioprep.webapp import app, main  # noqa: E402,F401

# Re-export the module's configuration as well, not just the Flask object, so
# anything that read app.RESULTS_DIR, app.SESSION_LIMIT or app.sessions still
# finds them here.
globals().update(
    {name: value for name, value in vars(webapp).items()
     if not name.startswith('_')})

if __name__ == '__main__':
    main()
