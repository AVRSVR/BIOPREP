FROM python:3.13-slim

# OpenBabel gives PDBQT/AutoDock/Vina export. Not required for the app to run -
# convert_to_pdbqt() already fails cleanly with a clear message when the
# `obabel` binary isn't on PATH - but installing it here means that feature
# actually works on the deployed instance instead of always reporting missing.
RUN apt-get update && apt-get install -y --no-install-recommends \
    openbabel \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# The real package lives at bioprep/bioprep/ (protein prep/bioprep/bioprep/).
# Running with this directory as the working directory - not the repo root -
# is what makes `import bioprep` resolve to it instead of shadowing it with
# the outer directory as an empty namespace package; see conftest.py's own
# note on this for local test runs.
COPY bioprep/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt gunicorn

COPY bioprep/ ./

# Ephemeral: fine, since BIOPREP_PUBLIC_DEMO=1 turns off the only features
# (History, Templates) that would otherwise need this to survive a restart.
ENV BIOPREP_DATA_DIR=/app/data
ENV BIOPREP_PUBLIC_DEMO=1
RUN mkdir -p /app/data

EXPOSE 8000
# A single sync worker was tried first and found the hard way: one visitor
# hitting docking export (an OpenBabel subprocess call) wedged that worker,
# and since it's the *only* one, every other visitor's request - including
# the plain homepage - queued behind it with nothing else able to answer.
# gthread keeps the memory footprint close to one process (threads share the
# already-loaded numpy/scipy/OpenMM import, unlike separate worker
# processes) while making sure a slow or stuck request no longer blocks
# everyone else's. One long timeout because a real minimisation run
# legitimately takes minutes, not the framework's usual few seconds. Shell
# form so $PORT (set by the hosting platform) actually expands; falls back
# to 8000 for a local docker run.
CMD gunicorn --worker-class gthread --workers 1 --threads 4 --timeout 600 --bind 0.0.0.0:${PORT:-8000} bioprep.webapp:app
