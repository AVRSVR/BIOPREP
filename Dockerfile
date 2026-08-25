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
# Single worker: this app's own request-time compute (minimisation) is
# CPU/GPU-bound already: more workers would just contend for the same core
# rather than add real throughput, and cost more memory doing it. One long
# timeout because a real minimisation run legitimately takes minutes, not
# the framework's usual few seconds. Shell form so $PORT (set by the hosting
# platform) actually expands; falls back to 8000 for a local docker run.
CMD gunicorn --workers 1 --timeout 600 --bind 0.0.0.0:${PORT:-8000} bioprep.webapp:app
