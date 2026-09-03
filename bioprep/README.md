# BioPrep (package)

This is the installable package. See the [project README](../README.md) for what
BioPrep does, how to run it, and the engineering notes.

## Install

```bash
pip install -e .
```

Exposes two entry points:

- `bioprep` — the CLI (`bioprep --input in.pdb --output out.pdb --chain A --minimize`)
- `bioprep-web` — the Flask web app

Mutable state (job history, saved templates, processed structures) is written to
`BIOPREP_DATA_DIR`, defaulting to the working directory — deliberately not
inside the package, since `site-packages` is frequently read-only and
reinstalling would otherwise delete it.
