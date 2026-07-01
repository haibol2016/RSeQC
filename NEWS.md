# NEWS

## 2026-07-01 

By haibol2017@gmail.com

### Packaging

- Simplified setup.py into a minimal compatibility shim that calls setuptools.setup().
- Moved script registration to pyproject.toml via [tool.setuptools].script-files.
- Added MPIRE-powered chromosome-level parallel execution in scripts/tin.py
  using WorkerPool, controlled by `-p/--processes`.
- Added and pinned runtime dependencies in pyproject.toml:
  - pysam>=0.15.0
  - bx-python>=0.8.13
  - numpy>=1.17.0
  - pyBigWig>=0.3.18
  - logomaker>=0.8
  - pandas>=0.25.0
  - matplotlib>=3.0.0
  - mpire>=2.0.0
- Updated build-system requirements to setuptools>=61 and wheel.
- Updated requires-python to >=3.7 to align with dependency constraints.
- Added explicit Python classifiers for 3.7 through 3.12 and Python 3 only.

### Docker

- Updated Dockerfile to copy pyproject.toml and LICENSE into the build context.
- Removed manual pip dependency list from Dockerfile.
- Switched Docker install step to `pip install .` so dependency resolution is driven by pyproject.toml.
- Built image successfully as rseqc:5.0.4.
- Verified runtime smoke test:
  - python3 -u scripts/tin.py --help

### Documentation

- Reworked top-level README.md into a canonical project guide with:
  - requirements
  - PyPI install/upgrade/uninstall
  - source and editable install
  - Docker build and smoke test commands
- Updated doc/README.md to a short pointer to top-level README.md to avoid duplicate maintenance.
- Updated README Docker run example to include --user "$(id -u):$(id -g)".

### Repository Hygiene

- Added .gitignore with macOS, Python cache, build artifact, virtual environment, test cache, and IDE ignore patterns.
