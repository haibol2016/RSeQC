# RSeQC
**This is an improved version of RSeQC, downloaded from [PyPI](https://pypi.org/project/RSeQC/). I added MPIRE-based parallel computing to tin.py, update the documentation and add Docker support.**

RSeQC provides quality control tools for RNA-seq data.

## Requirements

- Python 3.7+
- pip

Runtime dependencies are declared in `pyproject.toml`.

## Install From PyPI

```bash
python3 -m pip install rseqc
```

Upgrade:

```bash
python3 -m pip install --upgrade rseqc
```

Uninstall:

```bash
python3 -m pip uninstall -y rseqc
```

## Install From Source

From the repository root:

```bash
python3 -m pip install .
```

Editable install for development:

```bash
python3 -m pip install -e .
```

## Docker

Build image:

```bash
docker build -t rseqc:5.0.4 .
```

Run a quick smoke test:

```bash
docker run --rm --user "$(id -u):$(id -g)" rseqc:5.0.4 python3 -u scripts/tin.py --help
```

## Documentation

- Project site: http://rseqc.sourceforge.net/
- Additional files: `doc/`

## Notes

- The top-level `README.md` is the canonical project readme.
- `doc/README.md` is a short pointer to this file to avoid duplicated maintenance.
