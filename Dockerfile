FROM python:3.11-slim

LABEL org.opencontainers.image.title="RSeQC" \
    org.opencontainers.image.description="RSeQC 5.0.4 runtime (Python 3)" \
    org.opencontainers.image.version="5.0.4"

ENV PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1 \
    DEBIAN_FRONTEND=noninteractive \
    PYTHONPATH=/opt/rseqc/src

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        procps \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/rseqc

# Copy metadata and source in separate steps for better build cache reuse.
COPY pyproject.toml setup.py setup.cfg README.md MANIFEST.in LICENSE /opt/rseqc/
COPY src /opt/rseqc/src
COPY scripts /opt/rseqc/scripts

# Install from project metadata so dependency versions stay in one place.
RUN python -m pip install --no-cache-dir --upgrade pip setuptools wheel \
    && python -m pip install --no-cache-dir .

# Create non-root user for safer container execution.
RUN groupadd --system rseqc \
    && useradd --system --gid rseqc --create-home --shell /bin/bash rseqc \
    && chown -R rseqc:rseqc /opt/rseqc

USER rseqc

# Keep runtime command overridable for workflow engines like Nextflow.
CMD ["python3", "-u", "scripts/tin.py", "--help"]
