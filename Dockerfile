FROM kbase/sdkpython:3.8.10
LABEL maintainer="ac.shahnam"
USER root

# -----------------------------------------------------------------------------
# Layer 1 (slow, rarely changes): micromamba, with a fallback mirror.
# RUN curl rather than ADD: buildkit's ADD may not honour the build proxy.
# -----------------------------------------------------------------------------
RUN ( curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xvj -C /usr/local bin/micromamba ) \
    || ( curl -fsSL https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64.tar.bz2 \
        | tar -xvj -C /usr/local bin/micromamba ) \
 && /usr/local/bin/micromamba --version

# -----------------------------------------------------------------------------
# Layer 2 (slow): isolated RepeatMasker env. NOT on the image PATH - it would
# shadow the KBase server's Python. The Runner injects it per subprocess.
# -----------------------------------------------------------------------------
COPY env-repeatmasker.yml /tmp/env-repeatmasker.yml
RUN /usr/local/bin/micromamba create -y -r /opt/conda -n repeatmasker -f /tmp/env-repeatmasker.yml \
 && MAMBA_ROOT_PREFIX=/opt/conda /usr/local/bin/micromamba clean -a -y

# -----------------------------------------------------------------------------
# Layer 3: point FamDB at the KBase reference data directory (/data/famdb).
# famdb.py reads FAMDB_DATA_DIR from famdb.conf in its install dir, otherwise
# Libraries/famdb beside famdb.py - set both so either lookup lands on /data.
# -----------------------------------------------------------------------------
RUN set -e; \
    FAMDB_SHARE=$(ls -d /opt/conda/envs/repeatmasker/share/famdb-* | head -1); \
    test -x "$FAMDB_SHARE/famdb.py"; \
    rm -rf "$FAMDB_SHARE/Libraries/famdb"; \
    mkdir -p "$FAMDB_SHARE/Libraries"; \
    ln -s /data/famdb "$FAMDB_SHARE/Libraries/famdb"; \
    printf '[famdb]\nFAMDB_DATA_DIR = /data/famdb\n' > "$FAMDB_SHARE/famdb.conf"; \
    echo "FamDB wired to /data/famdb in $FAMDB_SHARE"

# -----------------------------------------------------------------------------
# Layer 4: build-time check with the RUNTIME PATH (scoped to this RUN only).
# Every helper RepeatMasker calls by bare name must resolve.
# -----------------------------------------------------------------------------
RUN export PATH=/opt/conda/envs/repeatmasker/bin:$PATH \
 && RepeatMasker -v \
 && for t in RepeatMasker ProcessRepeats famdb.py rmblastn makeblastdb trf nhmmer perl python3 curl md5sum gunzip; do \
        command -v "$t" >/dev/null || { echo "MISSING helper: $t"; exit 1; }; done \
 && python3 -c "import h5py; print('h5py', h5py.__version__)" \
 && famdb.py -h | head -3 \
 && echo REPEATMASKER_ENV_OK

# -----------------------------------------------------------------------------
# Layer 5: KBase server-side Python deps (new layer AFTER the conda env so
# editing it never forces a conda re-solve).
# -----------------------------------------------------------------------------
RUN pip install --no-cache-dir jsonrpcbase pytest pytest-cov coverage
RUN python -c "import sys; print('server python:', sys.executable)" \
 && python -c "import pytest, jsonrpcbase; print('server deps OK')"

# -----------------------------------------------------------------------------
# Volatile layers: module code.
# -----------------------------------------------------------------------------
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
