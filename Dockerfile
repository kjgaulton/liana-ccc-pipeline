# LIANA+ ligand-receptor + metabolite-sensor CCC pipeline
#
# Combines R (Seurat + sceasy, for .rds -> .h5ad conversion) and Python
# (liana-py, for the actual signaling inference) in one conda/micromamba
# environment, since sceasy's reticulate bridge needs both languages to
# share the same Python interpreter.
#
# Build:  docker build -t liana-ccc-pipeline .
# Run:    docker run --rm liana-ccc-pipeline --help

FROM mambaorg/micromamba:1.5.9

# ---- Combined R + Python environment ----
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# sceasy (Seurat -> AnnData converter) has no conda/CRAN package;
# install from GitHub. Point reticulate at this same env's Python,
# which already has anndata installed via environment.yml.
ENV RETICULATE_PYTHON=/opt/conda/bin/python
RUN R -e "remotes::install_github('cellgeni/sceasy', upgrade = 'never')"

# decoupler (a dependency pulled in for LIANA+'s metabolite-estimation step,
# li.mt.fun.estimate_metalinks) JIT-compiles some functions with numba's
# on-disk caching (@nb.njit(cache=True)). Numba's cache locator can't find a
# writable location for the installed package path in some container setups,
# raising "RuntimeError: cannot cache function ...: no locator available".
# Pointing NUMBA_CACHE_DIR at a plain writable directory avoids that lookup.
ENV NUMBA_CACHE_DIR=/tmp/numba_cache
RUN mkdir -p /tmp/numba_cache && chmod -R 777 /tmp/numba_cache

# ---- Pipeline scripts ----
# Runs as root (see below) so it can always read/execute these regardless of
# host-specific permission quirks on some Docker setups (NFS-backed storage,
# SELinux, userns-remap, etc.) that caused "Permission denied" even when the
# numeric UID appeared to match the file owner.
USER root
COPY 01_convert_seurat_to_h5ad.R 02_run_liana_ccc_pipeline.py entrypoint.sh /pipeline/
RUN chmod -R 755 /pipeline

# Mount your data here, e.g.: docker run -v $PWD/data:/data ...
# Do NOT pass --user to `docker run` for this image; entrypoint.sh chowns
# anything it writes under /data back to match /data's own owner, so you
# still end up with normal (non-root) file ownership on the host side.
WORKDIR /data

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "bash", "/pipeline/entrypoint.sh"]
CMD ["--help"]
