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

# ---- Pipeline scripts ----
COPY --chown=$MAMBA_USER:$MAMBA_USER 01_convert_seurat_to_h5ad.R 02_run_liana_ccc_pipeline.py entrypoint.sh /pipeline/

USER root
RUN chmod +x /pipeline/entrypoint.sh
USER $MAMBA_USER

# Mount your data here, e.g.: docker run -v $PWD/data:/data ...
WORKDIR /data

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/pipeline/entrypoint.sh"]
CMD ["--help"]
