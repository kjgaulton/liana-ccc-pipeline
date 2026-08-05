# Cell-cell communication from single cell data 

Predicts ligand-receptor and metabolite-sensor signaling from a single-cell
Seurat object using LIANA+. LIANA+ is Python-only (scverse ecosystem,
AnnData/MuData), so a Seurat `.rds` is converted to `.h5ad` first.

```
Seurat .rds ──(R: sceasy)──> .h5ad ──(Python: liana-py)──> LR results + metabolite-sensor results
```

## Files

- `Dockerfile`, `environment.yml`, `entrypoint.sh` — containerized version of
  the whole pipeline (recommended way to run this reproducibly).
- `01_convert_seurat_to_h5ad.R` — converts the Seurat object to AnnData.
- `02_run_liana_ccc_pipeline.py` — runs LIANA+ rank_aggregate for
  ligand-receptor signaling, and LIANA+'s MetalinksDB-based module for
  metabolite-sensor signaling.

## Run with Docker (recommended)

Build once:
```bash
docker build -t liana-ccc-pipeline .
```
This installs R + Seurat + sceasy and Python + liana-py into one image, so
build time is significant (Seurat compiles several packages) — expect
15–30 minutes on first build, cached afterward.

Put your `.rds` file in a local folder, e.g. `./data/my_data.rds`, then:

```bash
# Convert only
docker run --rm -v "$PWD/data:/data" liana-ccc-pipeline \
    convert /data/my_data.rds /data/my_data.h5ad RNA cell_type

# Analyze only (if you already have a .h5ad)
docker run --rm -v "$PWD/data:/data" liana-ccc-pipeline \
    analyze --h5ad /data/my_data.h5ad --groupby cell_type --outdir /data/liana_results

# Both steps in one shot
docker run --rm -v "$PWD/data:/data" liana-ccc-pipeline \
    all /data/my_data.rds /data/my_data.h5ad cell_type --outdir /data/liana_results
```

`--outdir` and other flags accepted by `02_run_liana_ccc_pipeline.py` (e.g.
`--resource`, `--expr_prop`, `--biospecimen`) can be appended after the
required arguments in `analyze` / `all`. Run with `--help` for the full
subcommand list.

Everything you mount to `/data` is visible inside the container and results
are written back out to that same host folder, since it's a bind mount.

## Run without Docker

**R side:**
```r
install.packages(c("Seurat", "remotes"))
remotes::install_github("cellgeni/sceasy")
reticulate::py_install("anndata")   # sceasy needs a Python env with anndata
```

**Python side:**
```bash
pip install liana scanpy mudata plotnine --break-system-packages
```

**Run:**
```bash
Rscript 01_convert_seurat_to_h5ad.R my_data.rds my_data.h5ad RNA cell_type

python 02_run_liana_ccc_pipeline.py \
    --h5ad my_data.h5ad \
    --groupby cell_type \
    --outdir liana_results/
```

## Outputs

Written to `liana_results/` (or wherever `--outdir` points):
- `lr_results.csv` — ligand-receptor interactions per source/target cell type
  pair, with `magnitude_rank` and `specificity_rank` (lower = stronger/more
  specific evidence, aggregated across CellPhoneDB, Connectome, log2FC, NATMI,
  SingleCellSignalR).
- `metabolite_sensor_results.csv` — same columns, but `source` = metabolite
  (e.g. Prostaglandin J2) and `target`/`receptor` = the sensing cell type/gene.
  Metabolite abundance per cell is estimated from production/degradation
  enzyme expression (MetalinksDB), then treated as the "ligand" in the same
  rank_aggregate scoring used for LR.

## Key parameters to check before running

- `--groupby`: must match an actual column in your Seurat metadata (cell
  type/cluster annotation).
- `--resource`: `consensus` for human, `mouseconsensus` for mouse LR pairs.
  MetalinksDB (the metabolite-sensor resource) is human-curated — running on
  mouse data requires ortholog mapping first (see LIANA+'s "Prior Knowledge"
  tutorial, `li.resource.get_hcop_orthologs`), not included here since this
  build targets human data.
- `--biospecimen`: MetalinksDB filter, e.g. `Blood`, `Tissue` — pick what's
  closest to your sample type; run `li.resource.get_metalinks_values()` in
  Python to see all valid options.
- `--expr_prop` (default 0.1): minimum fraction of cells per group expressing
  a ligand/receptor for an interaction to count as detected.

## Caveats (worth keeping in mind / reporting)

- Metabolite abundance is *inferred* from enzyme gene expression via linear
  regression — it is not measured. Treat metabolite-sensor hits as
  hypotheses to validate, not ground truth (this is explicitly flagged by
  the LIANA+/MetalinksDB authors too).
- `rank_aggregate` ranks are relative within your dataset; they are not
  p-values in the traditional sense, though `cellphone_pvals` (one of the
  underlying method columns) is permutation-based.
- Results depend heavily on cell-type granularity in `--groupby` — coarse
  annotations can mask communication between subtypes.
- The Dockerfile hasn't been build-tested in this environment (no Docker
  daemon available here) — the R/Python package names, Dockerfile syntax,
  YAML, and entrypoint script were all validated individually, but please
  run the actual `docker build` yourself and flag anything that breaks.

## Sources consulted

- [LIANA+ steady-state LR inference tutorial](https://liana-py.readthedocs.io/en/latest/notebooks/basic_usage.html)
- [LIANA+ multi-modal / metabolite-mediated CCC tutorial](https://liana-py.readthedocs.io/en/latest/notebooks/sc_multi.html)
- [LIANA+ installation guide](https://liana-py.readthedocs.io/en/latest/installation.html)
- [LIANA+ paper, Nature Cell Biology](https://www.nature.com/articles/s41556-024-01469-w)
- [saezlab/liana-py GitHub](https://github.com/saezlab/liana-py/)
