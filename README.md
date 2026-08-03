# Predicting cell-cell communication from single cell data 

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
- `lr_results.csv` / `metabolite_sensor_results.csv` — the full, unfiltered
  results for every source/target cell-type pair, with `magnitude_rank` and
  `specificity_rank` (lower = stronger/more specific evidence, aggregated
  across CellPhoneDB, Connectome, log2FC, NATMI, SingleCellSignalR). In the
  metabolite file, `source` = metabolite (e.g. Prostaglandin J2) and
  `target`/`receptor` = the sensing cell type/gene; metabolite abundance per
  cell is estimated from production/degradation enzyme expression
  (MetalinksDB), then treated as the "ligand" in the same rank_aggregate
  scoring used for LR.
- `lr_results_filtered.csv` / `metabolite_sensor_results_filtered.csv` — the
  same results restricted to `specificity_rank <= --specificity_cutoff`
  (default 0.05) and sorted by `magnitude_rank` ascending (strongest first).
  This is a reasonable starting point for "which of these are worth trusting
  as bona fide," but see Caveats below — it's not sufficient on its own.
  Skip this file with `--no_filtered_output`, or pass `--specificity_cutoff
  1.0` to keep it but effectively disable the filter.
- `outgoing_ligand_sets.gmt` / `incoming_receptor_sets.gmt` /
  `combined_signal_sets.gmt` — per-cell-type gene sets (standard `.gmt`
  format: `name<TAB>description<TAB>gene1<TAB>gene2...`), one gene set per
  cell type per file, built from `lr_results_filtered.csv` only (not the
  metabolite step, since a metabolite name isn't a gene). Outgoing = ligand
  genes from interactions where that cell type is the source; incoming =
  receptor genes where it's the target; combined = the union of the two for
  that cell type. Complex subunits (e.g. `ITGAV_ITGB3`) are split into
  individual gene symbols. Skip with `--no_gmt_output`.

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
- `--min_cells` (default 100): minimum cells per cell-type group required to
  compute stats for that group at all. Set lower if you have genuinely small
  but biologically important populations you don't want dropped.
- `--use_raw`: only pass this if your `.h5ad` actually has `.raw` populated;
  `01_convert_seurat_to_h5ad.R` writes log-normalized data into `.X`, not
  `.raw`, so the default is `False`.
- `--specificity_cutoff` (default 0.05) / `--no_filtered_output`: control the
  `*_filtered.csv` files described above.
- `--no_gmt_output`: skip the `.gmt` gene-set files (also gated on the same
  `--specificity_cutoff`).

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
- The `specificity_rank <= 0.05` + magnitude-sort filter is a floor, not a
  finish line. Treat computational output as hypothesis-generating: also
  check consistency across donors/samples (not automated here), spatial or
  known-anatomy plausibility (dissociated scRNA-seq has no spatial context),
  and whether hits are dominated by ubiquitously-expressed genes rather than
  truly cell-type-specific ones.
- Container permissions: some Docker hosts (NFS-backed storage, SELinux,
  userns-remap) reject writes to bind-mounted `/data` even when the UID
  looks like it should match. To sidestep this, the container runs as root
  internally and `entrypoint.sh` chowns everything under `/data` back to
  match `/data`'s own owner afterward — don't pass `--user` to `docker run`
  for this image.

## Sources consulted

- [LIANA+ steady-state LR inference tutorial](https://liana-py.readthedocs.io/en/latest/notebooks/basic_usage.html)
- [LIANA+ multi-modal / metabolite-mediated CCC tutorial](https://liana-py.readthedocs.io/en/latest/notebooks/sc_multi.html)
- [LIANA+ installation guide](https://liana-py.readthedocs.io/en/latest/installation.html)
- [LIANA+ paper, Nature Cell Biology](https://www.nature.com/articles/s41556-024-01469-w)
- [saezlab/liana-py GitHub](https://github.com/saezlab/liana-py/)
