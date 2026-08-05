#!/usr/bin/env python3
"""
02_run_liana_ccc_pipeline.py

Predicts, from a single-cell AnnData (.h5ad) object:
  (A) Ligand-receptor (LR) signaling            -- LIANA+ rank_aggregate consensus
  (B) Metabolite-sensor (metabolite-receptor)   -- LIANA+ + MetalinksDB
      signaling

Input:
  A .h5ad file produced by 01_convert_seurat_to_h5ad.R (or any AnnData with
  log1p-normalized expression in .X and a categorical cell-type/cluster
  column in .obs).

Install (once):
  pip install liana scanpy mudata plotnine --break-system-packages

Usage:
  python 02_run_liana_ccc_pipeline.py \
      --h5ad my_data.h5ad \
      --groupby cell_type \
      --outdir liana_results/

Human data uses LIANA+'s default 'consensus' LR resource and MetalinksDB
(human-curated). For mouse, pass --resource mouseconsensus for the LR step;
for the metabolite step you'll need to map MetalinksDB gene symbols to mouse
orthologs first (see li.resource.get_hcop_orthologs / translate_resource in
the LIANA+ "Prior Knowledge" tutorial) -- not included here since this run
is configured for human data.

Filtering to a "bona fide" subset:
  In addition to the full, unfiltered results (lr_results.csv /
  metabolite_sensor_results.csv), this script writes a second, filtered +
  sorted file (lr_results_filtered.csv / metabolite_sensor_results_filtered.csv)
  keeping only interactions with specificity_rank <= --specificity_cutoff
  (default 0.05), ordered by magnitude_rank ascending (strongest first).
  This mirrors the common convention of filtering on specificity, then
  ranking by magnitude, rather than trusting either statistic alone.
  Use --no_filtered_output to skip writing the filtered file, or set
  --specificity_cutoff 1.0 to effectively disable the filter while still
  getting a magnitude-sorted file.

Gene set (.gmt) outputs:
  Built from the same specificity-filtered LR subset (not the metabolite
  step -- metabolite "ligands" are metabolite names, not genes), three .gmt
  files are written, each with one gene set per cell type:
    - outgoing_ligand_sets.gmt    ligand genes from interactions where the
                                   cell type is the source (sender)
    - incoming_receptor_sets.gmt  receptor genes from interactions where the
                                   cell type is the target (receiver)
    - combined_signal_sets.gmt    union of that cell type's outgoing ligand
                                   genes and incoming receptor genes
  Complex subunits (e.g. ITGAV_ITGB3) are split into individual gene
  symbols. Skip these with --no_gmt_output.
"""

import argparse
import os

import scanpy as sc
import liana as li


# ----------------------------------------------------------------------------
# CONFIG - sensible defaults; override via CLI flags below
# ----------------------------------------------------------------------------
DEFAULT_GROUPBY    = "cell_type"   # obs column with cell type / cluster labels
DEFAULT_LR_RESOURCE = "consensus"  # 'consensus' (human) or 'mouseconsensus' (mouse)
DEFAULT_EXPR_PROP  = 0.10          # min fraction of cells expressing L/R per group
DEFAULT_MIN_CELLS  = 100           # min cells per group; 100 is conservative for large atlases
                                   # (LIANA's own default is 5, which is low for datasets this size)
DEFAULT_N_PERMS    = 1000
DEFAULT_SPECIFICITY_CUTOFF = 0.05  # keep interactions with specificity_rank <= this value
DEFAULT_WRITE_FILTERED = True      # also write a specificity-filtered, magnitude-sorted CSV
DEFAULT_WRITE_GMT = True           # also write per-cell-type outgoing/incoming/combined .gmt files
DEFAULT_USE_RAW    = False         # sceasy (01_convert_seurat_to_h5ad.R) writes log-normalized
                                   # data straight into .X, not .raw -- set True only if your
                                   # .h5ad has .raw populated with the expression you want to use
DEFAULT_BIOSPECIMEN = "Blood"      # MetalinksDB tissue/biofluid filter (e.g. Blood, Tissue)
DEFAULT_MET_SOURCES = [
    "CellPhoneDB", "Cellinker", "scConnect", "NeuronChat",  # metabolite-receptor sub-resources
                                                             # (NeuronChat adds neurotransmitter/
                                                             # neuropeptide-receptor pairs)
    "recon", "hmr", "rhea", "hmdb", "Stich",                # production/degradation + broad
                                                             # chemical-protein sources. NOTE: 'Stich'
                                                             # (not 'STITCH') is the literal value used
                                                             # in MetalinksDB's own source column -- this
                                                             # isn't a typo on our part, it's theirs.
                                                             # STITCH is text-mined/computationally
                                                             # predicted, not literature-curated, so it
                                                             # adds much broader but lower-confidence
                                                             # coverage than the other sources above.
]


def filter_and_sort(res_df, specificity_cutoff, specificity_col="specificity_rank", magnitude_col="magnitude_rank"):
    """Keep interactions with specificity_col <= specificity_cutoff, sorted by
    magnitude_col ascending (strongest first). Returns None if the expected
    columns aren't present (e.g. a custom rank_aggregate method subset was
    used that doesn't produce them) so callers can skip writing the file."""
    if specificity_col not in res_df.columns or magnitude_col not in res_df.columns:
        print(
            f"  [!] '{specificity_col}' and/or '{magnitude_col}' not found in results "
            f"(columns present: {list(res_df.columns)}) -- skipping filtered output."
        )
        return None
    filtered = res_df[res_df[specificity_col] <= specificity_cutoff].copy()
    filtered = filtered.sort_values(magnitude_col, ascending=True)
    return filtered


def _split_complex(gene_str):
    """LIANA represents heteromeric complexes as underscore-joined subunit
    symbols (e.g. 'ITGAV_ITGB3'); split those into individual gene symbols
    for gene-set purposes. Plain single genes just pass through unchanged."""
    if gene_str is None:
        return []
    return [g for g in str(gene_str).split("_") if g]


def _sanitize_gmt_name(name):
    # GMT is tab-delimited; strip anything that would break that.
    return str(name).replace("\t", " ").replace("\n", " ").strip()


def _write_gmt(path, sets_by_celltype, suffix, description_label):
    """Write one gene set per cell type as a line in a .gmt file:
    <name>\t<description>\t<gene1>\t<gene2>\t...
    Cell types with an empty gene set are skipped (nothing to write)."""
    n_written = 0
    with open(path, "w") as fh:
        for cell_type in sorted(sets_by_celltype):
            genes = sorted(sets_by_celltype[cell_type])
            if not genes:
                continue
            name = f"{_sanitize_gmt_name(cell_type)}_{suffix}"
            description = f"{description_label} (n_genes={len(genes)})"
            fh.write("\t".join([name, description] + genes) + "\n")
            n_written += 1
    return n_written


def write_gmt_outputs(lr_res_filtered, specificity_cutoff, outdir):
    """From the specificity-filtered LR results, build three .gmt files:
    outgoing ligand genes per cell type, incoming receptor genes per cell
    type, and their per-cell-type union (combined)."""
    if lr_res_filtered is None or lr_res_filtered.empty:
        print("  [!] No filtered LR interactions available -- skipping .gmt outputs.")
        return

    outgoing = {}  # cell type (as source) -> set of ligand genes
    incoming = {}  # cell type (as target) -> set of receptor genes

    for _, row in lr_res_filtered.iterrows():
        outgoing.setdefault(row["source"], set()).update(_split_complex(row.get("ligand_complex")))
        incoming.setdefault(row["target"], set()).update(_split_complex(row.get("receptor_complex")))

    combined = {}
    for cell_type in set(outgoing) | set(incoming):
        combined[cell_type] = outgoing.get(cell_type, set()) | incoming.get(cell_type, set())

    desc = f"specificity_rank<={specificity_cutoff}"

    out_path = os.path.join(outdir, "outgoing_ligand_sets.gmt")
    n = _write_gmt(out_path, outgoing, "outgoing_ligands", f"Significant outgoing ligand genes; {desc}")
    print(f"Saved outgoing ligand gene sets ({n} cell types) -> {out_path}")

    in_path = os.path.join(outdir, "incoming_receptor_sets.gmt")
    n = _write_gmt(in_path, incoming, "incoming_receptors", f"Significant incoming receptor genes; {desc}")
    print(f"Saved incoming receptor gene sets ({n} cell types) -> {in_path}")

    comb_path = os.path.join(outdir, "combined_signal_sets.gmt")
    n = _write_gmt(comb_path, combined, "combined_signals", f"Significant outgoing + incoming genes; {desc}")
    print(f"Saved combined signal gene sets ({n} cell types) -> {comb_path}")


def run_lr_analysis(adata, groupby, resource_name, expr_prop, min_cells, n_perms, use_raw,
                     specificity_cutoff, write_filtered, write_gmt, outdir):
    print("\n=== [Step A] Ligand-receptor signaling (LIANA+ rank_aggregate) ===")
    li.mt.rank_aggregate(
        adata,
        groupby=groupby,
        resource_name=resource_name,
        expr_prop=expr_prop,
        min_cells=min_cells,
        n_perms=n_perms,
        use_raw=use_raw,
        key_added="liana_lr_res",
        verbose=True,
    )
    lr_res = adata.uns["liana_lr_res"]
    out_path = os.path.join(outdir, "lr_results.csv")
    lr_res.to_csv(out_path, index=False)
    print(f"Saved LR results ({lr_res.shape[0]} interactions) -> {out_path}")

    # Computed whenever either the filtered CSV or the .gmt gene sets are
    # wanted, since both are derived from the same significant subset.
    filtered = None
    if write_filtered or write_gmt:
        filtered = filter_and_sort(lr_res, specificity_cutoff)

    if write_filtered and filtered is not None:
        filt_path = os.path.join(outdir, "lr_results_filtered.csv")
        filtered.to_csv(filt_path, index=False)
        print(
            f"Saved filtered LR results (specificity_rank <= {specificity_cutoff}, "
            f"sorted by magnitude_rank): {filtered.shape[0]}/{lr_res.shape[0]} interactions -> {filt_path}"
        )

    if write_gmt and filtered is not None:
        write_gmt_outputs(filtered, specificity_cutoff, outdir)

    return lr_res


def run_metabolite_sensor_analysis(adata, groupby, biospecimen_location, met_sources,
                                    specificity_cutoff, write_filtered, outdir):
    print("\n=== [Step B] Metabolite-sensor signaling (MetalinksDB) ===")

    print("  Fetching MetalinksDB prior knowledge...")
    metalinks = li.resource.get_metalinks(
        biospecimen_location=biospecimen_location,
        source=met_sources,
        types=["pd", "lr"],
    )

    # Metabolite -> receptor ("sensor") resource
    resource = metalinks[metalinks["type"] == "lr"].copy()
    resource = (
        resource[["metabolite", "gene_symbol"]]
        .rename(columns={"metabolite": "source", "gene_symbol": "receptor"})
        .drop_duplicates()
    )

    # Production-degradation enzyme network, used to estimate per-cell
    # metabolite abundance from enzyme gene expression.
    pd_net = metalinks[metalinks["type"] == "pd"]
    pd_net = (
        pd_net[["metabolite", "gene_symbol", "mor"]]
        .groupby(["metabolite", "gene_symbol"])
        .agg("mean")
        .reset_index()
        .rename(columns={"metabolite": "source", "gene_symbol": "target", "mor": "weight"})
    )

    # Transporter network (optional refinement: export=+1, import=-1)
    t_net = metalinks[metalinks["type"] == "pd"]
    t_net = t_net[["metabolite", "gene_symbol", "transport_direction"]].dropna()
    t_net["mor"] = t_net["transport_direction"].apply(
        lambda x: 1 if x == "out" else -1 if x == "in" else None
    )
    t_net = (
        t_net[["metabolite", "gene_symbol", "mor"]]
        .dropna()
        .groupby(["metabolite", "gene_symbol"])
        .agg("mean")
        .reset_index()
    )
    t_net = t_net[t_net["mor"] != 0]
    t_net = t_net.rename(columns={"metabolite": "source", "gene_symbol": "target", "mor": "weight"})

    print("  Estimating per-cell metabolite abundances from enzyme expression...")
    meta = li.mt.fun.estimate_metalinks(
        adata,
        resource,
        pd_net=pd_net,
        t_net=t_net,
        use_raw=False,
        tmin=3,
    )
    meta.obs[groupby] = adata.obs[groupby]

    print("  Inferring metabolite-receptor (sensor) interactions with rank_aggregate...")
    li.mt.rank_aggregate(
        adata=meta,
        groupby=groupby,
        resource=resource.rename(columns={"source": "ligand"}),
        mdata_kwargs={
            "x_mod": "metabolite",
            "y_mod": "receptor",
            "x_use_raw": False,
            "y_use_raw": False,
            "x_transform": li.ut.zi_minmax,
            "y_transform": li.ut.zi_minmax,
        },
        key_added="liana_met_res",
        verbose=True,
    )

    met_res = meta.uns["liana_met_res"]
    out_path = os.path.join(outdir, "metabolite_sensor_results.csv")
    met_res.to_csv(out_path, index=False)
    print(f"Saved metabolite-sensor results ({met_res.shape[0]} interactions) -> {out_path}")

    if write_filtered:
        filtered = filter_and_sort(met_res, specificity_cutoff)
        if filtered is not None:
            filt_path = os.path.join(outdir, "metabolite_sensor_results_filtered.csv")
            filtered.to_csv(filt_path, index=False)
            print(
                f"Saved filtered metabolite-sensor results (specificity_rank <= {specificity_cutoff}, "
                f"sorted by magnitude_rank): {filtered.shape[0]}/{met_res.shape[0]} interactions -> {filt_path}"
            )

    return met_res, meta


def main():
    parser = argparse.ArgumentParser(description="LIANA+ ligand-receptor + metabolite-sensor CCC pipeline")
    parser.add_argument("--h5ad", required=True, help="Path to input .h5ad (from Seurat conversion)")
    parser.add_argument("--groupby", default=DEFAULT_GROUPBY, help="obs column with cell type/cluster labels")
    parser.add_argument("--resource", default=DEFAULT_LR_RESOURCE, help="'consensus' (human) or 'mouseconsensus' (mouse)")
    parser.add_argument("--expr_prop", type=float, default=DEFAULT_EXPR_PROP)
    parser.add_argument("--min_cells", type=int, default=DEFAULT_MIN_CELLS)
    parser.add_argument("--n_perms", type=int, default=DEFAULT_N_PERMS)
    parser.add_argument("--use_raw", action="store_true", default=DEFAULT_USE_RAW,
                         help="Use adata.raw instead of adata.X (only if adata.raw is populated)")
    parser.add_argument("--specificity_cutoff", type=float, default=DEFAULT_SPECIFICITY_CUTOFF,
                         help="Keep interactions with specificity_rank <= this value in the "
                              "*_filtered.csv output (default 0.05). Set to 1.0 to keep everything.")
    parser.add_argument("--no_filtered_output", action="store_true",
                         help="Skip writing the specificity-filtered, magnitude-sorted *_filtered.csv files")
    parser.add_argument("--no_gmt_output", action="store_true",
                         help="Skip writing outgoing_ligand_sets.gmt / incoming_receptor_sets.gmt / "
                              "combined_signal_sets.gmt (built from the LR step only)")
    parser.add_argument("--biospecimen", default=DEFAULT_BIOSPECIMEN, help="MetalinksDB tissue/biofluid filter")
    parser.add_argument("--outdir", default="liana_results")
    parser.add_argument("--skip_lr", action="store_true", help="Skip the ligand-receptor step")
    parser.add_argument("--skip_metabolite", action="store_true", help="Skip the metabolite-sensor step")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading AnnData: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)

    if args.groupby not in adata.obs.columns:
        raise ValueError(
            f"'{args.groupby}' not found in adata.obs. Available columns: {list(adata.obs.columns)}"
        )
    adata.obs[args.groupby] = adata.obs[args.groupby].astype("category")

    write_filtered = not args.no_filtered_output
    write_gmt = not args.no_gmt_output

    if not args.skip_lr:
        run_lr_analysis(
            adata, args.groupby, args.resource, args.expr_prop,
            args.min_cells, args.n_perms, args.use_raw,
            args.specificity_cutoff, write_filtered, write_gmt, args.outdir,
        )

    if not args.skip_metabolite:
        run_metabolite_sensor_analysis(
            adata, args.groupby, args.biospecimen, DEFAULT_MET_SOURCES,
            args.specificity_cutoff, write_filtered, args.outdir,
        )

    print("\nDone. Results written to:", os.path.abspath(args.outdir))


if __name__ == "__main__":
    main()
