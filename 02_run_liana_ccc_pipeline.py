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
DEFAULT_MIN_CELLS  = 5
DEFAULT_N_PERMS    = 1000
DEFAULT_BIOSPECIMEN = "Blood"      # MetalinksDB tissue/biofluid filter (e.g. Blood, Tissue)
DEFAULT_MET_SOURCES = [
    "CellPhoneDB", "Cellinker", "scConnect",   # metabolite-receptor sub-resources
    "recon", "hmr", "rhea", "hmdb",            # production/degradation enzyme sources
]


def run_lr_analysis(adata, groupby, resource_name, expr_prop, min_cells, n_perms, outdir):
    print("\n=== [Step A] Ligand-receptor signaling (LIANA+ rank_aggregate) ===")
    li.mt.rank_aggregate(
        adata,
        groupby=groupby,
        resource_name=resource_name,
        expr_prop=expr_prop,
        min_cells=min_cells,
        n_perms=n_perms,
        key_added="liana_lr_res",
        verbose=True,
    )
    lr_res = adata.uns["liana_lr_res"]
    out_path = os.path.join(outdir, "lr_results.csv")
    lr_res.to_csv(out_path, index=False)
    print(f"Saved LR results ({lr_res.shape[0]} interactions) -> {out_path}")
    return lr_res


def run_metabolite_sensor_analysis(adata, groupby, biospecimen_location, met_sources, outdir):
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
    return met_res, meta


def main():
    parser = argparse.ArgumentParser(description="LIANA+ ligand-receptor + metabolite-sensor CCC pipeline")
    parser.add_argument("--h5ad", required=True, help="Path to input .h5ad (from Seurat conversion)")
    parser.add_argument("--groupby", default=DEFAULT_GROUPBY, help="obs column with cell type/cluster labels")
    parser.add_argument("--resource", default=DEFAULT_LR_RESOURCE, help="'consensus' (human) or 'mouseconsensus' (mouse)")
    parser.add_argument("--expr_prop", type=float, default=DEFAULT_EXPR_PROP)
    parser.add_argument("--min_cells", type=int, default=DEFAULT_MIN_CELLS)
    parser.add_argument("--n_perms", type=int, default=DEFAULT_N_PERMS)
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

    if not args.skip_lr:
        run_lr_analysis(
            adata, args.groupby, args.resource, args.expr_prop,
            args.min_cells, args.n_perms, args.outdir,
        )

    if not args.skip_metabolite:
        run_metabolite_sensor_analysis(
            adata, args.groupby, args.biospecimen, DEFAULT_MET_SOURCES, args.outdir,
        )

    print("\nDone. Results written to:", os.path.abspath(args.outdir))


if __name__ == "__main__":
    main()
