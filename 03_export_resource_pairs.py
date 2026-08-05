#!/usr/bin/env python3
"""
03_export_resource_pairs.py

Exports the raw prior-knowledge ligand-receptor and metabolite-sensor
resources that 02_run_liana_ccc_pipeline.py draws on -- independent of any
particular dataset. Useful for auditing exactly which pairs were eligible
for a given run, or for browsing/reusing the underlying databases on their
own (e.g. as input to other tools).

No .h5ad / Seurat input is needed here: these are prior-knowledge databases,
queried before any expression data is ever considered.

Install (once):
  pip install liana --break-system-packages

Usage:
  python 03_export_resource_pairs.py --outdir resource_pairs/

Match this to a specific 02_run_liana_ccc_pipeline.py run by passing the
same --resource / --biospecimen values you used there.

Outputs (written to --outdir):
  lr_resource_pairs.csv                  -- full LR resource (ligand, receptor,
                                             plus any annotation columns the
                                             resource itself carries)
  metabolite_receptor_pairs.csv          -- MetalinksDB metabolite -> receptor
                                             pairs ("lr" type rows)
  metabolite_production_degradation.csv  -- MetalinksDB metabolite production/
                                             degradation enzyme pairs ("pd"
                                             type rows) -- these are what
                                             02_run_liana_ccc_pipeline.py uses
                                             to estimate metabolite abundance
                                             from enzyme expression, not
                                             receptor interactions themselves
"""

import argparse
import os

import liana as li


DEFAULT_LR_RESOURCE = "consensus"  # 'consensus' (human) or 'mouseconsensus' (mouse)
DEFAULT_BIOSPECIMEN = "Blood"      # MetalinksDB tissue/biofluid filter (e.g. Blood, Tissue)
DEFAULT_MET_SOURCES = [
    "CellPhoneDB", "Cellinker", "scConnect", "NeuronChat",  # metabolite-receptor sub-resources
    "recon", "hmr", "rhea", "hmdb", "Stich",                # production/degradation + broad
                                                             # chemical-protein sources. NOTE: 'Stich'
                                                             # (not 'STITCH') is the literal value used
                                                             # in MetalinksDB's own source column.
]


def export_lr_resource(resource_name, outdir):
    print(f"Fetching LR resource: '{resource_name}' (see li.resource.show_resources() for the full list)")
    resource = li.resource.select_resource(resource_name)
    out_path = os.path.join(outdir, "lr_resource_pairs.csv")
    resource.to_csv(out_path, index=False)
    print(f"Saved {resource.shape[0]} ligand-receptor pairs -> {out_path}")
    return resource


def export_metalinks_resource(biospecimen_location, met_sources, outdir):
    print(f"Fetching MetalinksDB (biospecimen_location={biospecimen_location!r}, sources={met_sources})")
    metalinks = li.resource.get_metalinks(
        biospecimen_location=biospecimen_location,
        source=met_sources,
        types=["pd", "lr"],
    )

    lr = metalinks[metalinks["type"] == "lr"].copy()
    lr_path = os.path.join(outdir, "metabolite_receptor_pairs.csv")
    lr.to_csv(lr_path, index=False)
    print(
        f"Saved {lr.shape[0]} metabolite-receptor pairs "
        f"({lr['metabolite'].nunique()} metabolites, {lr['gene_symbol'].nunique()} receptor genes) -> {lr_path}"
    )

    pd_ = metalinks[metalinks["type"] == "pd"].copy()
    pd_path = os.path.join(outdir, "metabolite_production_degradation.csv")
    pd_.to_csv(pd_path, index=False)
    print(
        f"Saved {pd_.shape[0]} production/degradation enzyme pairs "
        f"({pd_['metabolite'].nunique()} metabolites, {pd_['gene_symbol'].nunique()} enzyme genes) -> {pd_path}"
    )

    return lr, pd_


def main():
    parser = argparse.ArgumentParser(
        description="Export the raw LR + metabolite-sensor prior-knowledge resources (no dataset required)"
    )
    parser.add_argument("--resource", default=DEFAULT_LR_RESOURCE,
                         help="LIANA+ LR resource name, e.g. 'consensus' (human) or 'mouseconsensus' (mouse)")
    parser.add_argument("--biospecimen", default=DEFAULT_BIOSPECIMEN,
                         help="MetalinksDB tissue/biofluid filter, e.g. Blood, Tissue")
    parser.add_argument("--outdir", default="resource_pairs")
    parser.add_argument("--skip_lr", action="store_true", help="Skip exporting the LR resource")
    parser.add_argument("--skip_metabolite", action="store_true", help="Skip exporting the metabolite-sensor resource")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if not args.skip_lr:
        export_lr_resource(args.resource, args.outdir)

    if not args.skip_metabolite:
        export_metalinks_resource(args.biospecimen, DEFAULT_MET_SOURCES, args.outdir)

    print("\nDone. Resource pairs written to:", os.path.abspath(args.outdir))


if __name__ == "__main__":
    main()
