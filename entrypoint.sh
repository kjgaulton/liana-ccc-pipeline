#!/usr/bin/env bash
# Dispatcher for the LIANA+ CCC pipeline container.
set -euo pipefail

usage() {
  cat <<'EOF'
LIANA+ ligand-receptor + metabolite-sensor CCC pipeline

Usage:
  docker run -v <host_data_dir>:/data <image> convert <in.rds> <out.h5ad> [assay] [celltype_col]
  docker run -v <host_data_dir>:/data <image> analyze --h5ad <file.h5ad> --groupby <col> [more python args...]
  docker run -v <host_data_dir>:/data <image> all <in.rds> <out.h5ad> <celltype_col> [more python args...]

Subcommands:
  convert   Seurat .rds -> .h5ad          (runs 01_convert_seurat_to_h5ad.R)
  analyze   LR + metabolite-sensor scoring (runs 02_run_liana_ccc_pipeline.py)
  all       convert, then analyze, sharing the same celltype/groupby column
  --help    Show this message

All paths are resolved inside the container, so mount your data directory to
/data and reference files as /data/<filename>.

Examples:
  docker run --rm -v "$PWD/data:/data" liana-ccc-pipeline \
      convert /data/my_data.rds /data/my_data.h5ad RNA cell_type

  docker run --rm -v "$PWD/data:/data" liana-ccc-pipeline \
      analyze --h5ad /data/my_data.h5ad --groupby cell_type --outdir /data/liana_results

  docker run --rm -v "$PWD/data:/data" liana-ccc-pipeline \
      all /data/my_data.rds /data/my_data.h5ad cell_type --outdir /data/liana_results
EOF
}

if [[ $# -eq 0 || "$1" == "--help" || "$1" == "-h" ]]; then
  usage
  exit 0
fi

cmd="$1"; shift

case "$cmd" in
  convert)
    exec Rscript /pipeline/01_convert_seurat_to_h5ad.R "$@"
    ;;
  analyze)
    exec python /pipeline/02_run_liana_ccc_pipeline.py "$@"
    ;;
  all)
    if [[ $# -lt 3 ]]; then
      echo "Usage: all <input.rds> <output.h5ad> <celltype_col> [extra python args...]" >&2
      exit 1
    fi
    input_rds="$1"; output_h5ad="$2"; celltype_col="$3"; shift 3
    Rscript /pipeline/01_convert_seurat_to_h5ad.R "$input_rds" "$output_h5ad" RNA "$celltype_col"
    exec python /pipeline/02_run_liana_ccc_pipeline.py --h5ad "$output_h5ad" --groupby "$celltype_col" "$@"
    ;;
  *)
    echo "Unknown subcommand: $cmd" >&2
    usage
    exit 1
    ;;
esac
