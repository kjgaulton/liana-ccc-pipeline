#!/usr/bin/env bash
# Dispatcher for the LIANA+ CCC pipeline container.
#
# The container runs as root by default (see Dockerfile) so it can always
# read/execute the baked-in pipeline scripts and write anywhere under /data,
# regardless of how the host's Docker setup maps UIDs/permissions on bind
# mounts (this sidesteps permission mismatches seen with --user tricks on
# some hosts, e.g. NFS-backed mounts, SELinux, or userns-remap). Any files
# created under /data are chowned back to match /data's own owner at the
# end, so you don't end up with root-owned output files.
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
/data and reference files as /data/<filename>. Do NOT pass --user to `docker
run` for this image -- it needs to run as root to avoid host-specific
bind-mount permission issues; output file ownership is fixed up automatically.

Examples:
  docker run --rm -v "$PWD/data:/data" liana-ccc-pipeline \
      convert /data/my_data.rds /data/my_data.h5ad RNA cell_type

  docker run --rm -v "$PWD/data:/data" liana-ccc-pipeline \
      analyze --h5ad /data/my_data.h5ad --groupby cell_type --outdir /data/liana_results

  docker run --rm -v "$PWD/data:/data" liana-ccc-pipeline \
      all /data/my_data.rds /data/my_data.h5ad cell_type --outdir /data/liana_results
EOF
}

# Reassign ownership of everything under /data to match /data's own owner,
# so files created as root inside the container are usable by the host user
# who owns the mounted directory.
fix_ownership() {
  if [[ -d /data ]]; then
    local uid gid
    uid=$(stat -c '%u' /data 2>/dev/null || true)
    gid=$(stat -c '%g' /data 2>/dev/null || true)
    if [[ -n "${uid:-}" && -n "${gid:-}" ]]; then
      chown -R "${uid}:${gid}" /data 2>/dev/null || true
    fi
  fi
}

if [[ $# -eq 0 || "$1" == "--help" || "$1" == "-h" ]]; then
  usage
  exit 0
fi

cmd="$1"; shift
status=0

case "$cmd" in
  convert)
    Rscript /pipeline/01_convert_seurat_to_h5ad.R "$@" || status=$?
    ;;
  analyze)
    python /pipeline/02_run_liana_ccc_pipeline.py "$@" || status=$?
    ;;
  all)
    if [[ $# -lt 3 ]]; then
      echo "Usage: all <input.rds> <output.h5ad> <celltype_col> [extra python args...]" >&2
      fix_ownership
      exit 1
    fi
    input_rds="$1"; output_h5ad="$2"; celltype_col="$3"; shift 3
    Rscript /pipeline/01_convert_seurat_to_h5ad.R "$input_rds" "$output_h5ad" RNA "$celltype_col" || status=$?
    if [[ $status -eq 0 ]]; then
      python /pipeline/02_run_liana_ccc_pipeline.py --h5ad "$output_h5ad" --groupby "$celltype_col" "$@" || status=$?
    fi
    ;;
  *)
    echo "Unknown subcommand: $cmd" >&2
    usage
    fix_ownership
    exit 1
    ;;
esac

fix_ownership
exit "$status"
