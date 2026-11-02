#!/bin/bash
set -euo pipefail

ROOT_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/all_sample"
R_SCRIPT="/public/users/wuqingqin/wuqingqin/mapping_assemble/leafcutter/scripts/leafcutter_ds.R"
EXON_FILE="/public/users/wuqingqin/wuqingqin/leafcutter/hg38_UCSC_exon.txt.gz"
OUT_ROOT="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/result"
N_THREADS=30

mkdir -p "$OUT_ROOT"

sub_dirs=("$ROOT_DIR"/*)

run_one() {
  local dir=$1
  local base=$(basename "$dir")
  local counts="$dir/vaccine_perind_numers.filtered.counts.gz"
  local group_file="$dir/${base}.txt"   # 每个子文件夹对应的 GROUP_FILE
  local prefix="$OUT_ROOT/${base}_vac7_vs_vac0"

  if [ ! -f "$group_file" ]; then
    echo ">>> WARNING: GROUP_FILE not found for $base, skipping..."
    return
  fi

  echo ">>> Processing $base ..."
  /usr/local/bin/Rscript "$R_SCRIPT" \
    --num_threads "$N_THREADS" \
    --exon_file="$EXON_FILE" \
    --output_prefix="$prefix" \
    "$counts" \
    "$group_file"
}

export -f run_one
export R_SCRIPT EXON_FILE OUT_ROOT N_THREADS

printf '%s\n' "${sub_dirs[@]}" | parallel -j 6 run_one {}
