#!/bin/bash

set -euo pipefail


BASE="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/all_sample/matched_depth_celltypeQ20"

OUT_BASE="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/matched_depth_DAS_result_celltypeQ20"

R_SCRIPT="/public/users/wuqingqin/wuqingqin/mapping_assemble/leafcutter/scripts/leafcutter_ds.R"

EXON_FILE="/public/users/wuqingqin/wuqingqin/leafcutter/hg38_UCSC_exon.txt.gz"

N_THREADS=30

N_JOBS=6

mkdir -p "$OUT_BASE"


# ============================================================
# cell type
# ============================================================

run_one() {

    local dir="$1"
    local seed="$2"

    local base
    base=$(basename "$dir")

    local counts="${dir}/vaccine_perind_numers.filtered.downsampled.counts.gz"

    local group_file="${dir}/${base}.txt"

    local outdir="${OUT_BASE}/seed_${seed}"

    mkdir -p "$outdir"

    local prefix="${outdir}/${base}_vac7_vs_vac0"


    # --------------------------------------------------------

    # --------------------------------------------------------

    if [ ! -f "$counts" ]; then

        echo "Missing counts:"
        echo "$counts"

        return 1
    fi


    if [ ! -f "$group_file" ]; then

        echo "Missing group file:"
        echo "$group_file"

        return 1
    fi


    # --------------------------------------------------------
    # paired sample
    # --------------------------------------------------------

    local n_e1
    local n_e2

    n_e1=$(awk '$2=="e1"{n++} END{print n+0}' "$group_file")

    n_e2=$(awk '$2=="e2"{n++} END{print n+0}' "$group_file")


    if [ "$n_e1" -ne "$n_e2" ]; then

        echo "ERROR: unequal e1/e2 samples:"
        echo "${base}: e1=${n_e1}, e2=${n_e2}"

        return 1
    fi


    # --------------------------------------------------------
    # --------------------------------------------------------

    if [ -f "${prefix}_effect_sizes.txt" ] && \
       [ -f "${prefix}_cluster_significance.txt" ]; then

        echo "Already finished:"
        echo "seed_${seed} ${base}"

        return 0
    fi


    echo
    echo "============================================================"
    echo "Seed      : ${seed}"
    echo "Cell type : ${base}"
    echo "e1        : ${n_e1}"
    echo "e2        : ${n_e2}"
    echo "============================================================"


    # ========================================================
    # LeafCutter DAS
    #
    # ========================================================

    /usr/local/bin/Rscript "$R_SCRIPT" \
        --num_threads "$N_THREADS" \
        --exon_file="$EXON_FILE" \
        --output_prefix="$prefix" \
        "$counts" \
        "$group_file"
}


export -f run_one

export BASE
export OUT_BASE
export R_SCRIPT
export EXON_FILE
export N_THREADS


# ============================================================
# seed 1–5
# ============================================================

for seed in {1..5}
do

    echo
    echo "############################################################"
    echo "START seed_${seed}"
    echo "############################################################"

    SEED_DIR="${BASE}/seed_${seed}"

    find "$SEED_DIR" \
        -maxdepth 1 \
        -type d \
        -name "leaftcutter_*" \
        -print0 \
    | parallel \
        -0 \
        -j "$N_JOBS" \
        run_one {} "$seed"

    echo
    echo "FINISHED seed_${seed}"

done


echo
echo "============================================================"
echo "ALL DAS ANALYSES FINISHED"
echo "Output: ${OUT_BASE}"
echo "============================================================"