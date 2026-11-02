#!/bin/bash

# >>> conda initialize >>>
__conda_setup="$('/public/users/wuqingqin/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/public/users/wuqingqin/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/public/users/wuqingqin/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/public/users/wuqingqin/miniconda3/bin/:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate LDSC

#!/bin/bash


SUMSTATS_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/hg_trans/gwas_summary_gz"
ANNOT_ROOT="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/hg_trans/ldscore/anno_file"
OUT_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/hg_trans/ldscore/enrichment"

WEIGHT_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/iqtl_coloc/ldsc/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC."
BASELINE_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/iqtl_coloc/ldsc/baseline/baselineLD."
FRQ_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/iqtl_coloc/ldsc/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC."

mkdir -p "$OUT_DIR"


for sumstats in ${SUMSTATS_DIR}/*.sumstats.gz; do

    trait=$(basename "$sumstats" .sumstats.gz)
    echo "==== Trait: $trait ===="

    for cell_dir in ${ANNOT_ROOT}/*; do

        cell=$(basename "$cell_dir")
        echo "  --> Cell: $cell"

        ref_ld="${cell_dir}/${cell}_chr,${BASELINE_DIR}"

        out_prefix="${OUT_DIR}/${trait}_${cell}_PIP_enrichment"

        ldsc.py \
          --h2 "$sumstats" \
          --ref-ld-chr "$ref_ld" \
          --w-ld-chr "$WEIGHT_DIR" \
          --overlap-annot \
          --frqfile-chr "$FRQ_DIR" \
          --out "$out_prefix" \
          --print-coefficients

    done
done