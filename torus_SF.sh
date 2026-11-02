#!/bin/bash

WORK_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/torus_RBP"
OUTPUT_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/torus_RBP/reported_qtl_prepare"
TORUS_PATH="/public/users/qinwanting/Softwares/torus-master/src/torus"
SCRIPT_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/torus_RBP/torus_scripts"

mkdir -p "$SCRIPT_DIR"
cd "$WORK_DIR" || exit 1

for summ_file in *.merged.summ.gz; do
    base_name="${summ_file%.merged.summ.gz}"

    gene_pos_file="${base_name}.merged.gene.pos.gz"
    snp_pos_file="${base_name}.merged.snp.pos.gz"

    if [[ -f "$gene_pos_file" && -f "$snp_pos_file" ]]; then
        
        for annot_file in ${base_name}.merged_part*.an.gz; do
            if [[ -f "$annot_file" ]]; then
                part_name=$(basename "$annot_file" .an.gz)
                part_name=${part_name/_filter_annotated_final/}
                script_file="${SCRIPT_DIR}/run_torus_${part_name}.sh"
                est_file="${OUTPUT_DIR}/${part_name}.enrichment.est"

                cat > "$script_file" << EOF
#!/bin/bash
# Torus analysis script for $part_name
# Generated on: $(date)

echo "=== Running TORUS for ${part_name} ==="
cd "$WORK_DIR" || exit 1
mkdir -p "$OUTPUT_DIR"

$TORUS_PATH \\
    -d "$summ_file" \\
    -gmap "$gene_pos_file" \\
    -smap "$snp_pos_file" \\
    -annot "$annot_file" \\
    -est > "$est_file"

if [ \$? -eq 0 ]; then
    echo "Finished successfully: ${part_name}"
    echo "  Output: ${est_file}"
else
    echo "Failed: ${part_name}"
    exit 1
fi
EOF

                chmod +x "$script_file"
                echo "Generated script: $script_file"
            fi
        done
    else
        echo "Warning: Missing gene/snp file for $base_name, skipping."
    fi
done

echo "✅ All scripts generated in: $SCRIPT_DIR"
echo "You can run all with:"
echo "  for script in $SCRIPT_DIR/*.sh; do sbatch \$script; done"