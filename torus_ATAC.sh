#!/bin/bash

# Set paths
WORK_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/torus"
OUTPUT_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/torus/reported_qtl_prepare"
TORUS_PATH="/public/users/qinwanting/Softwares/torus-master/src/torus"
SCRIPT_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/torus/torus_scripts"  # Script output directory

# Create script directory
mkdir -p "$SCRIPT_DIR"

# Change to working directory
cd "$WORK_DIR"

# Find all summ.gz files and generate corresponding scripts
for summ_file in *.merged.summ.gz; do
    # Extract base filename (remove .summ.gz suffix)
    base_name="${summ_file%.merged.summ.gz}"
    
    # Construct other filenames
    gene_pos_file="${base_name}.merged.gene.pos.gz"
    snp_pos_file="${base_name}.merged.snp.pos.gz"
    annot_file="${base_name}.merged.an.gz"
    
    # Check if all required files exist
    if [[ -f "$gene_pos_file" && -f "$snp_pos_file" && -f "$annot_file" ]]; then
        # Script filename
        script_file="${SCRIPT_DIR}/run_torus_${base_name}.sh"
        
        # Construct output filenames
        prior_file="${OUTPUT_DIR}/${base_name}.enrichment.prior"
        est_file="${OUTPUT_DIR}/${base_name}.enrichment.est"
        
        # Generate script content
        cat > "$script_file" << EOF
#!/bin/bash
# Torus analysis script - $base_name
# Generated at: $(date)

echo "Starting processing cell type: $base_name"
echo "Working directory: $WORK_DIR"
echo "Output directory: $OUTPUT_DIR"

# Change to working directory
cd "$WORK_DIR"
mkdir -p "$OUTPUT_DIR"

# Run torus command
$TORUS_PATH \\
    -d "$summ_file" \\
    -gmap "$gene_pos_file" \\
    -smap "$snp_pos_file" \\
    -annot "$annot_file" \\
    -est \\
    -dump_prior "$prior_file" \\
    > "$est_file"

# Check command execution status
if [ \$? -eq 0 ]; then
    echo "Successfully completed: $base_name"
    echo "Output files:"
    echo "  - $prior_file"
    echo "  - $est_file"
else
    echo "Processing failed: $base_name"
    exit 1
fi

echo "Completion time: $(date)"
EOF
        
        # Add execute permission to the script
        chmod +x "$script_file"
        
        echo "Generated script: $script_file"
    else
        echo "Warning: Incomplete files, skipping $base_name"
    fi
done

echo "Script generation completed!"
echo "All scripts saved in: $SCRIPT_DIR"
echo "You can use the following command to run in batch:"
echo "  for script in $SCRIPT_DIR/*.sh; do sbatch \$script; done"
