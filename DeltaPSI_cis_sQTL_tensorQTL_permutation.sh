#!/bin/bash

bed_dir="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/delta_bed/bed"

cov_dir="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/delta_bed/cov"

vcf_dir="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/cell_type/vcf/annotation_e1_27celltypes/plink2"

out_base="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/delta_bed/permute_SQTL"

script_dir="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/delta_bed/permute_script"

mkdir -p "$script_dir"

for gpu in {0..7}
do
    script="$script_dir/run_gpu${gpu}.sh"

    cat > "$script" <<EOF
#!/bin/bash

__conda_setup="\$('/public/users/wuqingqin/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ \$? -eq 0 ]; then
    eval "\$__conda_setup"
else
    if [ -f "/public/users/wuqingqin/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/public/users/wuqingqin/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/public/users/wuqingqin/miniconda3/bin/:\$PATH"
    fi
fi
unset __conda_setup

conda activate tensorqtl-gpu

EOF

    chmod +x "$script"
done


i=0

for bed_file in "$bed_dir"/leaftcutter_*.merged.bed.gz
do

    [ -e "$bed_file" ] || continue

    filename=$(basename "$bed_file")

    base=${filename%.bed.gz}

    gpu=$((i % 8))

    cov_file="$cov_dir/${base}.txt"
    vcf_prefix="$vcf_dir/${base}"

    output_dir="$out_base/${base}"
    output_prefix="$output_dir/${base}"

    script="$script_dir/run_gpu${gpu}.sh"

    cat >> "$script" <<EOF

echo "============================================================"
echo "--> Processing ${base}"
echo "--> GPU ${gpu}"
echo "============================================================"

cov_file="${cov_file}"
bed_file="${bed_file}"
vcf_prefix="${vcf_prefix}"
output_dir="${output_dir}"
output_prefix="${output_prefix}"

if [ ! -f "\$cov_file" ]; then
    echo "ERROR: covariate file not found: \$cov_file"
    exit 1
fi

if [ ! -f "\${vcf_prefix}.pgen" ]; then
    echo "ERROR: PGEN file not found: \${vcf_prefix}.pgen"
    exit 1
fi

if [ ! -f "\${vcf_prefix}.pvar" ]; then
    echo "ERROR: PVAR file not found: \${vcf_prefix}.pvar"
    exit 1
fi

if [ ! -f "\${vcf_prefix}.psam" ]; then
    echo "ERROR: PSAM file not found: \${vcf_prefix}.psam"
    exit 1
fi

mkdir -p "\$output_dir"

CUDA_VISIBLE_DEVICES=${gpu} python3 -m tensorqtl \
    "\$vcf_prefix" \
    "\$bed_file" \
    "\$output_prefix" \
    --covariates "\$cov_file" \
    --mode cis

echo "--> Finished ${base}"

EOF

    echo "$base -> GPU $gpu"

    i=$((i + 1))

done


echo ""
echo "Generated scripts:"
for gpu in {0..7}
do
    echo "$script_dir/run_gpu${gpu}.sh"
done