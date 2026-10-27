#!/bin/bash

cov_root="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/snaptron/cov"

script_root="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/snaptron/tqtl_scripts"
mkdir -p "$script_root"

for cov_file in "$cov_root"/*.txt; do
        base=$(basename "$cov_file" .txt)
    script_file="$script_root/run_tensorqtl_interaction_${base}.sh"

    cat <<EOL > "$script_file"

#!/bin/bash

# >>> conda initialize >>>
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

cov_file="$cov_file"
bed_file="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/snaptron/bed/${base}.bed.gz"
vcf_prefix="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/cell_type/vcf/annotation/plink2/${base}"
output_dir="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/snaptron/permute_qtl/${base}"
mkdir -p "\$output_dir"

output_prefix="\$output_dir/${base}"

echo "--> Processing ${base}"

CUDA_VISIBLE_DEVICES=7 python3 -m tensorqtl "\$vcf_prefix" "\$bed_file" "\$output_prefix" \\
    --covariates "\$cov_file" \\
    --mode cis

EOL

    chmod +x "$script_file"
    echo "Generated script: $script_file"
done

echo "All scripts generated."