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

conda activate plink2

input_dir="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/merged/annotation/plink2"

output_dir="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/GP_PCA"
mkdir -p "$output_dir"

threads=4

ls ${input_dir}/*.pgen | sed 's/\.pgen$//' > prefix_list.txt

cat prefix_list.txt | parallel -j 42 "
    out_prefix=\$(basename {})
    echo '[INFO] Running PCA on' \${out_prefix}
    plink2 \
      --pfile {} \
      --pca 20 \
      --threads ${threads} \
      --out ${output_dir}/\${out_prefix}_pca