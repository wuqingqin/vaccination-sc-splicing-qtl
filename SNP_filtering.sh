#!/bin/bash

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
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


conda activate vcftools



cd /public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/merged

mkdir -p filter

ls *.vcf.gz | parallel -j 42 '
  fname=$(basename {})
  prefix=${fname%.vcf.gz}
  vcftools --gzvcf {} \
           --max-missing 0.99 \
           --maf 0.05 \
           --hwe 0.00001 \
           --recode --stdout | gzip -c > filter/${prefix}_filter.vcf.gz'
