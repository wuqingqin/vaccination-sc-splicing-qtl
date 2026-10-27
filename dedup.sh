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

conda activate star

base_dir="/public/users/wuqingqin/wuqingqin/STAR_CELL"

for project_dir in $base_dir/*/ ; do
    if [ -d "$project_dir" ]; then
        for bam_file in $project_dir/*GeneAligned.sortedByCoord.out.bam ; do
            if [ -f "$bam_file" ]; then
                base_name=$(basename "$bam_file" -GeneAligned.sortedByCoord.out.bam)
                
                output_dir="$project_dir"
                output_bam="$output_dir/${base_name}_dedup.bam"
                metrics_file="$output_dir/marked_dup_metrics.${base_name}.txt"

                picard MarkDuplicates -Xmx64g \
                    REMOVE_DUPLICATES=true \
                    BARCODE_TAG=CB \
                    MOLECULAR_IDENTIFIER_TAG=UB \
                    I="$bam_file" \
                    O="$output_bam" \
                    M="$metrics_file"
            fi
        done
    fi
done