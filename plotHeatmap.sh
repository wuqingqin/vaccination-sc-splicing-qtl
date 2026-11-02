
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

# 激活 conda 环境
conda activate deeptools

cd  /public/users/wuqingqin/wuqingqin/STAR_CELL/P1


plotHeatmap -m /public/users/wuqingqin/wuqingqin/STAR_CELL/P1/P1_R2_R1_geneBody_matrix.gz \
            --samplesLabel Read1 Read2 \
            --heatmapWidth 4 --heatmapHeight 9 \
            --colorMap Blues Blues \
            -out mergeAll001.all_celltypes.heatmap.final.pdf 