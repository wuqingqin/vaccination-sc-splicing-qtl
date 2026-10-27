
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

# 激活 conda 环境
conda activate peer_env

EXPR_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/splicing"
COV_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/cov"
OUT_ROOT="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/peer"
mkdir -p "$OUT_ROOT"

export COV_DIR OUT_ROOT 

ls ${EXPR_DIR}/*.csv | parallel -j 42 '
    expr_file={}
    base=$(basename "$expr_file" .csv)
    cov_file="${COV_DIR}/${base}.csv"
    out_dir="${OUT_ROOT}/${base}"
    mkdir -p "$out_dir"
    peertool -f "$expr_file" \
             -c "$cov_file" \
             -n 15 \
             --has_rownames \
             --has_header \
             --add_mean \
             -o "$out_dir"
'