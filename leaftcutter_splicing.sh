#!/bin/bash
# Step 1: Generate per-cell-type regtools junction extraction scripts from BAM files.
TEMPLATE=$(cat <<'EOF'
#!/bin/bash
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

conda activate leaftcutter

TARGET_DIR="__TARGET__"
FOLDER_NAME=$(basename "$TARGET_DIR")

JUNC_LIST_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/list_file"
JUNC_LIST_FILE="${JUNC_LIST_DIR}/${FOLDER_NAME}_juncfiles.txt"

mkdir -p "$JUNC_LIST_DIR"

if ! command -v regtools &> /dev/null; then
    echo "? regtools 未找到，请确保已正确安装并加入 PATH。" >&2
    exit 1
fi

if [ ! -d "$TARGET_DIR" ]; then
    echo "? 错误: 目标目录 $TARGET_DIR 不存在" >&2
    exit 1
fi

cd "$TARGET_DIR" || { echo "? 错误: 无法进入目录 $TARGET_DIR"; exit 1; }

echo "================================================================"
echo "?? 开始处理目录: $TARGET_DIR (允许覆盖已存在文件)"
echo "?? BAM文件数量: $(find . -maxdepth 1 -type f -name "*.bam" | wc -l) 个"
echo "?? JUNC列表文件: $JUNC_LIST_FILE"
echo "================================================================"

> "$JUNC_LIST_FILE"

extract_junctions() {
    BAM_FILE=$1
    CLEAN_BAM_FILE=$(echo "$BAM_FILE" | sed 's|^\./||')
    JUNC_FILE="${CLEAN_BAM_FILE}.junc"
    ABS_JUNC_PATH="${TARGET_DIR}/${JUNC_FILE}"

    [ -f "$JUNC_FILE" ] && rm -f "$JUNC_FILE"

    echo "?? Extracting junctions from $CLEAN_BAM_FILE..." >&2
    regtools junctions extract -a 8 -m 50 -M 500000 -s FR "$CLEAN_BAM_FILE" -o "$JUNC_FILE" >/dev/null 2>&1

    if [ $? -eq 0 ] && [ -f "$JUNC_FILE" ]; then
        echo "? Junctions extracted successfully from $CLEAN_BAM_FILE." >&2
        echo "$ABS_JUNC_PATH"
        return 0
    else
        echo "??  Failed to extract junctions from $CLEAN_BAM_FILE." >&2
        return 1
    fi
}

export -f extract_junctions
export TARGET_DIR

TEMP_LIST=$(mktemp)

find . -maxdepth 1 -type f -name "*.bam" -print0 | \
    parallel -0 --jobs 15 --env extract_junctions --env TARGET_DIR 'extract_junctions {}' > "$TEMP_LIST"

grep -v '^$' "$TEMP_LIST" | sort -u > "$JUNC_LIST_FILE"
rm -f "$TEMP_LIST"

echo "================================================================"
echo "? 所有处理完成！"
echo "?? JUNC文件路径写入: $JUNC_LIST_FILE"
echo "?? 共记录文件数: $(wc -l < "$JUNC_LIST_FILE")"
echo "================================================================"
EOF
)

ROOT="/public/users/wuqingqin/wuqingqin/STAR_CELL/cell_type"
OUTDIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/job_scripts"
mkdir -p "$OUTDIR"

for subdir in "$ROOT"/*/; do
    absdir=$(realpath "$subdir")
    name=$(basename "$absdir")
    echo "$TEMPLATE" | sed "s|__TARGET__|$absdir|" > "$OUTDIR/run_regtools_${name}.sh"
    chmod +x "$OUTDIR/run_regtools_${name}.sh"
done

echo "? 已生成 $(ls -1 "$OUTDIR" | wc -l) 个脚本到 $OUTDIR"

#!/bin/bash
# Step 2: Generate LeafCutter clustering scripts for separate vaccination states.

JUNCFILES_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/list_file/time"
SCRIPT_DIR="$JUNCFILES_DIR/script"
LEAFCUTTER_SCRIPT="/public/users/wuqingqin/wuqingqin/software/leafcutter/clustering/leafcutter_cluster_regtools.py"

mkdir -p "$SCRIPT_DIR"

for juncfile in "$JUNCFILES_DIR"/*_juncfiles.txt; do
    if [ -f "$juncfile" ]; then
        filename=$(basename "$juncfile")
        prefix="${filename%%_juncfiles.txt}"
        output_dir="$JUNCFILES_DIR/leaftcutter_${prefix}"
        script_file="$SCRIPT_DIR/${prefix}.sh"

        cat << EOF > "$script_file"
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
conda activate leaftcutter

mkdir -p "$output_dir"
cd "$output_dir"

python "$LEAFCUTTER_SCRIPT" -j "$juncfile" -m 50 -o vaccine -l 500000
EOF

        chmod +x "$script_file"
    fi
done

echo "? 所有脚本已生成于：$SCRIPT_DIR"

#!/bin/bash
# Step 3: Prepare LeafCutter splicing phenotype tables for separate states.

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

conda activate leaftcutter

BASE_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/list_file/time"
SCRIPT_PATH="/public/users/wuqingqin/wuqingqin/software/leafcutter/scripts/prepare_phenotype_table.py"
COUNTS_FILE="vaccine_perind.counts.gz"

find "$BASE_DIR" -mindepth 1 -type d | while read -r subdir; do
    counts_file_path="$subdir/$COUNTS_FILE"

    if [ -f "$counts_file_path" ]; then
        subdirname=$(basename "$subdir")

        echo "Processing file: $counts_file_path in subdirectory: $subdir"

        python "$SCRIPT_PATH" "$counts_file_path" -p 100

        if [ $? -eq 0 ]; then
            echo "Successfully processed $counts_file_path."
        else
            echo "Failed to process $counts_file_path."
        fi
    else
        echo "File $counts_file_path does not exist or is not a regular file."
    fi
done

#!/bin/bash
# Step 4: Merge chromosome-wise normalized phenotype files for separate states.

ROOT="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/list_file/time"

CHRS=( {1..22} )

find "$ROOT" -maxdepth 1 -mindepth 1 -type d | while read -r subdir; do
    foldername=$(basename "$subdir")
    out_file="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/snaptron/${foldername}.merged.tsv"

    for chr in "${CHRS[@]}"; do
        infile="$subdir/vaccine_perind.counts.gz.qqnorm_chr${chr}"
        [[ -f $infile ]] && cat "$infile" >> "$out_file"
    done

    echo "Merged 22 chr files -> $out_file"
done

#!/bin/bash
# Step 5: Generate LeafCutter clustering scripts for merged vaccination states.

JUNCFILES_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/list_file"
SCRIPT_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/script"
LEAFCUTTER_SCRIPT="/public/users/wuqingqin/wuqingqin/software/leafcutter/clustering/leafcutter_cluster_regtools.py"

mkdir -p "$SCRIPT_DIR"

for juncfile in "$JUNCFILES_DIR"/*_juncfiles.txt; do
    if [ -f "$juncfile" ]; then
        filename=$(basename "$juncfile")
        prefix="${filename%%_juncfiles.txt}"
        output_dir="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/all_sample/leaftcutter_${prefix}"
        script_file="$SCRIPT_DIR/${prefix}.sh"

        cat << EOF > "$script_file"
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
conda activate leaftcutter

mkdir -p "$output_dir"
cd "$output_dir"

python "$LEAFCUTTER_SCRIPT" -j "$juncfile" -m 50 -o vaccine -l 500000
EOF

        chmod +x "$script_file"
    fi
done

echo "? 所有脚本已生成于：$SCRIPT_DIR"

#!/bin/bash
# Step 6: Prepare LeafCutter splicing phenotype tables for merged states.

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

conda activate leaftcutter

BASE_DIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/all_sample"
SCRIPT_PATH="/public/users/wuqingqin/wuqingqin/software/leafcutter/scripts/prepare_phenotype_table.py"
COUNTS_FILE="vaccine_perind.counts.gz"

find "$BASE_DIR" -mindepth 1 -type d | while read -r subdir; do
    counts_file_path="$subdir/$COUNTS_FILE"

    if [ -f "$counts_file_path" ]; then
        subdirname=$(basename "$subdir")

        echo "Processing file: $counts_file_path in subdirectory: $subdir"

        python "$SCRIPT_PATH" "$counts_file_path" -p 200

        if [ $? -eq 0 ]; then
            echo "Successfully processed $counts_file_path."
        else
            echo "Failed to process $counts_file_path."
        fi
    else
        echo "File $counts_file_path does not exist or is not a regular file."
    fi
done

#!/bin/bash
# Step 7: Merge chromosome-wise normalized phenotype files for merged states.

ROOT="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/all_sample"

CHRS=( {1..22} )

find "$ROOT" -maxdepth 1 -mindepth 1 -type d | while read -r subdir; do
    foldername=$(basename "$subdir")
    out_file="/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/${foldername}.merged.tsv"

    for chr in "${CHRS[@]}"; do
        infile="$subdir/vaccine_perind.counts.gz.qqnorm_chr${chr}"
        [[ -f $infile ]] && cat "$infile" >> "$out_file"
    done

    echo "Merged 22 chr files -> $out_file"

done
