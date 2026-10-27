
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

conda activate star

export TMPDIR="/public/users/wuqingqin/wuqingqin/STAR_CELL/tmp"
mkdir -p "$TMPDIR"

ulimit -n 4096

INPUT_BASE="/public/users/wuqingqin/data_pbmc_cite_241228"

OUTPUT_BASE="/public/users/wuqingqin/wuqingqin/STAR_CELL"

GENOME_DIR="/public/users/wuqingqin/wuqingqin/star_sp/STAR_index"
GTF_FILE="/public/users/wuqingqin/wuqingqin/star_sp/genes.gtf"

WHITELIST="/public/software/cellranger/cellranger-8.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt"


for P in {1..18}; do

    echo "=========================================="
    echo "Starting P${P}"
    echo "=========================================="

    
    MAIN_DIR="${INPUT_BASE}/PBMC-CITE-P${P}-Gene-total-rawdata"

    
    OUTPUT_BASE_DIR="${OUTPUT_BASE}/P${P}"

    
    if [ ! -d "$MAIN_DIR" ]; then
        echo "Input directory not found: $MAIN_DIR"
        echo "Skipping P${P}"
        continue
    fi

    mkdir -p "$OUTPUT_BASE_DIR"

    for SUBDIR in "$MAIN_DIR"/*/; do

        [ -d "$SUBDIR" ] || continue

        DIR_NAME=$(basename "$SUBDIR")

        OUTPUT_SUBDIR="$OUTPUT_BASE_DIR/$DIR_NAME"

        mkdir -p "$OUTPUT_SUBDIR"

        mapfile -d '' -t R1_FILES < <(
            find "$SUBDIR" -maxdepth 1 -type f \
                -name "*R1*.fastq.gz" -print0 | sort -z
        )

        mapfile -d '' -t R2_FILES < <(
            find "$SUBDIR" -maxdepth 1 -type f \
                -name "*R2*.fastq.gz" -print0 | sort -z
        )

        if [ ${#R1_FILES[@]} -eq 0 ] || [ ${#R2_FILES[@]} -eq 0 ]; then
            echo "No R1 or R2 files found: $SUBDIR"
            continue
        fi

        if [ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]; then
            echo "R1/R2 file count mismatch: $SUBDIR"
            continue
        fi

        PAIR_OK=1
        for i in "${!R1_FILES[@]}"; do
            R1_NAME=$(basename "${R1_FILES[$i]}")
            R2_NAME=$(basename "${R2_FILES[$i]}")
            if [ "${R1_NAME/R1/R2}" != "$R2_NAME" ]; then
                PAIR_OK=0
                break
            fi
        done

        if [ "$PAIR_OK" -eq 0 ]; then
            echo "R1/R2 filename pairing mismatch: $SUBDIR"
            continue
        fi

        R1_LIST=$(IFS=,; echo "${R1_FILES[*]}")
        R2_LIST=$(IFS=,; echo "${R2_FILES[*]}")

        echo "------------------------------------------"
        echo "Processing P${P}: $DIR_NAME"
        echo "R1 files: ${#R1_FILES[@]}"
        echo "R2 files: ${#R2_FILES[@]}"
        echo "Output: $OUTPUT_SUBDIR"
        echo "------------------------------------------"

        ##################################################
        # STAR
        ##################################################

        if STAR \
            --runThreadN 40 \
            --genomeDir "$GENOME_DIR" \
            --twopassMode Basic \
            --soloStrand Forward \
            --soloType CB_UMI_Simple \
            --soloCBstart 1 \
            --soloCBlen 16 \
            --soloUMIstart 17 \
            --soloUMIlen 10 \
            --soloBarcodeMate 1 \
            --clip5pNbases 39 0 \
            --soloCBwhitelist "$WHITELIST" \
            --soloCBmatchWLtype 1MM \
            --soloUMIdedup 1MM_Directional_UMItools \
            --readFilesIn "$R1_LIST" "$R2_LIST" \
            --readFilesCommand zcat \
            --outSAMattrRGline ID:RG1 SM:P${P} PL:Illumina LB:Library \
            --limitBAMsortRAM 158718282040 \
            --limitSjdbInsertNsj 2239227 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM XS \
            --sjdbGTFfile "$GTF_FILE" \
            --outFileNamePrefix "$OUTPUT_SUBDIR/${DIR_NAME}"; then

            echo "Finished P${P}: $DIR_NAME"

        else
            echo "ERROR: STAR failed for P${P}: $DIR_NAME"
        fi

    done

    echo "Completed P${P}"

done

echo "=========================================="
echo "All P1-P18 directories processed."
echo "=========================================="

