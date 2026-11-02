#!/bin/bash

input_dir="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/merged/annotation"
output_dir="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/anotation_C_all"
script_dir="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/scripts_C_all"

mkdir -p "$output_dir" "$script_dir"

vep_base="/public/users/qinwanting/Softwares/ensembl-vep/vep"
cache_dir="/public/users/qinwanting/Softwares/ensembl-vep/.vep"
gtf="/public/users/wuqingqin/wuqingqin/STAR/QTL/sqtl_e1_v2/step2/gencode.v43.annotation.gtf.gz"
fasta="/public/users/wuqingqin/wuqingqin/STAR/hg38/hg38_UCSC.fa"
rbp="/public/users/wuqingqin/wuqingqin/STAR/QTL/RBP_feature/merge_data_k562.sorted.bed.gz"

for vcf in "$input_dir"/*.vcf.gz; do
    name=$(basename "$vcf" .vcf.gz)
    script="$script_dir/run_vep_${name}.sh"

    cat <<EOF > "$script"
#!/bin/bash
cd /public/users/wuqingqin/wuqingqin/STAR/QTL/sqtl_e1_v2/step2

$vep_base \\
--tab --fields "Uploaded_variation,Consequence,RBPs" \\
-i $vcf \\
--cache --gtf $gtf \\
--fa $fasta \\
--dir_cache $cache_dir \\
-o $output_dir/${name}.txt \\
--regulatory \\
--custom $rbp,RBPs,bed,overlap,0
EOF

    chmod +x "$script"
    echo "? Created: $script"
done