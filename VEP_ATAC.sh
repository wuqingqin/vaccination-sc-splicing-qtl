#!/bin/bash

outdir="/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/scripts"
mkdir -p "$outdir"

bed_dir="/public/users/qinwanting/data/eQTL_enrichment/functional_annotation/open_chromatin"

peaks=(
"pDC"
"NK_CD56dim_CD16high"
"NK_CD56bright_CD16low"
"Memory_B"
"MAIT"
"CD8+_Tn"
"CD8+_Teff_Tm"
"CD4+_Treg"
"CD4+_Tn"
"CD4+_Teff_Tm"
"CD1C+_cDC"
"CD16+_Mono"
"CD14+_Mono"
)

for p in "${peaks[@]}"; do
    script_path="${outdir}/run_vep_${p}.sh"
    echo "Generating script: $script_path"

    cat > "$script_path" <<EOF
#!/bin/bash
cd ${bed_dir}

echo "Running VEP for ${p}..."

 /public/users/qinwanting/Softwares/ensembl-vep/vep \\
 --tab --fields "Uploaded_variation,Consequence,ATAC_Peaks" \\
 -i /public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/merged/annotation/leaftcutter_${p}.merged_filter_annotated_final.vcf.gz \\
 --cache \\
 --gtf /public/users/wuqingqin/wuqingqin/STAR/QTL/sqtl_e1_v2/step2/gencode.v43.annotation.gtf.gz \\
 --fa /public/users/wuqingqin/wuqingqin/STAR/hg38/hg38_UCSC.fa \\
 --dir_cache /public/users/qinwanting/Softwares/ensembl-vep/.vep \\
 -o /public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/anotation/leaftcutter_${p}.merged.txt \\
 --regulatory \\
 --custom ATAC_peaks_${p}.bed.gz,ATAC_Peaks,bed,overlap,0 

EOF

    chmod +x "$script_path"
done