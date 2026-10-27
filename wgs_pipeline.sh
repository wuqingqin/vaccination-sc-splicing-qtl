#!/bin/bash

# Set the source data path.
rawdata='/database/WGS/F22FTSCCKF10011-01_HUMbstuR_2_13JAN2023/'
# Set the output directory.
datadir='/database/WGS_bam/WT'

for sample_name in e3180A

do
# Create the directory with permissions 777 (read, write, and execute for everyone).
mkdir -m 777 ${datadir}/${sample_name}/

cd ${rawdata}/soapnuke/clean/${sample_name}
# Trim adapters.
trim_galore -q 20 --phred33 --stringency 3 --length 20 \
# q: Phred quality score threshold (default: Phred 20); trim sequence with quality scores below the specified threshold.
# Select -phred33 or -phred64 for ASCII+33 or ASCII+64 Phred quality score encoding.
# -s/--stringency sets the required adapter overlap in bases (default: 1).
--paired ${sample_name}_1.fq.gz ${sample_name}_2.fq.gz  \
# Input paired-end reads for the sample.
--gzip -o ${datadir}/${sample_name}
# Write and compress the output files.

# Alignment (fa -> bam).
cd ${datadir}/${sample_name}
bwa-mem2  mem -t 10 -Y -R '@RG\tID:ST-E00318:816:H3JWJCCX2\tPL:illumina\tLB:library\tSM:${sample_name}' \
# Align reads using bwa-mem2.
# Set the number of threads to 10.
# Define the header.
$reference_fasta \
# Specify the reference genome.
*_1.fq.gz *_2.fq.gz| \
# Input adapter-trimmed, quality-controlled paired-end reads.
samtools view -bS| \
# Convert SAM to BAM.
samtools sort -@ 10 -O bam \
# Sort the alignments and specify the number of additional worker threads, excluding the main thread.
# Set the output format to BAM.
 -o ${datadir}/${sample_name}/${sample_name}_sort.bam
# Specify the output filename.
sambamba markdup  -t 8 -r -p --tmpdir=~/mytmp  ${datadir}/${sample_name}/${sample_name}_sort.bam \
${datadir}/${sample_name}/${sample_name}_markdup.bam
# Mark and remove duplicate reads using sambamba.
# The input is a sorted BAM file; the output is a BAM file that overwrites the unsorted BAM file.
TMP_DIR=/home/calsx/Software/temp M=${sample_name}.markdup_metrics.txt
# Save duplicate-removal metrics.
### Modify the BAM header.
picard AddOrReplaceReadGroups \
    -I ${sample_name}.bam \
    -O ${sample_name}_id.bam \
    --RGID ${sample_name} \
    --RGLB lib1 \
    --RGPL illumina \
    --RGPU unit1 \
    --RGSM ${sample_name} \
    --VALIDATION_STRINGENCY LENIENT\
# Use a loop to process multiple samples in one run.
rm *.fq.gz_trimming_report.txt
rm *.markdup_metrics.txt
rm *.fq.gz
rm *_sort.bam
done

samtools index ${sample_name}.bam
# Create an index for the BAM file.
gatk HaplotypeCaller \
-R $reference_fasta  \
# Specify the reference genome.
-I ${sample_name}.bam -max-mnp-distance 0 -O ${datadir}/${sample_name}/${sample_name}.g.vcf.gz -ERC GVCF 
# Input: BAM; output: GVCF.
# MNP: Multiple Nucleotide Polymorphism. Two MNPs separated by 0 or fewer nucleotides will be merged into one MNP.
# Emit Reference Confidence Genotyping includes reference genotypes at variant sites and their associated confidence information.

gatk CombineGVCFs \
   -R $reference_fasta \
   -V ${sample_name}.g.vcf.gz\
   -V ${sample_name}.g.vcf.gz\
   -O cohort.g.vcf.gz
# Use the CombineGVCFs module to combine a small number of samples.

gatk GenomicsDBImport \
-V ${sample_name}.g.vcf.gz \
-V ${sample_name}.g.vcf.gz \
-L chr1 \
--genomicsdb-workspace-path mydata_chr1 \
--tmp-dir tmp\

# Use the GenomicsDBImport module to combine a large number of samples.
# For a large number of samples, split by chromosome and combine the data in a database (mydata_chr1); VCF.GZ information is stored in mydata_chr1, and L specifies the chromosome.

Jointly recalibrate Genotype Quality score	of all samples
gatk GenotypeGVCFs \
 -R $reference_fasta \
 -V gendb://my_database \
 -o $recalibrated_vcf


gatk SelectVariants \
     -R ${sample_name} \
     -V cohort.vcf.gz \
     --select-type-to-include SNP \
     -O cohort.snp.vcf.gz
# Extract SNPs or indels. If -V is a VCF combined through DBImport, supply the VCF.GZ file for the corresponding merged chromosome.



#Variant Quality Score Recalibration(VQSR)	
 gatk VariantRecalibrator /
 -R $reference_fasta /
 -nt 5 /
 -input $recalibrated_vcf /
 -mode SNP /
 -recalFile $vqsr_snp.recal /
 -tranchesFile $vqsr_snp.tranches /
 -rscriptFile $vqsr_snp_plots.R /
 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap
/
 -resource:omni,known=false,training=true,truth=true,prior=12.0 $kg_omni /
 -resource:1000G,known=false,training=true,truth=false,prior=10.0 $kg_snps /
 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp /
 -an QD /
 -an MQ /
 -an FS /
 -an MQRankSum /
 -an ReadPosRankSum /
 -an SOR /
 -an DP /
 -tranche 100.0 /
 -tranche 99.8 /
 -tranche 99.6 /
 -tranche 99.4 /
 -tranche 99.2 /
 -tranche 99.0 /
 -tranche 95.0 /
 -tranche 90.0

gatk VariantRecalibrator /
 -R $reference_fasta /
 -nt 5 /
 -input $recalibrated_vcf /
 -mode INDEL /
 -recalFile $recalibrate_indel.recal /
 -tranchesFile $recalibrate_indel.tranches /
 -rscriptFile $recalibrate_indel_plots.R /
 -resource:mills,known=true,training=true,truth=true,prior=12.0 $kg_mills /
 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp /
 -an QD /
 -an FS /
 -an ReadPosRankSum /
 -an MQRankSum /
 -an SOR /
 -an DP /
 -tranche 100.0 /
 -tranche 99.0 /
 -tranche 95.0 /
 -tranche 92.0 /
 -tranche 90.0 /
 --maxGaussians 4

# known: Whether the resource is a known-variant dataset used to annotate variants.
# training: Whether the resource is used as a training dataset for the VQSR model.
# truth: Whether the resource is a truth set for validating model training; it is also important for automatic parameter selection when VQSR trains the bad model.
# prior: The weight, or prior likelihood, of the resource in VQSR model training, expressed on the Phred scale (for example, 20 corresponds to 0.99).
# an: Annotation.
# nt: Set the number of threads to 5.

 gatk ApplyRecalibration /
 -R $reference_fasta /
 -nt 5 /
 -input $recalibrated_vcf /
 -mode SNP /
 --ts_filter_level 99.80 /
 -recalFile $recalibrate_SNP.recal /
 -tranchesFile $recalibrate_SNP.tranches /
 -o $vqsr_snp_vcf
 
gatk ApplyRecalibration /
 -R $reference_fasta /
 -nt 5 /
 -input $vqsr_snp_vcf /
 -mode INDEL /
 --ts_filter_level 99.0 /
 -recalFile $recalibrate_INDEL.recal /
 -tranchesFile $recalibrate_INDEL.tranches /
 -o $vqsr_snp_indel_vc

#Phase

shapeit4\
    phase_common --input $vqsr_snp_vcf --map b38.gmap.gz \
    --region chr1 --output $phased.vcf.gz --use-PS 0.0001 --thread 24

# b38.gmap.gz is the genetic map; human genetic maps are available at https://github.com/odelaneau/shapeit4/tree/master/maps.
