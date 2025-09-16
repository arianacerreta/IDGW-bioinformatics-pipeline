#!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J BQSR
#SBATCH --cpus-per-task=4
#SBATCH --mail-user email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

#Load the required modules
source /usr/modules/init/bash
module load gatk/4.1.7.0

#Define variables
BAMDIR="/path/to//RG_aligned_sorted/Clu10kTash"    #directory with RG and sorted .bam and .bai files
VCF="/path/to/reference_variants.vcf.gz"     #path to VCF file
GENOME="/path/to/genome/genome.fna"    #path to reference genome .fna
OUTDIR="/path/to/new/directory/recalibration/genomename" #will be made
OUTDIRBAM="/path/to/new/directory/recalibration/genomename/recal_bams" #will be made

mkdir -p "$OUTDIR" "$OUTDIRBAM"

#jobs to run in parallel
JOBS=8

cd "$BAMDIR" || exit 1      #set working directory to where BAM files are located

#Apply BQSR per sample
ls *.bam | parallel -j "$JOBS" '
    SAMPLE=$(basename {} .bam)  # Extract sample name from file

    gatk BaseRecalibrator \
        -R '"$GENOME"' \
        -I {} \
        --known-sites '"$VCF"' \
        -O '"$OUTDIR"'/${SAMPLE}.recal.table \
        > '"$OUTDIR"'/${SAMPLE}.recal.log 2>&1

    gatk ApplyBQSR \
    -R '"$GENOME"' \
    -I {} \
    --bqsr-recal-file '"$OUTDIR"'/${SAMPLE}.recal.table \
    -O '"$OUTDIRBAM"'/${SAMPLE}.recal.bam

    gatk BuildBamIndex -I '"$OUTDIRBAM"'/${SAMPLE}.recal.bam
'
