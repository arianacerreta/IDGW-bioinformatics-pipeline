#!/usr/bin/env bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J GATK
#SBATCH --cpus-per-task=48
#SBATCH --mem 192G
#SBATCH --mail-user email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

#Load the required modules
source /usr/modules/init/bash
module load gatk/4.1.7.0

#Variables and Paths
REF="/path/to/genome/genomic.fna"   #path ending in genome .fna
RECALBAMDIR="/path/to/dir/with/recal_bams"      #path to directory with recalibrated bams
OUTVCF="/path/to/GATK/CanFam3_1/indiv_vcfs_DPcap500"  #new: directory for individual vcfs; DPcap500 indicates that it will only use a max of 500 reads to make a call
BEDFILE="/path/to/neutral_200loci_CanFam3_1.bed"            #path to where .bed file stored
THREADS=8                                            # number of parallel samples to process

mkdir -p "$OUTVCF"

cd "$RECALBAMDIR" || exit 1      #set working directory to where BAM files are located

ls *.bam | parallel -j "$THREADS" '

SAMPLE=$(basename {} .bam)  # Extract sample name from file

OUTGVCF="$OUTVCF/${SAMPLE}.g.vcf.gz"

# Skip if already exists
if [[ -f "$OUTGVCF" ]]; then
    echo "[INFO] Skipping $SAMPLE, GVCF already exists: $OUTGVCF"
    else
        echo "[INFO] Processing $SAMPLE"
        gatk --java-options "-Xmx20G" HaplotypeCaller \
        -R '"$REF"' \
        -I {} \
        -O '"$OUTVCF"'/${SAMPLE}.g.vcf.gz \
        -ERC GVCF \
        -L '"$BEDFILE"' \
        --interval-padding 50 \
        --max-reads-per-alignment-start 500 \
        --native-pair-hmm-threads 6
fi
'

#after all samples have processed, make map file

MAPFILE="$OUTVCF/sample_map.txt"
ls "$OUTVCF"/*.g.vcf.gz | sed 's!.*/!!;s!.g.vcf.gz!!' | while read S; do
    echo -e "${S}\t$OUTVCF/${S}.g.vcf.gz"
done > "$MAPFILE"

echo "Mapfile made"
