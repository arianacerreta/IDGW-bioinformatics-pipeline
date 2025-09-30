#!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J summ_stats
#SBATCH --cpus-per-task=16
#SBATCH --mail-user example@email.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#Script to calculate summary statistics for panel performance

#Load the required modules
source /usr/modules/init/bash
module load bcftools
module load samtools

#Define variables
MERGED="/path/to/final/all.gatk.genotyped.unzipped.vcf" #unzipped
BAMDIR="/path/to/recalibrated/bams"
TARGETS_BED="/path/to/target.bed" #BED file with aplicon target regions
THREADS=16 #number of threads to use; should equal cpus-per-task in header

#Create output directory
OUTDIR="/path/to/summ_stats/$(date +%y%m%d)_genome_GATK" #where you want your output directory created

mkdir -p "$OUTDIR"

#1 calculate and output depth per locus
echo "Calculating depth per locus..."
bcftools query -f '%CHROM\t%POS[\t%DP]\n' "$MERGED" > "$OUTDIR/depth_per_locus.tsv"

#2 calculate genotype completness per sample
echo "Calulating genotype completeness per sample..."
bcftools query -f '[%SAMPLE\t%GT\n]' "$MERGED" | \
awk '$2 != "./." { count[$1]++ } END { for (s in count) print s, count[s] }' \
> "$OUTDIR/completeness_per_sample.tsv"

# 3. Percent reads mapping to target regions

cd "$BAMDIR" || exit 1

echo "Calculating % reads on target..."
echo -e "Sample\tTotal_Reads\tReads_on_Target\tPercent_on_Target" > "$OUTDIR/reads_on_target.tsv"

ls *.bam | parallel -j "$THREADS" '
    SAMPLE=$(basename {} .bam)
    TOTAL=$(samtools view -c {})
    ON_TARGET=$(samtools view -c -L '"$TARGETS_BED"' {})
    if [ "$TOTAL" -gt 0 ]; then
      PCT=$(echo "scale=4; $ON_TARGET/$TOTAL*100" | bc)
    else
      PCT="NA"
    fi
    echo -e "${SAMPLE}\t${TOTAL}\t${ON_TARGET}\t${PCT}"
' >> "$OUTDIR/reads_on_target.tsv"

echo "All summary statistics written to: $OUTDIR"
